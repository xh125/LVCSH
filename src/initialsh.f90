module initialsh
  use kinds     ,only   : dp,dpc
  use constants ,only   : maxlen,tpi,K_B_Ryd,ryd2eV
  implicit none
  
  contains
  
  subroutine set_subband(lelecsh,lholesh,ieband_min,ieband_max,ihband_min,ihband_max)
    use readepw,only : icbm
    use elph2  ,only : ibndmin,ibndmax 
    use io,     only : stdout
    implicit none
    logical , intent(in) :: lelecsh
    logical , intent(in) :: lholesh
    integer , intent(inout) :: ieband_min,&
                               ieband_max,&
                               ihband_min,&
                               ihband_max
    if(lelecsh) then
      if(ieband_min==0 .and. ieband_max==0) then
        ieband_min = icbm
        ieband_max = ibndmax
      endif
      
      write(stdout,"(/,5X,A)") "The electron are non-diabatic dynamica in valence band."
      WRITE(stdout,'(/14x,a,i5,2x,a,i5)') 'ieband_min = ', ieband_min, 'ieband_max = ', ieband_max
      
      if(ieband_min < ibndmin) then
        write(stdout,"(/5X,A)") "Error! Need to reset ieband_min,for the reason ieband_min<ibndmin."
        write(stdout,"(A)") "Check the ibndmin in epw.out"
      endif
      if(ieband_max > ibndmax) then
        write(stdout,"(/5X,A)") "Error! Need to reset ieband_max,for the reason ieband_max>ibndmax."
        write(stdout,"(A)") "Check the ibndmax in epw.out"
      endif
      if(ieband_min < icbm ) then
        write(stdout,"(/5X,A)") "Electron can dynamical to the conductor band."
        write(stdout,"(A)") "Calculation the electron-hole Carrier recombination. "
      endif
      
    endif
    
    if(lholesh) then
      if(ihband_min==0 .and. ihband_max==0) then
        ihband_min = ibndmin
        ihband_max = icbm - 1
      endif

      write(stdout,"(/,5X,A)") "The hole are non-diabatic dynamica in valance band."
      WRITE(stdout,'(/14x,a,i5,2x,a,i5)') 'ihband_min = ', ihband_min, 'ihband_max = ', ihband_max      
      
      if(ihband_min < ibndmin) then
        write(stdout,"(/5X,A)") "Error! Need to reset ihband_min,for the reason ihband_min<ibndmin."
        write(stdout,"(A)") "Check the ibndmin in epw.out"
      endif
      if(ihband_max > ibndmax) then
        write(stdout,"(/5X,A)") "Error! Need to reset ihband_max,for the reason ihband_max>ibndmax."
        write(stdout,"(A)") "Check the ibndmax in epw.out"
      endif
      if(ihband_max >= icbm ) then
        write(stdout,"(/5X,A)") "Hole can dynamical to the conductor band."
        write(stdout,"(A)") "Calculation the electron-hole Carrier recombination. "
      endif    
    endif
    
    
    
  end subroutine set_subband
  
  ! ref: 1 S. Fernandez-Alberti et al., The Journal of Chemical Physics 137 (2012) 
  ! ref: 2. HuangKun <固体物理> (9-29) (9-31)
  subroutine init_eh_KSstat(lelecsh,lholesh,llaser,init_ik,init_eband,init_hband)
    use parameters,only : init_ikx,init_iky,init_ikz,init_kx,init_ky,init_kz
    use epwcom,only : nkf1,nkf2,nkf3
    use readepw,only : etf,icbm
    use surfacecom,only : ieband_min,ieband_max,ihband_min,ihband_max,c_e_nk,c_h_nk
    use elph2,only  : vmef,ibndmin,ibndmax,nbndfst,nkf  !vmef(3,nbndsub,nbndsub,nkf)
    use getwcvk,only: W_cvk !(W_cvk(icbm:ieband_max,ihband_min:ivbm,nkf)
    use constants,only :ryd2eV
    use io,only : stdout
    implicit none
    
    logical , intent(in)  :: lelecsh
    logical , intent(in)  :: lholesh    
    logical,intent(in)    :: llaser
    integer,intent(out)   :: init_ik
    integer,intent(inout) :: init_eband,init_hband
    integer :: ik,ihband,ieband
    real(kind=dp) :: flagr,flagd,W_cvk_all
    real(kind=dp) :: obsorbEn
    integer :: ivbm
    
    ivbm = icbm - 1
    !allocate(W_cvk(icbm:ieband_max,ihband_min:ivbm,nkf) !光激发下的跃迁几率大小
    if (llaser) then
      W_cvk_all = SUM(W_cvk)
      call random_number(flagr)
      flagr = flagr*W_cvk_all
      flagd = 0.0
      outter:do ik=1,nkf
               do ihband=ihband_min,ivbm
                 do ieband=icbm,ieband_max
                   flagd = flagd+W_cvk(ihband,ieband,ik)
                   if (flagr <= flagd) then
                     init_ik    = ik
                     init_eband = ieband
                     init_hband = ihband
                     exit outter
                   endif
                 enddo
               enddo
      enddo outter
      obsorbEn = (etf(init_eband,2*init_ik-1)-etf(init_hband,2*init_ik-1))*ryd2eV
    else
      init_ikx = get_ik(init_kx,nkf1)
      init_iky = get_ik(init_ky,nkf2)
      init_ikz = get_ik(init_kz,nkf3)
      init_ik  =  (init_ikx - 1) * nkf2 * nkf3 + (init_iky - 1) * nkf3 + init_ikz
      if(init_eband < icbm .or. init_eband > ieband_max) then
        write(stdout,"(/5X,A29,I5)") "Wrong! The init_eband set as:",init_eband
        write(stdout,"(5X,A29,I5,A16,I5)") "The init_eband need to set : ",icbm,"<= init_eband <=",ieband_max
      endif
      if(init_hband < ihband_min .or. init_hband > ivbm) then
        write(stdout,"(/5X,A29,I5)") "Wrong! The init_hband set as:",init_hband
        write(stdout,"(5X,A29,I5,A16,I5)") "The init_hband need to set : ",ihband_min,"<= init_hband <=",ivbm
      endif      
    endif
    
    !if(lelecsh) then 
    !  init_eband = init_eband - ieband_min + 1
    !  c_e_nk = 0.0d0
    !  c_e_nk(init_eband,init_ik) = 1.0d0
    !endif
    !if(lholesh) then
    !  init_hband = init_hband - ihband_min + 1
    !  c_h_nk = 0.0d0
    !  c_h_nk(init_hband,init_ik) = 1.0d0
    !endif
    
  end subroutine init_eh_KSstat
  
  subroutine init_stat_diabatic(init_ik,init_band,iband_min,nband,nk,c_nk)
    implicit none
    integer,intent(in)    :: init_ik
    integer,intent(inout) :: init_band
    integer,intent(in)    :: iband_min
    integer,intent(in)    :: nband,nk
    complex(kind=dpc),intent(out) :: c_nk(nband,nk) 
    
    init_band = init_band - iband_min + 1
    c_nk = 0.0d0
    c_nk(init_band,init_ik) = 1.0d0
    
  end subroutine init_stat_diabatic
  
  function get_ik(kx,nkx)
    use kinds,only : dp
    implicit none
    real(kind=dp) :: kx
    integer       :: nkx
    integer       :: get_ik
    kx = MOD(kx,1.0)
    if(kx<=0.0) kx = kx + 1.0
    get_ik = Anint(kx*nkx)+1
    if(get_ik<1) then
      get_ik = get_ik +nkx
    elseif(get_ik>nkx) then
      get_ik = get_ik - nkx
    endif
    
  end function  


  !ref : 1 G. GRIMvall, <The electron-phonon interaction in metals by Goran Grimvall (z-lib.org).pdf> 1981),  
  !    : (3.24)  
  function bolziman(womiga,temp)
    use io,only :stdout
    use constants,only : K_B_Ryd
    implicit none
    real(kind=dp),intent(in)::womiga,temp
    real(kind=dp) :: bolziman
    if(womiga == 0) then
      write(stdout,*) "womiga == 0,phonon error"
      stop
    endif
    bolziman=1.0/(exp(womiga/(K_B_Ryd*temp))-1.0)
    !<nb>=1/(exp{hbar*w/kbT}-1)
  end function   

  !==============================================!
  != init Normal mode coordinate and  velocitie =!
  !==============================================!
  !ref : 1 G. GRIMvall, <The electron-phonon interaction in metals by Goran Grimvall (z-lib.org).pdf> 1981),  
  !    :  (3.17) (3.20) (3.24)
  !ref : 2 HuangKun 《固体物理》 (3-44) (3-45)
  subroutine init_normalmode_coordinate_velocity(nmodes,nq,w,T,ph_Q,ph_P)
    use kinds,only   : dp
    use randoms,only : gaussian_random_number

    implicit none
    integer,intent(in)       :: nmodes,nq
    real(kind=dp),intent(in) :: T
    real(kind=dp),intent(in) :: w(nmodes,nq)  
    real(kind=dp),intent(out):: ph_Q(nmodes,nq),ph_P(nmodes,nq)
    
    real(kind=dp) :: womiga
    real(kind=dp) :: E_ph_class,E_ph_quantum
    integer :: iq,imode
    
    do iq=1,nq
      do imode=1,nmodes
        womiga = w(imode,iq)
        E_ph_class   = K_B_Ryd*T  ! IN class 
        E_ph_quantum = (bolziman(womiga,T)+0.5)*womiga ! In Quantum
        ph_Q(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_quantum)/womiga)
        ph_P(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_quantum))      
        
        if(iq==1 .and. imode <=3) then
          ph_Q(imode,iq)=0.0
          ph_P(imode,iq)=0.0
        endif
        
      enddo
    enddo
    
  end subroutine init_normalmode_coordinate_velocity
  
  !=============================================!
  != init dynamical surface                    =!
  !=============================================!  
  subroutine init_surface(nfre,ww,isurface)
    implicit none
    
    integer , intent(in) :: nfre
    complex(kind=dpc),intent(in)    :: ww(nfre)
    integer , intent(out) :: isurface
    
    integer :: ifre
    real(kind=dp) :: flagr,flagd
   
    call random_number(flagr)
    flagd = 0.0d0
    do ifre = 1,nfre
      flagd = flagd + ww(ifre)*CONJG(ww(ifre))
      if(flagr <= flagd) then
        isurface = ifre
        exit
      endif
    enddo
    
    
    !en_eh = (ee(iesurface)-ee(ihsurface))*ryd2eV
    !write(stdout,"(/,5X,A)") "In the laser obsorbtion,the inital excited state as follow:"
    !write(stdout,"(5X,A22,F12.7,A3)")  "Laser centred energy :",w_laser*ryd2eV," eV"
    !
    !write(stdout,"(/,5X,A,I5)")"In diabatic base,electron excited :init_ik=",init_ik
    !write(stdout,"(5X,A,I5)")  "Electron in the conductor band:initi_init_eband=",init_init_eband
    !write(stdout,"(5X,A,I5)")  "Hole     in the valence   band:initi_init_hband=",init_init_hband
    !write(stdout,"(5X,A22,F12.7,A3)")  "The energy of elctron:",E_nk(init_init_eband,init_ik)*ryd2eV," eV"
    !write(stdout,"(5X,A22,F12.7,A3)")  "The energy of hole   :",E_nk(init_init_hband,init_ik)*ryd2eV," eV"
    !write(stdout,"(5X,A22,F12.7,A3)")  "The energy of exciton:",&
    !                                    (E_nk(init_init_eband,init_ik)-E_nk(init_init_hband,init_ik))*ryd2eV," eV"
    !
    !write(stdout,"(/,5X,A,I5)")"In adiabatic base,the elctron and hole state as follow"
    !write(stdout,"(5X,A,I5)")  "Electron in the energy surface:iesurface=",iesurface
    !write(stdout,"(5X,A,I5)")  "Hole in the energy surface    :ihsurface=",ihsurface
    !write(stdout,"(5X,A22,F12.7,A3)")  "The energy of elctron:",ee(iesurface)*ryd2eV," eV"
    !write(stdout,"(5X,A22,F12.7,A3)")  "The energy of hole   :",ee(ihsurface)*ryd2eV," eV"
    !write(stdout,"(5X,A22,F12.7,A3)")  "The energy of exciton:",&
    !                                    (ee(iesurface)-ee(ihsurface))*ryd2eV," eV"    
    
  end subroutine init_surface 
  
  
end module initialsh 