module initialsh
  use kinds     ,only   : dp,dpc
  use constants ,only   : maxlen,tpi,K_B_Ryd,ryd2eV
  use elph2     ,only   : epcq,nq=>nqtotf
  use modes     ,only   : nmodes
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
      if(ieband_min < ibndmin) then
        write(stdout,"(/5X,A)") "Error! Need to reset ieband_min,for the reason ieband_min<ibndmin."
        write(stdout,"(A)") "Check the ibndmin in epw.out"
      endif
      if(ieband_max > ibndmax) then
        write(stdout,"(/5X,A)") "Error! Need to reset ieband_max,for the reason ieband_max>ibndmax."
        write(stdout,"(A)") "Check the ibndmax in epw.out"
      endif
      if(ieband_min < icbm ) then
        write(stdout,"(/5X,A)") "Electron can dynamical to the valence band."
        write(stdout,"(A)") "Calculation the electron-hole Carrier recombination. "
      endif
    endif
    
    if(lholesh) then
      if(ihband_min==0 .and. ihband_max==0) then
        ihband_min = ibndmin
        ihband_max = icbm - 1
      endif
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
  subroutine init_eh_stat_diabatic(lelecsh,lholesh,llaser,init_ik,init_eband,init_hband)
    use parameters,only : init_ikx,init_iky,init_ikz,init_kx,init_ky,init_kz
    use epwcom,only : nkf1,nkf2,nkf3
    use readepw,only : etf,icbm
    use surfacecom,only : ieband_min,ieband_max,ihband_min,ihband_max,c_e_nk,c_h_nk
    use elph2,only  : vmef,ibndmin,ibndmax,nbndfst,nkf  !vmef(3,nbndsub,nbndsub,nkf)
    use getwcvk,only: W_cvk !W_cvk(nbndfst,nbndfst,nkf)
    use constants,only :ryd2eV
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
      !call get_Mcvk()
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
    endif
    
    if(lelecsh) then 
      init_eband = init_eband - ieband_min + 1
      c_e_nk = 0.0d0
      c_e_nk(init_eband,init_ik) = 1.0d0
    endif
    if(lholesh) then
      init_hband = init_hband - ihband_min + 1
      c_h_nk = 0.0d0
      c_h_nk(init_hband,init_ik) = 1.0d0
    endif
    
  end subroutine init_eh_stat_diabatic
  
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
  subroutine init_normalmode_coordinate_velocity(w,T,ph_Q,ph_P)
    use kinds,only   : dp
    use randoms,only : gaussian_random_number

    implicit none
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
  != init dynamical varibale                   =!
  !=============================================!  
  subroutine init_dynamical_variable(nband,nk,P_nk,c_nk,ww,isurface)
                  
    !use elph2,only       : nbndfst,nktotf
    !use epwcom,only      : kqmap
    use surfacehopping,only : convert_diabatic_adiabatic
    !use io,only : stdout
    !use readepw,only : E_nk,etf
    
    implicit none
    integer , intent(inout) :: nband,nk
    real(kind=dp) , intent(inout)    :: P_nk(nband,nk,nband*nk) 
    complex(kind=dpc),intent(inout)  :: c_nk(nband,nk)
    complex(kind=dpc),intent(out)    :: ww(nband*nk)
    integer , intent(out) :: isurface
    
    integer :: nfre
    integer :: ifre
    real(kind=dp) :: flagr,flagd
    real(kind=dp) :: en_eh
    
    nfre = nband*nk
    
    call convert_diabatic_adiabatic( nband,nk,p_nk,c_nk,ww )
   
    call random_number(flagr)
    flagd = 0.0d0
    do ifre = 1,nfre
      flagd = flagd + ww(ifre)**2
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
    
  end subroutine init_dynamical_variable  
  
  
end module initialsh 