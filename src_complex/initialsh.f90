module initialsh
  use kinds     ,only   : dp,dpc
  use constants ,only   : maxlen,tpi,K_B_Ryd,ryd2eV,ryd2mev,cone,czero,ci
  implicit none
  
  contains
  
  subroutine set_subband(lelecsh,lholesh,ieband_min,ieband_max,ihband_min,ihband_max)
    use readepw,only : icbm,nvbmax,ncbmin
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
        ieband_min = ncbmin
        ieband_max = ibndmax
      endif
      
      if(ieband_min < ibndmin) then
        write(stdout,"(A,I9,A3,I9)") "Error! The parameter: ieband_min need to set between ", ibndmin,"and",ibndmax
        stop
      endif
      if(ieband_max > ibndmax) then
        write(stdout,"(A,I9,A3,I9)") "Error! The parameter: ieband_max need to set between ", ibndmin,"and",ibndmax
        stop
      endif
      if(ieband_max<ieband_min) then
        write(stdout,"(A)") "Error! The parameter: ieband_max need to set larger than ieband_min "
        stop
      endif
      
      write(stdout,"(/,5X,A)") "The electron are non-diabatic dynamica in valence band."
      WRITE(stdout,'(/14x,a,i5,2x,a,i5)') 'ieband_min = ', ieband_min, 'ieband_max = ', ieband_max
      

      if(ieband_min < ncbmin ) then
        write(stdout,"(/5X,A)") "Electron can dynamical to the conductor band."
        write(stdout,"(A)") "Calculation the electron-hole Carrier recombination. "
      endif
      
    endif
    
    if(lholesh) then
      if(ihband_min==0 .and. ihband_max==0) then
        ihband_min = ibndmin
        ihband_max = nvbmax
      endif

      if(ihband_min < ibndmin) then
        write(stdout,"(A,I9,A3,I9)") "Error! The parameter: ihband_min need to set between ", ibndmin,"and",ibndmax
        stop
      endif
      if(ihband_max > ibndmax) then
        write(stdout,"(A,I9,A3,I9)") "Error! The parameter: ihband_max need to set between ", ibndmin,"and",ibndmax
        stop
      endif
      if(ihband_max<ihband_min) then
        write(stdout,"(A)") "Error! The parameter: ihband_max need to set larger than ihband_min "
        stop
      endif

      write(stdout,"(/,5X,A)") "The hole are non-diabatic dynamica in valance band."
      WRITE(stdout,'(/14x,a,i5,2x,a,i5)') 'ihband_min = ', ihband_min, 'ihband_max = ', ihband_max      
      
      if(ihband_max > nvbmax ) then
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
    use readepw,only : etf,icbm,xkf
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
    integer :: ivbm,ipol
    
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
                   flagd = flagd+W_cvk(ieband,ihband,ik)
                   if (flagr <= flagd) then
                     init_ik    = ik
                     init_eband = ieband
                     init_hband = ihband
                     exit outter
                   endif
                 enddo
               enddo
      enddo outter
      obsorbEn = (etf(init_eband,init_ik)-etf(init_hband,init_ik))*ryd2eV
      write(stdout,"(/5X,A)") "Initial eh_KSstate:  "
      write(stdout,"(5X,A3,I5,1X,A10,3(F12.6,1X),A2)") "ik=",init_ik, "coord.: ( ",(xkf(ipol,init_ik),ipol=1,3)," )"
      write(stdout,"(5X,A11,I5,A21,F12.5,A3)") "init_hband=",init_hband," Initial hole Energy:",&
                                                etf(init_hband,init_ik)*ryd2eV," eV"
      write(stdout,"(5X,A11,I5,A21,F12.5,A3)") "init_eband=",init_eband," Initial elec Energy:",&
                                                etf(init_eband,init_ik)*ryd2eV," eV"      
      write(stdout,"(5X,A18,F12.5,A3)")"elec-hole energy= ",obsorbEn," eV"
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
  subroutine init_normalmode_coordinate_velocity(nmodes,nq,w,T,l_ph_quantum,ph_Q,ph_P)
    use kinds,only   : dp,dpc
    use randoms,only : gaussian_random_number
    use parameters,only : lit_ephonon
    use surfacecom,only : E_ph_CA_sum,E_ph_QA_sum
    use elph2,only : iminusq
    implicit none
    integer,intent(in)           :: nmodes,nq
    real(kind=dp),intent(in)     :: T
    real(kind=dp),intent(in)     :: w(nmodes,nq)
    logical,intent(in)           :: l_ph_quantum
    complex(kind=dpc),intent(out):: ph_Q(nmodes,nq),ph_P(nmodes,nq)
    
    real(kind=dp) :: womiga,theta
    complex(kind=dpc) :: cplx_tmp
    real(kind=dp) :: E_ph_class,E_ph_quantum
    integer :: iq,imode
    
    
    E_ph_CA_sum   = 0.0
    E_ph_QA_sum   = 0.0
    do iq=1,nq
      if(iminusq(iq)>=iq) then
        ! ph_Q(v,q)=ph_Q(v,-q)*
        do imode=1,nmodes
          womiga = w(imode,iq)
          call random_number(theta)
          if(iq==iminusq(iq)) theta = 0.0
          theta = theta * tpi
          cplx_tmp = cos(theta)*cone+sin(theta)*ci
          
          if(womiga*ryd2mev <= lit_ephonon) then
            ph_Q(imode,iq)= czero
            ph_P(imode,iq)= czero
          else
            E_ph_class   = K_B_Ryd*T  ! IN class 
            E_ph_CA_sum  = E_ph_CA_sum + E_ph_class
            if(iminusq(iq)/=iq) E_ph_CA_sum  = E_ph_CA_sum + E_ph_class
            E_ph_quantum = (bolziman(womiga,T)+0.5)*womiga ! In Quantum
            E_ph_QA_sum  = E_ph_QA_sum + E_ph_quantum	
            if(iminusq(iq)/=iq) E_ph_QA_sum  = E_ph_QA_sum + E_ph_quantum
            if(l_ph_quantum) then
              ph_Q(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_quantum)/womiga)*cplx_tmp
              ph_P(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_quantum))*cplx_tmp
              if(iminusq(iq)/=iq) then
                ph_Q(imode,iminusq(iq)) = CONJG(ph_Q(imode,iq))
                ph_P(imode,iminusq(iq)) = CONJG(ph_P(imode,iq))
              endif
            else
              ph_Q(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_class)/womiga)*cplx_tmp
              ph_P(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_class))*cplx_tmp
              if(iminusq(iq)/=iq) then
                ph_Q(imode,iminusq(iq)) = CONJG(ph_Q(imode,iq))
                ph_P(imode,iminusq(iq)) = CONJG(ph_P(imode,iq))
              endif
            endif
          endif
          
        enddo
      endif
    enddo
    
  end subroutine init_normalmode_coordinate_velocity
  
  !=============================================!
  != init dynamical surface                    =!
  !=============================================!  
  subroutine init_surface(nfre,nfre_sh,ww,isurface)
    implicit none
    
    integer , intent(in) :: nfre,nfre_sh
    complex(kind=dpc),intent(in)    :: ww(nfre)
    integer , intent(out) :: isurface
    
    integer :: ifre
    real(kind=dp) :: flagr,flagd,flagsum
   
		flagsum = 0.0
		do ifre=1,nfre_sh
			flagsum = flagsum + ww(ifre)*CONJG(ww(ifre))
		enddo
	 
	 
    call random_number(flagr)
		flagr = flagr * flagsum
    flagd = 0.0d0
    do ifre = 1,nfre_sh
      flagd = flagd + ww(ifre)*CONJG(ww(ifre))
      if(flagr <= flagd) then
        isurface = ifre
        exit
      endif
    enddo
    
  end subroutine init_surface 
  
  
end module initialsh 