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
                               
    if(ieband_min==0 .and. ieband_max==0) then
      ieband_min = ncbmin
      ieband_max = ibndmax
    endif
    if(ihband_min==0 .and. ihband_max==0) then
      ihband_min = ibndmin
      ihband_max = nvbmax
    endif
    
    if(lelecsh) then
      if(ieband_min < ibndmin) then
        write(stdout,"(A,I9,1X,A3,1X,I9)") "Error! The parameter: ieband_min need to set between ", ibndmin,"and",ibndmax
        stop
      endif
      if(ieband_max > ibndmax) then
        write(stdout,"(A,I9,1X,A3,1X,I9)") "Error! The parameter: ieband_max need to set between ", ibndmin,"and",ibndmax
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
  subroutine init_eh_KSstat(lelecsh,lholesh,llaser,init_ik,init_eband,init_hband,e_en,h_en)
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
    real(kind=dp),intent(out) :: e_en,h_en
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

    obsorbEn = (etf(init_eband,init_ik)-etf(init_hband,init_ik))*ryd2eV
    e_en = etf(init_eband,init_ik)
    h_en = -1.0*etf(init_hband,init_ik)
    write(stdout,"(/5X,A)") "Initial eh_KSstate:  "
    write(stdout,"(5X,A3,I5,1X,A10,3(F12.6,1X),A2)") "ik=",init_ik, "coord.: ( ",(xkf(ipol,init_ik),ipol=1,3)," )"
    write(stdout,"(5X,A11,I5,A21,F12.5,A3)") "init_hband=",init_hband," Initial hole Energy:",&
                                              etf(init_hband,init_ik)*ryd2eV," eV"
    write(stdout,"(5X,A11,I5,A21,F12.5,A3)") "init_eband=",init_eband," Initial elec Energy:",&
                                              etf(init_eband,init_ik)*ryd2eV," eV"      
    write(stdout,"(5X,A18,F12.5,A3)")"elec-hole energy= ",obsorbEn," eV"
    
  end subroutine init_eh_KSstat
  
  subroutine init_stat_diabatic(init_ik,init_band,iband_min,nband,nk,c_nk)
    implicit none
    integer,intent(in)    :: init_ik
    integer,intent(inout) :: init_band
    integer,intent(in)    :: iband_min
    integer,intent(in)    :: nband,nk
    complex(kind=dpc),intent(out) :: c_nk(nband,nk) 
    
    init_band = init_band - iband_min + 1
    c_nk = czero
    c_nk(init_band,init_ik) = cone
    
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
    if(womiga == 0.0) then
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
  !ref : 3 "<Phonons Theory and Experiments I Lattice Dynamics and Models of Interatomic Forces by Dr. Peter Brüesch (auth.) (z-lib.org).pdf>."
  subroutine init_normalmode_coordinate_velocity(nmodes,nq,w,T,l_ph_quantum,ph_Q,ph_P)
    use kinds,only   : dp,dpc
    use randoms,only : gaussian_random_number
    !use parameters,only : lit_ephonon
    use surfacecom,only : nqv,E_ph_CA_sum,E_ph_QA_sum,eps_acustic
    use elph2,only : iminusq
    use dynamics,only : test_xv
    use io,only : stdout
    use constants,only : ryd2eV
    implicit none
    integer,intent(in)           :: nmodes,nq
    real(kind=dp),intent(in)     :: T
    real(kind=dp),intent(in)     :: w(nmodes,nq)
    logical,intent(in)           :: l_ph_quantum
    complex(kind=dpc),intent(out):: ph_Q(nmodes,nq),ph_P(nmodes,nq)
    ! complex normal coordinates ph_Q, and 
    ! ph_Q(v,-q)=ph_Q(v,q)*    (2.48)
    ! H = T + V =0.5*[|ph_P|^^2+w_qv^^2*|ph_Q|^^2]   (2.54)
    ! ph_P(v,-q)=ph_P(v,q)*    (2.57)
    ! ph_P(v,q) = d(T)/d(ph_P(v,q)*) (2.55)
    
    real(kind=dp) :: womiga,theta1,theta2
    complex(kind=dpc) :: cplx_tmp1,cplx_tmp2
    real(kind=dp) :: E_ph_class,E_ph_quantum
    integer :: iq,imode
    
    nqv = 0.0

    do iq=1,nq
      do imode=1,nmodes
        womiga = w(imode,iq)
        if(womiga >= eps_acustic) then
          nqv(imode,iq) = bolziman(womiga,T)+0.5
        endif
      enddo
    enddo
    
    E_ph_CA_sum   = 0.0
    E_ph_QA_sum   = 0.0
        
    do iq=1,nq
      if( iq <= iminusq(iq) ) then
        ! ph_Q(v,-q)=ph_Q(v,q)*
        ! ph_P(v,-q)=ph_P(v,q)*
        ! Ph_P(v,q) = d(ph_Q(v,q))/dt 
        do imode=1,nmodes
          womiga = w(imode,iq)
          
          call random_number(theta1)
          if(iq==iminusq(iq)) then
            ! for ph_Q and ph_P is real
            theta1 = 0.0
            theta2 = 0.0
          else
            theta1 = theta1 * tpi
            theta2 = theta1 - tpi/4.0
          endif
          
          cplx_tmp1 = cos(theta1)*cone+sin(theta1)*ci
          cplx_tmp2 = cos(theta2)*cone+sin(theta2)*ci
          
          if(womiga < eps_acustic) then
            ph_Q(imode,iq)= czero
            ph_P(imode,iq)= czero
          else
            E_ph_class   = K_B_Ryd*T  ! IN class 
            E_ph_CA_sum  = E_ph_CA_sum + E_ph_class
            if(iminusq(iq)/=iq) E_ph_CA_sum  = E_ph_CA_sum + E_ph_class
            E_ph_quantum = nqv(imode,iq)*womiga ! In Quantum
            E_ph_QA_sum  = E_ph_QA_sum + E_ph_quantum	
            if(iminusq(iq)/=iq) E_ph_QA_sum  = E_ph_QA_sum + E_ph_quantum
            if(l_ph_quantum) then
              ph_Q(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_quantum)/womiga)*cplx_tmp1
              ph_P(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_quantum))*cplx_tmp2
            else
              ph_Q(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_class)/womiga)*cplx_tmp1
              ph_P(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_class))*cplx_tmp2
            endif
          endif
          
        enddo
      else 
      !  ph_Q(v,-q)=ph_Q(v,q)*   ! ph_P(v,-q)=ph_P(v,q)*
        ph_Q(:,iq) = CONJG(ph_Q(:,iminusq(iq)))
        ph_P(:,iq) = CONJG(ph_P(:,iminusq(iq)))
      endif
    enddo
    
    !call test_xv(nmodes,nq,ph_Q,ph_P)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !% Write phonon energy information         %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    write(stdout,"(/5X,A51,F11.5,A2)") &
          "The temperature of the non-adiabatic dynamica is : ",T," K"
    write(stdout,"(5X,A49,F11.5,A4)") &
          "The average energy of phonon(quantum): <SUM_phE>=",E_ph_QA_sum*ryd2eV," eV."
    write(stdout,"(5X,A49,F11.5,A4)") &
          "The average energy of phonon(class)  : <SUM_phE>=",E_ph_CA_sum*ryd2eV," eV."
    
    if(l_ph_quantum) then
      write(stdout,"(/5X,A)") "The phonon dynamica set as quantum"
    else
      write(stdout,"(/5X,A)") "The phonon dynamica set as classical."
    endif
    write(stdout,"(5X,A38,F11.5,A4,A9,F11.5,A4)") &
    "The initial energy of phonon: SUM_phT=",0.5*SUM(ABS(ph_P)**2)*ryd2eV," eV",&
    " SUM_phU=",0.5*SUM(w**2*ABS(ph_Q)**2)*ryd2eV," eV"
    write(stdout,"(5X,A38,F11.5,A4)") &
    "The initial energy of phonon: SUM_phE=",0.5*SUM(ABS(ph_P)**2+w**2*ABS(ph_Q)**2)*ryd2eV," eV."    
    
    
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