module initialsh
  use kinds     ,only   : dp,dpc
  use constants ,only   : maxlen,tpi,K_B_Ryd
  use parameters
  implicit none
  
  contains
  
  subroutine init_eh_stat(laser,init_ik,cband,vband)
    use parameters,only : init_ikx,init_iky,init_ikz
    use epwcom,only : nkf1,nkf2,nkf3
    use elph2,only  : vmef,ibndmin,ibndmax,nbndfst  !vmef(3,nbndsub,nbndsub,nkf)
    use lasercom,only: W_cvk !W_cvk(nbndfst,nbndfst,nkf)
    implicit none
    
    logical,intent(in)    :: laser
    integer,intent(out)   :: init_ik
    integer,intent(inout) :: cband,vband
    !real(kind=dp),allocatable :: W_cvk(nbndfst,nbndfst,nkf) !光激发下的跃迁几率大小
    if (laser) then
      
      
      
      init_ik    = 1
      init_cband = 1
      init_vband = 2
    else
      init_ikx = get_ik(init_kx,nkf1)
      init_iky = get_ik(init_ky,nkf2)
      init_ikz = get_ik(init_kz,nkf3)
      init_ik  =  (init_ikx - 1) * nkf2 * nkf3 + (init_iky - 1) * nkf3 + init_ikz
    endif
    
  end subroutine init_eh_stat
  
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
  subroutine init_normalmode_coordinate_velocity(nq,nmodes,w,T,ph_Q,ph_P)
    use kinds,only   : dp
    use randoms,only : gaussian_random_number

    implicit none
    integer      ,intent(in) :: nq,nmodes
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
        ph_Q(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_quantum/womiga))
        ph_P(imode,iq) = gaussian_random_number(0.0d0,dsqrt(E_ph_quantum))      
      
      enddo
    enddo
    
  end subroutine init_normalmode_coordinate_velocity
  
  !=============================================!
  != init dynamical varibale                   =!
  !=============================================!  
  subroutine init_dynamical_variable(nq,nmodes,ph_Q,laser,c_nk,v_nk,ee,pp,ww)
    use parameters, only : init_cband,init_vband,init_ik
    use hamiltonian,only : H0_nk,H_nk,set_H_nk,calculate_eigen_energy_state
    implicit none
    integer,intent(in) :: nq,nmodes
    real(kind=dp),intent(in) :: ph_Q(nmodes,nq)
    logical,intent(in) :: laser
    real(kind=dp),intent(out) :: ee(nbndfst*nktotf),pp(nbndfst*nktotf,nbndfst*nktotf)
    complex(kind=dpc),intent(out) :: c_nk(nbndfst,nktotf),v_nk(nbndfst,nktotf),ww(nbndfst*nktotf)
    integer :: iefre
    real(kind=dp) :: flagr,flagd
    
    call init_eh_stat(laser,init_ik,init_cband,init_vband) 
    c_nk = 0.0d0
    v_nk = 0.0d0
    c_nk(init_cband,init_ik) = 1.0d0
    v_nk(init_vband,init_ik) = 1.0d0
      
    call set_H_nk(nq,nmodes,ph_l,nbndfst,nktotf,epcq,kqmap,H0_nk,H_nk)
    call calculate_eigen_energy_state(nktotf,nbndfst,H_nk,ee,pp)
    call convert_diabatic_adiabatic(pp,c_nk,ww)
    
    call random_number(flagr)
    flagd = 0.0d0
    do iefre = 1,nefre
      flagd = flagd + pp((init_ik-1)*nbndfst+init_ik,iefre)**2
      if(flagr <= flagd) then
        isurface = iefre
        exit
      endif
    enddo
    
  end subroutine init_dynamical_variable  
  
  
end module initialsh 