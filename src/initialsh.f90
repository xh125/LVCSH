module initialsh
  use kinds     ,only   : dp,dpc
  use constants ,only   : maxlen,tpi,K_B_Ryd,ryd2eV
  use elph2     ,only   : epcq,nq=>nqtotf
  use modes     ,only   : nmodes
  implicit none
  
  contains
  
  
  ! ref: 1 S. Fernandez-Alberti et al., The Journal of Chemical Physics 137 (2012) 
  ! ref: 2. HuangKun <固体物理> (9-29) (9-31)
  subroutine init_eh_stat(laser,init_ik,cband,vband)
    use parameters,only : init_ikx,init_iky,init_ikz,init_kx,init_ky,init_kz
    use epwcom,only : nkf1,nkf2,nkf3
    use readepw,only : E_nk
    use elph2,only  : vmef,ibndmin,ibndmax,nbndfst,nkf  !vmef(3,nbndsub,nbndsub,nkf)
    use getmcvk,only: W_cvk,get_Mcvk !W_cvk(nbndfst,nbndfst,nkf)
    use constants,only :ryd2eV
    implicit none
    
    logical,intent(in)    :: laser
    integer,intent(out)   :: init_ik
    integer,intent(inout) :: cband,vband
    integer :: ik,ivband,icband
    real(kind=dp) :: flagr,flagd,W_cvk_all
    real(kind=dp) :: obsorbEn
    !real(kind=dp),allocatable :: W_cvk(nbndfst,nbndfst,nkf) !光激发下的跃迁几率大小
    if (laser) then
      call get_Mcvk()
      W_cvk_all = SUM(W_cvk)
      call random_number(flagr)
      flagr = flagr*W_cvk_all
      flagd = 0.0
      outter:do ik=1,nkf
               do ivband=1,nbndfst-1
                 do icband=ivband+1,nbndfst
                   flagd = flagd+W_cvk(icband,ivband,ik)
                   if (flagr <= flagd) then
                     init_ik    = ik
                     cband = icband
                     vband = ivband
                     exit outter
                   endif
                 enddo
               enddo
      enddo outter
      obsorbEn = (E_nk(cband,init_ik)-E_nk(vband,init_ik))*ryd2eV
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
  subroutine init_dynamical_variable(ph_Q,laser,c_nk,v_nk,ee,pp,ww_e,ww_h)
    use parameters, only : init_cband,init_vband,init_ik
    use hamiltonian,only : H0_nk,H_nk,set_H_nk,calculate_eigen_energy_state,nefre
    use elph2,only       : nbndfst,nktotf
    use epwcom,only      : kqmap
    use surfacehopping,only : iesurface,ihsurface,convert_diabatic_adiabatic,p_nk
    use io,only : stdout
    use lasercom,only : w_laser
    use readepw,only : E_nk
    implicit none
    real(kind=dp),intent(in) :: ph_Q(nmodes,nq)
    logical,intent(in) :: laser
    real(kind=dp),intent(out) :: ee(nefre),pp(nefre,nefre)
    complex(kind=dpc),intent(out) :: c_nk(nbndfst,nktotf),v_nk(nbndfst,nktotf),ww_e(nefre),ww_h(nefre)
    integer :: iefre
    real(kind=dp) :: flagr,flagd
    real(kind=dp) :: en_eh
    
    call init_eh_stat(laser,init_ik,init_cband,init_vband) 
    c_nk = 0.0d0
    v_nk = 0.0d0
    c_nk(init_cband,init_ik) = 1.0d0
    v_nk(init_vband,init_ik) = 1.0d0
      
    call set_H_nk(ph_Q,H_nk)
    call calculate_eigen_energy_state(nktotf,nbndfst,H_nk,ee,pp)
    p_nk = reshape(pp,(/ nbndfst,nktotf,nefre /))
    call convert_diabatic_adiabatic(p_nk,c_nk,ww_e)
    call convert_diabatic_adiabatic(p_nk,v_nk,ww_h)
   
    call random_number(flagr)
    flagd = 0.0d0
    do iefre = 1,nefre
      flagd = flagd + p_nk(init_cband,init_ik,iefre)**2
      if(flagr <= flagd) then
        iesurface = iefre
        exit
      endif
    enddo
    
    call random_number(flagr)
    flagd = 0.0d0
    do iefre = 1,nefre
      flagd = flagd + p_nk(init_vband,init_ik,iefre)**2
      if(flagr <= flagd) then
        ihsurface = iefre
        exit
      endif
    enddo    
    
    en_eh = (ee(iesurface)-ee(ihsurface))*ryd2eV
    write(stdout,"(/,5X,A)") "In the laser obsorbtion,the inital excited state as follow:"
    write(stdout,"(5X,A22,F12.7,A3)")  "Laser centred energy :",w_laser*ryd2eV," eV"
    
    write(stdout,"(/,5X,A,I5)")"In diabatic base,electron excited :init_ik=",init_ik
    write(stdout,"(5X,A,I5)")  "Electron in the conductor band:initi_cband=",init_cband
    write(stdout,"(5X,A,I5)")  "Hole     in the valence   band:initi_vband=",init_vband
    write(stdout,"(5X,A22,F12.7,A3)")  "The energy of elctron:",E_nk(init_cband,init_ik)*ryd2eV," eV"
    write(stdout,"(5X,A22,F12.7,A3)")  "The energy of hole   :",E_nk(init_vband,init_ik)*ryd2eV," eV"
    write(stdout,"(5X,A22,F12.7,A3)")  "The energy of exciton:",&
                                        (E_nk(init_cband,init_ik)-E_nk(init_vband,init_ik))*ryd2eV," eV"
    
    write(stdout,"(/,5X,A,I5)")"In adiabatic base,the elctron and hole state as follow"
    write(stdout,"(5X,A,I5)")  "Electron in the energy surface:iesurface=",iesurface
    write(stdout,"(5X,A,I5)")  "Hole in the energy surface    :ihsurface=",ihsurface
    write(stdout,"(5X,A22,F12.7,A3)")  "The energy of elctron:",ee(iesurface)*ryd2eV," eV"
    write(stdout,"(5X,A22,F12.7,A3)")  "The energy of hole   :",ee(ihsurface)*ryd2eV," eV"
    write(stdout,"(5X,A22,F12.7,A3)")  "The energy of exciton:",&
                                        (ee(iesurface)-ee(ihsurface))*ryd2eV," eV"    
    
  end subroutine init_dynamical_variable  
  
  
end module initialsh 