!----------------------------------------------------------------------------!
! Copyright (C) 2020 WeiZhuang's group                                       !
! This file is distributed under the terms of the                            !
! GNU General Public License. See the file `License'                         !
! in the root directory of the present distribution,                         !
! or http://www.gnu.org/copyleft/gpl.txt .                                   !
!                                                                            !
!----------------------------------------------------------------------------!

program lvcsh
  !=====================================================================================================!
  != this program is used to simulate carrier transport with the self-consistent surface hopping method=!
  !=====================================================================================================!
  != of materials. We use a liner vibronic coupling model with electron-phonon couplings and           =!
  !=====================================================================================================!
  !=  system interaction with environment by system-bath interactions                                  =!     
  !=====================================================================================================!
  != Last updata 2021-04.14 Version:0.1.2                                                              =!   
  != Developed by XieHua at department of physic, USTC;xh125@mail.ustc.edu.cn                          =!
  !=====================================================================================================!
  !! Author: HuaXie
  !! Version: v0.1.1
  !! License: GNU
  !!=================================================================================================  
  use mkl_service
  use omp_lib
  use constants,only      : maxlen
  use environments,only   : environment_start,mkl_threads,&
                            set_mkl_threads,lsetthreads
  use readinput,only      : get_inputfile
  use readscf,only        : readpwscf_out
  use readphout,only      : readph_out
  use readepw,only        : readepwout
  use parameters, only    : lreadscfout,scfoutname,lreadphout,phoutname,epwoutname,temp,&
                            nsnap,nstep,dt,inputfilename,init_ik,llaser,methodsh
  use hamiltonian,only    : H_nk,set_H_nk,set_H0_nk,calculate_eigen_energy_state
  use randoms,only        : init_random_seed
  use initialsh,only      : init_normalmode_coordinate_velocity,init_dynamical_variable
  use surfacehopping,only : iaver,isnap,istep,naver,phQ,phP,phQ0,phP0,e,p,e0,p0,d,d0,ge,ge1,gh,gh1,&
                            allocatesh,celec_nk,chole_nk,w_e,w_h,w0_e,w0_h,iesurface,ihsurface,&
                            calculate_nonadiabatic_coupling,convert_diabatic_adiabatic,&
                            calculate_hopping_probability,calculate_sumg_pes,minde_e,minde_h,&
                            sumg0_e,sumg0_h,sumg1_e,sumg1_h,nonadiabatic_transition
  use elph2,only          : wf,nqtotf,nktotf,nbndfst
  use disp,only           : ph_configuration
  use modes,only          : nmodes
  use io      ,only       : stdout,io_time,time1,time2
  use dynamics,only       : rk4_nuclei,rk4_electron_diabatic,lelec,lhole,ADD_BATH_EFFECT
  !use dynamica
  implicit none
  !===============!
  != preparation =!
  !===============!
  real(kind=8) t0,t1
  call cpu_time(t0)  
  
  call environment_start( 'LVCSH' )
  call get_inputfile(inputfilename)
  if(lreadscfout) call readpwscf_out(scfoutname)
  if(lreadphout) call readph_out(phoutname)
  call readepwout(epwoutname)
  call set_H0_nk()
  call init_random_seed()
  if(lsetthreads) call set_mkl_threads(mkl_threads)
  call allocatesh(nmodes)
  
  
  !==========================!
  != loop over realizations =!
  !==========================!  
  do iaver=1,naver
    write(stdout,'(a,I4,a)') '###### iaver=',iaver,' ######'    
    !==================!
    != initialization =!
    !==================! 
    !!Get the initial normal mode coordinate phQ and versity phP
    call init_normalmode_coordinate_velocity(wf,temp,phQ,phP)
    !!得到初始电子和空穴的状态
    call init_dynamical_variable(phQ,llaser,celec_nk,chole_nk,e,p,w_e,w_h)
    call calculate_nonadiabatic_coupling(nmodes,e,p,d)
    phQ0=phQ; phP0=phP; e0=e; p0=p; d0=d; w0_e=w_e; w0_h=w_h
    

    !=======================!
    != loop over snapshots =!
    !=======================!

    do isnap=1,nsnap
      do istep=1,nstep
        
        time1   = io_time()  
        !==============================!
        != update phQ,phP,c,e,p,d,w,g =!
        !==============================!
        
        !use rk4 to calculate the dynamical of phonon normal modes
        !update phQ,phP in time t+dt
        call rk4_nuclei(p0,phQ,phP,dt)
        !electron and hole wave function is propagated in diabatic representation
        !update c_e,c_h in time t+dt
        call rk4_electron_diabatic(phQ0,celec_nk,dt,lelec)
        call rk4_electron_diabatic(phQ0,chole_nk,dt,lhole)
        
        ! hamiltonian in time t+dt
        ! update H_nk
        call set_H_nk(phQ,H_nk)
        ! update e,p in time t+dt
        call calculate_eigen_energy_state(nktotf,nbndfst,H_nk,e,p)
        ! update d in time t+dt
        call calculate_nonadiabatic_coupling(nmodes,e,p,d)
        call convert_diabatic_adiabatic(p,celec_nk,w_e)
        call convert_diabatic_adiabatic(p,chole_nk,w_h)
        call calculate_hopping_probability(iesurface,w0_e,phP0,d0,dt,ge,ge1)
        call calculate_hopping_probability(ihsurface,w0_h,phP0,d0,dt,gh,gh1)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !% CALCULATE SUMG0,SUMG1,MINDE and CHANGE POTENTIAL ENERGY SURFACE %!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!        
        call calculate_sumg_pes(sumg0_e,sumg1_e,w0_e,w_e,ge1,ge,iesurface,minde_e)
        call calculate_sumg_pes(sumg0_h,sumg1_h,w0_h,w_h,gh1,gh,ihsurface,minde_h)
        
        
        call nonadiabatic_transition(lelec,iesurface,E0,P0,d0,ge,w_e,phP)
        call nonadiabatic_transition(lhole,ihsurface,E0,P0,d0,ge,w_h,phP)        
   
        !===================!
        != add bath effect =!
        !===================!        
        call add_bath_effect(E0,P0,d0,dt,phQ,phP)

        !============================!
        != reset dynamical variable =!
        !============================!

        phQ0=phQ; phP0=phP; e0=e; p0=p; d0=d; w0_e=w_e; w0_h=w_h
        time2   = io_time()
        write(stdout,"(A5,1X,I10,1X,F15.6,1X,A10)") 'Step=',istep,(time2-time1),"seconds"        
      enddo
      
      !=====================!
      != store information =!
      !=====================!    

!      pes_elec(0,isnap,iaver)=E_e(isurface_e)
!      pes_hole(0,isnap,iaver)=E_h(isurface_h)
!      do ibasis=1,nbasis
!        csit_elec(ibasis,isnap)=csit_elec(ibasis,isnap)+abs(c_e(ibasis))**2
!        csit_hole(ibasis,isnap)=csit_hole(ibasis,isnap)+abs(c_h(ibasis))**2
!        wsit_elec(ibasis,isnap)=wsit_elec(ibasis,isnap)+abs(w_e(ibasis))**2
!        wsit_hole(ibasis,isnap)=wsit_hole(ibasis,isnap)+abs(w_h(ibasis))**2
!        psit_elec(ibasis,isnap)=psit_elec(ibasis,isnap)+p_e(ibasis,isurface_e)**2
!        psit_hole(ibasis,isnap)=psit_hole(ibasis,isnap)+p_h(ibasis,isurface_h)**2
!        pes_elec(ibasis,isnap,iaver)=e_e(ibasis)
!        pes_hole(ibasis,isnap,iaver)=e_h(ibasis)
!      enddo      
!      
!      do ifreem =1 , nfreem
!        !!平均核构型
!        xsit(ifreem,isnap)=xsit(ifreem,isnap)+phQ(ifreem)
!        !!平均动能
!        ksit(ifreem,isnap)=ksit(ifreem,isnap)+0.5d0*phP(ifreem)**2
!           
    enddo
  enddo
!  
!  csit_elec=csit_elec/naver
!  csit_hole=csit_hole/naver
!  wsit_elec=wsit_elec/naver
!  wsit_hole=wsit_hole/naver
!  psit_elec=psit_elec/naver
!  psit_hole=psit_hole/naver
!  xsit=xsit/naver*Au2ang*dsqrt(au2amu)
!  ksit=ksit/naver*Au2eV
!
!  !====================!
!  != save information =!
!  !====================!
!  call saveresult()
  !call cpu_time(t1)
  !write(6,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'
  !write(stdout,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'  
  
  close(stdout)
  
  stop
  
end program lvcsh


