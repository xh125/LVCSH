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
  != Last updata 2020-12.30 Version:0.1.1                                                              =!   
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
  use readinput,only      : get_inputfile,treat_parameters
  use readscf,only        : readpwscf_out
  use readphout,only      : readph_out
  use readepw,only        : readepwout
  use parameters, only    : scfoutname,phoutname,epwoutname,temp
  use hamiltonian,only    : set_H0_nk,e,p
  use randoms,only        : init_random_seed
  use surfacehopping,only : iaver,naver,phQ,phP,&
                            allocatesh,init_normalmode_coordinate_velocity,c_nk,w,&
                            init_dynamical_variable
  use elph2,only          : wf,nqtotf
  use phdisp,only         : ph_configuration,ph_lqv,ph_l,ph_p
  use modes,only          : nmodes
  use io      ,only       : stdout
  implicit none
  !===============!
  != preparation =!
  !===============!
  
  call environment_start( 'LVCSH' )
  call get_inputfile()
  call readpwscf_out(scfoutname)
  call readph_out(phoutname)
  call readepwout(epwoutname)
  call treat_parameters()
  call set_H0_nk()
  call ph_configuration(nqtotf,nmodes,wf,temp)
  call init_random_seed()
  if(lsetthreads) call set_mkl_threads(mkl_threads)
  call allocatesh()
  
  
  !==========================!
  != loop over realizations =!
  !==========================!  
  do iaver=1,naver
    write(stdout,'(a,I4,a)') '###### iaver=',iaver,' ######'    
    !==================!
    != initialization =!
    !==================! 
    !!得到简正坐标的初始位置ph_l=phQ/sqrt(hbar/(2*wqv))和速度ph_p=phP /sqrt(hbar/(2*wqv))
    call init_normalmode_coordinate_velocity(nqtotf,nmodes,ph_l,ph_p,wf,temp)
    call init_dynamical_variable(nqtotf,nmodes,ph_l,c_nk,e,p,w)
    
    
    
!    call dia_syH(nbasis,H0,E,P)
!    call set_HH(nfreem,nbasis,phQ,H0,Hep,HH)
!    if(lelecsh) then
!      HH_e = HH
!      call dia_syH(nbasis,HH_e,E_e,P_e)
!    endif

!    !!得到初始的HH_e,c_e,w_e,n_e,isurface_e,HH_h,c_h,w_h,n_h,isurface_h
!    call get_initialsurface(initehstat,representation)

!    if(lelecsh) then
!      call calculate_nonadiabatic_coupling(nbasis,nfreem,E_e,P_e,Hep,d_e,dE_dQ_e)
!      E0_e = E_e; P0_e=P_e;d0_e=d_e;w0_e=w_e
!    endif
!    if(lholesh) then
!      call calculate_nonadiabatic_coupling(nbasis,nfreem,E_h,P_h,Hep,d_h,dE_dQ_h)
!      d_h = -d_h
!      dE_dQ_h = -dE_dQ_h
!      E0_h = E_h; P0_h=P_h;d0_h=d_h;w0_h=w_h
!    endif
!    phQ0 = phQ; phP0=phP
!
!    !=======================!
!    != loop over snapshots =!
!    !=======================!
!
!    do isnap=1,nsnap
!      do istep=1,nstep
!        
!        time1   = io_time()  
!        !==========================!
!        != update x,v,c,e,p,d,w,g =!
!        !==========================!
!        !rk4方法数组计算核的动力学，得到新的Q
!        call rk4_nuclei(nfreem,phQ,phP,dt)
!        !rk3方法计算电子空穴在透热表象下的演化
!        if(lelecsh) call rk4_electron_diabatic(nbasis,C_e,n_e,dt,HH_e)      
!        if(lholesh) call rk4_electron_diabatic(nbasis,C_h,n_h,dt,HH_h)
!        call set_HH(nfreem,nbasis,phQ,H0,Hep,HH)
!        if(lelecsh) HH_e = HH
!        if(lholesh) HH_h = -HH
!        if(trim(adjustl(shtype)) == "exciton") then
!          call set_H_with_Coulomb(nbasis,n_e,n_h,HH_e,HH_h)
!        endif
!        if(lelecsh) then
!          call dia_syH(nbasis,HH_e,E_e,P_e)
!          call calculate_nonadiabatic_coupling(nbasis,nfreem,E_e,P_e,Hep,d_e,dE_dQ_e)
!          call convert_diabatic_adiabatic(nbasis,P_e,C_e,W_e)
!          call calculate_hopping_probability(nbasis,nfreem,w0_e,phP0,d0_e,dt,g_e,g1_e,isurface_e)
!          call correct_hopping_probability(MSH,nbasis,g_e,g1_e,E0_e,w0_e,w_e,P0_e,P_e,isurface_e)
!          call nonadiabatic_transition(MSH,nbasis,nfreem,E0_e,P0_e,P_e,d_e,isurface_e,g_e,w0_e,phP,c_e)
!        endif
!        if(lholesh) then
!          call dia_syH(nbasis,HH_h,E_h,P_h)
!          call calculate_nonadiabatic_coupling(nbasis,nfreem,E_h,P_h,Hep,d_h,dE_dQ_h)
!          d_h     = -d_h
!          dE_dQ_h = -dE_dQ_h
!          call convert_diabatic_adiabatic(nbasis,P_h,C_h,W_h)
!          call calculate_hopping_probability(nbasis,nfreem,w0_h,phP0,d0_h,dt,g_h,g1_h,isurface_h)
!          call correct_hopping_probability(MSH,nbasis,g_h,g1_h,E0_h,w0_h,w_h,P0_h,P_h,isurface_h)
!          call nonadiabatic_transition(MSH,nbasis,nfreem,E0_h,P0_h,P_h,d_h,isurface_h,g_h,w0_h,phP,c_e)          
!        endif
!        
!        !===================!
!        != add bath effect =!
!        !===================!        
!        call add_bath_effect(nbasis,nfreem,d0_e,p0_e,d0_h,p0_h,womiga,dt,phQ,phP)
!
!        !============================!
!        != reset dynamical variable =!
!        !============================!
!
!        phQ0=phQ; phQ0=phQ; E0_e=E_e; p0_e=p_e; d0_e=d_e; w0_e=w_e;w0_h=w_h
!        E0_h=E_h; p0_h=p_h; d0_h=d_h
!        time2   = io_time()
!        write(stdout,"(A5,1X,I10,1X,F15.6,1X,A10)") 'Step=',istep,(time2-time1),"seconds"        
!      enddo
!      
!      !=====================!
!      != store information =!
!      !=====================!    
!
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
!      enddo
!           
!    enddo
!    !!
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
!  call cpu_time(t1)
!  write(6,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'
!  write(stdout,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'  
!  
!  close(stdout)
!  
!  stop
  
end program lvcsh


