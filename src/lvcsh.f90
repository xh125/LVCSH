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
  != Last updata 2021-05.18 Version:0.1.3                                                              =!   
  != Developed by XieHua at department of physic, USTC;xh125@mail.ustc.edu.cn                          =!
  !=====================================================================================================!
  !! Author: HuaXie
  !! Version: v0.1.1
  !! License: GNU
  !!=================================================================================================  
  use mkl_service
  use omp_lib
  use constants,only      : maxlen,ryd2eV,ry_to_fs
  use environments,only   : environment_start,mkl_threads,&
                            set_mkl_threads,lsetthreads
  use readinput,only      : get_inputfile
  use readscf,only        : readpwscf_out
  use readphout,only      : readph_out
  use readepw,only        : readepwout
  use parameters, only    : lreadscfout,scfoutname,lreadphout,phoutname,epwoutname,inputfilename,&
                            nsnap,nstep,dt,gamma,temp,init_ik,init_eband,init_hband,&
                            llaser,methodsh
  use hamiltonian,only    : nefre,neband,H_e,H_e_nk,E_e,P_e,P_e_nk,P0_e_nk,epcq_e,H0_e_nk,E0_e,P0_e,&
                            nhfre,nhband,H_h,H_h_nk,E_h,P_h,P_h_nk,P0_h_nk,epcq_h,H0_h_nk,E0_h,P0_h,&
                            allocate_hamiltonian,set_H_nk,set_H0_nk,&
                            calculate_eigen_energy_state
  use randoms,only        : init_random_seed
  use lasercom,only       : fwhm,w_laser
  use getwcvk,only        : get_Wcvk
  use initialsh,only      : set_subband,init_normalmode_coordinate_velocity,init_eh_KSstat,&
                            init_stat_diabatic,init_surface
  use surfacecom,only     : iaver,isnap,istep,naver,iesurface,ihsurface,iesurface_j,ihsurface_j,&
                            c_e,c_e_nk,d_e,d0_e,&
                            c_h,c_h_nk,d_h,d0_h,&
                            dEa_dQ,dEa_dQ_e,dEa_dQ_h,dEa2_dQ2,dEa2_dQ2_e,dEa2_dQ2_h
  use fssh,only           : nonadiabatic_transition_fssh
  use sc_fssh,only        : get_G_SC_FSSH,nonadiabatic_transition_scfssh
  use cc_fssh,only        : S_ai_e,S_ai_h,S_bi_e,S_bi_h,&
                            get_G_CC_FSSH,nonadiabatic_transition_ccfssh
  use surfacehopping,only : phQ,phP,phQ0,phP0,ph_T,ph_U,SUM_ph_E,SUM_ph_T,SUM_ph_U,&
                            w_e,w0_e,g_e,g1_e,esurface_type,cc0_e,dc1_e,dc2_e,dc3_e,dc4_e,&
                            w_h,w0_h,g_h,g1_h,hsurface_type,cc0_h,dc1_h,dc2_h,dc3_h,dc4_h,&
                            allocatesh,&
                            calculate_nonadiabatic_coupling,convert_diabatic_adiabatic,&
                            calculate_hopping_probability
  use surfacecom,only     : lelecsh,lholesh,ieband_min,ieband_max,ihband_min,ihband_max
  use elph2,only          : wf,nqtotf,nktotf,nbndfst
  use disp,only           : ph_configuration
  use modes,only          : nmodes
  use io      ,only       : stdout,io_time,time1,time2
  use dynamics,only       : get_dEa_dQ,get_dEa2_dQ2,rk4_nuclei,rk4_electron_diabatic,ADD_BATH_EFFECT
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
  call set_subband(lelecsh,lholesh,ieband_min,ieband_max,ihband_min,ihband_max)
  !get ieband_min,ieband_max,ihband_min,ihband_max
  call allocate_hamiltonian(lelecsh,lholesh,ieband_min,ieband_max,ihband_min,ihband_max)
  if(lelecsh) then
    call set_H0_nk(nktotf,neband,H0_e_nk,ieband_min,epcq_e)
  endif
  if(lholesh) then
    call set_H0_nk(nktotf,nhband,H0_h_nk,ihband_min,epcq_h)
    !H0_h_nk = -1.0 * H0_h_nk
    !epcq_h  = -1.0 * H0_h_nk
  endif
  !get H0_e_nk(neband,nktotf,neband,nktotf),epcq_e(neband,neband,nktotf,nmodes,nqtotf)
  !get H0_h_nk(nhband,nktotf,nhband,nktotf),epcq_h(nhband,nhband,nktotf,nmodes,nqtotf)
  
  call allocatesh(methodsh,lelecsh,lholesh,nmodes,nqtotf)
  
  if(llaser) then
    call get_Wcvk(ihband_min,ieband_max,fwhm,w_laser)
    write(stdout,"(/,5X,A)") "In the laser obsorbtion,the Pump laser as follow:"
    write(stdout,"(5X,A22,F12.7,A4)")  "Laser centred energy :",w_laser*ryd2eV," eV."
    write(stdout,"(5X,A38,F12.7,A4)")  "The full width at half-maximum:fwhm = ",fwhm*ry_to_fs," fs."
    !get W_cvk(icband,ivband,ik)
  endif
  
  call init_random_seed()
  if(lsetthreads) call set_mkl_threads(mkl_threads)
  
  
  !==========================!
  != loop over realizations =!
  !==========================!  
  do iaver=1,naver
    write(stdout,'(/,a,I4,a)') '###### iaver=',iaver,' ######'    
    
    !==================!
    != initialization =!
    !==================! 
    
    !!Get the initial normal mode coordinate phQ and versity phP
    call init_normalmode_coordinate_velocity(nmodes,nqtotf,wf,temp,phQ,phP)
    !应该先跑平衡后，再做电子空穴动力学计算
    
    !!得到初始电子和空穴的初始的KS状态 init_ik,init_eband,init_hband
    call init_eh_KSstat(lelecsh,lholesh,llaser,init_ik,init_eband,init_hband)
    
    
    if(lelecsh) then
      call init_stat_diabatic(init_ik,init_eband,ieband_min,neband,nktotf,c_e_nk)
      c_e = reshape(c_e_nk,(/nefre/))
      call set_H_nk(neband,nktotf,nmodes,nqtotf,phQ,wf,epcq_e,H0_e_nk,H_e_nk)
      H_e = reshape(H_e_nk,(/ nefre,nefre /))
      call calculate_eigen_energy_state(nefre,H_e,E_e,P_e)
      P_e_nk = reshape(P_e,(/ neband,nktotf,nefre /))
      call convert_diabatic_adiabatic(nefre,P_e,c_e,w_e)
      call init_surface(nefre,w_e,iesurface)
      call calculate_nonadiabatic_coupling(nmodes,nqtotf,neband,nktotf,wf,E_e,P_e_nk,epcq_e,d_e)
      E0_e = E_e;P0_e=P_e;P0_e_nk=P_e_nk;d0_e=d_e;w0_e=w_e
    endif
    
    if(lholesh) then
      call init_stat_diabatic(init_ik,init_hband,ihband_min,nhband,nktotf,c_h_nk)
      c_h = reshape(c_h_nk,(/nhfre/))
      call set_H_nk(nhband,nktotf,nmodes,nqtotf,phQ,wf,epcq_h,H0_h_nk,H_h_nk)
      H_h = reshape(H_h_nk,(/ nhfre,nhfre /))    
      call calculate_eigen_energy_state(nhfre,H_h,E_h,P_h)
      P_h_nk = reshape(P_h,(/ nhband,nktotf,nhfre /))
      call convert_diabatic_adiabatic(nhfre,P_h,c_h,w_h)
      call init_surface(nhfre,w_h,ihsurface)                
      call calculate_nonadiabatic_coupling(nmodes,nqtotf,nhband,nktotf,wf,E_h,P_h_nk,epcq_h,d_h)
      E0_h = E_h;P0_h=P_h;P0_h_nk=P_h_nk;d0_h=d_h;w0_h=w_h
    endif
    
    write(stdout,"(/,5X,A)") "In adiabatic base,the elctron-hole state as follow:"
    if(lholesh) then
      write(stdout,"(5X,A14,I5,1X,A20,F12.7,A3)") &
      "Init_hsurface=",ihsurface+ihband_min*nktotf,"Initial hole Energy:",E_h(ihsurface)*ryd2eV," eV"             
    endif
    if(lelecsh) then
      write(stdout,"(5X,A14,I5,1X,A20,F12.7,A3)") &
      "Init_esurface=",iesurface+ieband_min*nktotf,"Initial elec Energy:",E_e(iesurface)*ryd2eV," eV"       
    endif
    if(lelecsh .and. lholesh) then
      write(stdout,"(5X,A17,F12.7,A3)")  "elec-hole energy=",(E_e(iesurface)-E_h(ihsurface))*ryd2eV," eV"  
    endif
    
    
    
    phQ0=phQ; phP0=phP    


    if(lholesh) then
      if(lelecsh) then
        write(stdout,"(/,A)") "   time(fs) rt(s) hsur esur&
        &  E_h(eV)  E_e(eV) E_eh(eV) T_ph(eV) U_ph(eV) E_ph(eV) E_tot(eV)"         
      else 
        write(stdout,"(/,A)") "   time(fs) rt(s) hsur&
        &  E_h(eV) T_ph(eV) U_ph(eV) E_ph(eV) E_tot(eV)"       
      endif
    else
      if(lelecsh) then
        write(stdout,"(/,A)") "   time(fs) rt(s) esur&
        &  E_e(eV) T_ph(eV) U_ph(eV) E_ph(eV) E_tot(eV)"           
      else
        write(stdout,"(/,A)") "Error!! lelecsh and lholesh must have one need to be set TRUE."
      endif
    endif


    ph_U = 0.5*(wf**2)*(phQ**2)
    ph_T = 0.5*phP**2        
    SUM_ph_T = SUM(ph_T)
    SUM_ph_U = SUM(ph_U)
    SUM_ph_E = SUM_ph_T+SUM_ph_U    
    write(stdout,"(F11.2,F6.2,I5,I5,7(1X,F8.4))") 0.00,0.00,&
    ihsurface+ihband_min*nktotf,iesurface+ieband_min*nktotf,&
    e_h(ihsurface)*ryd2eV,e_e(iesurface)*ryd2eV,(e_e(iesurface)-e_h(ihsurface))*ryd2eV,&
    SUM_ph_T*ryd2eV,SUM_ph_U*ryd2eV,SUM_ph_E*ryd2eV,(E_e(iesurface)-E_h(ihsurface)+SUM_ph_E)*ryd2eV          
    !=======================!
    != loop over snapshots =!
    !=======================!

    do isnap=1,nsnap
      do istep=1,nstep
        
        time1   = io_time()  
        
        !==============================!
        != update phQ,phP             =!
        !==============================!
        
        dEa_dQ = 0.0
        !dEa_dQ in time t0
        if(lelecsh) then
          call get_dEa_dQ(nmodes,nqtotf,neband,nktotf,wf,P0_e_nk,epcq_e,iesurface,dEa_dQ_e)
          dEa_dQ = dEa_dQ + dEa_dQ_e
        endif
        if(lholesh) then
          call get_dEa_dQ(nmodes,nqtotf,nhband,nktotf,wf,P0_h_nk,epcq_h,ihsurface,dEa_dQ_h)
          dEa_dQ = dEa_dQ - dEa_dQ_h
        endif
        
        !use rk4 to calculate the dynamical of phonon normal modes
        !update phQ,phP to time t0+dt
        call rk4_nuclei(nmodes,nqtotf,dEa_dQ,gamma,wf,phQ,phP,dt)
        

        if(lelecsh) then
        
          !==============================!
          != update c,e,p,d,w,g         =!
          !==============================!
          
          if(methodsh == "SC-FSSH" .OR. methodsh == "CC-FSSH") then
            !electron wave function is propagated in diabatic representation
            !hamiltonian in time t0
            call set_H_nk(neband,nktotf,nmodes,nqtotf,phQ0,wf,epcq_e,H0_e_nk,H_e_nk)
            H_e = reshape(H_e_nk,(/ nefre,nefre /))
            !update c_e to time t0+dt
            call rk4_electron_diabatic(nefre,H_e,c_e,cc0_e,dc1_e,dc2_e,dc3_e,dc4_e,dt)
          endif
          
          ! update H_nk in time t0+dt
          call set_H_nk(neband,nktotf,nmodes,nqtotf,phQ,wf,epcq_e,H0_e_nk,H_e_nk)
          H_e = reshape(H_e_nk,(/ nefre,nefre /))          
          ! update E_e,P_e to time t0+dt
          call calculate_eigen_energy_state(nefre,H_e,E_e,P_e)
          P_e_nk = reshape(P_e,(/ neband,nktotf,nefre /))
          
          ! Calculate non-adiabatic coupling vectors with the Hellmann-Feynman theorem.
          ! update d_e in time t0+dt
          call calculate_nonadiabatic_coupling(nmodes,nqtotf,neband,nktotf,wf,E_e,p_e,epcq_e,d_e)
          
          if(methodsh == "SC-FSSH" .OR. methodsh == "CC-FSSH") then
            ! use p_e in time t0+dt, to convert c_e(t0+dt) to w_e(t0+dt) 
            call convert_diabatic_adiabatic(nefre,p_e,c_e,w_e)
          endif          
          
          ! use FSSH calculation hopping probability in adiabatic representation,get g_e,g1_e
          call calculate_hopping_probability(iesurface,nefre,nmodes,nqtotf,w0_e,phP0,d0_e,dt,g_e,g1_e)          
          
          
          !dealwith trilvial crossing,fixed ge
          if(methodsh == "SC-FSSH") then
            !use SC-FSSH method to fixed ge
            call get_G_SC_FSSH(iesurface,nefre,E0_e,w0_e,w_e,g1_e,g_e)
          elseif(methodsh == "CC-FSSH") then
            !use CC-FSSH method to fixed ge
            call get_G_CC_FSSH(nefre,iesurface,iesurface_j,p0_e,p_e,w0_e,w_e,S_ai_e,g1_e,g_e)
          endif
          
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
          !% change potential energy surface %!
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!          
          if(methodsh == "FSSH") then
            call nonadiabatic_transition_fssh(nefre,nqtotf,nmodes,iesurface,E0_e,P0_e,d0_e,g_e,phP)
          elseif( methodsh == "SC-FSSH") then
            call nonadiabatic_transition_scfssh(nefre,nqtotf,nmodes,iesurface,-E0_e,P0_e,d0_e,g_e,phP)              
          elseif(methodsh == "CC-FSSH") then
            call nonadiabatic_transition_ccfssh(nefre,nqtotf,nmodes,iesurface,iesurface_j,&
                                                esurface_type,E0_e,P0_e,P_e,d0_e,S_bi_e,g_e,phP)
          endif
          
        endif
        
        if(lholesh) then

          !==============================!
          != update c,e,p,d,w,g         =!
          !==============================!        
        
          if(methodsh == "SC-FSSH" .OR. methodsh == "CC-FSSH") then
            !hole wave function is propagated in diabatic representation
            !hamiltonian in time t0
            call set_H_nk(nhband,nktotf,nmodes,nqtotf,phQ0,wf,epcq_h,H0_h_nk,H_h_nk)
            H_h = reshape(H_h_nk,(/ nhfre,nhfre /))
            H_h = -1.0*H_h
            !update c_h to time t0+dt
            call rk4_electron_diabatic(nhfre,H_h,c_h,cc0_h,dc1_h,dc2_h,dc3_h,dc4_h,dt)        
          endif
          
          ! update H_nk in time t0+dt
          call set_H_nk(nhband,nktotf,nmodes,nqtotf,phQ,wf,epcq_h,H0_h_nk,H_h_nk)
          H_h = reshape(H_h_nk,(/ nhfre,nhfre /))          
          ! update E_h,P_h in time t0+dt
          call calculate_eigen_energy_state(nhfre,H_h,E_h,P_h)
          P_h_nk = reshape(P_h,(/ nhband,nktotf,nhfre /))
          
          ! Calculate non-adiabatic coupling vectors with the Hellmann-Feynman theorem.
          ! update d_h in time t0+dt
          call calculate_nonadiabatic_coupling(nmodes,nqtotf,neband,nktotf,wf,E_h,p_h,epcq_h,d_h)          


          if(methodsh == "SC-FSSH" .OR. methodsh == "CC-FSSH") then
            ! use p_h in time t0+dt, to convert c_h(t0+dt) to w_h(t0+dt) 
            if(lholesh) call convert_diabatic_adiabatic(nhfre,p_h,c_h,w_h)
          endif          

          ! use FSSH calculation hopping probability in adiabatic representation
          call calculate_hopping_probability(ihsurface,nhfre,nmodes,nqtotf,w0_h,phP0,d0_h,dt,g_h,g1_h)          
          
          !dealwith trilvial crossing,fixed ge
          if(methodsh == "SC-FSSH") then
            !use SC-FSSH method to fixed ge
            call get_G_SC_FSSH(ihsurface,nhfre,E0_h,w0_h,w_h,g1_h,g_h)
          elseif(methodsh == "CC-FSSH") then
            !use CC-FSSH method to fixed ge
            call get_G_CC_FSSH(nhfre,ihsurface,ihsurface_j,p0_h,p_h,w0_h,w_h,S_ai_h,g1_h,g_h)
          endif
          
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
          !% change potential energy surface %!
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!          
          if(methodsh == "FSSH") then
            call nonadiabatic_transition_fssh(nhfre,nqtotf,nmodes,ihsurface,-1.0*E0_h,P0_h,d0_h,g_h,phP)
          elseif( methodsh == "SC-FSSH") then
            call nonadiabatic_transition_scfssh(nhfre,nqtotf,nmodes,ihsurface,-1.0*E0_h,P0_h,d0_h,g_h,phP)              
          elseif(methodsh == "CC-FSSH") then
            call nonadiabatic_transition_ccfssh(nhfre,nqtotf,nmodes,ihsurface,ihsurface_j,hsurface_type,-1.0*E0_h,P0_h,P_h,d0_h,S_bi_h,g_h,phP)
          endif          
          
        endif 
   
        !===================!
        != add bath effect =!
        !===================! 
        dEa2_dQ2 = 0.0
        !dEa2_dQ2 in time t0
        if(lelecsh) then
          call get_dEa2_dQ2(nmodes,nqtotf,nefre,iesurface,E0_e,d0_e,dEa_dQ_e)
          dEa2_dQ2 = dEa2_dQ2 + dEa2_dQ2_e
        endif
        if(lholesh) then
          call get_dEa2_dQ2(nmodes,nqtotf,nhfre,ihsurface,E0_h,d0_h,dEa_dQ_h)
          dEa2_dQ2 = dEa2_dQ2 - dEa2_dQ2_h
        endif        
        call add_bath_effect(nmodes,nqtotf,gamma,temp,dEa2_dQ2,dt,phQ,phP)

        !============================!
        != reset dynamical variable =!
        !============================!

        phQ0=phQ; phP0=phP
        
        
        if(lelecsh) then
          E0_e = E_e;P0_e = P_e;P0_e_nk = P_e_nk; d0_e = d_e;w0_e = w_e
        endif
        if(lholesh) then
          E0_h = E_h;P0_h = P_h;P0_h_nk = P_h_nk; d0_h = d_h;w0_h = w_h
        endif
        time2   = io_time()

        
        ph_U = 0.5*(wf**2)*(phQ**2)
        ph_T = 0.5*phP**2        
        SUM_ph_T = SUM(ph_T)
        SUM_ph_U = SUM(ph_U)
        SUM_ph_E = SUM_ph_T+SUM_ph_U
        
        if(lholesh) then
          if(lelecsh) then
            !write(stdout,"(/,A)") "isnap istep runtime iesur ihsur  &
            !&en_e(eV)  en_h(eV)  en_eh(eV)  T_ph(eV)  U_ph(eV)  E_ph(eV)  E_tot(eV)" 
            write(stdout,"(F11.2,F6.2,I5,I5,7(1X,F8.4))") ((isnap-1)*nstep+istep)*dt*ry_to_fs,(time2-time1),&
            ihsurface+ihband_min*nktotf,iesurface+ieband_min*nktotf,&
            e_h(ihsurface)*ryd2eV,e_e(iesurface)*ryd2eV,(e_e(iesurface)-e_h(ihsurface))*ryd2eV,&
            SUM_ph_T*ryd2eV,SUM_ph_U*ryd2eV,SUM_ph_E*ryd2eV,(E_e(iesurface)-E_h(ihsurface)+SUM_ph_E)*ryd2eV            
          else 
            !write(stdout,"(/,A)") "isnap istep runtime ihsur  &
            !& en_h(eV)  T_ph(eV)  U_ph(eV)  E_ph(eV)  E_tot(eV)"  
            write(stdout,"(F11.2,F6.2,I5,I5,7(1X,F8.4))") ((isnap-1)*nstep+istep)*dt*ry_to_fs,(time2-time1),ihsurface,&
            e_h(ihsurface)*ryd2eV,&
            SUM_ph_T*ryd2eV,SUM_ph_U*ryd2eV,SUM_ph_E*ryd2eV,(-e_h(ihsurface)+SUM_ph_E)*ryd2eV            
          endif
        else
          if(lelecsh) then
            !write(stdout,"(/,A)") "isnap istep runtime iesur  &
            !&en_e(eV)  T_ph(eV)  U_ph(eV)  E_ph(eV)  E_tot(eV)"
            write(stdout,"(F11.2,F6.2,I5,I5,7(1X,F8.4))") ((isnap-1)*nstep+istep)*dt*ry_to_fs,(time2-time1),iesurface,&
            E_e(iesurface)*ryd2eV,(E_e(iesurface))*ryd2eV,&
            SUM_ph_T*ryd2eV,SUM_ph_U*ryd2eV,SUM_ph_E*ryd2eV,(E_e(iesurface)+SUM_ph_E)*ryd2eV            
          else
            write(stdout,"(/,A)") "Error!! lelecsh and lholesh must have one need to be set TRUE."
          endif
        endif        
        

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
  !====================!
  != save information =!
  !====================!
!  call saveresult()
  !call cpu_time(t1)
  !write(6,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'
  !write(stdout,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'  
  
  close(stdout)
  
  stop
  
end program lvcsh


