!----------------------------------------------------------------------------!
! Copyright (C) 2021 WeiZhuang's group                                       !
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
  != Last updata 2021-08.27 Version:0.3.0                                                              =!   
  != Developed by XieHua at department of physic, USTC;xh125@mail.ustc.edu.cn                          =!
  !=====================================================================================================!
  !! Author: HuaXie
  !! Version: v0.3.0
  !! License: GNU
  !!=================================================================================================  
  use mkl_service
  use omp_lib
  use f95_precision
  use blas95	
  use constants,only      : maxlen,ryd2eV,ry_to_fs,ryd2meV
  use environments,only   : environment_start,mkl_threads,&
                            set_mkl_threads,lsetthreads
  use readinput,only      : get_inputfile
  use readscf,only        : readpwscf_out
  use readphout,only      : readph_out
  use readepw,only        : readepwout
	use types,only          : gmnvkq_n0
  use parameters, only    : lreadscfout,scfoutname,lreadphout,phoutname,epwoutname,inputfilename,&
                            llaser,init_ik,init_eband,init_hband,mix_thr,lsortpes,&
														calculation,verbosity,nnode,ncore,naver_sum,savedsnap
  use hamiltonian,only    : nefre,neband,H_e,H_e_nk,E_e,P_e,P_e_nk,P0_e_nk,gmnvkq_e,Enk_e,H0_e_nk,E0_e,P0_e,&
                            nhfre,nhband,H_h,H_h_nk,E_h,P_h,P_h_nk,P0_h_nk,gmnvkq_h,Enk_h,H0_h_nk,E0_h,P0_h,&
                            allocate_hamiltonian,set_H_nk,set_H0_nk,&
                            calculate_eigen_energy_state,&
														ngfre_e,ngfre_h,gmnvkq_n0_e,gmnvkq_n0_h,get_gmnvkq_n0
	use sortting,only       : resort_eigen_energy_stat
  use randoms,only        : init_random_seed
  use lasercom,only       : fwhm,w_laser
  use getwcvk,only        : get_Wcvk
  use initialsh,only      : set_subband,init_normalmode_coordinate_velocity,init_eh_KSstat,&
                            init_stat_diabatic,init_surface
  use surfacecom,only     : methodsh,lfeedback,lit_gmnvkq,naver,nsnap,nstep,dt,pre_nstep,pre_dt,l_ph_quantum,&
                            gamma,temp,iaver,isnap,istep,ldecoherence,Cdecoherence,&
                            l_gamma_energy,gamma_min,gamma_max,ld_fric_min,ld_fric_max,n_gamma,&
                            lelecsh,lholesh,ieband_min,ieband_max,ihband_min,ihband_max,nefre_sh,nhfre_sh,&
                            iesurface,ihsurface,iesurface_j,ihsurface_j,&
                            c_e,c_e_nk,d_e,d0_e,&
                            c_h,c_h_nk,d_h,d0_h,&
                            dEa_dQ,dEa_dQ_e,dEa_dQ_h,dEa2_dQ2,dEa2_dQ2_e,dEa2_dQ2_h,&
                            csit_e,wsit_e,pes_one_e,pes_e,psit_e,&
                            csit_h,wsit_h,pes_one_h,pes_h,psit_h,&
                            E_ph_CA_sum,E_ph_QA_sum,ld_fric,ld_gamma
  use fssh,only           : nonadiabatic_transition_fssh
  use sc_fssh,only        : get_G_SC_FSSH,nonadiabatic_transition_scfssh
  use cc_fssh,only        : S_ai_e,S_ai_h,S_bi_e,S_bi_h,&
                            get_G_CC_FSSH,nonadiabatic_transition_ccfssh
  use surfacehopping,only : phQ,phP,phQ0,phP0,phK,phU,SUM_phE,SUM_phK,SUM_phK0,SUM_phU,&
                            phQsit,phPsit,phKsit,phUsit,&
                            w_e,w0_e,g_e,g1_e,esurface_type,cc0_e,dc1_e,dc2_e,dc3_e,dc4_e,&
                            w_h,w0_h,g_h,g1_h,hsurface_type,cc0_h,dc1_h,dc2_h,dc3_h,dc4_h,&
                            allocatesh,&
                            calculate_nonadiabatic_coupling,convert_diabatic_adiabatic,&
                            calculate_hopping_probability,convert_adiabatic_diabatic,&
														add_decoherence
  use elph2,only          : wf,xkf,nqtotf,nktotf,nbndfst,gmnvkq,etf
  use modes,only          : nmodes
	use cell_base,only      : bg
  use date_and_times,only : get_date_and_time
  use io      ,only       : stdout,io_time,time1,time2,time_,msg,io_error
  use dynamics,only       : set_gamma,get_dEa_dQ,get_dEa2_dQ2,rk4_nuclei,rk4_electron_diabatic,&
                            ADD_BATH_EFFECT
	use memory_report,only : MB,GB,complex_size, real_size,int_size,ram,print_memory  
  use saveinf,only : pes_e_file,pes_h_file,&
									   save_pes,save_phQ,save_phP,save_phK,save_phU,&
                     read_pes,read_phQ,read_phP,read_phK,read_phU,&
										 plot_pes,plot_phQ,plot_phP,plot_phK,plot_phU,&
										 csit_e_file,csit_h_file,save_csit,read_csit,plot_csit,&
                     wsit_e_file,wsit_h_file,save_wsit,read_wsit,plot_wsit,&
                     psit_e_file,psit_h_file,save_psit,read_psit,plot_psit,&
										 band_e_file,band_h_file,plot_band_occupatin_withtime
  implicit none
  
  !===============!
  != preparation =!
  !===============!
  integer :: iq,imode,ifre,igamma,ik,iband,inode,icore,iaver_i,iaver_f,ierr
  real(kind=8) :: t0,t1,time
	real(kind=dp) :: flagd,xkg_(3),xk_(3)
  character(len=9) :: cdate,ctime
  character(len=2) :: ctimeunit
  call cpu_time(t0)  

  
  call environment_start( 'LVCSH' )
  call get_inputfile(inputfilename)
  if(lreadscfout) call readpwscf_out(scfoutname)
  if(lreadphout)  call readph_out(phoutname)
  call readepwout(epwoutname)
  call set_subband(lelecsh,lholesh,ieband_min,ieband_max,ihband_min,ihband_max)
  !get ieband_min,ieband_max,ihband_min,ihband_max
  call allocate_hamiltonian(lelecsh,lholesh,ieband_min,ieband_max,ihband_min,ihband_max)
  
	if(lelecsh) then
    call set_H0_nk(nktotf,ieband_min,ieband_max,Enk_e,H0_e_nk,gmnvkq_e)
		!set H0_e_nk and gmnvkq_e
		if(nefre_sh == 0 .or. nefre_sh > nefre) nefre_sh = nefre
	endif
	
  if(lholesh) then
    call set_H0_nk(nktotf,ihband_min,ihband_max,Enk_h,H0_h_nk,gmnvkq_h)
    H0_h_nk = -1.0 * H0_h_nk
    gmnvkq_h  = -1.0 * gmnvkq_h		
		!set H0_h_nk and gmnvkq_h	
		if(nhfre_sh == 0 .or. nhfre_sh > nhfre) nhfre_sh = nhfre
  endif
	deallocate(gmnvkq)
  
  call allocatesh(methodsh,lelecsh,lholesh,nmodes,nqtotf)
  
	if(trim(calculation) == "lvcsh") then
  if(l_gamma_energy) then
	  !在BO势能面进行郎之万动力学计算，设置不同的gamma，测试摩擦系数对于动力学的影响。
    
    call init_normalmode_coordinate_velocity(nmodes,nqtotf,wf,temp,l_ph_quantum,phQ,phP)
    write(stdout,"(/5X,A51,F11.5,A2)") "The temperature of the non-adiabatic dynamica is : ",temp," K"
    write(stdout,"(5X,A49,F11.5,A4)") "The average energy of phonon(quantum): <SUM_phE>=",E_ph_QA_sum*ryd2eV," eV."
    write(stdout,"(5X,A49,F11.5,A4)") "The average energy of phonon(class)  : <SUM_phE>=",E_ph_CA_sum*ryd2eV," eV."
    if(l_ph_quantum) then
      write(stdout,"(/5X,A)") "The phonon dynamica set as quantum"
    else
      write(stdout,"(/5X,A)") "The phonon dynamica set as classical."
    endif
    
    write(stdout,"(5X,A)") repeat("=",67)
    write(stdout,"(5X,A)") "   gamma  ld_fric  phK(eV)  phU(eV)  phE(eV)"
    
    do igamma=0,n_gamma
      gamma   = gamma_min+igamma*(gamma_max-gamma_min)/n_gamma
      ld_fric = ld_fric_min+igamma*(ld_fric_max-ld_fric_min)/n_gamma
      call set_gamma(nmodes,nqtotf,gamma,ld_fric,wf,ld_gamma)
      dEa_dQ = 0.0
      dEa2_dQ2 = 0.0
      call init_normalmode_coordinate_velocity(nmodes,nqtotf,wf,temp,l_ph_quantum,phQ,phP)
      do istep=1,pre_nstep
        call rk4_nuclei(nmodes,nqtotf,dEa_dQ,ld_gamma,wf,phQ,phP,pre_dt)
        call add_bath_effect(nmodes,nqtotf,wf,ld_gamma,temp,dEa2_dQ2,pre_dt,l_ph_quantum,phQ,phP)
      enddo
      SUM_phK = 0.0
      SUM_phU = 0.0
      do istep=1,pre_nstep
        call rk4_nuclei(nmodes,nqtotf,dEa_dQ,ld_gamma,wf,phQ,phP,pre_dt)
        call add_bath_effect(nmodes,nqtotf,wf,ld_gamma,temp,dEa2_dQ2,pre_dt,l_ph_quantum,phQ,phP)
        SUM_phK = SUM_phK+ 0.5*SUM(phP**2)
        SUM_phU = SUM_phU+ 0.5*SUM((wf**2)*(phQ**2))
      enddo
      SUM_phK = SUM_phK/pre_nstep
      SUM_phU = SUM_phU/pre_nstep
      write(stdout,"(5X,5(F8.3,1X))") gamma,ld_fric,SUM_phK*ryd2eV,SUM_phU*ryd2eV,(SUM_phK+SUM_phU)*ryd2eV
    enddo
    
    stop
  endif
  
  ! set the friction coefficient of Langevin dynamica of all phonon modes.
  call set_gamma(nmodes,nqtotf,gamma,ld_fric,wf,ld_gamma)
  
  
  
  if(llaser) then
    call get_Wcvk(ihband_min,ieband_max,fwhm,w_laser)
    !get W_cvk(icband,ivband,ik)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !% Write laser information            %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    write(stdout,"(/,5X,A)") "In the laser obsorbtion,the Pump laser as follow:"
    write(stdout,"(5X,A22,F12.7,A4)")  "Laser centred energy :",w_laser*ryd2eV," eV."
    write(stdout,"(5X,A38,F12.7,A4)")  "The full width at half-maximum:fwhm = ",fwhm*ry_to_fs," fs."
    
    
  endif
  
  call init_random_seed()
  if(lsetthreads) call set_mkl_threads(mkl_threads)
  
  
  !==========================!
  != loop over realizations =!
  !==========================!  
  do iaver=1,naver
    write(stdout,'(/,a,I4,a)') '###### iaver=',iaver,' ######'    
		call get_date_and_time(cdate,ctime)
		write(stdout,'(1X,"This trajectory start on ",A9," at ",A9)') cdate,ctime      
    !==================!
    != initialization =!
    !==================!
    
    !!Get the initial normal mode coordinate phQ and versity phP
    call init_normalmode_coordinate_velocity(nmodes,nqtotf,wf,temp,l_ph_quantum,phQ,phP)
    !应该先跑平衡后，再做电子空穴动力学计算   
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !% Write phonon energy information         %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    write(stdout,"(/5X,A51,F11.5,A2)") "The temperature of the non-adiabatic dynamica is : ",temp," K"
    write(stdout,"(5X,A49,F11.5,A4)") "The average energy of phonon(quantum): <SUM_phE>=",E_ph_QA_sum*ryd2eV," eV."
    write(stdout,"(5X,A49,F11.5,A4)") "The average energy of phonon(class)  : <SUM_phE>=",E_ph_CA_sum*ryd2eV," eV."
    
    if(l_ph_quantum) then
      write(stdout,"(/5X,A)") "The phonon dynamica set as quantum"
    else
      write(stdout,"(/5X,A)") "The phonon dynamica set as classical."
    endif
    write(stdout,"(5X,A38,F11.5,A4,A9,F11.5,A4)") &
    "The initial energy of phonon: SUM_phT=",0.5*SUM(phP**2)*ryd2eV," eV",&
    " SUM_phU=",0.5*SUM(wf**2*phQ**2)*ryd2eV," eV"
    write(stdout,"(5X,A38,F11.5,A4)") &
    "The initial energy of phonon: SUM_phE=",0.5*SUM(phP**2+wf**2*phQ**2)*ryd2eV," eV."
    
    dEa_dQ = 0.0
    dEa2_dQ2 = 0.0
    do istep=1,pre_nstep
      call rk4_nuclei(nmodes,nqtotf,dEa_dQ,ld_gamma,wf,phQ,phP,pre_dt)
      call add_bath_effect(nmodes,nqtotf,wf,ld_gamma,temp,dEa2_dQ2,pre_dt,l_ph_quantum,phQ,phP)
    enddo

    time = pre_nstep*pre_dt*ry_to_fs
    if(time<1.0E3) then
      ctimeunit = 'fs'
    elseif(time<1.0E6) then
      time = time/1.0E3
      ctimeunit = 'ps'
    elseif(time<1.0E9) then
      time = time/1.0E6
      ctimeunit = 'ns'
    endif
    
    write(stdout,"(5X,A23,F6.2,A2,A19,F11.5,A4,A9,F11.5,A4)") &
    "Energy of phonon after ", time,ctimeunit," dynamica: SUM_phT=",0.5*SUM(phP**2)*ryd2eV," eV",&
    " SUM_phU=",0.5*SUM(wf**2*phQ**2)*ryd2eV," eV"
    write(stdout,"(5X,A23,F6.2,A2,A19,F11.5,A4)") &
    "Energy of phonon after ", time,ctimeunit," dynamica: SUM_phE=",0.5*SUM(phP**2+wf**2*phQ**2)*ryd2eV," eV."    
    
    !!得到初始电子和空穴的初始的KS状态 init_ik,init_eband,init_hband(in the diabatic states)
    call init_eh_KSstat(lelecsh,lholesh,llaser,init_ik,init_eband,init_hband)
    
    if(lelecsh) then
      call init_stat_diabatic(init_ik,init_eband,ieband_min,neband,nktotf,c_e_nk)
      c_e = reshape(c_e_nk,(/nefre/))
      call set_H_nk(neband,nktotf,nmodes,nqtotf,phQ,gmnvkq_e,H0_e_nk,H_e_nk)
      H_e = reshape(H_e_nk,(/ nefre,nefre /))
      call calculate_eigen_energy_state(nefre,H_e,E_e,P_e)
      P_e_nk = reshape(P_e,(/ neband,nktotf,nefre /))
      call convert_diabatic_adiabatic(nefre,P_e,c_e,w_e)
      call init_surface(nefre,nefre_sh,w_e,iesurface)

      call calculate_nonadiabatic_coupling(nmodes,nqtotf,neband,nktotf,E_e,P_e_nk,gmnvkq_e,lit_gmnvkq,nefre_sh,d_e)
      E0_e = E_e;P0_e=P_e;P0_e_nk=P_e_nk;d0_e=d_e;w0_e=w_e
    endif
    
    if(lholesh) then
      call init_stat_diabatic(init_ik,init_hband,ihband_min,nhband,nktotf,c_h_nk)
      c_h = reshape(c_h_nk,(/nhfre/))
      call set_H_nk(nhband,nktotf,nmodes,nqtotf,phQ,gmnvkq_h,H0_h_nk,H_h_nk)
      H_h = reshape(H_h_nk,(/ nhfre,nhfre /))    
      call calculate_eigen_energy_state(nhfre,H_h,E_h,P_h)
      P_h_nk = reshape(P_h,(/ nhband,nktotf,nhfre /))
      call convert_diabatic_adiabatic(nhfre,P_h,c_h,w_h)
      call init_surface(nhfre,nhfre_sh,w_h,ihsurface)
		
      call calculate_nonadiabatic_coupling(nmodes,nqtotf,nhband,nktotf,E_h,P_h_nk,gmnvkq_h,lit_gmnvkq,nhfre_sh,d_h)
      E0_h = E_h;P0_h=P_h;P0_h_nk=P_h_nk;d0_h=d_h;w0_h=w_h
    endif

    
    phQ0=phQ; phP0=phP 
    
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !% Write the elctron-hole information in adiabatic base   %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    write(stdout,"(/,5X,A)") "In adiabatic base,the elctron-hole state as follow:"
    if(lholesh) then
      write(stdout,"(5X,A14,I5,1X,A20,F12.7,A3)") &
      "Init_hsurface=",ihsurface,"Initial hole Energy:",-E_h(ihsurface)*ryd2eV," eV" 
    endif
    if(lelecsh) then
      write(stdout,"(5X,A14,I5,1X,A20,F12.7,A3)") &
      "Init_esurface=",iesurface,"Initial elec Energy:",E_e(iesurface)*ryd2eV," eV"       
    endif
    if(lelecsh .and. lholesh) then
      write(stdout,"(5X,A17,F12.7,A3)")  "elec-hole energy=",(E_e(iesurface)+E_h(ihsurface))*ryd2eV," eV"  
    endif
    
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !% calculate phonon energy                                %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    phU = 0.5*(wf**2)*(phQ**2)
    phK = 0.5*phP**2        
    SUM_phK = SUM(phK)
    SUM_phU = SUM(phU)
    SUM_phE = SUM_phK+SUM_phU
    
		SUM_phK0 = SUM_phK
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !% Write initial non-adiabatic dynamica information        %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    if(lholesh) then
      if(lelecsh) then
        write(stdout,"(/,A)") "   time(fs) rt(s) hsur esur&
        &  E_h(eV)  E_e(eV) E_eh(eV) T_ph(eV) U_ph(eV) E_ph(eV) E_tot(eV)"         
        write(stdout,"(F11.2,F6.2,I5,I5,7(1X,F8.4))") 0.00,0.00,&
        ihsurface,iesurface,&
        -e_h(ihsurface)*ryd2eV,e_e(iesurface)*ryd2eV,(e_e(iesurface)+e_h(ihsurface))*ryd2eV,&
        SUM_phK*ryd2eV,SUM_phU*ryd2eV,SUM_phE*ryd2eV,(E_e(iesurface)+E_h(ihsurface)+SUM_phE)*ryd2eV       
      else 
        write(stdout,"(/,A)") "   time(fs) rt(s) hsur&
        &  E_h(eV) T_ph(eV) U_ph(eV) E_ph(eV) E_tot(eV)"       
        write(stdout,"(F11.2,F6.2,I5,5(1X,F8.4))") 0.00,0.00,&
        ihsurface,-e_h(ihsurface)*ryd2eV,&
        SUM_phK*ryd2eV,SUM_phU*ryd2eV,SUM_phE*ryd2eV,(E_h(ihsurface)+SUM_phE)*ryd2eV     
      endif
    else
      if(lelecsh) then
        write(stdout,"(/,A)") "   time(fs) rt(s) esur&
        &  E_e(eV) T_ph(eV) U_ph(eV) E_ph(eV) E_tot(eV)" 
        write(stdout,"(F11.2,F6.2,I5,5(1X,F8.4))") 0.00,0.00,&
        iesurface,e_e(iesurface)*ryd2eV,&
        SUM_phK*ryd2eV,SUM_phU*ryd2eV,SUM_phE*ryd2eV,(E_e(iesurface)+SUM_phE)*ryd2eV     
      else
        write(stdout,"(/,A)") "Error!! lelecsh and lholesh must have one need to be set TRUE."
      endif
    endif



    !=============================!
    != store initial information =!
    !=============================!    
    
			isnap = 0
			do iq=1,nqtotf
        do imode=1,nmodes
          phQsit(imode,iq,isnap) = phQsit(imode,iq,isnap)+phQ(imode,iq)
          phPsit(imode,iq,isnap) = phPsit(imode,iq,isnap)+phP(imode,iq)
          phKsit(imode,iq,isnap) = phKsit(imode,iq,isnap)+phK(imode,iq)
          phUsit(imode,iq,isnap) = phUsit(imode,iq,isnap)+phU(imode,iq)
        enddo
      enddo
      
      if(lelecsh) then
        pes_e(0,isnap) = pes_e(0,isnap)+E_e(iesurface)
				if(iaver == 1) pes_one_e(0,isnap) = E_e(iesurface)
				do ifre=1,nefre	
					csit_e(ifre,isnap) = csit_e(ifre,isnap)+REAL(c_e(ifre)*CONJG(c_e(ifre)))
					wsit_e(ifre,isnap) = wsit_e(ifre,isnap)+REAL(w_e(ifre)*CONJG(w_e(ifre)))
					psit_e(ifre,isnap) = psit_e(ifre,isnap)+P_e(ifre,iesurface)**2
					pes_e( ifre,isnap) = pes_e( ifre,isnap) + E_e(ifre)
					if(iaver == 1) pes_one_e(ifre,isnap) = E_e(ifre)
				enddo
      endif
      
      if(lholesh) then
        pes_h(0,isnap) = pes_h(0,isnap)-E_h(ihsurface)
				if(iaver == 1) pes_one_h(0,isnap) = -E_h(ihsurface)
				do ifre=1,nhfre
					csit_h(ifre,isnap) = csit_h(ifre,isnap)+REAL(c_h(ifre)*CONJG(c_h(ifre)))
					wsit_h(ifre,isnap) = wsit_h(ifre,isnap)+REAL(w_h(ifre)*CONJG(w_h(ifre)))
					psit_h(ifre,isnap) = psit_h(ifre,isnap)+P_h(ifre,iesurface)**2
					pes_h( ifre,isnap) = pes_h( ifre,isnap)-E_h(ifre)
					if(iaver == 1) pes_one_h(ifre,isnap) = -E_h(ifre)
        enddo		
      endif
    
		!=================================!
    != End store initial information =!
    !=================================!    
			
		 
    !=======================!
    != loop over snapshots =!
    !=======================!

    do isnap=1,nsnap
      do istep=1,nstep
        
        time1   = io_time()  
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !% calculate dEa_dQ                     %!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        dEa_dQ = 0.0
        !dEa_dQ in time t0
        if(lelecsh) then
          call get_dEa_dQ(nmodes,nqtotf,neband,nktotf,P0_e_nk,gmnvkq_e,iesurface,dEa_dQ_e)
          dEa_dQ = dEa_dQ + dEa_dQ_e
        endif
        if(lholesh) then
          call get_dEa_dQ(nmodes,nqtotf,nhband,nktotf,P0_h_nk,gmnvkq_h,ihsurface,dEa_dQ_h)
          dEa_dQ = dEa_dQ + dEa_dQ_h
        endif
        
        
        !==============================!
        != update phQ,phP             =!
        !==============================!        
        !use rk4 to calculate the dynamical of phonon normal modes
        !update phQ,phP to time t0+dt
        call rk4_nuclei(nmodes,nqtotf,dEa_dQ,ld_gamma,wf,phQ,phP,dt)
        
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !% update c,e,p,d,w,g and change potential energy surface        %!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        if(lelecsh) then
      
          !electron wave function is propagated in diabatic representation
          !hamiltonian in time t0
          call set_H_nk(neband,nktotf,nmodes,nqtotf,phQ0,gmnvkq_e,H0_e_nk,H_e_nk)
          H_e = reshape(H_e_nk,(/ nefre,nefre /))
          !update c_e to time t0+dt
          call rk4_electron_diabatic(nefre,H_e,c_e,cc0_e,dc1_e,dc2_e,dc3_e,dc4_e,dt)
         
          
          ! update H_nk in time t0+dt
          call set_H_nk(neband,nktotf,nmodes,nqtotf,phQ,gmnvkq_e,H0_e_nk,H_e_nk)
          H_e = reshape(H_e_nk,(/ nefre,nefre /))          
          ! update E_e,P_e to time t0+dt
          call calculate_eigen_energy_state(nefre,H_e,E_e,P_e)
					!call resort_eigen_energy_stat(nefre,E_e,P_e,E_e_eq,P_e_eq)
					if(lsortpes) call resort_eigen_energy_stat(nefre,E_e,P_e,E0_e,P0_e,mix_thr)
          P_e_nk = reshape(P_e,(/ neband,nktotf,nefre /))
          
          ! Calculate non-adiabatic coupling vectors with the Hellmann-Feynman theorem.
          ! update d_e in time t0+dt
          call calculate_nonadiabatic_coupling(nmodes,nqtotf,neband,nktotf,E_e,p_e,gmnvkq_e,lit_gmnvkq,nefre_sh,d_e)
          
   
          ! use p_e in time t0+dt, to convert c_e(t0+dt) to w_e(t0+dt) 
          call convert_diabatic_adiabatic(nefre,p_e,c_e,w_e)
					!add decoherence induced change in adiabatic representation
					if(ldecoherence) then
						call add_decoherence(Cdecoherence,SUM_phK0,dt,nefre,iesurface,E0_e,w_e)
						! convert w_e(t0+dt) to c_e(t0+dt)
						call convert_adiabatic_diabatic(nefre,p_e,w_e,c_e)
					endif
						
            
          
          ! use FSSH calculation hopping probability in adiabatic representation,get g_e,g1_e
          call calculate_hopping_probability(iesurface,nefre,nefre_sh,nmodes,nqtotf,w0_e,phP0,d0_e,dt,g_e,g1_e)          
          
          !dealwith trilvial crossing,by fixed ge
          if(methodsh == "SC-FSSH") then
            !use SC-FSSH method to fixed ge
            call get_G_SC_FSSH(iesurface,nefre,nefre_sh,E0_e,w0_e,w_e,g1_e,g_e)
          elseif(methodsh == "CC-FSSH") then
            !use CC-FSSH method to fixed ge
            call get_G_CC_FSSH(nefre,nefre_sh,iesurface,iesurface_j,p0_e,p_e,w0_e,w_e,S_ai_e,g1_e,g_e)
          endif
          
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
          !% change potential energy surface %!
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!          
          if(methodsh == "FSSH") then
            call nonadiabatic_transition_fssh(lfeedback,nefre,nefre_sh,nqtotf,nmodes,iesurface,E0_e,P0_e,d0_e,g_e,phP)
          elseif( methodsh == "SC-FSSH") then
            call nonadiabatic_transition_scfssh(lfeedback,nefre,nefre_sh,nqtotf,nmodes,iesurface,E0_e,P0_e,d0_e,g_e,phP)              
          elseif(methodsh == "CC-FSSH") then
            call nonadiabatic_transition_ccfssh(lfeedback,nefre,nefre_sh,nqtotf,nmodes,iesurface,iesurface_j,&
                                                &esurface_type,E0_e,P0_e,P_e,d0_e,S_bi_e,g_e,phP)
          endif
          
        endif
        
        if(lholesh) then    
          
          !hole wave function is propagated in diabatic representation
          !hamiltonian in time t0
          call set_H_nk(nhband,nktotf,nmodes,nqtotf,phQ0,gmnvkq_h,H0_h_nk,H_h_nk)
          H_h = reshape(H_h_nk,(/ nhfre,nhfre /))
          !update c_h to time t0+dt
          call rk4_electron_diabatic(nhfre,H_h,c_h,cc0_h,dc1_h,dc2_h,dc3_h,dc4_h,dt)        
          
          
          ! update H_nk in time t0+dt
          call set_H_nk(nhband,nktotf,nmodes,nqtotf,phQ,gmnvkq_h,H0_h_nk,H_h_nk)
          H_h = reshape(H_h_nk,(/ nhfre,nhfre /))          
          ! update E_h,P_h in time t0+dt
          call calculate_eigen_energy_state(nhfre,H_h,E_h,P_h)
					!call resort_eigen_energy_stat(nhfre,E_h,P_h,E_h_eq,P_h_eq)
					if(lsortpes) call resort_eigen_energy_stat(nhfre,E_h,P_h,E0_h,P0_h,mix_thr)
          P_h_nk = reshape(P_h,(/ nhband,nktotf,nhfre /))
          
          ! Calculate non-adiabatic coupling vectors with the Hellmann-Feynman theorem.
          ! update d_h in time t0+dt
          call calculate_nonadiabatic_coupling(nmodes,nqtotf,nhband,nktotf,E_h,p_h,gmnvkq_h,lit_gmnvkq,nhfre_sh,d_h)          

          
          ! use p_h in time t0+dt, to convert c_h(t0+dt) to w_h(t0+dt) 
          call convert_diabatic_adiabatic(nhfre,p_h,c_h,w_h)
					if(ldecoherence) then
						call add_decoherence(Cdecoherence,SUM_phK0,dt,nhfre,ihsurface,E0_h,w_h)
						call convert_adiabatic_diabatic(nhfre,p_h,w_h,c_h)
					endif						
						
                   

          ! use FSSH calculation hopping probability in adiabatic representation
          call calculate_hopping_probability(ihsurface,nhfre,nhfre_sh,nmodes,nqtotf,w0_h,phP0,d0_h,dt,g_h,g1_h)          
          
          !dealwith trilvial crossing,fixed ge
          if(methodsh == "SC-FSSH") then
            !use SC-FSSH method to fixed ge
            call get_G_SC_FSSH(ihsurface,nhfre,nhfre_sh,E0_h,w0_h,w_h,g1_h,g_h)
          elseif(methodsh == "CC-FSSH") then
            !use CC-FSSH method to fixed ge
            call get_G_CC_FSSH(nhfre,nhfre_sh,ihsurface,ihsurface_j,p0_h,p_h,w0_h,w_h,S_ai_h,g1_h,g_h)
          endif
          
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
          !% change potential energy surface %!
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!          
          if(methodsh == "FSSH") then
            call nonadiabatic_transition_fssh(lfeedback,nhfre,nhfre_sh,nqtotf,nmodes,ihsurface,E0_h,P0_h,d0_h,g_h,phP)
          elseif( methodsh == "SC-FSSH") then
            call nonadiabatic_transition_scfssh(lfeedback,nhfre,nhfre_sh,nqtotf,nmodes,ihsurface,E0_h,P0_h,d0_h,g_h,phP)              
          elseif(methodsh == "CC-FSSH") then
            call nonadiabatic_transition_ccfssh(lfeedback,nhfre,nhfre_sh,nqtotf,nmodes,ihsurface,ihsurface_j,&
                                          &hsurface_type,E0_h,P0_h,P_h,d0_h,S_bi_h,g_h,phP)
          endif          
          
        endif 
   
   
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !% calculate dEa2_dQ2                   %!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        dEa2_dQ2 = 0.0
        !dEa2_dQ2 in time t0
        if(lelecsh) then
          call get_dEa2_dQ2(nmodes,nqtotf,nefre,nefre_sh,iesurface,E0_e,d0_e,dEa_dQ_e)
          dEa2_dQ2 = dEa2_dQ2 + dEa2_dQ2_e
        endif
        if(lholesh) then
          call get_dEa2_dQ2(nmodes,nqtotf,nhfre,nhfre_sh,ihsurface,E0_h,d0_h,dEa_dQ_h)
          dEa2_dQ2 = dEa2_dQ2 + dEa2_dQ2_h
        endif


        !===================!
        != add bath effect =!
        !===================!         
        call add_bath_effect(nmodes,nqtotf,wf,ld_gamma,temp,dEa2_dQ2,dt,l_ph_quantum,phQ,phP)


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
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        !% Write non-adiabatic dynamica energy information every step       %!
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
        phU = 0.5*(wf**2)*(phQ**2)
        phK = 0.5*phP**2        
        SUM_phK = SUM(phK)
        SUM_phU = SUM(phU)
        SUM_phE = SUM_phK+SUM_phU 
				SUM_phK0= SUM_phK
        time_ = ((isnap-1)*nstep+istep)*dt*ry_to_fs
        if(lholesh) then
          if(lelecsh) then
            !write(stdout,"(/,A)") "isnap istep runtime iesur ihsur  &
            !&en_e(eV)  en_h(eV)  en_eh(eV)  T_ph(eV)  U_ph(eV)  E_ph(eV)  E_tot(eV)" 
           if(trim(verbosity)=="high" .or. isnap == 1 .or. isnap == nsnap) then
						write(stdout,"(F11.2,F6.2,I5,I5,7(1X,F8.4))") time_,(time2-time1),ihsurface,iesurface,&
            -e_h(ihsurface)*ryd2eV,e_e(iesurface)*ryd2eV,(e_e(iesurface)+e_h(ihsurface))*ryd2eV,&
            SUM_phK*ryd2eV,SUM_phU*ryd2eV,SUM_phE*ryd2eV,(E_e(iesurface)+E_h(ihsurface)+SUM_phE)*ryd2eV
					 endif
					else 
            !write(stdout,"(/,A)") "isnap istep runtime ihsur  &
            !& en_h(eV)  T_ph(eV)  U_ph(eV)  E_ph(eV)  E_tot(eV)"  
           if(trim(verbosity)=="high" .or. isnap == 1 .or. isnap == nsnap) then
						write(stdout,"(F11.2,F6.2,I5,5(1X,F8.4))") time_,(time2-time1),ihsurface,-e_h(ihsurface)*ryd2eV,&
            SUM_phK*ryd2eV,SUM_phU*ryd2eV,SUM_phE*ryd2eV,(E_h(ihsurface)+SUM_phE)*ryd2eV
					 endif
					endif
        else
          if(lelecsh) then
            !write(stdout,"(/,A)") "isnap istep runtime iesur  &
            !&en_e(eV)  T_ph(eV)  U_ph(eV)  E_ph(eV)  E_tot(eV)"
           if(trim(verbosity)=="high" .or. isnap == 1 .or. isnap == nsnap) then
						write(stdout,"(F11.2,F6.2,I5,5(1X,F8.4))") time_,(time2-time1),iesurface,E_e(iesurface)*ryd2eV,&
            SUM_phK*ryd2eV,SUM_phU*ryd2eV,SUM_phE*ryd2eV,(E_e(iesurface)+SUM_phE)*ryd2eV            
           endif
					else
            write(stdout,"(/,A)") "Error!! lelecsh and lholesh must have one need to be set TRUE."
          endif
        endif        
      
      
      enddo
      
      !=====================!
      != store information =!
      !=====================!    
			do iq=1,nqtotf
        do imode=1,nmodes
          phQsit(imode,iq,isnap) = phQsit(imode,iq,isnap)+phQ(imode,iq)
          phPsit(imode,iq,isnap) = phPsit(imode,iq,isnap)+phP(imode,iq)
          phKsit(imode,iq,isnap) = phKsit(imode,iq,isnap)+phK(imode,iq)
          phUsit(imode,iq,isnap) = phUsit(imode,iq,isnap)+phU(imode,iq)
        enddo
      enddo
      
      if(lelecsh) then
        pes_e(0,isnap) = pes_e(0,isnap)+E_e(iesurface)
				if(iaver == 1) pes_one_e(0,isnap) = E_e(iesurface)
				do ifre=1,nefre	
					csit_e(ifre,isnap) = csit_e(ifre,isnap)+REAL(c_e(ifre)*CONJG(c_e(ifre)))
					wsit_e(ifre,isnap) = wsit_e(ifre,isnap)+REAL(w_e(ifre)*CONJG(w_e(ifre)))
					psit_e(ifre,isnap) = psit_e(ifre,isnap)+P_e(ifre,iesurface)**2
					pes_e( ifre,isnap) = pes_e( ifre,isnap) + E_e(ifre)
					if(iaver == 1) pes_one_e(ifre,isnap) = E_e(ifre)
				enddo
      endif
      
      if(lholesh) then
        pes_h(0,isnap) = pes_h(0,isnap)-E_h(ihsurface)
				if(iaver == 1) pes_one_h(0,isnap) = -E_h(ihsurface)
				do ifre=1,nhfre
					csit_h(ifre,isnap) = csit_h(ifre,isnap)+REAL(c_h(ifre)*CONJG(c_h(ifre)))
					wsit_h(ifre,isnap) = wsit_h(ifre,isnap)+REAL(w_h(ifre)*CONJG(w_h(ifre)))
					psit_h(ifre,isnap) = psit_h(ifre,isnap)+P_h(ifre,iesurface)**2
					pes_h( ifre,isnap) = pes_h( ifre,isnap)-E_h(ifre)
					if(iaver == 1) pes_one_h(ifre,isnap) = -E_h(ifre)
        enddo		
      endif
    
			!=========================!
			!= End store information =!
			!=========================!
			
    enddo
		
		call get_date_and_time(cdate,ctime)
		write(stdout,'(1X,"This trajectory end on ",A9," at ",A9)') cdate,ctime      
		
  enddo
  
  phQsit = phQsit / naver
  phPsit = phPsit / naver
  phKsit = phKsit / naver
  phUsit = phUsit / naver
  if(lelecsh) then
    csit_e = csit_e /naver
    wsit_e = wsit_e /naver
    psit_e = psit_e /naver
		pes_e  = pes_e  /naver
		!pes_one_e
  endif
  
  if(lholesh) then
    csit_h = csit_h /naver
    wsit_h = wsit_h /naver
    psit_h = psit_h /naver
		pes_h  = pes_h  /naver
		!pes_one_h
  endif


  !====================!
  != save information =!
  !====================!
  call save_phQ(nmodes,nqtotf,nsnap,phQsit)
  call save_phP(nmodes,nqtotf,nsnap,phPsit)
  call save_phK(nmodes,nqtotf,nsnap,phKsit)
  call save_phU(nmodes,nqtotf,nsnap,phUsit)
  
  if(lelecsh) then
    call save_pes(nefre,nsnap,naver,pes_one_e,pes_e,pes_e_file)
    call save_csit(nefre,nsnap,naver,csit_e,csit_e_file)
    call save_wsit(nefre,nsnap,naver,wsit_e,wsit_e_file)
    call save_psit(nefre,nsnap,naver,psit_e,psit_e_file)
		call plot_band_occupatin_withtime(neband,nktotf,Enk_e,xkf,nsnap,psit_e,csit_e,savedsnap,band_e_file)
  endif
  
  if(lholesh) then
    call save_pes(nhfre,nsnap,naver,pes_one_h,pes_h,pes_h_file)
    call save_csit(nhfre,nsnap,naver,csit_h,csit_h_file)
    call save_wsit(nhfre,nsnap,naver,wsit_h,wsit_h_file)
    call save_psit(nhfre,nsnap,naver,psit_h,psit_h_file)
		call plot_band_occupatin_withtime(nhband,nktotf,Enk_h,xkf,nsnap,psit_h,csit_h,savedsnap,band_h_file)
  endif
  
  
	elseif(trim(calculation)=="plot") then
		write(stdout,"(/A)") "Start to write files for plotting!"
		write(stdout,"(A,I4)") "Number of nodes for non-adiabatic calculation:",nnode
		write(stdout,"(A,I4)") "Number of samples for each node:",ncore		
		naver_sum = naver*nnode*ncore
		phQsit = 0.0
		phPsit = 0.0
		phKsit = 0.0
		phUsit = 0.0
		
		do inode=1,nnode
			write(stdout,"(A,I4)") "Read information in node:",inode
			do icore=1,ncore
				call read_phQ(inode,icore,nmodes,nqtotf,nsnap,phQsit)
				call read_phP(inode,icore,nmodes,nqtotf,nsnap,phPsit)
				call read_phK(inode,icore,nmodes,nqtotf,nsnap,phKsit)
				call read_phU(inode,icore,nmodes,nqtotf,nsnap,phUsit)
				
				if(lelecsh) then
					call read_pes( inode,icore,nefre,nsnap,naver,pes_one_e,pes_e,pes_e_file)
					call read_csit(inode,icore,nefre,nsnap,naver,csit_e,csit_e_file)
					call read_wsit(inode,icore,nefre,nsnap,naver,wsit_e,wsit_e_file)
					call read_psit(inode,icore,nefre,nsnap,naver,psit_e,psit_e_file)
				endif
				
				if(lholesh) then
					call read_pes( inode,icore,nhfre,nsnap,naver,pes_one_h,pes_h,pes_h_file)
					call read_csit(inode,icore,nhfre,nsnap,naver,csit_h,csit_h_file)
					call read_wsit(inode,icore,nhfre,nsnap,naver,wsit_h,wsit_h_file)
					call read_psit(inode,icore,nhfre,nsnap,naver,psit_h,psit_h_file)
				endif

			enddo
	  enddo
		write(stdout,"(A)") "Read all resut in different nodes Success."

		
		phQsit = phQsit /(nnode*ncore)
		phPsit = phPsit /(nnode*ncore)
		phKsit = phKsit /(nnode*ncore)
		phUsit = phUsit /(nnode*ncore)		
		if(lelecsh) then
			csit_e = csit_e /(nnode*ncore)
			wsit_e = wsit_e /(nnode*ncore)
			psit_e = psit_e /(nnode*ncore)
			pes_e  = pes_e  /(nnode*ncore)

		endif
		
		if(lholesh) then
			csit_h = csit_h /(nnode*ncore)
			wsit_h = wsit_h /(nnode*ncore)
			psit_h = psit_h /(nnode*ncore)
			pes_h  = pes_h  /(nnode*ncore)
		endif		
	

		!!plot !!!!	
		call plot_phQ(nmodes,nqtotf,nsnap,phQsit)
		call plot_phP(nmodes,nqtotf,nsnap,phPsit)
		call plot_phK(nmodes,nqtotf,nsnap,phKsit)
		call plot_phU(nmodes,nqtotf,nsnap,phUsit)
		deallocate(phPsit,phQsit,phKsit,phUsit)
		
		
		if(lelecsh) then
			write(stdout,"(/,A)") "Plotting electron non-adiabatic dynamica Information."
			call plot_pes(nefre,nsnap,pes_one_e,pes_e,wsit_e,pes_e_file)
			call plot_csit(nefre,nsnap,naver,csit_e,csit_e_file)
			call plot_wsit(nefre,nsnap,naver,wsit_e,wsit_e_file)
			call plot_psit(nefre,nsnap,naver,psit_e,psit_e_file)
			call plot_band_occupatin_withtime(neband,nktotf,Enk_e,xkf,nsnap,psit_e,csit_e,savedsnap,band_e_file)
			deallocate(pes_e,pes_one_e,csit_e,wsit_e,psit_e)
			
		endif
		
		if(lholesh) then
			write(stdout,"(/,A)") "Plotting hole non-adiabatic dynamica Information."
			call plot_pes(nhfre,nsnap,pes_one_h,pes_h,wsit_h,pes_h_file)
			call plot_csit(nhfre,nsnap,naver,csit_h,csit_h_file)
			call plot_wsit(nhfre,nsnap,naver,wsit_h,wsit_h_file)
			call plot_psit(nhfre,nsnap,naver,psit_h,psit_h_file)
			call plot_band_occupatin_withtime(nhband,nktotf,Enk_h,xkf,nsnap,psit_h,csit_h,savedsnap,band_h_file)
			deallocate(pes_h,pes_one_h,csit_h,wsit_h,psit_h)
		endif
		
	endif
  
	
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !% Write End information          %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  call get_date_and_time(cdate,ctime)
  write(stdout,"(/,1X,A77)") repeat("=",77)
  write(stdout,'(1X,"Program LVCSH stoped on ",A9," at ",A9)') cdate,ctime  
  
  call cpu_time(t1)
  !write(6,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'
  write(stdout,'(a,f10.2,a)') 'total time is',(t1-t0)/3600,'hours'  
  
  close(stdout)
  
  stop
  
  
end program lvcsh


