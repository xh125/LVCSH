module write_sh_information
  use kinds,only : dp,dpc
  implicit none
  
  contains
  
  subroutine write_initial_information(iaver,nmodes,nqtotf,wf,phQ,phP)
    use io,only : stdout
    use constants , only : ryd2eV
    use surfacecom, only : lholesh,lelecsh,ihsurface,iesurface, &
                           csit_e,wsit_e,pes_one_e,pes_e,psit_e, &
                           csit_h,wsit_h,pes_one_h,pes_h,psit_h, &
                           c_e,c_h,w_e,w_h
    use hamiltonian, only : E_h,E_e,P_e,P_h,nefre,nhfre
    use surfacehopping,only : PhU,phK,SUM_phK,SUM_phU,SUM_phE,SUM_phK0,&
                              phQsit,phPsit,phKsit,phUsit
    implicit none
    integer, intent(in) :: iaver,nmodes,nqtotf
    real(kind=dp), intent(in) :: wf(nmodes,nqtotf)
    complex(kind=dpc),intent(in) :: phQ(nmodes,nqtotf),phP(nmodes,nqtotf)
    
    integer :: isnap,iq,imode
    
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
      phU = 0.5*(wf**2)*(phQ*CONJG(phQ))
      phK = 0.5*phP*CONJG(phP)        
      SUM_phK = SUM(phK)
      SUM_phU = SUM(phU)
      SUM_phE = SUM_phK+SUM_phU
      
      SUM_phK0 = SUM_phK
      
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      !% Write initial non-adiabatic dynamica information        %!
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      if(lholesh) then
        if(lelecsh) then
          write(stdout,"(/,A)") " time(fs)    rt(s) hsur esur&
          &     E_h(eV)     E_e(eV)    E_eh(eV)    T_ph(eV)    U_ph(eV)    E_ph(eV)   E_tot(eV)"         
          write(stdout,"(F9.2,F9.2,I5,I5,7(1X,F11.4))") 0.00,0.00,&
          ihsurface,iesurface,&
          -e_h(ihsurface)*ryd2eV,e_e(iesurface)*ryd2eV,(e_e(iesurface)+e_h(ihsurface))*ryd2eV,&
          SUM_phK*ryd2eV,SUM_phU*ryd2eV,SUM_phE*ryd2eV,(E_e(iesurface)+E_h(ihsurface)+SUM_phE)*ryd2eV       
        else 
          write(stdout,"(/,A)") " time(fs)    rt(s) hsur&
          &     E_h(eV)    T_ph(eV)    U_ph(eV)    E_ph(eV)   E_tot(eV)"       
          write(stdout,"(F9.2,F9.2,I5,5(1X,F11.4))") 0.00,0.00,&
          ihsurface,-e_h(ihsurface)*ryd2eV,&
          SUM_phK*ryd2eV,SUM_phU*ryd2eV,SUM_phE*ryd2eV,(E_h(ihsurface)+SUM_phE)*ryd2eV     
        endif
      else
        if(lelecsh) then
          write(stdout,"(/,A)") " time(fs)    rt(s) esur&
          &     E_e(eV)    T_ph(eV)    U_ph(eV)     E_ph(eV)   E_tot(eV)" 
          write(stdout,"(F9.2,F9.2,I5,5(1X,F11.4))") 0.00,0.00,&
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
          !do ifre=1,nefre  
            csit_e(:,isnap) = csit_e(:,isnap)+REAL(c_e(:)*CONJG(c_e(:)))
            wsit_e(:,isnap) = wsit_e(:,isnap)+REAL(w_e(:)*CONJG(w_e(:)))
            psit_e(:,isnap) = psit_e(:,isnap)+ABS(P_e(:,iesurface))**2
            pes_e(1:nefre,isnap) = pes_e( 1:nefre,isnap) + E_e(1:nhfre)
            if(iaver == 1) pes_one_e(1:nefre,isnap) = E_e(1:nefre)
          !enddo
        endif
        
        if(lholesh) then
          pes_h(0,isnap) = pes_h(0,isnap)-E_h(ihsurface)
          if(iaver == 1) pes_one_h(0,isnap) = -E_h(ihsurface)
          !do ifre=1,nhfre
            csit_h(:,isnap) = csit_h(:,isnap)+REAL(c_h(:)*CONJG(c_h(:)))
            wsit_h(:,isnap) = wsit_h(:,isnap)+REAL(w_h(:)*CONJG(w_h(:)))
            psit_h(:,isnap) = psit_h(:,isnap)+P_h(:,iesurface)**2
            pes_h(1:nhfre,isnap) = pes_h(1:nhfre,isnap)-E_h(1:nhfre)
            if(iaver == 1) pes_one_h(1:nhfre,isnap) = -E_h(1:nhfre)
          !enddo    
        endif
      
      !=================================!
      != End store initial information =!
      !=================================!      
  end subroutine write_initial_information  

end module write_sh_information