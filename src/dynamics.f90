module dynamics
  use constants
  use io
  use parameters
  use hamiltonian
  use randoms
  use surfacehopping
  
  implicit none
  contains
  
  !=======================================================================!
  != rk4 method to obtain coordinate and velocitie after a time interval =!
  !=======================================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods               =!
  !=======================================================================!

  subroutine rk4_nuclei(nfreem,xx,vv,tt)
    implicit none
    integer,intent(in)         :: nfreem
    real(kind=dp),intent(inout):: xx(1:nfreem),vv(1:nfreem)
    real(kind=dp),intent(in)   :: tt
    real(kind=dp):: tt2,tt6
    real(kind=dp):: xx0(1:nfreem),dx1(1:nfreem),dx2(1:nfreem),dx3(1:nfreem),dx4(1:nfreem)
    real(kind=dp):: vv0(1:nfreem),dv1(1:nfreem),dv2(1:nfreem),dv3(1:nfreem),dv4(1:nfreem)

    tt2=tt/2.0d0; tt6=tt/6.0d0
    call derivs_nuclei(xx,vv,dx1,dv1)
    xx0=xx+tt2*dx1; vv0=vv+tt2*dv1
    call derivs_nuclei(xx0,vv0,dx2,dv2)
    xx0=xx+tt2*dx2; vv0=vv+tt2*dv2
    call derivs_nuclei(xx0,vv0,dx3,dv3)
    xx0=xx+tt*dx3; vv0=vv+tt*dv3
    call derivs_nuclei(xx0,vv0,dx4,dv4)
    xx=xx+tt6*(dx1+2.0d0*dx2+2.0d0*dx3+dx4)
    vv=vv+tt6*(dv1+2.0d0*dv2+2.0d0*dv3+dv4)
  endsubroutine
  
  !====================================================!
  != calculate derivative of coordinate and velocitie =!
  !====================================================!
  != ref: notebook page 630                           =!
  !====================================================!

  subroutine derivs_nuclei(xx,vv,dx,dv)
    use readphonon
    implicit none

    !integer ibasis
    !real(kind=dp)::dd_elec(1:nbasis,1:nbasis,nfreem),dd_hole(1:nbasis,1:nbasis,nfreem)
    real(kind=dp),intent(in)  ::  xx(1:nfreem),vv(1:nfreem)
    real(kind=dp),intent(out) ::  dx(1:nfreem),dv(1:nfreem)
    !integer :: ifreem
    real(kind=dp)::w2x,dEa_dQ,gv
    do ifreem=1,nfreem      
      w2x = (womiga(ifreem)**2)*xx(ifreem)      
      dEa_dQ = 0.0d0
      if(lelecsh) dEa_dQ = dEa_dQ + dE_dQ_e(isurface_e,ifreem)
      if(lholesh) dEa_dQ = dEa_dQ + dE_dQ_h(isurface_h,ifreem)

      gv = gamma*vv(ifreem)
      !dv(ifreem)=(-k(ifreem)*xx(ifreem)-dEa_dpf)/mass(ifreem)-gamma*vv(ifreem)
      dv(ifreem)= -w2x - dEa_dQ - gv
      dx(ifreem)= vv(ifreem)
      
    enddo
    
  endsubroutine derivs_nuclei
  
  subroutine derivs_electron_diabatic(nbasis,cc,dc,HH)
    use f95_precision
    use blas95
    implicit none
    integer,intent(in)   :: nbasis
    complex(kind=dpc),intent(in) :: cc(nbasis)
    complex(kind=dpc),intent(out):: dc(nbasis)
    real(kind=dp),intent(in)     :: HH(nbasis,nbasis)
    complex(kind=dpc),allocatable:: HH_cmplx(:,:)
    allocate(HH_cmplx(nbasis,nbasis))
    HH_cmplx = HH*cmplx_1
    dc=cmplx_0
    !ih dc = Hcc  !use MKL BLAS
    call gemv(HH_cmplx,cc,dc)
    dc = -dc*cmplx_i
    
    deallocate(HH_cmplx)
    
  endsubroutine derivs_electron_diabatic

  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!

  subroutine rk4_electron_diabatic(nbasis,cc,nn,tt,hh)
    implicit none
    integer,intent(in)              :: nbasis
    complex(kind=dpc),intent(inout) :: cc(nbasis)
    real(kind=dp),intent(out)       :: nn(nbasis)
    real(kind=dp),intent(in)        :: tt
    real(kind=dp),intent(in)        :: hh(nbasis,nbasis)
    real(kind=dp):: tt2,tt6
    real(kind=dp):: sum_nn
    complex(kind=dpc):: cc0(1:nbasis),dc1(1:nbasis),&
                        dc2(1:nbasis),dc3(1:nbasis),dc4(1:nbasis)
    nn=CONJG(cc)*cc
    sum_nn = SUM(nn)
    tt2=tt/2.0d0; tt6=tt/6.0d0
    call derivs_electron_diabatic(nbasis,cc,dc1,hh)
    cc0=cc+tt2*dc1
    nn=CONJG(cc0)*cc0
    sum_nn = SUM(nn)    
    call derivs_electron_diabatic(nbasis,cc0,dc2,hh)
    cc0=cc+tt2*dc2
    nn=CONJG(cc0)*cc0
    sum_nn = SUM(nn)       
    call derivs_electron_diabatic(nbasis,cc0,dc3,hh)
    cc0=cc+tt*dc3
    nn=CONJG(cc0)*cc0
    sum_nn = SUM(nn)       
    call derivs_electron_diabatic(nbasis,cc0,dc4,hh)
    cc=cc+tt6*(dc1+2.0d0*dc2+2.0d0*dc3+dc4)

    nn=CONJG(cc)*cc
    sum_nn = SUM(nn)
    cc = cc/dsqrt(sum_nn)
    nn=CONJG(cc)*cc    
  endsubroutine rk4_electron_diabatic  
  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!

  !subroutine rk4_electron_diabatic(xx,cc,nn,tt,hh)
  !  implicit none
  !
  !  real(kind=dp):: tt,tt2,tt6
  !  real(kind=dp):: xx(1:nfreem),hh(1:nbasis,1:nbasis)
  !  complex(kind=dpc):: cc(1:nbasis),cc0(1:nbasis),dc1(1:nbasis),&
  !                      dc2(1:nbasis),dc3(1:nbasis),dc4(1:nbasis)
  !  real(kind=dp)::nn(num_wann,na1site,na2site)
  !  logical ::  llhole
  !  integer :: ia1,ia2,ia3,iw,n_wann
  !  integer :: ibasis
  !  tt2=tt/2.0d0; tt6=tt/6.0d0
  !  call derivs_electron_diabatic(xx,cc,dc1,hh)
  !  cc0=cc+tt2*dc1
  !  call derivs_electron_diabatic(xx,cc0,dc2,hh)
  !  cc0=cc+tt2*dc2
  !  call derivs_electron_diabatic(xx,cc0,dc3,hh)
  !  cc0=cc+tt*dc3
  !  call derivs_electron_diabatic(xx,cc0,dc4,hh)
  !  cc=cc+tt6*(dc1+2.0d0*dc2+2.0d0*dc3+dc4)
  !  do ia2=0,na2
  !    do ia1=0,na1
  !      do iw=1,num_wann
  !        ibasis=ia2site*na1site*num_wann+ia1site*num_wann+n_wann
  !        nn(n_wann,ia1site+1,ia2site+1)=CONJG(cc(ibasis))*cc(ibasis)
  !      enddo
  !    enddo
  !  enddo
  !  
  !endsubroutine rk4_electron_diabatic
  
  !===============================================!
  != add bath effect to coordinate and velocitie =!
  !===============================================!
  != ref: notebook page 462 and 638              =!
  !===============================================!

  subroutine add_bath_effect(nbasis,nfreem,dd_elec,pp_elec,dd_hole,pp_hole,womiga,tt,xx,vv)
    implicit none
    integer,intent(in)       :: nbasis,nfreem
    real(kind=dp),intent(in) :: dd_elec(1:nbasis,1:nbasis,1:nfreem),&
                                dd_hole(1:nbasis,1:nbasis,1:nfreem)
    real(kind=dp),intent(in) :: pp_elec(1:nbasis,1:nbasis),pp_hole(1:nbasis,1:nbasis)
    real(kind=dp),intent(in) :: womiga(nfreem)
    real(kind=dp),intent(in) :: tt
    real(kind=dp),intent(inout) :: xx(1:nfreem),vv(1:nfreem)    
    
    integer       :: ik1,ik2
    integer       :: ifreem,jbasis
    real(kind=dp):: sigmar
    real(kind=dp):: kk,r1,r2,r3,r4,z1,z2,z3,z4

    !sigmar=dsqrt(2.0d0*kb*temp*gamma*tt/mass) !gamma -> 1/t   dangwei v-> dx/dt
    sigmar=dsqrt(2.0d0*kb*temp*gamma*tt)  !dangwei yu dQ/dt   (mv) yiyang.
    do ifreem=1,nfreem
      kk=womiga(ifreem)**2
      if(lelecsh) then
        do jbasis=1,nbasis
          !calculate d2Ei/dx(ifreem)**2 
          if(jbasis /= isurface_e) then
            do ik1=1,nbasis
              do ik2=1,nbasis
                  kk=kk+2.0d0*dd_elec(jbasis,isurface_e,ifreem)*pp_elec(ik1,isurface_e)*&
                  pp_elec(ik2,jbasis)*Hep(ik1,ik2,ifreem)
              enddo
            enddo
          endif
        enddo
      endif
      if(lholesh) then
        do jbasis=1,nbasis
          !calculate d2Ei/dx(ifreem)**2 
          if(jbasis /= isurface_h) then
            do ik1=1,nbasis
              do ik2=1,nbasis
                  kk=kk-2.0d0*dd_hole(jbasis,isurface_h,ifreem)*pp_hole(ik1,isurface_h)*&
                  pp_hole(ik2,jbasis)*Hep(ik1,ik2,ifreem)
              enddo
            enddo
          endif
        enddo
      endif      
      
      r1=gaussian_random_number_fast(0.0d0,sigmar)
      r2=gaussian_random_number_fast(0.0d0,sigmar)
      r3=gaussian_random_number_fast(0.0d0,sigmar)
      r4=gaussian_random_number_fast(0.0d0,sigmar)
      z1=r1                                                     !dQ/dt
      z2=tt*(r1/2.0d0+r2/sqrt3/2.0d0)                           !dQ
      z3=tt**2*(r1/6.0d0+r2*sqrt3/12.0d0+r3/sqrt5/12.0d0)       !dQdt
      z4=tt**3*(r1/24.0d0+r2*sqrt3/40.0d0+r3/sqrt5/24.0d0+r4/sqrt7/120.0d0) !dQ*(dt**2)
      xx(ifreem)=xx(ifreem)+(z2-gamma*z3+(-kk+gamma**2)*z4)
      vv(ifreem)=vv(ifreem)+(z1-gamma*z2+(-kk+gamma**2)*z3+(2.0d0*gamma*kk-gamma**3)*z4)
    
    enddo
    
  endsubroutine add_bath_effect
  
  subroutine saveresult()
  implicit none
    character(len=maxlen) ::  pes_name,csit_name,wsit_name,xsit_name,psit_name,ksit_name
    integer               ::  pes_unit,csit_unit,wsit_unit,xsit_unit,psit_unit,ksit_unit
    integer :: ibasis,ifreem
    inquire(directory = './result',exist=lexist)
    if (.not. lexist) call system('mkdir ./result')
    
    pes_unit = io_file_unit()
    pes_name = './result/pes_elec.out'
    call open_file(pes_name,pes_unit)
    write(pes_unit,"(2(A12,1X))") " time(fs) ","E_isurface"
    do iaver=1,1
      do isnap=1,nsnap
        write(pes_unit,'(999999(e12.5,1X))') dt*nstep*isnap*Au2fs,(pes_elec(ibasis,isnap,iaver),ibasis=0,nbasis)
      enddo
    enddo
    call close_file(pes_name,pes_unit)

    pes_unit = io_file_unit()
    pes_name = './result/pes_hole.out'
    call open_file(pes_name,pes_unit)
    write(pes_unit,"(2(A12,1X))") " time(fs) ","E_isurface"
    do iaver=1,1
      do isnap=1,nsnap
        write(pes_unit,'(999999(e12.5,1X))') dt*nstep*isnap*Au2fs,(pes_hole(ibasis,isnap,iaver),ibasis=0,nbasis)
      enddo
    enddo
    call close_file(pes_name,pes_unit)
    
    csit_unit = io_file_unit()
    csit_name = './result/csit_elec.out'
    call open_file(csit_name,csit_unit)
    do isnap=1,nsnap
      write(csit_unit,'(999999(e12.5,1X))') dt*nstep*isnap*Au2fs,(csit_elec(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(csit_name,csit_unit)
    
    csit_unit = io_file_unit()
    csit_name = './result/csit_hole.out'
    call open_file(csit_name,csit_unit)
    do isnap=1,nsnap
      write(csit_unit,'(999999(e12.5,1X))') dt*nstep*isnap*Au2fs,(csit_hole(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(csit_name,csit_unit)
    
    wsit_unit = io_file_unit()
    wsit_name = './result/wsit_elec.out'
    call open_file(wsit_name,wsit_unit)
    do isnap=1,nsnap
      write(wsit_unit,'(999999(e12.5,1X))') dt*nstep*isnap*Au2fs,(wsit_elec(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(wsit_name,wsit_unit)

    wsit_unit = io_file_unit()
    wsit_name = './result/wsit_hole.out'
    call open_file(wsit_name,wsit_unit)
    do isnap=1,nsnap
      write(wsit_unit,'(999999(e12.5,1X))') dt*nstep*isnap*Au2fs,(wsit_hole(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(wsit_name,wsit_unit)
    
    psit_unit = io_file_unit()
    psit_name = './result/psit_elec.out'
    call open_file(psit_name,psit_unit)
    do isnap=1,nsnap
      write(psit_unit,'(999999(e12.5,1X))') dt*nstep*isnap*Au2fs,(psit_elec(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(psit_name,psit_unit)
    
    psit_unit = io_file_unit()
    psit_name = './result/psit_hole.out'
    call open_file(psit_name,psit_unit)
    do isnap=1,nsnap
      write(psit_unit,'(999999(e12.5,1X))') dt*nstep*isnap*Au2fs,(psit_hole(ibasis,isnap),ibasis=1,nbasis)
    enddo
    call close_file(psit_name,psit_unit)
    
    xsit_unit = io_file_unit()
    xsit_name = './result/xsit.out'
    call open_file(xsit_name,xsit_unit)
    do isnap=1,nsnap
      write(xsit_unit,'(999999(e12.5,1X))') dt*nstep*isnap*Au2fs,(xsit(ifreem,isnap),ifreem=1,nfreem)
    enddo
    call close_file(xsit_name,xsit_unit)
    
    ksit_unit = io_file_unit()
    ksit_name = './result/ksit.out'
    call open_file(ksit_name,ksit_unit)
    do isnap=1,nsnap
      write(ksit_unit,'(999999(e12.5,1X))') dt*nstep*isnap*Au2fs,(ksit(ifreem,isnap),ifreem=1,nfreem)
    enddo
    call  close_file(ksit_name,ksit_unit)
  
  end subroutine saveresult  
  
  
end module dynamics