module dynamics
  use kinds,only : dp,dpc
  use io,only :
  use parameters,only : gamma
  use hamiltonian,only : nphfre,nefre
  use elph2,only          : nktotf,nqtotf,wf,nbndfst,epcq
  use modes,only          : nmodes
  use epwcom,only         : kqmap
  use surfacehopping,only : iesurface,ihsurface
  
  implicit none
  contains
  
  !=======================================================================!
  != rk4 method to obtain coordinate and velocitie after a time interval =!
  !=======================================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods               =!
  !=======================================================================!

  subroutine rk4_nuclei(pp_nk,xx,vv,tt)
    implicit none
    real(kind=dp),intent(in)   :: pp_nk(nbndfst,nktotf,nefre)
    real(kind=dp),intent(inout):: xx(nmodes,nqtotf),vv(nmodes,nqtotf)
    real(kind=dp),intent(in)   :: tt
    real(kind=dp):: tt2,tt6
    real(kind=dp):: xx0(nmodes,nqtotf),dx1(nmodes,nqtotf),dx2(nmodes,nqtotf),dx3(nmodes,nqtotf),dx4(nmodes,nqtotf)
    real(kind=dp):: vv0(nmodes,nqtotf),dv1(nmodes,nqtotf),dv2(nmodes,nqtotf),dv3(nmodes,nqtotf),dv4(nmodes,nqtotf)

    tt2=tt/2.0d0; tt6=tt/6.0d0
    call derivs_nuclei(pp_nk,xx,vv,dx1,dv1)
    xx0=xx+tt2*dx1; vv0=vv+tt2*dv1
    call derivs_nuclei(pp_nk,xx0,vv0,dx2,dv2)
    xx0=xx+tt2*dx2; vv0=vv+tt2*dv2
    call derivs_nuclei(pp_nk,xx0,vv0,dx3,dv3)
    xx0=xx+tt*dx3; vv0=vv+tt*dv3
    call derivs_nuclei(pp_nk,xx0,vv0,dx4,dv4)
    xx=xx+tt6*(dx1+2.0d0*dx2+2.0d0*dx3+dx4)
    vv=vv+tt6*(dv1+2.0d0*dv2+2.0d0*dv3+dv4)
  endsubroutine
  
  !====================================================!
  != calculate derivative of coordinate and velocitie =!
  !====================================================!
  != ref: notebook page 630                           =!
  !====================================================!

  subroutine derivs_nuclei(pp_nk,xx,vv,dx,dv)
    implicit none
    real(kind=dp),intent(in)   :: pp_nk(nbndfst,nktotf,nefre)
    real(kind=dp),intent(in)  ::  xx(nmodes,nqtotf),vv(nmodes,nqtotf)
    real(kind=dp),intent(out) ::  dx(nmodes,nqtotf),dv(nmodes,nqtotf)
    integer :: iq,imode,ik,iband1,iband2,ikq
    real(kind=dp)::dEa_dQ
    do iq=1,nqtotf
      do imode=1,nmodes
        dx(imode,iq) = vv(imode,iq)
        dv(imode,iq) = (-wf(imode,iq)**2*xx(imode,iq)-gamma*vv(imode,iq))
        
        dEa_dQ = 0.0
        
        do ik=1,nktotf
          ikq = kqmap(ik,iq)
          do iband1=1,nbndfst
            do iband2=1,nbndfst
              ! 电子能量随简正模的变化
              dEa_dQ = dEa_dQ + &
              pp_nk(iband1,ik,iesurface)*pp_nk(iband2,ikq,iesurface)*epcq(iband1,iband2,ik,imode,iq)
              ! 空穴能量随简正模的变化
              dEa_dQ = dEa_dQ - &
              pp_nk(iband1,ik,ihsurface)*pp_nk(iband2,ikq,iesurface)*epcq(iband1,iband2,ik,imode,iq)                    
            enddo
          enddo
        enddo
        dEa_dQ = sqrt(2.0*wf(imode,iq)/nqtotf) * dEa_dQ
        dv(imode,iq) = dv(imode,iq) - dEa_dQ
      enddo
    enddo
    
  endsubroutine derivs_nuclei
  
  !subroutine derivs_electron_diabatic(nbasis,cc,dc,HH)
  !  use f95_precision
  !  use blas95
  !  implicit none
  !  integer,intent(in)   :: nbasis
  !  complex(kind=dpc),intent(in) :: cc(nbasis)
  !  complex(kind=dpc),intent(out):: dc(nbasis)
  !  real(kind=dp),intent(in)     :: HH(nbasis,nbasis)
  !  complex(kind=dpc),allocatable:: HH_cmplx(:,:)
  !  allocate(HH_cmplx(nbasis,nbasis))
  !  HH_cmplx = HH*cmplx_1
  !  dc=cmplx_0
  !  !ih dc = Hcc  !use MKL BLAS
  !  call gemv(HH_cmplx,cc,dc)
  !  dc = -dc*cmplx_i
  !  
  !  deallocate(HH_cmplx)
  !  
  !endsubroutine derivs_electron_diabatic

  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!

  !subroutine rk4_electron_diabatic(nbasis,cc,nn,tt,hh)
  !  implicit none
  !  integer,intent(in)              :: nbasis
  !  complex(kind=dpc),intent(inout) :: cc(nbasis)
  !  real(kind=dp),intent(out)       :: nn(nbasis)
  !  real(kind=dp),intent(in)        :: tt
  !  real(kind=dp),intent(in)        :: hh(nbasis,nbasis)
  !  real(kind=dp):: tt2,tt6
  !  real(kind=dp):: sum_nn
  !  complex(kind=dpc):: cc0(1:nbasis),dc1(1:nbasis),&
  !                      dc2(1:nbasis),dc3(1:nbasis),dc4(1:nbasis)
  !  nn=CONJG(cc)*cc
  !  sum_nn = SUM(nn)
  !  tt2=tt/2.0d0; tt6=tt/6.0d0
  !  call derivs_electron_diabatic(nbasis,cc,dc1,hh)
  !  cc0=cc+tt2*dc1
  !  nn=CONJG(cc0)*cc0
  !  sum_nn = SUM(nn)    
  !  call derivs_electron_diabatic(nbasis,cc0,dc2,hh)
  !  cc0=cc+tt2*dc2
  !  nn=CONJG(cc0)*cc0
  !  sum_nn = SUM(nn)       
  !  call derivs_electron_diabatic(nbasis,cc0,dc3,hh)
  !  cc0=cc+tt*dc3
  !  nn=CONJG(cc0)*cc0
  !  sum_nn = SUM(nn)       
  !  call derivs_electron_diabatic(nbasis,cc0,dc4,hh)
  !  cc=cc+tt6*(dc1+2.0d0*dc2+2.0d0*dc3+dc4)
  !
  !  nn=CONJG(cc)*cc
  !  sum_nn = SUM(nn)
  !  cc = cc/dsqrt(sum_nn)
  !  nn=CONJG(cc)*cc    
  !endsubroutine rk4_electron_diabatic  
  !===========================================================!
  != rk4 method to obtain wavefunction after a time interval =!
  !===========================================================!
  != ref: http://en.wikipedia.org/wiki/runge_kutta_methods   =!
  !===========================================================!

  !subroutine rk4_electron_diabatic(xx,cc,nn,tt,hh)
  !  implicit none
  !
  !  real(kind=dp):: tt,tt2,tt6
  !  real(kind=dp):: xx(nmodes,nqtotf),hh(1:nbasis,1:nbasis)
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

  
  
end module dynamics