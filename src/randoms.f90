module randoms 
  implicit none
  contains
  !===================================!
  != a better version of random_seed =!
  !===================================!

  subroutine init_random_seed()
    implicit none

    integer             :: ii,nn,value_tmp(1:8)
    integer,allocatable :: seed(:)
    real(kind=8)        :: flagd

    call random_seed(size=nn)
    allocate(seed(nn))
    call date_and_time(values=value_tmp)
    seed = value_tmp(8)+37*(/(ii-1,ii=1,nn)/)
    call random_seed(put=seed)
    deallocate(seed)

    do ii=1,value_tmp(6)*3600+value_tmp(7)*60+value_tmp(8)
      call random_number(flagd)
    enddo
  end subroutine init_random_seed

  !==================================!
  != make random number more random =!
  !==================================!

  subroutine more_random()
    implicit none

    integer     :: ii,value_tmp(1:8)
    real(kind=8):: flagd

    call date_and_time(values=value_tmp)
    do ii=1,value_tmp(8)/100
      call random_number(flagd)
    enddo
  end subroutine more_random

  !==============================================================!
  != gaussian random number generator using box-muller method   =!
  !==============================================================!
  != ref: http://en.wikipedia.org/wiki/gaussian_random_variable =!
  !==============================================================!

  function gaussian_random_number(mean,sigma)
    implicit none

    real(kind=8) :: gaussian_random_number,mean,sigma,pi,r1,r2

    pi=4.0d0*datan(1.0d0)
    call more_random()
    call random_number(r1)
    call more_random()
    call random_number(r2)
    gaussian_random_number=mean+sigma*dsqrt(-2.0d0*dlog(r1))*dcos(2.0d0*pi*r2)
  endfunction

  function gaussian_random_number_fast(mean,sigma)
    implicit none

    real(kind=8) :: gaussian_random_number_fast,mean,sigma,pi,r1,r2

    pi=4.0d0*datan(1.0d0)
    call random_number(r1)
    call random_number(r2)
    gaussian_random_number_fast=mean+sigma*dsqrt(-2.0d0*dlog(r1))*dcos(2.0d0*pi*r2)
  endfunction


end module randoms