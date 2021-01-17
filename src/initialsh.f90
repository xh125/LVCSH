module initialsh
  use kinds     ,only   : dp,dpc
  use constants !,only   : maxlen,twopi
  use parameters
  implicit none
  
  contains
  
  function bolziman(wwomiga,ttemp)
    use io,only :stdout
    implicit none
    real(kind=dp)::wwomiga,ttemp,bolziman
    if(wwomiga == 0) then
      write(stdout,*) "womiga == 0,phonon error"
      stop
    endif
    bolziman=1.0/(exp(wwomiga/(ttemp))-1.0)
    !<nb>=1/(exp{hw/kbT}-1)
  end function   
  
  !=============================================!
  != init coordinate and Normal mode velocitie =!
  !=============================================!
  subroutine init_coordinate_velocity(phQ,phP,wf,temp)
    use randoms
    implicit none
    integer,intent(in)        :: nfreem
    real(kind=dp),intent(in)  :: temp
    real(kind=dp),intent(in)  :: womiga(nfreem)
    logical,intent(in)        :: lrandomsita
    real(kind=dp),intent(out) :: xx(1:nfreem),vv(1:nfreem),xxl(1:nfreem)
    character(len=maxlen),intent(in) :: initnmstat
    integer :: ifreem
    !real(kind=dp),external::gaussian_random_number
    real(kind=dp) :: nb   ! bolziman fenbuzhi
    real(kind=dp) :: sita

    do ifreem=1,nfreem
      nb=bolziman(womiga(ifreem),temp)
      if(trim(adjustl(initnmstat))=="quantum") then
        if(Lrandomsita) then
          xxl(ifreem) = 2*dsqrt((nb+0.5)/womiga(ifreem))
          call random_number(sita)
          sita=sita*twopi
          xx(ifreem) = xxl(ifreem)*cos(sita)
          vv(ifreem) = xxl(ifreem)*womiga(ifreem)*sin(sita)
        else
          xx(ifreem)=gaussian_random_number(0.0d0,dsqrt((nb+0.5)*womiga(ifreem))/womiga(ifreem) )
          vv(ifreem)=gaussian_random_number(0.0d0,dsqrt((nb+0.5)*womiga(ifreem)) )
        endif
      elseif(trim(adjustl(initnmstat))=="class") then
        if(Lrandomsita) then
          xxl(ifreem) = 2*dsqrt(temp)/womiga(ifreem)
          call random_number(sita)
          sita=sita*twopi   
          xx(ifreem) = xxl(ifreem)*cos(sita)
          vv(ifreem) = xxl(ifreem)*womiga(ifreem)*sin(sita)   
        else
          xx(ifreem)=gaussian_random_number(0.0d0,dsqrt(temp)/womiga(ifreem) )
          vv(ifreem)=gaussian_random_number(0.0d0,dsqrt(temp) )
        endif
      endif
      !in the unit of au
    enddo
    
  end subroutine init_coordinate_velocity
  
  subroutine get_initialsurface(iinitehstat,represtation)
    use surfacehopping
    use hamiltonian
    implicit none
    character(len=maxlen),intent(in) :: iinitehstat,represtation
    !integer,intent(out)              :: isurface
    integer :: ibasis
    real(kind=dp) :: flagd,flagr
    real(kind=dp) :: sum_ne,sum_nh
    
    if(trim(adjustl(iinitehstat))== "bk" ) then
      if(trim(adjustl(represtation))== "wfstat") then
        if(lelecsh) then
          lallocate = allocated(C_e)
          if(.not. lallocate) allocate(C_e(nbasis))
          lallocate = allocated(w_e)
          if(.not. lallocate) allocate(w_e(nbasis))          
          lallocate = allocated(n_e)
          if(.not. lallocate) allocate(n_e(nbasis))          
          call bk2cc(initeb,initek,nbasis,C_e)
          N_e = CONJG(C_e)*C_e
          sum_ne = SUM(N_e)
          call CC2isurface(nbasis,P_e,C_e,w_e,isurface_e)
        endif
        
        if(lholesh) then
          lallocate = allocated(C_h)
          if(.not. lallocate) allocate(C_h(nbasis))
          lallocate = allocated(w_h)
          if(.not. lallocate) allocate(w_h(nbasis))          
          lallocate = allocated(n_h)
          if(.not. lallocate) allocate(n_h(nbasis))
          call bk2cc(inithb,inithk,nbasis,C_h)
          n_h = abs(C_h)**2
          sum_nh = SUM(N_h)
          call CC2isurface(nbasis,P_h,C_h,w_h,isurface_h)
        endif        
      elseif(trim(adjustl(represtation))== "blochstat") then
        !call blochstat2cc(init_b,init_k,nbasis,CC)
      endif
      
    elseif(trim(adjustl(iinitehstat))== "wf") then
      !call wf2cc(init_wf,nbasis,CC)
    elseif(trim(adjustl(iinitehstat))== "en") then
      if(lelecsh) then
        lallocate = allocated(C_e)
        if(.not. lallocate) allocate(C_e(nbasis))
        lallocate = allocated(w_e)
        if(.not. lallocate) allocate(w_e(nbasis))          
        lallocate = allocated(n_e)
        if(.not. lallocate) allocate(n_e(nbasis)) 
        call en2isurface(initeEN,nbasis,E_e,P_e,Nocc,isurface_e,w_e,C_e,n_e)
      endif
      if(lholesh) then
        lallocate = allocated(C_h)
        if(.not. lallocate) allocate(C_h(nbasis))
        lallocate = allocated(w_h)
        if(.not. lallocate) allocate(w_h(nbasis))          
        lallocate = allocated(n_h)
        if(.not. lallocate) allocate(n_h(nbasis))
        call en2isurface(inithEN,nbasis,E_h,P_h,Nunocc,isurface_h,w_h,C_h,n_h)
        
      endif              
      sum_ne = SUM(N_e)
      sum_nh = SUM(N_h)
      !call en2isurface(initen,isurface)
    elseif(trim(adjustl(iinitehstat))== "ES") then
      !call es2isurface(inites,isurface)
    endif
  
  end subroutine   

  subroutine bk2cc(init_b,init_k,nnbasis,CC)
    use hamiltonian
    use readwannierfile
    implicit none
    integer,intent(in) :: init_b,init_k,nnbasis
    complex(kind=dpc),intent(out) :: CC(nnbasis)
    character(len=maxlen) :: H0parafile
    call readbandkpoints()
    kpoint = plot_kpoint(:,init_k)
    H0parafile = trim(adjustl(SHROOT_dir))//'wannier/wannier90_hr.dat'
    call getHam(H0parafile)
    call Ham2Hamkprm(kpoint)
    call U_int2CC(init_b,kpoint,num_wann,U_int,nnbasis,CC)
  end subroutine
  
  subroutine U_int2CC(init_b,kpoint,num_wann,UU_int,nnbasis,CC)
    implicit none
    integer,intent(in):: init_b,num_wann,nnbasis
    real(kind=dp),intent(in) :: kpoint(3)
    complex(kind=dpc),intent(in) :: UU_int(num_wann,num_wann)
    complex(kind=dpc),intent(out):: CC(nnbasis)
    real(kind=dp) :: sum2c
    integer:: ia1,ia2,ia3
    integer:: nf
    real(kind=dp)     :: rdotk
    complex(kind=dpc) :: fac 
    real(kind=dp)     :: rvec(3)
    do ia3=0,na3-1
      rvec(3) = ia3
      do ia2=0,na2-1
        rvec(2) = ia2
        do ia1=0,na1-1
          rvec(1) = ia1
          nf=(ia3*na2*na1+ia2*na1+ia1)*num_wann
          rdotk=twopi*dot_product(kpoint(:),rvec(:))
          fac=cmplx(cos(rdotk),sin(rdotk),dp)
          CC(nf+1:nf+num_wann) = UU_int(:,init_b)*fac
        enddo
      enddo
    enddo
    sum2c = SUM(CONJG(CC)*CC)
    CC = CC / sqrt(sum2c)
    
  end subroutine U_int2CC
          
  subroutine CC2isurface(nnbasis,PP,CC,WW,isurface)
    use hamiltonian
    implicit none
    integer,intent(in) :: nnbasis
    real(kind=dp),intent(in) :: PP(nnbasis,nnbasis)
    complex(kind=dpc),intent(in) :: CC(nnbasis)
    complex(kind=dpc),intent(out):: WW(nnbasis)
    integer,intent(out) :: isurface
    real(kind=dp),allocatable    :: n_ww(:)
    integer :: ibasis
    real(kind=dp) :: flagd,flagr
    allocate(n_ww(nnbasis))
    call convert_diabatic_adiabatic(nnbasis,PP,CC,WW)
    call random_number(flagr)
    flagd=0.0d0
    n_ww = CONJG(WW)*WW
    do ibasis=1,nnbasis
      flagd=flagd+n_ww(ibasis)
      if(flagr < flagd) then
        isurface=ibasis
        deallocate(n_ww)
        exit
      endif
    enddo
    
  end subroutine CC2isurface
  
  subroutine en2isurface(initEN,nnbasis,EE,PP,Numexit,isurface,WW,CC,nn)
    use hamiltonian
    implicit none
    integer,intent(in)       :: nnbasis
    real(kind=dp),intent(in) :: initEN
    real(kind=dp),intent(in) :: EE(nnbasis),PP(nnbasis,nnbasis)
    integer,intent(in)       :: Numexit
    integer,intent(out)      :: isurface
    complex(kind=dpc),intent(out) ::WW(nnbasis),CC(nnbasis)
    real(kind=dp),intent(out)     :: nn(nnbasis)
    integer :: ibasis
    do ibasis = Numexit+1,nnbasis
      if(EE(ibasis)-EE(Numexit+1)>=initEN) then
        if(abs(EE(ibasis)-EE(Numexit+1))<abs(EE(ibasis-1)-EE(Numexit+1))) then
          isurface = ibasis
        else
          isurface = ibasis-1
        endif
        exit
      endif
    enddo
    WW = cmplx_0
    WW(isurface) = cmplx_1
    CC = PP(:,isurface)
    nn = CONJG(CC)*CC

  end subroutine en2isurface
  
  
end module initialsh 