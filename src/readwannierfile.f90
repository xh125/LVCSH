module readwannierfile
  use kinds,only      : dp,dpc
  use io
  use constants
  use parameters !,only : pj_adjust,na1,na2,na3
  implicit none
  !character(len=maxlen) :: Hr_name
  integer               :: Hr_unit    
  character(len=maxlen) :: header
  integer :: num_wann
  integer :: nrpts,irpt,ir1,ir2,ir3  
  integer,allocatable   :: ndegen(:)
  integer,allocatable   :: irvec(:,:)
  real(kind=dp)         :: ReH, ImH  
  complex(kind=dpc),allocatable :: Ham_r(:,:,:),Ham_r_0(:,:,:)
  complex(kind=dpc),allocatable :: Ham_r_dQ(:,:,:,:)
  
  complex(kind=dpc),allocatable :: Ham_kprm(:,:),U_int(:,:)
  real(kind=dp),allocatable     :: E_k(:)
  !通过对实空间紧束缚模型哈密顿量进行傅里叶变换得到的
  !某一k点处的哈密顿量，以及相应k点的本征值与本征矢
                                   

  real(kind=dp),allocatable     :: Rwann(:,:)
  !logical :: lallocate = .False.
  integer :: total_pts
  real(kind=dp),allocatable::  plot_kpoint(:,:)
  real(kind=dp)            ::  kpoint(3)
  real(kind=dp),allocatable::  bands_projs(:,:,:)
  !用于能带计算的高对称k点
  contains
  
  subroutine set_num_wann(Hr_name,nnum_wann)
    implicit none
  	
    character(len=maxlen),intent(in) :: Hr_name
    integer,intent(out)              :: nnum_wann
    logical            :: lexist
    character(len=255) :: ctmp
    inquire(directory = trim(adjustl(SHROOT_dir))//'wannier',exist=lexist)
    if(lexist) then
      Hr_unit = io_file_unit()
      call open_file(Hr_name,Hr_unit)
      read(Hr_unit, *) ctmp
      read(Hr_unit, *) nnum_wann
      !read(Hr_unit, *) nrpts
      !allocate(ndegen(nrpts),irvec(3,nrpts))    
      !allocate(Ham_r(num_wann,num_wann,nrpts))      
      rewind(Hr_unit)
      call close_file(Hr_name,Hr_unit)
    endif
    !nbasis = na1site*na2site*num_wann
    write(stdout,*) "NUM_wann=",num_wann
    
  end subroutine set_num_wann
  
  !====================================================!
  != read wannier90_hr.dat 文件 得到 irvec(3,nrpts) ===!
  != 以及 ham_r(num_wann,num_wann,nrpts)            ===!
  !====================================================!  
  subroutine readwannhr(Hr_name)
    implicit none
    character(len=maxlen),intent(in) :: Hr_name
    integer :: n_wann,m_wann,m,n
    
    Hr_unit  = io_file_unit()
    call open_file(Hr_name,Hr_unit)    
  
    read(Hr_unit, *) header
    read(Hr_unit, *) num_wann
    read(Hr_unit, *) nrpts
    
    lallocate = ALLOCATED(ndegen)
    if(.not. lallocate) allocate(ndegen(nrpts))
    lallocate = Allocated(irvec)
    if(.not. lallocate) allocate(irvec(3,nrpts))
    lallocate = Allocated(Ham_r)
    if(.not. lallocate) allocate(Ham_r(num_wann,num_wann,nrpts))
    ndegen = 0
    irvec  = 0
    Ham_r  = cmplx_0

    read(Hr_unit, *) (ndegen(irpt), irpt=1, nrpts)
    !read <m0|H|nR>
    !in wannier90_hr.dat R,m,n,HR,HI
    do irpt=1, nrpts       !R
      do n_wann=1, num_wann   !n
        do m_wann=1, num_wann  !m
          read(Hr_unit,*) ir1, ir2, ir3, m, n, ReH, ImH
          irvec(1,irpt) = ir1
          irvec(2,irpt) = ir2
          irvec(3,irpt) = ir3
          Ham_r(m,n,irpt) = dcmplx(ReH,ImH)
        enddo
      enddo
    enddo
    call close_file(Hr_name,Hr_unit)
    
  end subroutine readwannhr
  
  !==============================================================!
  != write wannier90_hr.dat 文件 根据 num_wann,irvec(3,nrpts) ===!
  != 以及 ham_r(num_wann,num_wann,nrpts) 写入文件Hr_name中    ===!
  !==============================================================!  
  subroutine writewannhr(Hr_name,num_wann,nrpts,ndegen,irvec,ham_r)
    implicit none
    character(len=maxlen),intent(in) :: Hr_name
    integer,intent(in)               :: num_wann
    integer,intent(in)               :: nrpts
    integer,intent(in)               :: ndegen(nrpts)
    integer,intent(in)               :: irvec(3,nrpts)
    complex(kind=dp),intent(in)      :: ham_r(num_wann,num_wann,nrpts)
    character(len=9)  :: cdate,ctime 
    integer :: i,j
    logical :: lexist
    Hr_unit = io_file_unit()
    call open_file(Hr_name,Hr_unit)
    call io_date(cdate,ctime)
    header = 'written on '//cdate//' at '//ctime
    write(Hr_unit,*) header
    write(Hr_unit,*) num_wann
    write(Hr_unit,*) nrpts
    write(Hr_unit,'(15I5)') (ndegen(irpt),irpt=1,nrpts)
    do irpt=1,nrpts
      do i=1,num_wann
        do j=1,num_wann
          write(Hr_unit,'(5I5,2F12.6)') irvec(:,irpt),j,i,ham_r(j,i,irpt)
        enddo
      enddo
    enddo
    call close_file(Hr_name,Hr_unit)
    
  end subroutine
  
  !================================================!
  != read "wannier90.wout" setting Rwann(3,nwann) =!
  ! in atoms unit                                  !
  !================================================!
  subroutine readWFcentre(wannier_centres_name)
    use kinds,only:dp
    use io 
    implicit none
    character(len=maxlen),intent(in) :: wannier_centres_name
    integer::i,wannier_centres_unit,iwann
    character(len=maxlen)::ctmp
    !wannier_centres_name = './wannier/wannier90_centres.xyz'    
    wannier_centres_unit = io_file_unit()
    call open_file(wannier_centres_name,wannier_centres_unit)
    allocate(Rwann(3,num_wann))    
    read(wannier_centres_unit,"(/,A)") ctmp
    do iwann=1,num_wann
      read(wannier_centres_unit,*) ctmp,(Rwann(i,iwann),i=1,3)
    enddo
    call close_file(wannier_centres_name,wannier_centres_unit)
    
    Rwann = Rwann/au2ang
    
  end subroutine readWFcentre  
  

  subroutine readbandkpoints()
    implicit none
    integer :: bndunit,loop_spts,i
    character(len=maxlen) :: bndname
    bndname = trim(adjustl(SHROOT_dir))//"wannier/wannier90_band.kpt"
    bndunit = io_file_unit()
    call open_file(bndname,bndunit)
    read(bndunit,*) total_pts
    lallocate = allocated(plot_kpoint)
    if (.not. lallocate) allocate(plot_kpoint(3,total_pts))
    do loop_spts=1,total_pts
      read(bndunit,"(3F12.6)") (plot_kpoint(i,loop_spts),i=1,3)
    end do  
  end subroutine readbandkpoints
  
  !================================================!
  !=get the project of line band on wannier orbital!
  !================================================!
  subroutine readBandProjs()
    use constants
    use io
    implicit none
    integer::total_pts,ipt   !num of kpointer in band line
    real(kind=dp),allocatable ::  xval(:)   !kpointer lable in x value
    real(kind=dp),allocatable ::  eig_int(:,:) !line band
    integer:: i , j
    integer:: bandproj_unit
    character(len=maxlen):: bandproj_name
    
    bandproj_unit = io_file_unit()
    bandproj_name = trim(adjustl(SHROOT_dir))//"wannier/wannier90_bandproj.dat"
    call open_file(bandproj_name,bandproj_unit)
    read(bandproj_unit,*) ctmp,itmp,ctmp,total_pts
    
    allocate(xval(total_pts))
    allocate(eig_int(total_pts,num_wann))
    allocate(bands_projs(total_pts,num_wann,num_wann))
    ! benzheng tai zai wannier orbital shang de tou ying.
    
    write(ctmp,*) num_wann + 2
    ctmp = '('//trim(adjustl(ctmp))//'E16.8)'
    
    do i=1,num_wann
      do ipt=1,total_pts
        read(bandproj_unit,ctmp) xval(ipt),eig_int(ipt,i),(bands_projs(ipt,i,j),j=1,num_wann)
      enddo
    enddo
    
    call close_file(bandproj_name,bandproj_unit)
    deallocate(xval,eig_int)
  
  end subroutine readBandProjs  

  
end module readwannierfile