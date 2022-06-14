module grid
  use kinds,only : dp
  implicit none
  contains
  
  subroutine loadkmesh()
    use epwcom,   only : filkf, nkf1,nkf2,nkf3
    use elph2 ,   only : nkqtotf, nkqf, xkf, wkf, nkf,xkfd,xkf_irr,wkf_irr
    use cell_base,only : at,bg
    use symm_base,only : s, t_rev, time_reversal, nrot,nsym
    use io ,      only : io_file_unit 
    implicit none
    !
    character(len=10) :: coordinate_type
    !! filkf coordinate type (crystal or cartesian)
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: idir
    !! Crystal direction (G-vector) 
    integer :: i,j,k
    integer :: iunkf,iunqf
    integer :: ierr
    real(kind=dp), allocatable :: xkf_(:,:)
    !! coordinates k-points
    real(kind=dp), allocatable :: wkf_(:)
    !! weights k-points
    !
    if(filkf /='') then ! load from file
      iunkf = io_file_unit()
      open(unit=iunkf,file=filkf)
      read(iunkf,*) nkqtotf, coordinate_type
      allocate(xkf_(3,2*nkqtotf))
      allocate(wkf_(2*nkqtotf))
      do ik=1,nkqtotf
        ikk = 2* ik -1
        ikq = ikk+1
        read(iunkf,*) xkf_(:,ikk),wkf_(ikk)
        !
        ! SP: This is so we can input a weight of 1 to random file
        !     This way you can feed the same file for the k and q grid
        wkf_(ikk) = wkf_(ikk) * 2.d0
        xkf_(:, ikq) = xkf_(:, ikk)
        wkf_(ikq) = 0.d0
        !
      enddo
      close(iunkf)
      !
      ! redefine nkqtotf to include the k+q points
      !
      nkqtotf = 2 * nkqtotf       
    elseif((nkf1 /= 0) .and. (nkf2 /= 0) .and. (nkf3 /= 0)) then ! generate grid
      ! mp_mesh_k
      !WRITE(stdout, '(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3        
      !read(unitepwout,"(27X,3i4)") nqf1,nqf2,nqf3
      nkqtotf = 2 * nkf1 * nkf2 * nkf3
      ALLOCATE(xkf_(3, nkqtotf), STAT = ierr)
      !IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating xkf_ ', 1)
      ALLOCATE(wkf_(nkqtotf), STAT = ierr)
      !IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_', 1)
      wkf_(:) = 0.d0
      DO ik = 1, nkf1 * nkf2 * nkf3
        wkf_(2 * ik - 1) = 2.d0 / DBLE(nkqtotf / 2)
      ENDDO
      DO i = 1, nkf1
        DO j = 1, nkf2
          DO k = 1, nkf3
            ik = (i - 1) * nkf2 * nkf3 + (j - 1) * nkf3 + k
            ikk = 2 * ik - 1
            ikq = ikk + 1
            xkf_(1, ikk) = DBLE(i - 1) / DBLE(nkf1)
            xkf_(2, ikk) = DBLE(j - 1) / DBLE(nkf2)
            xkf_(3, ikk) = DBLE(k - 1) / DBLE(nkf3)
            xkf_(1, ikq) = xkf_(1, ikk)
            xkf_(2, ikq) = xkf_(2, ikk)
            xkf_(3, ikq) = xkf_(3, ikk)
          ENDDO
        ENDDO
      ENDDO
    ENDIF !mp_mesh_k        
    
    nkqf = nkqtotf
    nkf = nkqf / 2
    
    allocate(xkf(3,nkqf))
    allocate(wkf(nkqf))
    xkf = xkf_
    
    deallocate(xkf_,wkf_)
    
  end subroutine loadkmesh
  
  subroutine loadkmesh_fullBZ()
  !-----------------------------------------------------------------------
  !!
  !!  Create a k-mesh for fine grid without symmetries on the full grid
  !!
  !-----------------------------------------------------------------------
    use kinds,    only: dp
    use epwcom,   only: nkf1,nkf2,nkf3
    use elph2,    only: xkf_bz
    implicit none
    
    integer :: ik,i,j,k
    
    allocate(xkf_bz(3,nkf1*nkf2*nkf3))
    xkf_bz = 0.0
    IF ((nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0)) THEN
      DO i = 1, nkf1
        DO j = 1, nkf2
          DO k = 1, nkf3
            ik = (i - 1) * nkf2 * nkf3 + (j - 1) * nkf3 + k
            xkf_bz(1, ik) = DBLE(i - 1) / DBLE(nkf1)
            xkf_bz(2, ik) = DBLE(j - 1) / DBLE(nkf2)
            xkf_bz(3, ik) = DBLE(k - 1) / DBLE(nkf3)
          ENDDO
        ENDDO
      ENDDO
    ENDIF    
  end subroutine loadkmesh_fullBZ
  
  subroutine loadqmesh()
  use kinds,    only : dp
  use epwcom,   only : filqf,nqf1,nqf2,nqf3
  use elph2,    only : xqf,wqf,nqf,nqtotf
  use io ,      only : io_file_unit
  implicit none
  integer :: ierr
  integer :: iunqf
  character(len=10) :: coordinate_type
  integer :: iq
  real(kind=dp),allocatable :: xqf_(:,:)
  real(kind=dp),allocatable :: wqf_(:)
  integer :: i,j ,k
  
  filqf = ''
  if(filqf /='') then ! load from file
    iunqf = io_file_unit()
    open(unit=iunqf,file=filqf)
    read(iunqf,*) nqtotf, coordinate_type
    allocate(xqf_(3,nqtotf))
    allocate(wqf_(nqtotf))
    do iq=1,nqtotf
      read(iunqf,*) xqf_(:,iq),wqf_(iq)
    enddo
    close(iunqf)     
  elseif((nqf1 /= 0) .and. (nqf2 /= 0) .and. (nqf3 /= 0)) then ! generate grid
    ! mp_mesh_k
    !WRITE(stdout, '(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3        
    !read(unitepwout,"(27X,3i4)") nqf1,nqf2,nqf3
    nqtotf = nqf1*nqf2*nqf3
    ALLOCATE(xqf_(3, nqtotf), STAT = ierr)
    ALLOCATE(wqf_(nqtotf), STAT = ierr)
    !IF (ierr /= 0) CALL errore('loadkmesh_para', 'Error allocating wkf_', 1)
    wqf_(:) = 0.d0
    wqf_ = 1.0 / dble(nqtotf)
    DO i = 1, nqf1
      DO j = 1, nqf2
        DO k = 1, nqf3
          iq = (i - 1) * nqf2 * nqf3 + (j - 1) * nqf3 + k
          xqf_(1, iq) = DBLE(i - 1) / DBLE(nqf1)
          xqf_(2, iq) = DBLE(j - 1) / DBLE(nqf2)
          xqf_(3, iq) = DBLE(k - 1) / DBLE(nqf3)             
        ENDDO
      ENDDO
    ENDDO
  ENDIF !mp_mesh_k        
    
  nqf = nqtotf
  
  allocate(xqf(3,nqf))
  allocate(wqf(nqf))
  xqf = xqf_
  
  deallocate(xqf_,wqf_)    
  
  end subroutine loadqmesh
  
  integer function get_ikq(xk,xq)
    use kinds,        only : dp
    use epwcom,       only : nkf1,nkf2,nkf3
    use constants,    only : eps3
    use io,           only : stdout
    implicit none
    !integer :: get_ikq
    !! the index of k+q
    real(kind=dp),intent(in) :: xk(3)
    real(kind=dp),intent(in) :: xq(3)
    real(kind=dp) :: xx, yy, zz
    real(kind=dp) :: xkq(3)
    logical :: in_the_list 
    !! Check if k point is in the list      
    
    xkq = xk + xq
    xx = xkq(1)*nkf1
    yy = xkq(2)*nkf2
    zz = xkq(3)*nkf3
    in_the_list = ABS(xx - NINT(xx)) <= eps3 .AND. &
                  ABS(yy - NINT(yy)) <= eps3 .AND. &
                  ABS(zz - NINT(zz)) <= eps3      
    if(.not. in_the_list) then
      write(stdout,*) 'k+q does not fall on k-grid'
    endif
    CALL backtoBZ(xx, yy, zz, nkf1, nkf2, nkf3)
    !
    ! since k- and q- meshes are commensurate, nkq can be easily found
    !
    get_ikq = NINT(xx) * nkf2 * nkf3 + NINT(yy) * nkf3 + NINT(zz) + 1
    !
    !  Now nkq represents the index of k+sign*q on the fine k-grid.
    !
    RETURN 
    !      
    
    
    
  end function get_ikq
  
  !---------------------------------
  SUBROUTINE backtoBZ(xx, yy, zz, n1, n2, n3)
  !---------------------------------
  !!
  !!  Brings xx, yy, and zz  into first BZ
  !!
  !---------------------------------
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: n1, n2, n3
  !! cell size
  REAL(KIND = DP), INTENT(inout) :: xx, yy, zz
  !! kgrid
  !
  ! Local variables
  INTEGER :: ib
  !! Size of replicas
  !
  ! More translations are needed to go back to the first BZ when the unit cell
  ! is far from cubic
  !
  DO ib = -2, 0
    IF (NINT(xx) < ib * n1) xx = xx + (-ib + 1) * n1
    IF (NINT(yy) < ib * n2) yy = yy + (-ib + 1) * n2
    IF (NINT(zz) < ib * n3) zz = zz + (-ib + 1) * n3
  ENDDO
  DO ib = 2, 1, -1
    IF (NINT(xx) >= ib * n1) xx = xx - ib * n1
    IF (NINT(yy) >= ib * n2) yy = yy - ib * n2
    IF (NINT(zz) >= ib * n3) zz = zz - ib * n3
  ENDDO
  !
  !-------------------------------------------
  END SUBROUTINE backtoBZ
  !-------------------------------------------    
  
  subroutine kq2k_map()
    use kinds, only :dp
    use epwcom, only : kqmap,nkf1,nkf2,nkf3,nqf1,nqf2,nqf3
    implicit none
    integer :: iqx,iqy,iqz
    integer :: iq
    integer :: ikx,iky,ikz
    integer :: ik
    integer :: ikqx,ikqy,ikqz
    integer :: ikq
    if (.not. allocated(kqmap)) allocate(kqmap(nkf1*nkf2*nkf3,nqf1*nqf2*nqf3))
    do iqx=1,nqf1
      do iqy=1,nqf2
        do iqz=1,nqf3
          iq= (iqx-1)*nqf2*nqf3+(iqy-1)*nkf3+iqz
          do ikx=1,nkf1
            do iky=1,nkf2
              do ikz=1,nkf3
                ik = (ikx-1)*nkf2*nkf3+(iky-1)*nkf3+ikz
                ikqx = ikx + iqx -1
                ikqy = iky + iqy -1
                ikqz = ikz + iqz -1
                if (ikqx > nkf1) ikqx = ikqx - nkf1
                if (ikqy > nkf2) ikqy = ikqy - nkf2
                if (ikqz > nkf3) ikqz = ikqz - nkf3
                ikq = (ikqx-1)*nkf2*nkf3+(ikqy-1)*nkf3+ikqz
                kqmap(ik,iq) = ikq
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
              
  end subroutine kq2k_map
  
  subroutine mq2q_map()
    use kinds,only : dp
    use epwcom, only : mqmap,nqf1,nqf2,nqf3
    implicit none
    integer :: iqx,iqy,iqz,iqx_,iqy_,iqz_
    integer :: iq,imq
    if(.not. allocated(mqmap)) allocate(mqmap(nqf1*nqf2*nqf3))
    do iqx = 0,nqf1-1
      iqx_ = -1 * iqx
      if(iqx_<0) iqx_ = iqx_ + nqf1
      do iqy = 0,nqf2-1
        iqy_= -1 * iqy
        if(iqy_<0) iqy_ = iqy_ + nqf2
        do iqz = 0,nqf3-1
          iqz_ = -1*iqz
          if(iqz_<0) iqz_ = iqz_ + nqf3
          iq=iqx*nqf2*nqf3+iqy*nqf3+iqz+1
          imq=iqx_*nqf2*nqf3+iqy_*nqf3+iqz_+1
          mqmap(iq)=imq
        enddo
      enddo
    enddo
    
  end subroutine mq2q_map
  
  
end module grid
    