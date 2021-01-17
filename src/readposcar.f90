module readposcar
  use kinds     ,only :dp
  use constants ,only : maxlen,au2ang
  implicit none
    character(len=maxlen)           :: poscarcomment,Charac
    real(kind=dp),     public, save :: real_lattice(3,3) 
    real(kind=dp),     public, save :: a_lattice,b_lattice,c_lattice 
    real(kind=dp),     public, save :: cell_volume
    real(kind=dp),     public, save :: recip_lattice(3,3)
    integer                          ,public, save :: num_species
    real(kind=dp)        ,allocatable,public, save :: atoms_pos_frac(:,:,:)
    real(kind=dp)        ,allocatable,public, save :: atoms_pos_cart(:,:,:)
    integer              ,allocatable,public, save :: atoms_species_num(:)  
    character(len=maxlen),allocatable,public, save :: atoms_label(:)
    character(len=2)     ,allocatable,public, save :: atoms_symbol(:)
    integer                          ,public, save :: num_atoms,iatom  
    
    contains
    
  subroutine getPOSCAR(dir)
  use io,        only : stdout,io_file_unit,io_error,open_file,close_file
  use utility,   only : utility_recip_lattice
  implicit none      
    character(len=maxlen),intent(in) :: dir
    real(kind=dp)::scaling
    integer::i,j    
    integer::poscar_unit
    character(len=maxlen)::poscar_name
    character(len=maxlen) :: ctmp
    integer :: ierr
    integer :: maxatomnum(1)
    
    poscar_name  = trim(adjustl(dir))//'wannier/POSCAR'
    poscar_unit  = io_file_unit()
    call open_file(poscar_name,poscar_unit)
    read(poscar_unit,*) poscarcomment
    read(poscar_unit,"(F16.8)") scaling
    read(poscar_unit,'(3F20.12)') ((real_lattice(I,J),J=1,3),I=1,3)
    real_lattice = real_lattice*scaling
    a_lattice = sqrt(Sum(real_lattice(1,:)**2))
    b_lattice = sqrt(Sum(real_lattice(2,:)**2))
    c_lattice = sqrt(Sum(real_lattice(3,:)**2))
    call utility_recip_lattice (real_lattice,recip_lattice,cell_volume)
    
    num_atoms    = 0
    num_species  = 1
    read(poscar_unit,*)
    do while(ierr <= 0)
      allocate(atoms_species_num(num_species))
      read(poscar_unit,FMT=*,iostat=ierr) (atoms_species_num(i),i=1,num_species)
      num_species = num_species + 1
      deallocate (atoms_species_num)
      backspace(poscar_unit)
    end do
    num_species = num_species - 2
    
    allocate(atoms_symbol(num_species))
    allocate(atoms_species_num(num_species))
    backspace(poscar_unit)
    backspace(poscar_unit)
    read(poscar_unit,FMT=*,iostat=ierr) (atoms_symbol(i),i=1,num_species)
    read(poscar_unit,FMT=*,iostat=ierr) (atoms_species_num(i),i=1,num_species)
    num_atoms = sum(atoms_species_num)
    maxatomnum = MAXVAL(atoms_species_num)
    allocate(atoms_pos_cart(3,maxatomnum(1),num_species))
    read(poscar_unit,*) Charac
    do i=1,num_species
      do j=1,atoms_species_num(i)
        read(poscar_unit,*) atoms_pos_cart(:,j,i)
      enddo
    enddo
    
    call close_file(poscar_name,poscar_unit)
    write(stdout,"(/,1X,A67)") repeat("=",67)
    write(stdout,"(T5,A)")     "POSCAR information"
    write(stdout,"(1X,A67)")   repeat("=",67)
    write(stdout,*) "Lattice :"
    write(stdout,"(3F20.12)") ((real_lattice(I,J),J=1,3),I=1,3)
    write(stdout,*) "NUM_atoms=",num_atoms
    write(stdout,*) "Atoms_symbol:",(atoms_symbol(i),i=1,num_species)
    write(stdout,*) "Atoms of each symbol:",(atoms_species_num(i),i=1,num_species)
    write(stdout,"(/,1X,A)") "Atoms positions:"
    do i=1,num_species
      do j=1,atoms_species_num(i)
        write(stdout,*) atoms_symbol(i),atoms_pos_cart(:,j,i)
      enddo
    enddo    
    write(stdout,"(1X,A67)")   repeat("=",67)   
    
    real_lattice = real_lattice/au2ang
    a_lattice    = a_lattice/au2ang
    b_lattice    = b_lattice/au2ang
    c_lattice    = c_lattice/au2ang
    
  end subroutine getPOSCAR
  
end module readposcar