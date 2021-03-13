module ephwann
  use kinds,only : dp
  implicit none
  complex(kind=dp) ,allocatable :: vmef(:,:,:,:)
  !! velocity matrix elements on the fine
  complex(kind=dp) ,allocatable :: dmef(:,:,:,:)
  !! dipole matrix elements on the fine mesh
  logical :: vme
  !! If .true. then calculate the velocity matrix elements beyond the local approximation.
  contains
  subroutine readvmef()
    vme = .false.
    if(vme) then
      allocate(vmef(3,nbndsub,nbndsub,nkqtotf),stat=ierr)
      if(ierr/=0) call errore('ephwann','Error allocating vmef',1)
    else
      allocate(dmef(3,nbndsub,nbndsub,nkqtotf),stat=ierr)
      if(ierr/=0) call errore('ephwann','Error allocating dmef',1)
    endif
    
  end subroutine readvmef
  
end module ephwann