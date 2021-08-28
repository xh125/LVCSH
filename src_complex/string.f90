module string
! this module offers the following routines and functions for string manipulation:

! subroutine cut(string,delimiter,first,second)
!	in	string		inputstring
!	in	delimiter	character at which is cut
!	out	first		part of string before delimiter
!	out	second		part of string after delimiter

! subroutine compact(string,delimiter)
!	inout	string		string to be modified
!	in	delimiter	character to be collapsed

! subroutine tab2blank(string)
!	inout	string		string to be detabbed

! subroutine lowercase(string)
!	inout	string		string to be lowercased

! subroutine uppercase(string)
!	inout	string		string to be uppercased

! subroutine split(instring,delimiter,substrings,n)
!	in	instring	string to be split
!	in	delimiter	character at which is to be split
!	out	substrings	CURRENTLY NOT ALLOCATED array of strings
!	out	n		length of substrings on exit
integer :: maxlen=256

  contains

  ! =================================================================== !
  
  subroutine cut(string,delimiter,first,second)
  ! This routine takes an inputstring and a single character delimiter 
  ! and returns the substrings before and after the first occurence of this 
  ! delimiter in first and second
  implicit none
  character(len=maxlen), intent(out) :: first, second
  character(len=maxlen), intent(in) :: string
  character(len=maxlen) :: dummystring
  character, intent(in) :: delimiter
  integer :: i
  
  dummystring=string
  call tab2blank(dummystring)
  i=index(dummystring,delimiter)
  ! do
  !   if (i==0) exit
  !   if (string(i-1:i-1)=='\') then
  !     i=i+index(string(i+1:),delimiter)
  !   else
  !     exit
  !   endif
  ! enddo
  
  if (i==0) then
    first=dummystring
    second=''
  else
    first=dummystring(1:i-1)
    second=dummystring(i+1:)
  endif
  return
  
  endsubroutine
  
  ! =================================================================== !
  
  subroutine compact(string,delimiter)
  ! this subroutine removes leading delimiters and collapses multiple delimiters to a single one
  implicit none
  character(len=maxlen), intent(inout) :: string
  character, intent(in) :: delimiter
  integer :: i
  character(len=maxlen) :: tempstring,keepstring
  
  i=0
  keepstring=''
  tempstring=string
  
  do
    call generaladjustl(tempstring,delimiter)
    i=index(trim(tempstring),delimiter)
  
    if (i==0) then
      keepstring=trim(keepstring)//delimiter//tempstring
      exit
    endif
  
    keepstring=trim(keepstring)//delimiter//tempstring(1:i-1)
    tempstring=tempstring(i+1:len(tempstring))
  enddo
  
  call generaladjustl(keepstring,delimiter)
  string=keepstring
  
  return
  endsubroutine
  
  ! =================================================================== !
  
  ! recursive subroutine compact(string,delimiter)
  ! ! this subroutine removes leading delimiters and collapses multiple delimiters to a single one
  ! implicit none
  ! character(len=maxlen), intent(inout) :: string
  ! character, intent(in) :: delimiter
  ! integer :: i
  ! character(len=maxlen) :: tempstring,keepstring
  ! 
  ! tempstring=string
  ! call generaladjustl(tempstring,delimiter)
  ! i=index(trim(tempstring),delimiter)
  ! 
  ! if (i==0) then
  !   string=tempstring
  ! else
  !   keepstring=tempstring(1:i-1)
  !   tempstring=tempstring(i+1:len(tempstring))
  !   call compact(tempstring,delimiter)
  !   string=keepstring(1:i-1)//delimiter//tempstring
  ! endif
  ! 
  ! return
  ! endsubroutine
  
  ! =================================================================== !
  
  subroutine tab2blank(string)
  ! this routine converts tab characters to blank characters
  implicit none
  character(len=maxlen), intent(inout) :: string
  integer :: i
  character :: ch
  
  do i=1,len(string)
    ch=string(i:i)
    if (iachar(ch)==9) string(i:i)=' '
  enddo
  return
  
  endsubroutine
  
  ! =================================================================== !
  
  subroutine lowercase(string)
  ! this subroutine converts all uppercase letters A-Z to lowercase a-z
  implicit none
  character(len=maxlen), intent(inout) :: string
  integer :: diff=iachar('A')-iachar('a'), i,ia
  
  do i=1,len(string)
    ia=iachar(string(i:i))
    if ( (iachar('A')<=ia).and.(ia<=iachar('Z')) ) string(i:i)=achar(ia-diff)
  enddo
  return
  
  endsubroutine
  
  ! =================================================================== !
  
  subroutine uppercase(string)
  ! this subroutine converts all lowercase letters a-z to uppercase A-Z
  implicit none
  character(len=maxlen), intent(inout) :: string
  integer :: diff=iachar('A')-iachar('a'), i,ia
  
  do i=1,len(string)
    ia=iachar(string(i:i))
    if ( (iachar('a')<=ia).and.(ia<=iachar('z')) ) string(i:i)=achar(ia+diff)
  enddo
  return
  
  endsubroutine
  
  ! =================================================================== !
  
  subroutine split(instring,delimiter,substrings,n)
  ! this subroutine takes a string and splits it into an array of substrings
  ! it allocates substrings and returns the parts of the string in this array
  implicit none
  character(len=maxlen), intent(out), allocatable :: substrings(:)
  character(len=maxlen), intent(in) :: instring
  character, intent(in) :: delimiter
  integer, intent(out) :: n
  character(len=maxlen) :: string,tempstring
  integer :: i,numdel
  
  string=instring
  call tab2blank(string)
  call compact(string,delimiter)
  ! number of substrings
  numdel=1
  do i=2,strlen(string)
    if (string(i:i)==delimiter) numdel=numdel+1
  enddo
  allocate( substrings(numdel) )
  ! substrings
  call cut(string,delimiter,substrings(1),tempstring)
  do i=2,numdel
    string=tempstring
    call cut(string,delimiter,substrings(i),tempstring)
  enddo
  n=numdel
  return
  
  endsubroutine
  
  ! =================================================================== !
  
  integer function strlen(string)
  implicit none
  character(len=maxlen), intent(in) :: string
  
  strlen=len(trim(string))
  return
  
  endfunction
  
  ! =================================================================== !
  
  subroutine generaladjustl(string,delimiter)
  implicit none
  character(len=maxlen), intent(inout) :: string
  character, intent(in) :: delimiter
  character(len=maxlen) :: tempstring
  integer :: i
  
  do i=1,len(string)
    if (string(i:i)/=delimiter) exit
  enddo
  string=string(i:)
  return
  
  endsubroutine
  
  ! =================================================================== !
  
  ! recursive subroutine generaladjustl(string,delimiter)
  ! implicit none
  ! character(len=maxlen), intent(inout) :: string
  ! character, intent(in) :: delimiter
  ! character(len=maxlen) :: tempstring
  ! 
  ! if (string=='') return
  ! if (string(1:1)==delimiter) then
  !   tempstring=string(2:)
  !   call generaladjustl(tempstring,delimiter)
  !   string=tempstring
  ! endif
  ! return
  ! 
  ! endsubroutine
  
  ! =================================================================== !
  
  subroutine get_quoted(line,stringinquote)
  ! extracts the part of the string between quotes
  implicit none
  character(len=maxlen),intent(in) :: line
  character(len=maxlen),intent(out) :: stringinquote
  character :: cha
  integer :: i,q
  
  cha=' '
  stringinquote=''
  q=1
  do i=1,len(line)
    if (cha==' ') then
      if ((line(i:i)=='"').or.(line(i:i)=="'")) cha=line(i:i)
    else
      if (line(i:i)==cha) exit
      stringinquote(q:q)=line(i:i)
      q=q+1
    endif
  enddo
  return
  
  endsubroutine

! =================================================================== !












endmodule