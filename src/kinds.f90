module kinds
  implicit none
  save
  !...kind definitions
  integer ,parameter :: dp    = kind(1.0d0)
  integer ,parameter :: dpc   = kind((1.0d0,1.0d0))
  integer ,parameter :: i4b   = selected_int_kind(9)  !abs(n) < 10**9 (kind)
  integer ,parameter :: i8b   = selected_int_kind(18)
  integer ,parameter :: sgl   = selected_real_kind(6,30)
  integer ,parameter :: dpreal= selected_real_kind(14,200)
  
  private
  public :: dp, dpc, i4b, i8b, sgl, dpreal, print_kind_info
  
  contains
  
  !!   Print information about the used data types.
  subroutine print_kind_info (stdout)
    implicit none 
    integer ,intent(in):: stdout
    
    write(stdout,"(/,1X,A67)") repeat("=",67)
    write(stdout,'(1X,"=",T5,A,T68,"=")') 'DATA TYPE INFORMATION:'
    write(stdout,"(1X,A67)") repeat("=",67)

    write(stdout,'(/,T2,A,T78,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E15.8))') &
          'REAL: Data type name:', 'DP', '      Kind value:', kind(0.0_DP), &
          '      Precision:', precision(0.0_DP), &
          '      Smallest nonnegligible quantity relative to 1:', &
          epsilon(0.0_DP), '      Smallest positive number:', tiny(0.0_DP), &
          '      Largest representable number:', huge(0.0_DP)
    write( stdout,'(/,T2,A,T78,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E15.8))') &
          '      Data type name:', 'sgl', '      Kind value:', kind(0.0_sgl), &
          '      Precision:', precision(0.0_sgl), &
          '      Smallest nonnegligible quantity relative to 1:', &
          epsilon(0.0_sgl), '      Smallest positive number:', tiny(0.0_sgl), &
          '      Largest representable number:', huge(0.0_sgl)
    write( stdout,'(/,T2,A,T72,A,4(/,T2,A,T61,I20))') &
          'INTEGER: Data type name:', '(default)', '         Kind value:', &
          kind(0), '         Bit size:', bit_size(0), &
          '         Largest representable number:', huge(0)
    write( stdout,'(/,T2,A,T72,A,/,T2,A,T75,I6,/)') 'LOGICAL: Data type name:', &
          '(default)', '         Kind value:', kind(.TRUE.)
    write( stdout,'(/,T2,A,T72,A,/,T2,A,T75,I6,/)') &
          'CHARACTER: Data type name:', '(default)', '           Kind value:', &
          kind('C')
    write(stdout,*) repeat("=",67)
    
  END SUBROUTINE print_kind_info

END MODULE