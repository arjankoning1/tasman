subroutine readinput

!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
  use A1_error_handling_mod
!
! All global variables
!   numlines    ! maximum number of input lines
! Variables for reading TASMAN input lines
!   inline      ! input line
!   nlines      ! number of input lines
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=132)  :: line  ! input line
  integer             :: i     ! counter
  integer             :: istat ! logical for file access
!
! ************************** User Input ********************************
!
! We read the complete input file first as a set of character strings.
! The actual keywords will be read from these later on.
! For natural elements, the input file only needs to be read once.
!
  i = 1
  do
    read(*, '(a132)', iostat = istat) inline(i)
    if (istat ==  -1) exit
    if (istat /= 0) call read_error(inline(i), istat)
    i = i + 1
    call range_integer_error('inline', i, 1, numlines)
  enddo
  nlines = i - 1
!
! ************** Convert uppercase to lowercase characters *************
!
! For easy handling of all the input parameters, the whole input is converted to lowercase characters,
! with the exception of filenames or other character strings.
!
! convert: subroutine to convert input line from upper case to lowercase
!
  do i = 1, nlines
    line = inline(i)
    call convert(line)
    inline(i) = line
  enddo
  return
end subroutine readinput
! Copyright A.J. Koning 2021
