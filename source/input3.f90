subroutine input3
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for third set of variables
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
! Definition of single and double precision variables
!   sgl           ! single precision kind
! All global variables
!   nummt         ! maximum number of MT numbers
!  Variables for deviation
!   adevi         ! parameter for energy dependence of deviation
!   bdevi         ! parameter for energy dependence of deviation
!   cdevi         ! parameter for energy dependence of deviation
!   ddevi         ! parameter for energy dependence of deviation
!   devi          ! parameter for energy dependence of deviation
!   Ecent         ! parameter for energy dependence of deviation
! Variables for reading TASMAN input lines
!   inline        ! input line
!   nlines        ! number of input lines
! Variables for processing input
!   Atarget       ! mass number of target nucleus
!   changeline    ! line number at which TASMAN specific input starts
!   Ztarget       ! charge number of target nucleus
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: key         ! keyword
  character(len=132) :: value       ! value or string
  character(len=132) :: word(40)    ! words on input line
  character(len=132) :: line        ! input line
  integer            :: i           ! counter
  integer            :: imt         ! MT counter
  integer            :: istat       ! logical for file access
  real(sgl)          :: dripvar     ! function to set parameter uncertainties for nuclides
  real(sgl)          :: varunc      ! multiplication factor for nuclides closer to dripline
!
! ******* Default acceptable deviations for correlated covariances *****
!
! dripvar : function to set parameter uncertainties for nuclides closer to dripline
!
  devi = 1.
  Ecent = 6.
  adevi = 1.
  bdevi = 1.
  cdevi = 1.
  ddevi = 12.
  devi(1) = 0.08
  devi(2) = 0.12
  devi(3) = 0.12
  devi(4) = 0.12
  devi(16) = 0.24
  devi(17) = 0.40
  devi(18) = 0.50
  devi(22) = 0.40
  devi(28) = 0.40
  devi(37) = 0.60
  do imt = 51, 91
    devi(imt) = 0.40
  enddo
  devi(102) = 0.40
  devi(103) = 0.34
  devi(104) = 0.60
  devi(105) = 0.80
  devi(106) = 0.80
  devi(107) = 0.45
  varunc = dripvar(Ztarget, Atarget)
  devi = varunc * devi
  adevi(1) = 0.6
  bdevi(1) = 0.3
  cdevi(1) = 0.
  adevi(2) = 0.6
  bdevi(2) = 0.3
  cdevi(2) = 0.
  adevi(3) = 0.6
  bdevi(3) = 0.3
  cdevi(3) = 0.
  Ecent(4) = 5.
  ddevi(16) = 15.
  do imt = 51, 91
    bdevi(imt) = 1.
    Ecent(imt) = 3.
  enddo
  adevi(102) = 0.6
  bdevi(102) = 0.3
  cdevi(102) = 1.
  ddevi(102) = 20.
  Ecent(102) = 20.
  cdevi(104) = 2.
  cdevi(105) = 2.
  cdevi(106) = 2.
!
! ********************** Read input parameters *************************
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
  do i = changeline + 1, nlines
    line = inline(i)
    key = ''
    call getkeywords(line, word)
    key = word(1)
    value = word(2)
    if (key == '#dev') then
      read(value, * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imt, 1, nummt)
      read(word(3), * , iostat = istat) devi(imt)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#ecent') then
      read(value, * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imt, 1, nummt)
      read(word(3), * , iostat = istat) Ecent(imt)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#adev') then
      read(value, * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imt, 1, nummt)
      read(word(3), * , iostat = istat) adevi(imt)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#bdev') then
      read(value, * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imt, 1, nummt)
      read(word(3), * , iostat = istat) bdevi(imt)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#cdev') then
      read(value, * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imt, 1, nummt)
      read(word(3), * , iostat = istat) cdevi(imt)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#ddev') then
      read(value, * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, imt, 1, nummt)
      read(word(3), * , iostat = istat) ddevi(imt)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
  enddo
  return
end subroutine input3
! Copyright A.J. Koning 2021
