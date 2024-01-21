subroutine input1
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for first set of variables
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
!   numZ          ! maximum number of elements
! Variables for reading experimental data
!   Ein           ! incident energy
! Constants
!   nuc           ! symbol of nucleus
! Variables for reading TASMAN input lines
!   inline        ! input line
!   nlines        ! number of input lines
! Variables for processing input
!   Atarget       ! mass number of target nucleus
!   changeline    ! line number at which TASMAN specific input starts
!   energyfile    ! file with incident energies
!   Ltarget       ! excited level of target
!   numE          ! number of incident energies
!   ptype0        ! type of incident particle
!   Starget       ! symbol of target nucleus
!   Ztarget       ! charge number of target nucleus
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch          ! character
  character(len=132) :: key         ! keyword
  character(len=132) :: value       ! value or string
  character(len=132) :: word(40)    ! words on input line
  character(len=132) :: line        ! input line
  integer            :: i           ! counter
  integer            :: istat       ! logical for file access
  integer            :: iz          ! charge number of nucleus
!
! ************ Read input for first set of input variables *************
!
! An input line starting with '#change' has a special role:
! All input lines after it are either TASMAN input flags or TALYS keywords that will be varied.
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
  Ltarget = 0
  numE = 0
  do i = 1, nlines
    line = inline(i)
    key = ''
    call getkeywords(line, word)
    key = word(1)
    value = word(2)
    ch = word(2)(1:1)
    if (key == '#change') then
      changeline = i
      return
    endif
!
! TASMAN needs to know Z, A and the incident energies
!
    if (key == 'element') then
      if (ch >= '0' .and. ch <= '9') then
        read(value, * , iostat = istat) Ztarget
        if (istat /= 0) call read_error(line, istat)
        call range_integer_error(key, Ztarget, 1, numZ)
        cycle
      else
        read(value, '(a2)', iostat = istat) Starget
        if (istat /= 0) call read_error(line, istat)
        Starget(1:1) = char(ichar(Starget(1:1)) - 32)
        do iz = 1, numZ
          if (nuc(iz) == Starget) then
            Ztarget = iz
            cycle
          endif
        enddo
        cycle
      endif
    endif
    if (key == 'projectile') then
      ptype0 = ch
      cycle
    endif
    if (key == 'mass') then
      read(value, * , iostat = istat) Atarget
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'ltarget') then
      read(value, * , iostat = istat) Ltarget
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == 'energy') then
      if ((ch >= '0' .and. ch <= '9') .or. ch == '.') then
        read(value, * , iostat = istat) Ein(1)
        if (istat /= 0) call read_error(line, istat)
        numE = 1
        cycle
      else
        read(value, '(a73)', iostat = istat) energyfile
        if (istat /= 0) call read_error(line, istat)
      endif
    endif
  enddo
  return
end subroutine input1
! Copyright A.J. Koning 2021
