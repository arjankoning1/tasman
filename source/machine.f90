subroutine machine
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Machine dependent statements
!
! Author    : Arjan Koning
!
! 2023-12-29: Original code
! 2026-08-25: Runtime definition of TASMAN directory and user
! 2026-08-27: Resolve external executables through PATH
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
!
! Variables for path names
!   librariespath ! directory containing files to be read
!   exforpath     ! directory containing files to be read
!   binpath       ! optional prefix for external executables
!   tasmanpath    ! directory containing files to be read
!   psfpath       ! directory containing files to be read
!
! *** Declaration of local data
!
  implicit none
  logical             :: lexist       ! logical to determine existence
  character(len=1024) :: code_dir     ! TASMAN directory
  character(len=1024) :: base_dir     ! directory containing TASMAN and related packages
  character(len=1024) :: tasman_dir   ! TASMAN directory runtime defined
  character(len=1024) :: tasman_user  ! TASMAN user runtime defined
  integer             :: envstat
  integer             :: i            ! counter
  integer             :: n
  integer             :: values(8)
!
! ************************ Set directories *****************************
!
! The preferred option is to set an environment variable TASMAN_DIR,
! e.g. put in your ~/.profile or ~/.zshrc file.
!
! export TASMAN_DIR=/path/to/tasman
!
! If TASMAN_DIR is not set, get_environment_variable will simply return
! an empty string.
!
  call get_environment_variable('TASMAN_DIR', tasman_dir, length=n, status=envstat)
  if (envstat == 0 .and. n > 0) then
    code_dir = trim(tasman_dir)
  else
!
! If for some reason the above does not work, the code directory can be
! changed here manually.
!
    code_dir = '/path/to/tasman/'
  endif
!
! Remove a trailing slash, if present, and determine the parent directory.
! TASMAN expects libraries, exfortables and PSF to be sibling directories
! of the TASMAN directory.
!
  i = len_trim(code_dir)
  if (i > 1) then
    if (code_dir(i:i) == '/') code_dir = code_dir(:i - 1)
  endif
  i = scan(trim(code_dir), '/', back=.true.)
  if (i > 0) then
    base_dir = code_dir(:i)
  else
    base_dir = './'
  endif
!
  librariespath = trim(base_dir)//'libraries/'
  exforpath = trim(base_dir)//'exfortables/'
!
! External executables such as TALYS, TEFAL, TARES, TAFIS and TANES are
! resolved through the user's PATH.  input2.f90 appends the executable
! name to binpath, so an empty prefix lets the operating system find them.
!
  binpath = ''
  tasmanpath = trim(code_dir)//'/'
  psfpath = trim(base_dir)//'PSF/Photo/'
!
! Test to check accessibility of TASMAN data files.
!
  inquire (file=trim(tasmanpath)//'misc/score.tab', exist=lexist)
  if (.not. lexist) then
    write(*, '(a)') 'TASMAN error: misc database not found.'
    write(*, '(2a)') 'Expected file: ', trim(tasmanpath)//'misc/score.tab'
    write(*, '(a)') 'Set the TASMAN_DIR environment variable:'
    write(*, '(a)') '  export TASMAN_DIR=/path/to/tasman'
    write(*, '(a)') 'Alternatively, edit code_dir in source/machine.f90'
    write(*, '(a)') 'and rebuild TASMAN.'
    error stop 77
  endif
!
! Set date.
!
  call date_and_time(VALUES=values)
  year=values(1)
  month=values(2)
  day=values(3)
  date='xxxx-xx-xx'
  write(date(1:4),'(i4.4)') year
  write(date(6:7),'(i2.2)') month
  write(date(9:10),'(i2.2)') day
!
! Set user name for output files.
! The preferred option is to set an environment variable TASMAN_USER,
! e.g. put in your ~/.profile or ~/.zshrc file.
!
! export TASMAN_USER="Your Name"
!
  call get_environment_variable('TASMAN_USER', tasman_user, length=n, status=envstat)
  if (envstat == 0 .and. n > 0) then
    user = trim(tasman_user)
  else
    user = 'Unknown User'
  endif
  return
end subroutine machine
! Copyright A.J. Koning 2026
