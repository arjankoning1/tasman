subroutine machine
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Machine dependent statements
!
! Author    : Arjan Koning
!
! 2023-12-29: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
!
! Variables for path names
!   librariespath ! directory containing files to be read
!   exforpath     ! directory containing files to be read
!   binpath       ! directory containing files to be read
!   tasmanpath    ! directory containing files to be read
!   psfpath       ! directory containing files to be read
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: codedir   ! code directory
  character(len=132) :: basedir   ! base directory
  character(len=132) :: homedir   ! home directory
  character(len=132) :: envvar    ! environment variable
  integer            :: values(8)
  integer            :: ix
  integer            :: istat     ! I/O status
!
! ************************ Set directories *****************************
!
! Try to get TASMAN directory from environment variable first
  call get_environment_variable('TASMAN_DIR', envvar, status=istat)
  if (istat == 0 .and. len_trim(envvar) > 0) then
    codedir = trim(envvar)
    if (codedir(len_trim(codedir):len_trim(codedir)) /= '/') then
      codedir = trim(codedir) // '/'
    endif
  else
    ! Try to get from HOME directory
    call get_environment_variable('HOME', homedir, status=istat)
    if (istat == 0 .and. len_trim(homedir) > 0) then
      codedir = trim(homedir) // '/tasman/'
    else
      ! Fallback to current working directory
      codedir = './tasman/'
    endif
  endif
  ix = index(codedir,'/tasman/')
  basedir = codedir(1:ix)
  librariespath = trim(basedir) // 'libraries/'
  exforpath = trim(basedir) // 'exfortables/'
  binpath = trim(basedir) // 'bin/'
  tasmanpath = codedir
  psfpath = trim(basedir) // 'PSF/Photo/'
  call date_and_time(VALUES=values)
  year=values(1)
  month=values(2)
  day=values(3)
  date='xxxx-xx-xx'
  write(date(1:4),'(i4.4)') year
  write(date(6:7),'(i2.2)') month
  write(date(9:10),'(i2.2)') day
  user = 'Arjan Koning'
  return
end subroutine machine
! Copyright A.J. Koning 2021
