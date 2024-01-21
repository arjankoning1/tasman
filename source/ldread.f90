subroutine ldread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read level density GOF data from TALYS files
!
! Author    : Arjan Koning
!
! 2023-10-01: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl ! single precision kind
! Variables for automatic parameter variation
!   ACN            ! A of compound nucleus
!   ZCN            ! Z of compound nucleus
! Variables for writing TALYS input
!   italys          ! TALYS run counter
! Variables for reading TALYS output
!   Nchanall        ! total number of channels
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist    ! logical to determine existence
  character(len=12) :: xsf       ! filename
  character(len=132) :: line
  character(len=132) :: key
  integer           :: isamp     ! variable for 0th or random run
  integer           :: istat     ! logical for file access
  integer           :: keyix
  real(sgl)         :: Dexp
  real(sgl)         :: dDexp
  real(sgl)         :: Dth
  real(sgl)         :: Frms
!
! Read level densities
!
  Nchanall = 1
  Nenexp(1, 1) = 2
  Nsets(1) = 1
  MT(1) = 1
  MTexp(1) = 1
  NMTexp = 1
  xsf = 'ld000000.tot'
  write(xsf(3:5), '(i3.3)') ZCN
  write(xsf(6:8), '(i3.3)') ACN
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  inquire (file = xsf, exist = lexist)
  if (lexist) then
    open (unit = 2, file = xsf, status = 'old')
    read(2, '()', iostat = istat)
    if (istat /= 0) call read_error(xsf, istat)
    do
      read(2,'(a)', iostat = istat) line
      if (istat == -1) exit
      key='experimental D0 [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Dexp
      if (istat /= 0) call read_error(xsf, istat)
      key='experimental D0 unc. [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) dDexp
      if (istat /= 0) call read_error(xsf, istat)
      key='theoretical D0 [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Dth
      if (istat /= 0) call read_error(xsf, istat)
      key='Frms per level'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Frms
      if (istat /= 0) call read_error(xsf, istat)
    enddo
    close (2)
    xsexp(1, 1, 1) = Dexp
    dxsexp(1, 1, 1) = dDexp
    xsth(1, 1, 1) = Dth
    xsexp(1, 1, 2) = 1.
    dxsexp(1, 1, 2) = 0.1
    xsth(1, 1, 2) = Frms
  endif
  return
end subroutine ldread
! Copyright A.J. Koning 2023
