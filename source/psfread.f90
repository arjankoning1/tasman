subroutine psfread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read photon strength functions from TALYS files
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
!   sgl ! single precision kind
! Variables for automatic parameter variation
!   ACN            ! A of compound nucleus
!   ZCN            ! Z of compound nucleus
! Variables for writing TALYS input
!   italys          ! TALYS run counter
! Variables for reading TALYS output
!   Nchanall        ! total number of channels
! Variables for reading cross sections
!   E           ! incident energy
!   Efracmax    ! starting energy for maximum cross section
!   MT          ! MT number
!   Nen         ! number of incident energies
!   Nen0        ! number of incident energies
!   xseval      ! evaluated cross section
!   xsmax       ! maximum cross section per channel
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
  integer           :: k         ! counter
  integer           :: keyix
  real(sgl)         :: ee        ! energy
  real(sgl)         :: psf       ! PSF
  real(sgl)         :: xsnew     ! new cross section
  real(sgl)         :: xstop     ! maximum cross section
!
! Read PSF
!
  Nchanall = 1
  xsf = 'psf000000.E1'
  write(xsf(4:6), '(i3.3)') ZCN
  write(xsf(7:9), '(i3.3)') ACN
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  inquire (file = xsf, exist = lexist)
  if (lexist) then
    open (unit = 2, file = xsf, status = 'old')
    E(1, 0) = 0.
    read(2, '()', iostat = istat)
    if (istat /= 0) call read_error(xsf, istat)
    do
      read(2,'(a)', iostat = istat) line
      if (istat == -1) exit
      key='entries'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Nen(1)
        if (istat /= 0) call read_error(xsf, istat)
        read(2,'(/)')
        do k = 1, Nen(1)
          read(2, * , iostat = istat) ee, psf
          if (istat == -1) exit
          if (istat /= 0) call read_error(xsf, istat)
          E(1, k) = ee
          xseval(isamp, 1, k) = psf
        enddo
        exit
      endif
    enddo
    Nen0(1) = Nen(1)
    close (2)
  endif
  if (italys <= 0) then
    xstop = 0.
    do k = 1, Nen(1)
      if (E(1, k) <= Efracmax(MT(1))) cycle
      xsnew = xseval(0, 1, k)
      if (xsnew > xstop) xstop = xsnew
    enddo
    xsmax(1) = xstop
  endif
  return
end subroutine psfread
! Copyright A.J. Koning 2021
