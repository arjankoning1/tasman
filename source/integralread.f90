subroutine integralread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read integral activation data from TALYS files
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
!   sgl          ! single precision kind
! All global variables
!   numpar       ! maximum number of parameters
! Variables for writing TALYS input files
!   italys       ! TALYS run counter
!   mode         ! TASMAN mode
! Variables for reading integral cross sections
!   Nxseff       ! number of effective cross sections
!   xseffchan    ! channel for effective cross section
!   xseffexp     ! experimental effective cross section
!   xseffflux    ! flux for effective cross section
!   xseffrat     ! C / E of effective cross section
!   xseffsave    ! effective cross section from TALYS
!   xsefftalys   ! effective cross section from TALYS
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist    ! logical to determine existence
  character(len=12) :: xsf       ! filename
  character(len=15) :: ch        ! character
  character(len=15) :: fl        ! flux
  integer           :: isamp     ! variable for 0th or random run
  integer           :: istat     ! logical for file access
  integer           :: k         ! Legendre order
  real(sgl)         :: rat       ! C/E ratio
  real(sgl)         :: xsc       ! interpolated cross section
  real(sgl)         :: xse       ! experimental cross section
!
! ************************ Read integral results ***********************
!
  xsf = 'integral.dat'
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  Nxseff = 0
  inquire (file = xsf, exist = lexist)
  if (lexist) then
    open (unit = 2, file = xsf, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(xsf, istat)
    read(2, '(/)', iostat = istat)
    if (istat /= 0) call read_error(xsf, istat)
    k = 0
    do
      read(2, '(2a15, 2e12.5, f15.5)', iostat = istat) ch, fl, xsc, xse, rat
      if (istat /= -1) exit
      k = k + 1
      xseffchan(k) = ch
      xseffflux(k) = fl
      xsefftalys(isamp, k) = xsc
      xseffexp(k) = xse
      xseffrat(isamp, k) = rat
      if (mode == 2 .and. italys >= 0 .and. italys <= numpar) xseffsave(italys, k) = xsc
    enddo
    Nxseff = k
    close (2)
  endif
  return
end subroutine integralread
! Copyright A.J. Koning 2021
