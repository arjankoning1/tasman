subroutine macsread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read Maxwellian averaged cross sections from TALYS files
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
! Variables for reading MACS
!   macssave     ! average Maxwellian averaged cross section
!   macstalys    ! Maxwellian averaged cross section from TALYS
!   ratmacs      ! Maxwellian rate
!   Tmacs        ! Maxwellian temperature
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist    ! logical to determine existence
  character(len=11) :: xsf       ! filename
  integer           :: isamp     ! variable for 0th or random run
  integer           :: istat     ! logical for file access
!
! ************************ Read integral results ***********************
!
  xsf = 'macs.g'
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  inquire (file = xsf, exist = lexist)
  if (lexist) then
    open (unit = 2, file = xsf, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(xsf, istat)
    do
      read(2, '(/8x, 3e12.5)', iostat = istat) thmacs, expmacs, dexpmacs
      if (istat /= 0) exit
      macstalys(isamp) = thmacs
    enddo
    close (2)
    Nchanall = Nchanall + 1
    Nsets(Nchanall) = 1
    MT(Nchanall) = 102
    MTexp(Nchanall) = 102
    NMTexp = 1
    Nenexp(Nchanall,1) = 1
    xsexp(Nchanall,1,1) = expmacs
    dxsexp(Nchanall,1,1) = dexpmacs
    xsth(Nchanall,1,1) = thmacs
  endif
  return
end subroutine macsread
! Copyright A.J. Koning 2021
