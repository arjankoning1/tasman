subroutine gamgamread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read average radiatve widths from TALYS files
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
!   mode         ! TASMAN mode
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist     ! logical to determine existence
  character(len=3)  :: Zstring    ! help variable
  character(len=3)  :: Astring    ! help variable
  character(len=12) :: psffile    ! filename
  character(len=132):: line
  character(len=132):: key
  integer           :: istat      ! logical for file access
  integer           :: keyix
  real              :: expgamgam  ! experimental gamma gamma
  real              :: dexpgamgam ! d experimental gamma gamma
  real              :: thgamgam   ! theoretical gamma gamma
!
! ************************ Read integral results ***********************
!
  Zstring = '   '
  Astring = '   '
  write(Zstring,'(i3.3)') ZCN
  write(Astring,'(i3.3)') ACN
  psffile = 'psf'//Zstring//Astring//'.E1'
  inquire (file = psffile, exist = lexist)
  if (lexist) then
    open (unit = 2, file = psffile, status = 'old', iostat = istat)
    if (istat /= 0) call read_error(psffile, istat)
    do
      read(2,'(a)',iostat = istat) line
      if (istat == -1) exit
      key='experimental Gamma_gamma [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(line(keyix+len_trim(key)+2:80),*, iostat = istat) expgamgam
        if (istat /= 0) call read_error(psffile, istat)
      endif
      key='experimental Gamma_gamma unc. [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(line(keyix+len_trim(key)+2:80),*, iostat = istat) dexpgamgam
        if (istat /= 0) call read_error(psffile, istat)
      endif
      key='theoretical Gamma_gamma [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(line(keyix+len_trim(key)+2:80),*, iostat = istat) thgamgam
        if (istat /= 0) call read_error(psffile, istat)
      endif
    enddo
    close (2) 
    Nchanall = Nchanall + 1
    Nsets(Nchanall) = 1
    MT(Nchanall) = 102
    MTexp(Nchanall) = 102
    NMTexp = 1
    Nenexp(Nchanall,1) = 1
    xsexp(Nchanall,1,1) = expgamgam
    dxsexp(Nchanall,1,1) = dexpgamgam
    xsth(Nchanall,1,1) = thgamgam
  endif
  return
end subroutine gamgamread
! Copyright A.J. Koning 2022
