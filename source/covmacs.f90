subroutine covmacs
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for Maxwellian averaged cross sections
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
!
! Definition of single and double precision variables
!   sgl          ! single precision kind
! Variables for writing TALYS input files
!   italys       ! TALYS run counter
! Variables for cross section covariances
!   coveps       ! limit for covariance
! Variables for covariances
!   Sw0          ! sum of weights
!   Swp0         ! sum of weights
!   Sws0         ! sum of weights
! Variables for MACS
!   errmacs      ! Maxwellian averaged cross section uncertainty
!   macsav       ! average Maxwellian averaged cross section
!   RmacsD       ! diagonal of covariance matrix for Maxwellian averaged cross sections
! Variables for reading MACS
!   macstalys    ! Maxwellian averaged cross section from TALYS
!   ratmacs      ! Maxwellian rate
!   Tmacs        ! Maxwellian temperature
!
! *** Declaration of local data
!
  implicit none
  character(len=13) :: ofile     ! output file
  real(sgl)         :: err       ! error
  real(sgl)         :: term      ! help variable
  real(sgl)         :: term0     ! help variable
  real(sgl)         :: xs2       ! help variable
  real(sgl)         :: xsdifi    ! difference in cross section
  real(sgl)         :: xsi0      ! cross section of run 0
  real(sgl)         :: xsi1      ! cross section of random run
!
! Maxwellian averaged cross sections and uncertainties
!
  xsi0 = macstalys(0)
  xsi1 = min(macstalys(1), 10. * xsi0)
  xsdifi = xsi0 - xsi1
  term0 = xsdifi * xsdifi
  if (abs(term0) > coveps .and. Sws0 > 0.) then
    xs2 = xsi0 * xsi0
    if (xs2 > 0.) then
      term = term0 / xs2
      RmacsD = (RmacsD * Swp0 + Sw0 * term) / Sws0
      errmacs = sqrt(RmacsD)
    endif
  endif
  macsav = (macsav * Swp0 + Sw0 * xsi1) / Sws0
!
! Output
!
  ofile = 'astrorate.ave'
  open (unit = 3, file = ofile, status = 'replace')
  write(3, '("# astrorate.g Number of runs: ", i6)') italys
  write(3, '("#    T       Rate       MACS     uncertainty")')
  err = macsav * errmacs
  write(3, '(f8.4, 3es12.5)') Tmacs, ratmacs, macsav, err
  close (3)
  return
end subroutine covmacs
! Copyright A.J. Koning 2021
