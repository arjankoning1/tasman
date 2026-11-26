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
  real(dbl)         :: err       ! error
  real(dbl)         :: term      ! help variable
  real(dbl)         :: xsdifi    ! difference in cross section
  real(dbl)         :: xsi0      ! cross section of run 0
  real(dbl)         :: xsi1      ! cross section of random run
!
! Maxwellian averaged cross sections and uncertainties
!
  xsi0 = macstalys(0)
  xsi1 = min(real(macstalys(1), dbl), 10. * xsi0)
  xsdifi = xsi0 - xsi1
  if (Sws0 > 0.) then
    term = 0.
    if (xsi0 /= 0.) term = (xsdifi / xsi0) * (xsdifi / xsi0)
    RmacsD = (RmacsD * Swp0 + Sw0 * term) / Sws0
    errmacs = sqrt(RmacsD)
    macsav = (macsav * Swp0 + Sw0 * xsi1) / Sws0
  endif
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
