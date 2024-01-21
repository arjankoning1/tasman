subroutine covintegral
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for integral activation data
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
! Variables for integral data covariances
!   erreff       ! effective cross section uncertainty
!   Reff         ! covariance matrix for effective cross sections
!   ReffD        ! diagonal of covariance matrix for effective cross sections
!   xseffav      ! average effective cross section
! Variables for cross section covariances
!   coveps       ! limit for covariance
! Variables for covariances
!   Sw0          ! sum of weights
!   Swp0         ! sum of weights
!   Sws0         ! sum of weights
! Variables for reading integral cross sections
!   Nxseff       ! number of effective cross sections
!   xseffchan    ! channel for effective cross section
!   xseffexp     ! experimental effective cross section
!   xseffflux    ! flux for effective cross section
!   xseffrat     ! C / E of effective cross section
!   xsefftalys   ! effective cross section from TALYS
!
! *** Declaration of local data
!
  implicit none
  character(len=12) :: ofile     ! output file
  integer           :: i         ! counter
  integer           :: j         ! counter
  real(sgl)         :: err       ! error
  real(sgl)         :: errrat    ! ratio of error
  real(sgl)         :: term      ! help variable
  real(sgl)         :: term0     ! help variable
  real(sgl)         :: xs2       ! help variable
  real(sgl)         :: xsdifi    ! difference in cross section
  real(sgl)         :: xsdifk    ! difference in cross section
  real(sgl)         :: xsi0      ! cross section of run 0
  real(sgl)         :: xsi1      ! cross section of random run
  real(sgl)         :: xsk0      ! cross section of run 0
  real(sgl)         :: xsk1      ! cross section of random run
!
! Average effective cross sections and covariances
!
  do i = 1, Nxseff
    xsi0 = xsefftalys(0, i)
    xsi1 = min(xsefftalys(1, i), 10. * xsi0)
    xsdifi = xsi0 - xsi1
    term0 = xsdifi * xsdifi
    if (abs(term0) > coveps .and. Sws0 > 0.) then
      xs2 = xsi0 * xsi0
      if (xs2 > 0.) then
        term = term0 / xs2
        ReffD(i) = (ReffD(i) * Swp0 + Sw0 * term) / Sws0
        erreff(i) = sqrt(ReffD(i))
      endif
    endif
    if (Sws0 > 0.) xseffav(i) = (xseffav(i) * Swp0 + Sw0 * xsi1) / Sws0
    do j = 1, Nxseff
      xsk0 = xsefftalys(0, j)
      xsk1 = min(xsefftalys(1, j), 10. * xsk0)
      xsdifk = xsk0 - xsk1
      term0 = xsdifi * xsdifk
      if (abs(term0) > coveps .and. Sws0 > 0.) then
        xs2 = xsi0 * xsk0
        if (xs2 > 0.) then
          term = term0 / xs2
          Reff(i, j) = (Reff(i, j) * Swp0 + Sw0 * term) / Sws0
        endif
      endif
    enddo
  enddo
!
! Output
!
  ofile = 'integral.ave'
  open (unit = 3, file = ofile, status = 'replace')
  write(3, '("# integral.dat Number of runs: ", i6)') italys
  write(3, '("# # Channel      Flux        Eff. c.s. (b)", " uncertainty Exp. c.s. (b) Ratio    uncertainty")')
  do j = 1, Nxseff
    err = xseffav(j) * erreff(j)
    errrat = xseffrat(0, j) * erreff(j)
    if (err /= 0.) err = min(err, xseffav(j) - 1.e-13)
    write(3, '(2a15, 5e12.5)') xseffchan(j), xseffflux(j), xseffav(j), err, xseffexp(j), xseffrat(0, j), errrat
  enddo
  close (3)
  return
end subroutine covintegral
! Copyright A.J. Koning 2021
