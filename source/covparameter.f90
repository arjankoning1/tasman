subroutine covparameter
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for parameters
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
!   sgl            ! single precision kind
! All global variables
!   numenin        ! maximum number of incident energies
! Variables for writing TALYS input files
!   italys         ! TALYS run counter
! Variables for reading experimental data
!   Ein            ! incident energy
! Variables for weights
!   flagweight     ! flag to use weights for random samples
! Variables for GOF function
!   Ntalys         ! number of TALYS runs
! Variables for fitting limits
!   parfile        ! file with parameters
! Variables for reading cross sections
!   Nen            ! number of incident energies
!   Qval           ! Q - value
! Variables for parameter variation
!   Npar           ! number of parameters
!   pardelta       ! uncertainty of parameter
!   parinp         ! input parameter value
!   parsave        ! parameter value
!   partalys       ! parameter value
! Variables for parameter covariances
!   parav          ! average parameters
!   parcor         ! correlation of parameters
!   parcov         ! covariance matrix for parameters
!   pardif         ! difference of parameters
!   parstring      ! parameter string
! Variables for covariances
!   Sweight        ! weight for TALYS run
!   Sweightprev    ! weight for previous TALYS run
!   Sweightsum     ! sum of weights for TALYS run
!
! *** Declaration of local data
!
  implicit none
  character(len=25) :: pfile      ! parameter file
  integer           :: j          ! counter
  integer           :: k          ! Legendre order
  integer           :: l          ! counter
  real(sgl)         :: Ei         ! incident energy
  real(sgl)         :: Eprev      ! previous energy
  real(sgl)         :: err        ! error
  real(sgl)         :: finperc    ! relative parameter uncertainty on output
  real(sgl)         :: inperr     ! parameter uncertainty on input
  real(sgl)         :: inpperc    ! relative parameter uncertainty on input
  real(sgl)         :: pdifk      ! Rdifference in parameter
  real(sgl)         :: pdifl      ! Rdifference in parameter
  real(sgl)         :: ratio      ! ratio of uncertainties
  real(sgl)         :: term       ! help variable
!
! Output of accepted and rejected parameters
!
! parameter covariance matrix
!
  pdifk = 0.
  pdifl = 0.
  do k = 1, Npar
    if (parsave(0, k) /= 0.) pdifk = (partalys(1, k) - parsave(0, k)) / parsave(0, k)
    if (Sweightsum(0) > 0.) then
      do l = 1, Npar
        if (parsave(0, l) /= 0.) pdifl = (partalys(1, l) - parsave(0, l)) / parsave(0, l)
        term = pdifk * pdifl
        parcov(k, l) = (parcov(k, l) * Sweightprev(0) + Sweight(0) * term) / Sweightsum(0)
      enddo
      pardif(k) = (pardif(k) * Sweightprev(0) + Sweight(0) * pdifk) / Sweightsum(0)
    endif
    do j = 0, numenin
      if (Sweightsum(j) > 0.) parav(k, j) = (parav(k, j) * Sweightprev(j) + &
        Sweight(j) * partalys(1, k)) / Sweightsum(j)
    enddo
  enddo
  do k = 1, Npar
    do l = 1, Npar
      if (pardif(k) /= 0..and.pardif(l) /= 0.) parcor(k, l) = parcov(k, l) / (sqrt(abs(parcov(k, k))) * sqrt(abs(parcov(l, l))))
    enddo
  enddo
!
! Output of covariance matrices for parameters
!
  open (unit = 3, file = 'cov_parameter.ave', status = 'replace')
  write(3, '(6x, "parameter", 17x, "parameter", 11x, "      Pij     ", "  correlation")')
  do k = 1, Npar
    do l = 1, Npar
      write(3, '(2a25, 2es15.6)') parstring(k), parstring(l), parcov(k, l), parcor(k, l)
    enddo
  enddo
  close (3)
!
! Output of average parameters
!
  open (unit = 3, file = 'parameter.ave', status = 'replace')
  write(3, '("Number of accepted runs: ", i4)') italys
  write(3, '(6x, "parameter", 11x, "    average       uncertainty      %         input        uncertainty     %      Ratio")')
  do k = 1, Npar
    err = parav(k, 0) * sqrt(abs(parcov(k, k)))
    inperr = parsave(0, k) * pardelta(k)
    inpperc = 0.
    finperc = 0.
    ratio = 0.
    if (parsave(0, k) /= 0.) inpperc = 100. * inperr / parsave(0, k)
    if (parav(k, 0) /= 0.) finperc = 100. * err / parav(k, 0)
    if (finperc /= 0.) ratio = inpperc / finperc
    write(3, '(a25, 2(2es15.6, f9.3), f9.3)') parstring(k), parav(k, 0), err, finperc, parsave(0, k), inperr, inpperc, ratio
  enddo
  close (3)
!
! Output of energy-dependent average parameters
!
  if (flagweight) then
    do k = 1, Npar
      if (italys == Ntalys) then
        pfile = parfile(k)
      else
        pfile = 'upd.'//trim(parfile(k))
      endif
      open (unit = 8, file = pfile, status = 'replace')
      Eprev = 0.
      do j = 1, Nen(1)
        if (parinp(k) /= 0.) then
          ratio = parav(k, j) / parinp(k)
          if (pfile == 'rspincut' .or. pfile == 'm2constant') then
            Ei = Ein(j) + Qval
          else
            Ei = Ein(j)
          endif
          if (Ei > Eprev) write(8, '(2es12.5)') Ei, ratio
          Eprev = Ei
        endif
      enddo
      close(8)
    enddo
  endif
  return
end subroutine covparameter
! Copyright A.J. Koning 2021
