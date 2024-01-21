subroutine covariance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for TALYS results and parameters
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
! All global variables
!   numenin         ! maximum number of incident energies
! Variables for reading TALYS output
!   flagcross       ! flag for covariances or optimization of cross sections
!   flaggamma       ! flag for covariances or optimization of discrete gamma ray transitions
!   flagintegral    ! flag for covariances or optimization of integral activation data
!   flagmacs        ! flag for covariances or optimization of Maxwellian averaged cross sections
!   flagprod        ! flag for covariances or optimization of particle production cross sections
!   flagproduct     ! flag for covariances or optimization of isotope production
!   flagresidual    ! flag for covariances or optimization of residual production cross sections
! Variables for writing TALYS input files
!   flagangle       ! flag for covariances or optimization of angular distributions
!   flagleg         ! flag for covariances or optimization of Legendre coefficients
!   flagspectra     ! flag for covariances or optimization of emission spectra
!   flagtmc         ! flag for Total Monte Carlo (ENDF - 6 data file at each run)
! Variables for covariances
!   Sw0             ! sum of weights
!   Sweight         ! weight for TALYS run
!   Sweightprev     ! weight for previous TALYS run
!   Sweightsum      ! sum of weights for TALYS run
!   Swp0            ! sum of weights
!   Sws0            ! sum of weights
!
! *** Declaration of local data
!
  implicit none
  integer :: i              ! counter
!
! Calculation of covariance matrices
!
! covparameter : subroutine for covariance matrices for parameters
! covinfo      : subroutine with covariance information for ENDF-6 files
! covcross     : subroutine for covariance matrices for cross sections production cross sections
! covresidual  : subroutine for covariance matrices for residual production cross sections
! covgamma     : subroutine for covariance matrices for gamma production cross sections
! covprod      : subroutine for covariance matrices for particle production cross sections
! covspectra   : subroutine for covariance matrices for emission spectra distributions
! covangle     : subroutine for covariance matrices for angular distributions
! covlegendre  : subroutine for covariance matrices for Legendre coefficients
! covintegral  : subroutine for Covariance matrices for integral activation data
! covproduction: subroutine for covariance matrices for isotope production
! covmacs      : subroutine for Covariance matrices for Maxwellian averaged cross sections
!
  do i = 0, numenin
    Sweightprev(i) = Sweightsum(i)
    Sweightsum(i) = Sweightsum(i) + Sweight(i)
  enddo
  Sws0 = Sweightsum(0)
  Swp0 = Sweightprev(0)
  Sw0 = Sweight(0)
  call covparameter
  call covinfo
  if (flagcross) call covcross
  if (flagresidual) call covresidual
  if (flaggamma) call covgamma
  if (flagprod) call covprod
  if (flagspectra) call covspectra
  if (flagangle) call covangle
  if (flagleg) call covlegendre
  if (flagintegral) call covintegral
  if (flagproduct) call covproduction
  if (flagmacs) call covmacs
!
! Jump to system for possible random ENDF-6 file generation
!
! tmc      : subroutine for creation of random ENDF-6 files
!
  if (flagtmc) call tmc
  return
end subroutine covariance
! Copyright A.J. Koning 2021
