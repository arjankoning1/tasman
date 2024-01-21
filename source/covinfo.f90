subroutine covinfo
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance information for ENDF-6 files
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
!   numleg         ! maximum number of Legendre coefficients
! Variables for writing TALYS input files
!   italys         ! TALYS run counter
! Variables for cross section covariances
!   Ecovindex      ! index for covariance energy on main energy grid
!   MTcov          ! MT number with covariance data
!   MTisocov       ! isomer of MT number with covariance data
!   Nchancovint    ! number of channels with covariance data
!   Nencov         ! number of covariance energies
! Variables for reading Legendre coefficients
!   Eleg           ! energy grid for Legendre coefficients
!   Nchanleg       ! total number of channels with Legendre coefficients
! Variables for reading cross sections
!   MT             ! MT number
!   MTiso          ! isomer of MT number
!   Nchanxs        ! total number of channels with cross sections
! Variables for reading residual production cross sections
!   Arp            ! mass number of residual product
!   Lrp            ! isomeric level of residual product
!   Nchanrp        ! total number of channels with residual production cross sections
!   Zrp            ! charge number of residual product
!
! *** Declaration of local data
!
  implicit none
  integer :: i              ! counter
!
! *************************** Output ***********************************
!
! Info file for TEFAL
!
  open (unit = 2, file = 'covariance.inf', status = 'replace')
  write(2, '(i4, " Number of random runs")') italys
  write(2, '(i4, " Number of energies for covariances")') Nencov
  do i = 1, Nencov
    write(2, '(i4)') Ecovindex(i)
  enddo
  write(2, '(i4, " Number of channels")') Nchanxs
  do i = 1, Nchanxs
    write(2, '(2i4)') MT(i), MTiso(i)
  enddo
  write(2, '(i4, " Number of channels with inter-channel covariances")') Nchancovint
  do i = 1, Nchancovint
    write(2, '(2i4)') MTcov(i), MTisocov(i)
  enddo
  write(2, '(2i4, " Number of channels with Legendre coeffients and number of Legendre coefficients")') Nchanleg, numleg
  do i = 1, Nchanleg
    write(2, '(es12.5)') Eleg(i)
  enddo
  write(2, '(i4, " Number of channels with residual production covariances")') Nchanrp
  do i = 1, Nchanrp
    write(2, '(3i4)') Zrp(i), Arp(i), Lrp(i)
  enddo
  close (2)
  return
end subroutine covinfo
! Copyright A.J. Koning 2021
