subroutine arbitrary(Plow, Nbin, hist, P)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Sampling from arbitrary distribution
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer, parameter :: dbl = selected_real_kind(15,307) ! double precision kind
  integer   :: j             ! counter
  integer   :: Nbin          ! number of bin
  real(sgl) :: cweight(Nbin) ! cumulative weight
  real(sgl) :: dP            ! parameter bin width
  real(sgl) :: hist(Nbin)    ! histogram value
  real(sgl) :: P             ! sampled parameter
  real(sgl) :: Plow(Nbin)    ! bottom of histogram
  real(sgl) :: R             ! random number
  real(sgl) :: sumweight     ! total cumulative weight
  real(sgl) :: w             ! cumulative weight
  real(dbl) :: mersenne      ! random number generator (Mersenne twister)
!
! Sampling from arbitrary distribution
!
  dP = (Plow(Nbin) - Plow(1)) / Nbin
  cweight = 0.
  sumweight = 0.
  do j = 1, Nbin
    sumweight = sumweight + hist(j)
    cweight(j) = sumweight
  enddo
  R = real(mersenne())
  w = cweight(1) + R * (sumweight - cweight(1))
  call locate_nr(cweight, Nbin, w, j)
  R = real(mersenne())
  P = Plow(j) + R * dP
  return
end subroutine arbitrary
! Copyright A.J. Koning 2021
