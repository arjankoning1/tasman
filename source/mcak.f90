subroutine mcak(ipar, iR, Plow, Nbin, hist, k, kmax, P, flagkd, flagequi)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Sampling
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
  logical            :: flagequi                ! flag to use equidistant excitation bins instead of  logarithmic bins
  logical            :: flagkd                  ! flag to use KD function
  integer, parameter :: numpar=500              ! maximum number of parameters
  integer, parameter :: numhis=100              ! number of histogram bins
  integer            :: i                       ! counter
  integer            :: ipar                    ! counter for parameter
  integer            :: iR                      ! parameter to read from existing distribution
  integer            :: j                       ! counter
  integer            :: k                       ! Legendre order
  integer            :: khalf                   ! half of cumulative histogram
  integer            :: kmax                    ! maximum of distribution
  integer            :: Nbin                    ! number of bin
  integer            :: npow                    ! power of distribution
  integer            :: Pdone(numpar, 0:numhis) ! flag to determine whether parameter has been visited
  real(sgl)          :: crit                    ! help variable
  real(sgl)          :: f                       ! E-Ef
  real(sgl)          :: fac                     ! factor
  real(sgl)          :: hist(Nbin)              ! histogram value
  real(sgl)          :: P                       ! sampled parameter
  real(sgl)          :: Plow(Nbin)              ! bottom of histogram
  real(sgl)          :: Prange                  ! parameter range
  real(sgl)          :: R                       ! random number
  real(sgl)          :: RR                      ! random number
  real(dbl)          :: mersenne                ! random number generator (Mersenne twister)
  save Pdone
!
! Read from existing parameter distribution if requested
!
  if (iR == 1) then
    call arbitrary(Plow, Nbin, hist, P)
    return
  endif
!
! Criterion for sampling from uniform or current distribution
!
  if (k == 1 .or. mod(k, Nbin) == 0) then
    do i = 0, Nbin
      Pdone(ipar, i) = 0
    enddo
    if (k == 1) Pdone(ipar, Nbin / 2) = 1
  endif
!
! KD function
! flagkd : flag to use KD function
!
  if (flagkd) then
    npow = 4
    khalf = (kmax + 1) / 2
    fac = 1. + 0.5 **npow
    f = fac * real(k) **npow / (real(k) **npow + real(khalf) **npow)
    crit = f
!
! Binary
!
  else
    if (k <= kmax) then
      crit = 0.
    else
      crit = 1.
    endif
  endif
  R = real(mersenne())
!
! Sample from uniform distribution
!
  if (R > crit) then
    Prange = Plow(Nbin) - Plow(1)
    i = 0
    do
      i = i + 1
      RR = real(mersenne())
      P = Plow(1) + RR * Prange
      call locate_nr(Plow, Nbin, P, j)
      if (.not. (flagequi .and. Pdone(ipar, j) == 1 .and. i <= 100 * Nbin)) exit
    enddo
    Pdone(ipar, j) = 1
  else
!
! Sample from current distribution
!
    call arbitrary(Plow, Nbin, hist, P)
  endif
  return
end subroutine mcak
! Copyright A.J. Koning 2021
