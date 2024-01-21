function gasdev()
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Gaussian random numbers
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
  use A0_tasman_mod
!
! Definition of single and double precision variables
!   dbl                 ! double precision kind
!   sgl                 ! single precision kind
!
! *** Declaration of local data
!
  real(sgl) :: gasdev              ! Gaussian random number
  integer   :: iset                ! help variable
  real(sgl) :: fac                 ! factor
  real(sgl) :: gset                ! factor
  real(sgl) :: rsq                 ! square of random numbers
  real(sgl) :: v1                  ! random number
  real(sgl) :: v2                  ! random number
  real(dbl) :: mersenne            ! random number generator (Mersenne twister)
  SAVE iset, gset
!
! ****************** Gaussian random numbers ***************************
!
  DATA iset / 0 /
  if (iset == 0) then
    v1 = 2. * real(mersenne()) - 1.
    v2 = 2. * real(mersenne()) - 1.
    do
      rsq = v1 **2 + v2 **2
      if(rsq >= 1..or.rsq == 0.)then
        v1 = 2. * real(mersenne()) - 1.
        v2 = 2. * real(mersenne()) - 1.
      else
        exit
      endif
    enddo
    fac = sqrt( - 2. * log(rsq) / rsq)
    gset = v1 * fac
    gasdev = v2 * fac
    iset = 1
  else
    gasdev = gset
    iset = 0
  endif
  return
end function gasdev
! Copyright A.J. Koning 2021
