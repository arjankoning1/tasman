subroutine gof
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Goodness-of-fit
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
!   sgl           ! single precision kind
! Variables for GOF function
!   isearch       ! number of trial run
! Variables for weights
!   flagweight    ! flag to use weights for random samples
! Variables for parameter variation
!   Npar          ! number of parameters
!   par           ! parameter value
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: Fdum            ! dummy value
!
! ****************** Determine goodness-of-fit and weights *************
!
  call fcn(Npar, par, Fdum)
  if (flagweight) then
    call weights
    if (isearch > 0) call histoweights
    call histogof
  endif
  return
end subroutine gof
! Copyright A.J. Koning 2021
