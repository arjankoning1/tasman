subroutine sensinput
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Input parameters for sensitivity matrix
!
! Author    : Arjan Koning
!
! 2025-12-25: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
  use A1_error_handling_mod
!
! Variables for parameter variation
!   Npar           ! number of parameters
!   par            ! parameter value
!   parsave        ! parameter value
!   paradjust      ! local adjustment parameters
!   parinp         ! input parameter value
!   pardelta       ! uncertainty of parameter
!   partype        ! type of parameter variation
!
! *** Declaration of local data
!
  implicit none
  integer   :: i        ! counter
  real(dbl) :: mersenne ! random number generator (Mersenne twister)
  real(sgl), save  :: parorg(numpar) 
!
! ******* One-parameter variation for linear sensitivity matrix ********
!
  if (italys == 0) then
    if (flagmorris) then
      do i = 1, Npar
        parorg(i) = parsave(italys,i)
      enddo
    endif
    return
  endif
  do i = 1, Npar
    if (flagmorris) then
      par(i) = parorg(i)
    else
      par(i) = parinp(i)
    endif
    parsave(italys, i) = par(i)
  enddo
  if (flagmorris) then
    par(italys) = parlow(italys) + real(mersenne()) * (parhigh(italys) - parlow(italys))
  else
    if (partype(italys) == 'shift ') then
      par(italys) = par(italys) + pardelta(italys)
    else
      par(italys) = par(italys) * (1. + pardelta(italys))
    endif
  endif
  parsave(italys, italys) = par(italys)
  return
end subroutine sensinput
! Copyright A.J. Koning 2025
