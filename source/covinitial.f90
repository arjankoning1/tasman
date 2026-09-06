subroutine covinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization for covariance and sensitivity calculations
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
! ************************** Initialization ****************************
!
  Sweightprev = 0.
  Sweight = 1.
  Sweightsum = 0.
  errmt = 0.
  SE = 0.
  Sall = 0.
  parcov = 0.
  parcor = 0.
  RactD = 0.
  erract = 0.
  actav = 0.
  RNisoD = 0.
  errNiso = 0.
  Nisoav = 0.
  RyieldD = 0.
  erryield = 0.
  yieldav = 0.
  RNisorelD = 0.
  errNisorel = 0.
  Nisorelav = 0.
  Schan = 0.
  xsdevav = 0.
  errmtC = 0.
  pardif = 0.
  parav = 0.
  return
end subroutine covinitial
! Copyright A.J. Koning 2021
