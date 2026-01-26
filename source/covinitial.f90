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
  Vmt = 0.
  Rmt = 0.
  Cmt = 0.
  Rrp = 0.
  Rgam = 0.
  Rprod = 0.
  RmtD = 0.
  errmt = 0.
  xsav = 0.
  SE = 0.
  RrpD = 0.
  errrp = 0.
  rpav = 0.
  RgamD = 0.
  errgam = 0.
  gamav = 0.
  RprodD = 0.
  errprod = 0.
  prodav = 0.
  Sall = 0.
  parcov = 0.
  parcor = 0.
  Vcov = 0.
  Rcov = 0.
  Ccov = 0.
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
  S = 0.
  Senum = 0.
  Sdenom = 0.
  Schan = 0.
  Pearson = 0.
  Pearson_enum = 0.
  Pearson_denom_par = 0.
  Pearson_denom_xs = 0.
  xsdevav = 0.
  xsavC = 0.
  errmtC = 0.
  pardif = 0.
  parav = 0.
  angav = 0.
  errang = 0.
  Rang = 0.
  errleg = 0.
  legav = 0.
  Rleg = 0.
  return
end subroutine covinitial
! Copyright A.J. Koning 2021
