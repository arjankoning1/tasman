function dripvar(iz, ia)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Set parameter uncertainties for nuclides closer to dripline
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
!   sgl      ! single precision kind
! Constants
!   heavy    ! heaviest isotope
!   light    ! lightest isotope
!
! *** Declaration of local data
!
  implicit none
  integer   :: Aheavy  ! heaviest isotope per element
  integer   :: Alight  ! lightest isotope per element
  integer   :: Amid    ! middle of mass range
  integer   :: distcon ! distance from main isotope for which parameter is not extra  varied
  integer   :: ia      ! mass number of nucleus
  integer   :: iz      ! charge number of nucleus
  real(sgl) :: b       ! A value where uncertainty is halfway between min and max
  real(sgl) :: D       ! distance from constant uncertainty area
  real(sgl) :: dripvar ! function to set parameter uncertainties for nuclides  closer to dripline
  real(sgl) :: varmax  ! multiplication factor for parameter uncertainty: maximum  value
!
! ********************* Set parameter limits ***************************
!
! A smooth functional form is used to increase the parameter uncertainties for nuclides away from the line of stability
!
  dripvar = 1.
  varmax = 3.
  Alight = light(iz)
  Aheavy = heavy(iz)
  distcon = 3
  if (Aheavy - Alight < 2 * distcon) then
    Amid = (Alight + Aheavy) / 2
    Alight = Amid - (distcon - 1)
    Aheavy = Amid + distcon
  endif
  if (ia >= Alight .and. ia <= Aheavy) return
  if (ia > Aheavy) then
    D = real(ia - Aheavy)
    b = 5.
  else
    D = real(Alight - ia)
    b = 10.
  endif
  dripvar = 1. + (varmax - 1.) * D * D / (D * D + b * b)
  return
end function dripvar
! Copyright A.J. Koning 2021
