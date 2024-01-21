subroutine locate_nr(xx, n, x, j)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Find value in ordered table
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tasman_mod
!
! Definition of single and double precision variables
!   sgl              ! single precision kind
!
! *** Declaration of local data
!
  integer   :: j                ! counter
  integer   :: n                ! counter
  real(sgl) :: x                ! initial parameter vector
  real(sgl) :: xx(n)            ! x value
  integer   :: jl               ! lower value
  integer   :: jm               ! middle value
  integer   :: ju               ! higher value
  jl = 0
  ju = n + 1
10    if(ju - jl > 1)then
    jm = (ju + jl) / 2
    if((xx(n) > xx(1)).eqv.(x > xx(jm)))then
      jl = jm
    else
      ju = jm
    endif
  goto 10
  endif
  j = jl
  return
end subroutine locate_nr
! Copyright A.J. Koning 2021
