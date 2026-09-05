subroutine deallocate_arrays
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Deallocate dynamic arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tasman_mod
  implicit none
!
  if (allocated(xssave)) deallocate(xssave)
  if (allocated(rpsave)) deallocate(rpsave)
  if (allocated(prodsave)) deallocate(prodsave)
  if (allocated(gamsave)) deallocate(gamsave)
! if (allocated(parsave)) deallocate(parsave)
  if (allocated(S)) deallocate(S)
  if (allocated(Sdenom)) deallocate(Sdenom)
  if (allocated(Senum)) deallocate(Senum)
  if (allocated(Pearson)) deallocate(Pearson)
  if (allocated(Pearson_enum)) deallocate(Pearson_enum)
  if (allocated(Pearson_denom_par)) deallocate(Pearson_denom_par)
  if (allocated(Pearson_denom_xs)) deallocate(Pearson_denom_xs)
!
  return
end subroutine deallocate_arrays
