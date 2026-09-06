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
  if (allocated(xsav)) deallocate(xsav)
  if (allocated(xsavC)) deallocate(xsavC)
  if (allocated(Vmt)) deallocate(Vmt)
  if (allocated(Rmt)) deallocate(Rmt)
  if (allocated(RmtD)) deallocate(RmtD)
  if (allocated(Rcov)) deallocate(Rcov)
  if (allocated(Vcov)) deallocate(Vcov)
  if (allocated(rpsave)) deallocate(rpsave)
  if (allocated(Rrp)) deallocate(Rrp)
  if (allocated(RrpD)) deallocate(RrpD)
  if (allocated(rpav)) deallocate(rpav)
  if (allocated(prodsave)) deallocate(prodsave)
  if (allocated(Rprod)) deallocate(Rprod)
  if (allocated(RprodD)) deallocate(RprodD)
  if (allocated(prodav)) deallocate(prodav)
  if (allocated(gamsave)) deallocate(gamsave)
  if (allocated(gamav)) deallocate(gamav)
  if (allocated(Rgam)) deallocate(Rgam)
  if (allocated(RgamD)) deallocate(RgamD)
  if (allocated(Rsp)) deallocate(Rsp)
  if (allocated(spav)) deallocate(spav)
  if (allocated(angav)) deallocate(angav)
  if (allocated(Rang)) deallocate(Rang)
  if (allocated(legav)) deallocate(legav)
  if (allocated(Rleg)) deallocate(Rleg)
! if (allocated(parsave)) deallocate(parsave)
  if (allocated(S)) deallocate(S)
  if (allocated(Sdenom)) deallocate(Sdenom)
  if (allocated(Senum)) deallocate(Senum)
  if (allocated(Pearson)) deallocate(Pearson)
  if (allocated(Pearson_enum)) deallocate(Pearson_enum)
  if (allocated(Pearson_denom_par)) deallocate(Pearson_denom_par)
  if (allocated(Pearson_denom_xs)) deallocate(Pearson_denom_xs)
  if (allocated(Gchannelsave)) deallocate(Gchannelsave)
  if (allocated(xseffsave)) deallocate(xseffsave)
!
  return
end subroutine deallocate_arrays
