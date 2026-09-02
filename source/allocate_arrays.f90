subroutine allocate_arrays
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Allocate dynamic arrays
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_tasman_mod
  implicit none
!
  if (flagcross) then
    allocate(xssave(0:Ntalys,numchanxs,0:numenS))
    xssave = 0.
  endif

  if (flagresidual) then
    allocate(rpsave(0:Ntalys,numchanrp,0:numenS))
    rpsave = 0.
  endif

  if (flagprod) then
    allocate(prodsave(0:Ntalys,numchanprod,0:numenS))
    prodsave = 0.
  endif

  if (flaggamma) then
    allocate(gamsave(0:Ntalys,numchangam,0:numenS))
    gamsave = 0.
  endif
  allocate(S(Npar,numchanxs,numenS))
  S = 0.
  allocate(Sdenom(Npar,numchanxs,numenS))
  Sdenom = 0.
  allocate(Senum(Npar,numchanxs,numenS))
  Senum = 0.
  allocate(Pearson(Npar,numchanxs,numenS))
  Pearson = 0.
  allocate(Pearson_enum(Npar,numchanxs,numenS))
  Pearson_enum = 0.
  allocate(Pearson_denom_par(Npar,numchanxs,numenS))
  Pearson_denom_par = 0.
  allocate(Pearson_denom_xs(Npar,numchanxs,numenS))
  Pearson_denom_xs = 0.
  return
end subroutine allocate_arrays
