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
    allocate(xsav(numchanxs,numenin))
    xsav = 0.
    allocate(xsavC(numchanxs,numencov))
    xsavC = 0.
    allocate(Vmt(numchanxs,numencov,numencov))
    Vmt = 0.
    allocate(Rmt(numchanxs,numencov,numencov))
    Rmt = 0.
    allocate(RmtD(numchanxs,numenin))
    RmtD = 0.
    allocate(Rcov(numchancov,numencov,numchancov,numencov))
    Rcov =0.
    allocate(Vcov(numchancov,numencov,numchancov,numencov))
    Vcov = 0.
  endif

  if (flagresidual) then
    allocate(rpsave(0:Ntalys,numchanrp,0:numenS))
    rpsave = 0.
    allocate(rpav(numchanrp,numenin))
    rpav = 0.
    allocate(Rrp(numchanrp,numenin,numencov))
    Rrp = 0.
    allocate(RrpD(numchanrp,numenin))
    RRpD = 0.
  endif

  if (flagspectra) then
    allocate(Rsp(numchansp,0:numen2, 0:numen2))
    Rsp = 0.
    allocate(spav(numchansp,0:numen2))
    spav = 0.
  endif

  if (flagprod) then
    allocate(prodsave(0:Ntalys,numchanprod,0:numenS))
    prodsave = 0.
    allocate(prodav(numchanprod,numenin))
    prodav = 0.
    allocate(Rprod(numchanprod,numencov,numencov))
    Rprod = 0.
    allocate(RprodD(numchanprod,numenin))
    RprodD = 0.
  endif

  if (flagangle) then
    allocate(angav(numchanang,0:numang))
    angav = 0.
    allocate(Rang(numchanang,0:numang,0:numang))
    Rang = 0.
  endif

  if (flagleg) then
    allocate(legav(numchanang,0:numleg))                 
    legav = 0.
    allocate(Rleg(numchanang,0:numleg,numchanang,0:numleg))
    Rleg = 0.
  endif

  if (flagintegral) then
    allocate(xseffsave(0:Ntalys,numchanxs))
    xseffsave = 0.
  endif

  if (flaggamma) then
    allocate(gamsave(0:Ntalys,numchangam,0:numenS))
    gamsave = 0.
    allocate(gamav(numchangam,numenin))
    gamav = 0.
    allocate(Rgam(numchangam,numencov,numencov))
    Rgam = 0.
    allocate(RgamD(numchangam,numenin))
    RgamD = 0.
  endif

  if (flagsens .or. flagreadsens .or. mode == 2) then
    allocate(S(Npar,numchanxs,numenS))
    S = 0.
  endif
  if (flagcross .and. flagsens .and. .not. flagreadsens) then
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
  endif
  if (flagexp .or. flaglib .or. flagtal) then
    allocate(Gchannelsave(0:Ntalys,numchanxs))
    Gchannelsave = 0.
  endif
  return
end subroutine allocate_arrays
