subroutine talread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read TALYS data
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
!   sgl            ! single precision kind
! All global variables
!   numenexp       ! maximum number of energies forexperimental data set
! Variables for reading TALYS output
!   Nchanall       ! total number of channels
! Variables for processing input
!   auto           ! flag to in - or exclude channel
! Variables for reading cross sections
!   E              ! incident energy
!   E1mb           ! energy of 1 mb value
!   MT             ! MT number
!   MTiso          ! isomer of MT number
!   Nen            ! number of incident energies
!   xstalys        ! cross section from TALYS
! Variables for deviation
!   adevi          ! parameter for energy dependence of deviation
!   bdevi          ! parameter for energy dependence of deviation
!   cdevi          ! parameter for energy dependence of deviation
!   ddevi          ! parameter for energy dependence of deviation
!   devi           ! parameter for energy dependence of deviation
!   Ecent          ! parameter for energy dependence of deviation
! Variables for reading experimental data
!   auth           ! author
!   dxsexp         ! uncertainty of experimental cross section
!   Eexp           ! incident energy
!   Esearch1       ! start energy of search
!   Esearch2       ! end energy of search
!   Eweight        ! incident energy
!   flagedep       ! flag for energy dependence of uncertainties
!   isoexp         ! isomer of experimental data
!   mttalweight    ! weight per TALYS channel
!   Nenexp         ! number of experimental data points
!   Nsets          ! number of experimental data sets
!   xsexp          ! experimental cross section
!   expyear        ! year
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                 ! counter
  integer   :: imt               ! MT counter
  integer   :: j                 ! counter
  integer   :: l                 ! counter
  integer   :: Nb                ! begin of energy loop
  integer   :: Ne                ! begin of energy loop
  real(sgl) :: Dhigh             ! high value of parameter function
  real(sgl) :: Dlow              ! low value of parameter function
  real(sgl) :: Ec                ! help variable for energy
  real(sgl) :: ee0               ! incident energy
  real(sgl) :: Eth               ! help variable
  real(sgl) :: Fdev              ! function value for deviation
  real(sgl) :: isofac            ! factor for isomer
!
! ************************** Read TALYS data ***************************
!
! Copy TALYS data per channel
!
  do i = 1, Nchanall
    imt = MT(i)
    if (imt == 0) cycle
    if (auto(imt) == 0) cycle
    isoexp(i, 1) = MTiso(i)
    if (MTiso(i) ==  -1) then
      isofac = 1.
    else
      isofac = 1.7
    endif
    j = 0
    Nb = max(1, Nen(i) - numenexp + 1)
    Ne = min(Nen(i), numenexp)
    do l = Nb, Ne
      ee0 = E(i, l)
      if (ee0 < Esearch1(imt) .or. ee0 > Esearch2(imt)) cycle
      j = j + 1
      Eexp(i, 1, j) = ee0
      xsexp(i, 1, j) = xstalys(0, i, l)
      if (flagedep) then
        Eth = max(ee0 - E1mb(i), 0.)
        Dlow = max(adevi(imt) - devi(imt), 0.) * exp( - bdevi(imt) * Eth)
        Ec = Eth - Ecent(imt)
        if (Ec > 0) then
          Dhigh = max(cdevi(imt) - devi(imt), 0.) * Ec * Ec / (Ec * Ec + ddevi(imt) * ddevi(imt))
        else
          Dhigh = 0.
        endif
        Fdev = devi(imt) + Dlow + Dhigh
      else
        Fdev = devi(imt)
      endif
      dxsexp(i, 1, j) = isofac * Fdev * xstalys(0, i, l)
    enddo
    auth(i, 1) = 'TALYS-global'
    expyear(i, 1) = 0
    Nenexp(i, 1) = j
    Nsets(i) = 1
    Eweight(i, 1) = mttalweight(imt)
  enddo
  return
end subroutine talread
! Copyright A.J. Koning 2021
