subroutine magnet
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Unify theoretical and experimental data
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
!   sgl             ! single precision kind
! All global variables
!   numenin         ! maximum number of incident energies
!   numsets         ! maximum number of experimental data sets
! Variables for writing TALYS input files
!   italys          ! TALYS run counter
! Variables for reading TALYS output
!   flagtal         ! flag to optimize to TALYS
!   Nchanall        ! total number of channels
! Variables for processing input
!   dsmooth         ! parameter for energy smoothing
!   Tweight         ! weight of TALYS vs experiment
! Variables for parameter variation
!   flagsample      ! flag for sampling experimental data
! Variables for GOF function
!   Eindex          ! index for experimental data
! Variables for reading cross sections
!   E               ! incident energy
!   MT              ! MT number
!   Nen             ! number of incident energies
!   xseval          ! evaluated cross section
!   xstalys         ! cross section from TALYS
!   xsth            ! theoretical cross section
! Variables for reading experimental data
!   Eexp            ! incident energy
!   Eweight         ! incident energy
!   Nenexp          ! number of experimental data points
!   Nsets           ! number of experimental data sets
!   xsexp           ! experimental cross section
!   xsexpsamp       ! sampled experimental cross section
! Variables for weights
!   Esamp           ! energy of sampled cross section
!   flagchmagnet    ! flag to unify theoretical and experimental data
!   Nxssamp         ! number of sampled cross sections
!   xssamp          ! sampled random cross section
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                                  ! counter
  integer   :: iexp                               ! index
  integer   :: isamp                              ! variable for 0th or random run
  integer   :: iset0                              ! start of counter for data sets
  integer   :: j                                  ! counter
  integer   :: k                                  ! Legendre order
  integer   :: m                                  ! counter
  integer   :: n                                  ! counter
  integer   :: Ndxs(numenin)                      ! number of experimental data points per energy bin
  real(sgl) :: d                                  ! parameter for energy smoothing
  real(sgl) :: dxs                                ! uncertainty of experimental cross section
  real(sgl) :: dxsE(numenin)                      ! average experimental cross section per energy bin
  real(sgl) :: dxsEs(numenin, numsets)            ! average experimental cross section per data set and energy bin
  real(sgl) :: dxseval                            ! uncertainty of evaluated cross section
  real(sgl) :: dxssum(numenin)                    ! normalization factor
  real(sgl) :: ee(0:numenin)                      ! incident energy
  real(sgl) :: Eexp0                              ! incident energy
  real(sgl) :: et                                 ! energy gain at scission with respect to ground state
  real(sgl) :: expo                               ! exponent
  real(sgl) :: fedis                              ! energy dependent deviation
  real(sgl) :: sum                                ! help variable
  real(sgl) :: sumw                               ! sum over weights
  real(sgl) :: Tw                                 ! weight of TALYS vs experiment
  real(sgl) :: w                                  ! cumulative weight
  real(sgl) :: ws                                 ! surface component of the imaginary potential
  real(sgl) :: xse                                ! experimental cross section
  real(sgl) :: xst                                ! theoretical cross section
  real(sgl) :: xst2                               ! theoretical cross section
!
! ******************* Final local adjustment ***************************
!
  if (italys == 0) then
    isamp = 0
  else
    isamp = 1
  endif
  if (flagtal) then
    iset0 = 2
  else
    iset0 = 1
  endif
  do i = 1, Nchanall
    if (Nsets(i) == 0) cycle
    if (MT(i) == 0) cycle
    flagchmagnet(i) = .true.
    ee(0) = 0.
!
! Set energies and energy bins for averaging
!
    do j = 1, Nen(i)
      ee(j) = E(i, j)
    enddo
!
! Loop over experimental data sets
!
    iexp = 0
    do n = iset0, Nsets(i)
      do j = 1, Nen(i)
        Ndxs(j) = 0
        dxssum(j) = 0.
        dxsEs(j, n) = 0.
      enddo
!
! Loop over experimental data points in set. Sample experimental data points and calculate deviation from model per energy bin.
!
      do k = 1, Nenexp(i, n)
        if (isamp == 0 .or. .not. flagsample) then
          xse = xsexp(i, n, k)
        else
          xse = xsexpsamp(i, n, k)
        endif
        Eexp0 = Eexp(i, n, k)
        iexp = iexp + 1
        Esamp(i, iexp) = Eexp0
        xssamp(i, iexp) = xse
        xst = xsth(i, n, k)
        m = Eindex(i, n, k)
        dxs = xse - xst
        if (m == 0) cycle
        Ndxs(m) = Ndxs(m) + 1
        dxssum(m) = dxssum(m) + dxs
      enddo
      do j = 1, Nen(i)
        if (Ndxs(j) > 0) dxsEs(j, n) = dxssum(j) / Ndxs(j)
      enddo
    enddo
    Nxssamp(i) = iexp
!
! Average deviation per energy point
!
    do j = 1, Nen(i)
      ws = 0
      dxsE(j) = 0.
      sum = 0.
      do n = iset0, Nsets(i)
        if (abs(dxsEs(j, n)) > 1.e-3) then
          ws = ws + Eweight(i, n)
          sum = sum + Eweight(i, n) * dxsEs(j, n)
        endif
      enddo
      if (ws > 0.) dxsE(j) = sum / ws
    enddo
!
! Calculate weight of experimental data points on all grid points
!
    Tw = Tweight(MT(i))
    d = dsmooth(MT(i))
    do j = 1, Nen(i)
      et = ee(j)
      dxseval = 0.
      sumw = 0.
      xst = xstalys(isamp, i, j)
      if (xst > 0.) then
        do k = 1, Nen(i)
          if (dxsE(k) == 0.) cycle
          xst2 = xstalys(isamp, i, k)
          if (xst2 == 0.) cycle
          Eexp0 = ee(k)
          expo = abs(et - Eexp0) / d
          if (expo < 80.) then
            fedis = exp( - expo)
          else
            fedis = 0.
          endif
          w = fedis * min(xst / xst2, xst2 / xst)
          dxs = w * w * dxsE(k)
          dxseval = dxseval + dxs
          sumw = sumw + w
        enddo
        if (sumw > 0.) dxseval = dxseval / sumw
      endif
      xse = xst + dxseval
      xseval(isamp, i, j) = Tw * xst + (1. - Tw) * xse
    enddo
  enddo
  return
end subroutine magnet
! Copyright A.J. Koning 2021
