subroutine checkvalue
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Check for errors in values
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
  use A1_error_handling_mod
!
! All global variables
!   numlib          ! maximum number of data libraries
!   nummass         ! maximum number of masses
!   nummt           ! maximum number of MT numbers
!   numsets         ! maximum number of experimental data sets
!   numtalys        ! maximum number of TALYS runs
!   numZ            ! maximum number of elements
! Variables for automatic parameter variation
!   Adeep           ! maximum depth of A
!   maxbar          ! maximal fission barrier
!   minbar          ! minimal fission barrier
!   Zdeep           ! maximum depth of Z
! Variables for deviation
!   adevi           ! parameter for energy dependence of deviation
!   bdevi           ! parameter for energy dependence of deviation
!   cdevi           ! parameter for energy dependence of deviation
!   ddevi           ! parameter for energy dependence of deviation
!   devi            ! parameter for energy dependence of deviation
!   Ecent           ! parameter for energy dependence of deviation
! Variables for cross section covariances
!   average         ! cross sections 1: average, 2: first run 3: optimum
! Variables for reading cross sections
!   Efracmax        ! starting energy for maximum cross section
!   fracmax         ! fraction of maximum cross section to be included
! Variables for uncertainty
!   cwidth          ! scaling constant for parameter variation
!   fiscor0         ! adjustment parameter for fission uncertainties
! Variables for reading TALYS output
!   flaggamma       ! flag for covariances or optimization of discrete gamma ray transitions
!   flagprod        ! flag for covariances or optimization of particle production cross sections
!   flagresidual    ! flag for covariances or optimization of residual production cross sections
! Variables for reading gamma production cross sections
!   maxgam          ! maximum number of discrete levels with gamma cross section
! Variables for reading experimental data
!   Esearch1        ! start energy of search
!   Esearch2        ! end energy of search
!   expclass        ! quality class of experimental data (1, 2 or 3)
!   maxexpsets      ! number of experimental data sets before weight reduction
!   mtexpweight     ! weight for experimental data set and channel
!   mtlibweight     ! weight for library and channel
! Variables for writing TALYS input files
!   Ehigh           ! energy cutoff of high energy runs
!   flagomponly     ! flag for covariances or optimization for basic optical model observables
!   mode            ! TASMAN mode
!   Nhigh           ! number of high energy runs
! Variables for GOF function
!   chi2max         ! maximal chi2 value per point taken into account
!   Fmax            ! maximal F value per point taken into account
!   ichi2           ! parameter to set goodness - of - fit estimator
!   Nburn           ! number of burn - in runs
!   Ntalys          ! number of TALYS runs
!   outsearch       ! level of output of search (1 (minimal) - 5 (all))
!   searchmode      ! search mode
! Variables for TMC
!   offset          ! offset for numbering or random files (TMC only)
!   tafislib        ! library name for adoption of nubar values (none[default], endfb8.0, jeff3.0, jendl5.0, brc09)
!   taneslib        ! library name for adoption of FNS (none[default], endfb8.0, jeff3.0, jendl5.0, brc09)
!   tareslib        ! library name for adoption of resonance parameters (default or endfb8.0)
!   tmcoffset       ! offset for starting creation of ENDF - 6 files (TMC only)
! Variables for processing input
!   dsmooth         ! parameter for energy smoothing
!   errlim          ! lower limit for experimental error
!   Liso            ! number of the isomeric state
!   seed            ! seed for random number generator
!   Tweight         ! weight of TALYS vs experiment
! Variables for weights
!   gofmode         ! gof mode
!   weightpower     ! power for weight of random run
! Variables for experimental data files
!   aamax           ! maximal A value for which exp. data are included in the search
!   aamin           ! minimal A value for which exp. data are included in the search
!   zzmax           ! maximal Z value for which exp. data are included in the search
!   zzmin           ! minimal Z value for which exp. data are included in the search
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   range_real_error    ! Test if real variable is out of range
!
! *** Declaration of local data
!
  implicit none
  character(len=16) :: library    ! data library
  integer           :: i          ! counter
  integer           :: imt        ! MT counter
  integer           :: j          ! counter
!
! All parameters need to fall within certain ranges.
! These ranges are specified in this subroutine and in the manual.
!
! ******************* Check for wrong input variables ******************
!
! 1. Check of values for main keywords.
!
  call range_integer_error('#mode', mode, 1, 4)
  call range_integer_error('#ntalys', Ntalys, 0, numtalys, default = -1)
  call range_integer_error('#nburn', Nburn, 0, numtalys, default = -1)
  call range_real_error('#Ehigh', Ehigh, 1., 1000.)
  call range_integer_error('#tmcoffset', tmcoffset, 0, numtalys)
  call range_integer_error('#nhigh', Nhigh, 0, numtalys)
  call range_integer_error('#liso', Liso, 0, 10)
  call range_integer_error('#maxgam', maxgam, 1, 4)
  call range_integer_error('#minbar', minbar, 0, 3)
  call range_integer_error('#maxbar', maxbar, 0, 3)
  call range_integer_error('#maxbar', maxbar, minbar, 3)
  call range_integer_error('#chi2', ichi2, 1, 3)
  call range_integer_error('#gofmode', gofmode, 1, 6)
  call range_integer_error('#searchmode', searchmode, 1, 4)
  call range_integer_error('#maxexpsets', maxexpsets, 1, 1000000)
  call range_integer_error('#average', average, 1, 3)
  call range_real_error('#weightpower', weightpower, 0.1, 10.)
  call range_integer_error('#seed', seed, 1, 1000000000)
  call range_real_error('#cwidth', cwidth, 0., 20.)
  call range_real_error('#fiscor', fiscor0, 0.1, 10.)
  call range_integer_error('#offset', offset, 0, 10000)
  call range_integer_error('#zmin', zzmin, 1, numZ)
  call range_integer_error('#zmax', zzmax, 1, numZ)
  call range_integer_error('#amin', aamin, 0, nummass)
  call range_integer_error('#amax', aamax, 0, nummass)
  call range_integer_error('#zmax', zzmax, zzmin, nummass)
  call range_integer_error('#amax', aamax, aamin, nummass, default = 0)
  call range_integer_error('#expclass', expclass, 1, 3)
  call range_integer_error('#outsearch', outsearch, 1, 5)
  call range_integer_error('#zdeep', Zdeep, 0, numZ, default = -1)
  call range_integer_error('#adeep', Adeep, 0, numZ, default = -1)
!
! Check unlogical input combinations
!
  if (flagomponly .and. flagresidual) then
    write(*, '(" TASMAN-error: if #omponly y then #residual n")')
    stop
  endif
  if (flagomponly .and. flagprod) then
    write(*, '(" TASMAN-error: if #omponly y then #production n")')
    stop
  endif
  if (flagomponly .and. flaggamma) then
    write(*, '(" TASMAN-error: if #omponly y then #gamma n")')
    stop
  endif
!
! Check deviations and weights
!
  do imt = 1, nummt
    call range_real_error('#dev', devi(imt), 0., 1.e38, index1 = imt, name1 = 'MT')
    call range_real_error('#adev', adevi(imt), 0., 1.e38, index1 = imt, name1 = 'MT')
    call range_real_error('#bdev', bdevi(imt), 0., 1.e38, index1 = imt, name1 = 'MT')
    call range_real_error('#cdev', cdevi(imt), 0., 1.e38, index1 = imt, name1 = 'MT')
    call range_real_error('#ddev', ddevi(imt), 0., 1.e38, index1 = imt, name1 = 'MT')
    call range_real_error('#Ecent', Ecent(imt), 0., 1.e38, index1 = imt, name1 = 'MT')
    call range_real_error('#Fmax', Fmax(imt), 1., 1.e38, index1 = imt, name1 = 'MT')
    call range_real_error('#chi2max', chi2max(imt), 1., 1.e38, index1 = imt, name1 = 'MT')
    call range_real_error('#errlim', errlim(imt), 0., 1.e10, index1 = imt, name1 = 'MT')
    call range_real_error('#efracmax', Efracmax(imt), 0., 10., index1 = imt, name1 = 'MT')
    call range_real_error('#fracmax', fracmax(imt), 0., 1., index1 = imt, name1 = 'MT')
    call range_real_error('#esearch1', Esearch1(imt), 0., 1000., index1 = imt, name1 = 'MT')
    call range_real_error('#esearch2', Esearch2(imt), 0., 1000., index1 = imt, name1 = 'MT')
    call range_real_error('#esearch2', Esearch1(imt), Esearch1(imt), 1000., index1 = imt, name1 = 'MT')
    call range_real_error('#tweight', Tweight(imt), 0., 1.e10, index1 = imt, name1 = 'MT')
    call range_real_error('#dsmooth', dsmooth(imt), 0.001, 10., index1 = imt, name1 = 'MT')
    do i = 1, numlib
      call range_real_error('#weight', mtlibweight(imt, i), 0., 1.e10, index1 = imt, name1 = 'MT')
    enddo
    do i = 1, numsets
      call range_real_error('#weight', mtexpweight(imt, i), 0., 1.e10, index1 = imt, name1 = 'MT')
    enddo
  enddo
!
! Check for correct name of libraries
!
  do j = 1, 3
    if (j == 1) library = tareslib(1:16)
    if (j == 2) library = tafislib(1:16)
    if (j == 3) library = taneslib(1:16)
    if (j == 1) then
      if (library(1:7) == 'default') cycle
      if (library(1:3) == 'urr') cycle
    else
      if (library(1:4) == 'none') cycle
    endif
    if (library(1:7) == 'jeff3.3') cycle
    if (library(1:8) == 'endfb8.0') cycle
    if (library(1:8) == 'jendl5.0') cycle
    if (library(1:5) == 'brc09') cycle
    write(*, '(" TASMAN-error: Wrong library name: ", a16)') library
    stop
  enddo
  return
end subroutine checkvalue
! Copyright A.J. Koning 2021
