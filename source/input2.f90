subroutine input2
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for second set of variables
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
! Definition of single and double precision variables
!   sgl               ! single precision kind
! All global variables
!   numlib            ! maximum number of data libraries
!   nummt             ! maximum number of MT numbers
!   numsets           ! maximum number of experimental data sets
!   numZ              ! maximum number of elements
! Variables for reading TASMAN input lines
!   inline            ! input line
!   nlines            ! number of input lines
! Variables for path names
!   binpath        ! directory containing files to be read
! Variables for automatic parameter variation
!   flagparvary    ! flag to automatically vary parameters
!   flagallpar     ! flag to vary all nuclear model parameters
!   partvary       ! particle for varying parameters
!   keyvary        ! keyword to vary
!   Adeep          ! maximum depth of A
!   Zdeep          ! maximum depth of Z
!   maxbar         ! maximum fission barrier
!   minbar         ! minimum fission barrier
!   Nkeyvary       ! number of keywords to vary
!   Npartvary      ! number of particles to vary keywords
!   NZAskip        ! number of Z,A pairs to skip for parameter variation
!   NZAvary        ! number of Z,A pairs for parameter variation
!   Askip          ! A value to skip
!   Zskip          ! Z value to skip
!   Avary          ! A value to vary
!   Zvary          ! Z value to vary
!   sigvary        ! sigma for parameter variation
! Variables for GOF function
!   chi2max           ! maximal chi2 value per point taken into account
!   Fmax              ! maximal F value per point taken into account
!   ichi2             ! parameter to set goodness - of - fit estimator
!   Nburn             ! number of burn - in runs
!   Ntalys            ! number of TALYS runs
!   outsearch         ! level of output of search (1 (minimal) - 5 (all))
!   searchmode        ! search mode
! Variables for uncertainty
!   cwidth            ! scaling constant for parameter variation
!   fiscor0           ! adjustment parameter for fission uncertainties
!   flagsave          ! flag to save all input and output files of each Talys run
! Variables for reading experimental data
!   Esearch1          ! start energy of search
!   Esearch1all       ! global start energy of search
!   Esearch2          ! end energy of search
!   Esearch2all       ! global end energy of search
!   expclass          ! quality class of experimental data (1, 2 or 3)
!   expexccount       ! counter for excluded experimental data sets
!   expexcfile        ! file with experimental data set to exclude for search
!   expinccount       ! counter for included experimental data sets
!   expincfile        ! file with experimental data set to include for search
!   flagedep          ! flag for energy dependence of uncertainties
!   flagscore         ! flag to use scores for inclusion or exclusion
!   flagerror         ! flag to include true experimental uncertainties
!   flagsort          ! flag to sort experimental data in increasing energies
!   lib               ! data library
!   libweight         ! weight for library
!   maxexpsets        ! number of experimental data sets before weight reduction
!   mtallexpweight    ! weight for all channels per experimental channel
!   mtalllibweight    ! weight for all channels per librar
!   mtexpweight       ! weight for experimental data set and channel
!   mtlibweight       ! weight for library and channel
!   mttalweight       ! weight per TALYS channel
! Variables for parameter variation
!   flagequisam       ! flag to use equal binning uniform sampling
!   flaggauss         ! flag to sample parameters from a Gaussian instead of uniform
!   flagkdburn        ! flag to use KD03 OMP function for burn - in
!   flagmcak          ! flag to use sampling from arbitrary distribution
!   flagmetro         ! flag for Metropolis sampling
!   flagsample        ! flag for sampling experimental data
!   flagextparvar     ! flag for external parameter file (by Erik Sunden)
!   getcalcscript     ! script to read pre-calculated TALYS results
! Variables for sensitivity
!   flagreadsens      ! flag to read sensitivities from tables
!   flagsens          ! flag to use sensitivity in parameter weighting
!   flagmorris        ! flag to use Morris screening
! Variables for GOF function
!   flagdeltaE        ! flag to weigh goodness - of - fit with dE of energy grid
! Variables for reading TALYS output
!   flagcross         ! flag for covariances or optimization of cross sections
!   flaggamma         ! flag for covariances or optimization of discrete gamma ray transitions
!   flagld            ! flag for covariances or optimization of level density parameters
!   flagmagnet        ! flag to push TALYS data to experimental or library data
!   flagintegral      ! flag for covariances or optimization of integral activation data
!   flagmacs          ! flag for covariances or optimization of Maxwellian averaged cross sections
!   flagprod          ! flag for covariances or optimization of particle production cross sections
!   flagproduct       ! flag for covariances or optimization of isotope production
!   flagpsf           ! flag for optimization of photon strength functions
!   flagresidual      ! flag for covariances or optimization of residual production cross sections
! Variables for weights
!   flagEweight       ! flag to use energy - dependent weights
!   flagweight        ! flag to use weights for random samples
! Variables for fitting limits
!   flagreadpar       ! flag to read in parameter distribution
! Variables for reading cross sections
!   Efracmax          ! starting energy for maximum cross section
!   flagdiscrete      ! flag for discrete level cross sections
!   flagisomer        ! flag for isomeric production cross sections
!   fracmax           ! fraction of maximum cross section to be included
! Variables for cross section covariances
!   average           ! cross sections 1: average, 2: first run 3: optimum
!   flagband          ! flag to represent results as error bands
! Variables for reading gamma production cross sections
!   maxgam            ! maximum number of discrete levels with gamma cross section
! Variables for writing TALYS input files
!   Ehigh             ! energy cutoff of high energy runs
!   flagangle         ! flag for covariances or optimization of angular distributions
!   flagblock         ! flag to block spectra, angle and gamma files
!   flagecis          ! flag to repeat ECIS calculations in a second run
!   flagglobal        ! flag to use global TALYS cross sections as prior
!   flaginponly       ! flag to produce only TALYS input files, no runs
!   flagleg           ! flag for covariances or optimization of Legendre coefficients
!   flagomponly       ! flag for covariances or optimization for basic optical model observables
!   flagspectra       ! flag for covariances or optimization of emission spectra
!   flagtmc           ! flag for Total Monte Carlo (ENDF - 6 data file at each run)
!   mode              ! TASMAN mode
!   Nhigh             ! number of high energy runs
! Variables for TMC
!   background        ! library name for adoption of back ground cross sections in MF3
!   flagacf           ! flag for creation of ENDF activation file
!   flagcovrand       ! flag for covariance ENDF - 6 file after every random run
!   flagdefonly       ! flag for creation of ENDF default file only
!   flageaf           ! flag for creation of EAF activation file
!   flagfns           ! flag for variation of fission neutron spectrum for ENDF - 6 general purpose file
!   flaggpf           ! flag for creation of ENDF general purpose file
!   flaglrf7          ! flag for resonance representation
!   flagmt            ! flag for creation of ENDF general purpose file without switch
!   flagnjoy          ! flag to run NJOY
!   flagnubar         ! flag for variation of nubar values for ENDF - 6 general purpose file
!   flagprepro        ! flag to run PREPRO codes
!   flagproconly      ! flag for processing of ENDF default file only
!   flagpurr          ! flag to run PURR in NJOY
!   flagres           ! flag for variation of resonance parameters for ENDF - 6 general purpose file
!   flagrunfns        ! flag to produce fns in TASMAN with TANES
!   flagrunnubar      ! flag to produce nubar in TASMAN with TAFIS
!   flagruntares      ! flag to produce resonances in TASMAN with TARES
!   flags20           ! flag for creation of ENDF general purpose file with switch at 20 MeV
!   flags30           ! flag for creation of ENDF general purpose file with switch at 30 MeV
!   flags60           ! flag for creation of ENDF general purpose file with switch at 60 MeV
!   flagsdefault      ! flag for creation of ENDF file with def. switch (30 MeV for neutrons, 0 for other particles)
!   flagselect        ! flag for varying only parts of an ENDF - 6 data file
!   offset            ! offset for numbering or random files (TMC only)
!   tafislib          ! library name for adoption of nubar values (none[default], endfb8.1, jeff4.0, jendl5.0)
!   taneslib          ! library name for adoption of FNS (none[default], endfb8.1, jeff4.0, jendl5.0)
!   tareslib          ! library name for adoption of resonance parameters (default or endfb8.1)
!   tmcoffset         ! offset for starting creation of ENDF - 6 files (TMC only)
! Variables for processing input
!   allexp            ! flag to include all experimental channels
!   allexpweight      ! weight for all experimental channels
!   alllib            ! flag to include all library channels
!   alllibweight      ! weight for all library channels
!   alltal            ! flag to include all TALYS channels
!   alltalweight      ! weight for all TALYS channels
!   changeline        ! line number at which TASMAN specific input starts
!   chi2maxall        ! global maximal chi2 value per point taken into account
!   dsmooth           ! parameter for energy smoothing
!   dsmoothall        ! global parameter for energy smoothing
!   Efracmaxall       ! global starting energy for maximum cross section
!   errlim            ! lower limit for experimental error
!   errlimall         ! global lower limit for experimental error
!   flagautoinc       ! flag to automatically in - or exclude important channels
!   Fmaxall           ! global maximal F value per point taken into account
!   fracmaxall        ! global fraction of maximum cross section to be included
!   liballmt          ! flag for assignment of all MTs per library
!   Liso              ! number of the isomeric state
!   mtallexp          ! weight for all experiments per channel
!   mtalllib          ! weight for all libraries per channel
!   mtlib             ! flag for library per channel
!   mttal             ! weight for TALYS per channel
!   noexp             ! flag to include no experimental channels
!   nolib             ! flag to include no library channels
!   notal             ! flag to include no TALYS channels
!   seed              ! seed for random number generator
!   tafisversion      ! version of TAFIS executable
!   talys             ! TALYS
!   talysversion      ! version of TALYS executable
!   tanesversion      ! version of TANES executable
!   taresversion      ! version of TARES executable
!   tefalversion      ! version of TEFAL executable
!   Tweight           ! weight of TALYS vs experiment
!   Tweightall        ! global weight of TALYS vs experiment
! Variables for weights
!   gofmode           ! gof mode
!   weightpower       ! power for weight of random run
! Variables for experimental data files
!   aamax             ! maximal A value for which exp. data are included in the search
!   aamin             ! minimal A value for which exp. data are included in the search
!   zzmax             ! maximal Z value for which exp. data are included in the search
!   zzmin             ! minimal Z value for which exp. data are included in the search
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical            :: flagset     ! flag for on-set
  logical            :: found_library
  character(len=1)   :: ch          ! character
  character(len=1)   :: w1          ! word
  character(len=3)   :: w3          ! word
  character(len=132) :: key         ! keyword
  character(len=132) :: value       ! value or string
  character(len=132) :: word(40)    ! words on input line
  character(len=132) :: line        ! input line
  integer            :: i           ! counter
  integer            :: ix          ! counter
  integer            :: ilib        ! counter for library
  integer            :: imt         ! MT counter
  integer            :: istat       ! logical for file access
  integer            :: k           ! Legendre order
  integer            :: klib        ! counter for library
  real(sgl)          :: chi2m       ! maximal chi2 value per point taken into account
  real(sgl)          :: e1          ! first energy
  real(sgl)          :: e2          ! second energy
  real(sgl)          :: fm          ! maximal F value per point taken into account
  real(sgl)          :: Tw          ! weight of TALYS vs experiment
  real(sgl)          :: ww          ! weight
!
! ********************** Read input parameters *************************
!
! Set default values
!
! mode=1: Monte Carlo covariances
! mode=2: sensitivity matrix with linear approximation
! mode=3: optimization with experimental data for one nucleus
! mode=4: optimization with experimental data for many nuclides
!
! gofmode 1       : w = exp(-chi2point/2)/w0
! gofmode 2       : w = exp(-(chi2point/chi2point0)/2)/w0
! gofmode 3       : w = exp(-(chi2point/chi2point0)**2)/w0
! gofmode 4       : w = exp(-chi2/2)/w0
! gofmode 5       : w = exp(-(chi2/chi20)/2)/w0
! gofmode 6       : w = exp(-(chi2/chi20)**2)/w0
!
  mode = 1
  Ntalys = -1
  Nburn = -1
  Nhigh = -1
  Ehigh = 1000.
  talysversion = 'talys'
  talys = trim(binpath)//talysversion
  tefalversion = 'tefal'
  tefal = trim(binpath)//tefalversion
  taresversion = 'tares'
  tares = trim(binpath)//taresversion
  tafisversion = 'tafis'
  tafis = trim(binpath)//tafisversion
  tanesversion = 'tanes'
  tanes = trim(binpath)//tanesversion
  flagextparvar = .false.
  getcalcscript = 'none'
  Liso = 0
  maxgam = 4
  minbar = 0
  maxbar = 3
  ichi2 = 1
  gofmode = 1
  searchmode = 4
  maxexpsets = 1000000
  average = 2
  weightpower = 2.
  flagreadsens = .false.
  flagsens = .false.
  flagmorris = .false.
  Nmorris = 20
  flagdeltaE = .false.
  flagdexp = .true.
  flagerf = .true.
  flagband = .false.
  flagscore = .true.
  flagerror = .true.
  flagsort = .true.
  flagcross = .true.
  flagangle = .false.
  flagleg = .false.
  flagspectra = .false.
  flaggamma = .false.
  flagresidual = .false.
  flagprod = .false.
  flagproduct = .false.
  flagdiscrete = .false.
  flagisomer = .false.
  flagintegral = .false.
  flagmacs = .false.
  flaggamgam = .false.
  flagpsf=.false.
  flagld=.false.
  flagecis = .true.
  flagblock = .true.
  flagE1vary = .true.
  flagM1vary = .false.
  flaggauss = .true.
  flagsave = .false.
  flagompinc = .true.
  flagomponly = .false.
  flagtmc = .false.
  flagselect = .false.
  flagcovrand = .false.
  flagsample = .false.
  flaggpf = .true.
  flagsdefault = .true.
  flagdefonly = .true.
  flags20 = .false.
  flags30 = .false.
  flags60 = .false.
  flagmt = .false.
  flagacf = .false.
  flageaf = .false.
  flagprepro = .true.
  flagnjoy = .true.
  flagres = .true.
  tareslib = 'default         '
  background = ''
  flaglrf7 = .false.
  flagpurr = .false.
  flagproconly = .true.
  flagnubar = .false.
  flagrunnubar = .false.
  flagruntares = .false.
  flagfns = .false.
  flagrunfns = .false.
  taneslib = 'none            '
  flaginponly = .false.
  flagmcak = .false.
  flagreadpar = .false.
  flagparvary=.false.
  flagallpar=.false.
  flagweight = .true.
  flagEweight = .false.
  flagmetro = .false.
  flagkdburn = .true.
  flagequisam = .true.
  flagmagnet = .false.
  flagglobal = .false.
  flagautoinc = .true.
  flagedep = .true.
  tafislib = 'none            '
  zzmin = 1
  zzmax = numZ
  aamin = 0
  aamax = 0
  cwidth = 1.
  fiscor0 = 2.
  offset = 0
  tmcoffset = 0
  seed = 20032703
  Fmaxall = 1.e30
  chi2maxall = 1.e30
  errlimall = 0.03
  Efracmaxall = 0.1
  fracmaxall = 0.2
  Esearch1all = -1.
  Esearch2all = 1000.
  Tweightall = 0.5
  dsmoothall = 1.
  alltal = .false.
  notal = .false.
  alllib = .false.
  nolib = .false.
  allexp = .false.
  noexp = .false.
  alltalweight = 1.
  alllibweight = 1.
  allexpweight = 1.
  libweight = 1.
  liballmt = -1
  e1 = 0.
  e2 = 0.
  mtlibweight = 1.
  mtlib = -1
  mtalllibweight = 1.
  mtallexpweight = 1.
  mttalweight = 1.
  mtalllib = -1
  mtallexp = -1
  mttal = -1
  expinccount = 0
  expexccount = 0
  expincfile = ''
  expexcfile = ''
  mtexpweight = 1.
  Fmax = -1.
  chi2max = -1.
  errlim = -1.
  Efracmax = -1.
  fracmax = -1.
  Esearch1 = -1.
  Esearch2 = -1.
  Tweight = -1.
  dsmooth = -1.
  Efracmax(102) = 0.01
  fracmax(1) = 0.1
  fracmax(2) = 0.1
  fracmax(102) = 0.001
  expclass = 3
  outsearch = 5
  source = 'TASMAN'
  oformat = 'YANDF-0.4'
!
! Read input
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
  Nkeyvary = 0
  NZAvary = 0
  NZAskip = 0
  Npartvary = 0
  do i = changeline + 1, nlines
    line = inline(i)
    key = ''
    call getkeywords(line, word)
    key = word(1)
    value = word(2)
    ch = word(2)(1:1)
!
! TASMAN input flags (see manual)
!
    if (key == '#mode') then
      read(value, * , iostat = istat) mode
      if (istat /= 0) call read_error(line, istat)
      call range_integer_error(key, mode, 1, 4)
      cycle
    endif
    if (key == '#ntalys') then
      read(value, * , iostat = istat) Ntalys
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#nmorris') then
      read(value, * , iostat = istat) Nmorris
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#nburn') then
      read(value, * , iostat = istat) Nburn
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#nhigh') then
      read(value, * , iostat = istat) Nhigh
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#ehigh') then
      read(value, * , iostat = istat) Ehigh
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#talysversion') then
      read(value, '(a)' , iostat = istat) talys
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#tefalversion') then
      read(value, '(a)' , iostat = istat) tefal
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#taresversion') then
      read(value, '(a)' , iostat = istat) tares
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#tafisversion') then
      read(value, '(a)' , iostat = istat) tafis
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#tanesversion') then
      read(value, '(a)' , iostat = istat) tanes
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#liso') then
      read(value, * , iostat = istat) Liso
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#maxbar') then
      read(value, * , iostat = istat) maxbar
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#minbar') then
      read(value, * , iostat = istat) minbar
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#maxgam') then
      read(value, * , iostat = istat) maxgam
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#chi2') then
      read(value, * , iostat = istat) ichi2
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#gofmode') then
      read(value, * , iostat = istat) gofmode
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#searchmode') then
      read(value, * , iostat = istat) searchmode
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#maxexpsets') then
      read(value, * , iostat = istat) maxexpsets
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#average') then
      read(value, * , iostat = istat) average
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#weightpower') then
      read(value, * , iostat = istat) weightpower
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#readsens') then
      if (ch == 'n') flagreadsens = .false.
      if (ch == 'y') flagreadsens = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#sens') then
      if (ch == 'n') flagsens = .false.
      if (ch == 'y') flagsens = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#morris') then
      if (ch == 'n') flagmorris = .false.
      if (ch == 'y') flagmorris = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#deltae') then
      if (ch == 'n') flagdeltaE = .false.
      if (ch == 'y') flagdeltaE = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#dexp') then
      if (ch == 'n') flagdexp = .false.
      if (ch == 'y') flagdexp = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#erf') then
      if (ch == 'n') flagerf = .false.
      if (ch == 'y') flagerf = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#band') then
      if (ch == 'n') flagband = .false.
      if (ch == 'y') flagband = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#score') then
      if (ch == 'n') flagscore = .false.
      if (ch == 'y') flagscore = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#error') then
      if (ch == 'n') flagerror = .false.
      if (ch == 'y') flagerror = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#sort') then
      if (ch == 'n') flagsort = .false.
      if (ch == 'y') flagsort = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#cross') then
      if (ch == 'n') flagcross = .false.
      if (ch == 'y') flagcross = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#angle') then
      if (ch == 'n') flagangle = .false.
      if (ch == 'y') flagangle = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#legendre') then
      if (ch == 'n') flagleg = .false.
      if (ch == 'y') flagleg = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#spectra') then
      if (ch == 'n') flagspectra = .false.
      if (ch == 'y') flagspectra = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#gamma') then
      if (ch == 'n') flaggamma = .false.
      if (ch == 'y') flaggamma = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#residual') then
      if (ch == 'n') flagresidual = .false.
      if (ch == 'y') flagresidual = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#prod') then
      if (ch == 'n') flagprod = .false.
      if (ch == 'y') flagprod = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#production') then
      if (ch == 'n') flagproduct = .false.
      if (ch == 'y') flagproduct = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#discrete') then
      if (ch == 'n') flagdiscrete = .false.
      if (ch == 'y') flagdiscrete = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#isomer') then
      if (ch == 'n') flagisomer = .false.
      if (ch == 'y') flagisomer = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#integral') then
      if (ch == 'n') flagintegral = .false.
      if (ch == 'y') flagintegral = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#macs') then
      if (ch == 'n') flagmacs = .false.
      if (ch == 'y') flagmacs = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#gamgam') then
      if (ch == 'n') flaggamgam = .false.
      if (ch == 'y') flaggamgam = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#psf') then
      if (ch == 'n') flagpsf = .false.
      if (ch == 'y') flagpsf = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#ld') then
      if (ch == 'n') flagld = .false.
      if (ch == 'y') flagld = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#ecis') then
      if (ch == 'n') flagecis = .false.
      if (ch == 'y') flagecis = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#block') then
      if (ch == 'n') flagblock = .false.
      if (ch == 'y') flagblock = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#m1vary') then
      if (ch == 'n') flagM1vary = .false.
      if (ch == 'y') flagM1vary = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#e1vary') then
      if (ch == 'n') flagE1vary = .false.
      if (ch == 'y') flagE1vary = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#gauss') then
      if (ch == 'n') flaggauss = .false.
      if (ch == 'y') flaggauss = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#save') then
      if (ch == 'n') flagsave = .false.
      if (ch == 'y') flagsave = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#ompinc') then
      if (ch == 'n') flagompinc = .false.
      if (ch == 'y') flagompinc = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#omponly') then
      if (ch == 'n') flagomponly = .false.
      if (ch == 'y') flagomponly = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#tmc') then
      if (ch == 'n') flagtmc = .false.
      if (ch == 'y') flagtmc = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#select') then
      if (ch == 'n') flagselect = .false.
      if (ch == 'y') flagselect = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#covrand') then
      if (ch == 'n') flagcovrand = .false.
      if (ch == 'y') flagcovrand = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#sample') then
      if (ch == 'n') flagsample = .false.
      if (ch == 'y') flagsample = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#gpf') then
      if (ch == 'n') flaggpf = .false.
      if (ch == 'y') flaggpf = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#sdefault') then
      if (ch == 'n') flagsdefault = .false.
      if (ch == 'y') flagsdefault = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#defonly') then
      if (ch == 'n') flagdefonly = .false.
      if (ch == 'y') flagdefonly = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#s20') then
      if (ch == 'n') flags20 = .false.
      if (ch == 'y') flags20 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#s30') then
      if (ch == 'n') flags30 = .false.
      if (ch == 'y') flags30 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#s60') then
      if (ch == 'n') flags60 = .false.
      if (ch == 'y') flags60 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#mt') then
      if (ch == 'n') flagmt = .false.
      if (ch == 'y') flagmt = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#acf') then
      if (ch == 'n') flagacf = .false.
      if (ch == 'y') flagacf = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#eaf') then
      if (ch == 'n') flageaf = .false.
      if (ch == 'y') flageaf = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#prepro') then
      if (ch == 'n') flagprepro = .false.
      if (ch == 'y') flagprepro = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#njoy') then
      if (ch == 'n') flagnjoy = .false.
      if (ch == 'y') flagnjoy = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#resonance') then
      if (ch == 'n') flagres = .false.
      if (ch == 'y') flagres = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#tareslib') then
      read(value, * , iostat = istat) tareslib
      cycle
    endif
    if (key == '#background') then
      read(value, * , iostat = istat) background
      cycle
    endif
    if (key.eq.'#extparvar') then
      if (ch.eq.'n') flagextparvar = .false.
      if (ch.eq.'y') flagextparvar = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#getcalcscript') then
      read(value, * , iostat = istat) getcalcscript
      cycle
    endif
    if (key == '#lrf7') then
      if (ch == 'n') flaglrf7 = .false.
      if (ch == 'y') flaglrf7 = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#purr') then
      if (ch == 'n') flagpurr = .false.
      if (ch == 'y') flagpurr = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#proconly') then
      if (ch == 'n') flagproconly = .false.
      if (ch == 'y') flagproconly = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#nubar') then
      if (ch == 'n') flagnubar = .false.
      if (ch == 'y') flagnubar = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#runnubar') then
      if (ch == 'n') flagrunnubar = .false.
      if (ch == 'y') flagrunnubar = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#runtares') then
      if (ch == 'n') flagruntares = .false.
      if (ch == 'y') flagruntares = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#tafislib') then
      read(value, * , iostat = istat) tafislib
      cycle
    endif
    if (key == '#fns') then
      if (ch == 'n') flagfns = .false.
      if (ch == 'y') flagfns = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#mcak') then
      if (ch == 'n') flagmcak = .false.
      if (ch == 'y') flagmcak = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#readpar') then
      if (ch == 'n') flagreadpar = .false.
      if (ch == 'y') flagreadpar = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#parameters') then
      if (ch == 'n') flagparvary = .false.
      if (ch == 'y') flagparvary = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#allpar') then
      if (ch == 'n') flagallpar = .false.
      if (ch == 'y') flagallpar = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#weight') then
      if (ch == 'n') flagweight = .false.
      if (ch == 'y') flagweight = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#eweight') then
      if (ch == 'n') flageweight = .false.
      if (ch == 'y') flageweight = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#metropolis') then
      if (ch == 'n') flagmetro = .false.
      if (ch == 'y') flagmetro = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#kdburn') then
      if (ch == 'n') flagkdburn = .false.
      if (ch == 'y') flagkdburn = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#equisample') then
      if (ch == 'n') flagequisam = .false.
      if (ch == 'y') flagequisam = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#magnet') then
      if (ch == 'n') flagmagnet = .false.
      if (ch == 'y') flagmagnet = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#global') then
      if (ch == 'n') flagglobal = .false.
      if (ch == 'y') flagglobal = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#autoinc') then
      if (ch == 'n') flagautoinc = .false.
      if (ch == 'y') flagautoinc = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#edependent') then
      if (ch == 'n') flagedep = .false.
      if (ch == 'y') flagedep = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#runfns') then
      if (ch == 'n') flagrunfns = .false.
      if (ch == 'y') flagrunfns = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#inponly') then
      if (ch == 'n') flaginponly = .false.
      if (ch == 'y') flaginponly = .true.
      if (ch /= 'y' .and. ch /= 'n') call read_error(line, istat)
      cycle
    endif
    if (key == '#taneslib') then
      read(value, * , iostat = istat) taneslib
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#zmin') then
      read(value, * , iostat = istat) zzmin
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#zmax') then
      read(value, * , iostat = istat) zzmax
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#amin') then
      read(value, * , iostat = istat) aamin
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#amax') then
      read(value, * , iostat = istat) aamax
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#cwidth') then
      read(value, * , iostat = istat) cwidth
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#fiscor') then
      read(value, * , iostat = istat) fiscor0
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#offset') then
      read(value, * , iostat = istat) offset
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#tmcoffset') then
      read(value, * , iostat = istat) tmcoffset
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#seed') then
      read(value, * , iostat = istat) seed
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#zavary') then
      NZAvary = NZAvary + 1
      read(word(2), *, iostat = istat) Zvary(NZAvary)
      if (istat /= 0) call read_error(line, istat)
      read(word(3), *, iostat = istat) Avary(NZAvary)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#zaskip') then
      NZAskip = NZAskip + 1
      read(word(2), *, iostat = istat) Zskip(NZAskip)
      if (istat /= 0) call read_error(line, istat)
      read(word(3), *, iostat = istat) Askip(NZAskip)
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#partvary') then
      Npartvary = Npartvary + 1
      partvary(Npartvary) = ch
      cycle
    endif
    if (key == '#keyvary') then
      Nkeyvary = Nkeyvary + 1
      read(value, *, iostat = istat) keyvary(Nkeyvary)
      if (istat /= 0) call read_error(line, istat)
      if (word(3)(1:1) /= ' ') read(word(3), *, iostat = istat) sigvary(Nkeyvary)
      if (istat /= 0) call read_error(line, istat)
      line = keyvary(Nkeyvary)
      call convert(line)
      keyvary(Nkeyvary) = line
      cycle
    endif
!
! ********* Include data libraries, MT numbers and data sets ***********
!
! 1. Include libraries
!
    if (key == '#libinclude') then
      imt = 0
      k = 2
!
! MT number
!
      if (ch >= '0' .and. ch <= '9') then
        read(value, * , iostat = istat) imt
        call range_integer_error(key, imt, 0, nummt)
        k = k + 1
      endif
!
! Libraries
!
      klib = 0
      w3 = word(k)(1:3)
      if (w3 == 'all' .or. w3(1:1) == ' ' .or. (w3(1:1) >= '0' .and. w3(1:1) <= '9')) then
        if (imt > 0) then
          mtalllib(imt) = 1
        else
          alllib = .true.
        endif
        if (w3 == 'all') k = k+1
        if (w3(1:1) == ' ') cycle
      else
!
! Library name
!
        do ilib = 1, numlib
          found_library = .false.
          if (trim(lib(ilib)) == trim(word(k))) then
            found_library = .true.
            klib = ilib
            if (imt > 0) then
              mtlib(imt, ilib) = 1
            else
              liballmt(ilib) = 1
            endif
            k = k + 1
            exit
          endif
        enddo
        if (.not. found_library) call read_error(line, istat)
      endif
!
! Weights
!
      if (word(k)(1:1) /= ' ') then
        read(word(k), * , iostat = istat) ww
        if (imt > 0) then
          if (klib > 0) then
            mtlibweight(imt, klib) = ww
          else
            mtalllibweight(imt) = ww
          endif
        else
          if (w3 == 'all') alllibweight = ww
          if (klib > 0) libweight(klib) = ww
        endif
      endif
      cycle
    endif
!
! 2. Exclude libraries
!
    if (key == '#libexclude') then
      imt = 0
      k = 2
!
! MT number
!
      if (ch >= '0' .and. ch <= '9') then
        read(value, * , iostat = istat) imt
        call range_integer_error(key, imt, 0, nummt)
        k = k + 1
      endif
!
! Libraries
!
      w3 = word(k)(1:3)
      if (w3 == 'all' .or. w3(1:1) == ' ') then
        if (imt > 0) then
          mtalllib(imt) = 0
        else
          nolib = .true.
        endif
        if (w3 == 'all') k = k+1
        if (w3(1:1) == ' ') cycle
        cycle
      endif
!
! Library name
!
      klib = 0
      found_library = .false.
      do ilib = 1, numlib
        if (trim(lib(ilib)) == trim(word(k))) then
          found_library = .true.
          klib = ilib
          if (imt > 0) then
            mtlib(imt, ilib) = 0
          else
            liballmt(ilib) = 0
          endif
          k = k + 1
          cycle
        endif
      enddo
      if (.not. found_library) call read_error(line, istat)
      cycle
    endif
!
! 3. Include experimental data
!
    if (key == '#expinclude') then
      imt = 0
      k = 2
!
! MT number
!
      if (ch >= '0' .and. ch <= '9') then
        read(value, * , iostat = istat) imt
        call range_integer_error(key, imt, 0, nummt)
        k = k + 1
      endif
!
! Experimental data
!
      flagset = .false.
      w3 = word(k)(1:3)
      if (w3 == 'all' .or. w3(1:1) == ' ' .or. (w3(1:1) >= '0' .and. w3(1:1) <= '9')) then
        if (imt > 0) then
          mtallexp(imt) = 1
        else
          allexp = .true.
        endif
        if (w3 == 'all') k = k+1
        if (w3(1:1) == ' ') cycle
      else
!
! Specific experimental data sets
!
        if (imt > 0) then
          w1 = word(k)(1:1)
          if ((w1 >= 'A' .and. w1 <= 'Z') .or. (w1 >= 'a' .and. w1 <= 'z')) then
            if (expinccount(imt) < numsets) then
              expinccount(imt) = expinccount(imt) + 1
              read(word(k), * , iostat = istat) expincfile(imt, expinccount(imt))
              flagset = .true.
              k = k + 1
            endif
          endif
        endif
      endif
!
! Weights
!
      if (word(k)(1:1) /= ' ') then
        read(word(k), * , iostat = istat) ww
        if (imt > 0) then
          if (flagset) then
            mtexpweight(imt, expinccount(imt)) = ww
          else
            mtallexpweight(imt) = ww
          endif
        else
          if (w3 == 'all') allexpweight = ww
        endif
      endif
      cycle
    endif
!
! 4. Exclude experimental data
!
    if (key == '#expexclude') then
      imt = 0
      k = 2
!
! MT number
!
      if (ch >= '0' .and. ch <= '9') then
        read(value, * , iostat = istat) imt
        call range_integer_error(key, imt, 0, nummt)
        k = k + 1
      endif
!
! Experimental data
!
      flagset = .false.
      w3 = word(k)(1:3)
      if (w3 == 'all' .or. w3(1:1) == ' ') then
        if (imt > 0) then
          mtallexp(imt) = 0
        else
          noexp = .true.
        endif
        if (w3 == 'all') k = k+1
        if (w3(1:1) == ' ') cycle
      endif
!
! Specific experimental data sets
!
      if (imt > 0) then
        w1 = word(k)(1:1)
        if ((w1 >= 'A' .and. w1 <= 'Z') .or. (w1 >= 'a' .and. w1 <= 'z')) then
          if (expexccount(imt) < numsets) then
            expexccount(imt) = expexccount(imt) + 1
            read(word(k), * , iostat = istat) expexcfile(imt, expexccount(imt))
            flagset = .true.
            k = k + 1
          endif
        endif
      endif
    endif
!
! 5. Include central values of TALYS
!
    if (key == '#talinclude') then
      imt = 0
      k = 2
!
! MT number
!
      if (ch >= '0' .and. ch <= '9') then
        read(value, * , iostat = istat) imt
        call range_integer_error(key, imt, 0, nummt)
        k = k + 1
      endif
!
! All or per MT number
!
      w3 = word(k)(1:3)
      if (w3 == 'all' .or. w3(1:1) == ' ' .or. (w3(1:1) >= '0' .and. w3(1:1) <= '9')) then
        if (imt > 0) then
          mttal(imt) = 1
        else
          alltal = .true.
        endif
        if (w3 == 'all') k = k+1
        if (w3(1:1) == ' ') cycle
      endif
!
! Weights
!
      if (word(k)(1:1) /= ' ') then
        read(word(k), * , iostat = istat) ww
        if (imt > 0) then
          mttalweight(imt) = ww
        else
          if (w3 == 'all') alltalweight = ww
        endif
      endif
      cycle
    endif
!
! 6. Exclude central values of TALYS
!
    if (key == '#talexclude') then
      imt = 0
      k = 2
!
! MT number
!
      if (ch >= '0' .and. ch <= '9') then
        read(value, * , iostat = istat) imt
        call range_integer_error(key, imt, 0, nummt)
        k = k + 1
      endif
!
! All or per MT number
!
      w3 = word(k)(1:3)
      if (w3 == 'all' .or. w3(1:1) == ' ') then
        if (imt > 0) then
          mttal(imt) = 0
        else
          notal = .true.
        endif
        if (w3 == 'all') k = k+1
      endif
      cycle
    endif
    if (key == '#zdeep') then
      read(value, * , iostat = istat) Zdeep
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#adeep') then
      read(value, * , iostat = istat) Adeep
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#expclass') then
      read(value, * , iostat = istat) expclass
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#outsearch') then
      read(value, * , iostat = istat) outsearch
      if (istat /= 0) call read_error(line, istat)
      cycle
    endif
    if (key == '#esearch') then
      imt = 0
      read(value, * , iostat = istat) e1
      if (istat > 0) call read_error(line, istat)
      read(word(3), * , iostat = istat) e2
      if (istat > 0) call read_error(line, istat)
      read(word(4), * , iostat = istat) imt
      if (istat > 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      if (imt == 0) then
        Esearch1all = e1
        Esearch2all = e2
      else
        Esearch1(imt) = e1
        Esearch2(imt) = e2
      endif
      cycle
    endif
    if (key == '#esearch1') then
      imt = 0
      read(value, * , iostat = istat) e1
      if (istat > 0) call read_error(line, istat)
      read(word(3), * , iostat = istat) imt
      if (istat > 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      if (imt == 0) then
        Esearch1all = e1
      else
        Esearch1(imt) = e1
      endif
      cycle
    endif
    if (key == '#esearch2') then
      imt = 0
      read(value, * , iostat = istat) e2
      if (istat > 0) call read_error(line, istat)
      read(word(3), * , iostat = istat) imt
      if (istat > 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      if (imt == 0) then
        Esearch2all = e2
      else
        Esearch2(imt) = e2
      endif
      cycle
    endif
    if (key == '#fmax') then
      read(value, * , iostat = istat) fm
      if (istat > 0) call read_error(line, istat)
      imt = 0
      read(word(3), * , iostat = istat) imt
      if (istat > 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      if (imt == 0) then
        Fmaxall = fm
      else
        Fmax(imt) = fm
      endif
      cycle
    endif
    if (key == '#chi2max') then
      read(value, * , iostat = istat) chi2m
      if (istat > 0) call read_error(line, istat)
      imt = 0
      read(word(3), * , iostat = istat) imt
      if (istat > 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      if (imt == 0) then
        chi2maxall = chi2m
      else
        chi2max(imt) = chi2m
      endif
      cycle
    endif
    if (key == '#errlim') then
      read(value, * , iostat = istat) fm
      if (istat > 0) call read_error(line, istat)
      imt = 0
      read(word(3), * , iostat = istat) imt
      if (istat > 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      if (imt == 0) then
        errlimall = fm
      else
        errlim(imt) = fm
      endif
      cycle
    endif
    if (key == '#efracmax') then
      read(value, * , iostat = istat) fm
      if (istat > 0) call read_error(line, istat)
      imt = 0
      read(word(3), * , iostat = istat) imt
      if (istat > 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      if (imt == 0) then
        Efracmaxall = fm
      else
        Efracmax(imt) = fm
      endif
      cycle
    endif
    if (key == '#fracmax') then
      read(value, * , iostat = istat) fm
      if (istat > 0) call read_error(line, istat)
      imt = 0
      read(word(3), * , iostat = istat) imt
      if (istat > 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      if (imt == 0) then
        fracmaxall = fm
      else
        fracmax(imt) = fm
      endif
      cycle
    endif
    if (key == '#tweight') then
      read(value, * , iostat = istat) Tw
      if (istat > 0) call read_error(line, istat)
      imt = 0
      read(word(3), * , iostat = istat) imt
      if (istat > 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      if (imt == 0) then
        Tweightall = Tw
      else
        Tweight(imt) = Tw
      endif
      cycle
    endif
    if (key == '#dsmooth') then
      read(value, * , iostat = istat) Tw
      if (istat > 0) call read_error(line, istat)
      imt = 0
      read(word(3), * , iostat = istat) imt
      if (istat > 0) call read_error(line, istat)
      call range_integer_error(key, imt, 0, nummt)
      if (imt == 0) then
        dsmoothall = Tw
      else
        dsmooth(imt) = Tw
      endif
      cycle
    endif
    if (key == '#source') then
      ix=index(line,'source')+7
      source=trim(adjustl(line(ix:132)))
      cycle
    endif
    if (key == '#user') then
      ix=index(line,'user')+5
      user=trim(adjustl(line(ix:132)))
      cycle
    endif
    if (key == '#format') then
      ix=index(line,'format')+7
      oformat=trim(adjustl(line(ix:132)))
      cycle
    endif
  enddo
  return
end subroutine input2
! Copyright A.J. Koning 2021
