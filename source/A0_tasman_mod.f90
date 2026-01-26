module A0_tasman_mod
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : General module with all global variables
!
! Author    : Arjan Koning
!
! 2023-12-29: Original code
! 2026-01-25: Current revision
!-----------------------------------------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Definition of single and double precision variables
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer, parameter :: dbl = selected_real_kind(15,307) ! double precision kind
!
!-----------------------------------------------------------------------------------------------------------------------------------
! All global dimensions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter :: numZ=124        ! maximum number of elements
  integer, parameter :: nump=6          ! maximum number of particles
  integer, parameter :: numchancov=14   ! maximum number of channels with covariances
  integer, parameter :: numlib=11       ! maximum number of data libraries
  integer, parameter :: nummt=850       ! maximum number of MT numbers
  integer, parameter :: numencov=52     ! maximum number of covariance energies
  integer, parameter :: numlines=5000   ! maximum number of input lines
  integer, parameter :: numenin=600     ! maximum number of incident energies
  integer, parameter :: numpar=500      ! maximum number of parameters
  integer, parameter :: numtalys=6000   ! maximum number of TALYS runs
  integer, parameter :: nummass=414     ! maximum number of masses
  integer, parameter :: numchanxs=200   ! maximum number of channels with cross sections
  integer, parameter :: numchanrp=200   ! maximum number of channels with residual production cross sections
  integer, parameter :: numchangam=100  ! maximum number of channels with gamma cross sections
  integer, parameter :: numchanprod=7   ! maximum number of channels with particle production cross sections
  integer, parameter :: numchansp=7     ! maximum number of channels with emission spectra
  integer, parameter :: numchanang=52   ! maximum number of channels with angular distributions
  integer, parameter :: numnuc=1000     ! maximum number of nuclides with experimental data
  integer, parameter :: numenendf=10000 ! maximum number of energies on ENDF grid
  integer, parameter :: numsets=100     ! maximum number of experimental data sets
  integer, parameter :: numenexp=100    ! maximum number of energies forexperimental data set
  integer, parameter :: numchanexp=40   ! maximum number of channels with experimental data
  integer, parameter :: numang=90       ! maximum number of angles
  integer, parameter :: numenS=60       ! maximum number of energies for sensitivities
  integer, parameter :: numchanY=200    ! maximum number of channels with isotope production
  integer, parameter :: numtime=100     ! maximum number of time points
  integer, parameter :: numen2=650      ! maximum number of energies for spectra
  integer, parameter :: numleg=6        ! maximum number of Legendre coefficients
  integer, parameter :: numhist=100     ! maximum number of points for histogram
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=1), dimension(0:nump) :: particle ! type of particle
  character(len=2), dimension(numZ)   :: nuc      ! symbol of nucleus
  integer, dimension(numZ)            :: heavy    ! heaviest isotope
  integer, dimension(numZ)            :: light    ! lightest isotope
  integer, dimension(numZ)            :: mainis   ! main isotope
  real(sgl)                           :: xseps    ! limit for cross sections
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for pathnames
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=132) :: librariespath ! directory containing files to be read
  character(len=132) :: exforpath     ! directory containing files to be read
  character(len=132) :: binpath       ! directory containing files to be read
  character(len=132) :: tasmanpath    ! directory containing files to be read
  character(len=132) :: psfpath       ! directory containing files to be read
  integer            :: year          ! year
  integer            :: month         ! month
  integer            :: day           ! day
  character(len=10)  :: date      ! date
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading TASMAN input lines
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=132), dimension(numlines) :: inline ! input line
  integer                                 :: nlines ! number of input lines
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for processing input
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                          :: allexp       ! flag to include all experimental channels
  logical                          :: alllib       ! flag to include all library channels
  logical                          :: alltal       ! flag to include all TALYS channels
  logical, dimension(numenin)      :: anginclude   ! flag to include energy point for angular grid
  logical                          :: flagautoinc  ! flag to automatically in- or exclude important channels
  logical                          :: noexp        ! flag to include no experimental channels
  logical                          :: nolib        ! flag to include no library channels
  logical                          :: notal        ! flag to include no TALYS channels
  character(len=1)                 :: isochar      ! symbol for isomer
  character(len=1)                 :: ptype0       ! type of incident particle
  character(len=2)                 :: Starget      ! symbol of target nucleus
  character(len=6)                 :: targetnuclide
  character(len=16), dimension(nummt) :: MTreac    ! reaction string
  character(len=132)               :: talysversion ! version of TALYS executable
  character(len=132)               :: tefalversion ! version of TEFAL executable
  character(len=132)               :: tafisversion ! version of TAFIS executable
  character(len=132)               :: tanesversion ! version of TANES executable
  character(len=132)               :: taresversion ! version of TARES executable
  character(len=132)               :: tafis        ! TAFIS executable
  character(len=132)               :: talys        ! TALYS executable
  character(len=132)               :: tanes        ! TANES executable
  character(len=132)               :: tares        ! TAREL executable
  character(len=132)               :: tefal        ! TEFAL executable
  character(len=132)               :: energyfile   ! file with incident energies
  character(len=132)               :: exppath      ! path for experimental data
  character(len=132)               :: source       ! source of data
  character(len=132)               :: oformat      ! format of data
  character(len=132)               :: user         ! user of data
  integer                          :: Atarget      ! mass number of target nucleus
  integer, dimension(nummt)        :: auto         ! flag to in- or exclude channel
  integer                          :: changeline   ! line number at which TASMAN specific input starts
  integer                          :: k0           ! index of incident particle
  integer                          :: lenexp       ! length of pathname for experimental data
  integer                          :: lennuc       ! length of nuclear symbol
  integer, dimension(numlib)       :: liballmt     ! flag for assignment of all MTs per library
  integer, dimension(numlib)       :: libinclude   ! flag to include data library for search
  integer                          :: Liso         ! number of the isomeric state
  integer                          :: Ltarget      ! excited level of target
  integer, dimension(nummt)        :: mtallexp     ! weight for all experiments per channel
  integer, dimension(nummt)        :: mtalllib     ! weight for all libraries per channel
  integer, dimension(nummt,numlib) :: mtlib        ! flag for library per channel
  integer, dimension(nummt)        :: mttal        ! weight for TALYS per channel
  integer                          :: numE         ! number of incident energies
  integer                          :: seed         ! seed for random number generator
  integer                          :: Ztarget      ! charge number of target nucleus
  real(sgl)                        :: allexpweight ! weight for all experimental channels
  real(sgl)                        :: alllibweight ! weight for all library channels
  real(sgl)                        :: alltalweight ! weight for all TALYS channels
  real(sgl)                        :: chi2maxall   ! global maximal chi2 value per point taken into account
  real(sgl), dimension(nummt)      :: dsmooth      ! parameter for energy smoothing
  real(sgl)                        :: dsmoothall   ! global parameter for energy smoothing
  real(sgl)                        :: Efracmaxall  ! global starting energy for maximum cross section
  real(sgl)                        :: errlimall    ! global lower limit for experimental error
  real(sgl)                        :: Fmaxall      ! global maximal F value per point taken into account
  real(sgl)                        :: fracmaxall   ! global fraction of maximum cross section to be included
  real(sgl), dimension(nummt)      :: errlim       ! lower limit for experimental error
  real(sgl), dimension(nummt)      :: Tweight      ! weight of TALYS vs experiment
  real(sgl)                        :: Tweightall   ! global weight of TALYS vs experiment
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for writing TALYS input files
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                              :: first       ! tracks whether a TALYS run is the first one
  logical                              :: flagangle   ! flag for covariances or optimization of angular distributions
  logical                              :: flagecis    ! flag to repeat ECIS calculations in a second run
  logical                              :: flagglobal  ! flag to use global TALYS cross sections as prior
  logical                              :: flaginponly ! flag to produce only TALYS input files, no runs
  logical                              :: flagleg     ! flag for covariances or optimization of Legendre coefficients
  logical                              :: flagompinc  ! flag to vary OMP parameters for incident channel
  logical                              :: flagomponly ! flag for covariances or optimization for basic OMP only
  logical                              :: flagompvar  ! flag for variation of OMP parameters
  logical                              :: flagspectra ! flag for covariances or optimization of emission spectra
  logical                              :: flagtmc     ! flag for Total Monte Carlo (ENDF-6 data file at each run)
  character(len=1), dimension(numpar)  :: parsym      ! particle symbol
  character(len=2), dimension(numpar)  :: pargam      ! gamma-ray related input string
  character(len=132), dimension(numpar):: string      ! line with parameter value
  integer, dimension(numnuc)           :: A           ! mass number of nucleus
  integer, dimension(numpar)           :: Apar        ! A-index for keyword
  integer, dimension(numpar)           :: barrier     ! fission barrier
  integer                              :: inuc        ! index of nuclide
  integer                              :: italys      ! TALYS run counter
  integer, dimension(numpar)           :: lval        ! l value
  integer                              :: Nhigh       ! number of high energy runs
  integer                              :: mode        ! TASMAN mode
  integer, dimension(numpar)           :: parmt       ! MT number
  integer, dimension(numnuc)           :: Z           ! charge number of nucleus
  integer, dimension(numpar)           :: Zpar        ! Z-index for keyword
  real(sgl)                            :: Ehigh       ! energy cutoff of high energy runs
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading TALYS output
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical              :: flagcross    ! flag for covariances or optimization of cross sections
  logical              :: flagexp      ! flag to optimize to experimental data
  logical              :: flaggamma    ! flag for covariances or optimization of discrete gamma ray transitions
  logical              :: flagintegral ! flag for covariances or optimization of integral activation data
  logical              :: flagld       ! flag for optimization of level density parameters
  logical              :: flaglib      ! flag to optimize to nuclear data libraries
  logical              :: flagmacs     ! flag for covariances or optimization of Maxwellian averaged cross sections
  logical              :: flaggamgam   ! flag for covariances or optimization of average radiative width
  logical              :: flagmagnet   ! flag to push TALYS data to experimental or library data
  logical              :: flagprod     ! flag for covariances or optimization of particle production cross sections
  logical              :: flagproduct  ! flag for covariances or optimization of isotope production
  logical              :: flagpsf      ! flag for optimization of photon strength functions
  logical              :: flagresidual ! flag for covariances or optimization of residual production cross sections
  logical              :: flagtal      ! flag to optimize to TALYS
  logical              :: flagtalys    ! flag to track successful TALYS calculation
  integer              :: iout         ! counter for output files
  integer              :: Nchanall     ! total number of channels
  integer              :: Ntalbeg      ! counter for first TALYS calculation
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for parameter variation
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                 :: flagequisam   ! flag to use equal binning uniform sampling
  logical                                 :: flaggauss     ! flag to sample parameters from a Gaussian instead of uniform
  logical                                 :: flagkdburn    ! flag to use KD03 OMP function for burn-in
  logical                                 :: flagmcak      ! flag to use sampling from arbitrary distribution
  logical                                 :: flagmetro     ! flag for Metropolis sampling
  logical                                 :: flagsample    ! flag for sampling experimental data
  logical                                 :: flagextparvar ! flag for external parameter file (by Erik Sunden)
  character(len=6), dimension(numpar)     :: partype       ! type of parameter variation
  character(len=12), dimension(numpar)    :: parinpname    ! input character parameter value
  character(len=16), dimension(numpar)    :: parname       ! character parameter value
  character(len=132), dimension(numpar)   :: parkey        ! TALYS keyword
  integer                                 :: Npar          ! number of parameters
  real(sgl), dimension(numpar)            :: par           ! parameter value
  real(sgl), dimension(numpar,4)          :: paradjust     ! local adjustment parameters
  real(sgl), dimension(numpar)            :: parD          ! random parameter value
  real(sgl), dimension(numpar)            :: pardelta      ! uncertainty of parameter
  real(sgl), dimension(numpar)            :: parinp        ! input parameter value
  real(sgl), dimension(0:numtalys,numpar) :: parsave       ! parameter value
  real(sgl), dimension(0:1,numpar)        :: partalys      ! parameter value
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for deviation
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(nummt) :: adevi ! parameter for energy dependence of deviation
  real(sgl), dimension(nummt) :: bdevi ! parameter for energy dependence of deviation
  real(sgl), dimension(nummt) :: cdevi ! parameter for energy dependence of deviation
  real(sgl), dimension(nummt) :: ddevi ! parameter for energy dependence of deviation
  real(sgl), dimension(nummt) :: devi  ! parameter for energy dependence of deviation
  real(sgl), dimension(nummt) :: Ecent ! parameter for energy dependence of deviation
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for automatic parameter variation
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                               :: flagparvary ! flag to automatically vary parameters
  logical                               :: flagallpar  ! flag to vary all nuclear model parameters
  logical                               :: flagE1vary  ! flag to automatically vary keywords for E1 radiation
  logical                               :: flagM1vary  ! flag to automatically vary keywords for M1 radiation
  character(len=1), dimension(numpar)   :: partvary    ! particle for varying parameters
  character(len=132), dimension(numpar) :: keyvary     ! keyword to vary
  integer                               :: ACN         ! A of compound nucleus
  integer                               :: Adeep       ! maximum depth of A
  integer                               :: ZCN         ! Z of compound nucleus
  integer                               :: Zdeep       ! maximum depth of Z
  integer                               :: Nkeyvary    ! number of keywords to vary
  integer                               :: Npartvary   ! number of particles to vary keywords
  integer                               :: NZAskip     ! number of Z,A pairs to skip for parameter variation
  integer                               :: NZAvary     ! number of Z,A pairs for parameter variation
  integer                               :: minbar      ! minimal fission barrier
  integer                               :: maxbar      ! maximal fission barrier
  integer, dimension(numpar)            :: Askip       ! A value to skip
  integer, dimension(numpar)            :: Zskip       ! Z value to skip
  integer, dimension(numpar)            :: Avary       ! A value to vary
  integer, dimension(numpar)            :: Zvary       ! Z value to vary
  real(sgl), dimension(numpar)          :: sigvary     ! sigma for parameter variation
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for uncertainty
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical            :: flagsave      ! flag to save all input and output files of each Talys run
  character(len=132) :: getcalcscript ! script to read pre-calculated TALYS results
  real(sgl)          :: cwidth        ! scaling constant for parameter variation
  real(sgl)          :: fiscor0       ! adjustment parameter for fission uncertainties
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for weights
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical, dimension(numchanxs)                    :: flagchmagnet ! flag to unify theoretical and experimental data
  logical                                          :: flagEweight  ! flag to use energy-dependent weights
  logical                                          :: flagweight   ! flag to use weights for random samples
  integer                                          :: gofmode      ! gof mode
  integer, dimension(numchanxs)                    :: Nxssamp      ! number of sampled cross sections
  real(sgl), dimension(numchanxs,numsets*numenexp) :: Esamp        ! energy of sampled cross section
  real(sgl), dimension(0:numtalys,numchanxs)       :: Gchannelsave ! GOF for channel
  real(sgl), dimension(0:numtalys)                 :: Gsave        ! GOF value
  real(sgl)                                        :: weightpower  ! power for weight of random run
  real(sgl), dimension(0:numtalys)                 :: weightsave   ! weight for TALYS run
  real(sgl), dimension(numchanxs, 0:numenin)       :: xsopt        ! optimal cross section
  real(sgl), dimension(numchanxs,numsets*numenexp) :: xssamp       ! sampled random cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading cross sections
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                             :: flagdiscrete ! flag for discrete level cross sections
  logical                                             :: flagisomer   ! flag for isomeric production cross sections
  character(len=20), dimension(numchanxs)             :: xsfile       ! file with cross sections
  character(len=16), dimension(numchanxs)             :: reaction_string  ! reaction string
  integer, dimension(numchanxs)                       :: MT           ! MT number
  integer, dimension(numchanxs)                       :: MF           ! MF number
  integer, dimension(numchanxs)                       :: MTiso        ! isomer of MT number
  integer, dimension(nummt)                           :: MTnum        ! MT number of reaction channel
  integer                                             :: Nchanxs      ! total number of channels with cross sections
  integer, dimension(numchanxs)                       :: Nen          ! number of incident energies
  integer, dimension(numchanxs)                       :: Nen0         ! number of incident energies
  integer                                             :: Nenendf      ! number of incident energies
  integer                                             :: Nenendf0     ! number of incident energies
  integer, dimension(numchanxs)                       :: Nenlow       ! number of incident energies
  integer                                             :: Nennon       ! number of incident energies
  integer, dimension(numchanxs,numenin)               :: Sindex       ! index for cross section
  real(sgl), dimension(numchanxs, 0:numenin)          :: dE           ! uncertainty of incident energy
  real(sgl), dimension(numchanxs,0:numenin)           :: E            ! incident energy
  real(sgl), dimension(0:numenendf)                   :: Eendf        ! incident energy
  real(sgl), dimension(nummt)                         :: Efracmax     ! starting energy for maximum cross section
  real(sgl), dimension(numchanxs)                     :: E1mb         ! energy of 1 mb value
  real(sgl), dimension(0:numenin)                     :: Enon         ! incident energy
  real(sgl), dimension(nummt)                         :: fracmax      ! fraction of maximum cross sec. to be included
  real(sgl)                                           :: Qval         ! Q-value
  real(sgl), dimension(0:1,3,0:numenendf)             :: xsendf       ! cross section from TALYS
  real(sgl), dimension(0:1,numchanxs,0:numenin)       :: xseval       ! evaluated cross section
  real(sgl), dimension(numchanxs)                     :: xsmax        ! maximum cross section per channel
  real(sgl), dimension (numchanxs,numsets,0:numenexp) :: xsnon        ! nonelastic cross section
  real(sgl), dimension(0:numenin)                     :: xsnon0       ! non-elastic crosssection
  real(sgl), dimension(0:numtalys,numchanxs,0:numenS) :: xssave       ! cross section from TALYS
  real(sgl), dimension(0:1,numchanxs,0:numenin)       :: xstalys      ! cross section from TALYS
  real(sgl), dimension(0:1,numchanxs,0:numenin)       :: xstalys2     ! cross section from TALYS
  real(sgl), dimension(0:1,numchanxs,0:numenin)       :: xstalys3     ! cross section from TALYS
  real(sgl), dimension (numchanxs,numsets,0:numenexp) :: xsth         ! theoretical cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading residual production cross sections
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=20), dimension(numchanrp)             :: rpfile   ! name of file with residual production cross section
  integer, dimension(numchanrp)                       :: Arp      ! mass number of residual product
  integer, dimension(numchanrp)                       :: Lrp      ! isomeric level of residual product
  integer                                             :: Nchanrp  ! total number of channels with res. prod. c.s.
  integer, dimension(numchanrp)                       :: Nenrp    ! number of incident energies
  integer, dimension(numchanrp)                       :: Nenrp0   ! number of incident energies
  integer, dimension(numchanxs,numenin)               :: Srpindex ! index for residual production cross section
  integer, dimension(numchanrp)                       :: Zrp      ! charge number of residual product
  real(sgl), dimension(numchanrp,0:numenin)           :: Erp      ! incident energy
  real(sgl), dimension(0:numtalys,numchanrp,0:numenS) :: rpsave   ! residual production cross section from TALYS
  real(sgl), dimension(0:1,numchanrp,0:numenin)       :: rptalys  ! residual production cross section from TALYS
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading production cross sections
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=20), dimension(numchanprod)             :: prodfile   ! name of file with particle prod. cross sections
  integer                                               :: Nchanprod  ! number of channels with particle prod. c.s.
  integer, dimension(numchanprod)                       :: Nenprod    ! number of covariance energies
  integer, dimension(numchanprod)                       :: Nenprod0   ! number of covariance energies
  integer, dimension(numchanxs,numenin)                 :: Sprodindex ! index of particle production cross sections
  real(sgl), dimension(numchanprod,0:numenin)           :: Eprod      ! incident energy
  real(sgl), dimension(0:numtalys,numchanprod,0:numenS) :: prodsave   ! particle production cross section from TALY
  real(sgl), dimension(0:1,numchanprod,0:numenin)       :: prodtalys  ! particle production cross section from TALYS
  real(sgl), dimension(0:1,numchanprod,0:numenin)       :: Ytalys     ! yield from TALYS
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading gamma production cross sections
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=20), dimension(numchangam)             :: gamfile   ! name of file with gamma production cross sections
  integer                                              :: maxgam    ! maximum number of discrete levels with gamma c.s.
  integer                                              :: Nchangam  ! total number of channels with gamma prod. c.s.
  integer, dimension(numchangam)                       :: Nengam    ! number of incident energies
  integer, dimension(numchangam)                       :: Nengam0   ! number of incident energies
  integer, dimension(numchanxs,numenin)                :: Sgamindex ! index for gamma cross section
  real(sgl), dimension(numchangam,0:numenin)           :: Egam      ! incident energy
  real(sgl), dimension(0:numtalys,numchangam,0:numenS) :: gamsave   ! gamma cross section from TALYS
  real(sgl), dimension(0:1,numchangam,0:numenin)       :: gamtalys  ! gamma cross section from TALYS
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading spectra
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                      :: flagblock ! flag to block spectra, angle and gamma files
  character(len=20), dimension(numchansp)      :: spfile    ! name of file with emission spectra
  integer                                      :: Nchansp   ! total number of channels with emission spectra
  integer, dimension(numchansp)                :: Nensp     ! number of covariance energies
  real(sgl), dimension(numchansp,0:numen2)     :: Eout      ! emission energy
  real(sgl), dimension(numchansp)              :: Espec     ! energy
  real(sgl), dimension(0:1,numchansp,0:numen2) :: sptalys   ! emission spectrum from TALYS
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading angular distributions
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=21), dimension(numchanang)      :: angfile  ! name of file with angular distributions
  integer, dimension(numchanang)                :: Nang     ! number of angles
  integer                                       :: Nchanang ! total number of channels with angular distributions
  real(sgl), dimension(numchanang)              :: Eang     ! energy
  real(sgl), dimension(numchanang,0:numang)     :: angle    ! angle
  real(sgl), dimension(0:1,numchanang,0:numang) :: angtalys ! angular distribution
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading Legendre coefficients
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=21), dimension(numchanang)      :: legfile  ! name of file with Legendre coefficients
  integer                                       :: Nchanleg ! total number of channels with Legendre coefficients
  real(sgl), dimension(numchanang)              :: Eleg     ! energy grid for Legendre coefficients
  real(sgl), dimension(0:1,numchanang,0:numleg) :: leg0     ! Legendre coefficients from TALYS
  real(sgl), dimension(0:1,numchanang,0:numleg) :: legtalys ! Legendre coefficients from TALYS
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading integral data
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=15), dimension(numchanxs)    :: xseffchan  ! channel for effective cross section
  character(len=15), dimension(numchanxs)    :: xseffflux  ! flux for effective cross section
  integer                                    :: Nxseff     ! number of effective cross sections
  real(sgl), dimension(numchanxs)            :: xseffexp   ! experimental effective cross section
  real(sgl), dimension(0:1,numchanxs)        :: xseffrat   ! C/E of effective cross section
  real(sgl), dimension(0:numtalys,numchanxs) :: xseffsave  ! effective cross section from TALYS
  real(sgl), dimension(0:1,numchanxs)        :: xsefftalys ! effective cross section from TALYS
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading isotope production
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=20), dimension(numchanY)        :: Yfile        ! name of file with isotope production
  integer                                       :: NchanY       ! total number of channels with isotope production
  integer, dimension(numchanY)                  :: Ntime        ! number of time points
  real(sgl), dimension(0:1,numchanY,0:numtime)  :: acttalys     ! activity of produced isotope in MBq
  real(sgl), dimension(0:1,numchanY,0:numtime)  :: Nisoreltalys ! fraction of number of produced isotopes per elemen
  real(sgl), dimension(0:1,numchanY,0:numtime)  :: Nisotalys    ! number of isotopes produced after irradiation
  real(sgl), dimension(numchanY,0:numtime)      :: timeY        ! time point
  real(sgl), dimension(0:1,numchanY,0:numtime)  :: yieldtalys   ! yield of produced isotope in MBq/(mA.h)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading MACS
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                          :: flagD0exp
  real(sgl), dimension(0:numtalys) :: macssave  ! average Maxwellian averaged cross section
  real(sgl), dimension(0:1)        :: macstalys ! Maxwellian averaged cross section from TALYS
  real(sgl)                        :: ratmacs   ! Maxwellian rate
  real(sgl)                        :: Tmacs     ! Maxwellian temperature
  real(sgl)                        :: thmacs    ! theoretical MACS at 30 keV
  real(sgl)                        :: expmacs   ! experimental MACS at 30 keV
  real(sgl)                        :: dexpmacs  ! experimental unc. MACS at 30 keV
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl)                       :: Sw0         ! sum of weights
  real(sgl), dimension(0:numenin) :: Sweight     ! weight for TALYS run
  real(sgl), dimension(0:numenin) :: Sweightprev ! weight for previous  TALYS run
  real(sgl), dimension(0:numenin) :: Sweightsum  ! sum of weights for TALYS run
  real(sgl)                       :: Swp0        ! sum of weights
  real(sgl)                       :: Sws0        ! sum of weights
  real(sgl), dimension(0:numenin) :: weight      ! weight for experimental data sets
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for parameter covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=30), dimension(numpar)   :: parstring ! parameter string
  real(sgl), dimension(numpar,0:numenin) :: parav     ! average parameters
  real(sgl), dimension(numpar, numpar)   :: parcor    ! correlation of parameters
  real(sgl), dimension(numpar, numpar)   :: parcov    ! covariance matrix for parameters
  real(sgl), dimension(numpar)           :: pardif    ! difference of parameters
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for cross section covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                           :: flagband    ! flag to represent results as error band
  integer                                           :: average     ! cross sections 1: average,2: first run,3: optimum
  integer, dimension(numencov)                      :: Ecovindex   ! index for covariance energy on main energy grid
  integer, dimension(numchancov)                    :: MTcov       ! MT number with covariance data
  integer, dimension(numchancov)                    :: MTisocov    ! isomer of MT number with covariance data
  integer                                           :: Nchancovint ! number of channels with covariance data
  integer                                           :: Nencov      ! number of covariance energies
  real(sgl), dimension(numchanxs,numencov,numencov) :: Cmt         ! intra-channel correlation matrix for cross secti
  real(sgl)                                         :: coveps      ! limit for covariance
  real(sgl), dimension(0:numencov)                  :: Ecov        ! covariance energy grid
  real(sgl), dimension(numchanxs,numenin)           :: errmt       ! cross section uncertainty
  real(sgl), dimension(numchanxs,numencov)          :: errmtC      ! cross section uncertainty (for cov. energy grid)
  real(sgl)                                         :: Rlimit      ! limit for covariance calculation
  real(sgl), dimension(numchanxs,numencov,numencov) :: Rmt         ! intra-channel rel. cov. matrix for cross section
  real(sgl), dimension(numchanxs,numenin)           :: RmtD        ! diagonal of covariance matrix for cross sections
  real(sgl), dimension(numpar,numchanxs,numenS)     :: S           ! sensitivity matrix
  real(sgl), dimension(numpar,numchanxs,numenS)     :: Pearson     ! Pearson correlation
  real(sgl), dimension(numpar,numchanxs,numenS)     :: Pearson_enum ! variable for Pearson correlation
  real(sgl), dimension(numpar,numchanxs,numenS)     :: Pearson_denom_par ! variable for Pearson correlation
  real(sgl), dimension(numpar,numchanxs,numenS)     :: Pearson_denom_xs ! variable for Pearson correlation
  real(sgl), dimension(numpar,numchanxs,numenS)     :: Sdenom      ! denominator of sensitivity matrix
  real(sgl), dimension(numpar,numchanxs,numenS)     :: Senum       ! enumerator of sensitivity matrix
  real(sgl), dimension(numchanxs,numencov,numencov) :: Vmt         ! intra-channel covariance matrix for cross sectio
  real(sgl), dimension(numchanxs,numenin)           :: xsav        ! average cross section
  real(sgl), dimension(numchanxs,numencov)          :: xsavC       ! average cross section (for covariance energy grid)
  real(sgl), dimension(numchancov,numencov,numchancov,numencov) :: Rcov        ! relative covariance matrix for cross s
  real(sgl), dimension(numchancov,numencov,numchancov,numencov) :: Ccov        ! correlation matrix for cross sections
  real(sgl), dimension(numchancov,numencov,numchancov,numencov) :: Vcov        ! covariance matrix for cross sections
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for residual production cross section covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numchanrp,numenin)          :: errrp ! cross section uncertainty
  real(sgl), dimension(numchanrp,numenin)          :: rpav  ! average residual production cross section
  real(sgl), dimension(numchanrp,numenin,numencov) :: Rrp   ! covariance matrix for residual production cross sections
  real(sgl), dimension(numchanrp,numenin)          :: RrpD  ! diagonal of covariance matrix for cross sections
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for gamma production covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numchangam,numenin)           :: errgam ! cross section uncertainty
  real(sgl), dimension(numchangam,numenin)           :: gamav  ! average gamma production cross section
  real(sgl), dimension(numchangam,numencov,numencov) :: Rgam   ! covariance matrix for gamma production cross sections
  real(sgl), dimension(numchangam,numenin)           :: RgamD  ! diagonal of covariance matrix for cross sections
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for production cross section covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numchanprod,numenin)           :: errprod ! cross section uncertainty
  real(sgl), dimension(numchanprod,numenin)           :: prodav  ! average particle production cross section
  real(sgl), dimension(numchanprod,numencov,numencov) :: Rprod   ! covariance matrix for particle production cross sect
  real(sgl), dimension(numchanprod,numenin)           :: RprodD  ! diagonal of covariance matrix for cross sections
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for spectra covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numchansp,0:numen2)           :: errsp ! emission spectrum uncertainty
  real(sgl), dimension(numchansp,0:numen2, 0:numen2) :: Rsp   ! covariance matrix for emission spectra
  real(sgl), dimension(numchansp,0:numen2)           :: spav  ! average emission spectra
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for angular distribution covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numchanang,0:numang)          :: angav  ! average angular distributions
  real(sgl), dimension(numchanang,0:numang)          :: errang ! angular distributions uncertainty
  real(sgl), dimension(numchanang,0:numang,0:numang) :: Rang   ! covariance matrix for angular distributions
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for Legendre coefficient covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numchanang,0:numleg)                     :: errleg ! Legendre coefficients uncertainty
  real(sgl), dimension(numchanang,0:numleg)                     :: legav  ! average Legendre coefficients
  real(sgl), dimension(numchanang,0:numleg,numchanang,0:numleg) :: Rleg   ! covariance matrix for Legendre coefficients
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for integral data covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numchanxs)            :: erreff  ! effective cross section uncertainty
  real(sgl), dimension(numchanxs, numchanxs) :: Reff    ! covariance matrix for effective cross sections
  real(sgl), dimension(numchanxs)            :: ReffD   ! diagonal of covariance matrix for effective cross sections
  real(sgl), dimension(numchanxs)            :: xseffav ! average effective cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for isotope production covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl), dimension(numchanY, numtime) :: actav      ! average activity
  real(sgl), dimension(numchanY, numtime) :: erract     ! activity uncertainty
  real(sgl), dimension(numchanY, numtime) :: errNiso    ! error in diagonal of average number of isotopes produced afte
  real(sgl), dimension(numchanY, numtime) :: errNisorel ! error in fraction of number of produced isotopes per element
  real(sgl), dimension(numchanY, numtime) :: erryield   ! error in yield
  real(sgl), dimension(numchanY, numtime) :: Nisoav     ! average number of isotopes produced after irradiation
  real(sgl), dimension(numchanY, numtime) :: Nisorelav  ! average fraction of number of produced isotopes per element
  real(sgl), dimension(numchanY, numtime) :: RactD      ! diagonal of covariance matrix for activity
  real(sgl), dimension(numchanY, numtime) :: RNisoD     ! diagonal of average number of isotopes produced after irradia
  real(sgl), dimension(numchanY, numtime) :: RNisorelD  ! diagonal in fraction of number of produced isotopes per eleme
  real(sgl), dimension(numchanY, numtime) :: RyieldD    ! diagonal of yield covariance matrix
  real(sgl), dimension(numchanY, numtime) :: yieldav    ! average yield of produced isotope in MBq/(mA.h)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for MACS covariances
!-----------------------------------------------------------------------------------------------------------------------------------
!
  real(sgl)  :: errmacs ! Maxwellian averaged cross section uncertainty
  real(sgl)  :: macsav  ! average Maxwellian averaged cross section
  real(sgl)  :: RmacsD  ! diagonal of covariance matrix for Maxwellian averaged cross sections
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for sensitivity
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                :: flagreadsens ! flag to read sensitivities from tables
  logical                                :: flagsens     ! flag to use sensitivity in parameter weighting
  logical                                :: flagmorris   ! flag to use Morris screening
  integer                                :: Nmorris      ! number of samples for Morris screening
  integer                                :: iloop        ! counter for Morris screening
  real(sgl), dimension(numpar)           :: Sall         ! sensitivity for all channels
  real(sgl), dimension(numpar,numchanxs) :: Schan        ! sensitivity for channel
  real(sgl), dimension(numpar,numenin)   :: SE           ! sensitivity for energy
  real(sgl), dimension(numpar,numchanxs) :: xsdevav      ! average deviation of cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for reading experimental data
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                              :: flagedep       ! flag for energy dependence of uncertainties
  logical                                              :: flagscore      ! flag to use scores for inclusion or exclusion
  logical                                              :: flagerror      ! flag to include true exp. uncertainties
  logical                                              :: flagsort       ! flag to sort exp. data in increasing energy
  character(len=2), dimension(numchanxs,numsets)       :: quality        ! quality class
  character(len=9), dimension(numchanxs,numsets)       :: subentry       ! subentry
  character(len=15), dimension(numlib)                 :: lib            ! data library
  character(len=30), dimension(numchanxs,numsets)      :: auth           ! author
  character(len=132), dimension(nummt,numsets)         :: expexcfile     ! file with exp. data set to exclude for searc
  character(len=132), dimension(nummt,numsets)         :: expincfile     ! file with exp. data set to include for searc
  integer                                              :: expclass       ! quality class of exp. data (1, 2 or 3)
  integer, dimension(nummt)                            :: expexccount    ! counter for excluded experimental data sets
  integer, dimension(nummt)                            :: expinccount    ! counter for included experimental data sets
  integer, dimension(numchanxs,numsets)                :: isoexp         ! isomer of experimental data
  integer, dimension(nummt)                            :: isolib         ! isomer of MT number of library
  integer, dimension(numlib)                           :: lenlib         ! length of library name
  integer                                              :: maxexpsets     ! number of data sets before weight reduction
  integer, dimension(numchanxs*numsets)                :: MTexp          ! MT number of experimental data
  integer, dimension(nummt)                            :: MTlibk         ! MT number of library
  integer, dimension(numchanxs,numsets)                :: Nenexp         ! number of experimental data points
  integer                                              :: NMTexp         ! number of MT numbers for experimental data
  integer                                              :: NMTlib         ! number of MT numbers of library
  integer, dimension(numchanxs)                        :: Nsets          ! number of experimental data sets
  integer, dimension(numchanxs,numsets)                :: expyear        ! year
  real(sgl), dimension(numchanxs,numsets,numenexp)     :: dEexp          ! uncertainty of incident energy
  real(sgl), dimension(numchanxs,numsets,numenexp)     :: dxsexp         ! uncertainty of experimental cross section
  real(sgl), dimension(numchanxs,numsets)              :: Ebeg           ! first energy of experimental energy grid
  real(sgl), dimension(numchanxs,numsets)              :: Eend           ! last energy of experimental energy grid
  real(sgl), dimension(numchanxs,numsets,0:numenexp+1) :: Eexp           ! incident energy
  real(sgl), dimension(numenin)                        :: Ein            ! incident energy
  real(sgl), dimension(nummt)                          :: Esearch1       ! start energy of search
  real(sgl)                                            :: Esearch1all    ! global start energy of search
  real(sgl), dimension(nummt)                          :: Esearch2       ! end energy of search
  real(sgl)                                            :: Esearch2all    ! global end energy of search
  real(sgl), dimension(numchanxs,numsets)              :: Eweight        ! incident energy
  real(sgl), dimension(numlib)                         :: libweight      ! weight for library
  real(sgl), dimension(nummt)                          :: mtallexpweight ! weight for all channels per exp. channel
  real(sgl), dimension(nummt)                          :: mtalllibweight ! weight for all channels per librar
  real(sgl), dimension(nummt,numsets)                  :: mtexpweight    ! weight for experimental data set and channel
  real(sgl), dimension(nummt,numlib)                   :: mtlibweight    ! weight for library and channel
  real(sgl), dimension(nummt)                          :: mttalweight    ! weight per TALYS channel
  real(sgl), dimension(numchanxs,numsets,numenexp)     :: xsexp          ! experimental cross section
  real(sgl), dimension(numchanxs,numsets,numenexp)     :: xsexpsamp      ! sampled experimental cross section
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for experimental data files
!-----------------------------------------------------------------------------------------------------------------------------------
!
  character(len=132), dimension(numnuc,numchanexp):: expfilemul ! file with experimental data
  integer                                         :: aamax      ! maximal A value for which exp. data are included
  integer                                         :: aamin      ! minimal A value for which exp. data are included
  integer, dimension(numnuc)                      :: Nchanmul   ! number of reaction channels per nuclide
  integer                                         :: Nnuc       ! number of nuclides for which experimental data exists
  integer                                         :: zzmax      ! maximal Z value for which exp. data are included
  integer                                         :: zzmin      ! minimal Z value for which exp. data are included
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for fitting limits
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                               :: flagreadpar ! flag to read in parameter distribution
  character(len=25), dimension(numpar)  :: parfile     ! file with parameters
  integer, dimension(numpar, numhist)   :: nhist       ! number of histogram bins
  integer, dimension(numpar)            :: Nhistbin    ! number of histogram bins
  integer                               :: Nhistbinall ! global number of histogram bins
  integer, dimension(numpar, numhist)   :: nhistsam    ! number of hits in bin
  integer, dimension(numpar)            :: parclass    ! class of keyword
  integer, dimension(numpar)            :: readpar     ! flag to read parameter distribution
  real(sgl), dimension(numpar, numhist) :: hist        ! value in histogram bin
  real(sgl), dimension(numpar, numhist) :: histsam     ! value in histogram bin
  real(sgl), dimension(numpar)          :: parhigh     ! upper value of parameter uncertainty
  real(sgl), dimension(numpar)          :: parlow      ! lower value of parameter uncertainty
  real(sgl), dimension(numpar, numhist) :: Pbin        ! central value of histogram bin
  real(sgl), dimension(numpar, numhist) :: Pbot        ! lower value of histogram bin
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for GOF function
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical                                           :: flagdeltaE ! flag to weigh goodness-of-fit with dE of E grid
  logical                                           :: flagdexp   ! flag to include exp. uncertainty in C/E
  logical                                           :: flagerf    ! flag to use error function for Frms
  integer, dimension (numchanxs,numsets,0:numenexp) :: Eindex     ! index for experimental data
  integer                                           :: ichi2      ! parameter to set goodness-of-fit estimator
  integer                                           :: isearch    ! number of trial run
  integer                                           :: Nburn      ! number of burn-in runs
  integer, dimension(numenin)                       :: Nindex     ! index for experimental data
  integer                                           :: Ntalys     ! number of TALYS runs
  integer                                           :: outsearch  ! level of output of search (1 (minimal) - 5 (all))
  integer                                           :: searchmode ! search mode
  real(sgl), dimension(nummt)                       :: chi2max    ! maximal chi2 value per point taken into account
  real(sgl), dimension(nummt)                       :: Fmax       ! maximal F value per point taken into account
  real(sgl), dimension(0:numenin)                   :: Gall       ! GOF value
  real(sgl), dimension(0:numenin)                   :: Gall0      ! GOF value
  real(sgl), dimension(numchanxs)                   :: Gchannel   ! GOF value of channel
  real(sgl), dimension(0:numenin)                   :: Grun       ! GOF value of run
  real(sgl), dimension(0:numenin)                   :: Grun0      ! GOF value of run
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Variables for TMC
!-----------------------------------------------------------------------------------------------------------------------------------
!
  logical           :: flagacf      ! flag for creation of ENDF activation file
  logical           :: flagcovrand  ! flag for covariance ENDF-6 file after every random run
  logical           :: flagdefonly  ! flag for creation of ENDF default file only
  logical           :: flageaf      ! flag for creation of EAF activation file
  logical           :: flagfns      ! flag for variation of fission neutron spectrum for ENDF-6 general purpose file
  logical           :: flaggpf      ! flag for creation of ENDF general purpose file
  logical           :: flaglrf7     ! flag for resonance representation
  logical           :: flagmt       ! flag for creation of ENDF general purpose file without switch
  logical           :: flagnjoy     ! flag to run NJOY
  logical           :: flagnubar    ! flag for variation of nubar values for ENDF-6 general purpose file
  logical           :: flagprepro   ! flag to run PREPRO codes
  logical           :: flagproconly ! flag for processing of ENDF default file only
  logical           :: flagpurr     ! flag to run PURR in NJOY
  logical           :: flagres      ! flag for variation of resonance parameters for ENDF-6 general purpose file
  logical           :: flagrunfns   ! flag to produce fns in TASMAN with TANES
  logical           :: flagrunnubar ! flag to produce nubar in TASMAN with TAFIS
  logical           :: flagruntares ! flag to produce resonances in TASMAN with TARES
  logical           :: flags20      ! flag for creation of ENDF general purpose file with switch at 20 MeV
  logical           :: flags30      ! flag for creation of ENDF general purpose file with switch at 30 MeV
  logical           :: flags60      ! flag for creation of ENDF general purpose file with switch at 60 MeV
  logical           :: flagsdefault ! flag for creation of ENDF file with def. switch (30 MeV for neutrons, 0 for other particles)
  logical           :: flagselect   ! flag for varying only parts of an ENDF-6 data file
  character(len=16) :: tareslib     ! library name for adoption of resonance parameters (default or endfb8.1)
  character(len=16) :: tafislib     ! library name for adoption of nubar values (none[default], endfb8.1, jeff4.0, jendl5.0)
  character(len=16) :: taneslib     ! library name for adoption of FNS (none[default], endfb8.1, jeff4.0, jendl5.0)
  character(len=132):: background   ! library name for adoption of back ground cross sections in MF3
  integer           :: fislim       ! mass above which nuclide fissions
  integer           :: offset       ! offset for numbering or random files (TMC only)
  integer           :: tmcoffset    ! offset for starting creation of ENDF-6 files (TMC only)
end module A0_tasman_mod
! Copyright A.J. Koning 2026
