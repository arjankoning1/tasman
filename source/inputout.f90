subroutine inputout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write input parameters
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
! Variables for reading TASMAN input lines
!   inline            ! input line
!   nlines            ! number of input lines
! Variables for GOF function
!   flagdeltaE      ! flag to weigh goodness - of - fit with dE of energy grid
!   Nburn           ! number of burn - in runs
!   Ntalys          ! number of TALYS runs
!   outsearch       ! level of output of search (1 (minimal) - 5 (all))
!   ichi2           ! parameter to set goodness - of - fit estimator
! Variables for automatic parameter variation
!   flagparvary     ! flag to automatically vary parameters
!   maxbar          ! maximal fission barrier
!   minbar          ! minimal fission barrier
!   flagallpar      ! flag to vary all nuclear model parameters
! Variables for uncertainty
!   cwidth          ! scaling constant for parameter variation
!   fiscor0         ! adjustment parameter for fission uncertainties
!   flagsave        ! flag to save all input and output files of each Talys run
! Variables for reading experimental data
!   expclass        ! quality class of experimental data (1, 2 or 3)
!   flagedep        ! flag for energy dependence of uncertainties
!   flagscore       ! flag to use scores for inclusion or exclusion
!   flagerror       ! flag to include true experimental uncertainties
!   flagsort        ! flag to sort experimental data in increasing energies
!   maxexpsets      ! number of experimental data sets before weight reduction
! Variables for parameter variation
!   flagequisam     ! flag to use equal binning uniform sampling
!   flaggauss       ! flag to sample parameters from a Gaussian instead of uniform
!   flagkdburn      ! flag to use KD03 OMP function for burn - in
!   flagmcak        ! flag to use sampling from arbitrary distribution
!   flagmetro       ! flag for Metropolis sampling
!   flagsample      ! flag for sampling experimental data
! Variables for reading TALYS output
!   flagcross       ! flag for covariances or optimization of cross sections
!   flagexp         ! flag to optimize to experimental data
!   flaggamma       ! flag for covariances or optimization of discrete gamma ray transitions
!   flagintegral    ! flag for covariances or optimization of integral activation data
!   flagld          ! flag for optimization of level density parameters
!   flaglib         ! flag to optimize to nuclear data libraries
!   flagmacs        ! flag for covariances or optimization of Maxwellian averaged cross sections
!   flagmagnet      ! flag to push TALYS data to experimental or library data
!   flagprod        ! flag for covariances or optimization of particle production cross sections
!   flagproduct     ! flag for covariances or optimization of isotope production
!   flagpsf         ! flag for optimization of photon strength functions
!   flagresidual    ! flag for covariances or optimization of residual production cross sections
! Variables for weights
!   flagEweight     ! flag to use energy - dependent weights
!   flagweight      ! flag to use weights for random samples
! Variables for fitting limits
!   flagreadpar     ! flag to read in parameter distribution
! Variables for reading cross sections
!   flagdiscrete    ! flag for discrete level cross sections
!   flagisomer      ! flag for isomeric production cross sections
! Variables for cross section covariances
!   average         ! cross sections 1: average, 2: first run 3: optimum
!   flagband        ! flag to represent results as error bands
! Variables for reading gamma production cross sections
!   maxgam          ! maximum number of discrete levels with gamma cross section
! Variables for writing TALYS input files
!   Ehigh           ! energy cutoff of high energy runs
!   flagangle       ! flag for covariances or optimization of angular distributions
!   flagecis        ! flag to repeat ECIS calculations in a second run
!   flagglobal      ! flag to use global TALYS cross sections as prior
!   flaginponly     ! flag to produce only TALYS input files, no runs
!   flagleg         ! flag for covariances or optimization of Legendre coefficients
!   flagomponly     ! flag for covariances or optimization for basic optical model observables
!   flagspectra     ! flag for covariances or optimization of emission spectra
!   flagtmc         ! flag for Total Monte Carlo (ENDF - 6 data file at each run)
!   mode            ! TASMAN mode
!   Nhigh           ! number of high energy runs
! Variables for TMC
!   background      ! library name for adoption of back ground cross sections in MF3
!   flagacf         ! flag for creation of ENDF activation file
!   flagcovrand     ! flag for covariance ENDF - 6 file after every random run
!   flagdefonly     ! flag for creation of ENDF default file only
!   flageaf         ! flag for creation of EAF activation file
!   flagfns         ! flag for variation of fission neutron spectrum for ENDF - 6 general purpose file
!   flaggpf         ! flag for creation of ENDF general purpose file
!   flaglrf7        ! flag for resonance representation
!   flagmt          ! flag for creation of ENDF general purpose file without switch
!   flagnjoy        ! flag to run NJOY
!   flagnubar       ! flag for variation of nubar values for ENDF - 6 general purpose file
!   flagprepro      ! flag to run PREPRO codes
!   flagproconly    ! flag for processing of ENDF default file only
!   flagpurr        ! flag to run PURR in NJOY
!   flagres         ! flag for variation of resonance parameters for ENDF - 6 general purpose file
!   flagrunnubar    ! flag to produce nubar in TASMAN with TAFIS
!   flagruntares    ! flag to produce resonances in TASMAN with TAFISRES
!   flags20         ! flag for creation of ENDF general purpose file with switch at 20 MeV
!   flags30         ! flag for creation of ENDF general purpose file with switch at 30 MeV
!   flags60         ! flag for creation of ENDF general purpose file with switch at 60 MeV
!   flagsdefault    ! flag for creation of ENDF file with def. switch (30 MeV for neutrons, 0 for other particles)
!   flagselect      ! flag for varying only parts of an ENDF - 6 data file
!   offset          ! offset for numbering or random files (TMC only)
!   tafislib        ! library name for adoption of nubar values (none[default], endfb8.0, jeff3.3, jendl4.0, brc09)
!   taneslib        ! library name for adoption of FNS (none[default], endfb8.0, jeff3.3, jendl4.0, brc09)
!   tareslib        ! library name for adoption of resonance parameters (default or endfb8.0)
!   tmcoffset       ! offset for starting creation of ENDF - 6 files (TMC only)
! Variables for processing input
!   flagautoinc     ! flag to automatically in - or exclude important channels
!   Liso            ! number of the isomeric state
!   seed            ! seed for random number generator
!   tafisversion    ! version of TAFIS executable
!   talysversion    ! version of TALYS executable
!   tanesversion    ! version of TANES executable
!   taresversion    ! version of TARES executable
!   tefalversion    ! version of TEFAL executable
! Variables for sensitivity
!   flagreadsens    ! flag to read sensitivities from tables
!   flagsens        ! flag to use sensitivity in parameter weighting
! Variables for weights
!   gofmode         ! gof mode
!   weightpower     ! power for weight of random run
!
! *** Declaration of local data
!
  implicit none
  character(len=1) :: yesno     ! y or n function
  integer          :: i         ! counter
!
! ************************** User input file ***************************
!
  write(*, '(/" ########## USER INPUT ##########")')
  write(*, '(/" USER INPUT FILE"/)')
  do i = 1, nlines
    write(*, '(1x, a)') trim(inline(i))
  enddo
!
! ********* All possible input parameters including defaults ***********
!
  write(*, '(/" USER INPUT FILE + DEFAULTS"/)')
  write(*, '(" Keyword           Value   Variable", "     Explanation"/)')
!
! 1. Main keywords
!
  write(*, '(" #"/" # Main TASMAN keywords"/" #")')
  if (mode == 1) write(*, '(" #mode               ", i1, "     mode         mode for TASMAN: Monte Carlo covariances ")') mode
  if (mode == 2) write(*, '(" #mode               ", i1, "     mode         mode for TASMAN: sensitivity matrix with linear ", &
 &  "approximation ")') mode
  if (mode == 3) write(*, '(" #mode               ", i1, "     mode         mode for TASMAN: optimization with experimental ", &
 &  "data for one nucleus ")') mode
  if (mode == 4) write(*, '(" #mode               ", i1, "     mode         mode for TASMAN: optimization with experimental ", &
 &  "data for many nuclides ")') mode
  write(*, '(" #Ntalys        ", i6, "     Ntalys", "       number of TALYS runs ")') Ntalys
  write(*, '(" #Nburn         ", i6, "     Nburn ", "       number of burn-in runs ")') Nburn
  write(*, '(" #Nhigh         ", i6, "     Nhigh ", "       number of high energy runs ")') Nhigh
  write(*, '(" #Ehigh       ", f8.3, "     Ehigh ", "       energy cutoff of high energy runs ")') Ehigh
  write(*, '(" #Liso              ", i2, "     Liso  ", "       number of the isomeric state ")') Liso
  write(*, '(" #maxgam            ", i2, "     maxgam", "       maximum number of discrete levels with gamma cross section")') &
 &  maxgam
  write(*, '(" #minbar            ", i2, "     minbar", "       minimum number of fission barriers")') minbar
  write(*, '(" #maxbar            ", i2, "     maxbar", "       maximum number of fission barriers")') maxbar
  write(*, '(" #cross              ", a1, "     flagcross    flag for covariances or optimization of cross sections")') &
 &  yesno(flagcross)
  write(*, '(" #angle              ", a1, "     flagangle    flag for covariances or optimization of angular distributions")') &
 &  yesno(flagangle)
  write(*, '(" #legendre           ", a1, "     flagleg      flag for covariances or optimization of Legendre coefficients")') &
 &  yesno(flagleg)
  write(*, '(" #spectra            ", a1, "     flagspectra  flag for covariances or optimization of emission spectra")') &
 &  yesno(flagspectra)
  write(*, '(" #gamma              ", a1, "     flaggamma    flag for covariances or optimization of discrete gamma ray ", &
 &  "transitions")')  yesno(flaggamma)
  write(*, '(" #residual           ", a1, "     flagresidual flag for covariances or optimization of residual production ", &
 &  "cross sections")')  yesno(flagresidual)
  write(*, '(" #prod               ", a1, "     flagprod     flag for covariances or optimization of particle production ", &
 &  "cross sections")')  yesno(flagprod)
  write(*, '(" #production         ", a1, "     flagproduct  flag for covariances or optimization of isotope production")') &
 &  yesno(flagproduct)
  write(*, '(" #discrete           ", a1, "     flagdiscrete flag for covariances or optimization of discrete level ", &
 &  "cross sections")')  yesno(flagdiscrete)
  write(*, '(" #isomer             ", a1, "     flagisomer   flag for covariances or optimization of isomeric production", &
 &  " cross sections")')  yesno(flagisomer)
  write(*, '(" #ompinc             ", a1, "     flagompinc   flag to vary OMP parameters for incident channel")') &
 &  yesno(flagompinc)
  write(*, '(" #omponly            ", a1, "     flagomponly  flag for covariances or optimization for basic optical ", &
 &  "model observables")') yesno(flagomponly)
  write(*, '(" #integral           ", a1, "     flagintegral flag for covariances or optimization of integral", &
 &  " activation data")') yesno(flagintegral)
  write(*, '(" #macs               ", a1, "     flagmacs     flag for covariances or optimization of Maxwellian", &
 &  " averaged cross sections")') yesno(flagmacs)
  write(*, '(" #gamgam             ", a1, "     flaggamgam   flag for covariances or optimization of average", &
 &  " radiative width")') yesno(flaggamgam)
  write(*, '(" #psf                ", a1, "     flagpsf      flag for optimization of photon strength functions")') yesno(flagpsf)
  write(*, '(" #ld                 ", a1, "     flagld       flag for optimization of level density parameters")') yesno(flagld)
  write(*, '(" #sample             ", a1, "     flagsample   flag for sampling experimental data ")') yesno(flagsample)
  write(*, '(" #talysversion ", a, t27, " talysversion", " version of TALYS")') trim(talysversion)
  write(*, '(" #user   ", a, t27, "user         user for this calculation")') trim(user)
  write(*, '(" #source ", a, t27, "source       source for this calculation")') trim(source)
  write(*, '(" #format ", a, t27, "oformat      format for output")') trim(oformat)
  write(*, '(" #"/" # Keywords for mode ", i1, /" #")') mode
!
! 2. Keywords for uncertainties and covariances
!
  if (mode <= 2) then
    write(*, '(" #offset        ", i6, "     offset       offset for numbering or random files (TMC only) ")') offset
    write(*, '(" #tmcoffset     ", i6, "     tmcoffset    offset for starting creation of ENDF-6 files (TMC only)")') tmcoffset
    write(*, '(" #seed       ", i9, "     seed         seed for random number generator ")') seed
    write(*, '(" #cwidth  ", es12.5, "     cwidth       scaling constant for parameter variation ")') cwidth
    write(*, '(" #fiscor  ", es12.5, "     fiscor0      adjustment parameter for fission uncertainties")') fiscor0
    write(*, '(" #ecis               ", a1, "     flagecis     flag to repeat ECIS calculations in a second run ")') &
 &    yesno(flagecis)
    write(*, '(" #block              ", a1, "     flagblock    flag to block spectra, angle and gamma files")') yesno(flagblock)
    write(*, '(" #gauss              ", a1, "     flaggauss    flag to sample parameters from a Gaussian", &
 &    " instead of a uniform distribution")') yesno(flaggauss)
    write(*, '(" #save               ", a1, "     flagsave     flag to save all input and output files of each ", &
 &    "TALYS run ")') yesno(flagsave)
    write(*, '(" #tmc                ", a1, "     flagtmc      flag for Total Monte Carlo (ENDF-6 data file at each run")') &
 &    yesno(flagtmc)
    write(*, '(" #select             ", a1, "     flagselect   flag for varying only parts of an ENDF-6 data file")') &
 &    yesno(flagselect)
    write(*, '(" #covrand            ", a1, "     flagcovrand  flag for covariance ENDF-6 file after every random run")') &
 &    yesno(flagcovrand)
    write(*, '(" #gpf                ", a1, "     flaggpf      flag for creation of ENDF general purpose file ")') yesno(flaggpf)
    write(*, '(" #sdefault           ", a1, "     flagsdefault flag for creation of ENDF general purpose file", &
 &    " with default switch (30 MeV for neutrons, 0 for", " other particles)")') yesno(flagsdefault)
    write(*, '(" #defonly            ", a1, "     defonly      flag for creation of ENDF default file only")') yesno(flagdefonly)
    write(*, '(" #s20                ", a1, "     flags20      flag for creation of ENDF general purpose file", &
 &    " with switch at 20 MeV")') yesno(flags20)
    write(*, '(" #s30                ", a1, "     flags30      flag for creation of ENDF general purpose file", &
 &    " with switch at 30 MeV")') yesno(flags30)
    write(*, '(" #s60                ", a1, "     flags60      flag for creation of ENDF general purpose file", &
 &    " with switch at 60 MeV")') yesno(flags60)
    write(*, '(" #mt                 ", a1, "     flagmt       flag for creation of ENDF general purpose file", &
 &    " without switch")') yesno(flagmt)
    write(*, '(" #act                ", a1, "     flagacf      flag for creation of ENDF activation file ")') yesno(flagacf)
    write(*, '(" #eaf                ", a1, "     flageaf      flag for creation of EAF activation file ")') yesno(flageaf)
    write(*, '(" #prepro             ", a1, "     flagprepro   flag to run PREPRO codes")') yesno(flagprepro)
    write(*, '(" #njoy               ", a1, "     flagnjoy     flag to run NJOY")') yesno(flagnjoy)
    write(*, '(" #proconly           ", a1, "     flagproconly flag for processing of ENDF default file only")') yesno(flagproconly)
    write(*, '(" #purr               ", a1, "     flagpurr     flag to run PURR in NJOY")') yesno(flagpurr)
    write(*, '(" #resonance          ", a1, "     flagres      flag for variation of resonance parameters for", &
 &    " ENDF-6 general purpose file")') yesno(flagres)
    write(*, '(" #runtares           ", a1, "     flagruntares flag to produce resonances in TASMAN with TARES")') &
 &    yesno(flagruntares)
    write(*, '(" #nubar              ", a1, "     flagnubar    flag for variation of nubar values for", &
 &    " ENDF-6 general purpose file")') yesno(flagnubar)
    write(*, '(" #runnubar           ", a1, "     flagrunnubar flag to produce nubar in TASMAN with TAFIS")') yesno(flagrunnubar)
    write(*, '(" #fns                ", a1, "     flagfns      flag or variation of fission neutron spectrum for", &
 &    " ENDF-6 general purpose file ")') yesno(flagfns)
    write(*, '(" #tefalversion ", a, t27, " tefalversion version of TEFAL ")') trim(tefalversion)
    write(*, '(" #taresversion ", a, t27, " taresversion version of TARES ")') trim(taresversion)
    write(*, '(" #tafisversion ", a, t27, " tafisversion version of TAFIS ")') trim(tafisversion)
    write(*, '(" #tanesversion ", a, t27, " tanesversion version of TANES ")') trim(tanesversion)
    write(*, '(" #tareslib ", a, t27, " tareslib     library name for adoption of resonance parameters")') trim(tareslib)
    write(*, '(" #background ", a, t27, " background   library name for adoption of background cross sections in MF3")') &
 &    trim(background)
    write(*, '(" #lrf7               ", a1, "     flaglrf7     flag for resonance representation")') yesno(flaglrf7)
    write(*, '(" #tafislib ", a, t27, " tafislib     library name for adoption of nubar values")') trim(tafislib)
    write(*, '(" #taneslib ", a, t27, " taneslib     library name for adoption of fission neutron spectrum")') trim(taneslib)
    write(*, '(" #mcak               ", a1, "     flagmcak     flag to use sampling from arbitrary distribution")') yesno(flagmcak)
    write(*, '(" #readpar            ", a1, "     flagreadpar  flag to read in parameter distribution")') yesno(flagreadpar)
    write(*, '(" #parameters         ", a1, "     flagparvary  flag to vary tabulated nuclear model parameters")') &
 &    yesno(flagparvary)
    write(*, '(" #allpar             ", a1, "     flagallpar   flag to vary all nuclear model parameters")') yesno(flagallpar)
    write(*, '(" #weight             ", a1, "     flagweight   flag to use weights for random samples")') yesno(flagweight)
    write(*, '(" #eweight            ", a1, "     flageweight  flag to use energy-dependent weights ")') yesno(flagEweight)
    write(*, '(" #metropolis         ", a1, "     flagmetro    flag for Metropolis sampling")') yesno(flagmetro)
    write(*, '(" #kdburn             ", a1, "     flagkdburn   flag to use KD03 OMP function for burn-in")') &
 &    yesno(flagkdburn)
    write(*, '(" #magnet             ", a1, "     flagmagnet   flag to push TALYS data to experimental or library data")') &
 &    yesno(flagmagnet)
    write(*, '(" #global             ", a1, "     flagglobal  ", " flag to use global TALYS cross sections as prior")') &
 &    yesno(flagglobal)
    write(*, '(" #autoinclude        ", a1, "     flagautoinc  flag to automatically in- or exclude important channels")') &
 &    yesno(flagautoinc)
    write(*, '(" #edependent         ", a1, "     flagedep     flag for energy dependence of uncertainties")') yesno(flagedep)
    write(*, '(" #equisample         ", a1, "     flagequisam  flag to use equal binning uniform sampling")') yesno(flagequisam)
    write(*, '(" #gofmode            ", i1, "     gofmode      1: exp(-0.5*G) 2: exp(-0.5*G/G0) ", &
 &    "3: exp(-(G/G0)**weightpower GOF weight")') gofmode
    write(*, '(" #maxexpsets   ", i7, "     maxexpsets   number of experimental data sets before weight reduction")') maxexpsets
    write(*, '(" #average            ", i1, "     average      average cross sections 1: average, 2: first run", &
 &    " 3: optimum")') average
    write(*, '(" #inponly            ", a1, "     flaginponly  flag to produce only TALYS input files, no runs")') &
 &    yesno(flaginponly)
  else
!
! 3. Keywords for optimization
!
    write(*, '(" #chi2               ", i1, "     ichi2        parameter to set goodness-of-fit estimator ")') ichi2
    write(*, '(" #weightpower      ", f7.3, " weightpower  power for weight of random run")') weightpower
    write(*, '(" #readsens           ", a1, "     flagreadsens flag to read sensitivities from tables")') yesno(flagreadsens)
    write(*, '(" #sens               ", a1, "     flagsens     flag to use sensitivity in parameter weighting")') yesno(flagsens)
    write(*, '(" #deltaE             ", a1, "     flagdeltaE   flag to weigh goodness-of-fit with dE of energy grid")') &
 &    yesno(flagdeltaE)
    write(*, '(" #dexp               ", a1, "     flagdexp     flag to include exp. uncertainty in C/E")') yesno(flagdexp)
    write(*, '(" #band               ", a1, "     flagband     flag to represent results as error bands ")') yesno(flagband)
    write(*, '(" #E1vary             ", a1, "     flagE1vary   flag to automatically vary keywords for E1 radiation")') &
 &    yesno(flagE1vary)
    write(*, '(" #M1vary             ", a1, "     flagM1vary   flag to automatically vary keywords for M1 radiation")') &
 &    yesno(flagM1vary)
    write(*, '(" #score              ", a1, "     flagscore    flag to use scores for inclusion or exclusion")') yesno(flagscore)
    write(*, '(" #error              ", a1, "     flagerror    flag to include true experimental uncertainties")') yesno(flagerror)
    write(*, '(" #expclass         ", i3, "     expclass     quality class of experimental data (1, 2 or 3)")') expclass
    write(*, '(" #sort               ", a1, "     flagsort     flag to sort experimental data in increasing energies in the ", &
 &    "output")') yesno(flagsort)
    write(*, '(" #searchexp          ", a1, "     flagexp      flag to optimize to experimental data")') yesno(flagexp)
    write(*, '(" #searchlib          ", a1, "     flaglib      flag to optimize to nuclear data libraries")') yesno(flaglib)
    write(*, '(" #outsearch        ", i3, "     outsearch    level of output of search (1 (minimal) - 5 (all))")') outsearch
  endif
  return
end subroutine inputout
! Copyright A.J. Koning 2021
