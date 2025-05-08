subroutine inputprocess
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Process input parameters
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
!   numchanang        ! maximum number of channels with angular distributions
!   numencov          ! maximum number of covariance energies
!   numenin           ! maximum number of incident energies
!   numlib            ! maximum number of data libraries
!   nummt             ! maximum number of MT numbers
!   numsets           ! maximum number of experimental data sets
!   numtalys          ! maximum number of TALYS runs
! Variables for GOF function
!   chi2max           ! maximal chi2 value per point taken into account
!   Fmax              ! maximal F value per point taken into account
!   Nburn             ! number of burn - in runs
!   Ntalys            ! number of TALYS runs
! Variables for cross section covariances
!   Ecov              ! covariance energy grid
!   Ecovindex         ! index for covariance energy on main energy grid
!   Nencov            ! number of covariance energies
! Variables for reading TALYS output
!   flagexp           ! flag to optimize to experimental data
!   flaglib           ! flag to optimize to nuclear data libraries
!   flagtal           ! flag to optimize to TALYS
!   Ntalbeg           ! counter for first TALYS calculation
!   flagcross         ! flag for covariances or optimization of cross sections
!   flaggamma         ! flag for covariances or optimization of discrete gamma ray transitions
!   flagld            ! flag for covariances or optimization of level density parameters
!   flagintegral      ! flag for covariances or optimization of integral activation data
!   flagmacs          ! flag for covariances or optimization of Maxwellian averaged cross sections
!   flagprod          ! flag for covariances or optimization of particle production cross sections
!   flagproduct       ! flag for covariances or optimization of isotope production
!   flagpsf           ! flag for optimization of photon strength functions
!   flagresidual      ! flag for covariances or optimization of residual production cross sections
! Constants
!   particle          ! type of particle
! Variables for path names
!   binpath           ! directory containing files to be read
!   librariespath     ! directory containing files to be read
!   tasmanpath        ! directory containing files to be read
! Variables for parameter variation
!   Npar              ! number of parameters
! Variables for reading cross sections
!   Efracmax          ! starting energy for maximum cross section
!   fracmax           ! fraction of maximum cross section to be included
! Variables for writing TALYS input files
!   flagglobal        ! flag to use global TALYS cross sections as prior
!   flagtmc           ! flag for Total Monte Carlo (ENDF - 6 data file at each run)
!   mode              ! TASMAN mode
!   Nhigh             ! number of high energy runs
! Variables for reading experimental data
!   Ein               ! incident energy
!   Esearch1          ! start energy of search
!   Esearch1all       ! global start energy of search
!   Esearch2          ! end energy of search
!   Esearch2all       ! global end energy of search
!   expinccount       ! counter for included experimental data sets
!   expincfile        ! file with experimental data set to include for search
!   libweight         ! weight for library
!   mtallexpweight    ! weight for all channels per experimental channel
!   mtalllibweight    ! weight for all channels per librar
!   mtexpweight       ! weight for experimental data set and channel
!   mtlibweight       ! weight for library and channel
!   mttalweight       ! weight per TALYS channel
! Variables for TMC
!   flagrunfns        ! flag to produce fns in TASMAN with TANES
!   flagrunnubar      ! flag to produce nubar in TASMAN with TAFIS
!   flagruntares      ! flag to produce resonances in TASMAN with TARES
!   flags30           ! flag for creation of ENDF general purpose file with switch at 30 MeV
!   flagsdefault      ! flag for creation of ENDF file with def. switch (30 MeV for neutrons, 0 for other particles)
! Variables for writing TALYS input
!   flagangle         ! flag for covariances or optimization of angular distributions
!   flagleg           ! flag for covariances or optimization of Legendre coefficients
!   flagspectra       ! flag for covariances or optimization of emission spectra
! Variables for processing input
!   allexp            ! flag to include all experimental channels
!   allexpweight      ! weight for all experimental channels
!   alllib            ! flag to include all library channels
!   alllibweight      ! weight for all library channels
!   alltal            ! flag to include all TALYS channels
!   alltalweight      ! weight for all TALYS channels
!   anginclude        ! flag to include energy point for angular grid
!   Atarget           ! mass number of target nucleus
!   auto              ! flag to in - or exclude channel
!   chi2maxall        ! global maximal chi2 value per point taken into account
!   dsmooth           ! parameter for energy smoothing
!   dsmoothall        ! global parameter for energy smoothing
!   Efracmaxall       ! global starting energy for maximum cross section
!   energyfile        ! file with incident energies
!   errlim            ! lower limit for experimental error
!   errlimall         ! global lower limit for experimental error
!   exppath           ! path for experimental data
!   flagautoinc       ! flag to automatically in - or exclude important channels
!   Fmaxall           ! global maximal F value per point taken into account
!   fracmaxall        ! global fraction of maximum cross section to be included
!   isochar           ! symbol for isomer
!   k0                ! index of incident particle
!   liballmt          ! flag for assignment of all MTs per library
!   libinclude        ! flag to include data library for search
!   Liso              ! number of the isomeric state
!   Ltarget           ! excited level of target
!   mtallexp          ! weight for all experiments per channel
!   mtalllib          ! weight for all libraries per channel
!   mtlib             ! flag for library per channel
!   mttal             ! weight for TALYS per channel
!   noexp             ! flag to include no experimental channels
!   nolib             ! flag to include no library channels
!   notal             ! flag to include no TALYS channels
!   numE              ! number of incident energies
!   ptype0            ! type of incident particle
!   seed              ! seed for random number generator
!   Starget           ! symbol of target nucleus
!   tafis             ! TAFIS executable
!   tafisversion      ! version of TAFIS executable
!   talys             ! TALYS executable
!   talysversion      ! version of TALYS executable
!   tanes             ! TANES executable
!   tanesversion      ! version of TANES executable
!   tares             ! TAREL executable
!   taresversion      ! version of TARES executable
!   tefal             ! TEFAL executable
!   tefalversion      ! version of TEFAL executable
!   Tweight           ! weight of TALYS vs experiment
!   Tweightall        ! global weight of TALYS vs experiment
! Error handling
!   range_integer_error    ! Test if integer variable is out of range
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: numsub=200000   ! number of subentries
  logical            :: lexist          ! logical to determine existence
  character(len=1)   :: ch              ! character
  character(len=4)   :: massstring      ! string for mass number
  character(len=132) :: exppath1        ! 
  character(len=132) :: exppath2        ! 
  character(len=132) :: aut(numsub)     ! author
  character(len=132) :: efile           ! file with incident energies
  character(len=132) :: sub             ! file with experimental data set to include for search
  character(len=132) :: suben(numsub)   ! subentry
  integer            :: i               ! counter
  integer            :: icov            ! index for covariances
  integer            :: icovprev        ! index for covariances
  integer            :: ii              ! counter
  integer            :: ilib            ! counter for library
  integer            :: imt             ! MT counter
  integer            :: istat           ! logical for file access
  integer            :: j               ! counter
  integer            :: nc              ! counter
  integer            :: nin             ! counter for incident energy
  integer            :: nsub            ! number of subentries
  real(sgl)          :: Edist           ! help variable
  real(sgl)          :: Edisttry        ! help variable
  real(sgl)          :: eee             ! incident energy
  real(sgl)          :: fac             ! factor
!
! ********************** Process input parameters **********************
!
  do i = 0, 6
    if (ptype0 == particle(i)) k0 = i
  enddo
  if (Ntalys ==  -1) then
    if (mode == 1) Ntalys = 10
    if (mode >= 3) Ntalys = numtalys
  endif
  if (mode == 2) Ntalys = Npar
  if (Nburn ==  -1) Nburn = min(Ntalys, numtalys)
  if (Nhigh ==  -1) Nhigh = min(Ntalys, numtalys)
  if (Ltarget /= 0 .and. Liso == 0) Liso = 1
  if (Liso == 0) isochar = ' '
  if (Liso == 1) isochar = 'm'
  if (Liso > 1) isochar = 'n'
  if (flagsdefault .and. k0 == 1) flags30 = .false.
  if (mode > 2) flaggauss = .false.
  do i=1,nummt
    MTreac(i)(2:2)=ptype0
  enddo
  massstring = '    '
  write(massstring,'(i3)') Atarget
  targetnuclide = trim(Starget)// trim(adjustl(massstring)) // trim(isochar)
!
! ********************** Avoid inconsistent flags **********************
!
  if (flagld .or. flagpsf) then
    flagcross = .false.
    flagresidual = .false.
    flagprod = .false.
    flaggamma = .false.
    flagspectra = .false.
    flagangle = .false.
    flagleg = .false.
    flagintegral = .false.
    flagproduct = .false.
    flagmacs = .false.
  endif
  if (k0 == 0) flagompinc = .true.
  if (flagmacs) flagcross = .false.
  if (flaggamgam) flagcross = .false.
  if (flagld .and. flagpsf) then
    write(*,'(" TASMAN-error: Only ld or psf can be equal to y")')
    stop
  endif
!
! Set TALYS version to be used and test availability of other codes
!
  inquire (file = talys, exist = lexist)
  if ( .not. lexist) then
    write(*, '(" TASMAN-error: TALYS can not be found: ", a)') trim(talys)
    stop
  endif
  if (flagtmc) then
    inquire (file = tefal, exist = lexist)
    if ( .not. lexist) then
      write(*, '(" TASMAN-error: TEFAL can not be found: ", a)') trim(tefal)
      stop
    endif
    if (flagruntares) then
      inquire (file = tares, exist = lexist)
      if ( .not. lexist) then
        write(*, '(" TASMAN-error: TARES can not be found: ", a)') trim(tares)
        stop
      endif
    endif
    if (flagrunnubar) then
      inquire (file = tafis, exist = lexist)
      if ( .not. lexist) then
        write(*, '(" TASMAN-error: TAFIS can not be found: ", a)') trim(tafis)
        stop
      endif
    endif
    if (flagrunfns) then
      inquire (file = tanes, exist = lexist)
      if ( .not. lexist) then
        write(*, '(" TASMAN-error: TANES can not be found: ", a)') trim(tanes)
        stop
      endif
    endif
    if (flaggpf .and. flagnjoy) then
      inquire (file = trim(binpath)//'autonjoy', exist = lexist)
      if (.not.lexist) then
        write(*,'(" TASMAN-error: autonjoy can not be found: ",a)') trim(binpath)//'autonjoy'
        stop
      endif
    endif
    if (flageaf .and. flagprepro) then
      inquire (file=trim(binpath)//'autoprepro', exist=lexist)
      if (.not.lexist) then
        write(*,'(" TASMAN-error: autoprepro can not be found: " ,a)') trim(binpath)//'autoprepro'
        stop
      endif
    endif
    if (flagselect) then
      inquire (file=trim(binpath)//'select', exist=lexist)
      if (.not.lexist) then
        write(*,'(" TASMAN-error: select can not be found: " ,a)') trim(binpath)//'select'
        stop
      endif
    endif
  endif
!
! Initialize random number generator
!
! init_genrand: subroutine for initialization of random numbers
!
  call init_genrand(seed)
!
! ********************** Read energies from energyfile *****************
!
  if (numE == 0) then
    efile = energyfile
    if (mode == 1) then
      inquire (file = 'energies.endf', exist = lexist)
      if (lexist) efile = 'energies.endf'
    endif
    nin = 0
    open (unit = 2, file = efile, status = 'old')
    do
      read(2, * , iostat = istat) eee
      if (istat == -1) exit
      if (istat /= 0) call read_error(efile, istat)
      if (eee /= 0.) then
        nin = nin + 1
        Ein(nin) = eee
        numE = nin
      endif
      call range_integer_error('energies.endf', nin, 1, numenin)
    enddo
  endif
!
! ********************** Set covariance energy grid ********************
!
  if (mode <= 2) then
    ii = 0
    icov = 0
    icovprev = 0
    do nc = 1, numencov
      Edist = 1000.
      do nin = 1, numE
        Edisttry = abs(Ecov(nc) - Ein(nin))
        if (Edisttry < Edist) then
          Edist = Edisttry
          icov = nin
        endif
      enddo
      if (icov /= icovprev) then
        ii = ii + 1
        Ecovindex(ii) = icov
        icovprev = icov
      endif
    enddo
    Nencov = ii
  endif
!
! ************* Determine energy grid for angles and spectra ***********
!
  if (numE <= numchanang) then
    do nin = 1, numE
      anginclude(nin) = .true.
    enddo
  else
    anginclude(1) = .true.
    anginclude(numE) = .true.
    fac = real(numE - 1) / (numchanang - 1)
    do i = 1, numchanang - 2
      nin = min(int(1. + fac * i), numchanang)
      anginclude(nin) = .true.
    enddo
  endif
!
! *** Include or exclude library or experimental data sets for search **
!
! In- or exclude data libraries
!
  flaglib = .false.
  i = -1
  if (alllib) i = 1
  if (nolib) i = 0
  do imt = 1, nummt
    if (mtalllib(imt) ==  -1) mtalllib(imt) = i
    do ilib = 1, numlib
      if (mtlib(imt, ilib) ==  -1) then
        mtlib(imt, ilib) = mtalllib(imt)
        if (liballmt(ilib) /=  -1) mtlib(imt, ilib) = liballmt(ilib)
      endif
    enddo
  enddo
  do imt = 1, nummt
    if (mtalllibweight(imt) == 1.) mtalllibweight(imt) = alllibweight
    do ilib = 1, numlib
       if (mtlibweight(imt, ilib) == 1.) mtlibweight(imt, ilib) = mtalllibweight(imt)
    enddo
  enddo
Loop1:  do ilib = 1, numlib
    do imt = 1, nummt
      if (mtlib(imt, ilib) == 1) then
        libinclude(ilib) = 1
        flaglib = .true.
        cycle Loop1
      endif
    enddo
  enddo Loop1
  do ilib = 1, numlib
    if (libweight(ilib) == 1.) libweight(ilib) = alllibweight
    do imt = 1, nummt
      if (mtlibweight(imt, ilib) == 1.) mtlibweight(imt, ilib) = libweight(ilib)
    enddo
  enddo
!
! In- or exclude experimental data
!
  flagexp = .false.
  if (allexp .and. noexp) then
    write(*, '(" TASMAN-error: do not include and exclude all ", " experimental data")')
    stop
  endif
  i = -1
  if (noexp) i = 0
  if (allexp) i = 1
  do imt = 1, nummt
    if (mtallexp(imt) ==  -1) mtallexp(imt) = i
  enddo
  do imt = 1, nummt
    if (mtallexpweight(imt) == 1.) mtallexpweight(imt) = allexpweight
    do i = 1, numsets
       if (mtexpweight(imt, i) == 1.) mtexpweight(imt, i) = mtallexpweight(imt)
    enddo
    if (mtallexp(imt) == 1 .or. expinccount(imt) > 0) flagexp = .true.
  enddo
!
! In- or exclude central values from TALYS
!
  flagtal = .false.
  if (alltal .and. notal) then
    write(*, '(" TASMAN-error: do not include and exclude all ", " TALYS data")')
    stop
  endif
  i = -1
  if (notal) i = 0
  if (alltal) i = 1
  do imt = 1, nummt
    if (mttal(imt) ==  -1) mttal(imt) = i
  enddo
  do imt = 1, nummt
    if (mttalweight(imt) == 1.) mttalweight(imt) = alltalweight
    if (mttal(imt) == 1) flagtal = .true.
  enddo
  if (alltal .and. flagglobal) then
    Ntalbeg = -1
  else
    Ntalbeg = 0
  endif
!
! Automatically in- or exclude important channels
!
!   channels
!
  if (flagautoinc) then
    do imt = 1, nummt
      auto(imt) = 0
    enddo
    do imt = 1, 4
      auto(imt) = 1
    enddo
    do imt = 16, 18
      auto(imt) = 1
    enddo
    auto(22) = 1
    auto(28) = 1
    do imt = 51, 54
      auto(imt) = 1
    enddo
    auto(102) = 1
    auto(103) = 1
    auto(107) = 1
  else
    do imt = 1, nummt
      auto(imt) = 1
    enddo
  endif
!
! ********************* Set path for experimental data *****************
!
  if (flagexp) then
    massstring = '000 '
    write(massstring(1:3), '(i3.3)') Atarget
    if (Liso > 0) massstring(4:4) = isochar
    exppath1 = trim(exforpath)//ptype0//'/'// trim(Starget)//trim(massstring)//'/xs/'
    exppath2 = trim(librariespath)//ptype0//'/'// trim(Starget)//trim(massstring)//'/exfor/xs/'
    inquire (file = trim(exppath1)//'xslist', exist = lexist)
    if (lexist) then
      exppath = exppath1
    else
      exppath = exppath2
    endif
  endif
!
! ****** Process weights for TALYS and experimental data ***************
!
  if (Esearch1all ==  -1.) Esearch1all = 27. * exp( - 0.7 * sqrt(real(Atarget)))
  do imt = 1, nummt
    if (Tweight(imt) ==  -1.) Tweight(imt) = Tweightall
    if (dsmooth(imt) ==  -1.) dsmooth(imt) = dsmoothall
    if (Fmax(imt) ==  -1.) Fmax(imt) = Fmaxall
    if (chi2max(imt) ==  -1.) chi2max(imt) = chi2maxall
    if (errlim(imt) ==  -1.) errlim(imt) = errlimall
    if (Efracmax(imt) ==  -1.) Efracmax(imt) = Efracmaxall
    if (fracmax(imt) ==  -1.) fracmax(imt) = fracmaxall
    if (Esearch1(imt) ==  -1.) Esearch1(imt) = Esearch1all
    if (Esearch2(imt) ==  -1.) Esearch2(imt) = Esearch2all
  enddo
!
! ****** Connect subentry number to experimental data filenames ********
!
  open (unit = 3, file = trim(tasmanpath)//'misc/subentries', status = 'old')
  i = 1
  do
    read(3, '(2a)', iostat = istat) suben(i), aut(i)
    if (istat == -1) exit
    if (istat /= 0) call read_error(trim(tasmanpath)//'misc/subentries', istat)
    i = i + 1
  enddo
  nsub = i - 1
  close (unit = 3)
Loop2:  do imt = 1, nummt
    do j = 1, expinccount(imt)
      ch = expincfile(imt, j)(1:1)
      if (ch >= '0' .and. ch <= '9') then
        sub = expincfile(imt, j)
        do i = 1, nsub
          if (trim(sub) == trim(suben(i))) then
            expincfile(imt, j) = aut(i)(1:40)
            cycle Loop2
          endif
        enddo
      endif
    enddo
  enddo Loop2
  return
end subroutine inputprocess
! Copyright A.J. Koning 2021
