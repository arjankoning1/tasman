subroutine koekel(N, x, xmin, xmax, xopt, Fopt, mode)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Optimization of parameters
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  logical   :: flagall                ! logical for variation of all parameters simultaneously
  logical   :: flaghistory            ! logical for using information from previous runs
  logical   :: flagquasi              ! logical for quasi-random search
  logical   :: flagsave               ! flag to save all input and output files of each  Talys run
  logical   :: flagsens               ! flag to use sensitivity in parameter weighting
  logical   :: flagsobol              ! logical for quasi-random search with Sobol sequence
  logical   :: flagwrite              ! logical for extensive output
  integer   :: mode                   ! mode for TASMAN
  integer   :: N                      ! number of parameters
  integer   :: Neps                   ! number of times tolerance is achieved before termination
  integer   :: NR                     ! number of runs with random parameters to determine  temperature
  integer   :: NS                     ! number of loops through all parameters before step  adjustment
  integer   :: NT                     ! number of loops before temperature reduction
  real(sgl) :: Cglob                  ! adjustment factor for step length for all parameters
  real(sgl) :: eps                    ! convergence criteria
  real(sgl) :: Fopt                   ! optimal function value
  real(sgl) :: RT                     ! factor for temperature reduction
  real(sgl) :: T                      ! initial temperature
  real(sgl) :: x(N)                   ! initial parameter vector
  real(sgl) :: xmax(N)                ! vector of upper boundary of parameter values
  real(sgl) :: xmin(N)                ! vector of lower boundary of parameter values
  real(sgl) :: xopt(N)                ! optimal parameters
!
! Optimization
!
! koekelinput : subroutine to set defaults and read input for KOEKEL optimization program
! koekelsearch: subroutine for simulated annealing optimization
! nelder      : subroutine for Nelder-Mead optimization
! rosenbrock  : subroutine for Rosenbrock optimization
! csendes     : subroutine for Csendes optimization
!
  call koekelinput(flagwrite, flagsens, flagall, flagquasi, flagsobol, flaghistory, flagsave, NR, NS, NT, Neps, RT, T, Cglob, eps)
!
! Check for errors in subroutine call.
!
! All number of parameters needs to fall within a certain range.
! If not, KOEKEL will stop.
!
  if (mode < 1 .or. mode > 4) then
    write(*,'(" KOEKEL-error: 1 <= mode <= 4")')
    stop
  endif
  if (N < 1 ) then
    write(*, '(" KOEKEL-error:  N >= 1")')
    stop
  endif
  if (mode == 4 .and. N > 100) then
    write(*,'(" KOEKEL-error: 1 <= N <= 100 for mode 4 or create file koekel.inp with line mode 1, 2 or 3")')
    stop
  endif
  if (mode == 1) call koekelsearch(N, x, xmin, xmax, xopt, Fopt, &
      flagwrite, flagsens, flagall, flagquasi, flagsobol, flaghistory, flagsave, NR, NS, NT, Neps, RT, T, Cglob, eps)
  if (mode == 2) call nelder(N, x, xmin, xmax, xopt, Fopt, eps)
  if (mode == 3) call rosenbrock(N, x, xmin, xmax, xopt)
  if (mode == 4) then
    call fcn(N, x, Fopt)
    call csendes(N, x, xmin, xmax, xopt, eps, NS)
  endif
!
! Final function evaluation with optimal values
!
  call fcn(N, xopt, Fopt)
  return
end subroutine koekel
! Copyright A.J. Koning 2021
subroutine koekelinput(flagwrite, flagsens, flagall, flagquasi, &
 &  flagsobol, flaghistory, flagsave, NR, NS, NT, Neps, RT, T, Cglob, eps)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Set defaults and read input for KOEKEL optimization program
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer, parameter :: numkey=16         ! number of keywords
  integer, parameter :: numlin=100        ! maximum number of input lines
  logical            :: flagall           ! logical for variation of all parameters simultaneously
  logical            :: flaghistory       ! logical for using information from previous runs
  logical            :: flagquasi         ! logical for quasi-random search
  logical            :: flagsave          ! flag to save all input and output files of each  Talys run
  logical            :: flagsens          ! flag to use sensitivity in parameter weighting
  logical            :: flagsobol         ! logical for quasi-random search with Sobol sequence
  logical            :: flagwrite         ! logical for extensive output
  character(len=1)   :: ch                ! character
  character(len=1)   :: chprev            ! character
  character(len=132) :: inline(numlin)    ! input line
  character(len=132) :: key               ! keyword
  character(len=132) :: keyword(numkey)   ! keyword
  character(len=132) :: param             ! parameter
  character(len=132) :: value             ! value or string
  integer            :: i                 ! counter
  integer            :: inp               ! unit for input
  integer            :: k                 ! Legendre order
  integer            :: kbeg              ! index to mark begin of word
  integer            :: kend              ! index to mark end of word
  integer            :: Neps              ! number of times tolerance is achieved before termination
  integer            :: nkey              ! word counter
  integer            :: nlines            ! number of input lines
  integer            :: NR                ! number of runs with random parameters to determine  temperature
  integer            :: NS                ! number of loops through all parameters before step  adjustment
  integer            :: NT                ! number of loops before temperature reduction
  integer            :: out               ! unit for main output
  integer            :: seed              ! seed for random number generator
  real(sgl)          :: Cglob             ! adjustment factor for step length for all parameters
  real(sgl)          :: eps               ! convergence criteria
  real(sgl)          :: RT                ! factor for temperature reduction
  real(sgl)          :: T                 ! initial temperature
!
! Data block: possible keywords.
!
  data (keyword(i), i = 1, numkey) / ' ', 'all', 'cglob', 'eps', 'history', 'quasi', 'neps', 'nr', &
    'ns', 'nt', 'write', 'seed', 'sens', 'sobol', 'rt', 't'/
!
! ******************** Defaults for input variables ********************
!
! Defaults are set. These can be overruled in the input file.
!
!   temperature
!   adjustment
!
  flagwrite = .false.
  flagsens = .false.
  flagall = .false.
  flagquasi = .false.
  flagsobol = .false.
  flaghistory = .false.
  flagsave = .false.
  NR = 10
  NT = 5
  NS = 5
  RT = 0.5
  T = 0.
  Cglob = 2.
  Neps = 4
  eps = 1.e-6
  inp = 69
  out = 70
  if (out /= 6) open(unit = out, status = 'unknown', file = 'koekel.out')
  seed = 20032703
!
! ********************** Read input variables **************************
!
! 1. User input file.
!
! We read the complete input file first as a set of character strings.
! The actual keywords will be read from these later on.
!
  open(unit = inp, status = 'unknown', file = 'koekel.inp')
  i = 1
   10 read(inp, '(a)', end = 100) inline(i)
  i = i + 1
  if (i > numlin) then
    write(out, '(" KOEKEL-error: Number of input lines", " exceeds ", i5)') numlin
    write(out, '(" numlin should be increased")')
    stop
  endif
  goto 10
100 close(inp)
  nlines = i - 1
!
! 2. Convert uppercase to lowercase characters.
!
! For easy handling of all the input parameters, everything is converted
! to lowercase characters.
!
  do i = 1, nlines
    do k = 1, 132
      if (inline(i)(k:k) >= 'A' .and. inline(i)(k:k) <= 'Z') inline(i)(k:k) = char(ichar(inline(i)(k:k)) + 32)
    enddo
  enddo
!
! 3. Read keywords and values from input line.
!
! From each input line we retrieve the keyword and its value.
!
  key = '                                                  '
  value = '                                                  '
  param = '                                                  '
  chprev = ' '
  do i = 1, nlines
    nkey = 0
    do k = 1, 132
      if (k > 1) chprev = inline(i)(k - 1:k - 1)
      ch = inline(i)(k:k)
      if (ch /= ' ' .and. chprev == ' ') then
        nkey = nkey + 1
        kbeg = k
      endif
      if (ch == ' ' .and. chprev /= ' ') then
        kend = k - 1
        if (nkey == 1) key = inline(i)(kbeg:kend)
        if (nkey == 2) value = inline(i)(kbeg:kend)
        if (nkey == 3) param = inline(i)(kbeg:kend)
      endif
      if (nkey > 3) goto 230
    enddo
  230   ch = value(1:1)
!
! 4. Check for errors in keywords.
!
! A keyword can be de-activated by putting a # in front of it.
! All first words of the input lines are checked against the list
! of keywords. KOEKEL will stop if a keyword is incorrect.
!
    if (key(1:1) == '#') cycle
    do k = 1, numkey
      if (keyword(k) == key) goto 250
    enddo
    write(out, '(/" KOEKEL-error: Wrong keyword: ", a)') trim(key)
    stop
!
! 5. Identify keywords.
!
  250   if (key == 'nr') then
      read(value, * , err = 400) NR
      cycle
    endif
    if (key == 'ns') then
      read(value, * , err = 400) NS
      cycle
    endif
    if (key == 'nt') then
      read(value, * , err = 400) NT
      cycle
    endif
    if (key == 'neps') then
      read(value, * , err = 400) Neps
      cycle
    endif
    if (key == 'seed') then
      read(value, * , err = 400) seed
      cycle
    endif
    if (key == 'rt') then
      read(value, * , err = 400) RT
      cycle
    endif
    if (key == 't') then
      read(value, * , err = 400) T
      cycle
    endif
    if (key == 'cglob') then
      read(value, * , err = 400) Cglob
      cycle
    endif
    if (key == 'eps') then
      read(value, * , err = 400) eps
      cycle
    endif
    if (key == 'write') then
      if (ch == 'n') flagwrite = .false.
      if (ch == 'y') flagwrite = .true.
      if (ch /= 'y' .and. ch /= 'n') goto 400
      cycle
    endif
    if (key == 'sens') then
      if (ch == 'n') flagsens = .false.
      if (ch == 'y') flagsens = .true.
      if (ch /= 'y' .and. ch /= 'n') goto 400
      cycle
    endif
    if (key == 'all') then
      if (ch == 'n') flagall = .false.
      if (ch == 'y') flagall = .true.
      if (ch /= 'y' .and. ch /= 'n') goto 400
      cycle
    endif
    if (key == 'quasi') then
      if (ch == 'n') flagquasi = .false.
      if (ch == 'y') flagquasi = .true.
      if (ch /= 'y' .and. ch /= 'n') goto 400
      cycle
    endif
    if (key == 'sobol') then
      if (ch == 'n') flagsobol = .false.
      if (ch == 'y') flagsobol = .true.
      if (ch /= 'y' .and. ch /= 'n') goto 400
      cycle
    endif
    if (key == 'history') then
      if (ch == 'n') flaghistory = .false.
      if (ch == 'y') flaghistory = .true.
      if (ch /= 'y' .and. ch /= 'n') goto 400
      cycle
    endif
    if (key == 'save') then
      if (ch == 'n') flagsave = .false.
      if (ch == 'y') flagsave = .true.
      if (ch /= 'y' .and. ch /= 'n') goto 400
      cycle
    endif
  enddo
!
! 6. Check for errors in values.
!
! All input parameters need to fall within certain ranges.
! If not, KOEKEL will stop.
!
  if (NR < 1 .or. NR > 100) then
    write(out, '(" KOEKEL-error: 1 < = NR < = 100")')
    stop
  endif
  if (NS < 1 .or. NS > 100) then
    write(out, '(" KOEKEL-error: 1 < = NS < = 100")')
    stop
  endif
  if (NT < 1 .or. NT > 10) then
    write(out, '(" KOEKEL-error: 1 < = NT < = 10")')
    stop
  endif
  if (Neps < 1 .or. Neps > 10) then
    write(out, '(" KOEKEL-error: 1 < = Neps < = 10")')
    stop
  endif
  if (seed < 1 .or. seed >= 1000000000) then
    write(out, '(" KOEKEL-error: 1 < = seed < 1000000000")')
    stop
  endif
  if (eps < 1.e-30 .or. eps > 1.e30) then
    write(out, '(" KOEKEL-error: 1.e-30 < = eps < = 1.e30")')
    stop
  endif
  if (RT < 0..or.RT >= 1.) then
    write(out, '(" KOEKEL-error: 0. < = RT < 1.")')
    stop
  endif
  if (T < 0..or.T >= 1.e20) then
    write(out, '(" KOEKEL-error: 0. < T < 1.e20")')
    stop
  endif
  if (Cglob < 0..or.Cglob > 100.) then
    write(out, '(" KOEKEL-error: 0. < = Cglob < = 100.")')
    stop
  endif
!
! 7. Main output
!
! Title
!
  write(out, '(/"     KOEKEL-0.3 Beta ", "(Version: March 24, 2008)"/)')
  write(out, '(" Copyright (C) 2008  Arjan Koning"//)')
  write(out, '("      Optimization with simulated annealing"/)')
!
! User input file
!
  write(out, '(/" ########## USER INPUT ##########")')
  write(out, '(/" USER INPUT FILE"/)')
  do i = 1, nlines
    write(out, '(1x, a)') trim(inline(i))
  enddo
!
! User input + defaults
!
  write(out, '(/" USER INPUT FILE + DEFAULTS"/)')
  write(out, '(" Keyword  Value   Explanation"/)')
  write(out, '(" write       ", l1, "       logical", " for extensive output")')  flagwrite
  write(out, '(" quasi       ", l1, "       logical", " for quasi-random search")') flagquasi
  write(out, '(" sobol       ", l1, "       logical", " for quasi-random search with Sobol sequence")') flagsobol
  write(out, '(" history     ", l1, "       logical", " for using information from previous runs")') flaghistory
  write(out, '(" save        ", l1, "       logical", " for saving information from previous runs")') flagsave
  write(out, '(" sens        ", l1, "       logical", " for sensitivity weighted random parameter")') flagsens
  write(out, '(" all         ", l1, "       logical", " for variation of all parameters simultaneously")') flagall
  write(out, '(" NR       ", i5, "      number of", " initial runs with random parameters")') NR
  write(out, '(" NS       ", i5, "      number of loops", " through all parameters before step adjustment")') NS
  write(out, '(" NT       ", i5, "      number of loops", " before temperature reduction")') NT
  write(out, '(" Neps     ", i5, "      number of times", " tolerance is achieved before termination")') Neps
  write(out, '(" seed  ", i8, "      seed for random", " number initialization")') seed
  write(out, '(" eps      ", 1p, e10.3, " convergence criteria")') eps
  write(out, '(" RT       ", f6.3, "     factor for", " temperature reduction")') RT
  write(out, '(" T        ", 1p, e10.3, " initial temperature")') T
  write(out, '(" Cglob    ", f6.3, "     adjustment", " factor for step length for all parameters")') Cglob
!
! Initialize random numbers. For pseudo-random numbers, we use the
! Mersenne twister by Matsumoto and Nakamura. It has a period
! of 2**19937-1. Even my code will have found a solution before
! that period is over.
!
  call init_genrand(seed)
  return
  400 write(out, '(" KOEKEL-error: Wrong input: ", a)') trim(inline(i))
  stop
end subroutine koekelinput
! Copyright A.J. Koning 2021
subroutine koekelsearch(Npar, P, Pmin, Pmax, Popt, Fopt, flagwrite, flagsens, flagall, flagquasi, flagsobol, flaghistory, &
 &  flagsave, NR, NS, NT, Neps, RT, T, Cglob, eps)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Simulated annealing optimization
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer, parameter :: dbl = selected_real_kind(15,307) ! double precision kind
  integer, parameter :: numruns=100000              !
  integer, parameter :: maxdim=1111                 !
  logical            :: flag(2)                     ! distribution
  logical            :: flagall                     ! logical for variation of all parameters simultaneously
  logical            :: flaghistory                 ! logical for using information from previous runs
  logical            :: flagquasi                   ! logical for quasi-random search
  logical            :: flagquit                    !
  logical            :: flagsave                    ! flag to save all input and output files of each  Talys run
  logical            :: flagsens                    ! flag to use sensitivity in parameter weighting
  logical            :: flagsobol                   ! logical for quasi-random search with Sobol sequence
  logical            :: flagwrite                   ! logical for extensive output
  integer            :: base                        ! help variable
  integer            :: i                           ! counter
  integer            :: Npar                        ! number of parameters
  integer            :: indexQ(Npar)                !
  integer            :: inew                        !
  integer            :: ipar                        ! counter for parameter
  integer            :: irun                        ! counter for total number of runs
  integer            :: irunbeg                     !
  integer            :: irunend                     !
  integer            :: irunTbeg                    !
  integer            :: iS                          !
  integer            :: iT                          !
  integer            :: j                           ! counter
  integer            :: jend                        !
  integer            :: k                           ! Legendre order
  integer            :: kacc                        !
  integer            :: krej                        !
  integer            :: ktot                        !
  integer            :: l                           ! counter
  integer            :: lbeg                        ! start of l summation
  integer            :: lend                        ! end of l summation
  integer            :: ll                          ! angular momentum
  integer            :: ltot                        !
  integer            :: Nacc(Npar)                  ! number of accepted values
  integer            :: NaccT(Npar)                 !
  integer            :: Ndown(Npar)                 !
  integer            :: NdownT(Npar)                !
  integer            :: Neps                        ! number of times tolerance is achieved before termination
  integer            :: Nopt(Npar)                  ! number of optimum run
  integer            :: NoptT(Npar)                 !
  integer            :: NR                          ! number of runs with random parameters to determine  temperature
  integer            :: Nrej(Npar)                  ! number of rejected values
  integer            :: NrejT(Npar)                 !
  integer            :: NS                          ! number of loops through all parameters before step  adjustment
  integer            :: NT                          ! number of loops before temperature reduction
  integer            :: Ntot(Npar)                  ! number of nucleon units in exit channel
  integer            :: Nup(Npar)                   !
  integer            :: NupT(Npar)                  !
  integer            :: out                         ! unit for main output
  integer            :: pbeg                        ! help variable
  integer            :: pend                        !
  integer            :: Pindexsave(0:numruns)       !
  integer            :: prime(Npar)                 !
  integer            :: taus                        !
  integer            :: verdict(0:numruns)          !
  integer            :: verloc(numruns / Npar)      !
  integer            :: vertmp                      !
  real(sgl)          :: alpha                       ! amplitude of the density-dependent term of the DDM3Y interaction
  real(sgl)          :: C(Npar)                     ! adjustment factor for step length
  real(sgl)          :: Cglob                       ! adjustment factor for step length for all parameters
  real(sgl)          :: corput                      !
  real(sgl)          :: denom(Npar)                 ! help variable
  real(sgl)          :: domain(Npar)                ! domain (minimal and maximal boundary) for parameters
  real(sgl)          :: enum(Npar)                  ! enumerator of Lorentzian
  real(sgl)          :: eps                         ! convergence criteria
  real(sgl)          :: expo                        ! exponent
  real(sgl)          :: F                           ! function value
  real(sgl)          :: Fdif                        !
  real(sgl)          :: Ffile(0:numruns)            !
  real(sgl)          :: Fnew                        !
  real(sgl)          :: Fnewsave(0:numruns)         !
  real(sgl)          :: Fopt                        ! optimal function value
  real(sgl)          :: Foptsave(0:numruns)         !
  real(sgl)          :: Fprev(0:10)                 ! optimal function value from previous temperature
  real(sgl)          :: Frand(Npar * NT * NS)       !
  real(sgl)          :: Fsave(0:numruns)            !
  real(sgl)          :: Ftmp                        !
  real(sgl)          :: omega                       ! particle-hole state density
  real(sgl)          :: P(Npar)                     ! sampled parameter
  real(sgl)          :: p0                          ! initial particle number
  real(sgl)          :: pardif                      ! difference of parameters
  real(sgl)          :: Pdif                        !
  real(sgl)          :: Pfile(0:numruns, Npar)      !
  real(sgl)          :: Phigh(0:Npar)               !
  real(sgl)          :: Ploc(numruns/Npar)          !
  real(sgl)          :: Plow(0:Npar)                ! bottom of histogram
  real(sgl)          :: Pmax(Npar)                  ! maximal parameters
  real(sgl)          :: Pmin(Npar)                  ! minimal parameters
  real(sgl)          :: Pnew(Npar)                  !
  real(sgl)          :: Pnewsave(0:numruns, Npar)   !
  real(sgl)          :: Poldsum                     !
  real(sgl)          :: Popt(Npar)                  ! optimal parameters
  real(sgl)          :: Poptsave(0:numruns, Npar)   !
  real(sgl)          :: PP                          !
  real(sgl)          :: prob                        ! probability
  real(sgl)          :: prob1                       !
  real(sgl)          :: probsave(0:numruns)         !
  real(sgl)          :: Psave(0:numruns, Npar)      !
  real(sgl)          :: Psum                        !
  real(sgl)          :: Ptmp                        !
  real(sgl)          :: R                           ! random number
  real(sgl)          :: ratio                       ! ratio of uncertainties
  real(sgl)          :: RR                          ! random number
  real(sgl)          :: RT                          ! factor for temperature reduction
  real(sgl)          :: S(Npar)                     ! sensitivity matrix
  real(sgl)          :: Sabs(Npar)                  ! absolute sensitivity
  real(sgl)          :: Sfrac                       !
  real(sgl)          :: stepsum                     !
  real(sgl)          :: Stot                        !
  real(sgl)          :: sum(Npar)                   ! help variable
  real(sgl)          :: T                           ! initial temperature
  real(sgl)          :: Tsave(0:numruns)            !
  real(sgl)          :: V(Npar)                     ! step length
  real(dbl)          :: mersenne                    ! random number generator (Mersenne twister)
  real(dbl)          :: quasi(maxdim)               !
!
! *************************** Initialization ***************************
!
! Fopt       : optimal function value
!   temperature
!   adjustment
!
! Output units
!
  out = 6
!
! For quasi-random numbers, we use the Halton sequence, for which
! we need the first Npar prime numbers, or the Sobol sequence.
!
  if (flagquasi) then
    if (flagsobol) then
      call insobl(flag, Npar, numruns, taus)
    else
      k = 1
      prime(1) = 2
      i = 1
   10     i = i + 2
      jend = int(sqrt(real(i)))
      do j = 2, jend
        if (mod(i, j) == 0) goto 10
      enddo
      k = k + 1
      prime(k) = i
      if (k < Npar) goto 10
    endif
  endif
!
! Initialize SA parameters
!
  do i = 1, Npar
    domain(i) = Pmax(i) - Pmin(i)
    V(i) = 0.5 * domain(i)
    C(i) = Cglob
  enddo
  if (flagsens) then
    Phigh(0) = 0.
    do i = 1, Npar
      Sfrac = 1. / Npar
      Plow(i) = Phigh(i - 1)
      Phigh(i) = Plow(i) + Sfrac
    enddo
  endif
!
! Function evaluation with the starting values. This is also the first
! optimal result.
!
! fcn     : function
! F       : function value
!
  call fcn(Npar, P, F)
  call writestart(Npar, P, Pmin, Pmax, F, out)
  do i = 1, Npar
    Popt(i) = P(i)
    Psave(0, i) = P(i)
    Poptsave(0, i) = P(i)
    indexQ(i) = 0
  enddo
  indexQ(1) = 1
  Fopt = F
  Fsave(0) = F
  Foptsave(0) = F
  verdict(0) = 1
  Tsave(0) = T
  do i = 1, Npar
    P(i) = Popt(i)
  enddo
  F = Fopt
  if (flagsave) then
    open(unit = 11, status = 'unknown', file = 'save')
    write(11, * ) (P(j), j = 1, Npar), F
    close(11)
  endif
!
! *************************** Simulated annealing **********************
!
! Fprev  : optimal function value from previous temperature
! writesa: subroutine for output of simulated annealing parameters
!
  Fprev(0) = Fopt
  do i = 1, Neps
    Fprev(i) = 1.e30
  enddo
  irun = 0
  100 call writesa(Npar, V, T, 0, C, "New     ", out)
  call writeval(Npar, P, F, "Optimal ", out)
!
! Possible action when new optimum is found
!
  if (irun == 0) call optimum
!
  do i = 1, Npar
    NaccT(i) = 0
    NrejT(i) = 0
    NdownT(i) = 0
    NupT(i) = 0
    NoptT(i) = 0
    indexQ(i) = 0
  enddo
  indexQ(1) = 1
!
  irunTbeg = irun
  do iT = 0, NT
    if (iT == 0) then
!
! Determine initial temperature.
! To estimate the variation of function values, given the initial
! boundaries of the parameters, we perform NR*Npar runs with all
! parameters taken random at every run. The initial temperature is
! the average function value. In the process, we also record a
! possible new optimal value.
!
! writeval: subroutine for output of parameters and function value
!
      if (T == 0..or.RT == 0.) then
        do iS = 1, NR
          if (flagquasi .and. flagsobol) call GOSOBL(quasi)
          do i = 1, Npar
            if (flagquasi) then
              if (flagsobol) then
                RR = real(quasi(i))
              else
                base = prime(i)
                indexQ(i) = indexQ(i) + 1
                k = indexQ(i)
                RR = corput(base, k)
              endif
              R = 2. * RR - 1.
            else
              R = 2. * real(mersenne()) - 1.
            endif
            P(i) = Popt(i) + R * V(i)
          enddo
          call fcn(Npar, P, F)
          if (flagwrite) call writeval(Npar, P, F, "Random  ", out)
          Frand(iS) = F
          if (F < Fopt) then
            do i = 1, Npar
              Popt(i) = P(i)
            enddo
            Fopt = F
            call writeval(Npar, Popt, Fopt, "Optimal ", out)
!
! Possible action when new optimum is found
!
            call optimum
!
          endif
        enddo
        do i = 1, NR
          do j = 1, i
            if (Frand(i) < Frand(j)) then
              Ftmp = Frand(j)
              Frand(j) = Frand(i)
              Frand(i) = Ftmp
            endif
    enddo
  enddo
        p0 = 0.5
        Fdif = Frand(NR / 2) - Fopt
        T = - Fdif / log(p0)
      else
        cycle
      endif
    endif
    call writesa(Npar, V, T, iT, C, "Current ", out)
    stepsum = 0.
    do j = 1, Npar
      Nacc(j) = 0
      Nrej(j) = 0
      Ntot(j) = 0
      Ndown(j) = 0
      Nup(j) = 0
      Nopt(j) = 0
      enum(j) = 0.
      denom(j) = 0.
      S(j) = 0.
      stepsum = stepsum + V(j) **2
    enddo
!
! 130: Loop over function evaluations before step adjustment
!
    irunbeg = irun
    alpha = 1. / (NS * NT)
    omega = 2.
    do iS = 1, NS
      do i = 1, Npar
!
! Stop if maximum number of function evaluations is reached
!
        irun = irun + 1
        if (irun > numruns) then
          write(out, '(" KOEKEL-error: Number of function ", &
 &          "evaluations exceeds ", i8, " Increase numruns")') numruns
          irun = irun - 1
          goto 600
        endif
!
! Keep all parameters at their current values, apart from P(i),
! which gets a new value, taken random between P(i)-V(i) and P(i)+V(i).
!
        if (flagsens) then
          R = 2. * real(mersenne()) - 1.
          do j = 1, Npar
            if (R > Plow(j) .and. R <= Phigh(j)) ipar = j
          enddo
        else
          ipar = i
        endif
        if (flagall) then
          pbeg = 1
          pend = Npar
        else
          pbeg = ipar
          pend = ipar
          do j = 1, Npar
            Pnew(j) = P(j)
          enddo
        endif
        if (flagquasi .and. flagsobol) call GOSOBL(quasi)
        do j = pbeg, pend
          if (flagquasi) then
            if (flagsobol) then
              RR = real(quasi(j))
            else
              base = prime(j)
              indexQ(j) = indexQ(j) + 1
              k = indexQ(j)
              RR = corput(base, k)
            endif
            R = 2. * RR - 1.
          else
            R = 2. * real(mersenne()) - 1.
          endif
          PP = P(j) + R * V(j)
          Pnew(j) = PP
          if (Pnew(j) < Pmin(j) .or. Pnew(j) > Pmax(j)) Pnew(j) = Pmin(j) + real(mersenne()) * domain(j)
          Ntot(j) = Ntot(j) + 1
        enddo
        if (flaghistory .and. .not. flagall .and. irun > Npar) then
          Psum = 0.
          do j = 1, Npar
            Psum = Psum + Pnew(j) **2
          enddo
          do j = 1, irun - 1
            Poldsum = 0.
            do k = 1, Npar
              Poldsum = Poldsum + Psave(j, k) **2
            enddo
            Pdif = abs(Psum - Poldsum)
            if (Pdif <= 0.01 * stepsum / Npar) then
              Fnew = 1.e30
              goto 250
            endif
          enddo
          l = 0
          do j = 1, irun - 1
            if (Pindexsave(j) == ipar) then
              l = l + 1
              Ploc(l) = Psave(j, ipar)
              verloc(l) = verdict(j)
            endif
          enddo
          ltot = l
          do l = 1, ltot
            do k = 1, l
              if (Ploc(l) < Ploc(k)) then
                Ptmp = Ploc(k)
                vertmp = verloc(k)
                Ploc(k) = Ploc(l)
                verloc(k) = verloc(l)
                Ploc(l) = Ptmp
                verloc(l) = vertmp
              endif
            enddo
          enddo
          do l = 1, ltot
            if (Pnew(ipar) > Ploc(l) .and. Pnew(ipar) <= Ploc(l + 1)) then
              ll = l
              lbeg = max(1, ll - 5)
              lend = min(ltot, ll + 5)
              krej = 0
              kacc = 0
              do k = lbeg, lend
                if (verdict(k) == 4) then
                  krej = krej + 1
                else
                  kacc = kacc + 1
                endif
              enddo
              ktot = krej + kacc
              if (real(krej) / ktot > 0.7) then
                Fnew = 1.e30
                goto 250
              endif
            endif
          enddo
        endif
        call fcn(Npar, Pnew, Fnew)
  250       if (flagwrite) then
          write(out, '(" -----"/" Run: ", i6)') irun
          call writeval(Npar, P, F, "Current ", out)
          call writeval(Npar, Pnew, Fnew, "Trial   ", out)
          if (Fnew == 1.e30) write(out, '(" Function skipped")')
        endif
!
! 1. If the new value is smaller than the current one, we accept it.
!
        if (Fnew < F) then
          do j = 1, Npar
            P(j) = Pnew(j)
          enddo
          F = Fnew
          verdict(irun) = 2
          probsave(irun) = 0.
!
! 1a. If the new value is smaller than the optimal one, we replace it.
!
          if (F < Fopt) then
            do j = 1, Npar
              Popt(j) = Pnew(j)
              P(j) = Pnew(j)
            enddo
            Fopt = F
            verdict(irun) = 1
            do j = pbeg, pend
              Nopt(j) = Nopt(j) + 1
            enddo
            call writeval(Npar, Popt, Fopt, "Optimal ", out)
!
! Possible action when new optimum is found
!
            call optimum
          endif
          do j = pbeg, pend
            Nacc(j) = Nacc(j) + 1
            Ndown(j) = Ndown(j) + 1
          enddo
          if (flagwrite) write(out, '("     Accept"/)')
        else
!
! 2. If the new value is larger than the current one, we apply the
!    Metropolis criteria.
!
          do j = pbeg, pend
            Nup(j) = Nup(j) + 1
          enddo
          expo = (F - Fnew) / T
          if (expo < 80.) then
            prob = exp(expo)
          else
            prob = 1.e38
          endif
          probsave(irun) = prob
          prob1 = real(mersenne())
!
! 2a. Accept if random number smaller than the temperature function.
!
          if (prob1 < prob) then
            do j = 1, Npar
              P(j) = Pnew(j)
            enddo
            F = Fnew
            do j = pbeg, pend
              Nacc(j) = Nacc(j) + 1
            enddo
            if (flagwrite) write(out, '("     Accept"/)')
            verdict(irun) = 3
          else
!
! 2b. Reject if random number larger than the temperature function.
!
            do j = pbeg, pend
              Nrej(j) = Nrej(j) + 1
            enddo
            if (flagwrite) write(out, '("     Reject"/)')
            verdict(irun) = 4
          endif
        endif
        if (Cglob == 0.) then
          if (verdict(irun) <= 3) then
            do j = pbeg, pend
              V(j) = (1. - alpha) * V(j) + alpha * omega * (Pnew(j) - Psave(irun - 1, j))
            enddo
          endif
        endif
        do j = 1, Npar
          Psave(irun, j) = P(j)
          Pnewsave(irun, j) = Pnew(j)
          Poptsave(irun, j) = Popt(j)
        enddo
        Pindexsave(irun) = ipar
        Tsave(irun) = T
        Fsave(irun) = F
        Fnewsave(irun) = Fnew
        Foptsave(irun) = Fopt
        if (flagsens) then
          pardif = Pnew(ipar) - Psave(irun - 1, ipar)
          Fdif = Fnew - Fsave(irun - 1)
          k = Ntot(ipar)
          if (Fdif /= 0.) S(ipar) = (S(ipar) * (k - 1) + pardif / Fdif) / k
          enum(ipar) = enum(ipar) + pardif * Fdif
          denom(ipar) = denom(ipar) + pardif **2
          sum(ipar) = sum(ipar) + Fdif * Fdif
        endif
        if (flagsave) then
          open(unit = 11, status = 'unknown', file = 'save')
          inew = 0
          do k = 1, irun
            read(11, * ) (Pfile(k, j), j = 1, Npar), Ffile(k)
            if (Fnew >= Ffile(k)) inew = k
          enddo
          close(11)
          open(unit = 11, status = 'unknown', file = 'save')
          do k = 1, inew
            write(11, * ) (Pfile(k, j), j = 1, Npar), Ffile(k)
          enddo
          write(11, * ) (Pnew(j), j = 1, Npar), Fnew
          do k = inew + 1, irun
            write(11, * ) (Pfile(k, j), j = 1, Npar), Ffile(k)
          enddo
          close(11)
        endif
      enddo
    enddo
    irunend = irun
    if (flagsens) then
      write(out,'(/" parameter Sensitivity Abs. sensitivity")')
      Stot = 0.
      do j = 1, Npar
!   if (denom(j).gt.0.) S(j)=enum(j)/denom(j)
!   Stot=Stot+sum(j)
!   Srel(j)=S(j)*Popt(j)/Fopt
        if (S(j) == 0.) S(j) = Sabs(j) / V(j)
        Sabs(j) = abs(S(j) * V(j))
        Stot = Stot + Sabs(j)
        write(out, '(5x, i3, 1p, 2e15.5)') j, S(j), Sabs(j)
      enddo
      Phigh(0) = 0.
      Sfrac = 1. / Npar
      do i = 1, Npar
        if (Stot > 0.) Sfrac = Sabs(i) / Stot
        Plow(i) = Phigh(i - 1)
        Phigh(i) = Plow(i) + Sfrac
      enddo
    endif
!
! Output of total number of acceptations and rejections and
! adjustment of the step length.
!
!   if (Cglob.eq.0.) then
!   do 442 k=1,Npar
!   do 443 i=irunbeg+1,irunend
!   Porder(i-irunbeg)=Pnewsave(i,k)
! 443       continue
!   do 444 i=1,irunend-irunbeg
!   do 444 j=1,i
!   if (Porder(i).lt.Porder(j)) then
!   Ptmp=Porder(j)
!   Porder(j)=Porder(i)
!   Porder(i)=Ptmp
!   endif
! 444       continue
! 442     continue
!   do 445 i=1,NS
! 445     continue
!   endif
    if (Cglob > 0.) then
      do i = 1, Npar
        ratio = real(Nacc(i)) / Ntot(i)
        if (ratio > 0.6) V(i) = V(i) * (1 + C(i) * (ratio - 0.6) / 0.4)
        if (ratio < 0.4) V(i) = V(i) / (1 + C(i) * (0.4 - ratio) / 0.4)
        V(i) = min(V(i), domain(i))
      enddo
    endif
    if (flagwrite) call writediag(Npar, T, iT, Nacc, Nrej, Ndown, Nup, Nopt, out)
    do i = 1, Npar
      NaccT(i) = NaccT(i) + Nacc(i)
      NrejT(i) = NrejT(i) + Nrej(i)
      NdownT(i) = NdownT(i) + Ndown(i)
      NupT(i) = NupT(i) + Nup(i)
      NoptT(i) = NoptT(i) + Nopt(i)
    enddo
  enddo
  call writediag(Npar, T, 0, NaccT, NrejT, NdownT, NupT, NoptT, out)
!
! end of all loops within the same temperature
!
! Decide whether to go back for another loop.
!
  do i = Neps - 1, 0, -1
    Fprev(i + 1) = Fprev(i)
  enddo
  Fprev(0) = Fopt
  flagquit = .true.
  do i = 1, Neps
    if (Fprev(i) - Fopt > eps) flagquit = .false.
  enddo
  if (flagquit) then
    goto 600
  else
    do j = 1, Npar
      P(j) = Popt(j)
      indexQ(j) = 0
    enddo
    indexQ(1) = 1
    F = Fopt
    if (RT > 0.) T = RT * T
    goto 100
  endif
  600 write(out, '(/" Total number of runs: ", i8)') irun
  call writeval(Npar, Popt, Fopt, "Final   ", out)
  close(out)
  return
end subroutine koekelsearch
! Copyright A.J. Koning 2021
subroutine writestart(Npar, P, Pmin, Pmax, F, unit)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of begin situation
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer   :: i                     ! counter
  integer   :: Npar                  ! number of parameters
  integer   :: unit                  !
  real(sgl) :: F                     ! function value
  real(sgl) :: P(Npar)               ! sampled parameter
  real(sgl) :: Pmax(Npar)            ! maximal parameters
  real(sgl) :: Pmin(Npar)            ! minimal parameters
!
! *************************** Output ***********************************
!
  write(unit, '(/" ########### Initial parameters ##########"/)')
  write(unit, '("   #    Lower bound   Starting values", " Upper bound"/)')
  do i = 1, Npar
    write(unit, '(1x, i3, 1p, 3e15.5)') i, Pmin(i), P(i), Pmax(i)
  enddo
  write(unit, '(/" Initial function value: ", 1p, e12.5)') F
  write(unit, '(/" #####################")')
  return
end subroutine writestart
! Copyright A.J. Koning 2021
subroutine writeval(Npar, P, F, type, unit)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of parameters and function value
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  character(len=8) :: type        ! particle type
  integer          :: i           ! counter
  integer          :: Npar        ! number of parameters
  integer          :: unit        !
  real(sgl)        :: F           ! function value
  real(sgl)        :: P(Npar)     ! sampled parameter
!
! *************************** Output ***********************************
!
  write(unit, '(/1x, a8, " parameters"/)') type
  write(unit, '("   #    Value")')
  do i = 1, Npar
    write(unit, '(1x, i3, 1p, e15.5)') i, P(i)
  enddo
  write(unit, '(/1x, a8, " function value:", 1p, e12.5)') type, F
  return
end subroutine writeval
! Copyright A.J. Koning 2021
subroutine writesa(Npar, V, T, iT, C, type, unit)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Output of simulated annealing parameters
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  character(len=8) :: type        ! particle type
  integer          :: i           ! counter
  integer          :: iT          !
  integer          :: Npar        ! number of parameters
  integer          :: unit        !
  real(sgl)        :: C(Npar)     ! adjustment factor for step length
  real(sgl)        :: T           ! initial temperature
  real(sgl)        :: V(Npar)     ! step length
!
! *************************** Output ***********************************
!
  write(unit, '(/" ########## ", a8, " Temperature ##########"/)') type
  write(unit, '(" Temperature:", 1p, e12.5)') T
  if (iT /= 0) write(unit, '(" iT loop    :", i5)') iT
  write(unit, '(/"  #    step length  adjustment factor")')
  do i = 1, Npar
    write(unit, '(i3, 1p, 2e15.5)') i, V(i), C(i)
  enddo
  write(unit, '(/" ####################")')
  return
end subroutine writesa
! Copyright A.J. Koning 2021
subroutine writediag(Npar, T, iT, Nacc, Nrej, Ndown, Nup, Nopt, unit)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Diagnosis of uphill and downhill moves, accepts and rejects
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer   :: i                      ! counter
  integer   :: iT                     !
  integer   :: Npar                   ! number of parameters
  integer   :: Nacc(Npar)             ! number of accepted values
  integer   :: Nacctot                ! total number of acceptations
  integer   :: Ndown(Npar)            !
  integer   :: Ndowntot               !
  integer   :: Nopt(Npar)             ! number of optimum run
  integer   :: Nopttot                ! number of optimum run for total
  integer   :: Nrej(Npar)             ! number of rejected values
  integer   :: Nrejtot                ! total number of rejections
  integer   :: Ntot                   ! number of nucleon units in exit channel
  integer   :: Nup(Npar)              !
  integer   :: Nuptot                 !
  integer   :: unit                   !
  real(sgl) :: T                      ! initial temperature
!
! *************************** Output ***********************************
!
  write(unit, '(/" ########## Diagnosis ##########"/)')
  write(unit, '(" Temperature:", 1p, e12.5)') T
  if (iT > 0) write(unit, '(" iT loop    :", i5)') iT
  write(unit, '(/"  #    Accept  Reject Downhill Uphill  Optimum")')
  Nacctot = 0
  Nrejtot = 0
  Ndowntot = 0
  Nuptot = 0
  Nopttot = 0
  do i = 1, Npar
    Nacctot = Nacctot + Nacc(i)
    Nrejtot = Nrejtot + Nrej(i)
    Ndowntot = Ndowntot + Ndown(i)
    Nuptot = Nuptot + Nup(i)
    Nopttot = Nopttot + Nopt(i)
    write(unit, '(i3, 5(2x, i6))') i, Nacc(i), Nrej(i), Ndown(i), Nup(i), Nopt(i)
  enddo
  Ntot = Nacctot + Nrejtot
  write(unit, '(/" Total number of moves         : ", i5)') Ntot
  write(unit, '(" Total number of accepted moves: ", i5)') Nacctot
  write(unit, '(" Total number of rejected moves: ", i5)') Nrejtot
  write(unit, '(" Total number of downhill moves: ", i5)') Ndowntot
  write(unit, '(" Total number of uphill moves  : ", i5)') Nuptot
  write(unit, '(" Total number of new optima    : ", i5)') Nopttot
  write(unit, '(/" ####################")')
  return
end subroutine writediag
! Copyright A.J. Koning 2021
