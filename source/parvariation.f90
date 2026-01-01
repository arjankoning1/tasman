subroutine parvariation
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Variation of model parameters
!
! Author    : Arjan Koning
!
! 2025-12-25: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl            ! single precision kind
!   dbl            ! double precision kind
! All global variables
!   numhist        ! maximum number of points for histogram
!   numpar         ! maximum number of parameters
! Variables for writing TALYS input files
!   italys         ! TALYS run counter
!   mode           ! TASMAN mode
!   parsym         ! particle symbol
! Variables for GOF function
!   Nburn          ! number of burn - in runs
! Variables for fitting limits
!   hist           ! value in histogram bin
!   Nhistbin       ! number of histogram bins
!   parhigh        ! upper value of parameter uncertainty
!   parlow         ! lower value of parameter uncertainty
!   Pbot           ! lower value of histogram bin
!   readpar        ! flag to read parameter distribution
! Variables for weights
!   Gsave          ! GOF value
! Variables for parameter variation
!   flagequisam    ! flag to use equal binning uniform sampling
!   flaggauss      ! flag to sample parameters from a Gaussian instead of uniform
!   flagkdburn     ! flag to use KD03 OMP function for burn - in
!   flagmcak       ! flag to use sampling from arbitrary distribution
!   flagmetro      ! flag for Metropolis sampling
!   Npar           ! number of parameters
!   par            ! parameter value
!   paradjust      ! local adjustment parameters
!   parD           ! random parameter value
!   pardelta       ! uncertainty of parameter
!   parinp         ! input parameter value
!   parinpname     ! input character parameter value
!   parkey         ! TALYS keyword
!   parname        ! character parameter value
!   parsave        ! parameter value
!   partalys       ! parameter value
!   partype        ! type of parameter variation
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical   :: flagaccept                 ! flag for acceptance of random sample
  integer   :: i                          ! counter
  integer   :: isamp                      ! variable for 0th or random run
  integer   :: istat                      ! logical for file access
  integer   :: j                          ! counter
  integer   :: k                          ! counter
  integer   :: jj                         ! counter
  integer   :: nfis                       ! parameter for fission
  integer   :: Nrepeat                    ! number of allowed repetitions for sampling
  integer   :: ns                         ! parameter for fission
  integer   :: pfis                       ! parameter for fission
  real(sgl) :: alpha                      ! amplitude of the density-dependent term of the DDM3Y interaction
  real(sgl) :: efis                       ! parameter for fission
  real(sgl) :: extpar(numpar)             ! external parameters
  real(sgl) :: gasdev                     ! Gaussian random number
  real(sgl) :: jfis                       ! parameter for fission
  real(sgl) :: parbase(numpar)            ! parameter
  real(sgl) :: phist(numhist)             ! parameter histogram value
  real(sgl) :: Plow(numhist)              ! bottom of histogram
  real(sgl) :: Rdist                      ! random number
  real(sgl) :: RR                         ! random number
! real(sgl) :: rnlow                      ! 
  real(dbl) :: mersenne                   ! random number generator (Mersenne twister)
!
! ******************* Monte Carlo variation of parameters **************
!
! Assign random numbers to parameters
!
! Metropolis sampling
!
  if (flagmetro) then
    flagaccept = .true.
    if (italys > 1) then
      alpha = Gsave(italys - 1) / Gsave(italys - 2)
      if (alpha < 1.) then
        RR = real(mersenne())
        if (RR > alpha) flagaccept = .false.
      endif
    endif
    if (flagaccept) then
      do i = 1, Npar
        parbase(i) = parsave(italys - 1, i)
      enddo
    else
      do i = 1, Npar
        parbase(i) = parsave(italys - 2, i)
      enddo
    endif
  else
    do i = 1, Npar
      parbase(i) = parinp(i)
    enddo
  endif
!
! New set of parameters
!
! Option for external parameters by Erik Sunden
!
  if (flagextparvar) read(80,*) (extpar(k),k=1,Npar)
  if (mode == 2 .and. flagmorris) then
    isamp = 0
  else
    isamp = 1
  endif
  do i = 1, Npar
    Nrepeat = 0
!
! 1. Normal parameters
!
    if (parkey(i) /= 'hbtransfile' .and. parkey(i) /= 'class2file') then
 Loop1: do
        do
          if (flaggauss) then
            RR = gasdev()
          else
            RR = 2. * (real(mersenne()) - 0.5)
          endif
          if (paradjust(i, 4) > 0.) then
            if (RR < 0.) then
              parD(i) = paradjust(i, 4) / (1. + abs(RR) * pardelta(i))
            else
              parD(i) = paradjust(i, 4) * (1. + RR * pardelta(i))
            endif
            par(i) = parbase(i)
            exit Loop1
          endif
          if (partype(i) == 'model ') then
            RR = real(mersenne())
            if (parkey(i) == 'massmodel' .or. parkey(i) == 'jlmmode' .or. parkey(i) == 'widthmode') then
              par(i) = RR * pardelta(i)
            else
              par(i) = 1. + RR * pardelta(i)
            endif
            exit Loop1
          endif
          if (flagmcak) then
            do j = 1, Nhistbin(i) + 1
              phist(j) = hist(i, j)
              Plow(j) = Pbot(i, j)
            enddo
            call mcak(i, readpar(i), Plow, Nhistbin(i) + 1, phist, italys, Nburn, Rdist, flagkdburn, flagequisam)
          else
            if (partype(i) == 'factor' .and. RR < 0.) then
              Rdist = parbase(i) / (1. + abs(RR) * pardelta(i))
            else
              Rdist = parbase(i) + RR * pardelta(i)
              if (partype(i) /= 'shift ' .and. Rdist < 0.) then
                if (Nrepeat <= 10) then
                  Nrepeat = Nrepeat + 1
                  cycle
                else
                  Rdist = parlow(i) + real(mersenne()) * (parhigh(i) - parlow(i))
                endif
              endif
            endif
          endif
          if (Rdist < parlow(i) .or. Rdist > parhigh(i)) then
            if (Nrepeat <= 10) then
              Nrepeat = Nrepeat + 1
              cycle
            else
              Rdist = parlow(i) + real(mersenne()) * (parhigh(i) - parlow(i))
            endif
          endif
          exit
        enddo
        par(i) = Rdist
        if (parkey(i) == 'a') par(i) = max(par(i), 1.)
        if (parkey(i) == 'e0') par(i) = min(max(par(i), -10.), 10.)
        if (parkey(i) == 't') par(i) = min(max(par(i), 0.001), 10.)
        if (parkey(i) == 'nlow') then
          par(i) = max(par(i), 2.)
          par(i) = min(par(i), 8.)
        endif
        if (parkey(i) == 'ntop') then
          par(i) = max(par(i), 20.)
        endif
        if (parkey(i) == 'fishw') par(i) = max(0.01, min(par(i), 10.))
        if (parkey(i) == 'gamgamadjust') par(i) = max(0.1, min(par(i), 10.))
        if (parkey(i) == 'deltaw') par(i) = min(max(par(i), -20.), 20.)
        exit
      enddo Loop1
      if (flagextparvar) then
        par(i) = extpar(i)
      endif
      partalys(isamp, i) = par(i)
      parsave(italys, i) = par(i)
    else
!
! 2. Randomize information in structure files
!
      open(unit = 11, file = parinpname(i), status = 'old')
      if (istat /= 0) call read_error(parinpname(i), istat)
      parname(i) = trim(parinpname(i))//'.new'
      open(unit = 12, file = parname(i), status = 'replace')
      do
        read(11, '(2i4)', iostat = istat) ns, nfis
        if (istat ==  -1) exit
        if (istat /= 0) call read_error(parinpname(i), istat)
        write(12, '(2i4)') ns, nfis
        do jj = 1, nfis
          read(11, '(4x, f11.6, f6.1, i5)', iostat = istat) efis, jfis, pfis
          if (istat ==  -1) exit
          if (istat /= 0) call read_error(parinpname(i), istat)
          if (flaggauss) then
            RR = gasdev()
          else
            RR = 2. * (real(mersenne()) - 0.5)
          endif
          Rdist = 1. + RR * pardelta(i)
          write(12, '(i4, f11.6, f6.1, i5)') jj, efis*Rdist, jfis, pfis
        enddo
      enddo
      close (11)
      close (12)
    endif
  enddo
!
! Let knock-out and stripping parameters vary in sync
!
  do i = 1, Npar
    if (parkey(i) == 'cstrip') then
      do j = 1, Npar
        if (parkey(j) == 'cknock' .and. parsym(j) == parsym(i)) then
          par(i) = par(j)
          partalys(isamp, i) = par(j)
          parsave(italys, i) = par(j)
        endif
      enddo
    endif
  enddo
  return
end subroutine parvariation
! Copyright A.J. Koning 2025
