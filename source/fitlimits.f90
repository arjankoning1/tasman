subroutine fitlimits
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Set parameter limits for optimization
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
!   sgl            ! single precision kind
! All global variables
!   numhist        ! maximum number of points for histogram
! Variables for writing TALYS input files
!   Apar           ! A - index for keyword
!   barrier        ! fission barrier
!   parsym         ! particle symbol
!   Zpar           ! Z - index for keyword
! Variables for path names
!   tasmanpath     ! directory containing files to be read
! Variables for processing input
!   Atarget        ! mass number of target nucleus
!   isochar        ! symbol for isomer
!   Liso           ! number of the isomeric state
!   Starget        ! symbol of target nucleus
! Variables for GOF function
!   Nburn          ! number of burn - in runs
! Variables for weights
!   flagweight     ! flag to use weights for random samples
! Variables for parameter variation
!   Npar           ! number of parameters
!   pardelta       ! uncertainty of parameter
!   parinp         ! input parameter value
!   parkey         ! TALYS keyword
!   partype        ! type of parameter variation
! Variables for fitting limits
!   flagreadpar    ! flag to read in parameter distribution
!   hist           ! value in histogram bin
!   histsam        ! value in histogram bin
!   nhist          ! number of histogram bins
!   Nhistbin       ! number of histogram bins
!   Nhistbinall    ! global number of histogram bins
!   nhistsam       ! number of hits in bin
!   parclass       ! class of keyword
!   parfile        ! file with parameters
!   parhigh        ! upper value of parameter uncertainty
!   parlow         ! lower value of parameter uncertainty
!   Pbin           ! central value of histogram bin
!   Pbot           ! lower value of histogram bin
!   readpar        ! flag to read parameter distribution
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=4)  :: massstring    ! string for mass number
  character(len=132):: phisfile      ! file with parameter histograms
  logical           :: lexist        ! logical to determine existence
  integer           :: i             ! counter
  integer           :: istat         ! logical for file access
  integer           :: j             ! counter
  integer           :: jj            ! counter
  integer           :: Nsamp         ! number of samples
  real(sgl)         :: dP            ! parameter bin width
  real(sgl)         :: p1            ! first parameter
  real(sgl)         :: p2            ! second parameter
  real(sgl)         :: Prange        ! parameter range
  real(sgl)         :: var           ! parameter variation
!
! ********************* Set parameter limits ***************************
!
! For parameters of the 'shift ' type a shift is applied.
! Otherwise parameter variation is done with relative ratios.
!
  if (flagweight) Nhistbinall = int(1.+2*Nburn**0.333333)
  do i = 1, Npar
    Nhistbin(i) = Nhistbinall
    if (partype(i) == 'shift ') then
      p1 = parinp(i) - pardelta(i)
      p2 = parinp(i) + pardelta(i)
    else
      var = 1. + pardelta(i)
      p1 = parinp(i) / var
      p2 = parinp(i) * var
    endif
!
! Make sure there is a valid parameter range
!
    parlow(i) = min(p1, p2)
    parhigh(i) = max(p1, p2)
    if (parlow(i) == parhigh(i)) parhigh(i) = parhigh(i) * 1.001
!
! Make histogram for parameter
!
    if (flagweight) then
      Prange = parhigh(i) - parlow(i)
      dP = Prange / Nhistbin(i)
      do j = 1, Nhistbin(i) + 1
        Pbot(i, j) = parlow(i) + (j - 1) * dP
        Pbin(i, j) = parlow(i) + (j - 0.5) * dP
      enddo
      do j = 1, numhist
        hist(i, j) = 0.
        nhist(i, j) = 0
        histsam(i, j) = 0.
        nhistsam(i, j) = 0
      enddo
!
! Create filenames for parameters
!
      jj = 11
      do j = 1, 12
        if (parkey(i)(j:j) == ' ') then
          jj = j - 1
          exit
        endif
      enddo
      parfile(i)(1:jj) = parkey(i)(1:jj)
      if (parclass(i) >= 2 .and. parclass(i) <= 6) then
        write(parfile(i)(jj+1:jj+3), '(i3.3)') Zpar(i)
        write(parfile(i)(jj+4:jj+6), '(i3.3)') Apar(i)
        jj = jj + 6
      endif
      if (parclass(i) == 4 .or. parclass(i) == 5) then
        write(parfile(i)(jj+1:jj+3), '("_", i2.2)') barrier(i)
        jj = jj + 3
      endif
      if (parclass(i) >= 7 .and. parclass(i) <= 9) then
        write(parfile(i)(jj+1:jj+1), '(a1)') parsym(i)
        jj = jj + 1
      endif
!
! If requested, read parameter distributions from file
!
      if (flagreadpar) then
        massstring = '000'
        write(massstring(1:3), '(i3.3)') Atarget
        if (Liso == 0) massstring = trim(massstring)//isochar
        phisfile = trim(tasmanpath) // 'parameters/' // trim(Starget) // trim(massstring) // '/his.' // parfile(i)
        inquire (file = phisfile, exist = lexist)
        if (lexist) then
          open (unit = 2, file = phisfile, status = 'old', iostat = istat)
          if (istat /= 0) call read_error(phisfile, istat)
          read(2, '(18x, i5)', iostat = istat) Nsamp
          if (istat /= 0) call read_error(phisfile, istat)
          if (Nsamp < 10) then
            close (unit = 2)
            cycle
          endif
          readpar(i) = 1
          read(2, '(18x, e12.5, 3x, e12.5/)', iostat = istat) parlow(i), parhigh(i)
          if (istat /= 0) call read_error(phisfile, istat)
          read(2, '(20x, i3, /)', iostat = istat) Nhistbin(i)
          if (istat /= 0) call read_error(phisfile, istat)
          do j = 1, Nhistbin(i)
            read(2, '(2e12.5, i5)', iostat = istat) Pbin(i, j), hist(i, j), nhist(i, j)
            if (istat /= 0) call read_error(phisfile, istat)
          enddo
          close (unit = 2)
          Prange = parhigh(i) - parlow(i)
          dP = Prange / Nhistbin(i)
          Pbin(i, Nhistbin(i) + 1) = Pbin(i, Nhistbin(i)) + dP
          do j = 1, Nhistbin(i) + 1
            Pbot(i, j) = Pbin(i, j) - 0.5 * dP
          enddo
        endif
      endif
    endif
  enddo
  return
end subroutine fitlimits
! Copyright A.J. Koning 2021
