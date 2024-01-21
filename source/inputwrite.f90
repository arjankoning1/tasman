subroutine inputwrite
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write parameters for TALYS input file
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
! Variables for reading TALYS output
!   flaggamma       ! flag for covariances or optimization of discrete gamma ray transitions
!   flagresidual    ! flag for covariances or optimization of residual production cross sections
! Variables for fitting limits
!   parclass        ! class of keyword
! Variables for GOF function
!   isearch         ! number of trial run
! Variables for reading TASMAN input lines
!   inline          ! input line
!   nlines          ! number of input lines
! Variables for parameter variation
!   Npar            ! number of parameters
!   par             ! parameter value
!   paradjust       ! local adjustment parameters
!   parD            ! random parameter value
!   parkey          ! TALYS keyword
!   parname         ! character parameter value
! Variables for processing input
!   k0              ! index of incident particle
! Variables for writing TALYS input files
!   A               ! mass number of nucleus
!   Apar            ! A - index for keyword
!   barrier         ! fission barrier
!   Ehigh           ! energy cutoff of high energy runs
!   first           ! tracks whether a TALYS run is the first one
!   flagangle       ! flag for covariances or optimization of angular distributions
!   flagecis        ! flag to repeat ECIS calculations in a second run
!   flagglobal      ! flag to use global TALYS cross sections as prior
!   flaginponly     ! flag to produce only TALYS input files, no runs
!   flagleg         ! flag for covariances or optimization of Legendre coefficients
!   flagomponly     ! flag for covariances or optimization for basic optical model observables
!   flagompvar      ! flag for variation of OMP parameters
!   flagspectra     ! flag for covariances or optimization of emission spectra
!   flagtmc         ! flag for Total Monte Carlo (ENDF - 6 data file at each run)
!   inuc            ! index of nuclide
!   italys          ! TALYS run counter
!   lval            ! l value
!   mode            ! TASMAN mode
!   Nhigh           ! number of high energy runs
!   pargam          ! gamma - ray related input string
!   parmt           ! MT number
!   parsym          ! particle symbol
!   string          ! line with parameter value
!   Z               ! charge number of nucleus
!   Zpar            ! Z - index for keyword
!
! *** Declaration of local data
!
  implicit none
  character(len=21)  :: tfile     ! TALYS input file
  character(len=132) :: line      ! input line
  character(len=132) :: rangestring 
  integer            :: i         ! counter
  integer            :: ii        ! counter
  integer            :: ilast     ! pdeng
  integer            :: j         ! counter
  integer            :: system    ! system call
!
! ************************** Write input file **************************
!
! Generally, inputwrite is first called for the original set of parameters, and in later runs for varied parameters
!
! Original input file
!
  open (unit = 1, file = 'talys.inp', status = 'replace')
  write(1, '("# Original input file")')
  do i = 1, nlines
    write(1, '(a)') trim(inline(i))
  enddo
  if (isearch == 0) open(unit=2, file='parameters_range.dat', status='unknown')
!
! Keywords added on the basis of TASMAN input
!
  write(1, '("# Current variables Run:", i6)') isearch
  if (flagomponly) then
    write(1, '("maxz 0")')
    write(1, '("maxn 0")')
    write(1, '("bins 10")')
    write(1, '("preequilibrium n")')
  endif
  if (flagspectra) then
    write(1, '("outspectra y")')
    write(1, '("filespectrum n p d t h a")')
  endif
  if (flagangle) then
    write(1, '("outangle y")')
    write(1, '("fileelastic y")')
  endif
  if (flagleg) write(1, '("outlegendre y")')
  if (flaggamma) write(1, '("filegamdis y")')
  if (flagresidual) write(1, '("fileresidual y")')
  write(1, '("#Ztarget ",i3)') Ztarget
  write(1, '("partable y")')
!
! Possible economization for OMP calculations when no OMP parameters are varied.
!
  write(1, '("ecissave y")')
  if (first .and. flagecis) then
    write(1, '("inccalc y")')
    write(1, '("eciscalc y")')
    if (flagtmc) write(1, '("endfecis y")')
    if ( .not. flagompvar) first = .false.
  else
    if (k0 == 0 .or. (flagompvar .and. flagecis)) then
      write(1, '("inccalc y")')
      write(1, '("eciscalc y")')
      if (flagtmc) write(1, '("endfecis y")')
    else
      if ( .not. first) write(1, '("inccalc n")')
      if ( .not. first) write(1, '("eciscalc n")')
      if (flagtmc) write(1, '("endfecis n")')
    endif
    if (.not.flagompinc) write(1, '("inccalc n")')
  endif
  if (mode == 4) then
    write(1, '("element ", i3)') Z(inuc)
    write(1, '("mass    ", i3)') A(inuc)
  endif
!
! Write altered parameters
!
  do i = 1, Npar
    line = '                                                         '
    line(1:12) = parkey(i)(1:12)
    ii = 14
    if (parclass(i) >= 2 .and. parclass(i) <= 6) then
      write(line(14:16), '(i3)') Zpar(i)
      write(line(18:20), '(i3)') Apar(i)
      ii = 22
    endif
    if (parclass(i) >= 7 .and. parclass(i) <= 9) then
      write(line(14:14), '(a1)') parsym(i)
      ii = 16
    endif
    if (parclass(i) == 10) then
      write(line(14:16), '(i3)') parmt(i)
      line(18:33) = parname(i)(1:16)
      ii = 35
    endif
    if (parclass(i) == 1 .or. parclass(i) == 2 .or. parclass(i) == 4 .or. parclass(i) >= 6) then
      if (par(i) >= 0.1) then
        write(line(ii:ii+9), '(f10.5)') par(i)
        ii = ii + 11
      else
        write(line(ii:ii+11), '(es12.5)') par(i)
        ii = ii + 13
      endif
    endif
    if (parclass(i) == 0 .or. parclass(i) == 5) then
      write(line(ii:ii+4), '(i5)') int(par(i))
      ii = ii + 6
    endif
    if (parclass(i) == 3) then
      line(ii:ii + 15) = parname(i)(1:16)
      ii = ii + 17
    endif
    if (parclass(i) == 4 .or. parclass(i) == 5) then
      write(line(ii:ii), '(i1)') barrier(i)
      ii = ii + 2
    endif
    if (parclass(i) == 9) then
      write(line(ii:ii+1), '(i2)') lval(i)
      ii = ii + 3
    endif
    if (parclass(i) == 6) then
      line(ii:ii + 1) = pargam(i)(1:2)
      ii = ii + 3
    endif
    if (parD(i) > 0.) then
      do j = 1, 3
        write(line(ii:ii+7), '(f8.3)') paradjust(i, j)
        ii = ii + 9
      enddo
      write(line(ii:ii+7), '(f8.3)') parD(i)
      ii = ii + 9
    endif
!
! Bug fix by Puran Deng: some keywords may be more than 12-character long
!
    if (len_trim(parkey(i)) > 12) then
      ilast = len_trim(line)
      if (ilast == 132) then
        write(*,'(" TASMAN-error: input line is too long, more than 132 characters")')
        stop
      endif
      line = trim(parkey(i))//line(13:ilast)
    endif
!
    write(1, '(a)') trim(line)
    string(i) = line
    if (isearch == 0) then
      if (parclass(i) <= 9) then
        rangestring = ' # Range:              -   '
        write(rangestring(11:22),'(f12.5)') parlow(i)
        write(rangestring(26:37),'(f12.5)') parhigh(i)
        write(2,'(2a)') trim(line),trim(rangestring)
      endif
    endif
  enddo
  if (isearch == 0) close(2)
  if (flagglobal .and. italys == -1) write(1, '("best n")')
  if (mode <= 2 .and. italys == 0) write(1, '("bestend y")')
  if (italys > Nhigh) write(1, '("Estop ", f12.5)') Ehigh
  close (1)
  if (Npar > 0) write(*, '(/1x, " Parameters for input file:", i6/)') italys
  do i = 1, Npar
    write(*, '(1x, a)') trim(string(i))
  enddo
!
! Save input file
!
  if (flaginponly) then
    tfile = 'talys.inp.0000'
    write(tfile(11:14), '(i4.4)') italys
    i = system('cat talys.inp > '//tfile)
  endif
  return
end subroutine inputwrite
! Copyright A.J. Koning 2021
