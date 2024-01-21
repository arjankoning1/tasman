subroutine libread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read library data
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
!   numenexp       ! maximum number of energies forexperimental data set
!   numlib         ! maximum number of data libraries
!   nummt          ! maximum number of MT numbers
!   numsets        ! maximum number of experimental data sets
! Variables for reading TALYS output
!   Nchanall       ! total number of channels
! Constants
!   xseps          ! limit for cross sections
! Variables for path names
!   librariespath  ! directory containing files to be read
! Variables for processing input
!   Atarget        ! mass number of target nucleus
!   auto           ! flag to in - or exclude channel
!   libinclude     ! flag to include data library for search
!   Liso           ! number of the isomeric state
!   mtlib          ! flag for library per channel
!   numE           ! number of incident energies
!   ptype0         ! type of incident particle
!   Starget        ! symbol of target nucleus
! Variables for deviation
!   devi           ! parameter for energy dependence of deviation
! Variables for reading cross sections
!   MT             ! MT number
!   MTiso          ! isomer of MT number
! Variables for reading experimental data
!   auth           ! author
!   lib            ! data library
!   dxsexp         ! uncertainty of experimental cross section
!   Eexp           ! incident energy
!   Ein            ! incident energy
!   Esearch1       ! start energy of search
!   Esearch2       ! end energy of search
!   Eweight        ! incident energy
!   isoexp         ! isomer of experimental data
!   isolib         ! isomer of MT number of library
!   MTlibk         ! MT number of library
!   mtlibweight    ! weight for library and channel
!   Nenexp         ! number of experimental data points
!   NMTlib         ! number of MT numbers of library
!   Nsets          ! number of experimental data sets
!   xsexp          ! experimental cross section
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: numexpen=150000   ! number of energies for experimental data
  logical            :: lexist            ! logical to determine existence
  character(len=1)   :: symiso            ! isomer symbol
  character(len=4)   :: massstring        ! string for mass number
  character(len=4)   :: MTstring          ! string for MT number
  character(len=15)  :: libname(nummt)    ! library name
  character(len=30)  :: mtfile            ! file for MT number
  character(len=30)  :: author            ! author
  character(len=30)  :: libf0(numsets)    ! library name
  character(len=132) :: filehead          ! head of filename
  character(len=132) :: filename          ! filename
  character(len=132) :: lfile             ! file with library
  character(len=132) :: libpath           ! directory containing nuclear data libraries
  character(len=132) :: line              ! 
  character(len=132) :: key
  integer            :: i                 ! counter
  integer            :: ilib              ! counter for library
  integer            :: imt               ! MT counter
  integer            :: iso               ! counter for isotope
  integer            :: istat             ! logical for file access
  integer            :: keyix
  integer            :: j                 ! counter
  integer            :: k                 ! Legendre order
  integer            :: l                 ! counter
  integer            :: n                 ! counter
  integer            :: nbeg              ! first outgoing energy
  integer            :: nE                ! number of energies
  integer            :: nend              ! help variable
  integer            :: nenlfile          ! number of incident energies in library data file
  integer            :: nexp              ! number of experimental data sets
  integer            :: nexp0             ! number of experimental data sets
  integer            :: nprev             ! help variable
  real(sgl)          :: deltaE            ! energy bin around outgoing energies
  real(sgl)          :: dxs(0:numexpen)   ! uncertainty of experimental cross section
  real(sgl)          :: Ea                ! parameters for local energy dependence
  real(sgl)          :: Eb                ! end energy of local adjustment
  real(sgl)          :: ee(0:numexpen)    ! incident energy
  real(sgl)          :: ee0               ! incident energy
  real(sgl)          :: En                ! outgoing energy  term03, conv03, enhance03 : help variable for breakup enhancement calcu
  real(sgl)          :: xs(0:numexpen)    ! experimental cross section
  real(sgl)          :: xs0               ! cross section
!
! ************************ Read library data ***************************
!
  massstring = '000 '
  if (Liso > 0) then
    if (Liso == 1) then
      massstring(4:4) = 'm'
    else
      massstring(4:4) = 'n'
    endif
  endif
  write(massstring(1:3), '(i3.3)') Atarget
  NMTlib = 0
  k = 0
  do ilib = 1, numlib
    if (libinclude(ilib) /= 1) cycle
    libpath = trim(librariespath)//ptype0//'/'// trim(Starget)//trim(massstring)//'/'//trim(lib(ilib))//'/tables/xs/'
!
! Read and select MT numbers of available library data
!
    filehead = trim(libpath)//ptype0//'-'//trim(Starget)// trim(massstring)//'-MT'
    do imt = 1, nummt
      MTstring = '000 '
      write(MTstring(1:3), '(i3.3)') imt
      do iso = -1, 2
        symiso = ' '
        if (iso == 0) symiso = 'g'
        if (iso == 1) symiso = 'm'
        if (iso == 2) symiso = 'n'
        if (iso >= 0) MTstring = MTstring(1:3) //symiso
        filename = trim(filehead)//trim(MTstring)//'.'//lib(ilib)
        inquire (file = filename, exist = lexist)
        if (lexist) then
          if (mtlib(imt, ilib) == 1) then
            k = k + 1
            MTlibk(k) = imt
            libname(k) = lib(ilib)
            isolib(k) = -1
            if (symiso == 'g') isolib(k) = 0
            if (symiso == 'm') isolib(k) = 1
            if (symiso == 'n') isolib(k) = 2
          endif
        endif
      enddo
    enddo
  enddo
  NMTlib = k
  if (NMTlib == 0) return
!
! Read library data per channel
!
  ee(0) = 0.
  xs(0) = 0.
  dxs(0) = 0.
  filehead = ptype0//'-'//trim(Starget)//trim(massstring)//'-MT'
  do i = 1, Nchanall
    imt = MT(i)
    if (imt == 0) cycle
    if (auto(imt) == 0) cycle
    MTstring = '000 '
    write(MTstring(1:3), '(i3.3)') imt
    if (MTiso(i) == 0) MTstring = MTstring(1:3)//'g'
    if (MTiso(i) == 1) MTstring = MTstring(1:3)//'m'
    if (MTiso(i) == 2) MTstring = MTstring(1:3)//'n'
    nexp = Nsets(i)
    nexp0 = Nsets(i)
    do k = 1, NMTlib
      if (MTlibk(k) /= imt .or. isolib(k) /= MTiso(i)) cycle
      mtfile = trim(filehead)//trim(MTstring)//'.'//libname(k)
      lfile = trim(librariespath)//ptype0//'/'// trim(Starget)//trim(massstring)//'/'//trim(libname(k))//'/tables/xs/'//mtfile
      inquire (file=lfile,exist=lexist)
      if (.not.lexist) then
        write(*,'(" TASMAN-error: ",a," can not be found: ")') trim(lfile)
        stop
      endif
      nexp = nexp + 1
      isoexp(i, nexp) = MTiso(i)
      libf0(nexp) = libname(k)
!
! Read and select library data from each set
!
      open (unit = 2, file = lfile, status = 'old')
      author = '                              '
      author(1:15) = libname(k)(1:15)
      do
        read(2,'(a)',iostat = istat) line
        if (istat == -1) exit
        j = 0
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(2,'(/)')
          do
            read(2, * , iostat = istat) ee0, xs0
            if (istat == -1) exit
            if (j > numexpen) exit
            if (ee0 < Esearch1(imt) .or. ee0 > Esearch2(imt)) cycle
            if (ee0 > 0.1 .and. abs(ee0 - ee(j)) <= 0.001) cycle
            j = j + 1
            ee(j) = ee0
            xs(j) = xs0
            dxs(j) = devi(imt) * xs0
          enddo
          exit
        endif
      enddo
      close (unit = 2)
      if (j == 0) cycle
      nenlfile = j
!
! Always take maximum of numexp (=1000) highest energy points
!
      nbeg = 0
      nend = 0
      do l = 1, nenlfile
        if (ee(l) >= Ein(1) .and. xs(l) > xseps) then
          nbeg = l
          exit
        endif
      enddo
      do l = nenlfile, 1, -1
        if (ee(l) <= Ein(numE) .and. xs(l) > xseps) then
          nend = l
          exit
        endif
      enddo
      if (nbeg == 0 .or. nend == 0) cycle
      j = 0
      if (nend - nbeg < numenexp) then
        do l = nbeg, nend
          j = j + 1
          Eexp(i, nexp, j) = ee(l)
          xsexp(i, nexp, j) = xs(l)
          dxsexp(i, nexp, j) = dxs(l)
        enddo
      else
        Eexp(i, nexp, 1) = ee(nbeg)
        xsexp(i, nexp, 1) = xs(nbeg)
        dxsexp(i, nexp, 1) = dxs(nbeg)
        Ea = ee(nbeg)
        Eb = ee(nend)
        deltaE = (log(Eb) - log(Ea)) / numenexp
        j = 1
        nprev = 0
        do n = 2, numenexp - 1
          En = exp(log(Ea) + n * deltaE)
          call locate(ee, nbeg, nend, En, nE)
          if (nE /= nprev) then
            j = j + 1
            Eexp(i, nexp, j) = ee(nE)
            xsexp(i, nexp, j) = xs(nE)
            dxsexp(i, nexp, j) = dxs(nE)
            nprev = nE
          endif
        enddo
        j = j + 1
        Eexp(i, nexp, j) = ee(nend)
        xsexp(i, nexp, j) = xs(nend)
        dxsexp(i, nexp, j) = dxs(nend)
      endif
      auth(i, nexp) = author
      Nenexp(i, nexp) = j
    enddo
    Nsets(i) = nexp
!
! Set weights for data library sets
!
Loop1:    do nexp = nexp0 + 1, Nsets(i)
      do ilib = 1, numlib
        if (lib(ilib) == libf0(nexp)) then
          Eweight(i, nexp) = mtlibweight(imt, ilib)
          cycle Loop1
        endif
      enddo
    enddo Loop1
  enddo
  return
end subroutine libread
! Copyright A.J. Koning 2021
