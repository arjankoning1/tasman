subroutine xsread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read cross sections from TALYS files
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
!   sgl             ! single precision kind
! All global variables
!   numchanxs       ! maximum number of channels with cross sections
!   numenin         ! maximum number of incident energies
!   numenS          ! maximum number of energies for sensitivities
!   nummt           ! maximum number of MT numbers
!   numpar          ! maximum number of parameters
! Variables for writing TALYS input
!   flagomponly     ! flag for covariances or optimization for basic optical model observables
!   italys          ! TALYS run counter
!   mode            ! TASMAN mode
!   Nhigh           ! number of high energy runs
! Variables for reading TALYS output
!   Ntalbeg         ! counter for first TALYS calculation
! Constants
!   particle        ! type of particle
! Variables for processing input
!   k0              ! index of incident particle
!   ptype0          ! type of incident particle
! Variables for weights
!   xsopt           ! optimal cross section
! Variables for sensitivity
!   flagsens        ! flag to use sensitivity in parameter weighting
! Variables for reading cross sections
!   dE              ! uncertainty of incident energy
!   E               ! incident energy
!   E1mb            ! energy of 1 mb value
!   Eendf           ! incident energy
!   Efracmax        ! starting energy for maximum cross section
!   Enon            ! incident energy
!   flagdiscrete    ! flag for discrete level cross sections
!   flagisomer      ! flag for isomeric production cross sections
!   MT              ! MT number
!   MTiso           ! isomer of MT number
!   MTnum           ! MT number of reaction channel
!   Nchanxs         ! total number of channels with cross sections
!   Nen             ! number of incident energies
!   Nen0            ! number of incident energies
!   Nenendf         ! number of incident energies
!   Nenendf0        ! number of incident energies
!   Nenlow          ! number of incident energies
!   Nennon          ! number of incident energies
!   Qval            ! Q - value
!   Sindex          ! index for cross section
!   xsendf          ! cross section from TALYS
!   xseval          ! evaluated cross section
!   xsfile          ! file with crosssections
!   xsmax           ! maximum cross section per channel
!   xsnon0          ! non - elastic crosssection
!   xssave          ! cross section from TALYS
!   xstalys         ! cross section from TALYS
!   xstalys2        ! cross section from TALYS
!   xstalys3        ! cross section from TALYS
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: numen6=12000       ! number of energies for ENDF-6 format
  logical            :: lexist             ! logical to determine existence
  character(len=12)  :: xsf                ! filename
  character(len=13)  :: xsfL               ! filename
  character(len=20)  :: xsf1               ! filename
  character(len=132) :: xsstring(numen6)   ! string for cross section
  character(len=132) :: line               ! line
  character(len=132) :: key
  integer            :: i                  ! counter
  integer            :: ia                 ! mass number of nucleus
  integer            :: ibeg               ! index to mark begin of word
  integer            :: ichan              ! counter for channels
  integer            :: id                 ! counter for deuterons
  integer            :: idc                ! help variable
  integer            :: ih                 ! hole number
  integer            :: imt                ! MT counter
  integer            :: imt0               ! MT counter
  integer            :: in                 ! counter for neutrons
  integer            :: ip                 ! particle number
  integer            :: isamp              ! variable for 0th or random run
  integer            :: iso                ! counter for isotope
  integer            :: istat              ! logical for file access
  integer            :: it                 ! counter for tritons
  integer            :: keyix
  integer            :: j                  ! counter
  integer            :: jheader            ! counter
  integer            :: k                  ! Legendre order
  integer            :: L                  ! counter for Legendre coefficients
  integer            :: n                  ! counter
  integer            :: nin                ! counter for incident energy
  integer            :: system             ! system call
  integer            :: type               ! particle type
  real(sgl)          :: cratE              ! cross section ratio
  real(sgl)          :: cratio             ! cross section ratio
  real(sgl)          :: cratio2            ! cross section ratio
  real(sgl)          :: cratio3            ! cross section ratio
  real(sgl)          :: denom              ! help variable
  real(sgl)          :: e1                 ! first energy
  real(sgl)          :: frac               ! help variable
  real(sgl)          :: Rs                 ! sensitivity factor
  real(sgl)          :: xs1                ! help variable
  real(sgl)          :: xs2                ! help variable
  real(sgl)          :: xs3                ! help variable
  real(sgl)          :: xsnew              ! new cross section
  real(sgl)          :: xsthresh           ! threshold cross section
  real(sgl)          :: xstop              ! maximum cross section
!
! ************************ Set channels ********************************
!
! Cross sections
!
  if (Nchanxs == 0) then
    ichan = 0
    inquire (file = 'total.tot', exist = lexist)
    if (lexist) then
      ichan = ichan + 1
      xsfile(ichan) = 'total.tot'
      MT(ichan) = 1
    endif
    if (flagomponly) then
      Nchanxs = 1
      goto 200
    else
      inquire (file = 'elastic.tot', exist = lexist)
      if (lexist) then
        ichan = ichan + 1
        xsfile(ichan) = 'elastic.tot'
        MT(ichan) = 2
      endif
      inquire (file = 'nonelastic.tot', exist = lexist)
      if (lexist) then
        ichan = ichan + 1
        xsfile(ichan) = 'nonelastic.tot'
        MT(ichan) = 3
      endif
      inquire (file = 'fission.tot', exist = lexist)
      if (lexist) then
        ichan = ichan + 1
        xsfile(ichan) = 'fission.tot'
        MT(ichan) = 18
      endif
      do ia = 0, 1
        do ih = 0, 1
          do it = 0, 1
            do id = 0, 1
              do ip = 0, 1
                do in = 0, 3
                  idc = 100000 * in + 10000 * ip + 1000 * id + 100 * it + 10 * ih + ia
                  do imt = 1, nummt
                    if (MTnum(imt) == idc) then
                      imt0 = imt
                      xsf = 'xs000000.tot'
                      write(xsf(3:8), '(i6.6)') idc
                      inquire (file = xsf, exist = lexist)
                      if (lexist) then
                        ichan = ichan + 1
                        if (ichan <= numchanxs) then
                          xsfile(ichan) = xsf
                          MT(ichan) = imt0
                        endif
                      endif
                      if (flagisomer) then
                        xsfL = xsf(1:9)//'L00'
                        iso = -1
                        do L = 0, 40
                          write(xsfL(11:12), '(i2.2)') L
                          inquire (file = xsfL, exist = lexist)
                          if (lexist) then
                            ichan = ichan + 1
                            if (ichan <= numchanxs) then
                              xsfile(ichan) = xsfL
                              iso = iso + 1
                              MT(ichan) = imt0
                              MTiso(ichan) = iso
                            endif
                          endif
                        enddo
                      endif
                      exit
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
! Discrete level scattering cross sections
!
      if (flagdiscrete) then
        do type = 0, 6
          ibeg = 0
          if (type == k0) ibeg = 1
          do i = ibeg, 40
            xsf = 'nn.L00      '
            write(xsf(1:1), '(a1)') ptype0
            write(xsf(2:2), '(a1)') particle(type)
            write(xsf(5:6), '(i2.2)') i
            xsf1 = xsf(1:6)
            inquire (file = xsf1, exist = lexist)
            if (lexist) then
              ichan = ichan + 1
              if (ichan <= numchanxs) then
                xsfile(ichan) = xsf1
                if (type == 1) then
                  MT(ichan) = 50 + i
                else
                  MT(ichan) = 500 + type * 50 + i
                endif
              endif
            endif
          enddo
          xsf = 'nn.con      '
          write(xsf(2:2), '(a1)') particle(type)
          xsf1 = xsf(1:6)
          inquire (file = xsf1, exist = lexist)
          if (lexist) then
            ichan = ichan + 1
            if (ichan <= numchanxs) then
              xsfile(ichan) = xsf1
              if (type == 1) then
                MT(ichan) = 91
              else
                MT(ichan) = 549 + type * 50
              endif
            endif
          endif
        enddo
      endif
      Nchanxs = min(ichan, numchanxs)
    endif
  endif
!
! ************************** Read cross sections ***********************
!
200 if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  do i = 1, Nchanxs
    if (italys > Nhigh) then
      do k = 1, Nen0(i)
        xstalys(isamp, i, k) = xstalys(0, i, k)
        xstalys2(isamp, i, k) = xstalys2(0, i, k)
        xstalys3(isamp, i, k) = xstalys3(0, i, k)
        xseval(isamp, i, k) = xseval(0, i, k)
      enddo
    endif
    inquire (file = xsfile(i), exist = lexist)
    if (lexist) then
      open (unit = 2, file = xsfile(i), status = 'old', iostat = istat)
      if (istat /= 0) call read_error(xsfile(i), istat)
      do
        read(2,'(a)',iostat = istat) line
        if (istat == -1) exit
        E(i, 0) = 0.
        k = 0
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(2,'(/)',iostat = istat)
          if (istat == -1) exit
          do
            read(2, '(4es15.6)', iostat = istat) e1, xs1, xs2, xs3
            if (istat == -1) exit
            if (istat /= 0) call read_error(xsfile(i), istat, eor = 'continue')
            k = k + 1
            if (k > numenin) then
              write(*, '(" TASMAN-error: Number of incident energies larger than numenin =",i6)') numenin
              stop
            endif
            E(i, k) = e1
            xstalys(isamp, i, k) = xs1
            xstalys2(isamp, i, k) = xs2
            xstalys3(isamp, i, k) = xs3
            xseval(isamp, i, k) = xs1
          enddo
          exit
        endif
      enddo
      close (2)
      Nen(i) = k
      if (italys <= Nhigh) then
        Nen0(i) = k
      else
        if (xstalys(0, i, k) > 0.) then
          cratio = xstalys(isamp, i, k) / xstalys(0, i, k)
        else
          cratio = 1.
        endif
        cratio = min(cratio, 10.)
        if (xstalys2(0, i, k) > 0.) then
          cratio2 = xstalys2(isamp, i, k) / xstalys2(0, i, k)
        else
          cratio2 = 1.
        endif
        cratio2 = min(cratio2, 10.)
        if (xstalys3(0, i, k) > 0.) then
          cratio3 = xstalys3(isamp, i, k) / xstalys3(0, i, k)
        else
          cratio3 = 1.
        endif
        cratio3 = min(cratio3, 10.)
        do j = k + 1, Nen0(i)
          xstalys(isamp, i, j) = cratio * xstalys(isamp, i, j)
          xstalys2(isamp, i, j) = cratio2 * xstalys2(isamp, i, j)
          xstalys3(isamp, i, j) = cratio3 * xstalys3(isamp, i, j)
          xseval(isamp, i, j) = cratio * xseval(isamp, i, j)
        enddo
        open (unit = 2, file = xsfile(i), status = 'old', iostat = istat)
        if (istat /= 0) call read_error(xsfile(i), istat)
        jheader = 0
        j = 1
        do
          read(2, '(a)', iostat = istat) xsstring(j)
          if (istat == -1) exit
          if (istat /= 0) call read_error(xsfile(i), istat)
          key='entries'
          keyix=index(xsstring(j),trim(key))
          if (keyix > 0) then
            write(xsstring(j)(keyix+len_trim(key)+2:80), '(i6)') Nen0(i)
            jheader = j + 2
          endif
          j = j + 1
        enddo
        close (2)
        open (unit = 2, file = xsfile(i), status = 'replace')
        do j = 1, jheader + Nen(i) 
          write(2, '(a)') trim(xsstring(j))
        enddo
        do j = Nen(i) + 1, Nen0(i)
          write(2, '(4es15.6)') E(i, j), xseval(isamp, i, j), xstalys2(isamp, i, j), xstalys3(isamp, i, j)
        enddo
        close (2)
        Nenlow(i) = Nen(i)
        Nen(i) = Nen0(i)
      endif
      dE(i, 0) = 0.5 * (E(i, 0) + E(i, 1))
      do k = 1, Nen(i) - 1
        dE(i, k) = 0.5 * (E(i, k + 1) - E(i, k - 1))
      enddo
      dE(i, Nen(i)) = 0.5 * (E(i, Nen(i)) - E(i, Nen(i) - 1))
    endif
  enddo
  if (italys > Nhigh) then
    do i = 1, 3
      do k = 1, Nenendf0
        xsendf(isamp, i, k) = xsendf(0, i, k)
      enddo
    enddo
  endif
  inquire (file = 'endf.tot', exist = lexist)
  if (lexist) then
    open (unit = 2, file = 'endf.tot', status = 'old', iostat = istat)
    if (istat /= 0) call read_error('endf.tot', istat)
    do
      read(2,'(a)',iostat = istat) line
      if (istat /= 0) call read_error('endf.tot', istat)
      Eendf(0) = 0.
      k = 0
      key='entries'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(2,'(/)',iostat = istat) 
        if (istat /= 0) call read_error('endf.tot', istat)
        do
          read(2, '(4es15.6)', iostat = istat) e1, xs1, xs2, xs3
          if (istat == -1) exit
          if (istat /= 0) call read_error('endf.tot', istat, eor = 'continue')
          k = k + 1
          Eendf(k) = e1
          xsendf(isamp, 1, k) = xs1
          xsendf(isamp, 2, k) = xs2
          xsendf(isamp, 3, k) = xs3
        enddo
        exit
      endif
    enddo
    close (2)
    Nenendf = k
    if (italys <= Nhigh) then
      Nenendf0 = k
    else
      do i = 1, 3
        if (xsendf(0, i, k) > 0.) then
          cratE = xsendf(isamp, i, k) / xsendf(0, i, k)
        else
          cratE = 1.
        endif
        cratE = min(cratE, 10.)
        do j = k + 1, Nenendf0
          xsendf(isamp, i, j) = cratE * xsendf(isamp, i, j)
        enddo
      enddo
      open (unit = 2, file = 'endf.tot', status = 'old', iostat = istat)
      if (istat /= 0) call read_error('endf.tot', istat)
      jheader = 0
      j = 1
      do
        read(2, '(a)', iostat = istat) xsstring(j)
        if (istat == -1) exit
        if (istat /= 0) call read_error('endf.tot', istat)
        key='entries'
        keyix=index(xsstring(j),trim(key))
        if (keyix > 0) then
          write(xsstring(j)(keyix+len_trim(key)+2:80), '(i6)') Nenendf0
          jheader = j + 2
        endif
        j = j + 1
      enddo
      close (2)
      open (unit = 2, file = 'endf.tot', status = 'replace')
      do j = 1, jheader + Nenendf 
        write(2, '(a)') trim(xsstring(j))
      enddo
      do j = Nenendf + 1, Nenendf0
        write(2, '(4es15.6)') Eendf(j), (xsendf(isamp, i, j), i = 1, 3)
      enddo
      close (unit = 2)
      Nenendf = Nenendf0
    endif
  endif
  if (italys > Nhigh) then
    inquire (file = 'all.tot', exist = lexist)
    if (lexist) then
      open (unit = 2, file = 'all.tot', status = 'old', iostat = istat)
      if (istat /= 0) call read_error('all.tot', istat)
      jheader = 0 
      j = 1 
      do 
        read(2, '(a)', iostat = istat) xsstring(j)
        if (istat == -1) exit
        if (istat /= 0) call read_error('all.tot', istat)
        key='entries'
        keyix=index(xsstring(j),trim(key))
        if (keyix > 0) then
          write(xsstring(j)(keyix+len_trim(key)+2:80), '(i6)') Nen0(1)
          jheader = j + 2
        endif
        j=j+1
      enddo
      close (unit = 2)
      open (unit = 2, file = 'all.tot', status = 'replace')
      do j = 1, jheader + Nenlow(1) 
        write(2, '(a)') trim(xsstring(j))
      enddo
      do j = Nenlow(1) + 1, Nen0(1)
        write(2, '(4es15.6)') E(1, j), (xseval(isamp, i, j), i = 3, 1, -1)
       enddo
      close (unit = 2)
    endif
  endif
  inquire (file = 'nonelastic.tot', exist = lexist)
  if (lexist) i = system('mv -f nonelastic.tot nonelastic.000')
!
! Store non-elastic cross section in special array
!
  if (italys == 0) then
    if (k0 == 1) then
      nin = 3
    else
      nin = 2
    endif
    Nennon = Nen(nin)
    do k = 1, Nennon
      Enon(k) = E(nin, k)
      xsnon0(k) = xstalys(0, nin, k)
    enddo
    Enon(0) = 0.
    xsnon0(0) = xsnon0(1)
!
! Determine maximum cross section for each channel (used for deviation criteria).
! Determine acceptable deviation for each channel.
!
! deviation : subroutine to set acceptable deviations from average
!
    do i = 1, Nchanxs
      xstop = 0.
      do k = 1, Nen(i)
        if (E(i, k) <= Efracmax(MT(i))) cycle
        xsnew = xstalys(0, i, k)
        if (xsnew > xstop) xstop = xsnew
      enddo
      xsmax(i) = xstop
    enddo
!
! Set optimum cross sections to initial set
!
    do i = 1, Nchanxs
      do k = 1, Nen(i)
        xsopt(i, k) = xstalys(0, i, k)
      enddo
    enddo
!
! Read Q-value
!
    Qval = 0.
    inquire (file = 'xs000000.tot', exist = lexist)
    if (lexist) then
      open (unit = 2, file = 'xs000000.tot', status = 'old', iostat = istat)
      if (istat /= 0) call read_error('xs000000.tot', istat)
      do
        read(2, '(a)', iostat = istat) line
        if (istat == -1) exit
        if (istat /= 0) call read_error('xs000000.tot', istat)
        key='Q-value [MeV]'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(line(keyix+len_trim(key)+2:80), *) Qval
          exit
        endif
      enddo
      close (2)
    endif
  endif
!
! Search for energy where cross section crosses the threshold value of 1 mb.
!
  if (italys == Ntalbeg) then
    do i = 1, Nchanxs
      E1mb(i) = 0.
      if (MT(i) > 3 .and. MT(i) /= 102) then
        xsthresh = 1.
        do k = 1, Nen(i) - 1
          if (E(i, k) > 1..and.xstalys(0, i, k) <= xsthresh .and. xstalys(0, i, k + 1) > xsthresh) then
            frac = (xsthresh - xstalys(0, i, k)) / (xstalys(0, i, k + 1) - xstalys(0, i, k))
            denom = E(i, k + 1) - E(i, k)
            if (denom > 0.) E1mb(i) = E(i, k) + frac * denom
          endif
        enddo
      endif
    enddo
  endif
!
! Create energy grid and cross sections for sensitivities
!
  if (flagsens .or. mode == 2) then
    if (italys == 0) then
      do i = 1, Nchanxs
        if (Nen(i) <= numenS) then
          do k = 1, Nen(i)
            Sindex(i, k) = k
          enddo
        else
          Rs = real(numenS - 1) / real(Nen(i)) + 0.001
          do k = 1, Nen(i)
            Sindex(i, k) = int(1 + k * Rs)
          enddo
        endif
      enddo
    endif
    if (italys >= 0 .and. italys <= numpar) then
      do i = 1, Nchanxs
        do k = 1, Nen(i)
          n = Sindex(i, k)
          xssave(italys, i, n) = xseval(isamp, i, k)
        enddo
      enddo
    endif
  endif
  do i = 1, Nchanxs
    reaction_string(i) = MTreac(MT(i))
    if (MTiso(i) == 0) reaction_string(i) = trim(reaction_string(i))//'g'
    if (MTiso(i) == 1) reaction_string(i) = trim(reaction_string(i))//'m'
    if (MTiso(i) == 2) reaction_string(i) = trim(reaction_string(i))//'n'
    MF(i) = 3
    if (MTiso(i) >= 0) then
      if (MT(i) == 102) then
        MF(i) = 9
      else
        MF(i) = 10
      endif
    endif
  enddo
  return
end subroutine xsread
! Copyright A.J. Koning 2021
