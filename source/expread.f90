subroutine expread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read experimental data
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
! Variables for path names
!   psfpath          ! directory containing files to be read
! All global variables
!   numchanxs         ! maximum number of channels with cross sections
!   numenexp          ! maximum number of energies forexperimental data set
!   nummt             ! maximum number of MT numbers
!   numsets           ! maximum number of experimental data sets
! Constants
!   xseps             ! limit for cross sections
! Variables for reading TALYS output
!   flagpsf           ! flag for optimization of photon strength functions
!   Nchanall          ! total number of channels
! Variables for processing input
!   Atarget           ! mass number of target nucleus
!   errlim            ! lower limit for experimental error
!   exppath           ! path for experimental data
!   mtallexp          ! weight for all experiments per channel
!   numE              ! number of incident energies
!   Ztarget           ! charge number of target nucleus
! Variables for reading cross sections
!   MT                ! MT number
!   MTiso             ! isomer of MT number
! Variables for reading experimental data
!   auth              ! author
!   dxsexp            ! uncertainty of experimental cross section
!   Eexp              ! incident energy
!   Ein               ! incident energy
!   Esearch1          ! start energy of search
!   Esearch2          ! end energy of search
!   Eweight           ! incident energy
!   expclass          ! quality class of experimental data (1, 2 or 3)
!   expexccount       ! counter for excluded experimental data sets
!   expexcfile        ! file with experimental data set to exclude for search
!   expinccount       ! counter for included experimental data sets
!   expincfile        ! file with experimental data set to include for search
!   flagscore         ! flag to use scores for inclusion or exclusion
!   flagerror         ! flag to include true experimental uncertainties
!   flagsort          ! flag to sort experimental data in increasing energies
!   isoexp            ! isomer of experimental data
!   maxexpsets        ! number of experimental data sets before weight reduction
!   mtallexpweight    ! weight for all channels per experimental channel
!   MTexp             ! MT number of experimental data
!   mtexpweight       ! weight for experimental data set and channel
!   Nenexp            ! number of experimental data points
!   NMTexp            ! number of MT numbers for experimental data
!   Nsets             ! number of experimental data sets
!   quality           ! quality class
!   subentry          ! subentry
!   xsexp             ! experimental cross section
!   expyear           ! year
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: numexpen=150000              ! number of energies for experimental data
  logical            :: adopt                        ! logical for existence of MF information (per MT)
  logical            :: lexist                       ! logical to determine existence
  character(len=2)   :: qual                         ! quality flag
  character(len=4)   :: mtstring                     ! string with reaction information
  character(len=7)   :: ZAstring                     !
  character(len=9)   :: suben                        ! subentry
  character(len=9)   :: sube                         ! subentry
  character(len=9)   :: subent(numchanxs*numsets)    ! subentry
  character(len=30)  :: author                       ! author
  character(len=132) :: word(40)                     ! words on input line
  character(len=132) :: expfi                        ! filename for experimental data
  character(len=132) :: expfile(numchanxs*numsets)   ! experimental data file
  character(len=132) :: expreac                      ! experimental data path name
  character(len=132) :: expf                         ! experimental data file
  character(len=132) :: expxs                        ! experimental cross section
  character(len=132) :: str                          ! input line
  character(len=132) :: psflist                      !
  character(len=132) :: scorefile                    !
  character(len=132) :: xslist                       ! file with cross sections
  integer            :: class                        ! input class
  integer            :: expiso(numchanxs*numsets)    ! isomeric cross section
  integer            :: i                            ! counter
  integer            :: ii                           ! counter
  integer            :: imt                          ! MT counter
  integer            :: iso                          ! counter for isotope
  integer            :: istat                        ! logical for file access
  integer            :: j                            ! counter
  integer            :: k                            ! Legendre order
  integer            :: kend                         ! index to mark end of word
  integer            :: l                            ! counter
  integer            :: n                            ! counter
  integer            :: nbeg                         ! first outgoing energy
  integer            :: nE                           ! number of energies
  integer            :: nend                         ! help variable
  integer            :: nenexpfile                   ! number of incident energies in experimental data file
  integer            :: nexp                         ! number of experimental data sets
  integer            :: nexpinc                      ! number opf experimental data sets
  integer            :: NMT(nummt)                   ! total number of MT sections
  integer            :: NMTexp0                      ! numner of MT sections
  integer            :: nprev                        ! help variable
  integer            :: yr                           ! help array for storing segment intersection points  with the bin sides
  integer            :: iw                           !
  real(sgl)          :: de0                          !
  real(sgl)          :: deltaE                       ! energy bin around outgoing energies
  real(sgl)          :: dxs(0:numexpen)              ! uncertainty of experimental cross section
  real(sgl)          :: dxs0                         ! uncertainty cross section
  real(sgl)          :: dxsexptmp                    ! uncertainty of experimental cross section
  real(sgl)          :: Ea                           ! parameters for local energy dependence
  real(sgl)          :: Eb                           ! end energy of local adjustment
  real(sgl)          :: ee(0:numexpen)               ! incident energy
  real(sgl)          :: ee0                          ! incident energy
  real(sgl)          :: Eexptmp                      ! energy for experimental cross section
  real(sgl)          :: En                           ! outgoing energy  term03, conv03, enhance03 : help variable for breakup enhanc
  real(sgl)          :: errsum                       ! sum of uncertainties
  real(sgl)          :: factor                       ! multiplication factor
  real(sgl)          :: relerr                       ! relative cross section uncertainty
  real(sgl)          :: ww                           ! weight
  real(sgl)          :: xs(0:numexpen)               ! experimental cross section
  real(sgl)          :: xs0                          ! cross section
  real(sgl)          :: xsexptmp                     ! experimental cross section
!
! ************************ Read experimental data **********************
!
  expxs = exppath
  NMT = 0
  scorefile = trim(tasmanpath)//'misc/score.tab'
!
! Read and select MT numbers of available experimental data
!
  if (flagcross) then
    i = 0
    xslist = trim(expxs)//'xslist'
    inquire (file = xslist, exist = lexist)
    if (lexist) then
      open (unit = 3, file = xslist, status = 'old', iostat = istat)
      if (istat /= 0) call read_error(xslist, istat)
      do
        read(3, '(a132)', iostat = istat) str
        if (istat == -1) exit
        if (istat /= 0) call read_error(xslist, istat)
        call getkeywords(str, word)
        i = i + 1
        read(word(1), * ) imt
        MTexp(i) = imt
        expfile(i) = word(2)
        subent(i) = word(3)(1:9)
        NMT(imt) = NMT(imt) + 1
      enddo
      close (unit = 3)
      NMTexp0 = i
    else
      write(*,'(" TASMAN-warning: ",a," does not exist")') trim(xslist)
    endif
    do imt = 1, nummt
Loop1: do j = 1, expinccount(imt)
        do ii = 1, NMTexp0
          if (expfile(ii) == expincfile(imt, j) .or. subent(ii)(1:9) == expincfile(imt, j)(1:9)) cycle Loop1
        enddo
        i = i + 1
        NMT(imt) = NMT(imt) + 1
        MTexp(i) = imt
        expfile(i) = expincfile(imt, j)
        subent(i) = '         '
      enddo Loop1
    enddo
    NMTexp = i
  endif
  if (flagpsf) then
    imt = 92
    ZAstring = '000_000'
    write(ZAstring(1:3), '(i3.3)') Ztarget
    write(ZAstring(5:7), '(i3.3)') Atarget
    psflist = trim(psfpath)//'psflist'
    inquire (file = psflist, exist = lexist)
    if (lexist) then
      open (unit = 3, file = psflist, status = 'old')
      i = 0
      do
        read(3, '(a132)', iostat = istat) str
        if (istat == -1) exit
        if (istat /= 0) call read_error(xslist, istat)
        k = index(str, ZAstring)
        l = index(str, 'photoneut')
        if (k > 0 .and. l > 0) then
          i = i + 1
          expfile(i) = trim(psfpath) //str
          NMT(imt) = NMT(imt) + 1
          MTexp(i) = 92
        endif
      enddo
    else 
      write(*,'(" TASMAN-error: ",a," does not exist")') trim(psflist)
      stop
    endif
    MT(1) = 92
    NMTexp = NMT(imt)
    NMTexp0 = NMT(imt)
  endif
Loop2:  do i = 1, NMTexp
    expiso(i) = -1
    do k = 2, 131
      if (expfile(i)(k:k+1) == 'MT') then
        kend = k + 5
        iso = -1
        if (expfile(i)(kend:kend) == 'g') iso = 0
        if (expfile(i)(kend:kend) == 'm') iso = 1
        if (expfile(i)(kend:kend) == 'n') iso = 2
        expiso(i) = iso
        cycle Loop2
      endif
    enddo
  enddo Loop2
  ee(0) = 0.
  xs(0) = 0.
  dxs(0) = 0.
  if (flagscore) open(unit=17,file='score.out',status='replace')
  do i = 1, Nchanall
    imt = MT(i)
    if (imt == 0) cycle
    if (NMT(imt) == 0) cycle
    mtstring = '000 '
    write(mtstring(1:3), '(i3.3)') imt
    if (MTiso(i) == 0) mtstring = mtstring(1:3)//'g'
    if (MTiso(i) == 1) mtstring = mtstring(1:3)//'m'
    if (MTiso(i) == 2) mtstring = mtstring(1:3)//'n'
    expreac = trim(expxs)//trim(mtstring)//'/'
    nexp = Nsets(i)
    nexpinc = nexp
    do ii = 1, NMTexp
      if (MTexp(ii) /= imt) cycle
      if (expiso(ii) /= MTiso(i)) cycle
      expf = trim(expfile(ii))
      if (ii <= NMTexp0) then
        expfi = trim(expreac) //expfile(ii)
      else
        expfi = expfile(ii)
      endif
      if (flagpsf) expfi = expfile(ii)
!
! Read and select experimental data sets for this reaction
!
      adopt = .false.
      if (mtallexp(imt) > 0) adopt = .true.
      ww = mtallexpweight(imt)
      do j = 1, expinccount(imt)
        if (trim(expf) == trim(expincfile(imt, j)) .or. subent(ii)(1:9) == expincfile(imt, j)(1:9)) then
          adopt = .true.
          ww = mtexpweight(imt, j)
        endif
      enddo
      do j = 1, expexccount(imt)
        if (trim(expf) == trim(expexcfile(imt, j))) adopt = .false.
      enddo
      if (flagpsf) adopt = .true.
      if ( .not. adopt) cycle
      inquire (file = expfi, exist = lexist)
      if ( .not. lexist) then
        write(*, '(" TASMAN-warning: Wrong data file: ", a)') trim(expfi)
      else
        open (unit = 2, file = expfi, status = 'old')
        errsum = 0.
        author = ''
        suben = '         '
        qual = '  '
        class = 0
        yr = 0
        k = 0
        if (flagpsf) then
          do
            read(2, '(a132)', iostat = istat) str
            if (istat == -1) exit
            if (istat /= 0) call read_error(xslist, istat)
!           if (str(1:8) == '# Author') author = str(11:40)
            if (str(1:10) == '#   author' .and. author == '') author = str(13:40)
            if (str(1:1) /= '#') then
              if (k > numexpen) cycle
              read(str, * , iostat = istat) ee0, de0, xs0, dxs0
              if (istat == -1) exit
              if (istat /= 0) call read_error(xslist, istat)
              k = k + 1
              ee(k) = ee0
              xs(k) = xs0
              dxs(k) = dxs0
            endif
          enddo
        else
          do
            read(2, '(a132)', iostat = istat) str
            if (istat == -1) exit
!           if (str(1:16) == '# X4 Subentry : ') suben = str(17:25)
            if (str(1:12) == '#   subentry') suben = str(15:27)
            if (ii <= NMTexp0) then
!             if (str(1:16) == '# Author      : ') author = str(17:45)
              if (str(1:10) == '#   author' .and. author == '') author = str(13:45)
            else
              author = trim(expf)
            endif
!           if (str(1:16) == '# Year        : ') read(str(17:20), '(i4)') yr
            if (str(1:8) == '#   year') read(str(11:20), '(i4)') yr
!           if (str(1:9) == '# Quality' .and. str(21:21) == ':') then
            if (str(1:11) == '#   quality') then
!             qual = str(23:24)
              qual = str(14:15)
              if (qual(2:2) >= '1' .and. qual(2:2) <= '3') read(qual(2:2), '(i1)') class
            endif
            if (str(1:1) /= '#') then
              if (k > numexpen) cycle
              read(str, * , iostat = istat) ee0, de0, xs0
              if (istat == -1) exit
              if (istat /= 0) cycle
              read(str, * , iostat = istat) ee0, de0, xs0, dxs0
              if (istat == 0) then
                if (ee0 < Esearch1(imt) .or. ee0 > Esearch2(imt)) cycle
                k = k + 1
                ee(k) = ee0
                xs(k) = xs0
              endif
              if (flagerror .and. xs0 > 0.) then
                relerr = dxs0 / xs0
                if (relerr < errlim(imt)) then
                  relerr = errlim(imt)
                  dxs0 = relerr * xs0
                endif
                errsum = errsum + relerr
                dxs(k) = dxs0
              else
                dxs(k) = 0.10 * xs0
              endif
            endif
          enddo
        endif
        close (unit = 2)
      endif
      if (k == 0) cycle
      nenexpfile = k
!
! Always take maximum of numexp (=1000) highest energy points
!
      if (flagpsf) then
        nbeg = 1
        nend = nenexpfile
      else
        nbeg = 0
        nend = 0
      endif
      do l = 1, nenexpfile
        if (ee(l) >= Ein(1) .and. xs(l) > xseps) then
          nbeg = l
          exit
        endif
      enddo
      do l = nenexpfile, 1, -1
        if (ee(l) <= Ein(numE) .and. xs(l) > xseps) then
          nend = l
          exit
        endif
      enddo
      if (nbeg == 0 .or. nend == 0) cycle
      j = 0
      nexp = nexp + 1
      if (nexp > numsets) exit
      Nsets(i) = nexp
      isoexp(i, nexp) = expiso(ii)
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
      expyear(i, nexp) = yr
      subentry(i, nexp) = suben
      quality(i, nexp) = qual
      Nenexp(i, nexp) = j
      if (class > expclass) then
        Eweight(i, nexp) = 0.
      else
        Eweight(i, nexp) = ww
        nexpinc = nexpinc + 1
      endif
      if (flagscore) then
        open (unit=18, file=scorefile, status='old')
        do
          read(18, '(a9,i2)', iostat = istat) sube, iw
          if (istat == -1) exit
          if (suben == sube) then
            Eweight(i, nexp) = real(iw)
            if (iw > 0) nexpinc = nexpinc + 1
            exit
          endif
        enddo
        close(18)
        write(17, '(a9,i2," # AK ",i4,"-",i2.2,"-",i2.2,1x,a1,2i4,1x,a4,1x,a,1x,i4)') suben,iw,year,month,day,ptype0, &
 &        Ztarget,Atarget,mtstring,trim(author),yr
      endif
      if (flagerror) then
        do k = 1, Nenexp(i, nexp)
          if (dxsexp(i, nexp, k) == 0.) then
            if (errsum == 0.) then
              dxsexp(i, nexp, k) = 0.10 * xsexp(i, nexp, k)
            else
              dxsexp(i, nexp, k) = errsum / Nenexp(i, nexp) * xsexp(i, nexp, k)
            endif
          endif
        enddo
      endif
!
! Sort experimental data in increasing energy, if necessary
!
      if (flagsort) then
        do k = 1, Nenexp(i, nexp)
          do l = 1, k
            if (Eexp(i, nexp, k) < Eexp(i, nexp, l)) then
              Eexptmp = Eexp(i, nexp, l)
              xsexptmp = xsexp(i, nexp, l)
              dxsexptmp = dxsexp(i, nexp, l)
              Eexp(i, nexp, l) = Eexp(i, nexp, k)
              xsexp(i, nexp, l) = xsexp(i, nexp, k)
              dxsexp(i, nexp, l) = dxsexp(i, nexp, k)
              Eexp(i, nexp, k) = Eexptmp
              xsexp(i, nexp, k) = xsexptmp
              dxsexp(i, nexp, k) = dxsexptmp
            endif
          enddo
        enddo
      endif
    enddo
!
! Option to limit the total weight of many experimental data sets
!
    if (nexpinc > maxexpsets) then
      factor = real(maxexpsets) / nexpinc
      do j = 1, nexp
        Eweight(i, j) = Eweight(i, j) * factor
      enddo
    endif
  enddo
  if (flagscore) close(unit=17)
  return
end subroutine expread
! Copyright A.J. Koning 2021
