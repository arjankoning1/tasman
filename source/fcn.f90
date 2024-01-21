subroutine fcn(N, P, G)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Goodness-of-fit estimator to be optimized
!
! Author    : Arjan Koning
!
! 2022-06-08: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
!
! Definition of single and double precision variables
!   sgl             ! single precision kind
!   dbl             ! double precision kind
! All global variables
!   numchanxs       ! maximum number of channels with cross sections
!   numenexp        ! maximum number of energies forexperimental data set
!   numenin         ! maximum number of incident energies
!   numsets         ! maximum number of experimental data sets
! Variables for reading TALYS output
!   flagcross       ! flag for covariances or optimization of cross sections
!   flagtal         ! flag to optimize to TALYS
!   flagtalys       ! flag to track successful TALYS calculation
!   Nchanall        ! total number of channels
! Variables for writing TALYS input files
!   A               ! mass number of nucleus
!   first           ! tracks whether a TALYS run is the first one
!   flagecis        ! flag to repeat ECIS calculations in a second run
!   flagompvar      ! flag for variation of OMP parameters
!   inuc            ! index of nuclide
!   italys          ! TALYS run counter
!   mode            ! TASMAN mode
!   Z               ! charge number of nucleus
! Variables for processing input
!   Atarget         ! mass number of target nucleus
!   isochar         ! symbol for isomer
!   talys           ! TALYS executable
!   Ztarget         ! charge number of target nucleus
! Variables for experimental data files
!   Nnuc            ! number of nuclides for which experimental data exists
! Variables for parameter variation
!   flagsample      ! flag for sampling experimental data
!   par             ! parameter value
! Variables for reading cross sections
!   fracmax         ! fraction of maximum cross section to be included
!   MT              ! MT number
!   MTiso           ! isomer of MT number
!   Nchanxs         ! total number of channels with cross sections
!   Nen             ! number of incident energies
!   xseval          ! evaluated cross section
!   xsfile          ! file with crosssections
!   xsmax           ! maximum cross section per channel
!   xsnon           ! nonelastic cross section
!   xsth            ! theoretical cross section
! Variables for weights
!   Gchannelsave    ! GOF for channel
!   Gsave           ! GOF value
!   xsopt           ! optimal cross section
! Variables for reading experimental data
!   auth            ! author
!   dEexp           ! uncertainty of incident energy
!   dxsexp          ! uncertainty of experimental cross section
!   Ebeg            ! first energy of experimental energy grid
!   Eend            ! last energy of experimental energy grid
!   Eexp            ! incident energy
!   Esearch1        ! start energy of search
!   Esearch2        ! end energy of search
!   Eweight         ! incident energy
!   isoexp          ! isomer of experimental data
!   isolib          ! isomer of MT number of library
!   MTexp           ! MT number of experimental data
!   MTlibk          ! MT number of library
!   NMTexp          ! number of MT numbers for experimental data
!   NMTlib          ! number of MT numbers of library
!   Nenexp          ! number of experimental data points
!   Nsets           ! number of experimental data sets
!   quality         ! quality class
!   subentry        ! subentry
!   xsexp           ! experimental cross section
!   xsexpsamp       ! sampled experimental cross section
!   expyear         ! year
! Variables for GOF function
!   chi2max         ! maximal chi2 value per point taken into account
!   Eindex          ! index for experimental data
!   flagdeltaE      ! flag to weigh goodness - of - fit with dE of energy grid
!   Fmax            ! maximal F value per point taken into account
!   Gall            ! GOF value
!   Gall0           ! GOF value
!   Gchannel        ! GOF value of channel
!   Grun            ! GOF value of run
!   Grun0           ! GOF value of run
!   ichi2           ! parameter to set goodness - of - fit estimator
!   isearch         ! number of trial run
!   Nindex          ! index for experimental data
!   Ntalys          ! number of TALYS runs
!   outsearch       ! level of output of search (1 (minimal) - 5 (all))
!   power           ! power for goodness - of - fit estimator K
!
! *** Declaration of local data
!
  implicit none

  integer, parameter :: numtype=3                                   ! number of GOF types
  character(len=1),save :: chinclude(numchanxs, numsets, 0:numenexp)! flag to include channel
  character(len=6)   :: ext                                         ! filename extension
  character(len=8)   :: goffile                                     ! file with GOF results
  character(len=10)  :: Gstring(numtype)                            ! string for G estimator
  character(len=132) :: cmd                                         ! command
  integer            :: i                                           ! counter
  integer            :: ia                                          ! mass number of nucleus
  integer            :: ifile                                       ! file
  integer            :: imt                                         ! MT counter
  integer            :: iz                                          ! charge number of nucleus
  integer            :: j                                           ! counter
  integer            :: k                                           ! Legendre order
  integer            :: m                                           ! counter
  integer            :: N                                           ! number of parameters
  integer            :: Nchannuc                                    ! number of channels per nucleus
  integer            :: Nchantot                                    ! total number of channels
  integer            :: nenchan(numchanxs)                          ! number of energies per channel
  integer            :: nenchanorg(numchanxs)                       ! number of energies per channel
  integer            :: nennuc                                      ! total number of energies per nucleus
  integer            :: nennucorg                                   ! total number of energies per nucleus
  integer            :: nenset(numsets)                             ! number of energies per set
  integer            :: nentot                                      ! total number of energies
  integer            :: Nnuctot                                     ! total number of channels for all nuclides
  integer, save      :: Nopt(numchanxs)                             ! number of optimum run
  integer, save      :: Nopttot(numchanxs)                          ! number of optimum run for total
  integer            :: Nsetchan(numchanxs)                         ! number of sets per channel
  integer            :: Nsetnuc                                     ! number of sets per nucleus
  integer            :: Nsettot                                     ! total number of sets
  integer            :: system                                      ! system call
  integer            :: ir(numtalys)                                ! 
  integer            :: Nr                                          ! 
  integer            :: istat                                       ! 
  real(sgl)          :: deltaE                                      ! energy bin around outgoing energies
  real(sgl)          :: expo                                        ! exponent
  real(sgl)          :: F0                                          ! GOF estimator
  real(sgl)          :: Fset(numsets, numtype)                      ! F value of data set
  real(sgl)          :: G                                           ! goodness-of-fit value  : G = 1: chi-2 : G = 2: squared distanc
  real(sgl)          :: Gchan(numtype)                              ! GOF for channel
  real(sgl)          :: Gchanav(numchanxs,numtype)                  ! average GOF for channel
  real(sgl), save    :: Gchanopt(numchanxs)                         ! optimal GOF for channel
  real(sgl), save    :: Gchanopttot(numchanxs)                      ! total optimal GOF for channel
  real(sgl), save    :: Gmin(0:numtype)                             ! minimal GOF
  real(sgl)          :: Gnuc(numtype)                               ! GOF for nucleus
  real(sgl)          :: Gnucav(numtype)                             ! average GOF for nucleus
  real(sgl)          :: Gp(numchanxs, numsets, 0:numenexp)          ! GOF per point
  real(sgl)          :: Gpoint(numtype)                             ! GOF per point
  real(sgl)          :: Gset(numsets, numtype)                      ! GOF for data set
  real(sgl)          :: Gsetav(numsets, numtype)                    ! average GOF for data set
  real(sgl)          :: Gtot(numtype)                               ! total GOF
  real(sgl)          :: Gtotav(numtype)                             ! average total GOF
  real(sgl)          :: P(N)                                        ! sampled parameter
  real(sgl)          :: relerr                                      ! relative cross section uncertainty
  real(sgl)          :: wchan                                       ! weight of channel
  real(sgl)          :: wnuc                                        ! weight of nucleus
  real(sgl)          :: wtot                                        ! total weight
  real(sgl)          :: ww                                          ! weight
  real(sgl)          :: xsde                                        ! uncertainty of cross section
  real(sgl)          :: xsdif                                       ! difference between experimental and theoretical cross  section
  real(sgl)          :: xse                                         ! experimental cross section
  real(sgl)          :: xslim                                       ! limit of cross section per channel
  real(sgl)          :: xst                                         ! theoretical cross section
  real(sgl)          :: x
  real(sgl)          :: pp
  real(sgl)          :: Ri                                          ! 
  real(sgl)          :: G0                                          ! 
  real(sgl)          :: term                                        ! 
  real(sgl)          :: CE                                          ! 
  real(sgl)          :: gofrel                                      ! 
  real(sgl)          :: gofR(numtalys)                              ! 
  real(dbl)          :: mersenne                                    ! random number generator (Mersenne twister)
!
! ******************************* Initialization ***********************
!
! libread    : subroutine to read library data
! expread    : subroutine to read experimental data
! expinterpol: subroutine to interpolate theoretical energy grid on experimental grid
! G          : goodness-of-fit value
!            : G = 1: chi-2
!            : G = 2: Frms
!            : G = 3: Erms
!
  Gstring(3) = '  Chi2    '
  Gstring(1) = '    Frms  '
  Gstring(2) = '    Erms  '
  isearch = isearch + 1
  do i = 1, N
    par(i) = P(i)
  enddo
  do k = 0, numenexp
    do j = 1, numsets
      do i = 1, numchanxs
        Gp(i, j, k) = 0.
        if (isearch == 0) chinclude(i, j, k) = 'y'
      enddo
    enddo
  enddo
  do i = 1, numchanxs
    Gchannel(i) = 0.
  enddo
  if (isearch == 0) then
    do m = 0, numtype
      Gmin(m) = 1.e38
    enddo
    do i = 1, numchanxs
      Gchanopt(i) = 1.e38
      Gchanopttot(i) = 1.e38
      Nopt(i) = 0
      Nopttot(i) = 0
    enddo
    open (unit=22,file='gofno',status='unknown')
    open (unit=23,file='gofnew',status='unknown')
    write(22, '("#  Run    GOF  ",a2,i4)') nuc(Ztarget), Atarget
    write(23, '("#  Run  New optimum  ",a2,i4)') nuc(Ztarget), Atarget
  endif
  if (flagcross) call expinterpol
  G = 0.
  do m = 1, numtype
    Gtot(m) = 0.
    Gtotav(m) = 0.
  enddo
  nentot = 0
  wtot = 0.
  Nsettot = 0
  Nnuctot = 0
  Nchantot = 0
  deltaE = 1.
  goffile = 'gof.0000'
  write(goffile(5:8), '(i4.4)') isearch
  open (unit = 7, file = goffile, status = 'unknown')
  write(*,'()')
  do ifile = 6, 7
    write(ifile,'("#                         Goodness-of-fit")')
    write(ifile,'("# TALYS run        :",i6)') isearch
    write(ifile,'("# GOF estimator    :",a)') Gstring(ichi2)
  enddo
!
! Write Talys input, run Talys, and read cross sections
!
! inputwrite: subroutine to write parameters for TALYS input file
!
! Loop over nuclides
!
  do inuc = 1, Nnuc
    if (mode <= 3) then
      iz = Ztarget
      ia = Atarget
    else
      iz = Z(inuc)
      ia = A(inuc)
      ext = '      '
      write(ext(1:6), '(2i3.3)') iz, ia
    endif
    if (mode >= 3) then
      call inputwrite
!
! Avoid new OMP calculation for every run if not required
!
      if (mode == 4 .and. (( .not. first .and. .not. flagompvar) .or. .not. flagecis)) then
        cmd = 'rename incident'//ext//' incident incident??????.*'
        i = system(cmd)
      endif
!
! Run Talys and read TALYS output
!
! talysread: subroutine to read results from TALYS files
! expcopy  : subroutine to remove and copy experimental data files for multi-nuclide case
!
      cmd = trim(talys)//' < talys.inp > talys.out'
      i = system(cmd)
      if (mode == 4) call expcopy
      call talysread
      if ( .not. flagtalys) then
        isearch = isearch - 1
        return
      endif
      if (flagcross) call expinterpol
      if (mode == 4 .and. .not. flagompvar) then
        cmd = 'rename incident incident'//ext//' incident.*'
        i = system(cmd)
        cmd = 'rename .tot .tot'//ext//' *.tot'
        i = system(cmd)
      endif
    endif
    if (inuc == Nnuc) first = .false.
    write(*,'()')
    do ifile = 6, 7
      write(ifile,'("# projectile       : ",a1)') ptype0
      write(ifile,'("# Z                :",i4)') iz
      write(ifile,'("# A                :",i4)') ia
      write(ifile,'("# Target isomer    : ",a1)') isochar
      write(ifile,'("# Nuclide          : ",a,i3.3,a1)') trim(nuc(iz)),ia,trim(isochar)
    enddo
!
! Calculate goodness-of-fit estimator
!
    nennucorg = 0
    nennuc = 0
    wnuc = 0.
    Nsetnuc = 0
    Nchannuc = 0
    Gnuc = 0.
    Gchanav = 0.
!
! Loop over reaction channels
!
Loop1: do i = 1, Nchanall
      if (Nsets(i) <= 0) cycle
      imt = MT(i)
      xslim = fracmax(imt) * xsmax(i)
      if (flagcross) then
        if (imt == 102) then
          xslim = max(xslim, 0.01)
        else
          xslim = max(xslim, 0.1)
        endif
      endif
Loop2: do
        do k = 1, NMTexp
          if (MTexp(k) == imt) exit Loop2
        enddo
        do k = 1, NMTlib
          if (MTlibk(k) == imt .and. isolib(k) == MTiso(i)) exit Loop2
        enddo
        if (flagtal) exit Loop2
        cycle Loop1
      enddo Loop2
      wchan = 0.
      Gchan = 0.
      nenchan(i) = 0
      nenchanorg(i) = 0
      Nsetchan(i) = 0
!
! Loop over incident energies
!
      nenset = 0
      Gset = 0.
      Gsetav = 0.
      Fset = 0.
      do j = 1, Nsets(i)
        if (isoexp(i, j) /= MTiso(i)) cycle
        if (Nenexp(i, j) <= 0) cycle
        if (outsearch >= 4) then
          if (j == 1) then
            write(7,'(/"# Channel          : ",a15)') xsfile(i)
            write(7,'("# Frms limit       :",es10.3)') Fmax(imt)
            write(7,'("# Chi2 limit       :",es10.3)') chi2max(imt)
            write(7,'("# Fracmax          :",es10.3)') fracmax(imt)
            write(7,'("# xs (mb) limit    :",es10.3)') xslim
          endif
          write(7,'(/"# Data set         :",i4,3x,a15)') j,xsfile(i)
          write(7,'("# Author           : ",a20)') auth(i,j)
          write(7,'("# Year             : ",i4)') expyear(i,j)
          write(7,'("# EXFOR subentry   : ",a9)') subentry(i,j)
          write(7,'("# Weight           : ",f12.5)') Eweight(i,j)
          write(7,'("# Quality          : ",a2)') quality(i,j)
        endif
        if (outsearch == 5) write(7,'("# energy     TALYS      exp       dexp  rel. unc. %   Chi2      C/E   include")')
        nenchanorg(i) = nenchanorg(i) + Nenexp(i, j)
        do k = 1, Nenexp(i, j)
          xst = xsth(i, j, k)
          if (isearch == 0) then
            if (xst < xslim) chinclude(i,j,k)=' '
            if (Eweight(i,j) == 0.) chinclude(i,j,k)=' '
          endif
          xsde = dxsexp(i, j, k)
          if (flagsample) then
            xse = xsexp(i, j, k) + (2. * real(mersenne()) - 1.) * xsde
            xsexpsamp(i, j, k) = xse
          else
            xse = xsexp(i, j, k)
          endif
          xse = max(xse, 0.)
          if (flagdeltaE) deltaE = dEexp(i, j, k)
          xsdif = abs(xst - xse)
          if (xse > 0.) then
            relerr = 100. * xsde / xse
          else
            relerr = 0.
          endif
!
! 1. F  (Frms)
! 2. Epsilon  (C/E) factor (asymmetry) per point)
!
! Cumulative probability derived as follows:
!
!          cdf=0.5*(1.+erf(x))
!          pp=2.*(cdf-0.5)
!
          if (xst > 0. .and. xse > 0.) then
            CE = xst/xse
            if (flagdexp .and. xsde > 0.) then
              x = (xst - xse)/(xsde*sqrt(2.))
              if (xst > xse) then
                pp = erf(x)
              else
                pp = -erf(x)
              endif
              Ri = 1. + (CE-1.)*pp
            else
              Ri = CE
            endif
            F0 = log(Ri)**2
            expo = min(F0, 80.**2)
            Gpoint(1) = exp(sqrt(expo))
            G0 = log(Ri)
            expo = min(G0, 80.)
            Gpoint(2) = exp(expo)
          else
            F0 = 80.**2
            Gpoint(1) = exp(sqrt(F0))
            G0 = 80.
            Gpoint(2) = exp(G0)
          endif
          if (isearch == 0 .and. Gpoint(1) > Fmax(imt)) chinclude(i, j, k) = ' '
!
! 3. Chi2
!
          if (xsde > 0.) then
            Gpoint(3) = (xsdif / xsde) **2
          else
            Gpoint(3) = 0.
          endif
          if (isearch == 0 .and. Gpoint(3) > chi2max(imt)) chinclude(i, j, k) = ' '
!
! Write goodness-of-fit results
!
          if (outsearch == 5) write(7,'(7es10.3,4x,a1)') Eexp(i,j,k),xst,xse,xsde,relerr, &
 &          (Gpoint(m),m=3,2,-1),chinclude(i,j,k)
!
! Summation. Outliers and too small cross sections are excluded
!
          if (chinclude(i, j, k) == 'y') then
            nenset(j) = nenset(j) + 1
            do m = 1, numtype
              Gset(j, m) = Gset(j, m) + Gpoint(m) * deltaE
            enddo
            Fset(j, 1) = Fset(j, 1) + F0 * deltaE
            Fset(j, 2) = Fset(j, 2) + G0 * deltaE
          else
            Gpoint(1) = 1.
            Gpoint(2) = 1.
            Gpoint(3) = 0.
          endif
          Gp(i, j, k) = Gpoint(ichi2)
        enddo
        if (nenset(j) == 0) cycle
        ww = Eweight(i, j)
!
! Average over incident energies
!
        if (flagdeltaE) then
          term = max(abs(Eend(i,j)-Ebeg(i,j)), 1.)
        else
          term = real(nenset(j))
        endif
        do m = 1, numtype
          if (m == 1) then
            expo = min(Fset(j,m)/term, 80.**2)
            Gsetav(j,m) = exp(sqrt(expo))
          endif
          if (m == 2) then
            expo = min(Fset(j,m)/term, 80.)
            Gsetav(j,m) = exp(expo)
          endif
          if (m == 3) Gsetav(j,m) = Gset(j,m)/term
          Gchan(m) = Gchan(m) + ww*Gsetav(j,m)
        enddo
        if (flagdeltaE) then
          nenchan(i) = nenchan(i) + 1
        else
          nenchan(i) = nenchan(i) + nenset(j)
        endif
        wchan = wchan + ww
        Nsetchan(i) = Nsetchan(i) + 1
        if (outsearch >= 4) then
          write(7,'("# Frms             :",es10.3)') Gsetav(j,1)
          write(7,'("# Erms             :",es10.3)') Gsetav(j,2)
          write(7,'("# Chi2             :",es10.3)') Gsetav(j,3)
          write(7,'("# Total points     :",i6)') Nenexp(i,j)
          write(7,'("# Included points  :",i6)') nenset(j)
          write(7,'("# End Data set     :",i4)') j
        endif
      enddo
      if (nenchan(i) == 0) cycle
      nennucorg = nennucorg + nenchanorg(i)
      nennuc = nennuc + nenchan(i)
      wnuc = wnuc + wchan
      do m = 1, numtype
        Gnuc(m) = Gnuc(m) + Gchan(m)
        Gchanav(i, m) = Gchan(m) / wchan
      enddo
      Nsetnuc = Nsetnuc + Nsetchan(i)
      Nchannuc = Nchannuc + 1
      if (outsearch >= 3) then
        write(7,'(/,"# Data set summary : ",a15)') xsfile(i)
        write(7,'("Set  Author              Year  Subentry Points Incl. Energy range      Weight Frms      Erms     Chi2")')
        do j = 1, Nsets(i)
          write(7, '(i4,1x,a20,1x,i4,1x,a9,2i5,f8.3," -",2f8.3,3es10.3)') j, auth(i,j), expyear(i,j), subentry(i,j), Nenexp(i,j), &
 &          nenset(j), Eexp(i,j,1), Eexp(i,j, Nenexp(i,j)), Eweight(i,j), (Gsetav(j,k), k=1,3)
        enddo
        write(7,'("# Channel Frms     : ",es10.3)') Gchanav(i, 1)
        write(7,'("# Channel Erms     : ",es10.3)') Gchanav(i, 2)
        write(7,'("# Channel Chi2     : ",es10.3)') Gchanav(i, 3)
        write(7,'("# End Channel      : ",a15)') xsfile(i)
      endif
      Gchannel(i) = Gchanav(i, ichi2)
      Gchannelsave(isearch, i) = Gchannel(i)
    enddo Loop1
!
! Average over reaction channels
!
    if (nennuc == 0) cycle
    nentot = nentot + nennuc
    wtot = wtot + wnuc
    do m = 1, numtype
      Gtot(m) = Gtot(m) + Gnuc(m)
      Gnucav(m) = Gnuc(m) / wnuc
    enddo
    Nsettot = Nsettot + Nsetnuc
    Nnuctot = Nnuctot + 1
    Nchantot = Nchantot + Nchannuc
    if (outsearch >= 3) then
      do ifile = 6, 7
        write(ifile, '(/,"# Channel summary  : ",a,i3.3,a1)') trim(nuc(iz)),ia,trim(isochar)
        write(ifile, '("# Channel        Sets  Points     Frms           Erms          Chi2           ", &
 &        "Total optimum         Channel optimum")')
        do i = 1,Nchanall
          if (Nsets(i) <= 0) cycle
          write(ifile, '(a15,2i6,1p,3g15.4,2(g15.4," (",i4,")"))') xsfile(i), Nsetchan(i), nenchan(i), (Gchanav(i, k), k=1,3), &
 &          Gchanopttot(i), Nopttot(i), Gchanopt(i), Nopt(i)
        enddo
      enddo
    endif
    if (outsearch >= 2) then
      write(7,'("# Nuclide channels :",i6)') Nchannuc
      write(7,'("# Nuclide data sets:",i6)') Nsetnuc
      write(7,'("# Nuclide points   :",i6)') nennucorg
      write(7,'("# Nuclide incl. pts:",i6)') nennuc
      write(7,'("# Nuclide Frms     : ",es10.3)') Gnucav(1)
      write(7,'("# Nuclide Erms     : ",es10.3)') Gnucav(2)
      write(7,'("# Nuclide Chi2     : ",es10.3)') Gnucav(3)
      write(7,'("# End Nuclide      : ",a,i3.3,a1)') trim(nuc(iz)), ia, trim(isochar)
    endif
  enddo
!
! Average over nuclides
!
  if (wtot > 0.) then
    do m = 1, numtype
      Gtotav(m) = Gtot(m) / wtot
    enddo
  endif
  G = Gtotav(ichi2)
  if (G <= 1.e-30) G = 1.e35
  Gall(0) = Gtot(ichi2)
  if (outsearch >= 1) then
    do ifile = 6, 7
      write(ifile,'()')
      write(ifile,'("# GOF no.           ",i6,1p,g12.4)') isearch,G
      write(ifile,'("# Nuclides         : ",i4)') Nnuctot
      write(ifile,'("# Channels         : ",i4)') Nchantot
      write(ifile,'("# Data sets        : ",i4)') Nsettot
      write(ifile,'("# Data points      : ",i6)') nentot
      write(ifile,'("# Frms             :",1p,g12.4," Optimum: ",g12.4)') Gtotav(1),Gmin(1)
      write(ifile,'("# Erms             :",1p,g12.4," Optimum: ",g12.4)') Gtotav(2),Gmin(2)
      write(ifile,'("# Chi2             :",1p,g12.4," Optimum: ",g12.4)') Gtotav(3),Gmin(3)
      write(ifile,'()')
    enddo
    write(22,'(i6,1p,g14.6)') isearch,G
  endif
  if (G < Gmin(0)) then
    Gmin(0) = G
    do m = 1, numtype
      Gmin(m) = Gtotav(m)
    enddo
    do i = 1, numchanxs
      Gchanopttot(i) = Gchannel(i)
      Nopttot(i) = isearch
    enddo
    if (isearch > 0) then
      do i = 1, Nchanxs
        do j = 1, Nen(i)
          xsopt(i, j) = xseval(1, i, j)
        enddo
      enddo
    endif
    do ifile = 6, 7
      write(ifile, '(2i4," Run:", i6, " New optimum (Frms Erms Chi-2): ",3es15.6)') Ztarget, Atarget, isearch, (Gmin(i),i=1,3)
    enddo
    write(23,'(i6,1p,g14.6)') isearch,Gmin(0)
    close(7)
    cmd = '\cp -f '//goffile//' gof.opt'
    i = system(cmd)
    call optimum
  else
    close(7)
  endif
  do i = 1, numchanxs
    if (Gchannel(i) < Gchanopt(i)) then
      Gchanopt(i) = Gchannel(i)
      Nopt(i) = isearch
    endif
  enddo
  Gsave(isearch) = G
  Grun(0) = G
  Gall(0) = Gall(ichi2)
  if (isearch == 0) then
    Grun0(0) = G
    Gall0(0) = Gall(ichi2)
  endif
!
! Energy dependent goodness-of-fit
!
  do i = 1, Nchanall
    do j = 1, Nsets(i)
      do k = 1, Nenexp(i, j)
        m = Eindex(i, j, k)
        Grun(m) = Grun(m) + Gp(i, j, k)
        Gall(m) = Gall(m) + Gp(i, j, k)
      enddo
    enddo
  enddo
  do m = 1, numenin
    if (Nindex(m) > 0) Grun(m) = Grun(m) / Nindex(m)
  enddo
  if (isearch == 0) then
    do m = 1, numenin
      Grun0(m) = Grun(m)
      Gall0(m) = Gall(m)
    enddo
  endif
!
! Record new optimum and make system call for saving the results belonging to the optimal parameter set.
!
! Gmin  : minimal function value
!
  if (mode >= 3) then
    if (isearch == Ntalys) then
      write( * , * ) "Maximum number of runs done: ", Ntalys
      close(unit=22)
      close(unit=23)
!
! Write improvement relative to the optimum
!
      open (unit=23,file='gofnew',status='unknown')
      read(23,'()')
      Nr = 0
      do
        Nr = Nr + 1
        read(23,*,iostat=istat) ir(Nr), gofR(Nr)
        if (istat == -1) exit
      enddo
      Nr = Nr - 1
      close(unit=23)
      open (unit=23,file='gofnew',status='unknown')
      write(23,'("#  Run  New optimum Rel. to final")')
      do i = 1, Nr
        gofrel = gofR(i)/gofR(Nr)
        write(23,'(i6,1p,2g14.6)') ir(i), gofR(i), gofrel
      enddo
      close(unit=23)

      call timer
      stop
    endif
    italys = italys + 1
  endif
  return
end subroutine fcn
! Copyright A.J. Koning 2021
