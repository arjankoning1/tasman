subroutine tasmaninitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Initialization
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
! Variables for cross section covariances
!   coveps          ! limit for covariance
!   Ecov            ! covariance energy grid
!   Rlimit          ! limit for covariance calculation
! Variables for parameter variation
!   Npar            ! number of parameters
! Variables for experimental data files
!   Nnuc            ! number of nuclides for which experimental data exists
! Variables for GOF function
!   isearch         ! number of trial run
! Variables for TMC
!   fislim          ! mass above which nuclide fissions
! Variables for reading TALYS output
!   flagtalys       ! flag to track successful TALYS calculation
!   iout            ! counter for output files
!   Nchanall        ! total number of channels
! Variables for writing TALYS input files
!   first           ! tracks whether a TALYS run is the first one
!   italys          ! TALYS run counter
! Variables for reading angular distributions
!   Nchanang        ! total number of channels with angular distributions
! Variables for reading gamma production cross sections
!   Nchangam        ! total number of channels with gamma production cross sections
! Variables for reading Legendre coefficients
!   Nchanleg        ! total number of channels with Legendre coefficients
! Variables for reading production cross sections
!   Nchanprod       ! number of channels with particle production cross sections
! Variables for reading residual production cross sections
!   Nchanrp         ! total number of channels with residual production cross sections
! Variables for reading spectra
!   Nchansp         ! total number of channels with emission spectra
! Variables for reading isotope production
!   NchanY          ! total number of channels with isotope production
! Variables for reading cross sections
!   MTnum           ! MT number of reaction channel
!   Nchanxs         ! total number of channels with cross sections
!   Nenendf         ! number of incident energies
!   Nenendf0        ! number of incident energies
! Variables for reading experimental data
!   lib             ! data library
!   NMTexp          ! number of MT numbers for experimental data
!   NMTlib          ! number of MT numbers of library
!
! *** Declaration of local data
!
  implicit none
  integer i
!
! ************************** Initialization ****************************
!
  flagtalys = .false.
  first = .true.
  italys = 0
  isearch = -1
  iout = 0
  NMTexp = 0
  NMTlib = 0
  Npar = 0
  Nnuc = 1
  Nchanxs = 0
  Nchanall = 0
  Nchanrp = 0
  NchanY = 0
  Nchangam = 0
  Nchanprod = 0
  Nchansp = 0
  Nchanang = 0
  Nchanleg = 0
  fislim = 215
  Nenendf = 0
  Nenendf0 = 0
  Ein = 0.
  Nindex = 0
  anginclude = .false.
  Ecovindex = 0
  MT = 0
  Nsets = 0
  flagchmagnet = .false.
  MTiso = -1
  xsmax = 0.
  Nen = 0
  Nen0 = 0
  Nenlow = 0
  E1mb = 0.
  xsfile = ''
  Nenexp = 0
  xsth = 0.
  Eindex = 0
  Eexp = 0.
  Ebeg = 0.
  Eend = 0.
  Eweight = 1.
  auth = ''
  quality = ''
  expyear = 0
  isoexp = -1
  subentry = ''
  xsexp = 0.
  xsexpsamp = 0.
  xsth = 0.
  Eindex = 0
  dxsexp = 0.
  Eexp = 0.
  dEexp = 0.
  Esamp = 0.
  xssamp = 0.
  Nxssamp = 0
  Nenrp = 0
  Nenrp0 = 0
  rpfile = ''
  Ntime = 0
  Yfile = ''
  Nengam = 0
  Nengam0 = 0
  gamfile = ''
  Nenprod = 0
  Nenprod0 = 0
  prodfile = ''
  Sindex = 0
  Srpindex = 0
  Sprodindex = 0
  Sgamindex = 0
  Espec = 0.
  Nensp = 0
  spfile = ''
  Nang = 0
  Eang = 0.
  Eleg = 0.
  angfile = ''
  legfile = ''
  parlow = 0
  parhigh = 0
  parsave = 0.
  pardelta = 0.
  par = 0.
  readpar = 0
  parD = 0.
  barrier = 0
  lval = 0
  Zpar = 0
  Apar = 0
  parclass = 0
  parmt = 0
  parinp = 0.
  parsym = ''
  pargam = ''
  parkey = ''
  partype = 'norm  '
  string = ''
  parinpname = ''
  parname = ''
  parstring = ''
  parfile = ''
  keyvary = ''
  sigvary = 0.
  Zvary = 0
  Avary = 0
  Zskip = 0
  Askip = 0
  parsave = 0.
  paradjust = 0.
  Pbot = 0.
  Pbin = 0.
  nhist = 0
  hist = 0.
  nhistsam = 0
  histsam = 0.
  Gsave = 1.e30
  weightsave = 1.e30
  Gchannelsave = 0.
  Z = 0
  A = 0
  Nchanmul = 0
  expfilemul = ''
  Grun0 = 0.
  Grun = 0.
  Gall0 = 0.
  Gall = 0.
  weight = 0.
  E = 0.
  dE = 0.
  xsopt = 0.
  xstalys = 0.
  xstalys2 = 0.
  xstalys3 = 0.
  xseval = 0.
  Erp = 0.
  rptalys = 0.
  Egam = 0.
  gamtalys = 0.
  Eprod = 0.
  prodtalys = 0.
  Ytalys = 0.
  weight = 1.
  Eendf = 0.
  xsendf = 0.
  timeY = 0.
  acttalys = 0.
  Nisotalys = 0.
  yieldtalys = 0.
  Nisoreltalys = 0.
  xssave = 0.
  rpsave = 0.
  gamsave = 0.
  prodsave = 0.
  angle = 0.
  angtalys = 0.
  Eout = 0.
  sptalys = 0.
  partalys = 0.
  legtalys = 0.
  leg0 = 0.
  libinclude = -1
  lib(1) = 'endfb8.0'
  lib(2) = 'jeff3.3'
  lib(3) = 'jendl5.0'
  lib(4) = 'cendl3.2'
  lib(5) = 'eaf.2010'
  lib(6) = 'irdff2.0'
  lib(7) = 'tendl.2023'
  lib(8) = 'iaea.2019'
  lib(9) = 'iaea.pd'
  lib(10) = 'ibandl'
!
! Relate MT numbers to reaction channels.
!
  MTnum = -1
  MTnum(4) = 100000
  MTnum(11) = 201000
  MTnum(16) = 200000
  MTnum(17) = 300000
  MTnum(18) = -99
  MTnum(22) = 100001
  MTnum(23) = 100003
  MTnum(24) = 200001
  MTnum(25) = 300001
  MTnum(28) = 110000
  MTnum(29) = 100002
  MTnum(30) = 200002
  MTnum(32) = 101000
  MTnum(33) = 100100
  MTnum(34) = 100010
  MTnum(35) = 101002
  MTnum(36) = 100102
  MTnum(37) = 400000
  MTnum(41) = 210000
  MTnum(42) = 310000
  MTnum(44) = 120000
  MTnum(45) = 110001
  MTnum(91) = 100000
  MTnum(102) = 000000
  MTnum(103) = 010000
  MTnum(104) = 001000
  MTnum(105) = 000100
  MTnum(106) = 000010
  MTnum(107) = 000001
  MTnum(108) = 000002
  MTnum(109) = 000003
  MTnum(111) = 020000
  MTnum(112) = 010001
  MTnum(113) = 000102
  MTnum(114) = 001002
  MTnum(115) = 011000
  MTnum(116) = 010100
  MTnum(117) = 001001
  MTnum(649) = 010000
  MTnum(699) = 001000
  MTnum(749) = 000100
  MTnum(799) = 000010
  MTnum(849) = 000001
!
! Special high-energy MT numbers for EAF-format.
!
  MTnum(152) = 500000
  MTnum(153) = 600000
  MTnum(154) = 200100
  MTnum(155) = 000101
  MTnum(156) = 410000
  MTnum(157) = 301000
  MTnum(158) = 101001
  MTnum(159) = 210001
  MTnum(160) = 700000
  MTnum(161) = 800000
  MTnum(162) = 510000
  MTnum(163) = 610000
  MTnum(164) = 710000
  MTnum(165) = 400001
  MTnum(166) = 500001
  MTnum(167) = 600001
  MTnum(168) = 700001
  MTnum(169) = 401000
  MTnum(170) = 501000
  MTnum(171) = 601000
  MTnum(172) = 300100
  MTnum(173) = 400100
  MTnum(174) = 500100
  MTnum(175) = 600100
  MTnum(176) = 200010
  MTnum(177) = 300010
  MTnum(178) = 400010
  MTnum(179) = 320000
  MTnum(180) = 300002
  MTnum(181) = 310001
  MTnum(182) = 001100
  MTnum(183) = 111000
  MTnum(184) = 110100
  MTnum(185) = 101100
  MTnum(186) = 110010
  MTnum(187) = 101010
  MTnum(188) = 100110
  MTnum(189) = 100101
  MTnum(190) = 220000
  MTnum(191) = 010010
  MTnum(192) = 001010
  MTnum(193) = 000011
  MTnum(194) = 420000
  MTnum(195) = 400002
  MTnum(196) = 410001
  MTnum(197) = 030000
  MTnum(198) = 130000
  MTnum(199) = 320001
  MTnum(200) = 520000
!
! Reaction strings
!
  MTreac=''
  MTreac(1)='(n,tot)'
  MTreac(2)='(n,el)'
  MTreac(3)='(n,non)'
  MTreac(4)="(n,n')"
  MTreac(5)="(n,x)"
  MTreac(11)='(n,2nd)'
  MTreac(16)='(n,2n)'
  MTreac(17)='(n,3n)'
  MTreac(18)='(n,f)'
  MTreac(22)='(n,na)'
  MTreac(23)='(n,n3a)'
  MTreac(24)='(n,2na)'
  MTreac(25)='(n,3na)'
  MTreac(28)='(n,np)'
  MTreac(29)='(n,n2a)'
  MTreac(30)='(n,2n2a)'
  MTreac(32)='(n,nd)'
  MTreac(33)='(n,nt)'
  MTreac(34)='(n,nh)'
  MTreac(35)='(n,nd2a)'
  MTreac(36)='(n,nt2a)'
  MTreac(37)='(n,4n)'
  MTreac(41)='(n,2np)'
  MTreac(42)='(n,3np)'
  MTreac(44)='(n,n2p)'
  MTreac(45)='(n,npa)'
  do i=0,40
    MTreac(50+i)="(n,n'_01) "
    write(MTreac(50+i)(7:8),'(i2.2)') i
  enddo
  MTreac(91)="(n,n'_con)"
  MTreac(101)='(n,abs)'
  MTreac(102)='(n,g)'
  MTreac(103)='(n,p)'
  MTreac(104)='(n,d)'
  MTreac(105)='(n,t)'
  MTreac(106)='(n,h)'
  MTreac(107)='(n,a)'
  MTreac(108)='(n,2a)'
  MTreac(109)='(n,3a)'
  MTreac(111)='(n,2p)'
  MTreac(112)='(n,pa)'
  MTreac(113)='(n,t2a)'
  MTreac(114)='(n,d2a)'
  MTreac(115)='(n,pd)'
  MTreac(116)='(n,pt)'
  MTreac(117)='(n,da)'
  MTreac(201)='(n,xn)'
  MTreac(202)='(n,xg)'
  MTreac(203)='(n,xp)'
  MTreac(204)='(n,xd)'
  MTreac(205)='(n,xt)'
  MTreac(206)='(n,xh)'
  MTreac(207)='(n,xa)'
  do i=0,40
    MTreac(600+i)="(n,p_00)"
    write(MTreac(600+i)(6:7),'(i2.2)') i
    MTreac(650+i)="(n,d_00)"
    write(MTreac(650+i)(6:7),'(i2.2)') i
    MTreac(700+i)="(n,t_00)"
    write(MTreac(700+i)(6:7),'(i2.2)') i
    MTreac(750+i)="(n,h_00)"
    write(MTreac(750+i)(6:7),'(i2.2)') i
    MTreac(800+i)="(n,a_00)"
    write(MTreac(800+i)(6:7),'(i2.2)') i
  enddo
  MTreac(649)='(n,p_con)'
  MTreac(699)='(n,d_con)'
  MTreac(749)='(n,t_con)'
  MTreac(799)='(n,h_con)'
  MTreac(849)='(n,a_con)'
!
! Create energy grid for covariance data
!
  Ecov(0) = 0.
  Ecov(1) = 1.e-11
  Ecov(2) = 2.53e-8
  Ecov(3) = 1.e-6
  Ecov(4) = 1.e-5
  Ecov(5) = 1.e-4
  Ecov(6) = 0.001
  Ecov(7) = 0.002
  Ecov(8) = 0.005
  Ecov(9) = 0.01
  Ecov(10) = 0.02
  Ecov(11) = 0.05
  Ecov(12) = 0.1
  Ecov(13) = 0.2
  Ecov(14) = 0.4
  Ecov(15) = 0.6
  Ecov(16) = 0.8
  Ecov(17) = 1.
  Ecov(18) = 1.4
  Ecov(19) = 1.8
  Ecov(20) = 2.2
  Ecov(21) = 2.6
  Ecov(22) = 3.
  Ecov(23) = 4.
  Ecov(24) = 5.
  Ecov(25) = 6.
  Ecov(26) = 7.
  Ecov(27) = 8.
  Ecov(28) = 10.
  Ecov(29) = 12.
  Ecov(30) = 14.
  Ecov(31) = 16.
  Ecov(32) = 18.
  Ecov(33) = 20.
  Ecov(34) = 22.
  Ecov(35) = 24.
  Ecov(36) = 26.
  Ecov(37) = 28.
  Ecov(38) = 30.
  Ecov(39) = 34.
  Ecov(40) = 38.
  Ecov(41) = 42.
  Ecov(42) = 46.
  Ecov(43) = 50.
  Ecov(44) = 55.
  Ecov(45) = 60.
  Ecov(46) = 70.
  Ecov(47) = 80.
  Ecov(48) = 100.
  Ecov(49) = 120.
  Ecov(50) = 140.
  Ecov(51) = 160.
  Ecov(52) = 200.
  coveps = 1.e-15
  Rlimit = 1.000001
  return
end subroutine tasmaninitial
! Copyright A.J. Koning 2021
