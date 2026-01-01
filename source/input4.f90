subroutine input4
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read input for fourth set of variables
!
! Author    : Arjan Koning
!
! 2025-12-26: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl           ! single precision kind
! All global variables
!   nummass       ! maximum number of masses
!   numpar        ! maximum number of parameters
! Variables for reading TASMAN input lines
!   inline        ! input line
!   nlines        ! number of input lines
! Variables for automatic parameter variation
!   flagparvary     ! flag to automatically vary parameters
!   ACN             ! A of compound nucleus
!   ZCN             ! Z of compound nucleus
! Variables for TMC
!   fislim        ! mass above which nuclide fissions
! Variables for parameter variation
!   Npar          ! number of parameters
!   par           ! parameter value
!   paradjust     ! local adjustment parameters
!   pardelta      ! uncertainty of parameter
!   parinp        ! input parameter value
!   parinpname    ! input character parameter value
!   parkey        ! TALYS keyword
!   parname       ! character parameter value
!   parsave       ! parameter value
!   partalys      ! parameter value
!   partype       ! type of parameter variation
! Variables for parameter covariances
!   parstring     ! parameter string
! Variables for fitting limits
!   parclass      ! class of keyword
! Variables for writing TALYS input files
!   Apar          ! A - index for keyword
!   barrier       ! fission barrier
!   flagompvar    ! flag for variation of OMP parameters
!   lval          ! l value
!   mode          ! TASMAN mode
!   pargam        ! gamma - ray related input string
!   parmt         ! MT number
!   parsym        ! particle symbol
!   Zpar          ! Z - index for keyword
! Variables for uncertainty
!   cwidth        ! scaling constant for parameter variation
!   fiscor0       ! adjustment parameter for fission uncertainties
! Variables for processing input
!   Atarget       ! mass number of target nucleus
!   changeline    ! line number at which TASMAN specific input starts
!   ptype0        ! type of incident particle
!   Ztarget       ! charge number of target nucleus
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: ch              ! character
  character(len=1)   :: ptype           ! particle type
  character(len=2)   :: gamtype         ! gamma type
  character(len=6)   :: partype0        ! particle type
  character(len=16)  :: cval            ! character value
  character(len=132) :: key             ! keyword
  character(len=132) :: word(40)        ! words on input line
  character(len=132) :: line            ! input line
  integer            :: class           ! input class
  integer            :: i               ! counter
  integer            :: ia              ! mass number of nucleus
  integer            :: ibar            ! fission barrier
  integer            :: ii              ! counter
  integer            :: imt             ! MT counter
  integer            :: istat           ! logical for file access
  integer            :: ival            ! integer value
  integer            :: iz              ! charge number of nucleus
  integer            :: j               ! counter
  integer            :: k               ! Legendre order
  integer            :: l               ! counter
  integer            :: lenpar          ! length of parameter string
  real(sgl)          :: D               ! distance from constant uncertainty area
  real(sgl)          :: delta           ! relative deviation
  real(sgl)          :: delta0          ! systematical pairing energy
  real(sgl)          :: dripvar         ! function to set parameter uncertainties for nuclides  closer to dripline
  real(sgl)          :: Ea              ! parameters for local energy dependence
  real(sgl)          :: Eb              ! end energy of local adjustment
  real(sgl)          :: Em              ! intermediate energy of local adjustment
  real(sgl)          :: fiscor(0:nummass) ! adjustment parameter for fission uncertainties
  real(sgl)          :: Gunc            ! multiplier for uncertainty
  real(sgl)          :: val             ! value
!
! **************************** Initialization **************************
!
! Set default values
!
  flagompvar = .false.
!
! For actinides, the relative uncertainty of 2nd, 3rd, etc fission is increased to get more realistic results.
! For this we use the fiscor parameter, which we multiply with some fission related parameters.
! Note that for non-fissile nuclides, fiscor=1.
!
  if (Atarget <= fislim .or. ptype0 /= 'n') then
    do k = 0, min(Atarget + 4, nummass)
      fiscor(k) = 1.
    enddo
  else
    do k = 0, min(Atarget + 4, nummass)
      fiscor(k) = fiscor0
    enddo
  endif
  fiscor(Atarget + 1) = 1.
!
! **************** Compound nucleus definition *************************
!
      if (ptype0 == 'g') then
        ZCN = Ztarget
        ACN = Atarget
      endif
      if (ptype0 == 'n') then
        ZCN = Ztarget
        ACN = Atarget + 1
      endif
      if (ptype0 == 'p') then
        ZCN = Ztarget + 1
        ACN = Atarget + 1
      endif
      if (ptype0 == 'd') then
        ZCN = Ztarget + 1
        ACN = Atarget + 2
      endif
      if (ptype0 == 't') then
        ZCN = Ztarget + 1
        ACN = Atarget + 3
      endif
      if (ptype0 == 'h') then
        ZCN = Ztarget + 2
        ACN = Atarget + 3
      endif
      if (ptype0 == 'a') then
        ZCN = Ztarget + 2
        ACN = Atarget + 4
      endif
!
! ********************** Read input parameters *************************
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
! Special option to read in and vary ALL TALYS parameters or a subset without having to give them explicitly in the
! TASMAN input file
!
  if (flagparvary) call parameters
  ii = 0
  Loop1: do i = changeline + 1, nlines
    line = inline(i)
    key = ''
    call getkeywords(line, word)
    key = word(1)
    ch = word(2)(1:1)
    if (ch == '#') cycle Loop1
!
! TALYS input parameters
!
    delta = 0.
    val = 0.
!
! For each parameter, we read the value using exactly the same input rules as for TALYS.
! However, at the end of each input line, now the variation range can be given.
! For this, we use a fraction of 1, i.e. delta=0.40 means 40%.
! If no variation is given in the input, the default variations are set below for each parameter separately.
! For nuclides around the stability line, dripvar=1.
!
! See the TALYS manual for a description of the input parameters
!
    iz = -1
    ia = -1
    ptype = ''
    ibar = 0
    l = -1
    imt = -1
    gamtype = ''
    cval = ''
    partype0 = 'norm  '
    Ea = 0.
    Eb = 0.
    Em = 0.
    D = 0.
!
! Different classes of keywords
!
! Class=0: keyword value(integer) delta
!
    do
      class = 0
      delta0 = 0.
      if (key == 'ldmodel') then
        delta0 = 6.
        partype0 = 'model '
        exit
      endif
      if (key == 'fismodel') then
        delta0 = 5.
        partype0 = 'model '
        exit
      endif
      if (key == 'widthmode') then
        delta0 = 4.
        partype0 = 'model '
        exit
      endif
      if (key == 'deuteronomp') then
        delta0 = 5.
        partype0 = 'model '
        exit
      endif
      if (key == 'alphaomp') then
        delta0 = 5.
        partype0 = 'model '
        exit
      endif
      if (key == 'kvibmodel') then
        delta0 = 2.
        partype0 = 'model '
        exit
      endif
      if (key == 'jlmmode') then
        delta0 = 4.
        partype0 = 'model '
        exit
      endif
      if (key == 'radialmodel') then
        delta0 = 2.
        partype0 = 'model '
        exit
      endif
      if (key == 'axtype') then
        delta0 = 5.
        partype0 = 'model '
        exit
      endif
      if (key == 'massmodel') then
        delta0 = 4.
        partype0 = 'model '
        exit
      endif
      if (key == 'shellmodel') then
        delta0 = 2.
        partype0 = 'model '
        exit
      endif
      if (key == 'spincutmodel') then
        delta0 = 2.
        partype0 = 'model '
        exit
      endif
      if (key == 'preeqmode') then
        delta0 = 4.
        partype0 = 'model '
        exit
      endif
      if (key == 'mpreeqmode') then
        delta0 = 2.
        partype0 = 'model '
        exit
      endif
      if (key == 'pairmodel') then
        delta0 = 2.
        partype0 = 'model '
        exit
      endif
      if (key == 'strength') then
        delta0 = 9.
        partype0 = 'model '
        exit
      endif
      if (key == 'strengthm1') then
        delta0 = 4.
        partype0 = 'model '
        exit
      endif
!
! Class=1: keyword value(real) delta
!
      class = 1
      if (key == 'rspincut') then
        delta0 = 0.30
        partype0 = 'factor'
        exit
      endif
      if (key == 'rpipi' .or. key == 'rpinu' .or. key == 'rnupi' .or. key == 'rnunu') then
        delta0 = 0.30
        partype0 = 'factor'
        exit
      endif
      if (key == 'alphald') then
        delta0 = 0.30
        partype0 = 'norm  '
        exit
      endif
      if (key == 'betald') then
        delta0 = 0.30
        partype0 = 'norm  '
        exit
      endif
      if (key == 'gammashell1' .or. key == 'gammashell2') then
        delta0 = 0.30
        partype0 = 'norm  '
        exit
      endif
      if (key == 'kph') then
        delta0 = 0.1125 - 3.125e-4 * Atarget
        partype0 = 'norm  '
        exit
      endif
      if (key == 'rgamma') then
        delta0 = 0.50
        partype0 = 'norm  '
        exit
      endif
      if (key == 'pairconstant') then
        delta0 = 0.20
        partype0 = 'norm  '
        exit
      endif
      if (key == 'm2constant') then
        delta0 = 0.20
        partype0 = 'factor'
        exit
      endif
      if (key == 'm2limit') then
        delta0 = 0.30
        partype0 = 'factor'
        exit
      endif
      if (key == 'm2shift') then
        delta0 = 0.30
        partype0 = 'factor'
        exit
      endif
      if (key == 'esurf') then
        delta0 = 0.20
        partype0 = 'norm  '
        exit
      endif
      if (key == 'aradialcor') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'adepthcor') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'astroe') then
        delta0 = 0.10
        partype0 = 'norm  '
        exit
      endif
      if (key == 'astrot') then
        delta0 = 0.10
        partype0 = 'norm  '
        exit
      endif
      if (key == 'cglobal') then
        delta0 = 0.30
        partype0 = 'shift '
        exit
      endif
      if (key == 'pglobal') then
        delta0 = 0.30
        partype0 = 'shift '
        exit
      endif
      if (key == 'elwidth') then
        delta0 = 0.40
        partype0 = 'norm  '
        exit
      endif
      if (key == 'lvadjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'lwadjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'lv1adjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'lw1adjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'lvsoadjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'lwsoadjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'rprime') then
        delta0 = 0.05
        partype0 = 'norm  '
        exit
      endif
      if (key == 'pshiftconstant') then
        delta0 = 0.3
        partype0 = 'shift '
        exit
      endif
      if (key == 'cbarrier') then
        delta0 = 0.80
        partype0 = 'factor'
        exit
      endif
      if (key == 'gmradjustd') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'gmradjuste') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'gmradjustg') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'gqradjustd') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'gqradjuste') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'gqradjustg') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'leoradjustd') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'leoradjuste') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'leoradjustg') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'heoradjustd') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'heoradjuste') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
      if (key == 'heoradjustg') then
        delta0 = 0.40
        partype0 = 'factor'
        exit
      endif
!
! Class=2: keyword Z A value(real) delta
!
      class = 2
      if (key == 'a' .or. key == 'alimit' .or. key == 'g' .or. key == 'gn' .or. key == 'gp') then
        delta0 = 0.1125 - 3.125e-4 * Atarget
        partype0 = 'norm  '
        exit
      endif
      if (key == 'gammald') then
        delta0 = 0.30
        partype0 = 'norm  '
        exit
      endif
      if (key == 'aadjust' .or. key == 'gnadjust' .or. key == 'gpadjust' .or. key == 'gadjust') then
        delta0 = 0.1125 - 3.125e-4 * Atarget
        partype0 = 'factor'
        exit
      endif
      if (key == 'betafiscor') then
        delta0 = 0.40
        partype0 = 'norm  '
        exit
      endif
      if (key == 'betafiscoradjust') then
        delta0 = 0.20
        partype0 = 'norm  '
        exit
      endif
      if (key == 'vfiscor') then
        delta0 = 0.30
        partype0 = 'norm  '
        exit
      endif
      if (key == 'vfiscoradjust') then
        delta0 = 0.15
        partype0 = 'norm  '
        exit
      endif
      if (key == 'rmiufiscor') then
        delta0 = 0.80
        partype0 = 'norm  '
        exit
      endif
      if (key == 'rmiufiscoradjust') then
        delta0 = 0.40
        partype0 = 'norm  '
        exit
      endif
      if (key == 'gamgam') then
        delta0 = 0.40
        partype0 = 'norm  '
        exit
      endif
      if (key == 'gamgamadjust') then
        delta0 = 0.20
        partype0 = 'factor'
        exit
      endif
      if (key == 'risomer') then
        delta0 = 0.80
        partype0 = 'factor'
        exit
      endif
      if (key == 'd0') then
        delta0 = 0.30
        partype0 = 'norm  '
        exit
      endif
      if (key == 'pair') then
        delta0 = 2.
        partype0 = 'shift '
        exit
      endif
      if (key == 'massnucleus') then
        delta0 = 0.01
        partype0 = 'norm  '
        exit
      endif
      if (key == 'massexcess') then
        delta0 = 50.
        partype0 = 'shift '
        exit
      endif
!
! Class=3: keyword Z A value(string) delta
!
      class = 3
      if (key == 'hbtransfile') then
        delta0 = 0.20
        partype0 = 'norm  '
        exit
      endif
      if (key == 'class2file') then
        class = 3
        delta0 = 0.20
        partype0 = 'norm  '
        exit
      endif
!
! Class=4: keyword Z A value(real) barrier delta
!
      class = 4
      if (key == 'deltaw') then
        delta0 = 1.
        partype0 = 'shift '
        exit
      endif
      if (key == 'rtransmom') then
        delta0 = 0.80
        partype0 = 'factor'
        exit
      endif
      if (key == 'rclass2mom') then
        delta0 = 0.80
        partype0 = 'factor'
        exit
      endif
      if (key == 'class2width') then
        delta0 = 0.50
        partype0 = 'norm  '
        exit
      endif
      if (key == 'krotconstant') then
        delta0 = 0.80
        partype0 = 'factor'
        exit
      endif
      if (key == 'ufermi' .or. key == 'cfermi' .or. key == 'ufermibf' .or. key == 'cfermibf') then
        delta0 = 1.
        partype0 = 'norm  '
        exit
      endif
      if (key == 's2adjust') then
        delta0 = 0.30
        partype0 = 'factor'
        exit
      endif
      if (key == 'ctable') then
        delta0 = 1.00
        partype0 = 'shift '
        exit
      endif
      if (key == 'ptable') then
        delta0 = 2.00
        partype0 = 'shift '
        exit
      endif
      if (key == 'ctableadjust') then
        delta0 = 0.50
        partype0 = 'shift '
        exit
      endif
      if (key == 'ptableadjust') then
        delta0 = 1.00
        partype0 = 'shift '
        exit
      endif
      if (key == 'fisbar') then
        delta0 = 0.10
        partype0 = 'norm  '
        exit
      endif
      if (key == 'fisbaradjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'beta2') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'fishw') then
        delta0 = 0.10
        partype0 = 'norm  '
        exit
      endif
      if (key == 'fishwadjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'bdamp') then
        delta0 = 0.10
        partype0 = 'norm  '
        exit
      endif
      if (key == 'bdampadjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'pshift') then
        delta0 = 1.
        partype0 = 'shift '
        exit
      endif
      if (key == 'e0') then
        delta0 = 0.15
        partype0 = 'shift '
        exit
      endif
      if (key == 'e0adjust') then
        delta0 = 0.075
        partype0 = 'shift '
        exit
      endif
      if (key == 'exmatch') then
        delta0 = 0.15
        partype0 = 'shift '
        exit
      endif
      if (key == 'exmatchadjust') then
        delta0 = 0.075
        partype0 = 'shift '
        exit
      endif
      if (key == 't') then
        delta0 = 0.10
        partype0 = 'norm  '
        exit
      endif
      if (key == 'tadjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
!
! Class=5: keyword Z A value(integer) barrier delta
!
      class = 5
      if (key == 'nlow') then
        delta0 = 3.
        partype0 = 'shift '
        exit
      endif
      if (key == 'ntop') then
        delta0 = 8.
        partype0 = 'shift '
        exit
      endif
!
! Class=6: keyword Z A value radiationtype delta
!
      class = 6
      if (key == 'sgr') then
        delta0 = 0.20
        partype0 = 'norm  '
        exit
      endif
      if (key == 'sgradjust') then
        delta0 = 0.20
        partype0 = 'factor'
        exit
      endif
      if (key == 'egr') then
        delta0 = 0.10
        partype0 = 'norm  '
        exit
      endif
      if (key == 'egradjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'ggr') then
        delta0 = 0.20
        partype0 = 'norm  '
        exit
      endif
      if (key == 'ggradjust') then
        delta0 = 0.20
        partype0 = 'factor'
        exit
      endif
      if (key == 'spr') then
        delta0 = 0.20
        partype0 = 'norm  '
        exit
      endif
      if (key == 'spradjust') then
        delta0 = 0.20
        partype0 = 'factor'
        exit
      endif
      if (key == 'epr') then
        delta0 = 0.10
        partype0 = 'norm  '
        exit
      endif
      if (key == 'epradjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'gpr') then
        delta0 = 0.20
        partype0 = 'norm  '
        exit
      endif
      if (key == 'gpradjust') then
        delta0 = 0.20
        partype0 = 'factor'
        exit
      endif
      if (key == 'etable') then
        delta0 = 0.80
        partype0 = 'shift '
        exit
      endif
      if (key == 'ftable') then
        delta0 = 0.50
        partype0 = 'factor'
        exit
      endif
      if (key == 'wtable') then
        delta0 = 0.30
        partype0 = 'factor'
        exit
      endif
      if (key == 'etableadjust') then
        delta0 = 0.40
        partype0 = 'shift '
        exit
      endif
      if (key == 'ftableadjust') then
        delta0 = 0.25
        partype0 = 'factor'
        exit
      endif
      if (key == 'wtableadjust') then
        delta0 = 0.15
        partype0 = 'factor'
        exit
      endif
      if (key == 'upbendc') then
        delta0 = 1.00
        partype0 = 'factor'
        exit
      endif
      if (key == 'upbende') then
        delta0 = 1.00
        partype0 = 'factor'
        exit
      endif
      if (key == 'upbendf') then
        delta0 = 1.00
        partype0 = 'factor'
        exit
      endif
      if (key == 'upbendcadjust') then
        delta0 = 0.50
        partype0 = 'factor'
        exit
      endif
      if (key == 'upbendeadjust') then
        delta0 = 0.50
        partype0 = 'factor'
        exit
      endif
      if (key == 'upbendfadjust') then
        delta0 = 0.50
        partype0 = 'factor'
        exit
      endif
!
! Class=7: keyword symbol value delta (optical model)
!
      class = 7
      if (key == 'v1adjust') then
        delta0 = 0.02
        partype0 = 'factor'
        exit
      endif
      if (key == 'v2adjust') then
        delta0 = 0.03
        partype0 = 'factor'
        exit
      endif
      if (key == 'v3adjust') then
        delta0 = 0.03
        partype0 = 'factor'
        exit
      endif
      if (key == 'v4adjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'rvadjust') then
        delta0 = 0.02
        partype0 = 'factor'
        exit
      endif
      if (key == 'avadjust') then
        delta0 = 0.02
        partype0 = 'factor'
        exit
      endif
      if (key == 'rvdadjust') then
        delta0 = 0.03
        partype0 = 'factor'
        exit
      endif
      if (key == 'avdadjust') then
        delta0 = 0.04
        partype0 = 'factor'
        exit
      endif
      if (key == 'rwadjust') then
        delta0 = 0.02
        partype0 = 'factor'
        exit
      endif
      if (key == 'awadjust') then
        delta0 = 0.02
        partype0 = 'factor'
        exit
      endif
      if (key == 'w1adjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'w2adjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'w3adjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'w4adjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'rwdadjust') then
        delta0 = 0.03
        partype0 = 'factor'
        exit
      endif
      if (key == 'awdadjust') then
        delta0 = 0.04
        partype0 = 'factor'
        exit
      endif
      if (key == 'd1adjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'd2adjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'd3adjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'rvsoadjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'avsoadjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'vso1adjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
      if (key == 'vso2adjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'rwsoadjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'awsoadjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'wso1adjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'wso2adjust') then
        delta0 = 0.10
        partype0 = 'factor'
        exit
      endif
      if (key == 'rcadjust') then
        delta0 = 0.05
        partype0 = 'factor'
        exit
      endif
!
! Class=8: keyword symbol value delta (pre-equilibrium)
!
      class = 8
      if (key == 'cstrip') then
        delta0 = 1.00
        partype0 = 'factor'
        exit
      endif
      if (key == 'cknock') then
        delta0 = 1.00
        partype0 = 'factor'
        exit
      endif
      if (key == 'cbreak') then
        delta0 = 1.00
        partype0 = 'factor'
        exit
      endif
!
! Class=9: keyword symbol value l-value delta (transmission coefficient)
!
      class = 9
      if (key == 'tljadjust') then
        delta0 = 0.80
        partype0 = 'factor'
        exit
      endif
!
! Class=10: keyword mtnumber filename value(real) delta
!
      class = 10
      if (key == 'rescuefile') then
        delta0 = 0.10
        partype0 = 'norm  '
        exit
      endif
      cycle Loop1
    enddo
!
! Read keywords
!
    k = 1
    if (class >= 2 .and. class <= 6) then
      k = k + 1
      read(word(k), * , iostat = istat) iz
      if (istat /= 0) call read_error(line, istat)
      k = k + 1
      read(word(k), * , iostat = istat) ia
      if (istat /= 0) call read_error(line, istat)
      if (iz <= 2 .and. ia <= 4) then
        iz = ZCN - iz
        ia = ACN - ia
      endif
    endif
    if (class >= 7 .and. class <= 9) then
      k = k + 1
      ptype = ch
    endif
    if (class == 10) then
      k = k + 1
      read(word(k), * , iostat = istat) imt
      if (istat /= 0) call read_error(line, istat)
      k = k + 1
      read(word(k), * , iostat = istat) cval
      if (istat /= 0) call read_error(line, istat)
    endif
    if (class == 0) then
      k = k + 1
      read(word(k), * , iostat = istat) ival
      if (istat /= 0) call read_error(line, istat)
      val = real(ival)
    endif
    if (class == 1 .or. class == 2 .or. class == 4 .or. class >= 6) then
      k = k + 1
      read(word(k), * , iostat = istat) val
      if (istat /= 0) call read_error(line, istat)
    endif
    if (class == 3) then
      k = k + 1
      read(word(k), * , iostat = istat) cval
      if (istat /= 0) call read_error(line, istat)
    endif
    if (class == 5) then
      k = k + 1
      read(word(k), * , iostat = istat) ival
      if (istat /= 0) call read_error(line, istat)
      val = real(ival)
    endif
    do
      if (class == 4 .or. class == 5) then
        k = k + 1
        read(word(k), '(i1)', iostat = istat) ibar
        if (istat == -1) exit
        if (istat /= 0) call read_error(line, istat)
      endif
      if (class == 6) then
        k = k + 1
        read(word(k), * , iostat = istat) gamtype
        if (istat == -1) exit
        if (istat /= 0) call read_error(line, istat)
      endif
      if (class == 9) then
        k = k + 1
        read(word(k), * , iostat = istat) l
        if (istat == -1) exit
        if (istat /= 0) call read_error(line, istat)
      endif
      if (word(k+1) /= ' ' .and. word(k+2) /= ' ' .and. word(k+3) /= ' ' .and. word(k+4) /= ' ') then
        k = k + 1
        read(word(k), * , iostat = istat) Ea
        if (istat /= 0) call read_error(line, istat)
        k = k + 1
        read(word(k), * , iostat = istat) Eb
        if (istat /= 0) call read_error(line, istat)
        k = k + 1
        read(word(k), * , iostat = istat) Em
        if (istat /= 0) call read_error(line, istat)
        k = k + 1
        read(word(k), * , iostat = istat) D
        if (istat /= 0) call read_error(line, istat)
      endif
      k = k + 1
      read(word(k), * , iostat = istat) delta
      if (istat == -1) exit
      if (istat /= 0) call read_error(line, istat)
      exit
    enddo
!
! Set default parameter uncertainties
!
! dripvar: function to set parameter uncertainties for nuclides closer to dripline
!
    if (delta == 0.) then
      if (class >= 2 .and. class <= 4) then
        Gunc = dripvar(iz, ia) * fiscor(ia)
      else
        Gunc = dripvar(Ztarget, Atarget) * fiscor(Atarget)
      endif
      Gunc = min(Gunc, 3.) 
      delta = Gunc * delta0
      if (class == 7) then
        if (ptype == 'p') delta = 2.*delta
        if (ptype == 'd') delta = 3.*delta
        if (ptype == 't') delta = 3.*delta
        if (ptype == 'h') delta = 4.*delta
        if (ptype == 'a') delta = 2.*delta
      endif
      if (class == 8) then
        if (ptype == 'd') delta = 2.*delta
        if (ptype == 't') delta = 2.*delta
        if (ptype == 'h') delta = 2.*delta
      endif
      if (key == 'pair' .and. mod(iz, 2) == 1 .and. mod(ia-iz, 2) == 1) delta = 0.0001
      if (class == 0) delta = delta0
    endif
!
! Record keywords, their values and uncertainties from TASMAN input file
!
    if (class == 7) flagompvar = .true.
    ii = ii + 1
    parkey(ii) = key
    parinp(ii) = val
    parsym(ii) = ptype
    Zpar(ii) = iz
    Apar(ii) = ia
    barrier(ii) = ibar
    lval(ii) = l
    partype(ii) = partype0
    pargam(ii) = gamtype
    parinpname(ii) = cval(1:12)
    parname(ii) = cval
    parmt(ii) = imt
    paradjust(ii, 1) = Ea
    paradjust(ii, 2) = Eb
    paradjust(ii, 3) = Em
    paradjust(ii, 4) = D
    parclass(ii) = class
!
! Create parameter string
!
    parstring(ii) = parkey(ii)(1:30)
    lenpar = 16
    do j = 1, 16
      if (parkey(ii)(j:j) == ' ') then
        lenpar = j
        exit
      endif
    enddo
    if (iz /=  -1) then
      write(parstring(ii)(lenpar+1:lenpar+3), '(i3)') iz
      lenpar = lenpar + 4
    endif
    if (ia /=  -1) then
      write(parstring(ii)(lenpar+1:lenpar+3), '(i3)') ia
      lenpar = lenpar + 4
    endif
    if (ptype /= ' ') then
      write(parstring(ii)(lenpar+1:lenpar+1), '(a1)') ptype
      lenpar = lenpar + 2
    endif
    if (ibar > 0) then
      write(parstring(ii)(lenpar+1:lenpar+1), '(i1)') ibar
      lenpar = lenpar + 2
    endif
    if (l >  -1) then
      write(parstring(ii)(lenpar+1:lenpar+1), '(i1)') l
      lenpar = lenpar + 2
    endif
    if (gamtype /= '  ') then
      write(parstring(ii)(lenpar+1:lenpar+2), '(a2)') gamtype
      lenpar = lenpar + 3
    endif
!
! The parameter variation can be scaled with a global parameter cwidth.
! For linear sensitivity matrices, we use a smaller variation.
!
    if (class > 0) then
      pardelta(ii) = cwidth * delta
!     if (mode == 2) then
!       if (flagmorris) then
!         pardelta(ii) = 0.3 * pardelta(ii)
!       else
!         pardelta(ii) = 0.1 * pardelta(ii)
!       endif
!     endif
    else
      pardelta(ii) = delta
    endif
    par(ii) = parinp(ii)
    partalys(0, ii) = parinp(ii)
    parsave(0, ii) = parinp(ii)
    if (ii == numpar) exit
  enddo Loop1
  Npar = ii
  return
end subroutine input4
! Copyright A.J. Koning 2021
