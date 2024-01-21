subroutine talysread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read results from TALYS files
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
! All global variables
!   numchanxs       ! maximum number of channels with cross sections
! Variables for reading TALYS output
!   flagcross       ! flag for covariances or optimization of cross sections
!   flagexp         ! flag to optimize to experimental data
!   flaggamma       ! flag for covariances or optimization of discrete gamma ray transitions
!   flagintegral    ! flag for covariances or optimization of integral activation data
!   flagld          ! flag for optimization of level density parameters
!   flaglib         ! flag to optimize to nuclear data libraries
!   flagmacs        ! flag for covariances or optimization of Maxwellian averaged cross sections
!   flagmagnet      ! flag to push TALYS data to experimental or library data
!   flagprod        ! flag for covariances or optimization of particle production cross sections
!   flagproduct     ! flag for covariances or optimization of isotope production
!   flagpsf         ! flag for optimization of photon strength functions
!   flagresidual    ! flag for covariances or optimization of residual production cross sections
!   flagtal         ! flag to optimize to TALYS
!   flagtalys       ! flag to track successful TALYS calculation
!   iout            ! counter for output files
!   Nchanall        ! total number of channels
!   Ntalbeg         ! counter for first TALYS calculation
! Variables for writing TALYS input files
!   flagangle       ! flag for covariances or optimization of angular distributions
!   flagleg         ! flag for covariances or optimization of Legendre coefficients
!   flagspectra     ! flag for covariances or optimization of emission spectra
!   italys          ! TALYS run counter
! Variables for reading gamma production cross sections
!   Nchangam        ! total number of channels with gamma production cross sections
! Variables for reading production cross sections
!   Nchanprod       ! number of channels with particle production cross sections
! Variables for reading residual production cross sections
!   Nchanrp         ! total number of channels with residual production cross sections
! Variables for reading cross sections
!   Nchanxs         ! total number of channels with cross sections
! Variables for GOF function
!   Ntalys          ! number of TALYS runs
!
! *** Declaration of local data
!
  implicit none
  logical            :: lexist      ! logical to determine existence
  character(len=14)  :: ioutfile    ! file with TALYS input
  character(len=132) :: cmd         ! command
  integer            :: isys        ! help variable
  integer            :: system      ! system call
!
! ****************** Check successful TALYS run ************************
!
! Check whether last TALYS calculation was successful, otherwise reset the counter and write the problematic TALYS input file to
! a separate file.
!
  inquire (file = 'nonelastic.tot', exist = lexist)
  if (lexist .or. flagmacs .or. flaggamgam  .or. flagld .or. flagpsf) then
    flagtalys = .true.
  else
    write( * , * ) " TALYS calculation ", italys, " unsuccessful"
    iout = iout + 1
    flagtalys = .false.
    if (iout <= 100) then
      ioutfile = 'talys.not.0000'
      write(ioutfile(11:14), '(i4.4)') iout
      cmd = 'cp talys.inp '//ioutfile
      isys = system(cmd)
    endif
    if (iout > Ntalys) write(*, '(" TASMAN-error: too many", " unsuccessful TALYS runs:", i4)') iout
    return
  endif
!
! ************ Determine open reaction channels and read results  ******
!
! xsread      : subroutine to read cross sections from TALYS production cross sections
! resread     : subroutine to read residual production cross sections from TALYS
! prodread    : subroutine to read particle production cross sections from TALYS
! gamread     : subroutine to read gamma production cross sections from TALYS
! specread    : subroutine to read emission spectra from TALYS distributions
! angread     : subroutine to read angular distributions from TALYS coefficients
! legread     : subroutine to read Legendre coefficients from TALYS activation data
! integralread: subroutine to read integral activation data from TALYS files
! productread : subroutine to read isotope production from TALYS averaged cross sections
! macsread    : subroutine to read Maxwellian averaged cross sections from TALYS files
!
  if (flagcross) call xsread
  if (flagresidual) call resread
  if (flagprod) call prodread
  if (flaggamma) call gamread
  Nchanall = Nchanxs + Nchanrp + Nchangam + Nchanprod
  Nchanall = min(Nchanxs + Nchanrp + Nchangam + Nchanprod, numchanxs)
  if (flagspectra) call specread
  if (flagangle) call angleread
  if (flagleg) call legread
  if (flagintegral) call integralread
  if (flagproduct) call productread
  if (flagmacs) call macsread
  if (flaggamgam) call gamgamread
  if (flagpsf) call psfread
  if (flagld) call ldread
  if (italys == Ntalbeg) then
    if (flagtal) call talread
    if (flaglib) call libread
    if (flagexp) call expread
  endif
  if (flagmagnet) call magnet
  return
end subroutine talysread
! Copyright A.J. Koning 2021
