subroutine uncertainty
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Uncertainty propagation and sensitivities
!
! Author    : Arjan Koning
!
! 2025-12-25: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
!
! Variables for uncertainty
!   flagsave        ! flag to save all input and output files of each Talys run
! Variables for reading TALYS output
!   flagcross       ! flag for covariances or optimization of cross sections
!   flagexp         ! flag to optimize to experimental data
!   flaglib         ! flag to optimize to nuclear data libraries
!   flagtal         ! flag to optimize to TALYS
!   flagtalys       ! flag to track successful TALYS calculation
!   Ntalbeg         ! counter for first TALYS calculation
! Variables for writing TALYS input
!   flaginponly     ! flag to produce only TALYS input files, no runs
!   italys          ! TALYS run counter
!   mode            ! TASMAN mode
! Variables for processing input
!   talys           ! TALYS executable
! Variables for weights
!   flagweight      ! flag to use weights for random samples
! Variables for GOF function
!   Ntalys          ! number of TALYS runs
! Variables for sensitivity
!   flagreadsens    ! flag to read sensitivities from tables
!   flagsens        ! flag to use sensitivity in parameter weighting
! Variables for TMC
!   tmcoffset       ! offset for starting creation of ENDF - 6 files (TMC only)
!
! *** Declaration of local data
!
  implicit none
  character(len=4)   :: calcnumstr! string for counter
  character(len=132) :: talcmd    ! TALYS command
  integer            :: i         ! counter
  integer            :: system    ! system call
  integer            :: Nloop     ! number of loops, > 1 for Morris screening
!
! ******************** Uncertainties and sensitivities *****************
!
! fitlimits   : subroutine to set parameter limits for optimization
! covinitial  : subroutine for initialization for covariance and sensitivity calculations
! readsens    : subroutine to read sensitivity matrix
! parvariation: subroutine for variation of model parameters
! inputwrite  : subroutine to write parameters for TALYS input file
!
! Initialization, parameter variation and running of TALYS
!
  call fitlimits
  call covinitial
  if (flagweight .and. flagreadsens) call readsens
  if (flagextparvar) open(unit=80,file='parvars.inp',status='old')
  if (flagmorris) then
    Nloop = Nmorris
  else
    Nloop = 1
  endif
  do iloop = 1, Nloop
    if (iloop > 1) write(*,'(/," Random loop :",i6)') iloop
    do italys = Ntalbeg, Ntalys
      if ((mode == 1 .and. italys > 0) .or. (mode == 2 .and. flagmorris .and. italys == 0 .and. iloop > 1)) call parvariation
      if (mode == 2) call sensinput
      call inputwrite
      if (flaginponly) cycle
!
! Patch by Georg Schnabel to read pre-calculated TALYS results
!
      write(calcnumstr, '(i4.4)') italys + 1
      if (trim(getcalcscript) /= 'none') then
        talcmd = trim(getcalcscript)// ' ' // calcnumstr
        write(*,*) trim(talcmd)
        i = system(talcmd)
      else
        talcmd = trim(talys) //' < talys.inp > talys.out'
        i = system(talcmd)
      endif
!
! Read and process output of TALYS
!
! talysread: subroutine to read results from TALYS files
!
      call talysread
      if ( .not. flagtalys) cycle
      if (flagsave) call output
      if (italys ==  -1) cycle
!
! Statistical processing
!
! gof        : subroutine for goodness-of-fit TALYS run
! output     : subroutine to write output files (TMC only)
! covariance : subroutine for creation of covariance matrices and average cross sections
! sensitivity: subroutine to create sensitivity matrix
! cleaner    : subroutine to clean up TALYS files before next random run
!
      if (italys >= 0) then
        if (flagexp .or. flaglib .or. flagtal) call gof
      endif
      if (italys > 0) then
        if (mode == 1 .and. italys > tmcoffset .and. flagcross) call covariance
        if ((flagsens .and. .not. flagreadsens) .or. mode == 2) call sensitivity
      endif
    enddo
  enddo
  if (flagextparvar) close(unit=80)
  return
end subroutine uncertainty
! Copyright A.J. Koning 2021
