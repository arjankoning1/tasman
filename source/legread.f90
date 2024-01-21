subroutine legread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read Legendre coefficients from TALYS files
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
!   sgl           ! single precision kind
! All global variables
!   numleg        ! maximum number of Legendre coefficients
! Variables for writing TALYS input files
!   italys        ! TALYS run counter
! Variables for processing input
!   anginclude    ! flag to include energy point for angular grid
!   numE          ! number of incident energies
! Variables for reading Legendre coefficients
!   Eleg          ! energy grid for Legendre coefficients
!   leg0          ! Legendre coefficients from TALYS
!   legfile       ! name of file with Legendre coefficients
!   legtalys      ! Legendre coefficients from TALYS
!   Nchanleg      ! total number of channels with Legendre coefficients
! Variables for reading experimental data
!   Ein           ! incident energy
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical            :: lexist   ! logical to determine existence
  character(len=17)  :: xsfa     ! filename
  character(len=132) :: line     ! input line
  character(len=132) :: key
  integer            :: i        ! level
  integer            :: keyix
  integer            :: ichan    ! counter for channels
  integer            :: isamp    ! variable for 0th or random run
  integer            :: istat    ! logical for file access
  integer            :: k        ! counter
  integer            :: kdum     ! help variable
  integer            :: nin      ! counter for incident energy
  integer            :: Nleg     ! number of Legendre coefficients
  real(sgl)          :: Efile    ! energy
  real(sgl)          :: leg1     ! first Legendre coefficient
  real(sgl)          :: leg2     ! second Legendre coefficient
  real(sgl)          :: xs1      ! help variable
!
! ************************ Set channels ********************************
!
! Legendre coefficients
!
  if (Nchanleg == 0) then
    ichan = 0
    do nin = 1, numE
      if (anginclude(nin)) then
        ichan = ichan + 1
        Eleg(ichan) = Ein(nin)
      endif
    enddo
    Nchanleg = ichan
  endif
!
! Read Legendre coefficients
!
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  if (flagblock) then
    xsfa = ptype0//ptype0//'leg.L00'
    inquire (file = xsfa, exist = lexist)
    if (lexist) then
      legfile(1) = xsfa
      do i = 1, Nchanleg
        open (unit = 2, file = xsfa, status = 'old')
        do
          read(2,'(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(legfile(1), istat)
          key='E-incident [MeV]'
          keyix=index(line,trim(key))
          if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
          if (istat /= 0) call read_error(legfile(1), istat)
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Nleg
            if (istat /= 0) call read_error(legfile(1), istat)
            read(2,'(/)')
            if (abs(Eleg(i) - Efile) < 1.e-5) then
              do k = 0, Nleg - 1
                read(2, * , iostat = istat) kdum, xs1, xs1, xs1, leg1, leg2
                if (istat ==  -1) exit
                if (istat /= 0) call read_error(legfile(i), istat)
                if (k < numleg) then
                  legtalys(isamp, i, k) = leg1
                  leg0(isamp, i, k) = leg2
                endif
              enddo
              exit
            else
              do k = 0, Nleg - 1
                read(2, '()')
              enddo
              cycle
            endif
          endif
        enddo
        close (2)
      enddo
    endif
  else
    do i = 1, Nchanleg
      xsfa = ptype0//ptype0//'0000.000leg.L00'
      write(xsfa(3:10), '(f8.3)') Eleg(i)
      write(xsfa(3:6), '(i4.4)') int(Eleg(i))
      legfile(i) = xsfa
      inquire (file = legfile(i), exist = lexist)
      if (lexist) then
        open (unit = 2, file = legfile(i), status = 'old')
        do
          read(1,'(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(legfile(i), istat)
          key='E-incident [MeV]'
          keyix=index(line,trim(key))
          if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
          if (istat /= 0) call read_error(legfile(i), istat)
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(1,'(/)')
            do
              read(2, * , iostat = istat) k, xs1, xs1, xs1, leg1, leg2
              if (istat == -1) exit
              if (istat /= 0) call read_error(legfile(i), istat)
              legtalys(isamp, i, k) = leg1
              leg0(isamp, i, k) = leg2
              if (k >= numleg) exit
            enddo
            exit
          endif
        enddo
        close (2)
      endif
    enddo
  endif
  return
end subroutine legread
! Copyright A.J. Koning 2021
