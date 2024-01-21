subroutine specread
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
!   sgl          ! single precision kind
! All global variables
!   numchansp    ! maximum number of channels with emission spectra
!   numen2       ! number of emission energies
! Constants
!   particle     ! type of particle
! Variables for writing TALYS input files
!   italys       ! TALYS run counter
! Variables for processing input
!   numE         ! number of incident energies
! Variables for reading spectra
!   Eout         ! emission energy
!   Nchansp      ! total number of channels with emission spectra
!   Nensp        ! number of covariance energies
!   spfile       ! name of file with emission spectra
!   sptalys      ! emission spectrum from TALYS
! Variables for reading experimental data
!   Ein          ! incident energy
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
  integer            :: k        ! designator for particle
  integer            :: nin      ! counter for incident energy
  real(sgl)          :: e1       ! first energy
  real(sgl)          :: Efile    !
  real(sgl)          :: xs1      ! help variable
!
! ************************ Set channels ********************************
!
! Emission spectra
!
  if (Nchansp == 0) then
    ichan = 0
    do nin = 1, numE
      ichan = ichan + 1
      if (ichan <= numchansp) then
        Espec(ichan) = Ein(nin)
        Nchansp = ichan
      endif
    enddo
  endif
!
! Read emission spectra
!
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  if (flagblock) then
    xsfa = ptype0//'spec.tot'
    inquire (file = xsfa, exist = lexist)
    if (lexist) then
      spfile(1) = xsfa
      do i = 1, Nchansp
        open (unit = 2, file = xsfa, status = 'old')
        do
          read(2,'(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(spfile(1), istat)
          key='E-incident [MeV]'
          keyix=index(line,trim(key))
          if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
          if (istat /= 0) call read_error(spfile(1), istat)
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Nensp(i)
            if (istat /= 0) call read_error(spfile(1), istat)
            read(2,'(/)')
            if (abs(Espec(i) - Efile) < 1.e-5) then
              do k = 1, Nensp(i)
                read(2, * , iostat = istat) e1, xs1
                if (istat ==  -1) exit
                if (istat /= 0) call read_error(xsfa, istat)
                Eout(i, k) = e1
                sptalys(isamp, i, k) = xs1
              enddo
              exit
            else
              do k = 0, Nensp(i) - 1
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
    do i = 1, Nchansp
      xsfa = ptype0//'spec0000.000.tot'
      write(xsfa(6:13), '(f8.3)') Espec(i)
      write(xsfa(6:9), '(i4.4)') int(Espec(i))
      spfile(i) = xsfa
      inquire (file = spfile(i), exist = lexist)
      if (lexist) then
        open (unit = 2, file = spfile(i), status = 'old')
        do
          read(1,'(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(spfile(i), istat)
          key='E-incident [MeV]'
          keyix=index(line,trim(key))
          if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
          if (istat /= 0) call read_error(spfile(i), istat)
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(1,'(/)')
            k = 0
            do
              read(2, * , iostat = istat) e1, xs1
              if (istat == -1) exit
              if (istat /= 0) call read_error(spfile(i), istat)
              k = k + 1
              if (k > numen2) then
                k = k - 1
                exit
              endif
              Eout(i, k) = e1
              sptalys(isamp, i, k) = xs1
            enddo
            Nensp(i) = k
            exit
          endif
        enddo
        close (2)
      endif
    enddo
  endif
  return
end subroutine specread
