subroutine angleread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read angular distributions from TALYS files
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
! Variables for writing TALYS input
!   italys        ! TALYS run counter
! Variables for processing input
!   anginclude    ! flag to include energy point for angular grid
!   numE          ! number of incident energies
! Variables for reading experimental data
!   Ein           ! incident energy
! Variables for reading angular distributions
!   angfile       ! name of file with angular distributions
!   angle         ! angle
!   angtalys      ! angular distribution
!   Nang          ! number of angles
!   Nchanang      ! total number of channels with angular distributions
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
  integer            :: Na       ! help variable
  integer            :: nin      ! counter for incident energy
  real(sgl)          :: ang      ! angle
  real(sgl)          :: Efile    ! help variable
  real(sgl)          :: xs1      ! help variable
!
! ************************ Set channels ********************************
!
! Angular distributions
!
  if (Nchanang == 0) then
    ichan = 0
    do nin = 1, numE
      if (anginclude(nin)) then
        ichan = ichan + 1
        Eang(ichan) = Ein(nin)
      endif
    enddo
    Nchanang = ichan
  endif
!
! Read angular distributions
!
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  if (flagblock) then
    xsfa = ptype0//ptype0//'ang.L00'
    inquire (file = xsfa, exist = lexist)
    if (lexist) then
      angfile(1) = xsfa
      do i = 1, Nchanang
        open (unit = 2, file = xsfa, status = 'old')
        do
          read(2,'(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(angfile(1), istat)
          key='E-incident [MeV]'
          keyix=index(line,trim(key))
          if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
          if (istat /= 0) call read_error(angfile(1), istat)
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Na
            Nang(i) = Na - 1
            if (istat /= 0) call read_error(angfile(1), istat)
            read(2,'(/)')
            if (abs(Eang(i) - Efile) < 1.e-5) then
              do k = 0, Nang(i)
                read(2, * , iostat = istat) ang, xs1
                if (istat ==  -1) exit
                if (istat /= 0) call read_error(angfile(i), istat, error = 'continue')
                if (k < numang) then
                  angle(i, k) = ang
                  angtalys(isamp, i, k) = xs1
                endif
              enddo
              exit
            else
              do k = 0, Nang(i)
                read(2, '()', iostat=istat)
                if (istat ==  -1) exit
              enddo
              cycle
            endif
          endif
        enddo
        close (2)
      enddo
    endif
  else
    do i = 1, Nchanang
      xsfa = ptype0//ptype0//'0000.000ang.L00'
      write(xsfa(3:10), '(f8.3)') Eang(i)
      write(xsfa(3:6), '(i4.4)') int(Eang(i))
      angfile(i) = xsfa
      inquire (file = angfile(i), exist = lexist)
      if (lexist) then
        open (unit = 2, file = angfile(i), status = 'old')
        do
          read(1,'(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(angfile(i), istat)
          key='E-incident [MeV]'
          keyix=index(line,trim(key))
          if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Efile
          if (istat /= 0) call read_error(angfile(i), istat)
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
            read(2,'(/)')
            k = -1
            do
              read(2, * , iostat = istat) ang, xs1
              if (istat == -1) exit
              if (istat /= 0) call read_error(angfile(i), istat, error = 'continue')
              k = k + 1
              angle(i, k) = ang
              angtalys(isamp, i, k) = xs1
            enddo
            exit
          endif
        enddo
        close (2)
      endif
    enddo
  endif
  return
end subroutine angleread
! Copyright A.J. Koning 2021
