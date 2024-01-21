subroutine gamread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read gamma cross sections from TALYS files
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
!   numchangam    ! maximum number of channels with gamma cross sections
!   numenin       ! maximum number of incident energies
!   numenS        ! maximum number of energies for sensitivities
!   numpar        ! maximum number of parameters
! Variables for writing TALYS input files
!   italys        ! TALYS run counter
!   mode          ! TASMAN mode
!   Nhigh         ! number of high energy runs
! Variables for processing input
!   Atarget       ! mass number of target nucleus
!   Ztarget       ! charge number of target nucleus
! Variables for sensitivity
!   flagsens      ! flag to use sensitivity in parameter weighting
! Variables for reading cross sections
!   E             ! incident energy
! Variables for reading gamma production cross sections
!   Egam          ! incident energy
!   maxgam        ! maximum number of discrete levels with gamma cross section
!   gamfile       ! name of file with gamma production cross sections
!   gamsave       ! gamma cross section from TALYS
!   gamtalys      ! gamma cross section from TALYS
!   Nchangam      ! total number of channels with gamma production cross sections
!   Nengam        ! number of incident energies
!   Nengam0       ! number of incident energies
!   Sgamindex     ! index for gamma cross section
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist                  ! logical to determine existence
  character(len=19) :: xsfg                    ! gamma file
  character(len=132):: gamstring(numenin+100)    ! gamma string
  character(len=132):: line
  character(len=132):: key
  integer           :: i                       ! counter
  integer           :: ia                      ! mass number of nucleus
  integer           :: ichan                   ! counter for channels
  integer           :: isamp                   ! variable for 0th or random run
  integer           :: istat                   ! logical for file access
  integer           :: iz                      ! charge number of nucleus
  integer           :: keyix
  integer           :: j                       ! counter
  integer           :: jheader                 ! counter
  integer           :: k                       ! Legendre order
  integer           :: L1                      ! integration limits
  integer           :: L2                      ! integration limits
  integer           :: n                       ! counter
  real(sgl)         :: cratio                  ! cross section ratio
  real(sgl)         :: e1                      ! first energy
  real(sgl)         :: Rs                      ! sensitivity factor
  real(sgl)         :: xs1                     ! help variable
!
! ************************ Set channels ********************************
!
! Gamma production
!
!
  if (Nchangam == 0) then
    ichan = 0
    xsfg = 'gam000000L00L00.tot'
    do iz = Ztarget - 2, Ztarget + 2
      do ia = Atarget - 4, Atarget + 4
        write(xsfg(4:6), '(i3.3)') iz
        write(xsfg(7:9), '(i3.3)') ia
        do L1 = 1, maxgam
          do L2 = 0, L1 - 1
            write(xsfg(11:12), '(i2.2)') L1
            write(xsfg(14:15), '(i2.2)') L2
            inquire (file = xsfg, exist = lexist)
            if (lexist) then
              ichan = ichan + 1
              if (ichan <= numchangam) then
                gamfile(ichan) = xsfg
                Nchangam = ichan
              endif
            endif
          enddo
        enddo
      enddo
    enddo
  endif
!
! Read gamma production cross sections
!
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  do i = 1, Nchangam
    gamstring = ''
    if (italys > Nhigh) then
      do k = 1, Nengam0(i)
        gamtalys(isamp, i, k) = gamtalys(0, i, k)
      enddo
    endif
    inquire (file = gamfile(i), exist = lexist)
    if (lexist) then
      open (unit = 2, file = gamfile(i), status = 'old', iostat = istat)
      if (istat /= 0) call read_error(gamfile(i), istat)
      Egam(i, 0) = 0.
      k = 0
      do
        read(2,'(a)',iostat = istat) line
        if (istat == -1) exit
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(2,'(/)')
          do
            read(2, * , iostat = istat) e1, xs1
            if (istat == -1) exit
            if (istat /= 0) call read_error(gamfile(i), istat, eor = 'continue')
            k = k + 1
            Egam(i, k) = e1
            gamtalys(isamp, i, k) = xs1
          enddo
          exit
        endif
      enddo
      Nengam(i) = k
      close (2)
      if (italys <= Nhigh) then
        Nengam0(i) = k
      else
        if (gamtalys(0, i, k) > 0.) then
          cratio = gamtalys(isamp, i, k) / gamtalys(0, i, k)
        else
          cratio = 1.
        endif
        cratio = min(cratio, 10.)
        do j = k + 1, Nengam0(i)
          gamtalys(isamp, i, j) = cratio * gamtalys(isamp, i, j)
        enddo
        open (unit = 2, file = gamfile(i), status = 'old', iostat = istat)
        if (istat /= 0) call read_error(gamfile(i), istat)
        j = 1
        do
          read(2, '(a)', iostat = istat) gamstring(j)
          if (istat == -1) exit
          if (istat /= 0) call read_error(gamfile(i), istat)
          key='entries'
          keyix=index(gamstring(j),trim(key))
          if (keyix > 0) then
            write(gamstring(j)(keyix+len_trim(key)+2:80),*, iostat = istat) Nengam0(i)
            jheader = j + 2
          endif
          j = j + 1
        enddo
        close (2)
        open (unit = 2, file = gamfile(i), status = 'replace')
        do j = 1, jheader + Nengam(i) 
          write(2, '(a)') trim(gamstring(j))
        enddo
        do j = Nengam(i) + 1, Nengam0(i)
          write(2, '(2es15.6)') E(i, j), gamtalys(isamp, i, j)
        enddo
        close (2)
        Nengam(i) = Nengam0(i)
      endif
    endif
  enddo
!
! Create energy grid and cross sections for sensitivities
!
  if (flagsens .or. mode == 2) then
    if (italys == 0) then
      do i = 1, Nchangam
        if (Nengam(i) <= numenS) then
          do k = 1, Nengam(i)
            Sgamindex(i, k) = k
          enddo
        else
          Rs = real(numenS) / real(Nengam(i)) + 0.001
          do k = 1, Nengam(i)
            Sgamindex(i, k) = int(k * Rs)
          enddo
        endif
      enddo
    endif
    if (italys <= numpar) then
      do i = 1, Nchangam
        do k = 1, Nengam(i)
          n = Sgamindex(i, k)
          gamsave(italys, i, n) = gamtalys(isamp, i, k)
        enddo
      enddo
    endif
  endif
  return
end subroutine gamread
! Copyright A.J. Koning 2021
