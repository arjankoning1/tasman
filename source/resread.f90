subroutine resread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read residual cross sections from TALYS files
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
!   numchanrp     ! maximum number of channels with residual production cross sections
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
! Variables for reading cross sections
!   flagisomer    ! flag for isomeric production cross sections
! Variables for reading residual production cross sections
!   Arp           ! mass number of residual product
!   Erp           ! incident energy
!   Lrp           ! isomeric level of residual product
!   Nchanrp       ! total number of channels with residual production cross sections
!   Nenrp         ! number of incident energies
!   Nenrp0        ! number of incident energies
!   rpfile        ! name of file with residual production cross sections
!   rpsave        ! residual production cross section from TALYS
!   rptalys       ! residual production cross section from TALYS
!   Srpindex      ! index for residual production cross section
!   Zrp           ! charge number of residual product
! Variables for sensitivity
!   flagsens      ! flag to use sensitivity in parameter weighting
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist                  ! logical to determine existence
  character(len=9)  :: xsf                     ! filename
  character(len=13) :: xsfT                    ! filename
  character(len=13) :: xsfL                    ! filename
  character(len=132):: resstring(numenin+100)    ! string for residual product
  character(len=132) :: line               ! line
  character(len=132) :: key
  integer           :: i                       ! counter
  integer           :: ia                      ! mass number of nucleus
  integer           :: ichan                   ! counter for channels
  integer           :: isamp                   ! variable for 0th or random run
  integer           :: istat                   ! logical for file access
  integer           :: iz                      ! charge number of nucleus
  integer           :: j                       ! counter
  integer           :: keyix
  integer           :: jheader                 ! counter
  integer           :: k                       ! Legendre order
  integer           :: L                       ! counter for Legendre coefficients
  integer           :: n                       ! counter
  real(sgl)         :: cratio                  ! cross section ratio
  real(sgl)         :: e1                      ! first energy
  real(sgl)         :: Rs                      ! sensitivity factor
  real(sgl)         :: xs1                     ! help variable
!
! ************************ Set channels ********************************
!
! Residual production
!
  if (Nchanrp == 0) then
    ichan = 0
    xsf = 'rp000000.'
    do ia = Atarget + 4, Atarget - 34, -1
      do iz = Ztarget + 2, Ztarget - 14, -1
        write(xsf(3:5), '(i3.3)') iz
        write(xsf(6:8), '(i3.3)') ia
        xsfT = xsf//'tot'
        inquire (file = xsfT, exist = lexist)
        if (lexist) then
          ichan = ichan + 1
          if (ichan <= numchanrp) then
            Zrp(ichan) = iz
            Arp(ichan) = ia
            Lrp(ichan) = 0
            rpfile(ichan) = xsfT
            Nchanrp = ichan
          endif
        endif
        if (flagisomer) then
          xsfL = xsf//'L00'
          do L = 0, 40
            write(xsfL(11:12), '(i2.2)') L
            inquire (file = xsfL, exist = lexist)
            if (lexist) then
              ichan = ichan + 1
              if (ichan <= numchanrp) then
                Zrp(ichan) = iz
                Arp(ichan) = ia
                Lrp(ichan) = L
                rpfile(ichan) = xsfL
                Nchanrp = ichan
              endif
            endif
          enddo
        endif
      enddo
    enddo
  endif
!
! Read residual production cross sections
!
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  do i = 1, Nchanrp
    resstring = ''
    if (italys > Nhigh) then
      do k = 1, Nenrp0(i)
        rptalys(isamp, i, k) = rptalys(0, i, k)
      enddo
    endif
    inquire (file = rpfile(i), exist = lexist)
    if (lexist) then
      open (unit = 2, file = rpfile(i), status = 'old', iostat = istat)
      if (istat /= 0) call read_error(rpfile(i), istat)
      k = 0
      Erp(i, 0) = 0.
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
            if (istat /= 0) call read_error(rpfile(i), istat, eor = 'continue')
            k = k + 1
            Erp(i, k) = e1
            rptalys(isamp, i, k) = xs1
          enddo
          exit
        endif
      enddo
      close (2)
      Nenrp(i) = k
      if (italys <= Nhigh) then
        Nenrp0(i) = k
      else
        if (rptalys(0, i, k) > 0.) then
          cratio = rptalys(isamp, i, k) / rptalys(0, i, k)
        else
          cratio = 1.
        endif
        cratio = min(cratio, 10.)
        do j = k + 1, Nenrp0(i)
          rptalys(isamp, i, j) = cratio * rptalys(isamp, i, j)
        enddo
        open (unit = 2, file = rpfile(i), status = 'old', iostat = istat)
        if (istat /= 0) call read_error(rpfile(i), istat)
        jheader = 0
        j = 1
        do
          read(2, '(a)', iostat = istat) resstring(j)
          if (istat == -1) exit
          if (istat /= 0) call read_error(rpfile(i), istat)
          key='entries'
          keyix=index(resstring(j),trim(key))
          if (keyix > 0) then
            read(resstring(j)(keyix+len_trim(key)+2:80),*, iostat = istat) Nenrp0(i)
            jheader = j + 2
          endif
          j = j + 1
        enddo
        close (2)
        open (unit = 2, file = rpfile(i), status = 'replace')
        do j = 1, jheader + Nenrp(i) 
          write(2, '(a)') trim(resstring(j))
        enddo
        do j = Nenrp(i) + 1, Nenrp0(i)
          write(2, '(2es15.6)') Erp(i, j), rptalys(isamp, i, j)
        enddo
        close (2)
        Nenrp(i) = Nenrp0(i)
      endif
    endif
  enddo
!
! Create energy grid and cross sections for sensitivities
!
  if (flagsens .or. mode == 2) then
    if (italys == 0) then
      do i = 1, Nchanrp
        if (Nenrp(i) <= numenS) then
          do k = 1, Nenrp(i)
            Srpindex(i, k) = k
          enddo
        else
          Rs = real(numenS) / real(Nenrp(i)) + 0.001
          do k = 1, Nenrp(i)
            Srpindex(i, k) = int(k * Rs)
          enddo
        endif
      enddo
    endif
    if (italys <= numpar) then
      do i = 1, Nchanrp
        do k = 1, Nenrp(i)
          n = Srpindex(i, k)
          rpsave(italys, i, n) = rptalys(isamp, i, k)
        enddo
      enddo
    endif
  endif
  return
end subroutine resread
! Copyright A.J. Koning 2021
