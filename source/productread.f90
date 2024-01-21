subroutine productread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read isotope production from TALYS files
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
!   sgl             ! single precision kind
! All global variables
!   numchanY        ! maximum number of channels with isotope production
! Variables for writing TALYS input files
!   italys          ! TALYS run counter
! Variables for processing input
!   Atarget         ! mass number of target nucleus
!   Ztarget         ! charge number of target nucleus
! Variables for reading isotope production
!   acttalys        ! activity of produced isotope in MBq
!   NchanY          ! total number of channels with isotope production
!   Nisoreltalys    ! fraction of number of produced isotopes per element
!   Nisotalys       ! number of isotopes produced after irradiation
!   Ntime           ! number of time points
!   timeY           ! time point
!   Yfile           ! name of file with isotope production
!   yieldtalys      ! yield of produced isotope in MBq / (mA.h)
! Variables for reading cross sections
!   flagisomer      ! flag for isomeric production cross sections
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist    ! logical to determine existence
  character(len=8)  :: xsf       ! filename
  character(len=12) :: xsfL      ! filename
  character(len=12) :: xsfT      ! filename
  character(len=132) :: line
  character(len=132) :: key
  integer           :: i         ! counter
  integer           :: ia        ! mass number of nucleus
  integer           :: ichan     ! counter for channels
  integer           :: isamp     ! variable for 0th or random run
  integer           :: istat     ! logical for file access
  integer           :: iz        ! charge number of nucleus
  integer           :: keyix
  integer           :: k         ! Legendre order
  integer           :: L         ! counter for Legendre coefficients
  real(sgl)         :: t1        ! help variable
  real(sgl)         :: x1        ! coordinates of intersection points inside the bin
  real(sgl)         :: x2        ! coordinates of the 2nd summit of the triangle
  real(sgl)         :: x3        ! coordinates of the 3rd summit of the triangle
  real(sgl)         :: x4        ! cosine of the angle
!
! ************************ Set channels ********************************
!
! Isotope production
!
!   cross sections
!
  if (NchanY == 0) then
    ichan = 0
    xsf = 'Y000000.'
    do ia = Atarget + 4, Atarget - 34, -1
      do iz = Ztarget + 2, Ztarget - 14, -1
        write(xsf(2:4), '(i3.3)') iz
        write(xsf(5:7), '(i3.3)') ia
        xsfT = xsf//'tot'
        inquire (file = xsfT, exist = lexist)
        if (lexist) then
          ichan = ichan + 1
          if (ichan <= numchanY) then
            Yfile(ichan) = xsfT
            NchanY = ichan
          endif
        endif
        if (flagisomer) then
          xsfL = xsf//'L00'
          do L = 0, 40
            write(xsfL(10:11), '(i2.2)') L
            inquire (file = xsfL, exist = lexist)
            if (lexist) then
              ichan = ichan + 1
              if (ichan <= numchanY) then
                Yfile(ichan) = xsfL
                NchanY = ichan
              endif
            endif
          enddo
        endif
      enddo
    enddo
  endif
!
! Read isotope production
!
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  do i = 1, NchanY
    inquire (file = Yfile(i), exist = lexist)
    if (lexist) then
      open (unit = 2, file = Yfile(i), status = 'old', iostat = istat)
      if (istat /= 0) call read_error(Yfile(i), istat)
      timeY(i, 0) = 0.
      k = 0
      do
        read(2,'(a)',iostat = istat) line
        if (istat == -1) exit
          key='entries'
          keyix=index(line,trim(key))
          if (keyix > 0) then
          read(2,'(/)')
          do
            read(2, * , iostat = istat) t1, x1, x2, x3, x4
            if (istat == -1) exit
            if (istat /= 0) call read_error(Yfile(i), istat)
            k = k + 1
            timeY(i, k) = t1
            acttalys(isamp, i, k) = x1
            Nisotalys(isamp, i, k) = x2
            yieldtalys(isamp, i, k) = x3
            Nisoreltalys(isamp, i, k) = x4
          enddo
          exit
        endif
      enddo
      Ntime(i) = k
      close (2)
    endif
  enddo
  return
end subroutine productread
! Copyright A.J. Koning 2021
