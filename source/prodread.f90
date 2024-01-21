 subroutine prodread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read particle production cross sections from TALYS files
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
!   sgl            ! single precision kind
! All global variables
!   numchanprod    ! maximum number of channels with particle production cross sections
!   numenin        ! maximum number of incident energies
!   numenS         ! maximum number of energies for sensitivities
!   numpar         ! maximum number of parameters
! Variables for writing TALYS input files
!   italys         ! TALYS run counter
!   mode           ! TASMAN mode
!   Nhigh          ! number of high energy runs
! Constants
!   particle       ! type of particle
! Variables for sensitivity
!   flagsens       ! flag to use sensitivity in parameter weighting
! Variables for reading production cross sections
!   Eprod          ! incident energy
!   Nchanprod      ! number of channels with particle production cross sections
!   Nenprod        ! number of covariance energies
!   Nenprod0       ! number of covariance energies
!   prodfile       ! name of file with particle production cross sections
!   prodsave       ! particle production cross section from TALYS
!   prodtalys      ! particle production cross section from TALYS
!   Sprodindex     ! index of particle production cross sections
!   Ytalys         ! yield from TALYS
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist                   ! logical to determine existence
  character(len=9)  :: xsf                      ! filename
  character(len=132):: prodstring(numenin+100)    ! string for file with production cross section
  character(len=132):: line
  character(len=132):: key
  integer           :: i                        ! counter
  integer           :: ichan                    ! counter for channels
  integer           :: isamp                    ! variable for 0th or random run
  integer           :: istat                    ! logical for file access
  integer           :: keyix
  integer           :: j                        ! counter
  integer           :: jheader                  ! counter
  integer           :: k                        ! Legendre order
  integer           :: n                        ! counter
  integer           :: type                     ! particle type
  real(sgl)         :: cratio                   ! cross section ratio
  real(sgl)         :: e1                       ! first energy
  real(sgl)         :: Rs                       ! sensitivity factor
  real(sgl)         :: xs1                      ! help variable
  real(sgl)         :: xs2                      ! help variable
  real(sgl)         :: yratio                   ! ratio of yields
!
! ************************ Set channels ********************************
!
! Particle production
!
  if (Nchanprod == 0) then
    ichan = 0
    do type = 0, 6
      xsf = ' prod.tot'
      write(xsf(1:1), '(a1)') particle(type)
      inquire (file = xsf, exist = lexist)
      if (lexist) then
        ichan = ichan + 1
        if (ichan <= numchanprod) then
          prodfile(ichan) = xsf
          Nchanprod = ichan
        endif
      endif
    enddo
  endif
!
! Read particle production cross sections
!
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  do i = 1, Nchanprod
    prodstring = ''
    if (italys > Nhigh) then
      do k = 1, Nenprod0(i)
        prodtalys(isamp, i, k) = prodtalys(0, i, k)
        Ytalys(isamp, i, k) = Ytalys(0, i, k)
      enddo
    endif
    inquire (file = prodfile(i), exist = lexist)
    if (lexist) then
      open (unit = 2, file = prodfile(i), status = 'old', iostat = istat)
      if (istat /= 0) call read_error(prodfile(i), istat)
      Eprod(i, 0) = 0.
      k = 0
      do
        read(2,'(a)',iostat = istat) line
        if (istat == -1) exit
        key='entries'
        keyix=index(line,trim(key))
        if (keyix > 0) then
          read(2,'(/)')
          do
            read(2, * , iostat = istat) e1, xs1, xs2
            if (istat == -1) exit
            if (istat /= 0) call read_error(prodfile(i), istat, error = 'continue')
            k = k + 1
            Eprod(i, k) = e1
            prodtalys(isamp, i, k) = xs1
            Ytalys(isamp, i, k) = xs2
          enddo
          exit
        endif
      enddo
      Nenprod(i) = k
      close (2)
      if (italys <= Nhigh) then
        Nenprod0(i) = k
      else
        if (prodtalys(0, i, k) > 0.) then
          cratio = prodtalys(isamp, i, k) / prodtalys(0, i, k)
        else
          cratio = 1.
        endif
        cratio = min(cratio, 10.)
        if (Ytalys(0, i, k) > 0.) then
          yratio = Ytalys(isamp, i, k) / Ytalys(0, i, k)
        else
          yratio = 1.
        endif
        do j = k + 1, Nenprod0(i)
          prodtalys(isamp, i, j) = cratio * prodtalys(isamp, i, j)
          Ytalys(isamp, i, j) = yratio * Ytalys(isamp, i, j)
        enddo
        open (unit = 2, file = prodfile(i), status = 'old', iostat = istat)
        if (istat /= 0) call read_error(prodfile(i), istat)
        jheader = 0
        j = 1
        do
          read(2, '(a)', iostat = istat) prodstring(j)
          if (istat == -1) exit
          if (istat /= 0) call read_error(prodfile(i), istat)
          key='entries'
          keyix=index(prodstring(j),trim(key))
          if (keyix > 0) then
            read(prodstring(j)(keyix+len_trim(key)+2:80),*, iostat = istat) Nenprod0(i)
            jheader = j + 2
          endif
          j = j + 1
        enddo
        close (2)
        open (unit = 2, file = prodfile(i), status = 'replace')
        do j = 1, Nenprod(i) + jheader
          write(2, '(a)') trim(prodstring(j))
        enddo
        do j = Nenprod(i) + 1, Nenprod0(i)
          write(2, '(3es15.6)') Eprod(i, j), prodtalys(isamp, i, j), Ytalys(isamp, i, j)
        enddo
        close (2)
        Nenprod(i) = Nenprod0(i)
      endif
    endif
  enddo
!
! Create energy grid and cross sections for sensitivities
!
  if (flagsens .or. mode == 2) then
    if (italys == 0) then
      do i = 1, Nchanprod
        if (Nenprod(i) <= numenS) then
          do k = 1, Nenprod(i)
            Sprodindex(i, k) = k
          enddo
        else
          Rs = real(numenS) / real(Nenprod(i)) + 0.001
          do k = 1, Nenprod(i)
            Sprodindex(i, k) = int(k * Rs)
          enddo
        endif
      enddo
    endif
    if (italys <= numpar) then
      do i = 1, Nchanprod
        do k = 1, Nenprod(i)
          n = Sprodindex(i, k)
          prodsave(italys, i, n) = prodtalys(isamp, i, k)
        enddo
      enddo
    endif
  endif
  return
end subroutine prodread
! Copyright A.J. Koning 2021
