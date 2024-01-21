subroutine output
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Write output files
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
! All global variables
!   numleg          ! maximum number of Legendre coefficients
! Variables for reading TALYS output
!   flagcross       ! flag for covariances or optimization of cross sections
!   flaggamma       ! flag for covariances or optimization of discrete gamma ray transitions
!   flagintegral    ! flag for covariances or optimization of integral activation data
!   flagmacs        ! flag for covariances or optimization of Maxwellian averaged cross sections
!   flagprod        ! flag for covariances or optimization of particle production cross sections
!   flagresidual    ! flag for covariances or optimization of residual production cross sections
! Variables for writing TALYS input files
!   flagangle       ! flag for covariances or optimization of angular distributions
!   flagleg         ! flag for covariances or optimization of Legendre coefficients
!   flagspectra     ! flag for covariances or optimization of emission spectra
!   flagtmc         ! flag for Total Monte Carlo (ENDF - 6 data file at each run)
!   italys          ! TALYS run counter
!   mode            ! TASMAN mode
! Variables for TMC
!   offset          ! offset for numbering or random files (TMC only)
! Variables for weights
!   Esamp           ! energy of sampled cross section
!   flagchmagnet    ! flag to unify theoretical and experimental data
!   Nxssamp         ! number of sampled cross sections
!   xssamp          ! sampled random cross section
! Variables for reading cross sections
!   E               ! incident energy
!   Nchanxs         ! total number of channels with cross sections
!   Nen             ! number of incident energies
!   xseval          ! evaluated cross section
!   xsfile          ! file with crosssections
! Variables for reading MACS
!   macstalys       ! Maxwellian averaged cross section from TALYS
!   ratmacs         ! Maxwellian rate
!   Tmacs           ! Maxwellian temperature
! Variables for reading residual production cross sections
!   Erp             ! incident energy
!   Nchanrp         ! total number of channels with residual production cross sections
!   Nenrp           ! number of incident energies
!   rpfile          ! name of file with residual production cross sections
!   rptalys         ! residual production cross section from TALYS
! Variables for reading spectra
!   Eout            ! emission energy
!   Nchansp         ! total number of channels with emission spectra
!   Nensp           ! number of covariance energies
!   spfile          ! name of file with emission spectra
!   sptalys         ! emission spectrum from TALYS
! Variables for reading gamma production cross sections
!   Egam            ! incident energy
!   gamfile         ! name of file with gamma production cross sections
!   gamtalys        ! gamma cross section from TALYS
!   Nchangam        ! total number of channels with gamma production cross sections
!   Nengam          ! number of incident energies
! Variables for reading angular distributions
!   angfile         ! name of file with angular distributions
!   angle           ! angle
!   angtalys        ! angular distribution
!   Nang            ! number of angles
!   Nchanang        ! total number of channels with angular distributions
! Variables for reading Legendre coefficients
!   legfile         ! name of file with Legendre coefficients
!   legtalys        ! Legendre coefficients from TALYS
!   Nchanleg        ! total number of channels with Legendre coefficients
! Variables for reading production cross sections
!   Eprod           ! incident energy
!   Nchanprod       ! number of channels with particle production cross sections
!   Nenprod         ! number of covariance energies
!   prodfile        ! name of file with particle production cross sections
!   prodtalys       ! particle production cross section from TALYS
!   Ytalys          ! yield from TALYS
! Variables for reading integral cross sections
!   Nxseff          ! number of effective cross sections
!   xseffchan       ! channel for effective cross section
!   xseffexp        ! experimental effective cross section
!   xseffflux       ! flux for effective cross section
!   xseffrat        ! C / E of effective cross section
!   xsefftalys      ! effective cross section from TALYS
!
! *** Declaration of local data
!
  implicit none
  logical            :: lexist    ! logical
  character(len=5)   :: ranstring ! output file
  character(len=132) :: str(50)   ! string
  character(len=132) :: key
  character(len=132) :: ifile     ! input file
  character(len=132) :: ofile     ! output file
  character(len=132) :: tfile     ! TALYS input file
  integer            :: i         ! counter
  integer            :: j         ! counter
  integer            :: keyix
  integer            :: istat 
  integer            :: isamp     ! variable for 0th or random run
  integer            :: ital      ! TALYS run counter
  integer            :: k         ! counter
  integer            :: kheader   ! Legendre order
  integer            :: system    ! system call
!
! **************************** Write output ****************************
!
! The results are written to files for each random run
!
  if (italys == 0) then
    isamp = 0
  else
    isamp = 1
  endif
  if (italys ==  -1) then
    ital = 9999
  else
    ital = italys + offset
  endif
  ranstring='.    '
  write(ranstring(2:5),'(i4.4)') ital
!
! Cross sections
!
  if (flagcross) then
    do i = 1, Nchanxs
      ifile = trim(xsfile(i))
      inquire (file=ifile,exist=lexist)
      if (.not.lexist) cycle
      ofile = trim(ifile)//ranstring
      if (ofile(1:3) == '   ') cycle
      open (unit=2, file=ifile, status='old')
      kheader = 0
      k = 1
      do
        read(2, '(a)', iostat = istat) str(k)
        if (istat == -1) exit
        if (istat /= 0) call read_error(str(k), istat)
        key='entries'
        keyix=index(str(k),trim(key))
        if (keyix > 0) then
          do j=1,2
            k=k+1
            read(2, '(a)', iostat = istat) str(k)
            if (istat /= 0) call read_error(str(k), istat)
          enddo
          kheader = k
          exit
        endif
        k = k + 1
      enddo
      close (2)
      open (unit=3, file=ofile, status='replace')
      str(2) = trim(str(2))//' Run '//ranstring(2:5)
      do k = 1, kheader
        write(3, '(a)') trim(str(k))
      enddo
      do k = 1, Nen(i)
        write(3, '(2es15.6)') E(i, k), xseval(isamp, i, k)
      enddo
      close (3)
      if (flagchmagnet(i)) then
        open (unit = 3, file = 'tal.'//ofile, status = 'replace')
        do k = 1, Nen(i)
          write(3, '(2es15.6)') E(i, k), xseval(isamp, i, k)
        enddo
        close (3)
        open (unit = 3, file = 'exp.'//ofile, status = 'replace')
        do k = 1, Nxssamp(i)
          write(3, '(2es15.6)') Esamp(i, k), xssamp(i, k)
        enddo
        close (3)
      endif
    enddo
  endif
!
! Residual production cross sections
!
  if (flagresidual) then
    do i = 1, Nchanrp
      ifile = trim(rpfile(i))
      inquire (file=ifile,exist=lexist)
      if (.not.lexist) cycle
      ofile = trim(ifile)//ranstring
      open (unit=2, file=ifile, status='old')
      kheader = 0
      k = 1
      do
        read(2, '(a)', iostat = istat) str(k)
        if (istat == -1) exit
        if (istat /= 0) call read_error(str(k), istat)
        key='entries'
        keyix=index(str(k),trim(key))
        if (keyix > 0) then
          do j=1,2
            k=k+1
            read(2, '(a)', iostat = istat) str(k)
            if (istat /= 0) call read_error(str(k), istat)
          enddo
          kheader = k
          exit
        endif
        k = k + 1
      enddo
      close (2)
      open (unit=3, file=ofile, status='replace')
      str(2) = trim(str(2))//' Run '//ranstring(2:5)
      do k = 1, kheader
        write(3, '(a)') trim(str(k))
      enddo
      do k = 1, Nenrp(i)
        write(3, '(2es15.6)') Erp(i, k), rptalys(isamp, i, k)
      enddo
      close (3)
    enddo
  endif
!
! Gamma production cross sections
!
  if (flaggamma) then
    do i = 1, Nchangam
      ifile = trim(gamfile(i))
      inquire (file=ifile,exist=lexist)
      if (.not.lexist) cycle
      ofile = trim(ifile)//ranstring
      open (unit=2, file=ifile, status='old')
      kheader = 0
      k = 1
      do
        read(2, '(a)', iostat = istat) str(k)
        if (istat == -1) exit
        if (istat /= 0) call read_error(str(k), istat)
        key='entries'
        keyix=index(str(k),trim(key))
        if (keyix > 0) then
          do j=1,2
            k=k+1
            read(2, '(a)', iostat = istat) str(k)
            if (istat /= 0) call read_error(str(k), istat)
          enddo
          kheader = k
          exit
        endif
        k = k + 1
      enddo
      close (2)
      open (unit=3, file=ofile, status='replace')
      str(2) = trim(str(2))//' Run '//ranstring(2:5)
      do k = 1, kheader
        write(3, '(a)') trim(str(k))
      enddo
      do k = 1, Nengam(i)
        write(3, '(2es15.6)') Egam(i, k), gamtalys(isamp, i, k)
      enddo
      close (3)
    enddo
  endif
!
! Particle production cross sections
!
  if (flagprod) then
    do i = 1, Nchanprod
      ifile = trim(prodfile(i))
      inquire (file=ifile,exist=lexist)
      if (.not.lexist) cycle
      ofile = trim(ifile)//ranstring
      open (unit=2, file=ifile, status='old')
      kheader = 0
      k = 1
      do
        read(2, '(a)', iostat = istat) str(k)
        if (istat == -1) exit
        if (istat /= 0) call read_error(str(k), istat)
        key='entries'
        keyix=index(str(k),trim(key))
        if (keyix > 0) then
          do j=1,2
            k=k+1
            read(2, '(a)', iostat = istat) str(k)
            if (istat /= 0) call read_error(str(k), istat)
          enddo
          kheader = k
          exit
        endif
        k = k + 1
      enddo
      close (2)
      open (unit=3, file=ofile, status='replace')
      str(2) = trim(str(2))//' Run '//ranstring(2:5)
      do k = 1, kheader
        write(3, '(a)') trim(str(k))
      enddo
      do k = 1, Nenprod(i)
        write(3, '(3es15.6)') Eprod(i, k), prodtalys(isamp, i, k), Ytalys(isamp, i, k)
      enddo
      close (3)
    enddo
  endif
!
! Particle spectra
!
  if (flagspectra.and..not.flagblock) then
    do i = 1, Nchansp
      ifile = trim(spfile(i))
      inquire (file=ifile,exist=lexist)
      if (.not.lexist) cycle
      ofile = trim(ifile)//ranstring
      open (unit=2, file=ifile, status='old')
      kheader = 0
      k = 1
      do
        read(2, '(a)', iostat = istat) str(k)
        if (istat == -1) exit
        if (istat /= 0) call read_error(str(k), istat)
        key='entries'
        keyix=index(str(k),trim(key))
        if (keyix > 0) then
          do j=1,2
            k=k+1
            read(2, '(a)', iostat = istat) str(k)
            if (istat /= 0) call read_error(str(k), istat)
          enddo
          kheader = k
          exit
        endif
        k = k + 1
      enddo
      close (2)
      open (unit=3, file=ofile, status='replace')
      str(2) = trim(str(2))//' Run '//ranstring(2:5)
      do k = 1, kheader
        write(3, '(a)') trim(str(k))
      enddo
      do k = 1, Nensp(i)
        write(3, '(2es15.6)') Eout(i, k), sptalys(isamp, i, k)
      enddo
      close (3)
    enddo
  endif
!
! Angular distributions
!
  if (flagangle.and..not.flagblock) then
    do i = 1, Nchanang
      ifile = trim(angfile(i))
      inquire (file=ifile,exist=lexist)
      if (.not.lexist) cycle
      ofile = trim(ifile)//ranstring
      open (unit=2, file=ifile, status='old')
      kheader = 0
      k = 1
      do
        read(2, '(a)', iostat = istat) str(k)
        if (istat == -1) exit
        if (istat /= 0) call read_error(str(k), istat)
        key='entries'
        keyix=index(str(k),trim(key))
        if (keyix > 0) then
          do j=1,2
            k=k+1
            read(2, '(a)', iostat = istat) str(k)
            if (istat /= 0) call read_error(str(k), istat)
          enddo
          kheader = k
          exit
        endif
        k = k + 1
      enddo
      close (2)
      open (unit=3, file=ofile, status='replace')
      str(2) = trim(str(2))//' Run '//ranstring(2:5)
      do k = 1, kheader
        write(3, '(a)') trim(str(k))
      enddo
      do k = 0, Nang(i)
        write(3, '(2es15.6)') angle(i, k), angtalys(isamp, i, k)
      enddo
      close (3)
    enddo
  endif
!
! Legendre coefficients
!
  if (flagleg.and..not.flagblock) then
    do i = 1, Nchanleg
      ifile = trim(legfile(i))
      inquire (file=ifile,exist=lexist)
      if (.not.lexist) cycle
      ofile = trim(ifile)//ranstring
      open (unit=2, file=ifile, status='old')
      kheader = 0
      k = 1
      do
        read(2, '(a)', iostat = istat) str(k)
        if (istat == -1) exit
        if (istat /= 0) call read_error(str(k), istat)
        key='entries'
        keyix=index(str(k),trim(key))
        if (keyix > 0) then
          do j=1,2
            k=k+1
            read(2, '(a)', iostat = istat) str(k)
            if (istat /= 0) call read_error(str(k), istat)
          enddo
          kheader = k
          exit
        endif
        k = k + 1
      enddo
      close (2)
      open (unit=3, file=ofile, status='replace')
      str(2) = trim(str(2))//' Run '//ranstring(2:5)
      do k = 1, kheader
        write(3, '(a)') trim(str(k))
      enddo
      do k = 0, numleg
        write(4, '(i3, es15.6)') k, legtalys(isamp, i, k)
      enddo
      close (4)
    enddo
  endif
!
! Integral activation data
!
  if (flagintegral) then
    ofile = 'integral.dat'//ranstring
    open (unit = 4, file = ofile, status = 'replace')
    do k = 1, Nxseff
      write(4, '(2a15, 3es15.6)') xseffchan(k), xseffflux(k), xsefftalys(isamp, k), xseffexp(k), xseffrat(isamp, k)
    enddo
    close (4)
  endif
!
! Maxwellian averaged cross sections
!
  if (flagmacs) then
    ofile = 'macs.g'//ranstring
    open (unit = 4, file = ofile, status = 'replace')
    write(4, '(3es15.6)') thmacs, expmacs, dexpmacs
    close (4)
  endif
!
! Save input and parameter file file
!
  if (mode <= 2) then
    tfile = 'talys.inp'//ranstring
    i = system('cat talys.inp > '//tfile)
    tfile = 'parameters'//ranstring
    i = system('cat parameters.dat > '//tfile)
  endif
!
! Jump to system for possible ENDF-6 file generation
!
! tmc    : subroutine for creation of random ENDF-6 files
!
  if ((italys == 0 .or. mode == 2) .and. flagtmc) call tmc
  return
end subroutine output
! Copyright A.J. Koning 2021
