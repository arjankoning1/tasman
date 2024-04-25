subroutine ldread
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read level density GOF data from TALYS files
!
! Author    : Arjan Koning
!
! 2023-10-01: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
  use A1_error_handling_mod
!
! Definition of single and double precision variables
!   sgl ! single precision kind
! Variables for automatic parameter variation
!   ACN            ! A of compound nucleus
!   ZCN            ! Z of compound nucleus
! Variables for writing TALYS input
!   italys          ! TALYS run counter
! Variables for reading TALYS output
!   Nchanall        ! total number of channels
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist    ! logical to determine existence
  character(len=12) :: xsf       ! filename
  character(len=132) :: line
  character(len=132) :: key
  integer           :: isamp     ! variable for 0th or random run
  integer           :: istat     ! logical for file access
  integer           :: keyix
  integer           :: k
  integer           :: l1
  integer           :: Nlow
  integer           :: Ntop
  integer           :: nlev
  integer           :: Nlimit
  real(sgl)         :: Dexp
  real(sgl)         :: dDexp
  real(sgl)         :: Dglobal
  real(sgl)         :: dDglobal
  real(sgl)         :: Dth
  real(sgl)         :: Frms
  real(sgl)         :: chi2
  real(sgl)         :: Nc1
!
! Read level densities
!
  Nchanall = 1
  Nenexp(1, 1) = 1
  Nsets(1) = 2
  MT(1) = 1
  MTexp(1) = 1
  NMTexp = 1
  xsf = 'ld000000.tot'
  write(xsf(3:5), '(i3.3)') ZCN
  write(xsf(6:8), '(i3.3)') ACN
  if (italys <= 0) then
    isamp = 0
  else
    isamp = 1
  endif
  Nlimit=15
  inquire (file = xsf, exist = lexist)
  if (lexist) then
    open (unit = 2, file = xsf, status = 'old')
    read(2, '()', iostat = istat)
    if (istat /= 0) call read_error(xsf, istat)
    do
      read(2,'(a)', iostat = istat) line
      if (istat == -1) exit
      key='number of excited levels'
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(line(keyix+len_trim(key)+2:80),*, iostat = istat) nlev
        if (istat /= 0) call read_error(xsf, istat)
        if (nlev <= Nlimit) then 
          write(*, '(" TASMAN-error: Not enough levels for level density optimization =",i6)') nlev
          stop  
        endif
      endif
      key='Nlow'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Nlow
      if (istat /= 0) call read_error(xsf, istat)
      key='Ntop'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Ntop
      if (istat /= 0) call read_error(xsf, istat)
      key='experimental D0 [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Dexp
      if (istat /= 0) call read_error(xsf, istat)
      key='experimental D0 unc. [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) dDexp
      if (istat /= 0) call read_error(xsf, istat)
      key='global D0 [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Dglobal
      if (istat /= 0) call read_error(xsf, istat)
      key='global D0 unc. [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) dDglobal
      if (istat /= 0) call read_error(xsf, istat)
      key='theoretical D0 [eV]'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Dth
      if (istat /= 0) call read_error(xsf, istat)
      key='Chi-2 per level'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) chi2
      if (istat /= 0) call read_error(xsf, istat)
      key='Frms per level'
      keyix=index(line,trim(key))
      if (keyix > 0) read(line(keyix+len_trim(key)+2:80),*, iostat = istat) Frms
      if (istat /= 0) call read_error(xsf, istat)
      k = 0
      key='entries' 
      keyix=index(line,trim(key))
      if (keyix > 0) then
        read(2,'(/)',iostat = istat)
        if (istat == -1) exit
        do
          read(2, '(15x,i6,9x,es15.6)', iostat = istat) l1, Nc1 
          if (istat == -1) exit
          if (istat /= 0) call read_error(xsf, istat, eor = 'continue')
          if (l1 < Nlow) cycle
          if (l1 > Ntop) exit
          k = k + 1
          if (k > numenin) then 
            write(*, '(" TASMAN-error: Number of incident energies larger than numenin =",i6)') numenin
            stop  
          endif
          xsexp(1, 2, k) = real(l1)
          dxsexp(1, 2, k) = sqrt(real(l1))
          xsth(1, 2, k) = Nc1
          if (k == numenexp) exit
        enddo     
        exit
      endif
    enddo
    close (2)
    Nenexp(1,2) = k
    if (Dexp > 0.) then
      flagD0exp=.true.
      xsexp(1, 1, 1) = Dexp
      dxsexp(1, 1, 1) = dDexp
    else
      flagD0exp=.false.
      xsexp(1, 1, 1) = Dglobal
      dxsexp(1, 1, 1) = dDglobal
    endif
    xsth(1, 1, 1) = Dth
  endif
  return
end subroutine ldread
! Copyright A.J. Koning 2023
