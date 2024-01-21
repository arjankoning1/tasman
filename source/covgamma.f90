subroutine covgamma
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for gamma cross sections
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
!
! Definition of single and double precision variables
!   sgl          ! single precision kind
! Variables for writing TALYS input files
!   italys       ! TALYS run counter
! Variables for cross section covariances
!   Ecovindex    ! index for covariance energy on main energy grid
!   flagband     ! flag to represent results as error bands
!   coveps       ! limit for covariance
!   Nencov       ! number of covariance energies
! Variables for reading gamma production cross sections
!   Egam         ! incident energy
!   gamfile      ! name of file with gamma production cross sections
!   gamtalys     ! gamma cross section from TALYS
!   Nchangam     ! total number of channels with gamma production cross sections
!   Nengam       ! number of incident energies
! Variables for gamma production cross section covariances
!   errgam       ! cross section uncertainty
!   gamav        ! average gamma production cross section
!   Rgam         ! covariance matrix for gamma production cross sections
!   RgamD        ! diagonal of covariance matrix for cross sections
! Variables for covariances
!   Sw0          ! sum of weights
!   Swp0         ! sum of weights
!   Sws0         ! sum of weights
!
! *** Declaration of local data
!
  implicit none
  character(len=20) :: ofile     ! output file
  character(len=15) :: col(5)                         ! header
  character(len=15) :: un(5)
  character(len=132) :: line
  character(len=132) :: key
  character(len=132) :: quantity
  character(len=132), dimension(100) :: headerline   ! header line of file
  integer           :: i         ! counter
  integer           :: j         ! counter
  integer           :: jj        ! counter
  integer           :: l         ! counter
  integer           :: Ncol
  integer           :: keyix
  integer           :: istat
  integer           :: jheader
  integer           :: ll        ! angular momentum
  real(sgl)         :: err       ! error
  real(sgl)         :: term      ! help variable
  real(sgl)         :: term0     ! help variable
  real(sgl)         :: xs2       ! help variable
  real(sgl)         :: xsdifi    ! difference in cross section
  real(sgl)         :: xsdifk    ! difference in cross section
  real(sgl)         :: xsi0      ! cross section of run 0
  real(sgl)         :: xsi1      ! cross section of random run
  real(sgl)         :: xsk0      ! cross section of run 0
  real(sgl)         :: xsk1      ! cross section of random run
!
! Average gamma production cross sections and covariances
!
  do i = 1, Nchangam
    do j = 1, Nencov
      jj = Ecovindex(j)
      xsi0 = gamtalys(0, i, jj)
      xsi1 = min(gamtalys(1, i, jj), 10. * xsi0)
      xsdifi = xsi0 - xsi1
      do l = 1, Nencov
        ll = Ecovindex(l)
        xsk0 = gamtalys(0, i, ll)
        xsk1 = min(gamtalys(1, i, ll), 10. * xsk0)
        xsdifk = xsk0 - xsk1
        term0 = xsdifi * xsdifk
        if (abs(term0) > coveps) then
          xs2 = xsi0 * xsk0
          if (xs2 > 0.) then
            term = term0 / xs2
            Rgam(i, j, l) = (Rgam(i, j, l) * Swp0 + Sw0 * term) / Sws0
          endif
        endif
      enddo
    enddo
    do j = 1, Nengam(i)
      xsi0 = gamtalys(0, i, j)
      xsi1 = min(gamtalys(1, i, j), 10. * xsi0)
      xsdifi = xsi0 - xsi1
      term0 = xsdifi * xsdifi
      if (abs(term0) > coveps .and. Sws0 > 0.) then
        xs2 = xsi0 * xsi0
        if (xs2 > 0.) then
          term = term0 / xs2
          RgamD(i, j) = (RgamD(i, j) * Swp0 + Sw0 * term) / Sws0
          errgam(i, j) = sqrt(RgamD(i, j))
        endif
      endif
      if (Sws0 > 0.) gamav(i, j) = (gamav(i, j) * Swp0 + Sw0 * xsi1) / Sws0
    enddo
  enddo
!
! Output
!
  un = 'mb'
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  col(3)='xslow'
  col(4)='xsup'
  headerline=''
  quantity='cross section covariance matrix'
  do i = 1, Nchangam
    ofile = gamfile(i)
    if (ofile(1:1) == ' ') cycle
    write(ofile(16:19), '(".ave")')
    if (flagband) then
      Ncol=4
    else
      col(3)='xserr'
      Ncol=3
    endif
    j=0
    jheader=0
    open (unit = 1, file = gamfile(i), status = 'old', iostat = istat)
    if (istat /= 0) cycle
    do
      j=j+1
      read(1,'(a)',iostat = istat) line
      if (istat == -1) exit
      if (line(1:1) /= '#') exit
      headerline(j)=line
      key='source:'
      keyix=index(line,trim(key))
      if (keyix > 0) write(headerline(j)(keyix+len_trim(key)+1:80),'("TASMAN")')
      key='datablock'
      keyix=index(headerline(j),trim(key))
      if (keyix > 0) jheader = j-1
      key='entries'
      keyix=index(headerline(j),trim(key))
      if (keyix > 0) then
        write(headerline(j)(keyix+len_trim(key)+2:80), '(i6)') Nengam(i)
        exit
      endif
    enddo
    close (1)
    open (unit = 1, file = ofile, status = 'replace')
    do j = 1, jheader
      write(1,'(a)') trim(headerline(j))
    enddo
    call write_datablock(quantity,Ncol,Nengam(i),col,un)
!   write(3, '("#", a20, " Number of runs: ", i6)') gamfile(i), italys
    do j = 1, Nengam(i)
      err = gamav(i, j) * errgam(i, j)
      if (err /= 0.) err = min(err, gamav(i, j) - 1.e-13)
      if (flagband) then
        write(1, '(4es15.6)') Egam(i, j), gamav(i, j), gamav(i, j) - err, gamav(i, j) + err
      else
        write(1, '(3es15.6)') Egam(i, j), gamav(i, j), err
      endif
    enddo
    close (1)
  enddo
  return
end subroutine covgamma
! Copyright A.J. Koning 2021
