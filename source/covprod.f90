subroutine covprod
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for particle production cross sections
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
!   coveps       ! limit for covariance
!   Ecovindex    ! index for covariance energy on main energy grid
!   flagband     ! flag to represent results as error bands
!   Nencov       ! number of covariance energies
! Variables for reading production cross sections
!   Eprod        ! incident energy
!   Nchanprod    ! number of channels with particle production cross sections
!   Nenprod      ! number of covariance energies
!   prodfile     ! name of file with particle production cross sections
!   prodtalys    ! particle production cross section from TALYS
! Variables for production cross section covariances
!   errprod      ! cross section uncertainty
!   prodav       ! average particle production cross section
!   Rprod        ! covariance matrix for particle production cross sections
!   RprodD       ! diagonal of covariance matrix for cross sections
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
  integer           :: ll        ! angular momentum
  integer           :: Ncol
  integer           :: keyix
  integer           :: istat
  integer           :: jheader
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
! Average particle production cross sections and covariances
!
!   sections
!
  do i = 1, Nchanprod
    do j = 1, Nencov
      jj = Ecovindex(j)
      xsi0 = prodtalys(0, i, jj)
      xsi1 = min(prodtalys(1, i, jj), 10. * xsi0)
      xsdifi = xsi0 - xsi1
      do l = 1, Nencov
        ll = Ecovindex(l)
        xsk0 = prodtalys(0, i, ll)
        xsk1 = min(prodtalys(1, i, ll), 10. * xsk0)
        xsdifk = xsk0 - xsk1
        term0 = xsdifi * xsdifk
        if (abs(term0) > coveps .and. Sws0 > 0.) then
          xs2 = xsi0 * xsk0
          if (xs2 > 0.) then
            term = term0 / xs2
            Rprod(i, j, l) = (Rprod(i, j, l) * Swp0 + Sw0 * term) / Sws0
          endif
        endif
      enddo
    enddo
    do j = 1, Nenprod(i)
      xsi0 = prodtalys(0, i, j)
      xsi1 = min(prodtalys(1, i, j), 10. * xsi0)
      xsdifi = xsi0 - xsi1
      if (Sws0 > 0.) prodav(i, j) = (prodav(i, j) * Swp0 + Sw0 * xsi1) / Sws0
      term0 = xsdifi * xsdifi
      if (abs(term0) > coveps .and. Sws0 > 0.) then
        xs2 = xsi0 * xsi0
        if (xs2 > 0.) then
          term = term0 / xs2
          RprodD(i, j) = (RprodD(i, j) * Swp0 + Sw0 * term) / Sws0
          errprod(i, j) = sqrt(RprodD(i, j))
        endif
      endif
    enddo
  enddo
!
! Output
!
  quantity='cross section covariance matrix'
  un = 'mb'
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  col(3)='xslow'
  col(4)='xsup'
  headerline=''
  do i = 1, Nchanprod
    ofile = prodfile(i)
    if (ofile(1:1) == ' ') cycle
    write(ofile(6:9), '(".ave")')
    if (flagband) then
      Ncol=4
    else
      col(3)='xserr'
      Ncol=3
    endif
    j=0
    jheader=0
    open (unit = 1, file = prodfile(i), status = 'old', iostat = istat)
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
        write(headerline(j)(keyix+len_trim(key)+2:80), '(i6)') Nenprod(i)
        exit
      endif
    enddo
    close (1)
    open (unit = 1, file = ofile, status = 'replace')
    do j = 1, jheader
      write(1,'(a)') trim(headerline(j))
    enddo
    call write_datablock(quantity,Ncol,Nenprod(i),col,un)
    do j = 1, Nenprod(i)
      err = prodav(i, j) * errprod(i, j)
      if (err /= 0.) err = min(err, prodav(i, j) - 1.e-13)
      if (flagband) then
        write(1, '(4es15.6)') Eprod(i, j), prodav(i, j), prodav(i, j) - err, prodav(i, j) + err
      else
        write(1, '(3es15.6)') Eprod(i, j), prodav(i, j), err
      endif
    enddo
    close (1)
  enddo
  return
end subroutine covprod
! Copyright A.J. Koning 2021
