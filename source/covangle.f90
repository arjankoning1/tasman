subroutine covangle
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for angular distributions
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
!   sgl         ! single precision kind
!  Variables for writing TALYS input files
!   italys      ! TALYS run counter
! Variables for cross section covariances
!   flagband    ! flag to represent results as error bands
!   coveps      ! limit for covariance
! Variables for reading angular distributions
!   angfile     ! name of file with angular distributions
!   angle       ! angle
!   angtalys    ! angular distribution
!   Nang        ! number of angles
!   Nchanang    ! total number of channels with angular distributions
! Variables for angular distribution covariances
!   angav       ! average angular distributions
!   errang      ! angular distributions uncertainty
!   Rang        ! covariance matrix for angular distributions
! Variables for covariances
!   Sw0         ! sum of weights
!   Swp0        ! sum of weights
!   Sws0        ! sum of weights
!
! *** Declaration of local data
!
  implicit none
  character(len=21) :: ofile      ! output file
  character(len=15) :: col(5)                         ! header
  character(len=15) :: un(5)
  character(len=132) :: line
  character(len=132) :: key
  character(len=16) :: reaction
  character(len=132) :: quantity
  character(len=132) :: topline    ! topline
  character(len=132), dimension(100) :: headerline   ! header line of file
  integer           :: i          ! counter
  integer           :: j          ! counter
  integer           :: l          ! counter
  integer           :: Ncol
  integer           :: istat
  integer           :: keyix
  integer           :: jheader
  real(sgl)         :: ang2       ! product of angular cross section
  real(sgl)         :: angdifi    ! difference in angular cross section
  real(sgl)         :: angdifk    ! difference in angular cross section
  real(sgl)         :: angi0      ! angular cross section of run 0
  real(sgl)         :: angi1      ! angular cross section of random run
  real(sgl)         :: angk0      ! angular cross section of run 0
  real(sgl)         :: angk1      ! angular cross section of random run
  real(sgl)         :: down       ! lower value
  real(sgl)         :: err        ! error
  real(sgl)         :: term       ! help variable
  real(sgl)         :: term0      ! help variable
  real(sgl)         :: up         ! upper value
  real(sgl)         :: xslimit    ! xs boundary for inclusion
!
! Average angular distributions
!
  do i = 1, Nchanang
    do j = 0, Nang(i)
      angi0 = angtalys(0, i, j)
      angi1 = min(angtalys(1, i, j), angi0 * xslimit)
      angi1 = max(angi1, angi0 / xslimit)
      if (Sws0 > 0.) angav(i, j) = (angav(i, j) * Swp0 + Sw0 * angi1) / Sws0
    enddo
  enddo
!
! Covariances
!
  xslimit = 10.
  do i = 1, Nchanang
    do j = 0, Nang(i)
      angi0 = angtalys(0, i, j)
      angi1 = min(angtalys(1, i, j), angi0 * xslimit)
      angi1 = max(angi1, angi0 / xslimit)
      angi0 = angav(i, j)
      angdifi = angi0 - angi1
      do l = 0, Nang(i)
        angk0 = angtalys(0, i, l)
        angk1 = min(angtalys(1, i, l), angi0 * xslimit)
        angk1 = max(angi1, angi0 / xslimit)
        angk0 = angav(i, l)
        angdifk = angk0 - angk1
        term0 = angdifi * angdifk
        if (abs(term0) > coveps .and. Sws0 > 0.) then
          ang2 = angi0 * angk0
          if (ang2 > 0.) then
            term = term0 / ang2
            Rang(i, j, l) = (Rang(i, j, l) * Swp0 + Sw0 * term) / Sws0
          endif
        endif
      enddo
      errang(i, j) = sqrt(Rang(i, j, j))
    enddo
  enddo
!
! Output
!
! Average angular distributions
!
  reaction='('//ptype0//',el)'
  quantity='angular distribution covariance matrix'
  if (.not.flagblock) then
    un = 'mb/sr'
    col(1)='Angle'
    un(1)='degrees'
    col(2)='xs'
    col(3)='xslow'
    col(4)='xsup'
    headerline=''
    do i = 1, Nchanang
      if (angfile(i)(1:1) == ' ') cycle
      ofile = trim(angfile(i)) // '.ave'
      if (flagband) then
        Ncol=4
      else
        col(3)='xserr'
        Ncol=3
      endif
      j=0
      jheader=0
      open (unit = 1, file = angfile(i), status = 'old', iostat = istat)
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
          write(headerline(j)(keyix+len_trim(key)+2:80), '(i6)') Nang(i)+1
          exit
        endif
      enddo
      close (1)
      open (unit = 1, file = ofile, status = 'replace')
!     write(1, '("#", a21, " Number of runs: ", i6)') ofile, italys
      do j = 1, jheader
        write(1,'(a)') trim(headerline(j))
      enddo
      call write_datablock(quantity,Ncol,Nang(i)+1,col,un)
      do j = 0, Nang(i)
        err = angav(i, j) * errang(i, j)
        if (flagband) then
          up = angav(i, j) * (1. + errang(i, j))
          down = angav(i, j) / (1. + errang(i, j))
          write(1, '(4es15.6)') angle(i, j), angav(i, j), down, up
        else
          write(1, '(3es15.6)') angle(i, j), angav(i, j), err
        endif
      enddo
      close (1)
    enddo
  endif
!
! Covariance matrices for angular distributions
!
  un = ''
  col(1)='Angle1'
  un(1)='degrees'
  col(2)='Angle2'
  un(2)='degrees'
  col(3)='Rangle'
  Ncol=3
  open (unit = 1, file = 'cov_angle.ave', status = 'replace')
  do i = 1, Nchanang
    reaction='('//ptype0//',el)'
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,0.D0,0.D0,4,2)
    call write_real(2,'E-incident [MeV]',Eang(i))
    call write_covariance(reaction,4,2,4,2)
    call write_real(4,'E-incident [MeV]',Eang(i))
    call write_datablock(quantity,Ncol,(Nang(i)+1)**2,col,un)
    do j = 0, Nang(i)
      do l = 0, Nang(i)
        write(1, '(2(3x,i6,6x),es15.6)') j, l, Rang(i, j, l)
      enddo
    enddo
  enddo
  close (1)
  return
end subroutine covangle
! Copyright A.J. Koning 2021
