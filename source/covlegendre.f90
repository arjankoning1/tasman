subroutine covlegendre
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for Legendre coefficients
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
! All global variables
!   numleg      ! maximum number of Legendre coefficients
! Variables for writing TALYS input files
!   italys      ! TALYS run counter
! Variables for Legendre coefficient covariances
!   errleg      ! Legendre coefficients uncertainty
!   legav       ! average Legendre coefficients
!   Rleg        ! covariance matrix for Legendre coefficients
! Variables for reading Legendre coefficients
!   Eleg        ! energy grid for Legendre coefficients
!   leg0        ! Legendre coefficients from TALYS
!   legfile     ! name of file with Legendre coefficients
!   Nchanleg    ! total number of channels with Legendre coefficients
! Variables for cross section covariances
!   coveps      ! limit for covariance
! Variables for covariances
!   Sw0         ! sum of weights
!   Swp0        ! sum of weights
!   Sws0        ! sum of weights
!
! *** Declaration of local data
!
  implicit none
  character(len=21) :: ofile      ! output file
  character(len=15) :: col(3)                         ! header
  character(len=15) :: un(3)
  character(len=16) :: reaction
  character(len=132) :: line
  character(len=132) :: key
  character(len=132) :: quantity
  character(len=132) :: topline    ! topline
  character(len=132), dimension(100) :: headerline   ! header line of file
  integer           :: i          ! counter
  integer           :: j          ! counter
  integer           :: k          ! Legendre order
  integer           :: l          ! counter
  integer           :: Ncol
  integer           :: istat
  integer           :: keyix
  integer           :: jheader
  real(sgl)         :: err        ! error
  real(sgl)         :: leg2       ! second Legendre coefficient
  real(sgl)         :: legdifi    ! difference in Lefgendre coefficients
  real(sgl)         :: legdifk    ! difference in Lefgendre coefficients
  real(sgl)         :: legi0      ! Legendre coefficient for run 0
  real(sgl)         :: legi1      ! Legendre coefficient for random run
  real(sgl)         :: legk0      ! Legendre coefficient for run 0
  real(sgl)         :: legk1      ! Legendre coefficient for random run
  real(sgl)         :: term       ! help variable
  real(sgl)         :: term0      ! help variable
!
! Average Legendre coefficients and covariances
!
  do i = 1, Nchanleg
    do j = 0, numleg
      legi0 = leg0(0, i, j)
      legi1 = leg0(1, i, j)
      legdifi = legi0 - legi1
      do k = 1, Nchanleg
        do l = 0, numleg
          legk0 = leg0(0, k, l)
          legk1 = leg0(1, k, l)
          legdifk = legk0 - legk1
          term0 = legdifi * legdifk
          if (abs(term0) > coveps .and. Sws0 > 0.) then
            leg2 = legi0 * legk0
            if (leg2 > 0.) then
              term = term0 / leg2
              Rleg(i, j, k, l) = (Rleg(i, j, k, l) * Swp0 + Sw0 * term) / Sws0
            endif
          endif
        enddo
      enddo
      if (Sws0 > 0.) legav(i, j) = (legav(i, j) * Swp0 + Sw0 * legi1) / Sws0
      errleg(i, j) = sqrt(Rleg(i, j, i, j))
    enddo
  enddo
!
! Average Legendre coefficients
!
  Ncol=3
  quantity='legendre coefficient covariance matrix'
  if (.not.flagblock) then
    un = ''
    col(1)='l'
    col(2)='Legendre coeff.'
    col(3)='uncertainty'
    headerline=''
    do i = 1, Nchanleg
      if (legfile(i)(1:1) == ' ') cycle
      ofile = trim(legfile(i)) // '.ave'
      j=0
      jheader=0
      open (unit = 1, file = legfile(i), status = 'old', iostat = istat)
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
          write(headerline(j)(keyix+len_trim(key)+2:80), '(i6)') numleg
          exit
        endif
      enddo
      close (1)
      open (unit = 1, file = ofile, status = 'replace')
!     write(1, '("#", a21, " Number of runs: ", i6)') ofile, italys
      do j = 1, jheader
        write(1,'(a)') trim(headerline(j))
      enddo
      call write_datablock(quantity,Ncol,numleg,col,un)
      do j = 0, numleg
        err = legav(i, j) * errleg(i, j)
        write(3, '(3x,i6,6x, 2es15.6)') j, legav(i, j), err
      enddo
      close (3)
    enddo
  endif
!
! Covariance matrices for Legendre coefficients
!
  un=''
  Ncol=3
  col(1)='l1'
  col(2)='l2'
  col(3)='Rleg'
  open (unit = 1, file = 'cov_legendre.ave', status = 'replace')
  do i = 1, Nchanleg
    do k = 1, Nchanleg
      reaction='('//ptype0//',el)'
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
      call write_header(topline,source,user,date,oformat)
      call write_target
      call write_reaction(reaction,0.D0,0.D0,4,2)
      call write_real(2,'E-incident [MeV]',Eleg(i))
      call write_covariance(reaction,4,2,4,2)
      call write_real(4,'E-incident [MeV]',Eleg(k))
      call write_datablock(quantity,Ncol,(numleg+1)**2,col,un)
      do j = 0, numleg
        do l = 0, numleg
          write(1, '(2(3x,i6,6x), es15.6)') j, l, Rleg(i, j, k, l)
        enddo
      enddo
    enddo
  enddo
  close (1)
  return
end subroutine covlegendre
! Copyright A.J. Koning 2021
