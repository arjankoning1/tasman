subroutine covresidual
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for residual production cross sections
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
! Variables for reading cross sections
!   E            ! incident energy
! Variables for cross section covariances
!   coveps       ! limit for covariance
!   Ecovindex    ! index for covariance energy on main energy grid
!   flagband     ! flag to represent results as error bands
!   Nencov       ! number of covariance energies
! Variables for covariances
!   Sw0          ! sum of weights
!   Swp0         ! sum of weights
!   Sws0         ! sum of weights
! Variables for reading residual production cross sections
!   Arp          ! mass number of residual product
!   Erp          ! incident energy
!   Lrp          ! isomeric level of residual product
!   Nchanrp      ! total number of channels with residual production cross sections
!   Nenrp        ! number of incident energies
!   rpfile       ! name of file with residual production cross sections
!   rptalys      ! residual production cross section from TALYS
!   Zrp          ! charge number of residual product
! Variables for residual production cross section variances
!   errrp        ! cross section uncertainty
!   Rrp          ! covariance matrix for residual production cross sections
!   rpav         ! average residual production cross section
!   RrpD         ! diagonal of covariance matrix for cross sections
!
! *** Declaration of local data
!
  implicit none
  character(len=20) :: ofile     ! output file
  character(len=1)  :: isoch
  character(len=3)  :: massstring
  character(len=6)  :: finalnuclide
  character(len=15) :: col(5)                         ! header
  character(len=15) :: un(5)
  character(len=16) :: reaction
  character(len=132) :: line
  character(len=132), dimension(100) :: headerline   ! header line of file
  character(len=132) :: key
  character(len=132) :: quantity
  character(len=132) :: topline    ! topline
  integer           :: i         ! counter
  integer           :: j         ! counter
  integer           :: jheader
  integer           :: istat
  integer           :: keyix
  integer           :: jj        ! counter
  integer           :: l         ! counter
  integer           :: ll        ! angular momentum
  integer           :: Ncol
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
! Average residual production cross sections and covariances
!
  col = ''
  un = ''
  quantity='cross section covariance matrix'
  do i = 1, Nchanrp
    do j = 1, Nencov
      jj = Ecovindex(j)
      xsi0 = rptalys(0, i, jj)
      xsi1 = min(rptalys(1, i, jj), 10. * xsi0)
      xsdifi = xsi0 - xsi1
      do l = 1, Nencov
        ll = Ecovindex(l)
        xsk0 = rptalys(0, i, ll)
        xsk1 = min(rptalys(1, i, ll), 10. * xsk0)
        xsdifk = xsk0 - xsk1
        term0 = xsdifi * xsdifk
        if (abs(term0) > coveps .and. Sws0 > 0.) then
          xs2 = xsi0 * xsk0
          if (xs2 > 0.) then
            term = term0 / xs2
            Rrp(i, j, l) = (Rrp(i, j, l) * Swp0 + Sw0 * term) / Sws0
          endif
        endif
      enddo
    enddo
    do j = 1, Nenrp(i)
      xsi0 = rptalys(0, i, j)
      xsi1 = min(rptalys(1, i, j), 10. * xsi0)
      xsdifi = xsi0 - xsi1
      term0 = xsdifi * xsdifi
      if (abs(term0) > coveps .and. Sws0 > 0.) then
        xs2 = xsi0 * xsi0
        if (xs2 > 0.) then
          term = term0 / xs2
          RrpD(i, j) = (RrpD(i, j) * Swp0 + Sw0 * term) / Sws0
          errrp(i, j) = sqrt(RrpD(i, j))
        endif
      endif
      if (Sws0 > 0.) rpav(i, j) = (rpav(i, j) * Swp0 + Sw0 * xsi1) / Sws0
    enddo
  enddo
!
! ******** Output of covariance matrix and average cross sections ******
!
! Covariance matrices for residual production cross sections
!
  un=''
  col(1)='E_a'
  un(1)='MeV'
  col(2)='E_b'
  un(2)='MeV'
  col(3)='Rcov'
  Ncol=3
!
! Intra-channel correlations
!
  open (unit = 1, file = 'cov_residual.ave', status = 'replace')
  do i = 1, Nchanrp
    massstring='   '
    write(massstring,'(i3)') Arp(i)
    if (Lrp(i) == 0) isoch = ' '
    if (Lrp(i) == 1) isoch = 'm'
    if (Lrp(i) > 1) isoch = 'n'
    finalnuclide=trim(nuc(Zrp(i)))//trim(adjustl(massstring))//isoch
    reaction='('//ptype0//',x)'
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,0.D0,0.D0,6,5)
    call write_residual(Zrp(i),Arp(i),finalnuclide)
    if (Lrp(i) >= 0) call write_level(2,Lrp(i),-1,0.,-1.,0,0.)
    call write_covariance(reaction,6,5,6,5)
    call write_datablock(quantity,Ncol,Nencov*Nencov,col,un)
    do j = 1, Nencov
      jj = Ecovindex(j)
      do l = 1, Nencov
        ll = Ecovindex(l)
!       write(4, '(2(4i4, es15.6), es15.6)') Zrp(i), Arp(i), Lrp(i), jj, e(i, jj), Zrp(i), Arp(i), Lrp(i), ll, e(i, ll), &
!&        Rrp(i, j, l)
        write(1, '(3es15.6)') e(i, jj),  e(i, ll), Rrp(i, j, l)
      enddo
    enddo
  enddo
  close (1)
!
! Average cross sections
!
  un = 'mb'
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  col(3)='xslow'
  col(4)='xsup'
  Ncol=4
  headerline = ''
  do i = 1, Nchanrp
    ofile = rpfile(i)
    if (ofile(1:1) == ' ') cycle
    if (ofile(10:10) == 'L') then
      write(ofile(13:16), '(".ave")')
    else
      write(ofile(10:12), '("ave")')
    endif
    if (flagband) then
      Ncol=4
    else
      col(3)='xserr'
      Ncol=3
    endif
    j=0
    jheader=0
    open (unit = 1, file = rpfile(i), status = 'old', iostat=istat)
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
        write(headerline(j)(keyix+len_trim(key)+2:80), '(i6)') Nenrp(i)
        exit
      endif
    enddo
    close (1)
!   write(1, '("#", a20, " Number of runs: ", i6)') rpfile(i), italys
    open (unit = 1, file = ofile, status = 'replace')
    do j = 1, jheader
      write(1,'(a)') trim(headerline(j))
    enddo
    call write_datablock(quantity,Ncol,Nen(i),col,un)
    do j = 1, Nenrp(i)
      err = rpav(i, j) * errrp(i, j)
      if (err /= 0.) err = min(err, rpav(i, j) - 1.e-13)
      if (flagband) then
        write(1, '(4es15.6)') Erp(i, j), rpav(i, j), rpav(i, j)-err, rpav(i, j) + err
      else
        write(1, '(3es15.6)') Erp(i, j), rpav(i, j), err
      endif
    enddo
    close (1)
  enddo
  return
end subroutine covresidual
! Copyright A.J. Koning 2021
