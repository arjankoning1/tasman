subroutine covspectra
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for emission spectra
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
! Variables for writing TALYS input files
!   italys      ! TALYS run counter
! Variables for reading experimental data
!   Ein         ! incident energy
! Variables for cross section covariances
!   coveps      ! limit for covariance
!   flagband    ! flag to represent results as error bands
! Variables for reading spectra
!   Eout        ! emission energy
!   Nchansp     ! total number of channels with emission spectra
!   Nensp       ! number of covariance energies
!   spfile      ! name of file with emission spectra
!   sptalys     ! emission spectrum from TALYS
! Variables for spectra covariances
!   errsp       ! emission spectrum uncertainty
!   Rsp         ! covariance matrix for emission spectra
!   spav        ! average emission spectra
! Variables for covariances
!   Sw0         ! sum of weights
!   Swp0        ! sum of weights
!   Sws0        ! sum of weights
!
! *** Declaration of local data
!
  implicit none
  character(len=24) :: ofile     ! output file
  character(len=15) :: col(5)                         ! header
  character(len=15) :: un(5)
  character(len=132) :: line
  character(len=132) :: key
  character(len=16) :: reaction
  character(len=132) :: topline    ! topline
  character(len=132) :: quantity
  character(len=132), dimension(100) :: headerline   ! header line of file
  integer           :: i         ! counter
  integer           :: j         ! counter
  integer           :: l         ! counter
  integer           :: Ncol
  integer           :: keyix
  integer           :: istat
  integer           :: jheader
  real(sgl)         :: down      ! lower value
  real(sgl)         :: err       ! error
  real(sgl)         :: sp2       ! product of spectra
  real(sgl)         :: spdifi    ! difference in spectrum
  real(sgl)         :: spdifk    ! difference in spectrum
  real(sgl)         :: spi0      ! spectrum of run 0
  real(sgl)         :: spi1      ! spectrum of random run
  real(sgl)         :: spk0      ! spectrum of run 0
  real(sgl)         :: spk1      ! spectrum of random run
  real(sgl)         :: term      ! help variable
  real(sgl)         :: term0     ! help variable
  real(sgl)         :: up        ! upper value
!
! Average emission spectra and covariances
!
  do i = 1, Nchansp
    do j = 0, Nensp(i)
      spi0 = sptalys(0, i, j)
      spi1 = min(sptalys(1, i, j), 10. * spi0)
      spdifi = spi0 - spi1
      do l = 1, Nensp(i)
        spk0 = sptalys(0, i, l)
        spk1 = min(sptalys(1, i, l), 10. * spk0)
        spdifk = spk0 - spk1
        term0 = spdifi * spdifk
        if (abs(term0) > coveps .and. Sws0 > 0.) then
          sp2 = spi0 * spk0
          if (sp2 > 0.) then
            term = term0 / sp2
            Rsp(i, j, l) = (Rsp(i, j, l) * Swp0 + Sw0 * term) / Sws0
          endif
        endif
      enddo
      if (Sws0 > 0.) spav(i, j) = (spav(i, j) * Swp0 + Sw0 * spi1) / Sws0
      errsp(i, j) = sqrt(Rsp(i, j, j))
    enddo
  enddo
!
! Output
!
  un = 'mb/MeV'
  col(1)='E'
  un(1)='MeV'
  col(2)='xs'
  col(3)='xslow'
  col(4)='xsup'
  headerline=''
  quantity='emission spectrum covariance matrix'
  if (.not.flagblock) then
    do i = 1, Nchansp
      if (spfile(i)(1:1) == ' ') cycle
      ofile = trim(spfile(i))//'.ave'
      if (flagband) then
        Ncol=4
      else
        col(3)='xserr'
        Ncol=3
      endif
      j=0
      jheader=0
      open (unit = 1, file = spfile(i), status = 'old', iostat = istat)
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
          write(headerline(j)(keyix+len_trim(key)+2:80), '(i6)') Nensp(i)
          exit
        endif
      enddo
      close (1)
      open (unit = 1, file = ofile, status = 'replace')
      do j = 1, jheader
        write(1,'(a)') trim(headerline(j))
      enddo
      call write_datablock(quantity,Ncol,Nensp(i),col,un)
!     write(3, '("#", a20, " Number of runs: ", i6)') spfile(i), italys
      do j = 1, Nensp(i)
        err = spav(i, j) * errsp(i, j)
        if (flagband) then
          up = spav(i, j) * (1. + errsp(i, j))
          down = spav(i, j) / (1. + errsp(i, j))
          write(1, '(4es15.6)') Eout(i, j), spav(i, j), down, up
        else
          write(1, '(3es15.6)') Eout(i, j), spav(i, j), err
        endif
      enddo
      close (1)
    enddo
  endif
!
! Covariance matrices for emission spectra
!
  un = ''
  col(1)='E-out1'
  un(1)='MeV'
  col(2)='E-out2'
  un(2)='MeV'
  col(3)='Rspec'
  Ncol=3
  reaction='('//ptype0//',x'//ptype0//')'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
  open (unit = 1, file = 'cov_spectra.ave', status = 'replace')
  do i = 1, Nchansp
    call write_header(topline,source,user,date,oformat)
    call write_target
    call write_reaction(reaction,0.D0,0.D0,6,5)
    call write_real(2,'E-incident [MeV]',Ein(i))
    call write_covariance(reaction,35,5,0,0)
    call write_real(4,'E-incident [MeV]',Ein(i))
    call write_datablock(quantity,Ncol,(Nensp(i)+1)**2,col,un)
    do j = 0, Nensp(i)
      do l = 0, Nensp(i)
        write(1, '(3es15.6)')  Eout(i, j), Eout(i, l), Rsp(i, j, l)
      enddo
    enddo
  enddo
  close (1)
  return
end subroutine covspectra
! Copyright A.J. Koning 2021
