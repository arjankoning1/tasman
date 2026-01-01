subroutine covcross
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for cross sections
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
!   numchanxs       ! maximum number of channels with cross sections
!   numenin         ! maximum number of incident energies
!   numpar          ! maximum number of parameters
! Variables for writing TALYS input files
!   italys          ! TALYS run counter
! Variables for processing input
!   Atarget         ! mass number of target nucleus
!   Ztarget         ! charge number of target nucleus
! Variables for GOF function
!   Eindex          ! index for experimental data
! Variables for reading experimental data
!   Eexp            ! incident energy
!   Nenexp          ! number of experimental data points
!   Nsets           ! number of experimental data sets
!   xsexp           ! experimental cross section
! Variables for parameter variation
!   Npar            ! number of parameters
!   parsave         ! parameter value
!   partalys        ! parameter value
! Variables for weights
!   xsopt           ! optimal cross section
! Variables for covariances
!   Sweight         ! weight for TALYS run
!   Sweightprev     ! weight for previous TALYS run
!   Sweightsum      ! sum of weights for TALYS run
! Variables for sensitivity
!   flagreadsens    ! flag to read sensitivities from tables
!   flagsens        ! flag to use sensitivity in parameter weighting
! Variables for reading cross sections
!   E               ! incident energy
!   MT              ! MT number
!   MTiso           ! isomer of MT number
!   Nchanxs         ! total number of channels with cross sections
!   Nen             ! number of incident energies
!   Sindex          ! index for cross section
!   xseval          ! evaluated cross section
!   xsfile          ! file with crosssections
! Variables for cross section covariances
!   average         ! cross sections 1: average, 2: first run 3: optimum
!   Ecovindex       ! index for covariance energy on main energy grid
!   Ccov            ! correlation matrix for cross sections
!   Cmt             ! intra - channel correlation matrix for cross sections
!   coveps          ! limit for covariance
!   errmt           ! cross section uncertainty
!   errmtC          ! cross section uncertainty (for cov. energy grid)
!   flagband        ! flag to represent results as error bands
!   MTcov           ! MT number with covariance data
!   MTisocov        ! isomer of MT number with covariance data
!   Nchancovint     ! number of channels with covariance data
!   Nencov          ! number of covariance energies
!   Rcov            ! relative covariance matrix for cross sections
!   Rlimit          ! limit for covariance calculation
!   Rmt             ! intra - channel rel. cov. matrix for cross sections
!   RmtD            ! diagonal of covariance matrix for cross sections
!   S               ! sensitivity matrix
!   Sdenom          ! denominator of sensitivity matrix
!   Senum           ! enumerator of sensitivity matrix
!   Vcov            ! covariance matrix for cross sections
!   Vmt             ! intra - channel covariance matrix for cross sections
!   xsav            ! average cross section
!   xsavC           ! average cross section (for covariance energy grid)
!
! *** Declaration of local data
!
  implicit none
  character(len=20) :: ofile                          ! output file
  character(len=15) :: col(5)                         ! header
  character(len=15) :: un(5)
  character(len=16) :: reaction
  character(len=16) :: reaction_cov
  character(len=132) :: line
  character(len=132) :: key
  character(len=132) :: quantity
  character(len=132) :: topline    ! topline
  character(len=132), dimension(100) :: headerline   ! header line of file
  integer           :: i                              ! counter
  integer           :: istat
  integer           :: ichan                          ! counter for channels
  integer           :: ifile                          ! file
  integer           :: j                              ! counter
  integer           :: jj                             ! counter
  integer           :: k                              ! Legendre order
  integer           :: keyix
  integer           :: Ncol
  integer           :: indent
  integer           :: id2
  integer           :: id4
  integer           :: kchan                          ! counter
  integer           :: l                              ! counter
  integer           :: ll                             ! angular momentum
  integer           :: m                              ! counter
  integer           :: n                              ! counter
  integer           :: NE                             ! number of incident energies
  integer           :: Ninside                        ! number of points inside error band
  integer           :: NinsideE(numenin)              ! number of points inside error band as function of energy
  integer           :: Ninsidetot                     ! total number of points inside error band
  integer           :: NinsidexsE(numenin)            ! number of points inside error band as function of energy
  integer           :: Noutside                       ! number of points outside error band
  integer           :: NoutsideE(numenin)             ! number of points outside error band as function of energy
  integer           :: Noutsidetot                    ! total number of points outside error band
  integer           :: NoutsidexsE(numenin)           ! number of points outside error band as function of energy
  integer           :: jheader
  integer           :: Np                             ! number of points
  integer           :: Npoints                        ! number of points
  integer           :: NpointsE(numenin)              ! number of points as function of energy
  integer           :: Npointstot                     ! total number of points
  integer           :: NpointsxsE(numenin)            ! number of points per channel as function of energy
  integer           :: Nsetstot                       ! total number of sets
  real(sgl)         :: e1                             ! first energy
  real(sgl)         :: e2                             ! second energy
  real(sgl)         :: ee(0:numenin)                  ! incident energy
  real(sgl)         :: Eex                            ! excitation energy
  real(sgl)         :: err                            ! error
  real(sgl)         :: fac                            ! factor
  real(sgl)         :: fracinside                     ! fraction of points inside error band
  real(sgl)         :: fracinsideE(numenin)           ! fraction of points inside error band as function of energy
  real(sgl)         :: fracinsidetot                  ! fraction of total points inside error band
  real(sgl)         :: fracinsidexsE(numenin)         ! fraction of points inside error band as function of  energy
  real(sgl)         :: pardifj                        ! difference between parameters
  real(sgl)         :: parj0                          ! parameter of run 0
  real(sgl)         :: parjk                          ! parameter of random run
  real(sgl)         :: Sw                             ! weight for TALYS run
  real(sgl)         :: Swp                            ! weight for TALYS run
  real(sgl)         :: Sws                            ! sum of weights for TALYS run
  real(sgl)         :: term                           ! help variable
  real(sgl)         :: term0                          ! help variable
  real(sgl)         :: V2                             ! square of covariance matrix element
  real(sgl)         :: xs2                            ! help variable
  real(sgl)         :: xsdifi                         ! difference in cross section
  real(sgl)         :: xsdifk                         ! difference in cross section
  real(sgl)         :: xsex                           ! help variable
  real(sgl)         :: Ssign                          ! help variable
  real(sgl)         :: xsA                            ! average cross section
  real(sgl)         :: xsi0                           ! cross section of run 0
  real(sgl)         :: xsi1                           ! cross section of random run
  real(sgl)         :: xsk0                           ! cross section of run 0
  real(sgl)         :: xsk1                           ! cross section of random run
  real(sgl)         :: xslow(numchanxs, 0:numenin)    ! lower cross section band
  real(sgl)         :: xslowint                       ! interpolated lower cross section band
  real(sgl)         :: xsupp(numchanxs, 0:numenin)    ! upper cross section band
  real(sgl)         :: xsuppint                       ! interpolated upper cross section band
  real(sgl)         :: xslimit                        ! xs boundary for inclusion
!
! ******** Create covariance matrix and average cross sections *********
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  col = ''
  un = ''
  quantity='cross section covariance matrix'
  xslimit = 10.
!
! Average cross sections
!
  do i = 1, Nchanxs
    do j = 1, Nen(i)
      Swp = Sweightprev(j)
      Sw = Sweight(j)
      Sws = Sweightsum(j)
      xsi0 = xseval(0, i, j)
      xsi1 = min(xseval(1, i, j), xsi0 * xslimit)
      xsi1 = max(xsi1, xsi0 / xslimit)
      if (Sws > 0.) xsav(i, j) = (xsav(i, j) * Swp + Sw * xsi1) / Sws
    enddo
  enddo
  do i = 1, Nchanxs
    do j = 1, Nencov
      jj = Ecovindex(j)
      Swp = Sweightprev(jj)
      Sw = Sweight(jj)
      Sws = Sweightsum(jj)
      xsi0 = xseval(0, i, jj)
      xsi1 = min(xseval(1, i, jj), xsi0 * xslimit)
      xsi1 = max(xsi1, xsi0 / xslimit)
      if (Sws > 0.) xsavC(i, j) = (xsavC(i, j) * Swp + Sw * xsi1) / Sws
    enddo
  enddo
!
! Cross section covariance matrix R (inter-channels)
!
  do ichan = 1, Nchancovint
    do i = 1, Nchanxs
      if (MT(i) == MTcov(ichan) .and. MTiso(i) == MTisocov(ichan)) then
        do j = 1, Nencov
          jj = Ecovindex(j)
          Swp = Sweightprev(jj)
          Sw = Sweight(jj)
          Sws = Sweightsum(jj)
          xsi0 = xseval(0, i, jj)
          xsi1 = min(xseval(1, i, jj), xsi0 * xslimit)
          xsi1 = max(xsi1, xsi0 / xslimit)
          xsi0 = xsavC(i, j)
          xsdifi = xsi0 - xsi1
          do kchan = 1, Nchancovint
            do k = 1, Nchanxs
              if (MT(k) == MTcov(kchan) .and. MTiso(k) == MTisocov(kchan)) then
                do l = 1, Nencov
                  ll = Ecovindex(l)
                  xsk0 = xseval(0, k, ll)
                  xsk1 = min(xseval(1, k, ll), xsk0 * xslimit)
                  xsk1 = max(xsk1, xsk0 / xslimit)
                  xsk0 = xsavC(k, l)
                  xsdifk = xsk0 - xsk1
                  term0 = xsdifi * xsdifk
                  if (abs(term0) > coveps .and. Sws > 0.) then
                    Vcov(ichan, j, kchan, l) = (Vcov(ichan, j, kchan, l) * Swp + Sw * term0) / Sws
                    xs2 = xsi0 * xsk0
                    if (xs2 > 0.) then
                      term = term0 / xs2
                      Rcov(ichan, j, kchan, l) = (Rcov(ichan, j, kchan, l) * Swp + Sw * term) / Sws
                    endif
                  endif
                enddo
              endif
            enddo
          enddo
        enddo
      endif
    enddo
  enddo
!
! Correlation matrix
!
  do ichan = 1, Nchancovint
    do i = 1, Nchanxs
      if (MT(i) == MTcov(ichan) .and. MTiso(i) == MTisocov(ichan)) then
        do j = 1, Nencov
          do kchan = 1, Nchancovint
            do k = 1, Nchanxs
              if (MT(k) == MTcov(kchan) .and. MTiso(k) == MTisocov(kchan)) then
                do l = 1, Nencov
                  V2 = Vcov(ichan, j, ichan, j) * Vcov(kchan, l, kchan, l)
                  if (V2 > 0.) then
                    Ccov(ichan, j, kchan, l) = Vcov(ichan, j, kchan, l) / sqrt(V2)
                    term = Ccov(ichan, j, kchan, l) **2
                    if (term > Rlimit) then
                      Vcov(ichan, j, kchan, l) = Vcov(ichan, j, kchan, l) / term
                      Rcov(ichan, j, kchan, l) = Rcov(ichan, j, kchan, l) / term
                    endif
                  endif
                  if (Vcov(ichan, j, ichan, j) == 0.) Vcov(ichan, j, kchan, l) = 0.
                  if (Rcov(ichan, j, ichan, j) == 0.) Rcov(ichan, j, kchan, l) = 0.
                enddo
              endif
            enddo
          enddo
        enddo
        do j = 1, Nencov
          do kchan = 1, Nchancovint
            do l = 1, Nencov
              if (Vcov(ichan, j, ichan, j) == 0.) then
                Vcov(ichan, j, kchan, l) = 0.
                Vcov(kchan, l, ichan, j) = 0.
              endif
              if (Rcov(ichan, j, ichan, j) == 0.) then
                Rcov(ichan, j, kchan, l) = 0.
                Rcov(kchan, l, ichan, j) = 0.
              endif
            enddo
          enddo
        enddo
      endif
    enddo
  enddo
!
! Diagonal elements of covariance matrix (intra-channels) and average cross sections
!
  do i = 1, Nchanxs
    do j = 1, Nencov
      jj = Ecovindex(j)
      Swp = Sweightprev(jj)
      Sw = Sweight(jj)
      Sws = Sweightsum(jj)
      xsi0 = xseval(0, i, jj)
      xsi1 = min(xseval(1, i, jj), xsi0 * xslimit)
      xsi1 = max(xsi1,  xsi0 / xslimit)
      xsi0 = xsavC(i, j)
      xsdifi = xsi0 - xsi1
      do l = 1, Nencov
        ll = Ecovindex(l)
        xsk0 = xseval(0, i, ll)
        xsk1 = min(xseval(1, i, ll), xsk0 * xslimit)
        xsk1 = max(xsk1, xsk0 / xslimit)
        xsk0 = xsavC(i, l)
        xsdifk = xsk0 - xsk1
        term0 = xsdifi * xsdifk
        if (abs(term0) > coveps .and. Sws > 0.) then
          Vmt(i, j, l) = (Vmt(i, j, l) * Swp + Sw * term0) / Sws
          xs2 = xsi0 * xsk0
          if (xs2 > 0.) then
            term = term0 / xs2
            Rmt(i, j, l) = (Rmt(i, j, l) * Swp + Sw * term) / Sws
          endif
        endif
      enddo
      errmtC(i, j) = sqrt(Rmt(i, j, j))
    enddo
    do j = 1, Nen(i)
      Swp = Sweightprev(j)
      Sw = Sweight(j)
      Sws = Sweightsum(j)
      xsi0 = xseval(0, i, j)
      xsi1 = min(xseval(1, i, j), xsi0 * xslimit)
      xsi1 = max(xsi1, xsi0 / xslimit)
      xsi0 = xsav(i, j)
      xsdifi = xsi0 - xsi1
      term0 = xsdifi * xsdifi
      xs2 = xsi0 * xsi0
      if (abs(term0) > coveps .and. xs2 > 0.) then
        term = term0 / xs2
        if (Sws > 0.) RmtD(i, j) = (RmtD(i, j) * Swp + Sw * term) / Sws
        errmt(i, j) = sqrt(RmtD(i, j))
!
! Monte Carlo sensitivity matrix
!
        if (flagsens .and. .not. flagreadsens) then
          Np = min(Npar, numpar)
          n = Sindex(i, j)
          do k = 1, Np
            if (parsave(0, k) /= 0.) then
              parj0 = parsave(0, k)
              parjk = partalys(1, k)
              Ssign = 1.
            else
              parj0 = partalys(1, k)
              parjk = parsave(0, k)
              Ssign = -1.
            endif
            pardifj = parj0 - parjk
            if (pardifj /= 0.) then
              term = pardifj * Ssign * xsdifi
              if (Sws > 0.) Senum(k, i, n) = (Senum(k, i, n) * Swp + Sw * term) / Sws
              term = pardifj **2
              if (term == 0.) term = 1.
              if (Sws > 0.) Sdenom(k, i, n) = (Sdenom(k, i, n) * Swp + Sw * term) / Sws
              S(k, i, n) = Senum(k, i, n) / Sdenom(k, i, n) * parj0 / xsi0
            endif
          enddo
        endif
      endif
    enddo
  enddo
!
! Correlation matrix and setting covariance to zero if variance is zero
!
  do i = 1, Nchanxs
    do j = 1, Nencov
      do l = 1, Nencov
        V2 = Vmt(i, j, j) * Vmt(i, l, l)
        if (V2 > 0.) then
          Cmt(i, j, l) = Vmt(i, j, l) / sqrt(V2)
          term = Cmt(i, j, l) **2
          if (term > Rlimit) then
            Vmt(i, j, l) = Vmt(i, j, l) / term
            Rmt(i, j, l) = Rmt(i, j, l) / term
          endif
        endif
      enddo
    enddo
    do j = 1, Nencov
      do l = 1, Nencov
        if (Vmt(i, j, j) == 0.) then
          Vmt(i, j, l) = 0.
          Vmt(i, l, j) = 0.
        endif
        if (Rmt(i, j, j) == 0.) then
          Rmt(i, j, l) = 0.
          Rmt(i, l, j) = 0.
        endif
      enddo
    enddo
  enddo
!
! ******** Output of covariance matrix and average cross sections ******
!
! Covariance matrices for cross sections
!
  un=''
  col(1)='E_a'
  un(1)='MeV'
  col(2)='E_b'
  un(2)='MeV'
  col(3)='Rel._covariance'
  col(4)='Correlation'
  col(5)='Covariance'
  Ncol=5
!
! Inter-channel correlations
!
  open (unit = 1, file = 'cov_inter.ave', status = 'replace')
  do i = 1, Nchancovint
    do k = i, Nchancovint
      reaction=MTreac(MTcov(i))
      topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
      call write_header(indent,topline,source,user,date,oformat)
      call write_target(indent)
      call write_reaction(indent,reaction,0.D0,0.D0,3,MTcov(i))
      if (MTisocov(i) >= 0) call write_level(id2,MTisocov(i),-1,0.,-1.,0,0.)
      reaction_cov=MTreac(MTcov(k))
      call write_covariance(id2,reaction_cov,3,MTcov(i),3,MTcov(k),italys)
      if (MTisocov(k) >= 0) call write_level(id4,MTisocov(k),-1,0.,-1.,0,0.)
      call write_quantity(id2,quantity)
      call write_datablock(id2,Ncol,Nencov*Nencov,col,un)
      do j = 1, Nencov
        jj = Ecovindex(j)
        do l = 1, Nencov
          ll = Ecovindex(l)
!         write(3, '(2(3i4, es15.6), 3es15.6)') MTcov(i), MTisocov(i), jj, e(i, jj), MTcov(k), MTisocov(k), ll, &
!&          e(k, ll), Rcov(i, j, k, l), Ccov(i, j, k, l), Vcov(i, j, k, l)
          write(1, '(5es15.6)') e(i, jj), e(k, ll), Rcov(i, j, k, l), Ccov(i, j, k, l), Vcov(i, j, k, l)
        enddo
      enddo
    enddo
  enddo
  close (1)
!
! Intra-channel correlations
!
  open (unit = 1, file = 'cov_intra.ave', status = 'replace')
  do i = 1, Nchanxs
    reaction=reaction_string(i)
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
    call write_header(indent,topline,source,user,date,oformat)
    call write_target(indent)
    call write_reaction(indent,reaction,0.D0,0.D0,MF(i),MT(i))
    if (MTiso(i) >= 0) call write_level(id2,MTiso(i),-1,0.,-1.,0,0.)
    call write_covariance(id2,reaction,MF(i),MT(i),MF(i),MT(i),italys)
    call write_quantity(id2,quantity)
    call write_datablock(id2,Ncol,Nencov*Nencov,col,un)
    do j = 1, Nencov
      jj = Ecovindex(j)
      do l = 1, Nencov
        ll = Ecovindex(l)
!       write(4, '(2(3i4, es15.6), 3es15.6)') MT(i), MTiso(i), jj, e(i, jj), MT(i), MTiso(i), ll, e(i, ll), Rmt(i, j, l), &
!&        Cmt(i, j, l), Vmt(i, j, l)
        write(1, '(5es15.6)')  e(i, jj), e(i, ll), Rmt(i, j, l), Cmt(i, j, l), Vmt(i, j, l)
      enddo
    enddo
  enddo
  close (1)
!
! Diagonal elements of covariance matrices
!
  un=''
  col(1)='E'
  un(1)='MeV'
  col(2)='Rel._variance'
  col(3)='Rel._std._dev.'
  col(4)='Average_xs'
  un(4)='mb'
  col(5)='Std._dev.'
  un(5)='mb'
  Ncol=5
  open (unit = 1, file = 'variance.ave', status = 'replace')
  do i = 1, Nchanxs
    reaction=reaction_string(i)
    topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
    call write_header(indent,topline,source,user,date,oformat)
    call write_target(indent)
    call write_reaction(indent,reaction,0.D0,0.D0,MF(i),MT(i))
    if (MTiso(i) >= 0) call write_level(id2,MTiso(i),-1,0.,-1.,0,0.)
    call write_covariance(id2,reaction,3,MT(i),0,0,italys)
    call write_quantity(id2,quantity)
    call write_datablock(id2,Ncol,Nencov,col,un)
    do j = 1, Nencov
      jj = Ecovindex(j)
      err = xsavC(i, j) * errmtC(i, j)
!     write(3, '(3i4, 5es15.6)') MT(i), MTiso(i), jj, e(i, jj), Rmt(i, j, j), errmtC(i, j), xsavC(i, j), err
      write(1, '(5es15.6)') e(i, jj), Rmt(i, j, j), errmtC(i, j), xsavC(i, j), err
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
  headerline=''
  quantity='average cross section'
  do i = 1, Nchanxs
    ofile = xsfile(i)
    if (ofile(1:1) == ' ') cycle
    if (ofile(1:7) == 'fission') write(ofile(9:11), '("ave")')
    if (ofile(1:5) == 'total') write(ofile(7:9), '("ave")')
    if (ofile(1:7) == 'elastic') write(ofile(9:11), '("ave")')
    if (ofile(1:10) == 'nonelastic') write(ofile(12:14), '("ave")')
    if (ofile(1:2) == 'xs' .and. ofile(10:10) /= 'L') write(ofile(10:12), '("ave")')
    if (ofile(4:4) == 'L') write(ofile(7:10), '(".ave")')
    if (ofile(10:10) == 'L') write(ofile(13:16), '(".ave")')
    if (ofile(4:6) == 'con') write(ofile(7:10), '(".ave")')
    if (flagband) then
      Ncol=4
    else
      col(3)='xserr'
      Ncol=3
    endif
    j=0
    jheader=0
    open (unit = 1, file = xsfile(i), status = 'old', iostat = istat)
    if (istat /= 0) cycle
    do
      j=j+1
      read(1,'(a)',iostat = istat) line
      if (istat == -1) exit
      if (line(1:1) /= '#') exit
      headerline(j)=line
      key='title:'
      keyix=index(line,trim(key))
      if (keyix > 0) headerline(j)=trim(line)//' - average'
      key='source:'
      keyix=index(line,trim(key))
      if (keyix > 0) write(headerline(j)(keyix+len_trim(key)+1:80),'("TASMAN")')
      key='datablock'
      keyix=index(headerline(j),trim(key))
      if (keyix > 0) jheader = j-1
      key='entries'
      keyix=index(headerline(j),trim(key))
      if (keyix > 0) then
        write(headerline(j)(keyix+len_trim(key)+2:80), '(i6)') Nen(i)
        exit
      endif
    enddo
    close (1)
    open (unit = 1, file = ofile, status = 'replace')
    do k = 1, jheader
      write(1,'(a)') trim(headerline(k))
    enddo
    call write_quantity(indent,quantity)
    call write_datablock(indent,Ncol,Nen(i),col,un)
    do j = 1, Nen(i)
      xsA = xsav(i, j)
      if (average == 2) xsA = xseval(0, i, j)
      if (average == 3) xsA = xsopt(i, j)
      err = xsA * errmt(i, j)
      if (err /= 0.) err = min(err, xsA - 1.e-13)
      xslow(i, j) = xsA - err
      xsupp(i, j) = xsA + err
      if (flagband) then
        write(1, '(4es15.6)') E(i, j), xsA, xslow(i, j), xsupp(i, j)
      else
        write(1, '(3es15.6)') E(i, j), xsA, err
      endif
    enddo
    xslow(i, 0) = xslow(i, 1)
    xsupp(i, 0) = xsupp(i, 1)
    close (1)
  enddo
!
! ************************** Statistics of outliers ********************
!
! Count number of experimental data points inside and outside 1-sigma error band. This is for global statistics purposes.
!
  Nsetstot = 0
  Npointstot = 0
  Ninsidetot = 0
  Noutsidetot = 0
  fracinsidetot = 1.
  open (unit = 1, file = 'sigmaband.xs', status = 'replace')
  open (unit = 2, file = 'sigmaband.xsE', status = 'replace')
  open (unit = 3, file = 'sigmaband.E', status = 'replace')
  write(1, '("# 1-sigma fractions per channel after ", i6, " runs")') italys
  write(2, '("# 1-sigma fractions per channel and energy after ", i6, " runs")') italys
  write(3, '("# 1-sigma fractions per energy after ", i6, " runs")') italys
  do ifile = 1, 3
    write(ifile, '("# Z=", i3, " A=", i3)') Ztarget, Atarget
  enddo
  write(1, '("#")')
  write(3, '("#", 30x, " Energies:", i4)') Nen(1)
  write(1, '("#  Channel         Sets Points Inside ", "Outside Inside fraction")')
  write(3, '("#  Energy (MeV)     Points Inside ", "Outside Inside fraction")')
  do k = 1, numenin
     NpointsE(k) = 0
     NinsideE(k) = 0
     NoutsideE(k) = 0
     NpointsE(k) = 0
     fracinsideE(k) = 1.
  enddo
  do i = 1, Nchanxs
    Npoints = 0
    Ninside = 0
    Noutside = 0
    fracinside = 1.
    if (Nsets(i) == 0) cycle
    write(2, '("# Channel: ", a20, " Energies:", i4)') xsfile(i), Nen(1)
    write(2, '("#  Energy (MeV)     Points Inside ", "Outside Inside fraction")')
    Nsetstot = Nsetstot + Nsets(i)
    ee(0) = 0.
    do k = 1, Nen(1)
      ee(k) = E(i, k)
      NpointsxsE(k) = 0
      NinsidexsE(k) = 0
      NoutsidexsE(k) = 0
      fracinsidexsE(k) = 1.
    enddo
    do j = 1, Nsets(i)
      do k = 1, Nenexp(i, j)
        Eex = Eexp(i, j, k)
        m = Eindex(i, j, k)
        xsex = xsexp(i, j, k)
        if (Eex > ee(Nen(i))) cycle
        call locate(ee, 0, Nen(i), Eex, NE)
        e1 = ee(NE)
        e2 = ee(NE + 1)
        fac = (Eex - e1) / (e2 - e1)
        xslowint = xslow(i, NE) + fac * (xslow(i, NE + 1) - xslow(i, NE))
        xsuppint = xsupp(i, NE) + fac * (xsupp(i, NE + 1) - xsupp(i, NE))
        if (xsex <= xsuppint .and. xsex >= xslowint) then
          Ninside = Ninside + 1
          NinsideE(m) = NinsideE(m) + 1
          NinsidexsE(m) = NinsidexsE(m) + 1
          Ninsidetot = Ninsidetot + 1
        else
          Noutside = Noutside + 1
          NoutsideE(m) = NoutsideE(m) + 1
          NoutsidexsE(m) = NoutsidexsE(m) + 1
          Noutsidetot = Noutsidetot + 1
        endif
        Npoints = Npoints + 1
        NpointsE(m) = NpointsE(m) + 1
        NpointsxsE(m) = NpointsxsE(m) + 1
        Npointstot = Npointstot + 1
      enddo
    enddo
    if (Npoints > 0) then
      fracinside = Ninside / real(Npoints)
      write(1, '("# ", a15, 4i6, f12.5)') xsfile(i), Nsets(i), Npoints, Ninside, Noutside, fracinside
    endif
    do k = 1, Nen(1)
      if (NpointsxsE(k) > 0) fracinsidexsE(k) = NinsidexsE(k) / real(NpointsxsE(k))
      write(2, '(2x, es12.5, 6x, 3i6, f12.5)') E(1, k), NpointsxsE(k), NinsidexsE(k), NoutsidexsE(k), fracinsidexsE(k)
    enddo
    write(2, '("#  Energy (MeV)     Points Inside ", "Outside Inside fraction")')
    write(2, '("#   Total     ", 6x, 3i6, f12.5)') Npoints, Ninside, Noutside, fracinside
  enddo
  do k = 1, Nen(1)
    if (NpointsE(k) > 0) fracinsideE(k) = NinsideE(k) / real(NpointsE(k))
    write(3, '(2x, es12.5, 6x, 3i6, f12.5)') E(1, k), NpointsE(k), NinsideE(k), NoutsideE(k), fracinsideE(k)
  enddo
  if (Npointstot > 0) fracinsidetot = Ninsidetot / real(Npointstot)
  write(1, '("#  Channel         Sets Points Inside ", "Outside Inside fraction")')
  write(3, '("#  Energy (MeV)     Points Inside ", "Outside Inside fraction")')
  write(1, '("#   Total        ", 4i6, f12.5)') Nsetstot, Npointstot, Ninsidetot, Noutsidetot, fracinsidetot
  write(3, '("#   Total     ", 6x, 3i6, f12.5)') Npointstot, Ninsidetot, Noutsidetot, fracinsidetot
  do ifile = 1, 3
    close(ifile)
  enddo
  return
end subroutine covcross
! Copyright A.J. Koning 2021
