subroutine sensitivity
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Create sensitivity matrix
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
!   sgl             ! single precision kind
! All global variables
!   numchanxs       ! maximum number of channels with cross sections
!   numenin         ! maximum number of incident energies
!   numenS          ! maximum number of energies for sensitivities
!   numpar          ! maximum number of parameters
! Variables for writing TALYS input files
!   italys          ! TALYS run counter
!   mode            ! TASMAN mode
! Variables for parameter covariances
!   parstring       ! parameter string
! Variables for cross section covariances
!   S               ! sensitivity matrix
! Variables for parameter variation
!   Npar            ! number of parameters
!   pardelta        ! uncertainty of parameter
!   parsave         ! parameter value
!   partalys        ! parameter value
! Variables for reading cross sections
!   dE              ! uncertainty of incident energy
!   E               ! incident energy
!   Nchanxs         ! total number of channels with cross sections
!   Nen             ! number of incident energies
!   Sindex          ! index for cross section
!   xssave          ! cross section from TALYS
!   xsfile          ! file with crosssections
!   xstalys         ! cross section from TALYS
! Variables for sensitivity
!   flagreadsens    ! flag to read sensitivities from tables
!   flagsens        ! flag to use sensitivity in parameter weighting
!   Sall            ! sensitivity for all channels
!   Schan           ! sensitivity for channel
!   SE              ! sensitivity for energy
!   xsdevav         ! average deviation of cross section
!
! *** Declaration of local data
!
  implicit none
  character(len=15) :: col(numpar)                         ! header
  character(len=15) :: un(numpar)
  character(len=16) :: reaction
  character(len=132) :: quantity
  character(len=132) :: topline    ! topline
  character(len=132) :: sensfile
  integer           :: i                                   ! counter
  integer           :: j                                   ! counter
  integer           :: Ncol
  integer           :: k                                   ! Legendre order
  integer           :: l                                   ! counter
  integer           :: n                                   ! counter
  integer           :: Np                                  ! number of points
  integer           :: indent
  integer           :: id2
  integer           :: id4
  integer           :: id6
  real(sgl)         :: dEtot                               ! total energy range
  real(sgl)         :: factor                              ! help variable
  real(sgl)         :: pardifj                             ! difference between parameters
  real(sgl)         :: parj0                               ! parameter of run 0
  real(sgl)         :: parjk                               ! parameter of random run
  real(sgl)         :: ploc(numpar)                        ! help variable
  real(sgl)         :: pearloc(numpar)                     ! help variable
  real(sgl)         :: ptmp                                ! help variable
  real(sgl)         :: peartmp                             ! help variable
  real(sgl)         :: Sabs                                ! absolute sensitivity
  real(sgl)         :: sens                                ! sensitivity
  real(sgl)         :: Sloc(numpar)                        ! help variable
  real(sgl)         :: Stmp                                ! help variable
  real(sgl)         :: sumdev                              ! summed deviation
  real(sgl)         :: sumE(numenin)                       ! summed sensitivity x cross section
  real(sgl)         :: sumS                                ! integrated sensitivity
  real(sgl)         :: sumSxs                              ! integrated sensitivity
  real(sgl)         :: sumxs(numenin)                      ! sum over emission channels
  real(sgl)         :: sumxsall                            ! sumxsE summed over channels
  real(sgl)         :: sumxsE                              ! summed cross section x energy
  real(sgl)         :: xsdev(numpar, numchanxs, numenS)    ! cross section deviation
  real(sgl)         :: xsdifi                              ! difference in cross section
  real(sgl)         :: xsi0                                ! cross section of run 0
  real(sgl)         :: xsik                                ! cross section of random run
  real(sgl)         :: xsloc(numpar)                       ! help variable
  real(sgl)         :: xstmp                               ! help variable
  real(sgl)         :: Sprev                               ! help variable
  character(len=33) :: format1                             ! format string
  character(len=26) :: format2                             ! format string
  character(len=132):: strloc(numpar)                      ! help variable
  character(len=132):: sttmp                               ! help variable
!
! ************************* Create S matrix ****************************
!
! Linear sensitivity matrix
!
  indent = 0
  id2 = indent + 2
  id4 = indent + 4
  id6 = indent + 6
  quantity = 'sensitivity coefficients'
  if ( .not. flagsens .or. flagreadsens) then
    Np = Npar
  else
    Np = min(italys, numpar)
  endif
  Np = min(Npar, numpar)
  if (mode == 2) then
    do k = 1, Np
      parj0 = partalys(0, k)
      parjk = parsave(k, k)
      pardifj = parjk - parj0
      if (pardifj /= 0.) then
        do i = 1, Nchanxs
          do j = 1, Nen(i)
            n = Sindex(i,j)
            xsi0 = xstalys(0,i,j)
            xsik = xssave(k,i,n)
            if (xsi0 /= 0. .and. xsik /= 0.) then
              xsdifi = xsik - xsi0
              if (parj0 == 0.) then
                factor = parjk / xsik
              else
                factor = parj0 / xsi0
              endif
              sens = factor * xsdifi / pardifj
              if (flagmorris .and. iloop > 1) then
                if (italys == k) then
                  Sprev = S(k,i,n)
                  S(k,i,n) = (Sprev*(iloop - 1)+sens) / iloop
                endif
              else
                S(k,i,n) = sens
                Sprev = sens
              endif
            endif
          enddo
        enddo
      endif
    enddo
  endif
!
! *************** Integrals over energy and channels *******************
!
  do k = 1, Np
    sumSxs = 0.
    sumxsall = 0.
    do j = 1, numenin
      sumE(j) = 0.
      sumxs(j) = 0.
    enddo
    do i = 1, Nchanxs
      sumS = 0.
      sumxsE = 0.
      sumdev = 0.
      do j = 1, Nen(i)
        n = Sindex(i, j)
        Sabs = abs(S(k, i, n))
        sumS = sumS + Sabs * dE(i, j)
        sumxsE = sumxsE + xstalys(0, i, j) * dE(i, j)
        sumxs(j) = sumxs(j) + xstalys(0, i, j)
        sumE(j) = sumE(j) + Sabs * xstalys(0, i, j)
        xsdev(k, i, n) = S(k, i, n) * xstalys(0, i, j) * pardelta(k)
        sumdev = sumdev + abs(xsdev(k, i, n)) * dE(i, j)
      enddo
      dEtot = E(i, nen(i)) - E(i, 0)
      if (dEtot > 0.) then
        Schan(k, i) = sumS / dEtot
        xsdevav(k, i) = sumdev / dEtot
      endif
      sumSxs = sumSxs + Schan(k, i) * sumxsE
      sumxsall = sumxsall + sumxsE
    enddo
    do j = 1, numenin
      if (sumxs(j) > 0.) SE(k, j) = sumE(j) / sumxs(j)
    enddo
    if (sumxsall > 0.) Sall(k) = sumSxs / sumxsall
  enddo
!
! *************************** Output ***********************************
!
! Sensitivity per parameter, cross section and energy
!
  topline=trim(targetnuclide)//' '//trim(quantity)//' per parameter, cross section and energy'
  sensfile = 'sensitivity.par'
! if (flagmorris) write(sensfile(len_trim(sensfile)+1:len_trim(sensfile)+4),'("_",i3.3)') iloop
  open (unit = 1, file = sensfile, status = 'replace')
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  un = ''
  col(1) = 'Energy'
  un(1) = 'MeV'
  col(2) = 'S'
  col(3) = 'max. deviation'
  un(3) = 'mb'
  col(4) = 'xs'
  un(4) = '%'
  Ncol = 4
  call write_char(indent,'parameters','')
  call write_integer(id2,'number of parameters',Np)
  do k = 1, Np
    call write_char(id4,'parameter',parstring(k))
    call write_real(id4,'parameter variation [%]',100. * pardelta(k))
    call write_integer(id4,'number of channels',Nchanxs)
    do i = 1, Nchanxs
      call write_char(id6,'channel',xsfile(i))
      call write_char(id6,'reaction',reaction_string(i))
      call write_quantity(indent,quantity)
      call write_datablock(indent,Ncol,Nen(i),col,un)
      do j = 1, Nen(i)
        n = Sindex(i, j)
        write(1, '(4es15.6)') E(i, j), S(k, i, n), xsdev(k, i, n), xstalys(0, i, j)
      enddo
    enddo
  enddo
  close (1)
!
! Sensitivity per cross section, energy and parameter, sorted
!
  topline=trim(targetnuclide)//' '//trim(quantity)//' per cross section, energy and parameter, sorted'
  sensfile = 'sensitivity.xs'
! if (flagmorris) write(sensfile(len_trim(sensfile)+1:len_trim(sensfile)+4),'("_",i3.3)') iloop
  open (unit = 1, file = sensfile, status = 'replace')
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  un = ''
  col(1) = 'Parameter'
  col(2) = ''
  col(3) = 'S'
  col(4) = 'max. deviation'
  un(4) = 'mb'
  col(5) = 'par. deviation'
  un(5) = '%'
  col(6) = 'Pearson corr.'
  un(6) = ''
  Ncol = 6
  call write_char(indent,'parameters','')
  call write_integer(id2,'number of channels',Nchanxs)
  do i = 1, Nchanxs
    call write_char(id4,'channel',xsfile(i))
    call write_char(id4,'reaction',reaction_string(i))
    call write_integer(id4,'number of energies',Nen(i))
    do j = 1, Nen(i)
      n = Sindex(i, j)
      call write_real(id6,'E-incident [MeV]',e(i, j))
      call write_real(id6,'Cross section [mb]',xstalys(0, i, j))
      do k = 1, Np
        strloc(k) = parstring(k)
        ploc(k) = pardelta(k)
        Sloc(k) = S(k, i, n)
        pearloc(k) = Pearson(k, i, n)
        xsloc(k) = xsdev(k, i, n)
      enddo
      do k = 1, Np
        do l = k, Np
           if (abs(Sloc(k)) >= abs(Sloc(l))) cycle
           ptmp = ploc(k)
           peartmp = pearloc(k)
           Stmp = Sloc(k)
           xstmp = xsloc(k)
           sttmp = strloc(k)
           ploc(k) = ploc(l)
           pearloc(k) = pearloc(l)
           Sloc(k) = Sloc(l)
           xsloc(k) = xsloc(l)
           strloc(k) = strloc(l)
           ploc(l) = ptmp
           pearloc(l) = peartmp
           Sloc(l) = Stmp
           xsloc(l) = xstmp
           strloc(l) = sttmp
        enddo
      enddo
      call write_quantity(indent,quantity)
      call write_datablock(indent,Ncol,Np,col,un)
      do k = 1, Np
        write(1, '(a30, 4es15.6)') strloc(k)(1:30), Sloc(k), xsloc(k), 100. * ploc(k), pearloc(k)
      enddo
    enddo
  enddo
  close (1)
!
! Sensitivity per energy, parameter and cross section
!
  topline=trim(targetnuclide)//' '//trim(quantity)//' per energy, parameter and cross section'
  sensfile = 'sensitivity.E'
! if (flagmorris) write(sensfile(len_trim(sensfile)+1:len_trim(sensfile)+4),'("_",i3.3)') iloop
  open (unit = 1, file = sensfile, status = 'replace')
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  un = ''
  col(1) = 'Channel'
  col(2) = 'File'
  col(3) = 'S'
  col(4) = 'max. deviation'
  un(4) = 'mb'
  col(5) = 'xs'
  un(5) = 'mb'
  Ncol = 5
  call write_char(indent,'parameters','')
  call write_integer(id2,'number of energies',Nen(1))
  do j = 1, Nen(1)
    call write_real(id4,'E-incident [MeV]',e(1, j))
    call write_integer(id4,'number of parameters',Np)
    do k = 1, Np
      call write_char(id6,'parameter',parstring(k))
      call write_real(id6,'parameter variation [%]',100. * pardelta(k))
      call write_quantity(indent,quantity)
      call write_datablock(indent,Ncol,Nchanxs,col,un)
      do i = 1, Nchanxs
        n = Sindex(i, j)
        write(1, '(2a15, 3es15.6)') reaction_string(i),xsfile(i), S(k, i, n), xsdev(k, i, n), xstalys(0, i, j)
      enddo
    enddo
  enddo
  close (1)
!
! Sensitivity per parameter integrated over energy, sorted
!
  topline=trim(targetnuclide)//' '//trim(quantity)//' per parameter integrated over energy, sorted'
  sensfile = 'sensitivity.parave'
! if (flagmorris) write(sensfile(len_trim(sensfile)+1:len_trim(sensfile)+4),'("_",i3.3)') iloop
  open (unit = 1, file = sensfile, status = 'replace')
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  un = ''
  col(1) = 'Channel'
  col(2) = 'File'
  col(3) = 'S'
  col(4) = 'max. deviation'
  un(4) = 'mb'
  Ncol = 4
  call write_char(indent,'parameters','')
  call write_integer(id2,'number of parameters',Np)
  do k = 1, Np
    call write_char(id4,'parameter',parstring(k))
    call write_real(id4,'parameter variation [%]',100. * pardelta(k))
    call write_quantity(indent,quantity)
    call write_datablock(indent,Ncol,Nchanxs,col,un)
    do i = 1, Nchanxs
      write(1, '(2a15, 2es15.6)') reaction_string(i), xsfile(i), Schan(k, i), xsdevav(k, i)
    enddo
  enddo
  close (1)
!
! Sensitivity per channel integrated over energy, sorted
!
  topline=trim(targetnuclide)//' '//trim(quantity)//' per channel integrated over energy, sorted'
  sensfile = 'sensitivity.xsave'
! if (flagmorris) write(sensfile(len_trim(sensfile)+1:len_trim(sensfile)+4),'("_",i3.3)') iloop
  open (unit = 1, file = sensfile, status = 'replace')
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  un = ''
  col(1) = 'Parameter'
  col(2) = ''
  col(3) = 'S'
  col(4) = 'max. deviation'
  un(4) = 'mb'
  col(5) = 'par. deviation'
  un(5) = '%'
  Ncol = 5
  call write_char(indent,'parameters','')
  call write_integer(id2,'number of channels',Nchanxs)
  do i = 1, Nchanxs
    call write_char(id4,'channel',xsfile(i))
    call write_char(id4,'reaction',reaction_string(i))
    do k = 1, Np
      strloc(k) = parstring(k)
      ploc(k) = pardelta(k)
      Sloc(k) = Schan(k, i)
      xsloc(k) = xsdevav(k, i)
    enddo
    do k = 1, Np
      do l = k, Np
        if (abs(Sloc(k)) >= abs(Sloc(l))) cycle
        ptmp = ploc(k)
        Stmp = Sloc(k)
        sttmp = strloc(k)
        xstmp = xsloc(k)
        ploc(k) = ploc(l)
        Sloc(k) = Sloc(l)
        xsloc(k) = xsloc(l)
        strloc(k) = strloc(l)
        ploc(l) = ptmp
        Sloc(l) = Stmp
        xsloc(l) = xstmp
        strloc(l) = sttmp
      enddo
    enddo
    call write_quantity(indent,quantity)
    call write_datablock(indent,Ncol,Np,col,un)
    do k = 1, Np
      write(1, '(a30, 3es15.6)') strloc(k)(1:30), Sloc(k), xsloc(k), 100. * ploc(k)
    enddo
  enddo
  close (1)
!
! Sensitivity per energy, summed over channels
!
  topline=trim(targetnuclide)//' '//trim(quantity)//' per energy, integrated over channels'
  sensfile = 'sensitivity.Eave'
! if (flagmorris) write(sensfile(len_trim(sensfile)+1:len_trim(sensfile)+4),'("_",i3.3)') iloop
  open (unit = 1, file = sensfile, status = 'replace')
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  un = ''
  col(1) = 'Parameter'
  col(2) = ''
  col(3) = 'S'
  Ncol = 3
  call write_char(indent,'parameters','')
  call write_integer(id2,'number of energies',Nen(1))
  do j = 1, Nen(1)
    call write_real(id4,'E-incident [MeV]',e(1, j))
    do k = 1, Np
      strloc(k) = parstring(k)
      Sloc(k) = SE(k, j)
    enddo
    do k = 1, Np
      do l = k, Np
        if (abs(Sloc(k)) >= abs(Sloc(l))) cycle
        Stmp = Sloc(k)
        sttmp = strloc(k)
        Sloc(k) = Sloc(l)
        strloc(k) = strloc(l)
        Sloc(l) = Stmp
        strloc(l) = sttmp
      enddo
    enddo
    call write_quantity(indent,quantity)
    call write_datablock(indent,Ncol,Np,col,un)
    do k = 1, Np
      write(1, '(a30, es15.6)') strloc(k)(1:30), Sloc(k)
    enddo
  enddo
  close (1)
!
! Sensitivity integrated over all channels
!
  topline=trim(targetnuclide)//' '//trim(quantity)//' integrated over all channels and energies'
  sensfile = 'sensitivity.all'
! if (flagmorris) write(sensfile(len_trim(sensfile)+1:len_trim(sensfile)+4),'("_",i3.3)') iloop
  open (unit = 1, file = sensfile, status = 'replace')
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  un = ''
  col(1) = 'Parameter'
  col(2) = ''
  col(3) = 'S'
  Ncol = 3
  call write_quantity(indent,quantity)
  call write_datablock(indent,Ncol,Np,col,un)
  do k = 1, Np
    strloc(k) = parstring(k)
    Sloc(k) = Sall(k)
  enddo
  do k = 1, Np
    do l = k, Np
      if (abs(Sloc(k)) >= abs(Sloc(l))) cycle
      Stmp = Sloc(k)
      sttmp = strloc(k)
      Sloc(k) = Sloc(l)
      strloc(k) = strloc(l)
      Sloc(l) = Stmp
      strloc(l) = sttmp
    enddo
  enddo
  do k = 1, Np
    write(1, '(a30, es15.6)') strloc(k)(1:30), Sloc(k)
  enddo
  close (1)
!
! Output of sensitivity profile per reaction channel
!
  if (Np == 0) return
  format1 = '("#E                   ",xxxa30))'
  write(format1(26:28), '(i3.3)') Np
  format2='(es15.5,xxxes15.6)'
  write(format2(9:11), '(i3.3)') Np
  do i = 1, Nchanxs
    reaction=reaction_string(i)
    topline=trim(targetnuclide)//' '//trim(reaction)//' '//trim(quantity)//' per reaction channel'
    sensfile = 'sens_'//xsfile(i)
!   if (flagmorris) write(sensfile(len_trim(sensfile)+1:len_trim(sensfile)+4),'("_",i3.3)') iloop
    open (unit = 1, file = sensfile, status = 'replace')
    call write_header(indent,topline,source,user,date,oformat)
    call write_target(indent)
    call write_reaction(indent,reaction,0.D0,0.D0,0,0)
    un = ''
    col(1) = 'E-incident'
    un(1) = 'MeV'
    do k = 1, Np
      col(k+1) = parstring(k)(1:15)
    enddo
    Ncol = Np + 1
    call write_quantity(indent,quantity)
    call write_datablock(indent,Ncol,Nen(i),col,un)
    do j = 1, Nen(i)
      n = Sindex(i, j)
      write(1, fmt = format2) E(i, j), (S(k, i, n), k = 1, Np)
    enddo
    close (1)
  enddo
  return
end subroutine sensitivity
! Copyright A.J. Koning 2021
