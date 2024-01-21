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
  integer           :: i                                   ! counter
  integer           :: j                                   ! counter
  integer           :: k                                   ! Legendre order
  integer           :: l                                   ! counter
  integer           :: n                                   ! counter
  integer           :: Np                                  ! number of points
  real(sgl)         :: dEtot                               ! total energy range
  real(sgl)         :: factor                              ! help variable
  real(sgl)         :: pardifj                             ! difference between parameters
  real(sgl)         :: parj0                               ! parameter of run 0
  real(sgl)         :: parjk                               ! parameter of random run
  real(sgl)         :: ploc(numpar)                        ! help variable
  real(sgl)         :: ptmp                                ! help variable
  real(sgl)         :: Sabs                                ! absolute sensitivity
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
  character(len=33) :: format1                             ! format string
  character(len=26) :: format2                             ! format string
  character(len=132):: strloc(numpar)                      ! help variable
  character(len=132):: sttmp                               ! help variable
!
! ************************* Create S matrix ****************************
!
! Linear sensitivity matrix
!
  if ( .not. flagsens .or. flagreadsens) then
    Np = Npar
  else
    Np = min(italys, numpar)
  endif
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
              S(k,i,n) = factor * xsdifi / pardifj
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
  open (unit = 3, file = 'sensitivity.par', status = 'replace')
  write(3, '("#Number of parameters :", i4)') Np
  do k = 1, Np
    write(3, '("#Parameter ", i4, "       : ", a30)') k, parstring(k)
    write(3, '("#Parameter uncertainty:", f12.5, " %")') 100. * pardelta(k)
    write(3, '("#Number of channels   :", i4)') Nchanxs
    do i = 1, Nchanxs
      write(3, '("#Parameter ", i4, "       : ", a30)') k, parstring(k)
      write(3, '("#Parameter uncertainty:", f12.5, " %")') 100. * pardelta(k)
      write(3, '("#Channel   ", i4, "       : ", a20)') i, xsfile(i)
      write(3, '("#Number of energies   :", i4)') Nen(i)
      write(3, '("#      E            S         ", "max. dev. (mb)  cross section (mb)")')
      do j = 1, Nen(i)
        n = Sindex(i, j)
        write(3, '(es12.5, 3(es15.6))') E(i, j), S(k, i, n), xsdev(k, i, n), xstalys(0, i, j)
      enddo
    enddo
  enddo
  close (3)
!
! Sensitivity per cross section, energy and parameter, sorted
!
  open (unit = 3, file = 'sensitivity.xs', status = 'replace')
  write(3, '("#Number of channels   :", i4)') Nchanxs
  do i = 1, Nchanxs
    write(3, '("#Channel   ", i4, "       : ", a20)') i, xsfile(i)
    write(3, '("#Number of energies   :", i4)') Nen(i)
    do j = 1, Nen(i)
      n = Sindex(i, j)
      write(3, '("#Channel   ", i4, "       : ", a20)') i, xsfile(i)
      write(3, '("#Energy               : ", es12.5 , " MeV")') E(i, j)
      write(3, '("#Cross section        : ", es12.5 , " mb")') xstalys(0, i, j)
      write(3, '("#Number of parameters :", i4)') Np
      write(3, '("#Parameter                            S         ", "max. dev. (mb)   pardev(%)")')
      do k = 1, Np
        strloc(k) = parstring(k)
        ploc(k) = pardelta(k)
        Sloc(k) = S(k, i, n)
        xsloc(k) = xsdev(k, i, n)
      enddo
      do k = 1, Np
        do l = k, Np
           if (abs(Sloc(k)) >= abs(Sloc(l))) cycle
           ptmp = ploc(k)
           Stmp = Sloc(k)
           xstmp = xsloc(k)
           sttmp = strloc(k)
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
      do k = 1, Np
        write(3, '(a30, 2es15.6, f12.5)') strloc(k)(1:30), Sloc(k), xsloc(k), 100. * ploc(k)
      enddo
    enddo
  enddo
  close (3)
!
! Sensitivity per energy, parameter and cross section
!
  open (unit = 3, file = 'sensitivity.E', status = 'replace')
  write(3, '("#Number of energies   :", i4)') Nen(1)
  do j = 1, Nen(1)
    write(3, '("#Energy               : ", es12.5 , " MeV")') e(1, j)
    write(3, '("#Number of parameters :", i4)') Np
    do k = 1, Np
      write(3, '("#Energy               : ", es12.5 , " MeV")') e(1, j)
      write(3, '("#Parameter ", i4, "       : ", a30)') k, parstring(k)
      write(3, '("#Parameter uncertainty:", f12.5, " %")') 100. * pardelta(k)
      write(3, '("#Number of channels   :", i4)') Nchanxs
      write(3, '("#Channel                    S         ", "max. dev. (mb) Cross section (mb)")')
      do i = 1, Nchanxs
        n = Sindex(i, j)
        write(3, '(a20, 3(es15.6))') xsfile(i), S(k, i, n), xsdev(k, i, n), xstalys(0, i, j)
      enddo
    enddo
  enddo
  close (3)
!
! Sensitivity per parameter integrated over energy, sorted
!
  open (unit = 1, file = 'sensitivity.parave', status = 'replace')
  write(1, '("#Number of parameters :", i4)') Np
  do k = 1, Np
    write(1, '("#Parameter ", i4, "       : ", a30)') k, parstring(k)
    write(1, '("#Parameter uncertainty:", f12.5, " %")') 100. * pardelta(k)
    write(1, '("#Number of channels   :", i4)') Nchanxs
    write(1, '("#Channel                    S         ", "max. dev. (mb)")')
    do i = 1, Nchanxs
      write(1, '(a20, 2(es15.6))') xsfile(i), Schan(k, i), xsdevav(k, i)
    enddo
  enddo
  close (1)
!
! Sensitivity per channel integrated over energy, sorted
!
  open (unit = 1, file = 'sensitivity.xsave', status = 'replace')
  write(1, '("#Number of channels   :", i4)') Nchanxs
  do i = 1, Nchanxs
    write(1, '("#Channel   ", i4, "       : ", a20)') i, xsfile(i)
    write(1, '("#Number of parameters :", i4)') Np
    write(1, '("#Parameter                            S         ", "max. dev. (mb)   pardev(%)")')
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
    do k = 1, Np
      write(1, '(a30, 2es15.6, f12.5)') strloc(k)(1:30), Sloc(k), xsloc(k), 100. * ploc(k)
    enddo
  enddo
  close (1)
!
! Sensitivity per energy, summed over channels
!
  open (unit = 3, file = 'sensitivity.Eave', status = 'replace')
  write(3, '("#Number of energies  :", i4)') Nen(1)
  do j = 1, Nen(1)
    write(3, '("#Energy               : ", es12.5 , " MeV")') e(1, j)
    write(3, '("#Number of parameters :", i4)') Np
    write(3, '("#Parameter                            S    ")')
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
    do k = 1, Np
      write(3, '(a30, es15.6)') strloc(k)(1:30), Sloc(k)
    enddo
  enddo
  close (3)
!
! Sensitivity integrated over all channels
!
  open (unit = 2, file = 'sensitivity.all', status = 'replace')
  write(2, '("#Number of parameters :", i4)') Np
  write(2, '("#Parameter                            S    ")')
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
    write(2, '(a30, es15.6)') strloc(k)(1:30), Sloc(k)
  enddo
  close (2)
!
! Output of sensitivity profile per reaction channel
!
  if (Np == 0) return
  format1 = '("#E                   ",xxxa30))'
  write(format1(26:28), '(i3.3)') Np
  format2='(es12.5,xxx(es20.5,10x))'
  write(format2(9:11), '(i3.3)') Np
  do i = 1, Nchanxs
    open (unit = 3, file = 'sens_'//xsfile(i), status = 'replace')
    write(3, fmt = format1) (parstring(k)(1:30), k = 1, Np)
    do j = 1, Nen(i)
      n = Sindex(i, j)
      write(3, fmt = format2) E(i, j), (S(k, i, n), k = 1, Np)
    enddo
    close (3)
  enddo
  return
end subroutine sensitivity
! Copyright A.J. Koning 2021
