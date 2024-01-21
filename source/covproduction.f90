subroutine covproduction
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for isotope production
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
! Variables for writing TALYS input files
!   italys          ! TALYS run counter
! Variables for cross section covariances
!   coveps          ! limit for covariance
!   flagband        ! flag to represent results as error bands
! Variables for reading isotope production
!   acttalys        ! activity of produced isotope in MBq
!   NchanY          ! total number of channels with isotope production
!   Nisoreltalys    ! fraction of number of produced isotopes per element
!   Nisotalys       ! number of isotopes produced after irradiation
!   Ntime           ! number of time points
!   timeY           ! time point
!   Yfile           ! name of file with isotope production
!   yieldtalys      ! yield of produced isotope in MBq / (mA.h)
! Variables for isotope production covariances
!   actav           ! average activity, ...
!   erract          ! activity uncertainty, ...
!   errNiso         ! error in diagonal of average number of isotopes produced after irradiation
!   errNisorel      ! error in fraction of number of produced isotopes per element
!   erryield        ! error in yield
!   Nisoav          ! average number of isotopes produced after irradiation
!   Nisorelav       ! average fraction of number of produced isotopes per element
!   RactD           ! diagonal of covariance matrix for activity, ...
!   RNisoD          ! diagonal of average number of isotopes produced after irradiation
!   RNisorelD       ! diagonal in fraction of number of produced isotopes per element
!   RyieldD         ! diagonal of yield covariance matrix
!   yieldav         ! average yield of produced isotope in MBq / (mA.h)
! Variables for covariances
!   Sw0             ! sum of weights
!   Swp0            ! sum of weights
!   Sws0            ! sum of weights
!
! *** Declaration of local data
!
  implicit none
  character(len=21) :: ofile     ! output file
  integer           :: i         ! counter
  integer           :: j         ! counter
  real(sgl)         :: err       ! error
  real(sgl)         :: term      ! help variable
  real(sgl)         :: term0     ! help variable
  real(sgl)         :: xs2       ! help variable
  real(sgl)         :: xsdifi    ! difference in cross section
  real(sgl)         :: xsi0      ! cross section of run 0
  real(sgl)         :: xsi1      ! cross section of random run
!
! Average isotope production and covariances
!
  do i = 1, NchanY
    do j = 1, Ntime(i)
      xsi0 = acttalys(0, i, j)
      xsi1 = min(acttalys(1, i, j), 10. * xsi0)
      xsdifi = xsi0 - xsi1
      term0 = xsdifi * xsdifi
      if (abs(term0) > coveps .and. Sws0 > 0.) then
        xs2 = xsi0 * xsi0
        if (xs2 > 0.) then
          term = term0 / xs2
          RactD(i, j) = (RactD(i, j) * Swp0 + Sw0 * term) / Sws0
          erract(i, j) = sqrt(RactD(i, j))
        endif
      endif
      if (Sws0 > 0.) actav(i, j) = (actav(i, j) * Swp0 + Sw0 * xsi1) / Sws0
    enddo
    do j = 1, Ntime(i)
      xsi0 = Nisotalys(0, i, j)
      xsi1 = min(Nisotalys(1, i, j), 10. * xsi0)
      xsdifi = xsi0 - xsi1
      term0 = xsdifi * xsdifi
      if (abs(term0) > coveps .and. Sws0 > 0.) then
        xs2 = xsi0 * xsi0
        if (xs2 > 0.) then
          term = term0 / xs2
          RNisoD(i, j) = (RNisoD(i, j) * Swp0 + Sw0 * term) / Sws0
          errNiso(i, j) = sqrt(RNisoD(i, j))
        endif
      endif
      if (Sws0 > 0.) Nisoav(i, j) = (Nisoav(i, j) * Swp0 + Sw0 * xsi1) / Sws0
    enddo
    do j = 1, Ntime(i)
      xsi0 = yieldtalys(0, i, j)
      xsi1 = min(yieldtalys(1, i, j), 10. * xsi0)
      xsdifi = xsi0 - xsi1
      term0 = xsdifi * xsdifi
      if (abs(term0) > coveps .and. Sws0 > 0.) then
        xs2 = xsi0 * xsi0
        if (xs2 > 0.) then
          term = term0 / xs2
          RyieldD(i, j) = (RyieldD(i, j) * Swp0 + Sw0 * term) / Sws0
          erryield(i, j) = sqrt(RyieldD(i, j))
        endif
      endif
      if (Sws0 > 0.) yieldav(i, j) = (yieldav(i, j) * Swp0 + Sw0 * xsi1) / Sws0
    enddo
    do j = 1, Ntime(i)
      xsi0 = Nisoreltalys(0, i, j)
      xsi1 = min(Nisoreltalys(1, i, j), 10. * xsi0)
      xsdifi = xsi0 - xsi1
      term0 = xsdifi * xsdifi
      if (abs(term0) > coveps .and. Sws0 > 0.) then
        xs2 = xsi0 * xsi0
        if (xs2 > 0.) then
          term = term0 / xs2
          RNisorelD(i, j) = (RNisorelD(i, j) * Swp0 + Sw0 * term) / Sws0
          errNisorel(i, j) = sqrt(RNisorelD(i, j))
        endif
      endif
      if (Sws0 > 0.) Nisorelav(i, j) = (Nisorelav(i, j) * Swp0 + Sw0 * xsi1) / Sws0
    enddo
  enddo
!
! ************** Output of average isotope production ******************
!
! Average isotope production
!
  do i = 1, NchanY
    ofile = 'act000000.xxx        '
    ofile(4:13) = Yfile(i)(2:11)
    if (Yfile(i)(9:9) == 'L') then
      write(ofile(14:17), '(".ave")')
    else
      write(ofile(11:13), '("ave")')
    endif
    open (unit = 3, file = ofile, status = 'replace')
    write(3, '("#", a20, " Number of runs: ", i6)') Yfile(i), italys
    if (flagband) then
      write(3, '("# Time [h] Activity [Mbq]  Upper [Mbq]", "   Lower [Mbq]")')
    else
      write(3, '("# Time [h] Activity [Mbq]  Error [Mbq]       ")')
    endif
    do j = 1, Ntime(i)
      err = actav(i, j) * erract(i, j)
      if (err /= 0.) err = min(err, actav(i, j) - 1.e-13)
      if (flagband) then
        write(3, '(4es12.5)') timeY(i, j), actav(i, j), actav(i, j)-err, actav(i, j) + err
      else
        write(3, '(3es12.5)') timeY(i, j), actav(i, j), err
      endif
    enddo
    close (3)
  enddo
  do i = 1, NchanY
    ofile = 'Niso000000.xxx       '
    ofile(5:14) = Yfile(i)(2:11)
    if (Yfile(i)(9:9) == 'L') then
      write(ofile(15:18), '(".ave")')
    else
      write(ofile(12:14), '("ave")')
    endif
    open (unit = 3, file = ofile, status = 'replace')
    write(3, '("#", a20, " Number of runs: ", i6)') Yfile(i), italys
    if (flagband) then
      write(3, '("# Time [h] # isotopes      Upper      ", "   Lower      ")')
    else
      write(3, '("# Time [h] # isotopes      Error             ")')
    endif
    do j = 1, Ntime(i)
      err = Nisoav(i, j) * errNiso(i, j)
      if (err /= 0.) err = min(err, Nisoav(i, j) - 1.e-13)
      if (flagband) then
        write(3, '(4es12.5)') timeY(i, j), Nisoav(i, j), Nisoav(i, j)- err, Nisoav(i, j) + err
      else
        write(3, '(3es12.5)') timeY(i, j), Nisoav(i, j), err
      endif
    enddo
    close (3)
  enddo
  do i = 1, NchanY
    ofile = 'yield000000.xxx      '
    ofile(6:15) = Yfile(i)(2:11)
    if (Yfile(i)(9:9) == 'L') then
      write(ofile(16:19), '(".ave")')
    else
      write(ofile(13:15), '("ave")')
    endif
    open (unit = 3, file = ofile, status = 'replace')
    write(3, '("#", a20, " Number of runs: ", i6)') Yfile(i), italys
    if (flagband) then
      write(3, '("# Time [h] Yield [Mbq/mAh] Upper      ", "   Lower      ")')
    else
      write(3, '("# Time [h] Yield [Mbq/mAh] Error             ")')
    endif
    do j = 1, Ntime(i)
      err = yieldav(i, j) * erryield(i, j)
      if (err /= 0.) err = min(err, yieldav(i, j) - 1.e-13)
      if (flagband) then
        write(3, '(4es12.5)') timeY(i, j), yieldav(i, j), yieldav(i, j) - err, yieldav(i, j) + err
      else
        write(3, '(3es12.5)') timeY(i, j), yieldav(i, j), err
      endif
    enddo
    close (3)
  enddo
  do i = 1, NchanY
    ofile = 'Nisorel000000.xxx    '
    ofile(8:17) = Yfile(i)(2:11)
    if (Yfile(i)(9:9) == 'L') then
      write(ofile(18:21), '(".ave")')
    else
      write(ofile(15:17), '("ave")')
    endif
    open (unit = 3, file = ofile, status = 'replace')
    write(3, '("#", a20, " Number of runs: ", i6)') Yfile(i), italys
    if (flagband) then
      write(3, '("# Time [h] Isotopic frac.  Upper      ", "   Lower      ")')
    else
      write(3, '("# Time [h] Isotopic frac.  Error             ")')
    endif
    do j = 1, Ntime(i)
      err = Nisorelav(i, j) * errNisorel(i, j)
      if (err /= 0.) err = min(err, Nisorelav(i, j) - 1.e-13)
      if (flagband) then
        write(3, '(4es12.5)') timeY(i, j), Nisorelav(i, j), Nisorelav(i, j) - err, Nisorelav(i, j) + err
      else
        write(3, '(3es12.5)') timeY(i, j), Nisorelav(i, j), err
      endif
    enddo
    close (3)
  enddo
  return
end subroutine covproduction
! Copyright A.J. Koning 2021
