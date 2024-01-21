subroutine tmc
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Creation of random ENDF-6 files
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
! Variables for writing TALYS input
!   italys          ! TALYS run counter
!   Nhigh           ! number of high energy runs
! Constants
!   nuc             ! symbol of nucleus
! Variables for path names
!   binpath       ! directory containing files to be read
!   librariespath ! directory containing files to be read
! Variables for processing input
!   Atarget         ! mass number of target nucleus
!   k0              ! index of incident particle
!   Liso            ! number of the isomeric state
!   Liso            ! isomeric number of target
!   ptype0          ! type of incident particle
!   tafis           ! TAFIS executable
!   tanes           ! TANES executable
!   tares           ! TAREL executable
!   tefal           ! TEFAL executable
!   Ztarget         ! charge number of target nucleus
! Variables for covariances
!   weight          ! weight for experimental data sets
! Variables for GOF function
!   Ntalys          ! number of TALYS runs
! Variables for TMC
!   background      ! library name for adoption of back ground cross sections in MF3
!   fislim          ! mass above which nuclide fissions
!   flagacf         ! flag for creation of ENDF activation file
!   flagcovrand     ! flag for covariance ENDF - 6 file after every random run
!   flagdefonly     ! flag for creation of ENDF default file only
!   flageaf         ! flag for creation of EAF activation file
!   flagfns         ! flag for variation of fission neutron spectrum for ENDF - 6 general purpose file
!   flaggpf         ! flag for creation of ENDF general purpose file
!   flagmt          ! flag for creation of ENDF general purpose file without switch
!   flagnjoy        ! flag to run NJOY
!   flagnubar       ! flag for variation of nubar values for ENDF - 6 general purpose file
!   flagprepro      ! flag to run PREPRO codes
!   flagproconly    ! flag for processing of ENDF default file only
!   flagpurr        ! flag to run PURR in NJOY
!   flagres         ! flag for variation of resonance parameters for ENDF - 6 general purpose file
!   flagrunfns      ! flag to produce fns in TASMAN with TANES
!   flagrunnubar    ! flag to produce nubar in TASMAN with TAFIS
!   flagruntares    ! flag to produce resonances in TASMAN with TARES
!   flags20         ! flag for creation of ENDF general purpose file with switch at 20 MeV
!   flags30         ! flag for creation of ENDF general purpose file with switch at 30 MeV
!   flags60         ! flag for creation of ENDF general purpose file with switch at 60 MeV
!   flagsdefault    ! flag for creation of ENDF file with def. switch (30 MeV for neutrons, 0 for other particles)
!   flagselect      ! flag for varying only parts of an ENDF - 6 data file
!   offset          ! offset for numbering or random files (TMC only)
!   tafislib        ! library name for adoption of nubar values (none[default], endfb8.0, jeff3.3, jendl4.0)
!   taneslib        ! library name for adoption of FNS (none[default], endfb8.0, jeff3.3, jendl4.0)
!   tareslib        ! library name for adoption of resonance parameters (default or endfb8.0)
!   tmcoffset       ! offset for starting creation of ENDF - 6 files (TMC only)
! Variables for reading cross sections
!   E               ! incident energy
!   Nen             ! number of incident energies
!   Nenlow          ! number of incident energies
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=1)   :: backsl         ! backslash
  logical            :: lexist         ! logical to determine existence
  character(len=3)   :: mcnpext        ! extension for ACE file
  character(len=3)   :: Astring        ! string for mass number
  character(len=4)   :: randstring     ! string for random number
  character(len=8)   :: nuclide        ! nuclide symbol
  character(len=8)   :: purrflag       ! flag for PURR
  character(len=10)  :: projnuc        ! projectile-nuclide combination
  character(len=10)  :: addfile        ! file for TEFAL additional keywords
  character(len=15)  :: tefalfile      ! input file for TEFAL
  character(len=15)  :: tefaloutfile   ! output file for TEFAL
  character(len=17)  :: select         ! command for select code
  character(len=22)  :: endffile       ! ENDF file
  character(len=132) :: gpffile        ! file for GPF
  character(len=132) :: nubarfile      ! file with nubar
  character(len=132) :: nudelfile      ! file with nudel
  character(len=132) :: Efisfile       ! file with fission parameters
  character(len=132) :: fnsfile        ! file with nudel
  character(len=132) :: delfile        ! file with delayed FNS
  character(len=132) :: mf2file        ! file with MF2
  character(len=28)  :: pendffile      ! file for PENDF
  character(len=132) :: str            ! input line
  character(len=132) :: base           ! help variable
  character(len=132) :: beststr(100)   ! string for best keywords
  character(len=132) :: tefalstr(20)   ! string with info for TEFAL
  character(len=132) :: inpfile        ! input file
  character(len=132) :: outfile        ! output file
  character(len=160) :: cmd            ! command
  integer            :: Afis           ! mass value for fission
  integer            :: i              ! counter
  integer            :: ibest          ! counter for best keywords
  integer            :: istat          ! logical for file access
  integer            :: ital           ! TALYS run counter
  integer            :: k              ! Legendre order
  integer            :: NE             ! number of incident energies
  integer            :: NEhigh         ! number of energies in high energy region
  integer            :: Nstr           ! number of lines
  integer            :: system         ! system call
  integer            :: Zfis           ! charge value for fission
!
! Jump to system for possible ENDF-6 file generation
!
! Resonance data
!
!   ENDF-6 general purpose file.
!
  backsl = '\'
  Astring = '000'
  write(Astring, '(i3.3)') Atarget
  nuclide = trim(nuc(Ztarget)) //Astring
  if (Liso > 0) then
    if (Liso == 1) nuclide = trim(nuclide)//'m'
    if (Liso >= 2) nuclide = trim(nuclide)//'n'
  endif
  projnuc = ptype0//'-'//trim(nuclide)
  ital = italys + offset
  if (italys > tmcoffset) ital = ital - tmcoffset
  randstring = '0000'
  write(randstring, '(i4.4)') ital
  write(*, '(" TASMAN: Create random ENDF-6 file ")')
  if (flagres .and. flaggpf) then
    if (italys == 0 .and. flagruntares) then
      inpfile = 'tares.inp'
      outfile = 'tares.out'
      open (unit = 3, file = trim(inpfile), status = 'replace')
      write(3, '("charge ", i3)') Ztarget
      write(3, '("mass ", i3)') Atarget
      write(3, '("isomer ", i3)') Liso
      write(3, '("output ", a)') trim(projnuc)
      if (Atarget <= 19) write(3, '("urr n")')
      write(3, '("mf32mf33 y")')
      write(3, '("nbrRandom ", i6)') Ntalys
      if (tareslib(1:7) /= 'default') write(3, '("#library ", a16)') tareslib
      close (unit = 3)
      write(*, '(" TASMAN: Run TARES ")')
      cmd = trim(tares) //" " //inpfile //" > " //outfile
      i = system(cmd)
    endif
  endif
!
! nubar data
!
  if ((flagrunnubar .or. flagrunfns) .and. flaggpf) then
    open (unit = 1, file = 'nubarfns', status = 'replace')
    write(1, '("#!/bin/bash")')
  endif
  if (flagrunnubar .and. flaggpf) then
    inpfile = 'tafis.'//randstring//'.inp'
    outfile = 'tafis.'//randstring//'.out'
    Zfis = max(Ztarget, 84)
    Afis = max(Atarget, fislim)
    open (unit = 3, file = trim(inpfile), status = 'replace')
    write(3, '("charge ", i3)') Zfis
    write(3, '("mass ", i3)') Afis
    write(3, '("projectile n")')
    write(3, '("normalization ", a16)') tafislib
    if (italys == 0) then
      write(3, '("random 0")')
    else
      write(3, '("random 1")')
    endif
    close (unit = 3)
    write(1, '("cp -f ", a, " tafis.inp")') trim(inpfile)
    write(1, '(a, " tafis.inp > tafis.out ")') trim(tafis)
    write(1, '("cp -f tafis.out ", a)') trim(outfile)
  endif
!
! Fission neutron spectrum
!
  if (flagrunfns .and. flaggpf) then
    inpfile = 'tanes.'//randstring//'.inp'
    outfile = 'tanes.'//randstring//'.out'
    Zfis = max(Ztarget, 84)
    Afis = max(Atarget, fislim)
    open (unit = 3, file = trim(inpfile), status = 'replace')
    write(3, '("charge ", i3)') Zfis
    write(3, '("mass ", i3)') Afis
    write(3, '("library ", a16)') taneslib
    if (italys == 0) then
      write(3, '("random 0")')
    else
      write(3, '("random 1")')
    endif
    close (unit = 3)
    write(1, '("cp -f ", a, " tanes.inp")') trim(inpfile)
    write(1, '(a, " tanes.inp > tanes.out ")') trim(tanes)
    write(1, '("cp -f tanes.out ", a)') trim(outfile)
  endif
  if ((flagrunnubar .or. flagrunfns) .and. flaggpf) then
    close (unit = 1)
    write(*, '(" TASMAN: Run nubarfns ")')
    i = system('chmod a+x nubarfns; ./nubarfns')
  endif
!
! Prepare TEFAL run
!
  tefalfile = '   0000.inp    '
  tefaloutfile = '   0000.out    '
  write(tefalfile(4:7), '(i4.4)') ital
  write(tefaloutfile(4:7), '(i4.4)') ital
!
! Overwrite file with incident energies in case of avoiding high-energy random calculations
!
  if (italys > Nhigh) then
    NEhigh = Nen(1) - Nenlow(1)
    inquire (file = 'tefal.inf', exist = lexist)
    if (lexist) then
      open (unit = 1, file = 'tefal.inf', status = 'old')
      i = 1
      do
        read(1, '(a)', iostat = istat) tefalstr(i)
        if (istat == -1) exit
        if (istat /= 0) call read_error('tefal.inf', istat)
        i = i + 1
      enddo
      close (unit = 1)
      Nstr = i - 1
      read(tefalstr(8), '(i4)') NE
      NE = NE + NEhigh
      write(tefalstr(8)(1:4), '(i4)') NE
      open (unit = 1, file = 'tefal.inf', status = 'replace')
      do i = 1, Nstr
        write(1, '(a)') trim(tefalstr(i))
      enddo
      close (unit = 1)
      open (unit = 1, file = 'energies.endf', status = 'old', position = 'append')
      do i = Nenlow(1) + 1, Nen(1)
        write(1, '(es12.5)') E(1, i)
      enddo
      close (unit = 1)
    endif
  endif
!
! If requested, perform 7 TEFAL runs: for general purpose file, activation file and EAF file
!
  do k = 1, 7
    if (k == 1 .and. .not. flagsdefault) cycle
    if (k > 1 .and. flagdefonly) cycle
    if (k == 2 .and. .not. flags20) cycle
    if (k == 3 .and. .not. flags30) cycle
    if (k == 4 .and. .not. flags60) cycle
    if (k == 5 .and. .not. flagmt) cycle
    if (k == 6 .and. .not. flagacf) cycle
    if (k == 7 .and. .not. flageaf) cycle
    if (k <= 5) tefalfile(1:3) = 'gpf'
    tefalfile(12:15) = '    '
    if (k == 2) tefalfile(12:15) = '.s20'
    if (k == 3) tefalfile(12:15) = '.s20'
    if (k == 4) tefalfile(12:15) = '.s60'
    if (k == 5) tefalfile(12:15) = '.MT '
    if (k == 6) tefalfile(1:3) = 'acf'
    if (k == 7) tefalfile(1:3) = 'eaf'
    tefaloutfile(1:3) = tefalfile(1:3)
    ibest = -1
    addfile = tefalfile(1:3)//'.inp'
    open (unit = 1, file = tefalfile, status = 'replace')
    inquire (file = addfile, exist = lexist)
    if (lexist) then
      open (unit = 2, file = addfile, status = 'old', iostat = istat)
      if (istat /= 0) call read_error(addfile, istat)
      do
        read(2, '(a)', iostat = istat) str
        if (istat /= 0) exit
        if (str(1:8) == 'endffile') cycle
        if (str(3:25) == 'Specific keywords added') ibest = 0
        if (ibest >= 0) then
          ibest = ibest + 1
          beststr(ibest) = str
        endif
        write(1, '(a)') trim(str)
      enddo
      close (unit = 2)
    endif
    gpffile = trim(projnuc)//'-rand-'//randstring
    if (k == 1) then
      if (k0 == 1) then
        write(1, '("Eswitch 30.")')
      else
        write(1, '("Eswitch 0.")')
      endif
    endif
    if (k == 2) then
      gpffile = trim(gpffile)//'.s20'
      write(1, '("Eswitch 20.")')
    endif
    if (k == 3) then
      gpffile = trim(gpffile)//'.s30'
      write(1, '("Eswitch 30.")')
    endif
    if (k == 4) then
      gpffile = trim(gpffile)//'.s60'
      write(1, '("Eswitch 60.")')
      if (k0 == 3) write(1, '("breakup y")')
    endif
    if (k == 5) then
      gpffile = trim(gpffile)//'.MT'
      if (k0 == 3) then
        write(1, '("Eswitch 0.")')
        write(1, '("breakup y")')
      else
        write(1, '("Eswitch 250.")')
        write(1, '("mtall n")')
      endif
    endif
    if (k <= 5) then
      endffile = trim(gpffile)
    endif
    if (k == 6) then
      endffile = trim(gpffile)//'.acf'
      write(1, '("gpf n")')
    endif
    if (k == 7) then
      endffile = trim(gpffile)//'.eaf'
      write(1, '("gpf n")')
      write(1, '("eaf y")')
    endif
    write(1, '("endffile ", a)') trim(endffile)
    write(1, '("diffweight ", es12.5)') weight(0)
    if (k <= 5) then
      if (flagres) then
        if (flagruntares) then
          mf2file = 'other.tares.'//trim(projnuc)//'.files/randomMF2/' //trim(projnuc)//'.mf2.'//randstring
        else
          base = trim(librariespath)//'resbase/'// trim(nuclide)//'/random/'//projnuc
          mf2file = trim(base)//'.mf2.'//randstring
        endif
        inquire (file = trim(mf2file), exist = lexist)
        if (.not.lexist) then
          write(*,'(" TASMAN-error: ",a," can not be found: ")') trim(mf2file)
          stop
        endif
        write(1, '("adopt 2 151 ", a)') trim(mf2file)
      endif
      if (background(1:1) /= ' ') write(1, '("background ", a)') trim(background)
      if (flagnubar) then
        if (flagrunnubar) then
          if (italys == 0) then
            nubarfile = 'fission.mt452.0000'
            nudelfile = 'fission.mt455.0000'
            Efisfile = 'fission.mt458.0000'
          else
            nubarfile = 'fission.mt452.0001'
            nudelfile = 'fission.mt455.0001'
            Efisfile = 'fission.mt458.0001'
          endif
        else
          base = trim(librariespath)//'nubarbase/'// trim(nuclide)//'/random/'//projnuc
          if (italys == 0) then
            nubarfile = trim(base)//'.mt452.0000'
            nudelfile = trim(base)//'.mt455.0000'
            Efisfile = trim(base)//'.mt458.0000'
          else
            nubarfile = trim(base)//'.mt452.'//randstring
            nudelfile = trim(base)//'.mt455.'//randstring
            Efisfile = trim(base)//'.mt458.'//randstring
          endif
        endif
        inquire (file = trim(nubarfile), exist = lexist)
        if (.not.lexist) then
          write(*,'(" TASMAN-error: ",a," can not be found: ")') trim(nubarfile)
          stop
        endif
        inquire (file = trim(nudelfile), exist = lexist)
        if (.not.lexist) then
          write(*,'(" TASMAN-error: ",a," can not be found: ")') trim(nudelfile)
          stop
        endif
        inquire (file = trim(Efisfile), exist = lexist)
        if (.not.lexist) then
          write(*,'(" TASMAN-error: ",a," can not be found: ")') trim(Efisfile)
          stop
        endif
        write(1, '("adopt 1 452 ", a)') trim(nubarfile)
        write(1, '("adopt 1 455 ", a)') trim(nudelfile)
        write(1, '("adopt 1 458 ", a)') trim(Efisfile)
      endif
      if (flagfns) then
        if (flagrunfns) then
          if (italys == 0) then
            fnsfile = 'fns.mf5.0000      '
            delfile = 'fns.mf5.delayed.0000      '
          else
            fnsfile = 'fns.mf5.0001      '
            delfile = 'fns.mf5.delayed.0001      '
          endif
        else
          base = trim(librariespath)//'fnsbase/'// trim(nuclide)//'/random/'//projnuc
          if (italys == 0) then
            fnsfile = trim(base)//'.mf5.0000'
            delfile = trim(base)//'.delayed.mf5.0000'
          else
            fnsfile = trim(base)//'.mf5.'//randstring
            delfile = trim(base)//'.delayed.mf5.'//randstring
          endif
        endif
        inquire (file = trim(fnsfile), exist = lexist)
        if (.not.lexist) then
          write(*,'(" TASMAN-error: ",a," can not be found: ")') trim(fnsfile)
          stop
        endif
        inquire (file = trim(delfile), exist = lexist)
        if (.not.lexist) then
          write(*,'(" TASMAN-error: ",a," can not be found: ")') trim(delfile)
          stop
        endif
        write(1, '("adopt 5 18 ", a)') trim(fnsfile)
        write(1, '("adopt 5 455 ", a)') trim(delfile)
      endif
    endif
!
! Create pointwise cross sections for EAF file
!
    if (k == 7) then
      if (flaggpf .and. flagprepro) then
        cmd = trim(binpath)//'autoprepro -file '//trim(gpffile)
        write(*, '(" TASMAN: Run PREPRO ")')
        i = system(cmd)
        pendffile = trim(gpffile)//'.pendf'
        write(1, '("adopt 3 102 ", a, " 1.e-5 2.e8")') trim(pendffile)
        inquire (file = 'fission.tot', exist = lexist)
        if (lexist) write(1, '("adopt 3 18 ", a, " 1.e-5 2.e8")') trim(pendffile)
      endif
    endif
    do i = 1, ibest
      write(1, '(a)') trim(beststr(i))
    enddo
    close (unit = 1)
    write(*, '(" TASMAN: Run TEFAL ")')
    i = system(tefal//'<'//tefalfile//'>'//tefaloutfile)
!
! Run NJOY
!
    if (k <= 5 .and. flaggpf .and. flagnjoy) then
      if (k == 1 .or. .not. flagproconly) then
        if (flagpurr) then
          purrflag = ' -purr'
        else
          purrflag = ' -nopurr'
        endif
        mcnpext = '.00'
        if (Liso == 1) mcnpext = '.10'
        if (Liso >= 2) mcnpext = '.20'
        cmd = trim(binpath)//'autonjoy -file '//trim(gpffile)// purrflag//' -outname 2 -mcnpext '//mcnpext
        write(*, '(" TASMAN: Run NJOY ")')
        write(*, '(a)') trim(cmd)
        i = system(cmd)
      endif
    endif
  enddo
!
! Create temporary directory with central values for the 'select' option to vary only certain parts of an ENDF-6 data file.
!
! For the TALYS nuclear data evaluation system, often a very large number of files need to be removed, copied or moved.
! A standard Linux operating system may not be able to handle this, which results in the message "Argument list too long".
! We circumvent this with the 'find' command.
!
  if ((flagselect .or. flagcovrand) .and. flaggpf) then
    if (italys == 0) then
      i = system('mkdir 0000')
      cmd = 'find . -maxdepth 1 -name "*" -exec cp -f {} 0000 '// backsl//';'
      i = system(cmd)
      i = system('cp -f nonelastic.000 0000/nonelastic.tot')
    else
      if (flagres) then
        open (unit = 3, file = 'taresela.inp', status = 'replace')
        write(3, '("charge ", i3)') Ztarget
        write(3, '("mass ", i3)') Atarget
        write(3, '("isomer ", i3)') Liso
        write(3, '("output mf2.ela")')
        if (Atarget <= 19) write(3, '("urr n")')
        write(3, '("random GnR")')
        if (tareslib(1:7) /= 'default') write(3, '("#library ", a16)') tareslib
        close (unit = 3)
        cmd = tares //" taresela.inp > taresela.out"
        i = system(cmd)
        open (unit = 3, file = 'tarescap.inp', status = 'replace')
        write(3, '("charge ", i3)') Ztarget
        write(3, '("mass ", i3)') Atarget
        write(3, '("isomer ", i3)') Liso
        write(3, '("output mf2.cap")')
        if (Atarget <= 19) write(3, '("urr n")')
        write(3, '("random Gg")')
        if (tareslib(1:7) /= 'default') write(3, '("#library ", a16)') tareslib
        close (unit = 3)
        cmd = tares //" tarescap.inp > tarescap.out"
        i = system(cmd)
        open (unit = 3, file = 'taresfis.inp', status = 'replace')
        write(3, '("charge ", i3)') Ztarget
        write(3, '("mass ", i3)') Atarget
        write(3, '("isomer ", i3)') Liso
        write(3, '("output mf2.fis")')
        if (Atarget <= 19) write(3, '("urr n")')
        write(3, '("random Gf")')
        if (tareslib(1:7) /= 'default') write(3, '("#library ", a16)') tareslib
        close (unit = 3)
        cmd = tares //" taresfis.inp > taresfis.out"
        i = system(cmd)
      endif
      select = '         0000 0 0'
      write(select(1:8), '(a8)') projnuc
      write(select(10:13), '(i4.4)') ital
      if (flagselect) write(select(15:15), '("1")')
      if (flagcovrand) write(select(17:17), '("1")')
      cmd = trim(tasmanpath)//'misc/select '//select
      write(*, '(" TASMAN: Run select: ", a)') trim(cmd)
      i = system(cmd)
    endif
  endif
  return
end subroutine tmc
! Copyright A.J. Koning 2021
