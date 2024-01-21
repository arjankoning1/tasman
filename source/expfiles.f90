subroutine expfiles
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read experimental cross sections for all nuclides
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
! All global variables
!   numchanexp    ! maximum number of channels with experimental data
!   numnuc        ! maximum number of nuclides with experimental data
! Variables for writing TALYS input files
!   A             ! mass number of nucleus
!   Z             ! charge number of nucleus
! Constants
!   mainis        ! main isotope
!   particle      ! type of particle
! Variables for experimental data files
!   aamax         ! maximal A value for which exp. data are included in the search
!   aamin         ! minimal A value for which exp. data are included in the search
!   expfilemul    ! file with experimental data
!   Nchanmul      ! number of reaction channels per nuclide
!   Nnuc          ! number of nuclides for which experimental data exists
!   zzmax         ! maximal Z value for which exp. data are included in the search
!   zzmin         ! minimal Z value for which exp. data are included in the search
!
! *** Declaration of local data
!
  implicit none
  logical           :: flag      ! distribution
  logical           :: lexist    ! logical to determine existence
  character(len=6)  :: ext       ! filename extension
  character(len=20) :: xsf       ! filename
  integer           :: aa        ! mass number
  integer           :: i         ! counter
  integer           :: ia        ! mass number of nucleus
  integer           :: ibeg      ! index to mark begin of word
  integer           :: ichan     ! counter for channels
  integer           :: id        ! counter for deuterons
  integer           :: ih        ! hole number
  integer           :: ii        ! counter
  integer           :: in        ! counter for neutrons
  integer           :: ip        ! particle number
  integer           :: it        ! counter for tritons
  integer           :: type      ! particle type
  integer           :: zz        ! charge number
!
! ********** Determine reaction channels with experimental data ********
!
! Note that for a multi-nuclide search, files with experimental data need to have an extension denoting the Z,A of the nucleus,
! e.g. xs200000.exp026056.
!
  ii = 0
Loop1: do zz = zzmin, zzmax
    if (aamin == 0) aamin = max(mainis(zz) - 10, 1)
    if (aamax == 0) aamax = mainis(zz) + 10
    do aa = aamin, aamax
      ichan = 0
      flag = .false.
      ext = '      '
      write(ext(1:6), '(2i3.3)') zz, aa
      xsf = 'total.exp'//ext
      inquire (file = xsf, exist = lexist)
      if (lexist) then
        ichan = ichan + 1
        if ( .not. flag) then
          ii = ii + 1
          flag = .true.
        endif
        expfilemul(ii, ichan) = xsf
      endif
      xsf = 'elastic.exp'//ext
      inquire (file = xsf, exist = lexist)
      if (lexist) then
        ichan = ichan + 1
        if ( .not. flag) then
          ii = ii + 1
          flag = .true.
        endif
        expfilemul(ii, ichan) = xsf
      endif
      xsf = 'nonelastic.exp'//ext
      inquire (file = xsf, exist = lexist)
      if (lexist) then
        ichan = ichan + 1
        if ( .not. flag) then
          ii = ii + 1
          flag = .true.
        endif
        expfilemul(ii, ichan) = xsf
      endif
      xsf = 'fission.exp'//ext
      inquire (file = xsf, exist = lexist)
      if (lexist) then
        ichan = ichan + 1
        if ( .not. flag) then
          ii = ii + 1
          flag = .true.
        endif
        expfilemul(ii, ichan) = xsf
      endif
      do ia = 0, 1
        do ih = 0, 1
          do it = 0, 1
            do id = 0, 1
              do ip = 0, 1
                do in = 0, 3
                  xsf = 'xs000000.exp'//ext
                  write(xsf(3:3), '(i1)') in
                  write(xsf(4:4), '(i1)') ip
                  write(xsf(5:5), '(i1)') id
                  write(xsf(6:6), '(i1)') it
                  write(xsf(7:7), '(i1)') ih
                  write(xsf(8:8), '(i1)') ia
                  inquire (file = xsf, exist = lexist)
                  if (lexist) then
                    ichan = ichan + 1
                    if (ichan <= numchanexp) then
                      if ( .not. flag) then
                        ii = ii + 1
                        flag = .true.
                      endif
                      expfilemul(ii, ichan) = xsf
                      write( * , * ) xsf
                    endif
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      do type = 1, 6
        ibeg = 0
        if (type == 1) ibeg = 1
        do i = ibeg, 40
          xsf = 'nn.L00exp'//ext
          write(xsf(2:2), '(a1)') particle(type)
          write(xsf(5:6), '(i2.2)') i
          inquire (file = xsf, exist = lexist)
          if (lexist) then
            ichan = ichan + 1
            if (ichan <= numchanexp) then
              if ( .not. flag) then
                ii = ii + 1
                flag = .true.
              endif
              expfilemul(ii, ichan) = xsf
            endif
          endif
        enddo
      enddo
      if (flag) then
        Z(ii) = zz
        A(ii) = aa
        Nchanmul(ii) = ichan
      endif
      if (ii == numnuc) exit Loop1
    enddo
  enddo Loop1
  Nnuc = ii
  return
end subroutine expfiles
! Copyright A.J. Koning 2021
