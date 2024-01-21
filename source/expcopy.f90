subroutine expcopy
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Remove and copy experimental data files for multi-nuclide
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
! Variables for writing TALYS input files
!   A           ! mass number of nucleus
!   inuc        ! index of nuclide
!   Z           ! charge number of nucleus
! Constants
!   particle    ! type of particle
!
! *** Declaration of local data
!
  implicit none
  logical           :: lexist    ! logical to determine existence
  character(len=6)  :: ext       ! filename extension
  character(len=20) :: xsf0      ! experimental data file
  character(len=40) :: cmd       ! command
  integer           :: i         ! counter
  integer           :: ia        ! mass number of nucleus
  integer           :: ibeg      ! index to mark begin of word
  integer           :: id        ! counter for deuterons
  integer           :: ih        ! hole number
  integer           :: in        ! counter for neutrons
  integer           :: ip        ! particle number
  integer           :: it        ! counter for tritons
  integer           :: system    ! system call
  integer           :: type      ! particle type
!
! *************** Remove and copy experimental data files **************
!
! For every new nuclide in the multi-nuclide search, the appropriate experimental data files are copied into the .exp files
!
  i = system('rm *.exp')
  ext = '      '
  write(ext(1:6), '(2i3.3)') Z(inuc), A(inuc)
  xsf0 = 'total.exp'//ext
  inquire (file = xsf0, exist = lexist)
  if (lexist) then
    cmd = 'cp '//xsf0//' total.exp'
    i = system(cmd)
  endif
  xsf0 = 'elastic.exp'//ext
  inquire (file = xsf0, exist = lexist)
  if (lexist) then
    cmd = 'cp '//xsf0//' elastic.exp'
    i = system(cmd)
  endif
  xsf0 = 'nonelastic.exp'//ext
  inquire (file = xsf0, exist = lexist)
  if (lexist) then
    cmd = 'cp '//xsf0//' nonelastic.exp'
    i = system(cmd)
  endif
  xsf0 = 'fission.exp'//ext
  inquire (file = xsf0, exist = lexist)
  if (lexist) then
    cmd = 'cp '//xsf0//' fission.exp'
    i = system(cmd)
  endif
  do ia = 0, 1
    do ih = 0, 1
      do it = 0, 1
        do id = 0, 1
          do ip = 0, 1
            do in = 0, 3
              xsf0 = 'xs000000.exp'//ext
              write(xsf0(3:3), '(i1)') in
              write(xsf0(4:4), '(i1)') ip
              write(xsf0(5:5), '(i1)') id
              write(xsf0(6:6), '(i1)') it
              write(xsf0(7:7), '(i1)') ih
              write(xsf0(8:8), '(i1)') ia
              inquire (file = xsf0, exist = lexist)
              if (lexist) then
                cmd = 'cp '//xsf0//' '//xsf0(1:12)
                i = system(cmd)
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
      xsf0 = 'nn.L00exp'//ext
      write(xsf0(2:2), '(a1)') particle(type)
      write(xsf0(5:6), '(i2.2)') i
      inquire (file = xsf0, exist = lexist)
      if (lexist) then
        cmd = 'cp '//xsf0//' '//xsf0(1:9)
        ia = system(cmd)
      endif
    enddo
  enddo
  return
end subroutine expcopy
! Copyright A.J. Koning 2021
