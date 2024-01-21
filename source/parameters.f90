   subroutine parameters
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Vary all tabulated parameters or a subset
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
! All global variables
!   numlines    ! maximum number of input lines
!   numpar      ! maximum number of parameters
! Variables for reading TASMAN input lines
!   inline      ! input line
!   nlines      ! number of input lines
! Variables for processing input
!   talys  ! TALYS executable
! Variables for automatic parameter variation
!   flagallpar     ! flag to vary all nuclear model parameters
!   maxbar         ! maximal fission barrier
!   minbar         ! minimal fission barrier
!   partvary       ! particle for varying parameters
!   keyvary        ! keyword to vary
!   ACN            ! A of compound nucleus
!   Adeep          ! maximum depth of A
!   ZCN            ! Z of compound nucleus
!   Zdeep          ! maximum depth of Z
!   Nkeyvary       ! number of keywords to vary
!   Npartvary      ! number of particles to vary keywords
!   NZAskip        ! number of Z,A pairs to skip for parameter variation
!   NZAvary        ! number of Z,A pairs for parameter variation
!   Askip          ! A value to skip
!   Zskip          ! Z value to skip
!   Avary          ! A value to vary
!   Zvary          ! Z value to vary
!   sigvary        ! sigma for parameter variation
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  character(len=3)   :: partstring          !
  character(len=8)   :: ZAstring            !
  character(len=10)  :: sigstring           !
  character(len=132) :: talcmd              ! TALYS command
  character(len=132) :: word(40)            ! words on input line
  character(len=132) :: line                ! input line
  character(len=132) :: parline(numlines)   !
  integer            :: Aend                ! last A to be included
  integer            :: i                   ! counter
  integer            :: istat               ! logical for file access
  integer            :: ia                  ! mass number of nucleus
  integer            :: iz                  ! charge number of nucleus
  integer            :: j                   ! counter
  integer            :: k                   ! counter
  integer            :: m                   ! counter
  integer            :: Mvary               !
  integer            :: ZZ                  !
  integer            :: AA                  !
  integer            :: N                   ! number of parameters
  integer            :: parlines            !
  integer            :: system              ! system call
  integer            :: Zend                ! maximal charge number
!
! **************** Include nuclides according to search depth **********
!
  if (Zdeep /= -1 .and. Adeep /= -1) then
    Zend = ZCN - Zdeep
    Aend = ACN - Adeep
    N = NZAvary
    do iz = Zend, ZCN
    ialoop: do ia = Aend, ACN
        do k = 1, NZAskip
          ZZ = ZCN - Zskip(k)
          AA = ACN - Askip(k)
          if (ZZ == iz .and. AA == ia) cycle ialoop
        enddo
        N = N + 1
        if (N > numpar) cycle
        Zvary(N) = iz
        Avary(N) = ia
      enddo ialoop
    enddo
    NZAvary = N
  endif
  NZAvary = min(NZAvary, numpar)
!
! ************** Run TALYS and read parameter file *********************
!
  write(*, '(" Performing first TALYS run to retrieve parameters")')
  open (unit = 1, file = 'talys.inp', status = 'unknown')
  do i = 1, nlines
    write(1, '(a)') trim(inline(i))
  enddo
  write(1, '("partable y")')
  write(1, '("spherical y")')
  close(1)
  talcmd = trim(talys)//' < talys.inp > talys.out'
  i = system(trim(talcmd))
  open (unit = 2, file = 'parameters.dat', status = 'old', iostat = istat)
  if (istat /= 0) then
    write(*, '(" TASMAN-error: Expected parameters.dat")')
    stop
  endif  
  i = 1
  do
    read(2, '(a)', iostat = istat) parline(i)
    if (istat ==  -1) exit
    i = i + 1
    if (i > numlines - nlines) then
      write(*, '(" TASMAN-error: Number of input lines exceeds ", i5)') numlines
      write(*, '(" numlines in tasman.cmb should be increased")')
      stop
    endif
  enddo
  parlines = i - 1
  close(2)
  if (parlines == 0) write(*, '(" TASMAN-warning: parameters.dat is empty")')
  N = nlines
  Mvary = 0
  lines: do i = 1, parlines
    if ( .not. flagallpar) then
      line = parline(i)
      call convert(line)
      parline(i) = line
      call getkeywords(parline(i), word)
      do k = 1, NZAvary
        ZAstring = '        '
        write(ZAstring(1:4), '(i4)') Zvary(k)
        write(ZAstring(5:8), '(i4)') Avary(k)
        j = index(parline(i), ZAstring)
        if (j > 0) then
          m = index(parline(i), ' 0 ')
          if (m > 0 .and. minbar > 0) cycle
          m = index(parline(i), ' 1 ')
          if (m > 0 .and. (minbar > 1 .or. maxbar < 1)) cycle
          m = index(parline(i), ' 2 ')
          if (m > 0 .and. (minbar > 2 .or. maxbar < 2)) cycle
          m = index(parline(i), ' 3 ')
          if (m > 0 .and. (minbar > 3 .or. maxbar < 3)) cycle
          m = index(parline(i), ' m1 ')
          if (m > 0 .and. .not. flagM1vary) cycle
          m = index(parline(i), ' e1 ')
          if (m > 0 .and. .not. flagE1vary) cycle
          do m = 1, Nkeyvary
            if (trim(word(1)) == trim(keyvary(m))) then
              Mvary = m
              goto 160
            endif
          enddo
        endif
      enddo
      do k = 1, Npartvary
        partstring = '   '
        write(partstring(2:2), '(a1)') partvary(k)
        j = index(parline(i), partstring)
        if (j > 0) then
          do m = 1, Nkeyvary
            if (trim(word(1)) == trim(keyvary(m))) then
              Mvary = m
              goto 160
            endif
          enddo
        endif
      enddo
      if (word(3)(1:20) == '                    ') then
        do m = 1, Nkeyvary
          if (trim(word(1)) == trim(keyvary(m))) then
            Mvary = m
            goto 160
          endif
        enddo
      endif
      cycle
    endif
  160   N = N + 1
    sigstring = '          '
    if (Mvary > 0) then
      if (sigvary(Mvary) > 0.) write(sigstring, '(f10.5)') sigvary(Mvary)
    endif
    inline(N) = trim(parline(i)) //trim(sigstring)
  enddo lines
  nlines = N
  return
end subroutine parameters
! Copyright A.J. Koning 2021
