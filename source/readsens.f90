subroutine readsens
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Read sensitivity matrix
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
!   sgl          ! single precision kind
! All global variables
!   numenin      ! maximum number of incident energies
! Variables for path names
!   tasmanpath   ! directory containing files to be read
! Variables for processing input
!   Atarget      ! mass number of target nucleus
!   isochar      ! symbol for isomer
!   Liso         ! number of the isomeric state
!   Starget      ! symbol of target nucleus
! Variables for cross section covariances
!   S            ! sensitivity matrix
! Variables for parameter covariances
!   parstring    ! parameter string
! Variables for reading cross sections
!   E            ! incident energy
!   Nchanxs      ! total number of channels with cross sections
!   Nen          ! number of incident energies
!   Sindex       ! index for cross section
!   xsfile       ! file with crosssections
! Error handling
!   read_error ! Message for file reading error
!
! *** Declaration of local data
!
  implicit none
  integer            :: i                 ! counter
  integer            :: ii                ! counter
  integer            :: istat             ! logical for file access
  integer            :: j                 ! counter
  integer            :: k                 ! Legendre order
  integer            :: kk                ! counter
  integer            :: l                 ! counter
  integer            :: n                 ! counter
  integer            :: Nchan             ! number of channels
  integer            :: NE                ! number of incident energies
  integer            :: nn                ! help variable
  integer            :: Np                ! number of points
  real(sgl)          :: e1                ! first energy
  real(sgl)          :: e2                ! second energy
  real(sgl)          :: ee(0:numenin)     ! incident energy
  real(sgl)          :: Einc              ! incident energy in MeV
  real(sgl)          :: fac               ! factor
  real(sgl)          :: s1                ! first sensitivity value
  real(sgl)          :: s2                ! second sensitivity value
  real(sgl)          :: sint              ! interpolated sensitivity value
  real(sgl)          :: Sval(0:numenin)   ! sensitivity value
  character(len=3)   :: massstring        ! string for mass number
  character(len=20)  :: chanfile          ! file with channel cross section
  character(len=30)  :: sensstring        ! string for sensitivities
  character(len=132) :: sensfile          ! file with sensitivities
!
! *************************** Read S matrix ****************************
!
! Linear sensitivity matrix
!
  massstring = '000'
  write(massstring(1:3), '(i3.3)') Atarget
  if (Liso == 0) then
    sensfile = trim(tasmanpath)//'sens/'// trim(Starget)//massstring//'/sensitivity.par'
  else
    sensfile = trim(tasmanpath)//'sens/'// trim(Starget)//massstring//isochar//'/sensitivity.par'
  endif
  open (unit = 3, file = sensfile, status = 'old', iostat = istat)
  if (istat /= 0) call read_error(sensfile, istat)
  read(3, '(23x, i4)', iostat = istat) Np
  if (istat /= 0) call read_error(sensfile, istat)
  do kk = 1, Np
    read(3, '(24x, a30)', iostat = istat) sensstring
    if (istat /= 0) call read_error(sensfile, istat)
    k = 0
    do l = 1, Np
      if (parstring(l) == sensstring) then
        k = l
        exit
      endif
    enddo
    read(3, '(/23x, i4)', iostat = istat) Nchan
    if (istat /= 0) call read_error(sensfile, istat)
    do ii = 1, Nchan
      read(3, '(24x, a20)', iostat = istat) chanfile
      if (istat /= 0) call read_error(sensfile, istat)
      i = 0
      do j = 1, Nchanxs
        if (chanfile == xsfile(j)) then
          i = j
          exit
        endif
      enddo
      read(3, '(24x, i4/)', iostat = istat) NE
      if (istat /= 0) call read_error(sensfile, istat)
      do n = 1, NE
        read(3, '(2e15.5)', iostat = istat) ee(n), Sval(n)
        if (istat /= 0) call read_error(sensfile, istat)
      enddo
      ee(0) = 0.
      Sval(0) = Sval(1)
      if (i == 0) cycle
!
! Interpolate on current energy grid
!
      do j = 1, Nen(i)
        Einc = E(i, j)
        n = Sindex(i, j)
        if (Einc <= ee(NE)) then
          call locate(ee, 0, NE, Einc, nn)
          e1 = ee(nn)
          e2 = ee(nn + 1)
          fac = (Einc - e1) / (e2 - e1)
          s1 = Sval(nn)
          s2 = Sval(nn + 1)
          sint = s1 + fac * (s2 - s1)
        else
          nn = Nen(i)
          sint = Sval(nn)
        endif
        S(k, i, n) = sint
      enddo
    enddo
  enddo
  close (3)
  return
end subroutine readsens
! Copyright A.J. Koning 2021
