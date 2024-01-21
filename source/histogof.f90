subroutine histogof
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make histograms of goodness-of-fit values
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
!   numhist         ! maximum number of points for histogram
!   numtalys        ! maximum number of TALYS runs
! Variables for reading TALYS output
!   Nchanall        ! total number of channels
! Variables for GOF function
!   isearch         ! number of trial run
! Variables for reading cross sections
!   xsfile          ! file with crosssections
! Variables for weights
!   Gchannelsave    ! GOF for channel
!   Gsave           ! GOF value
!   weightsave      ! weight for TALYS run
!
! *** Declaration of local data
!
  implicit none
  character(len=24) :: goffile               ! file with GOF results
  integer           :: gofhist(numhist)      ! GOF histogram value
  integer           :: i                     ! counter
  integer           :: iloc(0:numtalys)      ! local help variable
  integer           :: itmp                  ! help variable
  integer           :: j                     ! counter
  integer           :: k                     ! Legendre order
  integer           :: n                     ! counter
  integer           :: Nbeg                  ! begin of energy loop
  integer           :: Nbin                  ! number of bin
  integer           :: Nlast                 ! last discrete level
  real(sgl)         :: dbin                  ! width of bin
  real(sgl)         :: Gav(0:numtalys)       ! average GOF value
  real(sgl)         :: Gavdev(0:numtalys)    ! average GOF deviation
  real(sgl)         :: gbin(numhist)         ! GOF of bin
  real(sgl)         :: Gloc(0:numtalys)      ! local help variable
  real(sgl)         :: Gsum                  ! total GOF value
  real(sgl)         :: Gtmp                  ! help variable
  real(sgl)         :: maxgof                ! maximum GOF value
  real(sgl)         :: mingof                ! minimum GOF value
  real(sgl)         :: wloc(0:numtalys)      ! local help variable
  real(sgl)         :: wtmp                  ! help variable
  save    Gav, Gavdev
!
! Update distribution with previous result
!
! Running average and average deviation
!
  Gav(isearch) = 0.
  Gavdev(isearch) = 0.
  Nbeg = max(isearch - 100, 0)
  Nlast = isearch
  N = Nlast - Nbeg + 1
  if (N > 0) then
    Gsum = 0
    do i = Nbeg, Nlast
      Gsum = Gsum + Gsave(i)
    enddo
    Gav(isearch) = Gsum / N
    Gsum = 0
    do i = Nbeg, Nlast
      Gsum = Gsum + (Gsave(i) - Gav(isearch)) **2
    enddo
    Gavdev(isearch) = sqrt(Gsum / N)
  endif
!
! Make histogram
!
  mingof = 1.e38
  maxgof = 0.
  do i = 0, isearch
    if (Gsave(i) > 0.) mingof = min(Gsave(i), mingof)
    maxgof = max(Gsave(i), maxgof)
  enddo
  if (maxgof > mingof) then
    Nbin = numhist - 1
    dbin = (log(maxgof) - log(mingof)) / Nbin
    do k = 1, Nbin + 1
      gofhist(k) = 0
      gbin(k) = exp(log(mingof) + (k - 1) * dbin)
    enddo
    do i = 0, isearch
      if (Gsave(i) >= gbin(1) .and. Gsave(i) <= gbin(Nbin + 1)) then
        call locate_nr(gbin, Nbin + 1, Gsave(i), n)
        n = max(n, 1)
        gofhist(n) = gofhist(n) + 1
      endif
    enddo
    open (unit = 8, file = 'gof.hist', status = 'replace')
    write(8, '("# GOF bin   number")')
    do k = 1, Nbin
      write(8, '(f10.3, i5)') gbin(k), gofhist(k)
    enddo
    close(8)
    open (unit = 9, file = 'gof.all', status = 'replace')
    write(9, '("# number    GOF  Running av. Running std.  Weight")')
    do i = 0, isearch
      write(9, '(i5, 3f10.3, 5x, es12.5)') i, Gsave(i), Gav(i), Gavdev(i), weightsave(i)
    enddo
    close(9)
    do i = 1, Nchanall
      goffile = 'gof.'//xsfile(i)
      open (unit = 9, file = goffile, status = 'replace')
      write(9, '("# number    GOF for ", a20)') xsfile(i)
      do k = 0, isearch
        write(9, '(i5, f10.3)') k, Gchannelsave(k, i)
      enddo
    close(9)
    enddo
  endif
!
! Sort GOF values
!
  do i = 0, isearch
    iloc(i) = i
    Gloc(i) = Gsave(i)
    wloc(i) = weightsave(i)
  enddo
  do i = 0, isearch
    do j = i, isearch
      if (Gloc(i) > Gloc(j)) then
        Gtmp = Gloc(i)
        itmp = iloc(i)
        wtmp = wloc(i)
        Gloc(i) = Gloc(j)
        iloc(i) = iloc(j)
        wloc(i) = wloc(j)
        Gloc(j) = Gtmp
        iloc(j) = itmp
        wloc(j) = wtmp
      endif
    enddo
  enddo
  open (unit = 9, file = 'gof.sort', status = 'replace')
  write(9, '("# number    GOF   Weight")')
  do i = 0, isearch
    write(9, '(i5, f10.3, es12.5)') iloc(i), Gloc(i), wloc(i)
  enddo
  close(9)
  return
end subroutine histogof
! Copyright A.J. Koning 2021
