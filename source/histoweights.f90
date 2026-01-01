subroutine histoweights
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Make histograms per parameter based on weights
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
!   sgl         ! single precision kind
! All global variables
!   numhist     ! maximum number of points for histogram
! Variables for covariances
!   weight      ! weight for experimental data sets
! Variables for parameter variation
!   Npar        ! number of parameters
!   par         ! parameter value
! Variables for GOF function
!   isearch     ! number of trial run
! Variables for fitting limits
!   hist        ! value in histogram bin
!   histsam     ! value in histogram bin
!   nhist       ! number of histogram bins
!   Nhistbin    ! number of histogram bins
!   nhistsam    ! number of hits in bin
!   parfile     ! file with parameters
!   Pbin        ! central value of histogram bin
!   Pbot        ! lower value of histogram bin
!   readpar     ! flag to read parameter distribution
!
! *** Declaration of local data
!
  implicit none
  character(len=34) :: pfile                ! parameter file
  integer           :: i                    ! counter
  integer           :: ipar                 ! counter for parameter
  integer           :: j                    ! counter
  integer           :: k                    ! Legendre order
  real(sgl)         :: frac                 ! help variable
  real(sgl)         :: hweight              ! weight for histogram
  real(sgl)         :: ohist(numhist)       ! output histogram value
  real(sgl)         :: ohistsam(numhist)    ! sampled histogram value
  real(sgl)         :: P                    ! sampled parameter
  real(sgl)         :: Pav                  ! average parameter value
  real(sgl)         :: Pdev                 ! parameter deviation
  real(sgl)         :: Plow(numhist)        ! bottom of histogram
  real(sgl)         :: Pvar                 ! parameter variation
  real(sgl)         :: reldev               ! relative deviation
  real(sgl)         :: sum                  ! help variable
  real(sgl)         :: suma                 ! help variable
  real(sgl)         :: sumd                 ! help variable
!
! Update distribution with previous result
!
  do ipar = 1, Npar
    do i = 1, Nhistbin(ipar) + 1
      Plow(i) = Pbot(ipar, i)
    enddo
    P = par(ipar)
    hweight = weight(0)
    if (P >= Plow(1) .and. P <= Plow(Nhistbin(ipar) + 1)) then
      call locate_nr(Plow, Nhistbin(ipar) + 1, P, j)
      if (j > 0) then
        if (readpar(ipar) == 0) then
          hist(ipar, j) = (hist(ipar, j) * nhist(ipar, j) + hweight) / (nhist(ipar, j) + 1)
          nhist(ipar, j) = nhist(ipar, j) + 1
        endif
        histsam(ipar, j) = (histsam(ipar, j) * nhistsam(ipar, j) + hweight) / (nhistsam(ipar, j) + 1)
        nhistsam(ipar, j) = nhistsam(ipar, j) + 1
      endif
    else
      write(*, '(" TASMAN-warning: parameter ", a24, " with value ", es12.5, " outside boundaries", es12.5," - ", es12.5)') &
 &      parfile(ipar), P, Plow(1), Plow(Nhistbin(ipar) + 1)
    endif
!
! Normalize weights to unity for histogram output files
!
    sum = 0.
    do j = 1, Nhistbin(ipar) + 1
      sum = sum + histsam(ipar, j)
    enddo
    if (sum > 0.) then
      frac = 1. / sum
    else
      frac = 1.
    endif
    do j = 1, Nhistbin(ipar) + 1
      if (readpar(ipar) == 0) then
        ohist(j) = frac * hist(ipar, j)
      else
        ohist(j) = hist(ipar, j)
      endif
      ohistsam(j) = frac * histsam(ipar, j)
    enddo
!
! Calculate average and standard deviation
!
    sum = 0.
    suma = 0.
    do j = 1, Nhistbin(ipar) + 1
      sum = sum + ohistsam(j)
      suma = suma + Pbin(ipar, j) * ohistsam(j)
    enddo
    if (sum > 0.) then
      Pav = suma / sum
    else
      Pav = 0.
    endif
    sumd = 0.
    do j = 1, Nhistbin(ipar) + 1
      sumd = sumd + (Pbin(ipar, j) - Pav) **2 * ohistsam(j)
    enddo
    if (sum > 0.) then
      Pvar = sumd / sum
      Pdev = sqrt(Pvar)
    else
      Pvar = 0.
      Pdev = 0.
    endif
    if (Pav > 0.) then
      reldev = 100. * Pdev / Pav
    else
      reldev = 0.
    endif
!
! Write distribution to file
!
    pfile = '                                  '
    do k = 1, 2
      if (k == 1) then
        pfile = 'his.'//parfile(ipar)
      else
        if (mod(isearch, 100) /= 0) cycle
        pfile = 'his.0000.'//parfile(ipar)
        write(pfile(5:8), '(i4.4)') isearch
      endif
      open (unit = 2, file = pfile, status = 'replace')
      write(2, '("# Histogram after ", i5, " runs")') isearch
      write(2, '("# Range         : ", es12.5, " - ", es12.5)') Plow(1), Plow(Nhistbin(ipar) + 1)
      write(2, '("# Average       : ", es12.5, " Variance: ", es12.5," Stand. dev.: ", e12.5, " = ", f15.8, " %")') &
 &      Pav, Pvar, Pdev, reldev
      write(2, '("# Number of bins: ", i5)') Nhistbin(ipar)
      if (readpar(ipar) == 1) then
        write(2, '("#    Bin     Prior Weight Number", " Sampled weight Number")')
        do j = 1, Nhistbin(ipar)
          write(2, '(es12.5, 2(es12.5, i5, 5x))') Pbin(ipar, j), ohist(j), nhist(ipar, j), ohistsam(j), nhistsam(ipar, j)
        enddo
      else
        write(2, '("#    Bin   Sampled weight Number")')
        do j = 1, Nhistbin(ipar)
          write(2, '(2es12.5, i5)') Pbin(ipar, j), ohistsam(j), nhistsam(ipar, j)
        enddo
      endif
      close(2)
    enddo
  enddo
  return
end subroutine histoweights
! Copyright A.J. Koning 2021
