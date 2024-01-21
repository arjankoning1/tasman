subroutine weights
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Determine weights of random runs
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
!   sgl            ! single precision kind
! All global variables
!   numenin        ! maximum number of incident energies
! Variables for weights
!   gofmode        ! gof mode
! Variables for reading experimental data
!   Ein            ! incident energy
! Variables for reading cross sections
!   Nen            ! number of incident energies
! Variables for covariances
!   Sweight        ! weight for TALYS run
!   weight         ! weight for experimental data sets
! Variables for GOF function
!   Gall           ! GOF value
!   Gall0          ! GOF value
!   Grun           ! GOF value of run
!   Grun0          ! GOF value of run
!   ichi2          ! parameter to set goodness - of - fit estimator
!   isearch        ! number of trial run
! Variables for weights
!   flagEweight    ! flag to use energy - dependent weights
!   weightpower    ! power for weight of random run
!   weightsave     ! weight for TALYS run
!
! *** Declaration of local data
!
  implicit none
  integer   :: j                               ! counter
  integer   :: nend                            ! help variable
  real(sgl) :: expo                            ! exponent
  real(sgl) :: expo0                           ! exponent
  real(sgl) :: expoD                           ! exponent
  real(sgl) :: G                               ! goodness-of-fit value  : G = 1: chi-2 : G = 2: squared distance : G = 3: F factor :
  real(sgl) :: G0                              ! GOF for run 0
  real(sgl) :: GA                              ! total GOF value
  real(sgl) :: GA0                             ! total GOF value of run 0
  real(sgl) :: sumweight(0:numenin)            ! total cumulative weight
  real(sgl) :: w0                              ! weight of run 0
  save    sumweight
!
! ******************** Set total weight ********************************
!
  if (flagEweight) then
    nend = numenin
  else
    nend = 0
  endif
  if (isearch == 0) then
    open (unit = 19, file = 'weights', status = 'replace')
    write(19, '(" Run    Weight      Sum of weights")')
    do j = 0, nend
      sumweight(j) = 0.
    enddo
  endif
!
! gofmode 1: w = exp(-chi2point/2)/w0
! gofmode 2: w = exp(-(chi2point/chi2point0)/2)/w0
! gofmode 3: w = exp(-(chi2point/chi2point0)**2)/w0
! gofmode 4: w = exp(-chi2/2)/w0
! gofmode 5: w = exp(-(chi2/chi20)/2)/w0
! gofmode 6: w = exp(-(chi2/chi20)**2)/w0
!
  do j = 0, nend
    G0 = Grun0(j)
    if ((ichi2 == 2 .or. ichi2 == 3) .and. G0 <= 1.) G0=1.
    GA0 = Gall0(j)
    if ((ichi2 == 2 .or. ichi2 == 3) .and. G0 <= 1.) GA0=1.
    G = Grun(j)
    GA = Gall(j)
    if (G0 > 0.) then
      if (gofmode == 1 .or. gofmode == 4) then
        if (gofmode == 1) then
          expo = abs(G)
          expo0 = abs(G0)
        else
          expo = abs(GA)
          expo0 = abs(GA0)
        endif
        if (expo < 80..and.expo0 < 80.) then
          w0 = exp( - 0.5 * expo0)
          weight(j) = exp( - 0.5 * expo) / w0
        endif
      endif
      if (gofmode == 2 .or. gofmode == 5) then
        if (gofmode == 2) then
          expo = abs(G)
          expo0 = abs(G0)
        else
          expo = abs(GA)
          expo0 = abs(GA0)
        endif
        expoD = expo / expo0
        if (expoD < 80.) then
          w0 = exp( - 0.5)
          weight(j) = exp( - expoD) / w0
        endif
      endif
      if (gofmode == 3 .or. gofmode == 6) then
        if (gofmode == 3) then
          expo = abs(G / G0) **weightpower
        else
          expo = abs(GA / GA0) **weightpower
        endif
        if (expo < 80.) then
          w0 = exp( -1.)
          weight(j) = exp( - expo) / w0
        endif
      endif
      if (expo >= 80.) weight(j) = 0.
      sumweight(j) = sumweight(j) + weight(j)
      if (j == 0) then
        write(*, '("Weight of run ", i5, ":", es12.5, " Sum:", es12.5)') isearch, weight(0), sumweight(0)
        if (flagEweight) write(*, '(" Energy        Weight     Sum")')
      endif
      if (flagEweight .and. j > 0) then
        if (Ein(j) > 0.) write(*, '(3es12.5)') Ein(j), weight(j), sumweight(j)
      endif
    else
      weight(j) = 0.
    endif
    if (j == 0) write(19, '(i5, 2es12.5)') isearch, weight(0), sumweight(0)
  enddo
  weightsave(isearch) = weight(0)
  if (sumweight(0) > 0.) then
    if (flagEweight) then
      do j = 0, Nen(1)
        Sweight(j) = weight(j)
      enddo
    else
      do j = 0, Nen(1)
        Sweight(j) = weight(0)
      enddo
     endif
  endif
  return
end subroutine weights
! Copyright A.J. Koning 2021
