subroutine expinterpol
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Interpolate theoretical energy grid on experimental grid
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
!   sgl           ! single precision kind
! All global variables
!   numenin       ! maximum number of incident energies
! Variables for writing TALYS input files
!   italys        ! TALYS run counter
! Variables for reading TALYS output
!   Nchanall      ! total number of channels
! Variables for GOF function
!   Eindex        ! index for experimental data
!   flagdeltaE    ! flag to weigh goodness - of - fit with dE of energy grid
!   Nindex        ! index for experimental data
! Variables for reading cross sections
!   E             ! incident energy
!   Enon          ! incident energy
!   Nen           ! number of incident energies
!   Nennon        ! number of incident energies
!   xseval        ! evaluated cross section
!   xsnon         ! nonelastic cross section
!   xsnon0        ! non - elastic crosssection
!   xsth          ! theoretical cross section
! Variables for reading experimental data
!   dEexp         ! uncertainty of incident energy
!   Ebeg          ! first energy of experimental energy grid
!   Eend          ! last energy of experimental energy grid
!   Eexp          ! incident energy
!   Nenexp        ! number of experimental data points
!   Nsets         ! number of experimental data sets
!
! *** Declaration of local data
!
  implicit none
  integer   :: i                        ! counter
  integer   :: isamp                    ! variable for 0th or random run
  integer   :: j                        ! counter
  integer   :: k                        ! Legendre order
  integer   :: n                        ! counter
  real(sgl) :: e1                       ! first energy
  real(sgl) :: e2                       ! second energy
  real(sgl) :: ee(0:numenin)            ! incident energy
  real(sgl) :: Eexp0                    ! incident energy
  real(sgl) :: fac                      ! factor
  real(sgl) :: xs1                      ! help variable
  real(sgl) :: xs2                      ! help variable
!
! *** Provide theoretical cross sections on experimental energy grid ***
!
  if (italys == 0) then
    isamp = 0
  else
    isamp = 1
  endif
  do i = 1, Nchanall
    do j = 0, Nen(i)
      ee(j) = E(i, j)
    enddo
    do j = 1, Nsets(i)
      do k = 1, Nenexp(i, j)
        Eexp0 = Eexp(i, j, k)
        if (Eexp0 <= E(i, Nen(i))) then
          call locate(ee, 0, Nen(i), Eexp0, n)
          e1 = ee(n)
          e2 = ee(n + 1)
          fac = (Eexp0 - e1) / (e2 - e1)
          xs1 = xseval(isamp, i, n)
          xs2 = xseval(isamp, i, n + 1)
          xsth(i, j, k) = xs1 + fac * (xs2 - xs1)
        else
          n = Nen(i)
          xsth(i, j, k) = xseval(isamp, i, n)
        endif
        n = max(n, 1)
        Eindex(i, j, k) = n
        Nindex(n) = Nindex(n) + 1
        if (Eexp0 <= Enon(Nennon)) then
          call locate(Enon, 0, Nennon, Eexp0, n)
          e1 = Enon(n)
          e2 = Enon(n + 1)
          fac = (Eexp0 - e1) / (e2 - e1)
          xs1 = xsnon0(n)
          xs2 = xsnon0(n + 1)
          xsnon(i, j, k) = xs1 + fac * (xs2 - xs1)
        else
          n = Nennon
          xsnon(i, j, k) = xsnon0(n)
        endif
!
! Make deltaE energy grid for energy-weighted optimization
!
        if (flagdeltaE) then
          dEexp(i, j, k) = 0.5 * abs(Eexp(i, j, k + 1) - Eexp(i, j, k - 1))
          if (Eexp0 > 0..and.Eexp(i, j, k -1) == 0.) then
            Ebeg(i, j) = Eexp0
            dEexp(i, j, k) = abs(Eexp(i, j, k + 1) - Eexp0)
          endif
          if (Eexp0 > 0..and.Eexp(i, j, k + 1) == 0.) then
            Eend(i, j) = Eexp0
            dEexp(i, j, k) = abs(Eexp0 - Eexp(i, j, k - 1))
          endif
        endif
      enddo
    enddo
  enddo
  return
end subroutine expinterpol
! Copyright A.J. Koning 2021
