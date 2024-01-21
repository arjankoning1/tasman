subroutine optimization
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Optimization to experimental data
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
!   sgl        ! single precision kind
! All global variables
!   numpar     ! maximum number of parameters
! Variables for writing TALYS input files
!   mode       ! TASMAN mode
! Variables for reading TALYS output
!   iout       ! counter for output files
! Variables for fitting limits
!   parhigh    ! upper value of parameter uncertainty
!   parlow     ! lower value of parameter uncertainty
! Variables for GOF function
!   searchmode ! search mode
! Variables for parameter variation
!   Npar       ! number of parameters
!   parinp     ! input parameter value
!
! *** Declaration of local data
!
  implicit none
  real(sgl) :: Fopt                    ! optimal function value
  real(sgl) :: Popt(numpar)            ! optimal parameters
!
! ********************* Cross section fitting **************************
!
! expfiles  : subroutine to read experimental cross sections for all nuclides
! fitlimits : subroutine to set parameter limits for optimization
! inputwrite: subroutine to write parameters for TALYS input file
! Fopt      : optimal function value
!
  if (mode == 4) call expfiles
  call fitlimits
  if (Npar >= 1) then
    call koekel(Npar,parinp,parlow,parhigh,Popt,Fopt,searchmode)
    if (iout.ge.Npar.and.searchmode.ne.4) then
      searchmode=4
      call koekel(Npar,parinp,parlow,parhigh,Popt,Fopt,searchmode)
    endif
  else
    call fcn(Npar, parinp, Fopt)
  endif
  return
end subroutine optimization
! Copyright A.J. Koning 2021
