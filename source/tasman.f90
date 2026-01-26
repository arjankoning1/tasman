program tasman
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Uncertainty quantifiaction code for TALYS
!
! Author    : Arjan Koning
!
! 2023-12-29: Original code
! 2026-01-25: Current revision
!-----------------------------------------------------------------------------------------------------------------------------------
!
!   |-------------------------------------------------------|
!   |                 TASMAN-2.2                            |
!   |                 Arjan Koning                          |
!   |                                                       |
!   | Email: A.Koning@iaea.org                              |
!   |-------------------------------------------------------|
!
! MIT License
!
! Copyright (c) 2025 Arjan Koning
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
! *** Use data from other modules
!
  use A0_tasman_mod
!
! Variables for writing TALYS input files
!   mode    ! TASMAN mode
!
! ************* Input, initialization and reaction models **************
!
! machine      : subroutine for machine dependent statements
! constants    : subroutine for constants and basic properties of particles
! tasmaninitial: subroutine for initialization
! tasmaninput  : subroutine for user input and defaults
!
  call machine
  call constants
  call tasmaninitial
  call tasmaninput
!
! ******** Monte Carlo uncertainties, covariances and optimization *****
!
! uncertainty : subroutine for uncertainty propagation and sensitivities
! optimization: subroutine for optimization to experimental data
! timer       : subroutine for output of execution time
!
  if (mode <= 2) then
    call uncertainty
  else
    call optimization
  endif
  call timer
end program tasman
! Copyright A.J. Koning 2026
