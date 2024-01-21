subroutine tasmaninput
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : User input and defaults
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
! readinput   : subroutine to read input
! input1      : subroutine to read input for first set of variables
! input2      : subroutine to read input for second set of variables
! input3      : subroutine to read input for third set of variables
! input4      : subroutine to read input for fourth set of variables
! inputprocess: subroutine to process input parameters
! checkkeyword: subroutine to check for errors in keywords
! checkvalue  : subroutine to check for errors in values
! mainout     : subroutine for main output
!
  call readinput
  call input1
  call input2
  call input3
  call input4
  call inputprocess
  call checkkeyword
  call checkvalue
  call mainout
  return
end subroutine tasmaninput
! Copyright A.J. Koning 2021
