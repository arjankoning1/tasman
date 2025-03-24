subroutine mainout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Main output
!
! Author    : Arjan Koning
!
! 2023-12-29: Original code
! 2025-03-10: Current revision
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
!
! Main output
!
  write(*,'(/"    TASMAN-2.1 (Version: March 10, 2025)"/)')
  write(*, '(" Copyright (C) 2025  A.J. Koning"/)')
  write(*, '(" Dimensions - Cross sections: mb, Energies: MeV, ", "Angles: degrees")')
  write(*, '(/" User: ",a)') trim(user)
  write(*, '(" Date: ",a)') trim(date)
!
! ***************** Write input file and default parameters ************
!
! inputout: subroutine to write the input parameters
!
  call inputout
  return
end subroutine mainout
! Copyright A.J. Koning 2025
