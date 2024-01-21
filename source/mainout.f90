subroutine mainout
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Main output
!
! Author    : Arjan Koning
!
! 2023-12-29: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Declaration of local data
!
  write(*,'(/"    TASMAN-2.0 (Version: December 29, 2023)"/)')
  write(*, '(" Copyright (C) 2023  A.J. Koning"/)')
  write(*, '(" Dimensions - Cross sections: mb, Energies: MeV, ", "Angles: degrees")')
!
! ***************** Write input file and default parameters ************
!
! inputout: subroutine to write the input parameters
!
  call inputout
  return
end subroutine mainout
! Copyright A.J. Koning 2023
