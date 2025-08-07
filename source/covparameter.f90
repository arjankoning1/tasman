subroutine covparameter
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Covariance matrices for parameters
!
! Author    : Arjan Koning
!
! 2021-12-30: Original code
! 2024-02-01: Beautify output
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
! Variables for writing TALYS input files
!   italys         ! TALYS run counter
! Variables for reading experimental data
!   Ein            ! incident energy
! Variables for weights
!   flagweight     ! flag to use weights for random samples
! Variables for GOF function
!   Ntalys         ! number of TALYS runs
! Variables for fitting limits
!   parfile        ! file with parameters
! Variables for reading cross sections
!   Nen            ! number of incident energies
!   Qval           ! Q - value
! Variables for parameter variation
!   Npar           ! number of parameters
!   pardelta       ! uncertainty of parameter
!   parinp         ! input parameter value
!   parsave        ! parameter value
!   partalys       ! parameter value
! Variables for parameter covariances
!   parav          ! average parameters
!   parcor         ! correlation of parameters
!   parcov         ! covariance matrix for parameters
!   pardif         ! difference of parameters
!   parstring      ! parameter string
! Variables for covariances
!   Sweight        ! weight for TALYS run
!   Sweightprev    ! weight for previous TALYS run
!   Sweightsum     ! sum of weights for TALYS run
!
! *** Declaration of local data
!
  implicit none
  character(len=25) :: pfile      ! parameter file
  character(len=15) :: col(9)                         ! header 
  character(len=15) :: un(9)
  character(len=16) :: reaction                       
  character(len=132) :: quantity                      
  character(len=132) :: topline    ! topline          
  integer           :: j          ! counter
  integer           :: k          ! Legendre order
  integer           :: l          ! counter
  integer           :: Ncol                           
  integer           :: indent
  integer           :: id2
  real(sgl)         :: Ei         ! incident energy
  real(sgl)         :: Eprev      ! previous energy
  real(sgl)         :: err        ! error
  real(sgl)         :: finperc    ! relative parameter uncertainty on output
  real(sgl)         :: inperr     ! parameter uncertainty on input
  real(sgl)         :: inpperc    ! relative parameter uncertainty on input
  real(sgl)         :: pdifk      ! Rdifference in parameter
  real(sgl)         :: pdifl      ! Rdifference in parameter
  real(sgl)         :: ratio      ! ratio of uncertainties
  real(sgl)         :: term       ! help variable
!
! Output of accepted and rejected parameters
!
  indent = 0
  id2 = indent + 2
  col = ''
  un = ''
  reaction='('//ptype0//',x)'
  quantity='parameter covariance matrix'
  topline=trim(targetnuclide)//trim(reaction)//' '//trim(quantity)
!
! parameter covariance matrix
!
  pdifk = 0.
  pdifl = 0.
  do k = 1, Npar
    if (parsave(0, k) /= 0.) pdifk = (partalys(1, k) - parsave(0, k)) / parsave(0, k)
    if (Sweightsum(0) > 0.) then
      do l = 1, Npar
        if (parsave(0, l) /= 0.) pdifl = (partalys(1, l) - parsave(0, l)) / parsave(0, l)
        term = pdifk * pdifl
        parcov(k, l) = (parcov(k, l) * Sweightprev(0) + Sweight(0) * term) / Sweightsum(0)
      enddo
      pardif(k) = (pardif(k) * Sweightprev(0) + Sweight(0) * pdifk) / Sweightsum(0)
    endif
    do j = 0, numenin
      if (Sweightsum(j) > 0.) parav(k, j) = (parav(k, j) * Sweightprev(j) + &
        Sweight(j) * partalys(1, k)) / Sweightsum(j)
    enddo
  enddo
  do k = 1, Npar
    do l = 1, Npar
      if (pardif(k) /= 0..and.pardif(l) /= 0.) parcor(k, l) = parcov(k, l) / (sqrt(abs(parcov(k, k))) * sqrt(abs(parcov(l, l))))
    enddo
  enddo
!
! Output of covariance matrices for parameters
!
  col(1)='parameter_1'
  col(3)='parameter_2'
  col(5)='Covariance'
  col(6)='Correlation'
  Ncol=6
  open (unit = 1, file = 'cov_parameter.ave', status = 'replace')
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  call write_char(indent,'covariance','')
  call write_char(id2,'class: covariance','')
  call write_integer(id2,'number of accepted runs',italys)
  call write_quantity(id2,quantity)
  call write_datablock(id2,Ncol,Npar*Npar,col,un)
  do k = 1, Npar
    do l = 1, Npar
      write(1, '(2a30, 2es15.6)') parstring(k), parstring(l), parcov(k, l), parcor(k, l)
    enddo
  enddo
  close (1)
!
! Output of average parameters
!
  col(1)='parameter'
  col(3)='Average'
  col(4)='Std._dev.'
  col(5)='Rel._dev.'
  un(5)='%'
  col(6)='Input'
  col(7)='Std._dev.'
  col(8)='Rel._dev.'
  un(8)='%'
  col(9)='Unc._ratio'
  Ncol=9
  open (unit = 1, file = 'parameter.ave', status = 'replace')
  call write_header(indent,topline,source,user,date,oformat)
  call write_target(indent)
  call write_covariance(indent,reaction,0,0,0,0,italys)
  call write_quantity(indent,quantity)
  call write_datablock(indent,Ncol,Npar,col,un)
  do k = 1, Npar
    err = parav(k, 0) * sqrt(abs(parcov(k, k)))
    inperr = parsave(0, k) * pardelta(k)
    inpperc = 0.
    finperc = 0.
    ratio = 0.
    if (parsave(0, k) /= 0.) inpperc = 100. * inperr / parsave(0, k)
    if (parav(k, 0) /= 0.) finperc = 100. * err / parav(k, 0)
    if (finperc /= 0.) ratio = inpperc / finperc
    write(1, '(a30, 7es15.6)') parstring(k), parav(k, 0), err, finperc, parsave(0, k), inperr, inpperc, ratio
  enddo
  close (1)
!
! Output of energy-dependent average parameters
!
  if (flagweight) then
    do k = 1, Npar
      if (italys == Ntalys) then
        pfile = parfile(k)
      else
        pfile = 'upd.'//trim(parfile(k))
      endif
      open (unit = 1, file = pfile, status = 'replace')
      Eprev = 0.
      do j = 1, Nen(1)
        if (parinp(k) /= 0.) then
          ratio = parav(k, j) / parinp(k)
          if (pfile == 'rspincut' .or. pfile == 'm2constant') then
            Ei = Ein(j) + Qval
          else
            Ei = Ein(j)
          endif
          if (Ei > Eprev) write(1, '(2es12.5)') Ei, ratio
          Eprev = Ei
        endif
      enddo
      close(1)
    enddo
  endif
  return
end subroutine covparameter
! Copyright A.J. Koning 2021
