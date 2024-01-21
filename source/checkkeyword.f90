subroutine checkkeyword
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose   : Check for errors in keywords
!
! Author    : Arjan Koning
!
! 2022-05-09: Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_tasman_mod
!
! Variables for reading TASMAN input lines
!   inline          ! input line
!   nlines          ! number of input lines
! Variables for processing input
!   changeline      ! line number at which TASMAN specific input starts
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: numkey=136        ! number of keywords
  integer            :: i                 ! counter
  integer            :: j                 ! counter
  character(len=132) :: key               ! keyword
  character(len=132) :: keyword(numkey)   ! keyword
  character(len=132) :: word(40)          ! words on input line
!
! Although it is difficult to prevent the user from all possible input errors, we can check for the use of wrong keywords
! and for unphysical values for most of the input variables.
!
! *********************** Check for wrong keywords *********************
!
! TASMAN will stop if a keyword is incorrect
!
  data (keyword(i), i = 1, numkey) / '# ', '#acf', '#adeep', '#adevi', '#allpar', '#amax', '#amin', '#angle', &
    '#autoinclude', '#average', '#background', '#band', '#bdevi', '#block', &
    '#cdevi', '#chi2', '#chi2max', '#covrand', '#cross', '#cwidth', '#ddevi', '#defonly', '#deltae', '#dev', '#dexp', '#discrete', &
    '#dsmooth', '#e1vary', '#eaf', '#ecent', '#ecis', '#edependent', '#efracmax', '#ehigh', '#equisample', '#errlim', '#error', &
    '#esearch', '#esearch1', '#esearch2', '#eweight', '#expclass', '#expinclude', '#expexclude', '#extparvar', '#fiscor', '#fmax', &
    '#fns', '#format', '#fracmax', '#gamgam', '#gamma', '#gauss', &
    '#getcalcscript', '#global', '#gofmode', '#gpf', '#inponly', '#integral', '#isomer', '#kdburn', '#keyvary', '#ld', &
    '#legendre', '#libexclude', '#libinclude', '#liso', '#lrf7', '#m1vary', '#macs', '#magnet', '#maxbar', '#maxexpsets', &
    '#maxgam', '#mcak', '#metropolis', '#minbar', '#mode', '#mt', '#nburn', '#nhigh', '#njoy', '#ntalys', '#nubar', '#offset', &
    '#ompinc', '#omponly', '#outsearch', '#parameters', &
    '#partvary', '#proconly', '#prepro', '#prod', '#production', '#psf', '#purr', '#readpar', '#readsens', &
    '#residual', '#resonance', '#runfns', '#runnubar', '#runtares', '#s20', '#s30', '#s60', '#sample', '#save', &
    '#searchmode', '#seed', '#select', '#sens', '#sort', '#source', '#spectra', &
    '#tafislib', '#tafisversion', '#talexclude', '#talinclude', '#talysversion', '#taneslib', '#tanesversion', &
    '#tareslib', '#taresversion', '#tefalversion', '#tmc', '#tmcoffset', '#tweight', '#user', '#weight', '#weightpower', &
    '#zaskip', '#zavary', '#zdeep', '#zmax', '#zmin'/
!
! A keyword can be de-activated by putting a # in front of it.
! All first words of the input lines are checked against the list of keywords.
!
! getkeywords: subroutine to retrieve keywords and values from input line
!
! The keyword is identified.
!
Loop1:  do i = changeline + 1, nlines
    call getkeywords(inline(i), word)
    key = word(1)
    if (key(1:1) /= '#') cycle
    if (key(1:2) == '##') cycle
    do j = 1, numkey
      if (keyword(j) == key) cycle Loop1
    enddo
    write(*, '(/" TASMAN-error: Wrong keyword: ", a20)') key
    write(*, '(" Below #change a comment should start with ##")')
    stop
  enddo Loop1
  return
end subroutine checkkeyword
! Copyright A.J. Koning 2021
