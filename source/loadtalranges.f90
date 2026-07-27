subroutine loadtalranges
  use A0_tasman_mod
  use A1_error_handling_mod
  ! Global variables used
  !
  ! talys             ! name of talys executable
  ! parkey            ! array holding key names of parameters to be optimized
  ! parlow            ! array holding minimum acceptable value for parameters to be optimized
  ! parhigh           ! array holding maximum acceptable value for parameters to be optimized
  !
  ! Global subroutines used
  !
  ! read_error        ! Throw error for fileIO error and terminate program
  ! getkeywords       ! Get key words and associated values from input file line

  implicit none
  logical            :: lexist       ! Boolean for existence of file
  integer            :: istat        ! File IO state
  character(len=256) :: command      ! Shell command
  integer            :: cmdstat      ! Shell command status
  character(len=132) :: word(40)     ! To hold key-value pairs per line in range file
  character(len=132) :: line         ! To hold line from file temporarily
  character(len=256) :: rangeFile    ! Path to keyword range file
  integer            :: numkeys      ! number of keys stored in key.ranges of TALYS
  integer            :: i            ! Iterator
  integer            :: j            ! Iterator
  integer            :: ix           ! String index

  numkeys=0
  i=1
  ! Fix for now: pipe shell output to temp file to get the path to the ranges file
  ! Only works for linux based OSes at the moment (also need access to realpath, dirname)
  ! Based on minimal testing, should work for symlinks
  ! Probably not safe if TASMAN is called many times from the same directory
  command='echo $(dirname $(realpath '// talys
  command=trim(command)//')) > temp.txt'
  print *, "Length:", len(command), command
  cmdstat=system(trim(command))
  if (cmdstat/=0) then
    write(*, '(" TASMAN-error: Unable to create temporary file with command: ", &
    & a)') command
    stop
  end if
  open(13, file='temp.txt', status='unknown', iostat=istat)
  if (istat/=0) call read_error('temp.txt', istat)
  read (13, '(a)') rangeFile
  close(unit=13)
  print *, 'range file:', rangeFile
  cmdstat=system('rm temp.txt')
  ix=index(rangeFile,'/bin')
  rangeFile=trim(rangeFile(1:ix))//'misc/key.ranges'
  if (cmdstat/=0) then
    write(*, '(" TASMAN-error: Unable to remove temporary file with command: ", &
    & "rm temp.txt ")') 
    stop
  end if

  inquire(file=trim(rangeFile), exist=lexist)
  if (lexist) then
    open(3, file=trim(rangeFile), status='old', iostat=istat)
    if (istat/=0) call read_error(rangeFile, istat)
    do
      read(3, '(a132)', iostat=istat) line
      if (istat==-1) exit
      if (istat/=0) call read_error(rangeFile, istat)
      if (line(1:1)=='#') cycle
      numkeys= numkeys + 1
    end do
    ! Allocate space for talkeyranges based on number of keys stored &
    ! reopen file to obtain contents
    close(unit=3)
    allocate(talkeyranges(numkeys,3))
    open(3, file=trim(rangeFile), status='old', iostat=istat)
    do
      read(3, '(a132)', iostat=istat) line
      if (istat==-1) exit
      if (line(1:1)=='#') cycle
      call getkeywords(line, word)
      talkeyranges(i,1)=trim(word(1))
      talkeyranges(i,2)=trim(word(2))
      talkeyranges(i,3)=trim(word(3))
      i=i+1
    end do
    close(unit=3)
    ! Iterate over stored key values to vary and assign ranges
    do i=1, size(parkey)
      do j=1, numkeys
        if (parkey(i)==talkeyranges(j,1)) then
          ! If no explicit value given, pass over it idk ill figure this out later maybe
          ! As of 07/24/26, ~15 variables dont have explicit values (defined based on tal
          ! input or A0_talys_mod `constants`)
          if (talkeyranges(j,2)/='N/A') then
            read(talkeyranges(j,2), *) parlow(i)
          else
            continue
          end if
          if (talkeyranges(j,3)/='N/A') then
            read(talkeyranges(j,3), *) parhigh(i)
          else
            continue
          end if
        end if
      end do
    end do
    deallocate(talkeyranges)
    return
  end if
  write(*, '(" TASMAN-error: Unable to locate TALYS key word range file:", a)') trim(rangeFile)
  stop
end subroutine loadtalranges