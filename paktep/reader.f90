module readers
  implicit none 
  public :: reader

contains

subroutine reader(file_name,step,bp,x1,y1,z1,x2,y2,z2,reverse,circular,energy_out)

  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2
  integer :: i,bp,nlines,io,step
  logical :: energy_out,circular,reverse
  character :: file_name*20, search_string*4, top_line*20, energy_out_string*4
  real :: tol=0.0001

  search_string = 't = '
  open(1,file=file_name)
  nlines=0

  ! searches for string pattern "t = " to set logical energy_out and obatin timestep
  read(1,"(A20)") top_line
  read (top_line,'(A4)') energy_out_string
  close(1)
  open(1,file=file_name)
  if (energy_out_string==search_string) then
    read (1,"(A4,I8)") search_string,step
    energy_out=.True.
  else
    energy_out=.False.
  end if
  close(1)

  open(1,file=file_name)
  ! if energy_out = True, then there will be 3 lines to skip
  if (energy_out.eqv..False.) then 
    do
      read(1, *, iostat=io) 
      if (io/=0) exit
      nlines = nlines +1
    end do
  else
    do i = 1,3
      read(1,*, iostat=io)
    end do
    do 
      read(1,*,iostat=io)
      if (io/=0) exit
      nlines = nlines + 1
    end do
  end if
  close(1)

  bp = nlines/2

  ! now we load the base pair positions of both strands
  open(1, file=file_name)
  allocate(x1(bp))
  allocate(y1(bp))
  allocate(z1(bp))
  allocate(x2(bp))
  allocate(y2(bp))
  allocate(z2(bp))
  if (energy_out.eqv..False.) then
    do i = 1,bp
      read(1,*) x1(i), y1(i), z1(i)
    end do
    do i = 1, bp
      read(1,*) x2(i), y2(i), z2(i) 
    end do
  else
    do i = 1,3
      read(1,*,iostat=io)
    end do
    do i = 1,bp
      read(1,*,iostat=io) x1(i), y1(i), z1(i)
    end do
    do i = 1,bp
      if (reverse.eqv..True.) then
        read(1,*,iostat=io) x2(bp+1-i), y2(bp+1-i), z2(bp+1-i)
      else
        read(1,*,iostat=io) x2(i), y2(i), z2(i)
      end if
    end do
  end if
  close(1)

  if((abs(x1(1)-x1(bp)).le.tol).and.(abs(y1(1)-y1(bp)).le.tol).and.(abs(z1(1)-z1(bp)).le.tol)) then
    circular=.True.
  else
    circular=.False.
  end if



end subroutine reader


! This subroutine splits a trajectory file into separate configuration files
!subroutine splitter(file_name)
!
!  integer :: i,nlines,bp,steps,ui,io,step,linesperfile,neme_pos
!  character ::  init_conf_name*20,traj_file_name*20,final_conf_name*20, dir_name*4
!  character :: search_string*4, energy_out_string*4
!  character :: twist_writhe_file_name*20,plectoneme_file_name*20
!  real :: twist_integral,writhe_integral
!
!  search_string = 't = '
!
!
!!  read(*,*) init_conf_name
!  read(*,*) traj_file_name
!  read(*,*) final_conf_name
!
!  ! this script splits a trajectory file into separate configuration files so that scripts can be run on each config.
!  ! these conf.dat files are stored in a temporary directory call 'dump'
!  ! after the conf files are looped through, the dump folder is deleted
!  dir_name = 'dump'
!  write(cmd,"(A20)") 'mkdir '//dir_name
!  call system(cmd)
!  write(cmd,"(A20)") 'cd '//dir_name
!  call system(cmd)
!
!  ! call reader on the initical config to get bp,circular and energy_out
!  call reader(final_conf_name,step,bp,x1,y1,z1,x2,y2,z2,reverse.eqv..False.,circular,energy_out)
!  if (energy_out.eqv..True.) then
!    linesperfile=2*bp+3
!  else
!    linesperfile=2*bp
!  end if
!  write(*,*) linesperfile
!  write(lpf_string,"(I5)") linesperfile
!  lpf_string=adjustl(trim(lpf_string))
!  open(1,file=traj_file_name)
!  cmd = "gsplit --additional-suffix=.dat -d -l "//lpf_string//" -a 4 "//traj_file_name//" conf" ! this command is 'split' for linux and 'gsplit' for mac
!                                                                                                ! if you brew install coreutils
!  write(*,*) cmd
!  call system(cmd)
!  write(cmd,"(A20)") 'cd ../'
!  call system(cmd)
!
!
!end subroutine splitter

end module readers


