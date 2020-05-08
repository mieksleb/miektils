module readers
  implicit none 
  public :: reader

contains

subroutine reader(conf_file_name,top_file_name,step,bp,x1,y1,z1,x2,y2,z2,reverse,circular,energy_out)

  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2,x2_rev,y2_rev,z2_rev
  integer :: i,bp,nlines,io,step,bases,strands,strandid,circ_int,dum_int
  logical :: energy_out,circular,reverse
  character :: conf_file_name*20,top_file_name*20, search_string*4, top_line*20, energy_out_string*4
  character :: base*1, second_line*20
  real :: tol=0.01


  open(3,file=top_file_name)
  read(3,"(A20)") top_line
  read(top_line,*) bases,strands
  bp = int(bases/strands)
  read(3,"(A20)") second_line
  read(second_line,*) strandid,base,circ_int,dum_int
  if (circ_int.eq.-1) then
    circular=.False.
  else
    circular=.True.
  end if


  search_string = 't = '
  open(1,file=conf_file_name)
  nlines=0

  ! searches for string pattern "t = " to set logical energy_out and obatin timestep
  read(1,"(A20)") top_line
  read(top_line,'(A4)') energy_out_string
  close(1)
  open(1,file=conf_file_name)
  if (energy_out_string==search_string) then
    read (1,"(A4,I8)") search_string,step
    energy_out=.True.
  else
    energy_out=.False.
  end if
  close(1)


  ! now we load the base pair positions of both strands

  if (circular.eqv..False.) then
    open(1, file=conf_file_name)
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
      x1(bp+1)=x1(1)
      y1(bp+1)=y1(1)
      z1(bp+1)=z1(1)
      x2(bp+1)=x2(1)
      y2(bp+1)=y2(1)
      z2(bp+1)=z2(1)
    else
      do i = 1,3
        read(1,*,iostat=io)
      end do
      do i = 1,bp
        read(1,*,iostat=io) x1(i), y1(i), z1(i)
      end do
      do i = 1, bp
        read(1,*) x2(i), y2(i), z2(i)
      end do
    end if
    if (reverse.eqv..True.) then
    allocate(x2_rev(bp))
    allocate(y2_rev(bp))
    allocate(z2_rev(bp))
      do i=1,bp
        x2_rev(i)=x2(bp+1-i)
        y2_rev(i)=y2(bp+1-i)
        z2_rev(i)=z2(bp+1-i)
      end do
    x2(:)=x2_rev(:)
    y2(:)=y2_rev(:)
    z2(:)=z2_rev(:)
    else
    end if
    close(1)

  else if (circular.eqv..True.) then
    open(1, file=conf_file_name)
    allocate(x1(bp+1))
    allocate(y1(bp+1))
    allocate(z1(bp+1))
    allocate(x2(bp+1))
    allocate(y2(bp+1))
    allocate(z2(bp+1))
    if (energy_out.eqv..False.) then
      do i = 1,bp
        read(1,*) x1(i), y1(i), z1(i)
      end do
      do i = 1, bp
        read(1,*) x2(i), y2(i), z2(i)
      end do
      x1(bp+1)=x1(1)
      y1(bp+1)=y1(1)
      z1(bp+1)=z1(1)
      x2(bp+1)=x2(1)
      y2(bp+1)=y2(1)
      z2(bp+1)=z2(1)
    else
      do i = 1,3
        read(1,*,iostat=io)
      end do
      do i = 1,bp
        read(1,*,iostat=io) x1(i), y1(i), z1(i)
      end do
      do i = 1, bp
        read(1,*) x2(i), y2(i), z2(i)
      end do
      x1(bp+1)=x1(1)
      y1(bp+1)=y1(1)
      z1(bp+1)=z1(1)
      x2(bp+1)=x2(1)
      y2(bp+1)=y2(1)
      z2(bp+1)=z2(1)
    end if
    
    if (reverse.eqv..True.) then
    allocate(x2_rev(bp+1))
    allocate(y2_rev(bp+1))
    allocate(z2_rev(bp+1))
      do i=1,bp+1
        x2_rev(i)=x2(bp+2-i)
        y2_rev(i)=y2(bp+2-i)
        z2_rev(i)=z2(bp+2-i)
      end do
    x2(:)=x2_rev(:)
    y2(:)=y2_rev(:)
    z2(:)=z2_rev(:)
    else
    end if
    close(1)
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


