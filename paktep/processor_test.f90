program processor
  use readers
  implicit none
  

  integer :: i,nlines,bp,step,steps,ui
  character ::  init_conf_name*20,traj_file_name*20,final_conf_name*20, dir_name*4
  logical :: reverse,circular,energy_out
  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2
  character :: frmt="(I5.5)", steps_string*8, cmdfmt="(A20)", cmd*20
  dir_name = 'dump'
  write(cmd,"(A20)") 'mkdir '//dir_name
  write(*,*) cmd
  call system(cmd)
  read(*,*) init_conf_name
  read(*,*) traj_file_name
  read(*,*) final_conf_name

  ! call reader on the initical config to get bp,circular and energy_out
  call reader(init_conf_name,step,bp,x1,y1,z1,x2,y2,z2,reverse.eqv..False.,circular,energy_out)

!  steps = 0
!  open(1,file=traj_file_name)
!  do                       ! main loop is loop through timesteps
!    steps = steps + 1      ! number of timesteps 
!    ui = 10+steps          ! unit id of conf'steps'.dat
!    write(steps_string,frmt) steps
!    filename = 'conf'//trim(steps_string)//'.dat'
!    open(ui,file=filename)
!    if (energy_out.eqv..False.) then
!      do i=1,bp
!        read(1, *, iostat=io) x,y,z,bx,by,bz,nx,ny,nz,vx,vy,vz,lx,ly,lz
!        write(ui,*) x,y,z,bx,by,bz,nx,ny,nz,vx,vy,vz,lx,ly,lz
!      if (i0/=0) exit
!      end do
!    else
!      do i = 1,3
!        read(1,*, iostat=io)
!        if (io/=0) exit
!      end do
!      do i=1,bp
!        read(1, *, iostat=io) x,y,z,bx,by,bz,nx,ny,nz,vx,vy,vz,lx,ly,lz
!        write(ui,*) x,y,z,bx,by,bz,nx,ny,nz,vx,vy,vz,lx,ly,lz
!        if (io/=0) exit
!      end do
!    end if
!    close(ui)
!  end do
!  close(1)
! 




end program processor
