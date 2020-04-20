program processor
  use readers
  use twist_writhe_fortran
  use nemehunter
  implicit none
  

  integer :: i,nlines,bp,steps,ui,io,step,linesperfile,neme_pos,j,ier,k=3,npoints=1000
  character ::  init_conf_name*20,traj_file_name*20,final_conf_name*20, dir_name*4
  logical :: reverse=.True.,circular,energy_out
  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2,x,y,z,bx,by,bz,nx,ny,nz,vx,vy,vz,lx,ly,lz
  real, dimension(:), allocatable :: tx1,ty1,tz1,tx2,ty2,tz2,cx1,cy1,cz1,cx2,cy2,cz2,dum_x1
  integer :: nx1,ny1,nz1,nx2,ny2,nz2
  character :: frmt='(I5.5)', steps_string*10, cmdfmt='A', cmd*100, filename*50, lpf_string*10
  character :: search_string*4, energy_out_string*4, top_line*20
  character :: twist_writhe_file_name*20,plectoneme_file_name*20
  real :: twist_integral,writhe_integral

  search_string = 't = '



  ! read in the input files and create the output files
  twist_writhe_file_name = 'twist_writhe.dat'
  plectoneme_file_name = 'plectoneme.dat'

!  read(*,*) init_conf_name
  read(*,*) traj_file_name
  read(*,*) final_conf_name
 



  ! call reader on the initical config to get bp,circular and energy_out, this also allocates the memory
  ! for x1,...,c2z correctly
  call reader(final_conf_name,step,bp,x1,y1,z1,x2,y2,z2,reverse.eqv..False.,circular,energy_out)


  open(1,file=traj_file_name)
  open(2,file=twist_writhe_file_name)
  open(3,file=plectoneme_file_name)

  allocate(dum_x1(bp))
   do i=1,bp
    dum_x1(i) = i-1
  end do



  do                             ! main loop is loop through printed timesteps
    if (energy_out.eqv..False.) then
      do i = 1,bp
        read(1,*,iostat=io) x1(i), y1(i), z1(i)
        if (io/=0) exit
      end do
      do i = 1, bp
        if (reverse.eqv..True.) then
          read(1,*,iostat=io) x2(bp+1-i), y2(bp+1-i), z2(bp+1-i)
          if (io/=0) exit
        else
          read(1,*,iostat=io) x2(i), y2(i), z2(i)
          if (io/=0) exit
        end if
      end do
    else
      do i = 1,3
        read(1,"(A20)",iostat=io) top_line
        if (io/=0) exit
        read (top_line,'(A4)') energy_out_string
        if (energy_out_string==search_string) then
          read (top_line,"(A4,I8)") search_string,step
          write(*,*) step
        end if
      end do
      do i = 1,bp
        read(1,*,iostat=io) x1(i), y1(i), z1(i)
      end do
      do i = 1,bp
        if (reverse.eqv..True.) then
          read(1,*,iostat=io) x2(bp+1-i), y2(bp+1-i), z2(bp+1-i)
          if (io/=0) exit
        else
          read(1,*,iostat=io) x2(i), y2(i), z2(i)
          if (io/=0) exit
        end if
      end do
    end if
    

    ! now we have loaded in the coordimnates for the ith timestep interval, we can run any desired
    ! scripts
    ! calculate the splines of the x,y,z postions for both strands
    call get_spline(dum_x1,x1,tx1,cx1,k,nx1,bp,circular,ier)
    call get_spline(dum_x1,y1,ty1,cy1,k,ny1,bp,circular,ier)
    call get_spline(dum_x1,z1,tz1,cz1,k,nz1,bp,circular,ier)
    call get_spline(dum_x1,x2,tx2,cx2,k,nx2,bp,circular,ier)
    call get_spline(dum_x1,y2,ty2,cy2,k,ny2,bp,circular,ier)
    call get_spline(dum_x1,z2,tz2,cz2,k,nz2,bp,circular,ier)

    ! call main subroutines
    call twist_writhe(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2,&
                                            &cz2,nz2,circular,twist_integral,writhe_integral)
    call neme_hunter(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2&
                                              &,cz2,nz2,circular,neme_pos)

    write(2,*) step,twist_integral,writhe_integral
    write(3,*) step,neme_pos


    ! rather annoyingly we have to deallocate all the allocatbale spline objects, this is 
    ! unavoidable if we wish for the fitpack routines to be wrapped in the module spline
    deallocate(tx1)
    deallocate(ty1)
    deallocate(tz1)
    deallocate(tx2)
    deallocate(ty2)
    deallocate(tz2)
    deallocate(cx1)
    deallocate(cy1)
    deallocate(cz1)
    deallocate(cx2)
    deallocate(cy2)
    deallocate(cz2)
  end do

  close(1)
  close(2)
  close(3)


  return


end program processor
