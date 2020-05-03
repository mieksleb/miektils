program processor
  use readers
  use twist_writhe_fortran
  use nemehunter
  use persistence_length
  implicit none
  

  integer :: i,nlines,bp,steps,io,step,neme_pos,j,ier,k=3,npoints=1000,nsteps
  character ::  init_conf_name*20,traj_file_name*20,final_conf_name*20 
  logical :: reverse=.True.,circular,energy_out
  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2,x,ss,ss_sum
  real, dimension(:), allocatable :: tx1,ty1,tz1,tx2,ty2,tz2,cx1,cy1,cz1,cx2,cy2,cz2,dum_x1,corr,corr_sum
  integer :: nx1,ny1,nz1,nx2,ny2,nz2,tot_steps
  character :: frmt='(I5.5)', steps_string*10,filename*50
  character :: search_string*4, energy_out_string*4, top_line*20
  character :: twist_writhe_file_name*20,plectoneme_file_name*20,am_string*20
  real :: twist_integral,writhe_integral,a,am,neme_len
  logical :: twist_writhe_logical,plectoneme_position_logical,persistence_length_logical

  twist_writhe_logical=.True.
  plectoneme_position_logical=.True.
  persistence_length_logical=.True.

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


  ! we now open any required output files for writing and perform any other pre process tasks

  if (twist_writhe_logical.eqv..True.) then
    open(21,file=twist_writhe_file_name)
  end if
  
  if (plectoneme_position_logical.eqv..True.) then
    open(22,file=plectoneme_file_name)
  end if

  if (persistence_length_logical.eqv..True.) then
    allocate(corr_sum(npoints))
    allocate(corr(npoints))
    allocate(ss_sum(npoints))
    corr_sum=0
    ss_sum=0
  end if

  allocate(dum_x1(bp))
   do i=1,bp
    dum_x1(i) = i-1
  end do

  nsteps = 0
  ! we now perform the main loop which is over timesteps
  do j=1,100                      
    nsteps = nsteps + 1
    call progress_bar(nsteps,100)
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
    ! call main subroutines of desired quantities


    if (plectoneme_position_logical.eqv..True.) then
      call neme_hunter(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2&
                                             &,cz2,nz2,circular,neme_pos,neme_len)
      write(22,*) step,neme_pos,neme_len
    end if


    if (twist_writhe_logical.eqv..True.) then
      call twist_writhe(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2,&
                                            &cz2,nz2,circular,twist_integral,writhe_integral)
      write(21,*) step,twist_integral,writhe_integral
    end if

!      write(*,*) ""


    
    if (persistence_length_logical.eqv..True.) then
      call tangent_correlation(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2&
                                             &,cz2,nz2,circular,corr,ss)
    ! we need to average the tangent tangent correlations over all timesteps
      ss_sum(:) = ss_sum(:) + ss(:)
      corr_sum(:) = corr_sum(:)+corr(:)
    end if


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
    deallocate(ss)
  end do

  if (persistence_length_logical.eqv..True.) then
    corr_sum(:) = corr_sum(:)/nsteps
    ss_sum(:) = ss_sum(:)/nsteps
    a = exp_fit(ss_sum,corr_sum,npoints)
    a = -1/a
    write(*,*) a
    am = 8.518E-1*a
    write(*,*) am
    write(am_string,"(F2.7)") am
    write(*,*) am_string
    write(*,*) "Persistence Length = "//adjustl(trim(am_string))//"nm"
  end if



  close(1)
  close(21)
  close(22)


contains


subroutine progress_bar(iteration, maximum)
!
! Prints progress bar.
!
! Args: 
!     iteration - iteration number
!     maximum - total iterations
!
    implicit none
    integer :: iteration, maximum
    integer :: counter
    integer :: step, done

    step = nint(iteration * 100 / (1.0 * maximum))
    done = floor(step / 10.0)  ! mark every 10%

    do counter = 1, 36                    ! clear whole line - 36 chars
        write(6,'(a)',advance='no') '\b'  ! (\b - backslash)
    end do

    write(6,'(a)',advance='no') ' -> In progress... ['
    if (done .LE. 0) then
        do counter = 1, 10
            write(6,'(a)',advance='no') '='
        end do
    else if ((done .GT. 0) .and. (done .LT. 10)) then
        do counter = 1, done
            write(6,'(a)',advance='no') '#'
        end do
        do counter = done+1, 10
            write(6,'(a)',advance='no') '='
        end do 
    else
        do counter = 1, 10
            write(6,'(a)',advance='no') '#'
        end do
    end if
    write(6,'(a)',advance='no') '] '
    write(6,'(I3.1)',advance='no') step
    write(6,'(a)',advance='no') '%'
end


subroutine remaining_time(iteration, maximum)
!
! Prints remaining time.
!
! Args: 
!     iteration - iteration number
!     maximum - total iterations
!
    implicit none
    integer :: iteration, maximum
    integer :: counter
    integer :: hours, minutes, seconds
    real :: tarray(2), current, remains

    call etime(tarray, current)

    remains = 2 * current * (maximum / (1.0 * iteration) - 1)
    hours = floor(remains / 3600)
    minutes = floor((remains - hours * 3600) / 60)
    seconds = nint(remains - (hours * 3600 + minutes * 60))

    do counter = 1, 38                    ! clear whole line - 38 chars
        write(6,'(a)',advance='no') '\b'  ! (\b - backslash)
    end do

    write(6,'(a)',advance='no') ' -> Remaining time (h:m:s): '
    write(6,'(I4.1)',advance='no') hours
    write(6,'(a)',advance='no') ':'
    write(6,'(I2.2)',advance='no') minutes
    write(6,'(a)',advance='no') ':'
    write(6,'(I2.2)',advance='no') seconds
end


subroutine progress_bar_time(iteration, maximum)
!
! Prints progress bar with remaining time.
!
! Args: 
!     iteration - iteration number
!     maximum - total iterations
!
    implicit none
    integer :: iteration, maximum
    integer :: counter
    integer :: hours, minutes, seconds
    integer :: step, done
    real :: tarray(2), current, remains

    call etime(tarray, current)
    
    remains = 2 * current * (maximum / (1.0 * iteration) - 1)
    hours = floor(remains / 3600)
    minutes = floor((remains - hours * 3600) / 60)
    seconds = nint(remains - (hours * 3600 + minutes * 60))

    step = nint(100/(1.0*maximum)*iteration)
    done = floor(step / 10.0)  ! mark every 10%

    do counter = 1, 63                    ! clear whole line - 63 chars
        write(6,'(a)',advance='no') '\b'  ! (\b - backslash)
    end do  

    write(6,'(a)',advance='no') ' -> In progress... ['
    if (done .LE. 0) then
        do counter = 1, 10
            write(6,'(a)',advance='no') '='
        end do
    else if ((done .GT. 0) .and. (done .LT. 10)) then
        do counter = 1, done
            write(6,'(a)',advance='no') '#'
        end do
        do counter = done+1, 10
            write(6,'(a)',advance='no') '='
        end do  
    else
        do counter = 1, 10
            write(6,'(a)',advance='no') '#'
        end do
    end if
    write(6,'(a)',advance='no') '] '
    write(6,'(I3.1)',advance='no') x
    write(6,'(a)',advance='no') '% (remaining '
    write(6,'(I4.1)',advance='no') hours
    write(6,'(a)',advance='no') ':'
    write(6,'(I2.2)',advance='no') minutes
    write(6,'(a)',advance='no') "'"
    write(6,'(I2.2)',advance='no') seconds
    write(6,'(a)',advance='no') '")'
end


end program processor
