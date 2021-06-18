module processor
  use readers
  use twist_writhe_fortran
  use nemehunter
  use persistence_length
  implicit none
  
contains

subroutine process(sim_type,bp,conf_name,top_name,twist_writhe_logical,plectoneme_position_logical,persistence_length_logical,&
                    &circular,energy_out,reverse)

  integer :: i,nlines,bp,steps,io,step,neme_pos,j,ier,k=3,npoints=1000,nsteps,array_len,q,sim_type
  character ::  init_conf_name*40,traj_file_name*40,final_conf_name*40,conf_name*40,top_name*40
  logical :: reverse,circular,energy_out
  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2,x,ss,ss_sum,x2_rev,y2_rev,z2_rev
  real, dimension(:), allocatable :: tx1,ty1,tz1,tx2,ty2,tz2,cx1,cy1,cz1,cx2,cy2,cz2,dum_x1,corr,corr_sum
  integer :: nx1,ny1,nz1,nx2,ny2,nz2,tot_steps
  character :: frmt='(I5.5)', steps_string*10,filename*50
  character :: search_string*4, energy_out_string*4, top_line*20
  character :: twist_writhe_file_name*40,plectoneme_file_name*40,am_string*20,top_file_name*40
  real :: twist_integral,writhe_integral,a,am,neme_len
  logical :: twist_writhe_logical,plectoneme_position_logical,persistence_length_logical
  integer :: bases,circ_int,dum_int,m,n,nres,natoms,s
  integer(kind = 4) :: bpn
  character :: atom_name*4,a1*1
  character :: base*1, second_line*20,resname*10,atom*4,word*10,line*20,string*128,altloc,icode,chains,segid*4,charge*2,element*2
  real :: tol=0.01,occ,temp


  search_string = 't = '
  reverse=.False.


  ! read in the input files and create the output files
  twist_writhe_file_name = 'twist_writhe.dat'
  plectoneme_file_name = 'plectoneme.dat'

  

  ! length of position arrays is one extra for circular dna such that x(1)=x(array_len)
  if (circular.eqv..True.) then
    array_len = bp+1
  else 
    array_len = bp
  end if


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

  allocate(dum_x1(array_len))
   do i=1,array_len
    dum_x1(i) = i-1
  end do


  ! open the trajectory file and find the number of steps
  open(1,file=conf_name)
  if (sim_type.eq.0) then ! for oxdna we annoyingly have to obtain the energy_out logical here

    ! searches for string pattern "t = " to set logical energy_out and obatin timestep
    read(1,"(A20)") top_line
    read(top_line,'(A4)') energy_out_string
    close(1)
    open(1,file=conf_name)
    if (energy_out_string==search_string) then
      read (1,"(A4,I8)") search_string,step
      energy_out=.True.
    else
      energy_out=.False.
    end if
    close(1)
    write(*,*) energy_out



    io=0
    tot_steps=0
    do while(io.eq.0)
      tot_steps=tot_steps+1        
      if (energy_out.eqv..False.) then
        do i = 1,2*bp ! in future will be strands*bp
          read(1,*,iostat=io)
        end do
      else
        do i = 1,3
          read(1,*,iostat=io)
        end do
        do i = 1,2*bp
          read(1,*,iostat=io)
        end do
      end if
    end do
  else
    io=0
    tot_steps=0
    do while(io.eq.0)
      read(1,"(a5)",iostat=io) line
      if (line(1:5)=='MODEL') then
        tot_steps=tot_steps+1
      end if
    end do
  end if
  write(*,"(i6,a14)") tot_steps, ' total frames'

  close(1)

  allocate(x1(array_len))
  allocate(y1(array_len))
  allocate(z1(array_len))
  allocate(x2(array_len))
  allocate(y2(array_len))
  allocate(z2(array_len))

  ! we now perform the main loop which is over timesteps
  open(1,file=conf_name)
  nsteps = 0
  do j=1,tot_steps-1
    nsteps = nsteps + 1
!    write(*,"(a5,i4)") 'step ',nsteps
    call progress_bar(nsteps,tot_steps)
    if (sim_type.eq.1) then
      read(1, "(a12,i6)") string, step
      if (string(1:5)=='MODEL') then
        call amber_block_reader(1,step,bp,x1,y1,z1,x2,y2,z2,.True.,circular)
      end if
    end if
    ! now we can apply the routines
    if (twist_writhe_logical.eqv..True.) then
      call twist_writhe_xyz(bp,x1,y1,z1,x2,y2,z2,twist_integral,writhe_integral,circular)
      write(21,*) nsteps,twist_integral,writhe_integral
    end if
  end do


  
    ! now we have loaded in the coordinates for the ith timestep interval, we can run any desired
    ! scripts
    
    ! call main subroutines of desired quantities


    ! twist and writhe are calculated from the module TWIST_WRITHE_FORTRAN
!    if (twist_writhe_logical.eqv..True.) then
!      call twist_writhe_xyz(bp,x1,y1,z1,x2,y2,z2,twist_integral,writhe_integral,circular)
!      write(21,*) step,twist_integral,writhe_integral
!    end if

    ! plectoneme position is calculated from the module NEMEHUNTER
!    if (plectoneme_position_logical.eqv..True.) then
!      write(22,*) step,neme_pos,neme_len
!    end if

    ! persistence length is a property of the entire trajectory and correlations are appended to the global variable ss_sum   
!    if (persistence_length_logical.eqv..True.) then
    
    ! we need to average the tangent tangent correlations over all timesteps
!      ss_sum(:) = ss_sum(:) + ss(:)
!      corr_sum(:) = corr_sum(:)+corr(:)
!    end if

!  end do



  ! persistence length is then calculated by finding the time average of the tangent-tangent correlations
  ! and applying an exponential fit, routines are called from the module PERSISTENCE_LENGTH
!  if (persistence_length_logical.eqv..True.) then
!    corr_sum(:) = corr_sum(:)/nsteps
!    ss_sum(:) = ss_sum(:)/nsteps
!    a = exp_fit(ss_sum,corr_sum,npoints)
!    a = -1/a
!    write(*,*) a
!    am = 8.518E-1*a
!    write(*,*) am
!    write(am_string,"(F2.7)") am
!    write(*,*) am_string
!    write(*,*) "Persistence Length = "//adjustl(trim(am_string))//"nm"
!  end if


  ! close all files
  close(1)
  close(21)
  close(22)


end subroutine process

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
    character :: x*5

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


end module processor
