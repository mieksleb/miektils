program main
  use processor
  use readers
  use twist_writhe_fortran
!  use nemehunter
  use persistence_length
  implicit none

  character :: sim_type*20,twist__writhe*20,plectoneme_pos*20,persistence__length*30,input_file_name*50
  character :: top_file_name*50,top*50,conf_file_name*50,conf*50
  integer :: sim_binary,bp,step 
  logical :: twist_writhe_log,plectoneme_pos_log,persistence_length_log,circular,energy_out,reverse
  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2

  read(*,*) input_file_name
  open(1,file=input_file_name)
  read(1,"(a12)") sim_type
  read(1,"(a50)") top
  read(1,"(a50)") conf
  read(1,"(a20)") twist__writhe
  read(1,"(a20)") plectoneme_pos
  read(1,"(a50)") persistence__length

  
  top_file_name = top(12:50)
  write(*,"(a22,a30)") 'Reading toplogy file: ',top_file_name
  conf_file_name = conf(8:50)
  write(*,"(a28,a30)") 'Reading configuration file: ',conf_file_name
  if (sim_type(1:12)=='type = oxdna') then
    sim_binary = 0
  else 
    sim_binary = 1
  end if


  if (twist__writhe=='twist_writhe = 1') then
    twist_writhe_log=.True.
  else
    twist_writhe_log=.False.
  end if

  if (plectoneme_pos=='plectoneme_pos = 1') then
    plectoneme_pos_log=.True.
  else
    plectoneme_pos_log=.False.
  end if
  
  if (persistence__length =='persistence_length = 1') then
    persistence_length_log=.True.
  else
    persistence_length_log=.False.
  end if

  ! we call the topology readers to get the number of base pairs and open/closed toplogy
  ! extra routine for energy_out logical required if oxdna 
  if (sim_binary.eq.0) then
    call oxdna_top_reader(top_file_name,bp,circular)
  else
    !call amber_top_reader(top_file_name,bp,circular)
    call amber_init_reader(conf_file_name,bp,circular)
  end if
  
  write(*,"(i5,a21)") bp," base pairs detected"

  if (circular.eqv..True.) then
    write(*,*) 'Circular DNA detected'
  else
    write(*,*) 'Linear DNA detected'
  end if

  ! we now call the main subroutine: process from the processor module
  call process(sim_binary,bp,conf_file_name,top_file_name,twist_writhe_log,plectoneme_pos_log,persistence_length_log,&
                            &circular,energy_out,reverse)



end program main
