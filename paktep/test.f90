program test
  use readers
  use twist_writhe_fortran
  implicit none

  integer :: bp,k=3,step,stepa,bpa
  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2
  real, dimension(:), allocatable :: x1a,y1a,z1a,x2a,y2a,z2a
  character :: file_name*50, oxdna_conf_name*50, oxdna_top_name*50, amber_conf_name*50, amber_top_name*50
  logical :: circular,energy_out,reverse=.False.,circulara,reversea
  real :: twist_integral,writhe_integral,twist_integrala,writhe_integrala


  oxdna_conf_name = "test/circle_300-Lk7.dat"
  oxdna_top_name = "test/circle_300-Lk7.top"
  amber_conf_name = "test/dv150t13.pdb"
  amber_top_name = "test/dv150t13.top"
  
  write(*,*) "Performing oxDNA test on 300bp minicircle with deltaLk=7"
  write(*,*) ""
  write(*,*) "Lk0 = 300/10.36 = 29 and so Lk = Tw0 = 36. Writhe should be zero as minicircle is planar"
  call oxdna_reader(oxdna_conf_name,oxdna_top_name,step,bp,x1,y1,z1,x2,y2,z2,circular,energy_out,reverse=.True.)
  call twist_writhe_xyz(bp,x1,y1,z1,x2,y2,z2,twist_integral,writhe_integral,circular)
  write(*,"(a8,f8.5)") "Twist = ",twist_integral
  write(*,"(a9,f8.5)") "Writhe = ",writhe_integral

  write(*,*) "Performing amber test on 336bp minicircle with Lk=30"
  call amber_reader(amber_conf_name,amber_top_name,step,bpa,x1a,y1a,z1a,x2a,y2a,z2a,reversea,circulara)
  call twist_writhe_xyz(bpa,x1a,y1a,z1a,x2a,y2a,z2a,twist_integrala,writhe_integrala,circulara)
  write(*,"(a8,f8.5)") "Twist = ",twist_integrala
  write(*,"(a9,f8.5)") "Writhe = ",writhe_integrala

end program test
