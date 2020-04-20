program reader_test
  use readers
  implicit none

  integer :: bp, k=3,step
  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2
  real, dimension(:), allocatable :: tx, cx
  character :: file_name*20
  logical :: circular,energy_out,reverse


  call reader('conf.conf',step,bp, x1, y1, z1, x2, y2, z2, reverse.eqv..False.,circular,energy_out)
  write(*,*) x1(1), y1(1), z1(1)
  write(*,*) circular
  write(*,*) energy_out
  write(*,*) step
end program reader_test
