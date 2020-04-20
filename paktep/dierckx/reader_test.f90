program reader_test
  use readers
  implicit none

  integer :: bp, k=3, npoints=10
  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2
  real, dimension(:), allocatable :: tx, cx
  character :: file_name*20

  call reader('conf.conf',bp, x1, y1, z1, x2, y2, z2)
  write(*,*) x1(1), y1(1), z1(1)
 end program reader_test
