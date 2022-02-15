module bend
  implicit none

contains

subroutine get_angle(bp,x1,y1,z1,x2,y2,z2)

! finds the angle between two sets of points about an 

  integer :: hinge = 20,i,bp,num
  real, dimension(:), allocatable :: x_av1,y_av1,z_av1,x_av2,y_av2,z_av2,av1,av2,x1,y1,z1,x2,y2,z2
  real, dimension (:,:), allocatable :: cov1, cov2,mat1,mat2
  real :: xbar1,ybar1,zbar1,xbar2,ybar2,zbar2


  num = bp/2  

  allocate (x_av1(num))
  allocate (y_av1(num))
  allocate (z_av1(num))
  allocate (x_av2(num))
  allocate (y_av2(num))
  allocate (z_av2(num))
  allocate (av1(3))
  allocate (av2(3))

  allocate (mat1(num,3))
  allocate (mat2(num,3))
  allocate (cov1(num,3))
  allocate (cov2(num,3))

  ! get the midpoints positions
  do i = 1,bp
    if (i < hinge) then
      x_av1(i) = (x1(i) + x2(i))/2
      y_av1(i) = (y1(i) + y2(i))/2
      z_av1(i) = (z1(i) + z2(i))/2
    else
      x_av2(i) = (x1(i) + x2(i))/2
      y_av2(i) = (y1(i) + y2(i))/2
      z_av2(i) = (z1(i) + z2(i))/2
    end if
  end do



  xbar1 = 0
  ybar1 = 0
  zbar1 = 0
  xbar2 = 0
  ybar2 = 0
  zbar2 = 0
  do i=1,num
    xbar1 = xbar1 + x_av1(i)
    ybar1 = ybar1 + y_av1(i)
    zbar1 = zbar1 + z_av1(i)
    xbar2 = xbar2 + x_av2(i)
    ybar2 = ybar2 + y_av2(i)
    zbar2 = zbar2 + z_av2(i)
  end do

  xbar1 = xbar1/num
  ybar1 = ybar1/num
  zbar1 = zbar1/num
  xbar2 = xbar2/num
  ybar2 = ybar2/num
  zbar2 = zbar2/num
 
  av1(1) = xbar1
  av1(2) = ybar1
  av1(3) = zbar1
  av2(1) = xbar2
  av2(2) = ybar2
  av2(3) = zbar2


    

  do i=1,num
    mat1(i,1) = x_av1(i)
    mat1(i,2) = y_av1(i)
    mat1(i,3) = z_av1(i) 
    cov1(i,1) = mat1(i,1) - av1(1)
    cov1(i,2) = mat1(i,2) - av1(2)
    cov1(i,3) = mat1(i,3) - av1(3)
  end do

  !print *, cov1

 
  !do i = 1,num
  !  do j = 1,3
  !    cov1(i,j) = mat1(i,j)*mat1(j,i)







end subroutine get_angle


end module bend
