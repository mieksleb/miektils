module spline
  use readers
  implicit none

  integer :: bp, k=3, ier, n, i, m
  real, dimension(:), allocatable :: x,x1,y1,z1,x2,y2,z2
  real, dimension(:), allocatable :: tx, cx
  character :: file_name*20
  logical :: circular, nox
   
  
  real, dimension(:), allocatable :: x

!  file_name = 'conf.conf'
!  call reader(file_name,bp, x1, y1, z1, x2, y2, z2, circular)
!  write(*,*) x1(1), y1(1), z1(1)
!  m = bp
!  allocate(x(m))
!  call get_spline(x,x1,tx,cx,k,bp,n,ier,nox=.True.)
!  if (ier.eq.10) then
!     print *,"something went wrong"
!  else
!  end if
  

!  do i=1,bp
!    x(i)=i-1
!  end do
!  call evaluate_spline(tx,n,cx,k,x,x1,m,ier)


contains



! input m data points (x,y) and returns splien object t,c,k,n. If data is periodic then set circular=.True.
! set nox=.True. if there is no supplied x data
subroutine get_spline(x,y,t,c,k,m,n,circular,ier,nox)

  real :: xb,xe,s=0,fp
  integer :: iopt=0,m,k,nest,n,ier,i !nest = m/2 will suffice
  integer :: lwrk
  logical :: circular,nox
  real, dimension(:), allocatable :: x,y,w
  real, dimension(:), allocatable :: t,c
  real, dimension(:), allocatable :: wrk
  integer, dimension(:), allocatable :: iwrk

  nest=m+k+1
  n=nest
  lwrk=m*(k+1)+nest*(7+3*k)
  allocate(x(m))
  allocate(w(m))
  allocate(t(nest))
  allocate(c(nest))
  allocate(wrk(lwrk))
  allocate(iwrk(nest))
  
  
  do i=1,m
    w(i)=1
  end do

  if nox=.True. then
    do i=1,m
      x(i)=i-1
    end do
  else 
  end if

  xb=x(1)
  xe=x(m)

  if(circular.eq..False.) then 
    call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
  else 
    call purcur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
  end if

end subroutine get_spline



! evaluate a spline object t,c,k,n at positions x and get values y. m is length of evalutation points. Allocates memeory of output here, so no need to do so when called
subroutine evaluate_spline(t,n,c,k,x,y,m,circular,ier,der)

  integer :: k,n,m,ier,der
  real, dimension(:), allocatable :: x,y
  real, dimension(:), allocatable :: t,c
  real, dimension(:), allocatable :: wrk
  integer, dimension(:), allocatable :: iwrk
  logical :: circular
  allocate(y(m))
 


  if(der.eq.0) then
    call splev(t,n,c,k,x,y,m,ier)
  else
    call splder(t,n,c,k,der,x,y,m,wrk,ier)
  end if

end subroutine evaluate_spline


end module spline

