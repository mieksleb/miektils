program spline
  use readers
  implicit none

  integer :: bp, k=3, ier
  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2
  real, dimension(:), allocatable :: tx, cx
  character :: file_name*20
   
  

  file_name = 'conf.conf'
  call reader(file_name,bp, x1, y1, z1, x2, y2, z2)
  write(*,*) x1(1), y1(1), z1(1)
  call get_spline(x1,tx,cx,k,bp,ier)
  if (ier.eq.10) then
     print *,"something went wrong"
  else
    write(*,*) "wooopee"
  end if


contains

subroutine get_spline(y,t,c,k,m,ier)

  real :: xb,xe,s=0,fp
  integer :: iopt=0,m,k,nest,n,ier,i !nest = m/2 will suffice
  integer :: lwrk
  real, dimension(:), allocatable :: x,y,w
  real, dimension(:), allocatable :: t,c
  real, dimension(:), allocatable :: wrk
  integer, dimension(:), allocatable :: iwrk

! test parameters
!  integer :: ia,ib,ifp,ig,iq,iz,j,k1,k2,lwest,nmin,maxit,l,nk1,nk2,nk3,p
!  real :: tol,tj,tl


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

  do i=1,m
    x(i)=i-1
  end do

!  do i=1,nest
!    t(i)=i
!  end do

  xb=x(1)
  xe=x(m)
  call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)


!  before starting computations a data check is made. if the input data
!  are invalid, control is immediately repassed to the calling program.
!      ier = 10
!      if(k.le.0 .or. k.gt.5) go to 50
!      k1 = k+1
!      k2 = k1+1
!      write(*,*) "snoop"
!      if(iopt.lt.(-1) .or. iopt.gt.1) go to 50
!      nmin = 2*k1
!      write(*,*) "skidoosh"
!      if(m.lt.k1 .or. nest.lt.nmin) go to 50
!      lwest = m*k1+nest*(7+3*k)
!      write(*,*) "bleep"
!      if(lwrk.lt.lwest) go to 50
!      if(xb.gt.x(1) .or. xe.lt.x(m) .or. w(1).le.0.) go to 50
!      write(*,*) "bloop"
!      do 10 i=2,m
!         if(x(i-1).ge.x(i) .or. w(i).le.0.) go to 50
!  10  continue
!      if(iopt.ge.0) go to 30
!      if(n.lt.nmin .or. n.gt.nest) go to 50
!      j = n
!      do 20 i=1,k1
!         t(i) = xb
!         t(j) = xe
!         j = j-1
!      write(*,*) "sneep"
!  20  continue
!      call fpchec(x,m,t,n,k,ier)
!      if(ier) 40,50,40
!  30  if(s.lt.0.) go to 40
!      write(*,*) "a cauldron"
!      if(s.eq.0. .and. nest.lt.(m+k1)) go to 40
!      ier = 0
! we partition the working space and determine the spline approximation.
!  40  ifp = 1
!      iz = ifp+nest
!      ia = iz+nest
!      ib = ia+nest*k1
!      ig = ib+nest*k2
!      iq = ig+nest*k2
!      call fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,n,t,c,fp,wrk(ifp),wrk(iz),wrk(ia),wrk(ib),wrk(ig),wrk(iq),iwrk,ier)
!  50 continue
!      k1 = k+1
!      k2 = k1+1
!      nk1 = n-k1
!      nk2 = nk1+1
!      ier = 10
!  check condition no 1
!      write(*,*) "gammo"
!      if(nk1.lt.k1 .or. nk1.gt.m) go to 120
!      write(*,*) "bumblebess"
!  check condition no 2
!      j = n
!      do 60 p=1,k
!        if(t(p).gt.t(p+1)) go to 120
!          write(*,*) "pretty"
!        if(t(j).lt.t(j-1)) go to 120
!        write(*,*) "pony"
!        j = j-1
!  60  continue
!  check condition no 3
!      write(*,*) "humdinger"
!      do 70 p=k2,nk2
!        if(t(p).le.t(p-1)) go to 70
!  70  continue
!  check condition no 4
!      write(*,*) "big"
!      if(x(1).lt.t(k1) .or. x(m).gt.t(nk2)) go to 120
!      write(*,*) "great"
!  check condition no5
!      if(x(m).le.t(nk1)) go to 120
!      write(*,*) "fat"
!      write(*,*)  t(1), t(2), t(3), t(4), t(5), t(6)
!      write(*,*)  t(n-6), t(n-5), t(n-4), t(n-3), t(n-2), t(n-1), t(n)
!      if(x(1).ge.t(k2)) go to 120
!      write(*,*) "sin duda"
!      p = 1
!      l = k2
!      nk3 = nk1-1
!      write(*,*) "wowee"
!      if(nk3.lt.2) go to 110
!      do 100 j=2,nk3
!        write(*,*) j
!        tj = t(j)
!        l = l+1
!        tl = t(l)
!  80    p = p+1
!        write(*,*) "chow"
!        if(p.ge.m) go to 120
!        if(x(p).le.tj) go to 80
!        write(*,*) "chunder"
!        if(x(p).ge.tl) go to 120
!        write(*,*) "chunder baby"
!  100  continue
!  110  ier = 0
!  120  return
    
 
end subroutine get_spline


end program spline

