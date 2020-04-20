module nemehunter
  use readers
  use spline
  use geom
  implicit none

contains 

subroutine neme_hunter_main(file_name,neme_pos)
 
  integer :: neme_pos
  character :: file_name*20


  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2,tx1,ty1,tz1,tx2,ty2,tz2,cx1,cy1,cz1,cx2,cy2,cz2
  real, dimension(:), allocatable :: dum_x1
  integer :: bp,npoints=100,nx1,ny1,nz1,nx2,ny2,nz2,k=3,ier,i,step
  logical :: circular,reverse,energy_out
  
  ! call reader to load in positions
  call reader(file_name,step,bp,x1,y1,z1,x2,y2,z2,reverse.eqv..False.,circular,energy_out)

  ! generate a linear sequence to feed as independent variable in spline fitting procedure
  allocate(dum_x1(bp))
  do i=1,bp
    dum_x1(i)=i-1
  end do

  ! calculate the splines of the x,y,z postions for both strands
  call get_spline(dum_x1,x1,tx1,cx1,k,nx1,bp,circular,ier)
  call get_spline(dum_x1,y1,ty1,cy1,k,ny1,bp,circular,ier)
  call get_spline(dum_x1,z1,tz1,cz1,k,nz1,bp,circular,ier)
  call get_spline(dum_x1,x2,tx2,cx2,k,nx2,bp,circular,ier)
  call get_spline(dum_x1,y2,ty2,cy2,k,ny2,bp,circular,ier)
  call get_spline(dum_x1,z2,tz2,cz2,k,nz2,bp,circular,ier)
  
  call neme_hunter(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2&
                                              &,cz2,nz2,circular,neme_pos)

end subroutine neme_hunter_main

subroutine neme_hunter(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2&
                                              &,cz2,nz2,circular,neme_pos)

  real, dimension(:), allocatable   :: x1,y1,z1,x2,y2,z2,tx1,ty1,tz1,tx2,ty2,tz2,cx1,cy1,cz1,cx2,cy2,cz2
  integer                           :: ii,jj,srange,ier,k=3,nx1,ny1,nz1,nx2,ny2,nz2,bp,nuxx,nuyy,nuzz,intersection,ival
  integer                           :: npoints,m,nnuxx,nnuyy,nnuzz,nsx,nsy,nsz,neme_pos,i_end_loop,j_end_loop,diff_ij,i,j
  logical                           :: circular,reverse
  real                              :: bpinc,delta_s
  real, dimension(:), allocatable   :: msxx,msyy,mszz
  real, dimension(:,:), allocatable :: uu,tt,boyo
  real, dimension(:), allocatable   :: uxx_bpi,uyy_bpi,uzz_bpi
  real, dimension(:), allocatable   :: x,y,z,m1xx,m1yy,m1zz,dmxx,dmyy,dmzz,duxx,duyy,duzz,xx,yy,zz
  real, dimension(:), allocatable   :: diff,tuxx,tuyy,tuzz,cnuxx,cnuyy,cnuzz,csx,csy,csz,cuxx,cuyy,cuzz
  real, dimension(:), allocatable   :: tnuxx,tnuyy,tnuzz,tsx,tsy,tsz,ss,contour,bpi,mxx,myy,mzz
  real, dimension(:), allocatable   :: sx1,sy1,sz1,sx2,sy2,sz2,snuxx,snuyy,snuzz,dum_x2,curvature
  real :: x_spread,y_spread,z_spread, tol=0.01
  logical :: x_proj,y_proj,z_proj
  integer :: n, flag, max_curv,coarse_factor=5,npoints_coarse,mid,neme_pos2,cutoff=20
  real(kind=8) :: x5,y5,alpha=0.0


  ! we begin with a much coarser spline fit as the self intersection algorith scales quadratically
  npoints_coarse = int(npoints/coarse_factor)
  m = npoints_coarse

  allocate(contour(npoints_coarse))
  allocate(ss(npoints_coarse))
  allocate(bpi(npoints_coarse))
  allocate(sx1(npoints_coarse))
  allocate(sy1(npoints_coarse))
  allocate(sz1(npoints_coarse))
  allocate(sx2(npoints_coarse))
  allocate(sy2(npoints_coarse))
  allocate(sz2(npoints_coarse))

! calculate the base pair index (bpi) common to both strands
  bpinc = (real(bp,8)-1)/(real(npoints_coarse,8)-1)
  do ii=1,npoints_coarse
    bpi(ii)=(ii-1)*bpinc
  end do
  
  
 ! evaluate the splines in bpi, output is sij for i=x,y,z j=1,2
  call evaluate_spline(bpi,sx1,tx1,cx1,k,nx1,npoints_coarse,circular,ier,0)
  call evaluate_spline(bpi,sy1,ty1,cy1,k,ny1,npoints_coarse,circular,ier,0)
  call evaluate_spline(bpi,sz1,tz1,cz1,k,nz1,npoints_coarse,circular,ier,0)
  call evaluate_spline(bpi,sx2,tx2,cx2,k,nx2,npoints_coarse,circular,ier,0)
  call evaluate_spline(bpi,sy2,ty2,cy2,k,ny2,npoints_coarse,circular,ier,0)
  call evaluate_spline(bpi,sz2,tz2,cz2,k,nz2,npoints_coarse,circular,ier,0)
 ! calculate the midpoint splines in bpi
  m1xx = (sx1+sx2)/2
  m1yy = (sy1+sy2)/2
  m1zz = (sz1+sz2)/2

  ! now contruct contour_len, a sequence which passes through all base pairs of lenght npoints
  contour(1)=0
  do ii=2,npoints_coarse
    contour(ii) = contour(ii-1) + ((m1xx(ii)-m1xx(ii-1))**2+(m1yy(ii)-m1yy(ii-1))**2+(m1zz(ii)-m1zz(ii-1))**2)**0.5
  end do
  delta_s = contour(npoints_coarse)/(real(npoints_coarse,8)-1)

  ss(1)=0
  do ii=2,npoints_coarse
    ss(ii) = delta_s*(ii-1)
  end do

  allocate(dum_x2(npoints_coarse))
  do ii=1,npoints_coarse
    dum_x2(ii)=ii-1
  end do


  ! now we compute the spline objects t,c,k,n of the midpoint spline 
  call get_spline(contour,m1xx,tsx,csx,k,nsx,npoints_coarse,circular,ier)
  call get_spline(contour,m1yy,tsy,csy,k,nsy,npoints_coarse,circular,ier)
  call get_spline(contour,m1zz,tsz,csz,k,nsz,npoints_coarse,circular,ier)

  allocate(xx(npoints_coarse))
  allocate(yy(npoints_coarse))
  allocate(zz(npoints_coarse))
  allocate(dmxx(npoints_coarse))
  allocate(dmyy(npoints_coarse))
  allocate(dmzz(npoints_coarse))
  ! we now evaluate this midpoint spline and its derivatives
  call evaluate_spline(ss,xx,tsx,csx,k,nsx,npoints_coarse,circular,ier,0)
  call evaluate_spline(ss,yy,tsy,csy,k,nsy,npoints_coarse,circular,ier,0)
  call evaluate_spline(ss,zz,tsz,csz,k,nsz,npoints_coarse,circular,ier,0)
  call evaluate_spline(ss,dmxx,tsx,csx,k,nsx,npoints_coarse,circular,ier,1)
  call evaluate_spline(ss,dmyy,tsy,csy,k,nsy,npoints_coarse,circular,ier,1)
  call evaluate_spline(ss,dmzz,tsz,csz,k,nsz,npoints_coarse,circular,ier,1)


  ! we now find which axis we have the most spread and choose to project out this dimension to obtain a 2D spline
  x_spread = maxval(xx)-minval(xx)
  y_spread = maxval(yy)-minval(yy)
  z_spread = maxval(zz)-minval(zz)

  if (x_spread.le.y_spread.and.x_spread.le.z_spread) then
    x_proj = .True.
    y_proj = .False.
    z_proj = .False.
  else if (y_spread.le.x_spread.and.y_spread.le.z_spread) then
    x_proj = .False. 
    y_proj = .True.
    z_proj = .False.
    
  else if(z_spread.le.x_spread.and.z_spread.le.y_spread) then
    x_proj = .False.
    y_proj = .False.
    z_proj = .True.
  end if

    ! now we project out the dimension with the smallest range
    ! we now call a routine from DUTCH which determines if our polygon has any self intersections (scales quadratically)
  do i=1,npoints_coarse-1
    do j=1,npoints_coarse-1
      if (j>i+1) then
        
        if (x_proj.eqv..True.) then
          call LINES_SEG_INT_2D(real(yy(i),8), real(zz(i),8), real(yy(i+1),8),real(zz(i+1),8),&
                      & real(yy(j),8),real(zz(j),8),real(yy(j+1),8),real(zz(j+1),8), flag,x5,y5)
          if (flag==1) then
!            write(*,*) i,j
            diff_ij=npoints_coarse  ! need a large upper bound for diff_ij to begin iterations
            if ((j-i).le.diff_ij.and.(j-i).ge.cutoff) then ! if the differnece between i and j is smaller than the previous, then set i and j to be the new end loop indices
              diff_ij=j-i
              i_end_loop = i
              j_end_loop = j
            else
            end if
          else
          end if

        else if(y_proj.eqv..True.) then
          call LINES_SEG_INT_2D(real(xx(i),8), real(zz(i),8), real(xx(i+1),8),real(zz(i+1),8),&
                      & real(xx(j),8),real(zz(j),8),real(xx(j+1),8),real(zz(j+1),8), flag,x5,y5)
          if (flag==1) then
!            write(*,*) i,j
            diff_ij=npoints_coarse  ! need a large upper bound for diff_ij to begin iterations
            if ((j-i).le.diff_ij.and.(j-i).ge.cutoff) then ! if the differnece between i and j is smaller than the previous, then set i and j to be the new end loop indices
              diff_ij=j-i
              i_end_loop = i
              j_end_loop = j
            else
            end if
          else
          end if

        else if (z_proj.eqv..True.) then
          call LINES_SEG_INT_2D(real(xx(i),8), real(yy(i),8), real(xx(i+1),8),real(yy(i+1),8),&
                      & real(xx(j),8),real(yy(j),8),real(xx(j+1),8),real(yy(j+1),8), flag,x5,y5)
          if (flag==1) then
!            write(*,*) i,j
            diff_ij=npoints_coarse  ! need a large upper bound for diff_ij to begin iterations
            if ((j-i).le.diff_ij.and.(j-i).ge.cutoff) then ! if the differnece between i and j is smaller than the previous, then set i and j to be the new end loop indices
              diff_ij=j-i
              i_end_loop = i
              j_end_loop = j
            else
            end if
          else
          end if
        else
        end if
      else
      end if
    end do
  end do
   
!  we then convert these indices to
!  i_end_loop=i_end_loop*coarse_factor
!  j_end_loop=j_end_loop*coarse_factor
  mid = int((i_end_loop + j_end_loop)/2)
  allocate(curvature(npoints_coarse))
  do i=1,npoints_coarse
    if (i.le.i_end_loop) then
      curvature(i) = 0
    else if (i.ge.j_end_loop) then
      curvature(i)=0
    else
      curvature(i) = (dmxx(i)**2+dmyy(i)**2+dmzz(i)**2)*exp(alpha*abs(i-mid))
      curvature(1)=0
      curvature(npoints_coarse)=0
    end if
  end do
  max_curv = maxloc(curvature,dim=1)
  neme_pos=int(bpi(max_curv))
  neme_pos2=int(bpi(mid))
!  write(*,*) neme_pos
end subroutine neme_hunter







end module nemehunter
