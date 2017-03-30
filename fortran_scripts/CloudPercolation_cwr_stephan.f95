module m_percolations_float_stephan
implicit none
contains
subroutine perc(cluster_strength,L, nclouds, rA, rA_ring, cloud_field, cloud_centers, gauss_field_out)
  real,intent(in) :: rA(:), rA_ring(:)
  integer,intent(in) :: cluster_strength, L, nclouds

  real :: p0, A0, A1, dx, dy, distance_center2, distance_x, distance_y, value_func, value_field, val_center, offset
  real, dimension(2*nclouds) :: random
  integer,dimension(L,L) :: m0,m1,m12,m2
  integer,dimension(nclouds,2) :: center
  real,   dimension(nclouds,2) :: random_dx  
  real,   dimension(L,L) :: gauss_field
  real,   dimension(L,L), intent(out) :: gauss_field_out
  integer,dimension(nclouds-1,2),intent(out) :: cloud_centers
  integer,dimension(L,L),intent(out) :: cloud_field
  integer :: k1,k2,num,area_lim_1_max,area_lim_1_min,area_lim_2_max,area_lim_2_min,random_int, ini_ran, ini_ran2
  integer(kind=8) :: u,v,w,r,clock
  
  val_center=10.0 !value of gauss_field at cloud center (decreases linearly to offset at r=r_disc)
  offset=8.0
  
  !m fields indicate presence or absence of clouds or cloud rings

  m0(:,:) = 1   !m0 [x,y]=1 means no cloud at [x,y]
  m2(:,:) = 0   !m2 [x,y]=1 means cloud at [x,y]
  m12(:,:)= 0   !m12[x,y]=1 means cloud or cloud ring at [x,y]

  
  center(:,:)=1 !Initilize 2D array for cloud center coordinates
  gauss_field(:,:)=0 !Initilize cloud field which has a maximum at cloud center and decreases towards the edges
  
  call system_clock(count=clock)
  call initran(clock,u,v,w)
  
  !draw all random numbers necessary for drawing cloud center indices
  
  do ini_ran=1,2*nclouds
    call getran(r,u,v,w)
    random(ini_ran)=i2f(r)
  end do
  
  !draw all random numbers to get a position of the cloud center within a grid-box [dx, dy]
  !for example if [i,j] are indices of grid-box within which the cloud center is located, [i+dx, i+dy] is the position of the cloud center

  do ini_ran2=1,nclouds
    call getran(r,u,v,w)
    random_dx(ini_ran2,1)=i2f(r)
    call getran(r,u,v,w)
    random_dx(ini_ran2,2)=i2f(r)
  end do  
  
  !Distribute first center randomly
  
  center(1,1)=floor(L*random(1))+1 !choose one from interval {n,n+1,...,m} with equal probability
  center(1,2)=floor(L*random(2))+1
 
 
  do num=2,nclouds
  
    !Search around center to asign supressed, increased and neutral state
    
    !Start by defining a box around the cloud center within which the new cloud can have an impact: 
    !a box around the cloud center with side-lengths of 2*radius*fac_cw
    
    area_lim_1_min=max(center(num-1,1)-ceiling(rA_ring(num-1)),1)
    area_lim_1_max=min(center(num-1,1)+ceiling(rA_ring(num-1)),L)

    area_lim_2_min=max(center(num-1,2)-ceiling(rA_ring(num-1)),1)
    area_lim_2_max=min(center(num-1,2)+ceiling(rA_ring(num-1)),L)
    
    dx=random_dx(num-1,1)-0.5
    dy=random_dx(num-1,2)-0.5

    !identify within this box distances to the center smaller than radius --> cloud, smaller than fac_cw*radius --> cloud ring
    do k1=area_lim_1_min,area_lim_1_max
        do k2=area_lim_2_min, area_lim_2_max
        
        distance_x=(real(k1-center(num-1,1))-dx)
        distance_y=(real(k2-center(num-1,2))-dy)
        
        distance_center2=distance_x*distance_x+distance_y*distance_y

        if (distance_center2<=rA(num-1)*rA(num-1)) then
            m2(k1,k2)=1
            m12(k1,k2)=1
            value_func=f_r(sqrt(distance_center2), real(rA(num-1)), val_center, offset)
            value_field=gauss_field(k1,k2)
            gauss_field(k1,k2)=max(value_field,value_func)
            
        else if (distance_center2<=rA_ring(num-1)*rA_ring(num-1)) then
            m12(k1,k2)=1
        end if
        end do
    end do 
    
    !Determine masks for regions 0 (cloud free) and 1 (cloud ring)
    
    where(m12==1) m0 =0
    m1=m12-m2

    !Determine area of regions 0 and 1
    
    A0=sum(m0)
    A1=sum(m1)

    !Calculate probability of each site in region 0
    
    p0=1.0/(A1*cluster_strength+A0)
    
    !Choose region 0 or 1 randomly
    
    if (random((num-1)*2+1)<=p0*A0) then 
        !randomly pick one of the elements in A0
        random_int=floor(A0*random((num-1)*2+2))+1
        call find_new_center(L,num,random_int,m0,A0,center)
        
    else 
        !randomly pick one of the elements in A1
        random_int=floor(A1*random((num-1)*2+2))+1
        call find_new_center(L,num,random_int,m1,A1,center)
        
    end if
    
  end do

  cloud_field=m2

  cloud_centers=center(1:nclouds-1,:)
  gauss_field_out = gauss_field
! call writeMatrix(L,m2)
end subroutine

subroutine writeMatrix(L,field)
  integer, intent(in) :: L
  integer,dimension(L,L) :: field(:,:)
  integer :: i,j
  
  do i=1,L
    do j=1,L
      write(*, '(I1,X)',advance='no') field(i,j)
    end do
    write(*,*)''
  end do
  
end subroutine

real function f_r(r, r_disc, val_center, offset)
    real,intent(in) :: val_center, offset,r_disc
    real :: r
    f_r=val_center-(val_center-offset)*r/r_disc
    return
end

subroutine find_new_center(L,num,random_int,m_chosen,A_chosen,center)
    integer, intent(in) :: L,num,random_int
    integer,dimension(L,L) :: m_chosen(:,:)
    integer,dimension(L,2) :: center(:,:)
    integer :: l1, l2, cum_int 
    real :: A_chosen

    if (random_int<=nint(A_chosen/2.0)) then
    
        cum_int=0
    
        outer: do l2=1,L 
            do l1=1,L
                cum_int=m_chosen(l1,l2)+cum_int
                if (cum_int==random_int) then
                    center(num,1)=l1
                    center(num,2)=l2
                    exit outer
                end if
            end do
        end do outer   
        
    else
        cum_int=nint(A_chosen)
        
         outer2: do l2=L,1,-1 
            do l1=L,1,-1
                cum_int=cum_int-m_chosen(l1,l2)
                if (cum_int==random_int) then
                    center(num,1)=l1
                    center(num,2)=l2
                    exit outer2
                end if
            end do
        end do outer2         
        
    end if
    
end subroutine

subroutine getran(r,u,v,w)

    implicit none
    integer(kind=8), intent(out)   :: r
    integer(kind=8), intent(inout) :: u,v,w
    integer(kind=8) :: x

    u=u*2862933555777941757_8+7046029254386353087_8
    v=ieor(v,ishft(v,-17))
    v=ieor(v,ishft(v,31))
    v=ieor(v,ishft(v,-8))
    w=int(4294957665_8*iand(w,Z'ffffffff') + ishft(w,-32),8)
    x=ieor(u,ishft(u,21))
    x=ieor(x,ishft(x,-35))
    x=ieor(x,ishft(x,4))
    r=ieor(x+v,w)
    
end subroutine getran

function i2f(r) result(f)
    implicit none
    integer(kind=8), intent(in)   :: r
    real(kind=8) :: f
    f=5.42101086242752217e-20_8*r
    if (f<0) then
      f=1.0_8+f
    end if
end function i2f

subroutine initran(seed,u,v,w)
    implicit none
    integer(kind=8), intent(in)    :: seed
    integer(kind=8), intent(out)   :: u,v,w
    integer(kind=8) :: d
    
    v=4101842887655102017_8
    w=1_8
    u=ieor(seed,v)
    call getran(d,u,v,w)
    v=u
    call getran(d,u,v,w)
    w=v
    call getran(d,u,v,w)
    
end subroutine initran


end module

! program main
!   use m_percolations_float
!   integer,parameter :: L=40
!   integer,parameter :: nclouds=4
!   integer,parameter :: cluster_strength=1
!   real :: rA(nclouds), rA_ring(nclouds)
!   real :: start_time, stop_time
!   integer,dimension(L,L) :: cloud_field
!   integer,dimension(nclouds,2) :: cloud_centers
!   integer :: i
!   call cpu_time(start_time)
!   
!   do i=1,nclouds
!     rA(i)=i+2
!     rA_ring(i)=2*rA(i)
!   end do
!   
! 
! !   rA(:)=8
! !   rA_ring(:)=20
!   
!   call perc(cluster_strength,L,nclouds,rA,rA_ring,cloud_field,cloud_centers)
! !   call writeMatrix_real(L,cs_field)
!   call writeMatrix(L,cloud_field)
!   write (*,*) sum(cloud_field)
!   call cpu_time(stop_time)
!   print *, "Time:", &
!    stop_time - start_time, "seconds"  
! 
! end
