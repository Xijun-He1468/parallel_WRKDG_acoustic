!This is the main program. In this test, we use 50 processors to simulate the acoustic wave propagation in Section 6.2. 
!Figure 13: v1_03 (for Tecplot 360 drawing) [run t_f=0.3D0 to obtain the snapshot at T = 0.3 s. ]
!Figure 14: wave1.txt (for waveform plot)   [run t_f=0.5D0 to obtain the waveforms at T = 0.5 s.]
!The meaning of some variables:
!   Nloc: number of basis functions for each element. Here, we use P2, so nloc is 10;
!   Nface: number of faces for tetrahedron 
!   npr: the number of processors
!   nel: number of tetrahedrons used in this example
!   nv: number of vertices 
!   nfxy, the id of the element where the source is located;nRxy,the id of the element where the receiver is located
!   etov3D, matrix for element to vertices 
!       votex3D,the coordinates for eache vertex
!       epart, the processor to which each tetrahedron belongs
!        ftoe1, matrix for face to element 
!        ftof1, matrix for face to face 
!        vtoh1, matrix for vertex to h !comppute h based on ADER3D (Dumbser and Käser, 2006)
!        vsd1,!compute(s1,s2,s3,s4,volume,diam) for each element,s1~s4:the areas of four faces,
!                  volume of the tetrahedron, diameter of the inscribed sphere for each tetrahedron 
!        velo1, acoustic velocity for each tetrahedron
!        n1(npr),number of elements per processor 
!        n7(npr,npr),the j-th process has (n7)_ij elements adjacent to the i-th process

 
module meshinfo
include 'mpif.h'
integer,parameter :: Nloc=10,Nface=4,npr=50
integer,parameter :: nel=3072000,nv=531441,nfxy=85143,nRxy=2334354
double precision, PARAMETER :: Pi=3.1415926535897932D0
integer EToV(nel,Nface),epart(nel)
integer VToH(nel,Nface),FToE(nel,Nface),FToF(nel,Nface) !comppute h based on ADER3D
!nel_k:this rank's element;  nel_b,neighbor's element;   rkb, number of neighbors
integer NProc,Myrank,inttime,endtime,totaltime,i1,ele_f,ip,i,j,ele,k,ele_r
integer status(Mpi_status_size),Ierr,vibrate
integer n1(npr),n7(npr,npr),nel_r,rkb,nel_b,nel_g,n3(nel)
integer,allocatable :: n2(:),n8(:,:),n82(:,:),n9(:,:),n10(:,:),n2g(:),n2g2(:)
integer,allocatable :: blocklengths(:),displacements(:),blocklengths0(:),displacements0(:)
double precision vsd(nel,6)   !compute(s1,s2,s3,s4,volume,diam)
double precision vx(nv),vy(nv),vz(nv),cc(nel),cn(nel),vis
double precision rou,yita,dt,ff,f0,att

double precision,allocatable::type_send(:),type_recv(:)
integer,allocatable::type_send0(:),type_recv0(:)

contains
subroutine read_info()
implicit none
OPEN(file='etov3D.txt',UNIT=10)  
        DO I=1,Nel
		READ(10,*) etov(I,:)
        end do         
close(10)
OPEN(file='votex3D.txt',UNIT=12)  
        DO I=1,Nv
		READ(12,*) vx(i),vy(i),vz(i)
        end do
close(12)
OPEN(file='ftoe1.txt',UNIT=10)  
        DO I=1,Nel
		READ(10,*) ftoe(I,:)
        end do         
close(10)
OPEN(file='vtoh1.txt',UNIT=10)  
        DO I=1,Nel
		READ(10,*) vtoh(I,:)
        end do         
close(10)
OPEN(file='ftof1.txt',UNIT=10)  
        DO I=1,Nel
		READ(10,*) ftof(I,:)
        end do         
close(10)
OPEN(file='epart',UNIT=10)  
        DO I=1,Nel
		READ(10,*) epart(I)
        end do         
close(10)
OPEN(file='vsd1.txt',UNIT=10)  
        DO I=1,Nel
		READ(10,*) vsd(I,:)
        end do         
close(10)
OPEN(file='velo1.txt',UNIT=10)  
        DO I=1,Nel
		READ(10,*) cc(i)
        end do         
close(10)
OPEN(file='n1.txt',UNIT=10)  
        DO I=1,npr
		READ(10,*) n1(I)
        end do         
close(10)
OPEN(file='n7.txt',UNIT=10)  
        DO I=1,Npr
		READ(10,*) n7(i,:)
        end do         
close(10)
end subroutine read_info
end module meshinfo
    
subroutine comp_connect()
!Calculate some connection matrices between processors 
use meshinfo
implicit none
nel_r=n1(myrank+1);allocate(n2(nel_r));i1=0 !n2, in this processor, the local ID-->the global ID for the element
do i=1,nel
    if(epart(i)==myrank) then
        i1=i1+1
        n2(i1)=i
    endif
enddo

n3=0;i1=0; !n3
do i=1,nel
    if(epart(i)==epart(nfxy)) then
        i1=i1+1
        n3(i)=i1;
    endif
enddo
ele_f=n3(nfxy);  !The local ID of the source in the processor 
n3=0;i1=0; !
do i=1,nel
    if(epart(i)==epart(nRxy)) then
        i1=i1+1
        n3(i)=i1;
    endif
enddo
ele_r=n3(nRxy);  !The local ID of the receiver in the processor 

ip=myrank+1;k=0;
do i=1,nproc
    if(n7(ip,i)/=0) then
        k=k+1;
    endif
enddo
rkb=k;allocate(n8(rkb,2));n8=0;k=0; !n8(:,1): recv processor number;n8(:,2) recv processor ID
do i=1,nproc
    if(n7(ip,i)/=0) then
        k=k+1;
        n8(k,1)=n7(ip,i);
        n8(k,2)=i;
    endif
enddo
allocate(n82(rkb,2));n82=0;k=0;!n82 send processor number and ID
do i=1,nproc
    if(n7(i,ip)/=0) then
        k=k+1;
        n82(k,1)=n7(i,ip);
        n82(k,2)=i;
    endif
enddo
allocate(n9(nel_r,rkb));allocate(n10(nel_r,rkb));!n9,n10
n9=0;n10=0;
do j=1,rkb
do i=1,nel_r
        ele=n2(i); 
        if(any(epart(FToE(ele,:))==n82(j,2)-1)) then
            n9(i,j)=1;
            n10(i,j)=ele;
        endif
enddo
enddo
nel_g=nel_r+sum(n8(:,1))  !nel_g, the number of elements for each processor. 
!nel_r is the number of inner elements, sum(n8(:,1)) is the number of virtual meshes
end subroutine comp_connect

module DG_TRI_MATRICE1
! input the Mass and stiff matrices for quadrature-free DG method in the reference element
    implicit none
    integer,parameter :: nloc1=10
    double precision:: dgmass(nloc1,nloc1), dgmdx(nloc1,nloc1), dgmdy(nloc1,nloc1), dgmdz(nloc1,nloc1)
    double precision:: dginv_mass(nloc1,nloc1)
    double precision:: dgmedge_int(nloc1,nloc1,4)    
    double precision:: dgmedge(nloc1,nloc1,4,4,3)

contains
subroutine dg_tri_matrice_init()
    implicit none
    integer :: i,i1,i2,i3
    double precision:: dgm(nloc1,nloc1)
    
      open(file='Mass.txt',UNIT=12)  
        DO I=1,10
		READ(12,*) dgm(i,:)
        end do
      close(12)       
      dgmass=dgm

      open(file='inv_ma.txt',UNIT=12)  
        DO I=1,10
		READ(12,*) dgm(i,:)
        end do
      close(12)       
      dginv_mass=dgm
      
      OPEN(file='mdx.txt',UNIT=12)  
        DO I=1,10
		READ(12,*) dgm(i,:)
        end do
      close(12)
      dgmdx=dgm
 
      OPEN(file='mdy.txt',UNIT=12)  
        DO I=1,10
		READ(12,*) dgm(i,:)
        end do
      close(12)
      dgmdy=dgm
      
      OPEN(file='mdz.txt',UNIT=12)  
        DO I=1,10
		READ(12,*) dgm(i,:)
        end do
      close(12)
      dgmdz=dgm
      
      OPEN(file='MFace.txt',UNIT=12) 
      do i1=1,4
          do i2=1,4
              do i3=1,3
        DO I=1,10
		READ(12,*) dgm(i,:)
        end do
        dgmedge(:,:,i1,i2,i3)=dgm
              enddo
          enddo
      enddo

       OPEN(file='MintF.txt',UNIT=12) 
      do i1=1,4
        DO I=1,10
		READ(12,*) dgm(i,:)
        end do
        dgmedge_int(:,:,i1) =dgm
      enddo
      
    return
end subroutine dg_tri_matrice_init
end module DG_TRI_MATRICE1

!============================================main============================================    
program tetdriver
! main program
use meshinfo
use DG_TRI_MATRICE1
implicit none

integer Nt,iit
double precision, allocatable,dimension(:,:)::v1,p1,p2,p3
double precision, allocatable,dimension(:,:)::rv1,rp1,rp2,rp3
double precision, allocatable,dimension(:,:)::v11,p11,p21,p31
double precision, allocatable,dimension(:,:)::r2v1,r2p1,r2p2,r2p3
double precision, allocatable,dimension(:) :: recei,recei2
double precision cfl,t_f,dx0,t

!******************some parameters begin*************

rou=1.D0/2.D0-sqrt(3.D0)/6.D0;yita=0.46D0
call read_info()
call dg_tri_matrice_init()
do i=1,nel
    cn(i)=vsd(i,6)/cc(i)
enddo
!******************some parameters begin*************
dx0=minval(cn);cfl=0.43D0;dt=cfl*dx0;t_f=0.3D0; ! dx0=min(d/c),cfl numberm, dt=time step, t_f=end time
att=0.0D0;
!dt=6.25D0/(att+6.25D0/cfl/dx0)
Nt=floor(t_f/dt)+1;  ! loop number
F0=30.D0;!frequecy
allocate(recei(Nt));allocate(recei2(Nt))
!write(*,*) 'para'
!write(*,*) dx0,Nt,dt
!******************some parameters end *************

Call MPI_Init(Ierr)
Call Mpi_Comm_Size(Mpi_Comm_World,NProc,Ierr)  !
call Mpi_Comm_Rank(Mpi_Comm_World,Myrank,Ierr)

inttime=MPI_WTIME()
call comp_connect()

!******************define mpi type *************
allocate(type_send(rkb));allocate(type_recv(rkb));
allocate(blocklengths(nel_r));allocate(displacements(nel_r));
allocate(type_send0(rkb));allocate(type_recv0(rkb));
allocate(blocklengths0(nel_r));allocate(displacements0(nel_r));
! form type_send,message send index 
do i=1,rkb
    k=0;blocklengths=0;displacements=0;blocklengths0=0;displacements0=0;
    do i1=1,nel_r
        if(n9(i1,i)/=0) then
            k=k+1;          
            blocklengths(k)=10;   ! 10 is the nunmber of basis functions
            displacements(k)=10*(i1-1);
            blocklengths0(k)=1;
            displacements0(k)=(i1-1);
        endif  
    enddo
    call Mpi_type_indexed(n82(i,1),blocklengths,displacements,mpi_double_precision,type_send(i),ierr)
    call Mpi_type_commit(type_send(i),ierr)
    call Mpi_type_indexed(n82(i,1),blocklengths0,displacements0,mpi_integer,type_send0(i),ierr)
    call Mpi_type_commit(type_send0(i),ierr)
enddo
! form type_recv,message recv index 
k=0;k=nel_r;
do i=1,rkb
    blocklengths=0;displacements=0;blocklengths0=0;displacements0=0;
    do i1=1,n8(i,1)       
            blocklengths(i1)=10;
            displacements(i1)=10*(i1-1)+10*k; 
            blocklengths0(i1)=1;
            displacements0(i1)=(i1-1)+k; 
    enddo
    k=k+n8(i,1);
    call Mpi_type_indexed(n8(i,1),blocklengths,displacements,mpi_double_precision,type_recv(i),ierr)
    call Mpi_type_commit(type_recv(i),ierr)
    call Mpi_type_indexed(n8(i,1),blocklengths0,displacements0,mpi_integer,type_recv0(i),ierr)
    call Mpi_type_commit(type_recv0(i),ierr)
enddo
allocate(n2g(nel_g)); allocate(n2g2(nel));  
n2g(1:nel_r)=n2;

call Mpi_barrier(Mpi_Comm_World, ierr)

do i=1,rkb
        call mpi_send(n10(:,i),1,type_send0(i),n82(i,2)-1,vibrate+1,mpi_comm_world,ierr)    
        call mpi_recv(n2g,1,type_recv0(i),n8(i,2)-1,vibrate+1,mpi_comm_world,status,ierr)
enddo

n2g2=0;
do i=1,nel_g
    n2g2(n2g(i))=i;
enddo

    if(Myrank.eq.0) then
       write(*,*) myrank,iit,dt,nt,yita,cfl,dx0,att
    endif 

    ! nel_g, the number of elements for each processor. 
allocate(v1(nloc,nel_g),p1(nloc,nel_g),p2(nloc,nel_g),p3(nloc,nel_g))
allocate(rv1(nloc,nel_g),rp1(nloc,nel_g),rp2(nloc,nel_g),rp3(nloc,nel_g))
allocate(v11(nloc,nel_g),p11(nloc,nel_g),p21(nloc,nel_g),p31(nloc,nel_g))
allocate(r2v1(nloc,nel_g),r2p1(nloc,nel_g),r2p2(nloc,nel_g),r2p3(nloc,nel_g))
v1=0.0D0;p1=0.0D0;p2=0.0D0;p3=0.0D0;rv1=0.0D0;rp1=0.0D0;rp2=0.0D0;rp3=0.0D0;
r2v1=0.0D0;r2p1=0.0D0;r2p2=0.0D0;r2p3=0.0D0;
v11=0.0D0;p11=0.0D0;p21=0.0D0;p31=0.0D0;recei=0.0D0;recei2=0.0D0

call Mpi_barrier(Mpi_Comm_World, ierr)
! time iteration
DO iit=1,nt 
    t=(iit-1)*dt+rou*dt;  
    !ff=-9.6D0*F0*(0.6D0*F0*t-1.D0)*exp (-8.D0*(0.6D0*F0*t-1.D0)**2);
    !ff=(1.D0-2.D0*(pi*f0*(t-1.D0/f0))**2)*exp(-(pi*f0*(t-1.D0/f0))**2)
    ff=(t-1.D0/f0)*exp(-(pi*f0*(t-1.D0/f0))**2)  ! source function
    call eirk(v1,p1,p2,p3,rv1,rp1,rp2,rp3)
    call Mpi_barrier(Mpi_Comm_World, ierr)
    !deal with souce
    if(myrank==epart(nfxy))then
                rv1(:,ele_f)=rv1(:,ele_f)+ff*matmul(dginv_mass,(/1.D0,0.D0,0.D0,0.D0,&
            &-0.125D0,0.D0,-0.3125D0,0.D0,0.D0,-0.5625D0/))

    endif
    if(any(epart(ftoe(nfxy,:))/=epart(nfxy)))then 
    write(*,*) 'source is located at processor boundary '
        call exchange1(rv1,rp1,rp2,rp3)
    endif

    call Mpi_barrier(Mpi_Comm_World, ierr)
	v11=v1+(1.D0-2.D0*rou)*dt*rv1;   
	p11=p1+(1.D0-2.D0*rou)*dt*rp1;
	p21=p2+(1.D0-2.D0*rou)*dt*rp2;
	p31=p3+(1.D0-2.D0*rou)*dt*rp3; 

    t=(iit-1)*dt+(1.D0-rou)*dt;
    !ff=-9.6D0*F0*(0.6D0*F0*t-1.D0)*exp (-8.D0*(0.6D0*F0*t-1.D0)**2);
    !ff=(1.D0-2.D0*(pi*f0*(t-1.D0/f0))**2)*exp(-(pi*f0*(t-1.D0/f0))**2)
    ff=(t-1.D0/f0)*exp(-(pi*f0*(t-1.D0/f0))**2)


    call eirk(v11,p11,p21,p31,r2v1,r2p1,r2p2,r2p3)
    !deal with souce
    if(myrank==epart(nfxy))then
                r2v1(:,ele_f)=r2v1(:,ele_f)+ff*matmul(dginv_mass,(/1.D0,0.D0,0.D0,0.D0,&
            &-0.125D0,0.D0,-0.3125D0,0.D0,0.D0,-0.5625D0/))

    endif
    if(any(epart(ftoe(nfxy,:))/=epart(nfxy)))then 
    write(*,*) 'source is located at processor boundary'
    call exchange1(r2v1,r2p1,r2p2,r2p3)
    endif
    call Mpi_barrier(Mpi_Comm_World, ierr)  
 
	v1=v1+dt/2.*(rv1+r2v1);
	p1=p1+dt/2.*(rp1+r2p1);
	p2=p2+dt/2.*(rp2+r2p2);
    p3=p3+dt/2.*(rp3+r2p3);

    if(Myrank.eq.epart(nRxy)) then
        recei(iit)=dot_product(v1(:,ele_R),(/1.,-1.,-1.,-1.,1.,1.,1.,1.,1.,1./))
        ! record the waveforms at the receiver
    endif
    if(Myrank.eq.epart(nfxy)) then
       write(*,*) myrank,iit,nt,dt,'tf',v1(1,ele_f)
    endif    
end do !end the loop for iit

endtime=MPI_Wtime()
totaltime=endtime-inttime
call output1(v1)  
call output2(recei,Nt)
do i=1,rkb
call MPI_TYPE_FREE(type_send(i),ierr)
call MPI_TYPE_FREE(type_recv(i),ierr)
call MPI_TYPE_FREE(type_send0(i),ierr) 
call MPI_TYPE_FREE(type_recv0(i),ierr)
enddo
! running cputime 
 if(myrank==0) then
     write(*,*) totaltime
 endif
 
CALL Mpi_Finalize(Ierr)
    end program tetdriver
    
 subroutine output2(recei,Nt)
 ! output the waveforms at the receiver
  use meshinfo
 implicit none
 double precision recei(5000),recei2(5000)
 integer nt

     if(Myrank.eq.epart(nRxy)) then
          write(*,*) nt
     OPEN(18,FILE='wave1.txt',STATUS='unknown')
  !write(*,*)Nloc
     DO i=1,nt	  
		WRITE(18,*) recei(i)
     enddo
 CLOSE(18)
 endif
 end subroutine output2
    
subroutine output1(v1)
 ! output the wavefield to generate the snapshot
 use meshinfo
 implicit none
 double precision v1(Nloc,nel_g),x,y,z,f(Nloc),uvx(nel_r),uvg(nel),uvg0(nel)

 uvx=0.D0;uvg=0.D0;uvg0=0.D0;
 do i=1,nel_r
     x=0.25D0;y=0.25D0;z=0.25D0;
     f=[1.D0, -1.+2.*x+y+z,  -1.+3.*y+z, -1.+4.*z, 1.-6.*x+6.*x**2-2.*y+6.*x*y+y**2-2.*z+6.*x*z+2.*y*z+z**2,&
        1.-2.*x-6.*y+10.*x*y+5.*y**2-2.*z+2.*x*z+6.*y*z+z**2, 1.-8.*y+10.*y**2-2.*z+8.*y*z+z**2,&
        1.-2.*x-y-7.*z+12.*x*z+6.*y*z+6.*z*z, 1.-3.*y-7.*z+18.*y*z+6.*z**2, 1.-10.*z+15.*z**2];
     uvg(n2(i))=dot_product(f,v1(:,i))
 enddo
! pressure or displacement at each center of the element

! processor rank 0 to gather the pressure
 if(myrank/=0)then
     call mpi_send(uvg,nel,mpi_double_precision,0,vibrate+100,mpi_comm_world,ierr)
 else
     uvg0=uvg;
     do j=1,npr-1
         uvg=0.D0;
        call mpi_recv(uvg,nel,mpi_double_precision,j,vibrate+100,mpi_comm_world,status,ierr)
        uvg0=uvg0+uvg;
     enddo
 endif
 
 ! wavefield data "v1_03"(for Tecplot 360 drawing)
 if(myrank==0) then
 open(15,FILE='v1_03',STATUS='unknown')
 write(15,*) 'Title="tecplot_field 3D"'
 write(15,*) 'Variables="x","y","z","u"'
 write(15,*) 'Zone n=531441,e=3072000,datapacking=block,zonetype=fetetrahedron'
 write(15,*) 'Varlocation=([1-3]=Nodal,[4]=cellcentered)'
 do i=1,nv
     write(15,"(1X,E15.7)") vx(i)
 enddo
 do i=1,nv
     write(15,"(1X,E15.7)") vy(i)
 enddo
 do i=1,nv
     write(15,"(1X,E15.7)") vz(i)
 enddo
 
 do i=1,nel
     write(15,"(1X,E15.7)") uvg0(i)
 enddo
 do i=1,Nel
    write(15,*) etov(i,:)
 end do 
 CLOSE(15)
 endif
end subroutine output1
  
 subroutine eirk(v1,p1,p2,p3,rv1,rp1,rp2,rp3)
 ! weighted Runge-Kutta time discretization
 use meshinfo
 implicit none
 double precision,dimension(nloc,nel_g),intent(in) ::v1,p1,p2,p3
 double precision,dimension(nloc,nel_g),intent(out) ::rv1,rp1,rp2,rp3
 double precision,dimension(nloc,nel_g)::kv1,kp1,kp2,kp3
 double precision,dimension(nloc,nel_g)::kkv1,kkp1,kkp2,kkp3

!call exchange(v1,p1,p2,p3)
!write(*,*) nel_g
call egl_rhs(v1,p1,p2,p3,rv1,rp1,rp2,rp3)
!call Mpi_barrier(Mpi_Comm_World, ierr)
call exchange(rv1,rp1,rp2,rp3)
call egl_rhs(rv1,rp1,rp2,rp3,kv1,kp1,kp2,kp3)
call exchange(kv1,kp1,kp2,kp3)
kv1=rv1+rou*dt*kv1;
kp1=rp1+rou*dt*kp1;
kp2=rp2+rou*dt*kp2;
kp3=rp3+rou*dt*kp3;

call egl_rhs(kv1,kp1,kp2,kp3,kkv1,kkp1,kkp2,kkp3)
call exchange(kkv1,kkp1,kkp2,kkp3)
rv1=yita*rv1+(1-yita)*kv1+yita*rou*dt*kkv1
rp1=yita*rp1+(1-yita)*kp1+yita*rou*dt*kkp1
rp2=yita*rp2+(1-yita)*kp2+yita*rou*dt*kkp2
rp3=yita*rp3+(1-yita)*kp3+yita*rou*dt*kkp3

end subroutine eirk

subroutine exchange(rv1,rp1,rp2,rp3)
! message passing
use meshinfo
implicit none
double precision,dimension(nloc,nel_g),intent(out) ::rv1,rp1,rp2,rp3
call Mpi_barrier(Mpi_Comm_World, ierr)
	do i = 0, npr-1
        if ( myrank == i ) then
            do j = 1,rkb
				!dest = neighbors(j)
                !write(*,*) myrank,'begin send'
            call mpi_send(rv1,1,type_send(j),n82(j,2)-1,vibrate+100,mpi_comm_world,ierr)
            call mpi_send(rp1,1,type_send(j),n82(j,2)-1,vibrate+101,mpi_comm_world,ierr)
            call mpi_send(rp2,1,type_send(j),n82(j,2)-1,vibrate+102,mpi_comm_world,ierr)
            call mpi_send(rp3,1,type_send(j),n82(j,2)-1,vibrate+103,mpi_comm_world,ierr)
				! call MPI_SEND(vv3,1,type_out_sig(j),dest,tag,MPI_COMM_WORLD,ierr)
            end do
            !write(*,*) myrank,'end send'
        else
            do j=1,rkb
             if ( n8(j,2)-1==i  ) then
                  !write(*,*)myrank, 'begin recv'
            call mpi_recv(rv1,1,type_recv(j),i,vibrate+100,mpi_comm_world,status,ierr)
            call mpi_recv(rp1,1,type_recv(j),i,vibrate+101,mpi_comm_world,status,ierr)
            call mpi_recv(rp2,1,type_recv(j),i,vibrate+102,mpi_comm_world,status,ierr)
            call mpi_recv(rp3,1,type_recv(j),i,vibrate+103,mpi_comm_world,status,ierr) 
            endif
            end do
            !write(*,*) myrank,'end recv'
        end if
    end do
    
call Mpi_barrier(Mpi_Comm_World, ierr)
end subroutine exchange
    
subroutine exchange1(rv1,rp1,rp2,rp3)
!If the source is at the processor boundary, exchange the information of the source
use meshinfo
implicit none
integer rankk,i0,ir
double precision,dimension(nloc,nel_g),intent(out) ::rv1,rp1,rp2,rp3

call Mpi_barrier(Mpi_Comm_World, ierr)
rankk=epart(nfxy);

        if ( myrank == rankk ) then
            do j = 1,rkb
				!dest = neighbors(j)
                !write(*,*) myrank,'begin send'
            call mpi_send(rv1,1,type_send(j),n82(j,2)-1,vibrate+100,mpi_comm_world,ierr)
            call mpi_send(rp1,1,type_send(j),n82(j,2)-1,vibrate+101,mpi_comm_world,ierr)
            call mpi_send(rp2,1,type_send(j),n82(j,2)-1,vibrate+102,mpi_comm_world,ierr)
            call mpi_send(rp3,1,type_send(j),n82(j,2)-1,vibrate+103,mpi_comm_world,ierr)
				! call MPI_SEND(vv3,1,type_out_sig(j),dest,tag,MPI_COMM_WORLD,ierr)
            end do
            !write(*,*) myrank,'end send'
        else
            do j=1,rkb
             if ( n8(j,2)-1==rankk  ) then
                  !write(*,*)myrank, 'begin recv'
            call mpi_recv(rv1,1,type_recv(j),i,vibrate+100,mpi_comm_world,status,ierr)
            call mpi_recv(rp1,1,type_recv(j),i,vibrate+101,mpi_comm_world,status,ierr)
            call mpi_recv(rp2,1,type_recv(j),i,vibrate+102,mpi_comm_world,status,ierr)
            call mpi_recv(rp3,1,type_recv(j),i,vibrate+103,mpi_comm_world,status,ierr) 
            endif
            end do
            !write(*,*) myrank,'end recv'
        end if
    
call Mpi_barrier(Mpi_Comm_World, ierr)
end subroutine exchange1
  
    
subroutine egl_rhs(v1,p1,p2,p3,rv1,rp1,rp2,rp3)
 use meshinfo
 implicit none 
 double precision,dimension(nloc,nel_g),intent(in) ::v1,p1,p2,p3
 double precision,dimension(nloc,nel_g),intent(out) ::rv1,rp1,rp2,rp3
 double precision,dimension(nloc)::lv1,lp1,lp2,lp3
 integer iel
 
 do iel=1,nel_r
     call elo_rhs(iel,v1,p1,p2,p3,lv1,lp1,lp2,lp3)
     rv1(:,iel)=lv1   
     rp1(:,iel)=lp1
     rp2(:,iel)=lp2
     rp3(:,iel)=lp3
 end do
end subroutine egl_rhs

subroutine elo_rhs(iel,v1,p1,p2,p3,rv1,rp1,rp2,rp3)
 use meshinfo
 use DG_TRI_MATRICE1
 implicit none
 double precision,dimension(nloc,nel_g),intent(in) ::v1,p1,p2,p3
 double precision,dimension(nloc),intent(out)::rv1,rp1,rp2,rp3
 double precision,dimension(nloc,nloc):: mass,mdx,mdy,mdz,mass_edge1,mass_edge2,inv_mass
 integer iel,lel,iel2,lel2,iface,iface2,ih
 double precision determ,normm(4,3),no1,no2,no3,inv_be(3,3)
 lel=n2(iel); !lel is the global ID of the element; iel is the ID in this processor
 !****************body matrice in this element*************
 mass=0.0D0;mdy=0.0D0;
 mdx=0.0D0;inv_mass=0.0D0;mass_edge1=0.0D0;mass_edge2=0.0D0
 call elem_basis0(lel,determ,inv_be,normm)
	mdx = determ*inv_be(1,1)*dgmdx + determ*inv_be(2,1)*dgmdy + determ*inv_be(3,1)*dgmdz
	mdy = determ*inv_be(1,2)*dgmdx + determ*inv_be(2,2)*dgmdy + determ*inv_be(3,2)*dgmdz
    mdz = determ*inv_be(1,3)*dgmdx + determ*inv_be(2,3)*dgmdy + determ*inv_be(3,3)*dgmdz
	mass = determ*dgmass; inv_mass = dginv_mass / determ
 !****************body matrices end*************

 ! volume integral
 rv1=-1.0*cc(lel)*cc(lel)*matmul(mdx,p1(:,iel))-1.0*cc(lel)*cc(lel)*matmul(mdy,p2(:,iel)) &
     -1.0*cc(lel)*cc(lel)*matmul(mdz,p3(:,iel))-att*matmul(mass,v1(:,iel))
 rp1=-matmul(mdx,v1(:,iel))
 rp2=-matmul(mdy,v1(:,iel))
 rp3=-matmul(mdz,v1(:,iel))
 
 ! face integral
 do iface=1,4
     lel2=ftoe(lel,iface);iface2=ftof(lel,iface);ih=vtoh(lel,iface)
     vis=max(cc(lel),cc(lel2));!numerical viscosity constant in LLF flux
!write(*,*) vis
	 mass_edge1 = dgmedge_int(:,:,iface)*vsd(lel,iface)*2.0
     no1=normm(iface,1);no2=normm(iface,2);no3=normm(iface,3);

 if(lel2/=lel) then  !which means it's not a boundary face
     mass_edge2 = dgmedge(:,:,iface,iface2,ih)*vsd(lel,iface)*2.0
     iel2=n2g2(lel2);  !
if(iel2==0) then
write(*,*) 'error'
endif
        
 rv1=rv1-(-0.5*cc(lel)**2*no1*(matmul(mass_edge1,p1(:,iel))+ matmul(mass_edge2,p1(:,iel2)))&
         &-0.5*cc(lel)**2*no2*(matmul(mass_edge1,p2(:,iel))+ matmul(mass_edge2,p2(:,iel2)))&
         &-0.5*cc(lel)**2*no3*(matmul(mass_edge1,p3(:,iel))+ matmul(mass_edge2,p3(:,iel2)))&
         &-0.5*vis*(matmul(mass_edge2,v1(:,iel2))-matmul(mass_edge1,v1(:,iel)))  )
 rp1=rp1-(-0.5*no1*(matmul(mass_edge1,v1(:,iel))+ matmul(mass_edge2,v1(:,iel2)))&
         &-0.5*vis*(matmul(mass_edge2,p1(:,iel2))-matmul(mass_edge1,p1(:,iel)))  );
 rp2=rp2-(-0.5*no2*(matmul(mass_edge1,v1(:,iel))+ matmul(mass_edge2,v1(:,iel2)))&
         &-0.5*vis*(matmul(mass_edge2,p2(:,iel2))-matmul(mass_edge1,p2(:,iel)))  );
 rp3=rp3-(-0.5*no3*(matmul(mass_edge1,v1(:,iel))+ matmul(mass_edge2,v1(:,iel2)))&
         &-0.5*vis*(matmul(mass_edge2,p3(:,iel2))-matmul(mass_edge1,p3(:,iel)))  );


  !}if 
  else !which means it's a boundary face
   
 !using absorb boundary 
 mass_edge2 = 0.D0
 iel2=iel
 rv1=rv1-(-0.5*cc(lel)**2*no1*(matmul(mass_edge1,p1(:,iel))+ matmul(mass_edge2,p1(:,iel2)))&
         &-0.5*cc(lel)**2*no2*(matmul(mass_edge1,p2(:,iel))+ matmul(mass_edge2,p2(:,iel2)))&
         &-0.5*cc(lel)**2*no3*(matmul(mass_edge1,p3(:,iel))+ matmul(mass_edge2,p3(:,iel2)))&
         &-0.5*vis*(matmul(mass_edge2,v1(:,iel2))-matmul(mass_edge1,v1(:,iel)))  )
 rp1=rp1-(-0.5*no1*(matmul(mass_edge1,v1(:,iel))+ matmul(mass_edge2,v1(:,iel2)))&
         &-0.5*vis*(matmul(mass_edge2,p1(:,iel2))-matmul(mass_edge1,p1(:,iel)))  );
 rp2=rp2-(-0.5*no2*(matmul(mass_edge1,v1(:,iel))+ matmul(mass_edge2,v1(:,iel2)))&
         &-0.5*vis*(matmul(mass_edge2,p2(:,iel2))-matmul(mass_edge1,p2(:,iel)))  );
 rp3=rp3-(-0.5*no3*(matmul(mass_edge1,v1(:,iel))+ matmul(mass_edge2,v1(:,iel2)))&
         &-0.5*vis*(matmul(mass_edge2,p3(:,iel2))-matmul(mass_edge1,p3(:,iel)))  );
  !}ielse
  endif
  
 end do !end do ifaces

 !write(*,*) ff
 rv1=matmul(inv_mass,rv1)
 rp1=matmul(inv_mass,rp1)
 rp2=matmul(inv_mass,rp2)
 rp3=matmul(inv_mass,rp3)
end subroutine elo_rhs

  !*****************form the right-hand-side****************
    !----------------------------------------------------------------------------------------------------------------------
subroutine elem_basis0(e,determ,inv_be,normm)
    use meshinfo
	implicit none
	integer, intent(in) :: e
	double precision :: determ, inv_be(3,3), be(3,3),normm(4,3),nor(3)
	double precision :: node1(3), node2(3), node3(3), node4(3)
	integer :: nodecoor(4)

	nodecoor=EToV(e,:)
    ! coordinates of the four vertices
    node1(1)=vx(nodecoor(1));node1(2)=vy(nodecoor(1));node1(3)=vz(nodecoor(1));
    node2(1)=vx(nodecoor(2));node2(2)=vy(nodecoor(2));node2(3)=vz(nodecoor(2))
    node3(1)=vx(nodecoor(3));node3(2)=vy(nodecoor(3));node3(3)=vz(nodecoor(3))
    node4(1)=vx(nodecoor(4));node4(2)=vy(nodecoor(4));node4(3)=vz(nodecoor(4))

    be(:,1)=(node2-node1);be(:,2)=(node3-node1);be(:,3)=(node4-node1);
    determ=be(1,1)*be(2,2)*be(3,3)-be(1,1)*be(2,3)*be(3,2)-be(1,2)*be(2,1)*be(3,3)+&
          &be(1,2)*be(2,3)*be(3,1)+be(1,3)*be(2,1)*be(3,2)-be(1,3)*be(2,2)*be(3,1);
    inv_be(1,1)= (be(2,2)*be(3,3)- be(2,3)*be(3,2))/determ
    inv_be(1,2)=-(be(1,2)*be(3,3)- be(1,3)*be(3,2))/determ
    inv_be(1,3)= (be(1,2)*be(2,3)- be(1,3)*be(2,2))/determ
    inv_be(2,1)=-(be(2,1)*be(3,3)- be(2,3)*be(3,1))/determ
    inv_be(2,2)= (be(1,1)*be(3,3)- be(1,3)*be(3,1))/determ
    inv_be(2,3)=-(be(1,1)*be(2,3)- be(1,3)*be(2,1))/determ
    inv_be(3,1)= (be(2,1)*be(3,2)- be(2,2)*be(3,1))/determ
    inv_be(3,2)=-(be(1,1)*be(3,2)- be(1,2)*be(3,1))/determ
    inv_be(3,3)= (be(1,1)*be(2,2)- be(1,2)*be(2,1))/determ
    
    ! four unit normal vector
    nor=cros(be(:,2),be(:,1),3);normm(1,:)=nor/sqrt(dot_product(nor,nor));
    nor=cros(be(:,1),be(:,3),3);normm(2,:)=nor/sqrt(dot_product(nor,nor));
    nor=cros(be(:,3),be(:,2),3);normm(3,:)=nor/sqrt(dot_product(nor,nor));
    nor=cros(be(:,2)-be(:,1),be(:,3)-be(:,1),3);normm(4,:)=nor/sqrt(dot_product(nor,nor));
    
    contains
function cros(a1,a2,n)
!IMPLICIT NONE
!cross product
integer :: n
double precision cros(n),a1(n),a2(n)
 cros(:)=[a1(2)*a2(3) - a2(2)*a1(3), a2(1)*a1(3) - a1(1)*a2(3), a1(1)*a2(2) - a2(1)*a1(2)]
end function cros
end subroutine elem_basis0

  !*****************form the right-hand-side end****************

  

 

  
