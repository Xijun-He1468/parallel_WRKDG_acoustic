!Generate the files for the element connection matrices of the meshes, and the area of the four faces,
!           volume and diameter of the inscribed sphere for each tetrahedron, and two auxiliary files. 
!The meaning of some variables:
!Nloc: number of basis functions for each element. Here, we use P2, so nloc is 10;
!Nface: number of faces for tetrahedron 
!npr: the number of processors
!nel: number of tetrahedrons used in this example
!nv: number of vertices 
!input: etov3D.txt, matrix for element to vertices 
!       votex3D.txt,the coordinates for eache vertex
!       epart, the processor to which each tetrahedron belongs
!output: ftoe1.txt, matrix for face to element 
!        ftof1.txt, matrix for face to face 
!        vtoh1.txt, matrix for vertex to h !comppute h based on ADER3D (Dumbser and K?ser, 2006)
!        vsd1.txt,!compute(s1,s2,s3,s4,volume,diam) for each element,s1~s4:the areas of four faces,
!                  volume of the tetrahedron, diameter of the inscribed sphere for each tetrahedron 
!        velo1.txt, acoustic velocity for each tetrahedron
!        n1(npr),number of elements per processor 
!        n7(npr,npr),the j-th process has (n7)_ij elements adjacent to the i-th process 

module meshinfo
integer,parameter :: Nloc=10,Nface=4,npr=50
integer,parameter :: nel=3072000,nv=531441
integer EToV(nel,Nface),FToE(nel,Nface),FToF(nel,Nface),epart(nel)
integer VToH(nel,Nface) !comppute h based on ADER3D (Dumbser and K?ser, 2006)
double precision vsd(nel,6)   !compute(s1,s2,s3,s4,volume,diam)
double precision vx(nv),vy(nv),vz(nv),cc(nel),vis
end module meshinfo
    
subroutine compute_v_s_d()
!compute(s1,s2,s3,s4,volume,diam) for each element
 use meshinfo
 implicit none
    integer :: i,nodecoor(4)
    double precision :: node1(3), node2(3), node3(3), node4(3),be(3,3)
    double precision :: a1,a2,a3,a4,len12,len13,len14,len23,len24,len34,sper,det
  
     do i=1,nel
    nodecoor=EToV(i,:); 
    node1(1)=vx(nodecoor(1));node1(2)=vy(nodecoor(1));node1(3)=vz(nodecoor(1))
    node2(1)=vx(nodecoor(2));node2(2)=vy(nodecoor(2));node2(3)=vz(nodecoor(2))
    node3(1)=vx(nodecoor(3));node3(2)=vy(nodecoor(3));node3(3)=vz(nodecoor(3))
    node4(1)=vx(nodecoor(4));node4(2)=vy(nodecoor(4));node4(3)=vz(nodecoor(4))
    len12=sqrt((node2(1)-node1(1))**2.D0+(node2(2)-node1(2))**2.D0+(node2(3)-node1(3))**2.D0)
    len13=sqrt((node3(1)-node1(1))**2.D0+(node3(2)-node1(2))**2.D0+(node3(3)-node1(3))**2.D0)
    len14=sqrt((node4(1)-node1(1))**2.D0+(node4(2)-node1(2))**2.D0+(node4(3)-node1(3))**2.D0)
    len23=sqrt((node2(1)-node3(1))**2.D0+(node2(2)-node3(2))**2.D0+(node2(3)-node3(3))**2.D0)
    len24=sqrt((node2(1)-node4(1))**2.D0+(node2(2)-node4(2))**2.D0+(node2(3)-node4(3))**2.D0)
    len34=sqrt((node3(1)-node4(1))**2.D0+(node3(2)-node4(2))**2.D0+(node3(3)-node4(3))**2.D0)
    !len12,len13,len14,len23,len24,len34 , the lengths of six edges of this tetrahedron
!nodes 1,2,3
    sper=(len12+len13+len23)/2.D0
    a1=sqrt(sper*(sper-len12)*(sper-len13)*(sper-len23))
!nodes 1,2,4    
    sper=(len12+len14+len24)/2.D0
    a2=sqrt(sper*(sper-len12)*(sper-len14)*(sper-len24))
!nodes 1,3,4    
    sper=(len13+len14+len34)/2.D0
    a3=sqrt(sper*(sper-len13)*(sper-len14)*(sper-len34))
!nodes 2,3,4    
    sper=(len23+len24+len34)/2.D0
    a4=sqrt(sper*(sper-len23)*(sper-len24)*(sper-len34)) 
!volume
    be(1,1)=node2(1)-node1(1);be(2,1)=node3(1)-node1(1);be(3,1)=node4(1)-node1(1)
    be(1,2)=node2(2)-node1(2);be(2,2)=node3(2)-node1(2);be(3,2)=node4(2)-node1(2)
    be(1,3)=node2(3)-node1(3);be(2,3)=node3(3)-node1(3);be(3,3)=node4(3)-node1(3)
    det=be(1,1)*(be(2,2)*be(3,3)-be(2,3)*be(3,2))-be(1,2)*(be(2,1)*be(3,3)-be(2,3)*be(3,1))+&
        be(1,3)*(be(2,1)*be(3,2)-be(2,2)*be(3,1))
    
    vsd(i,1)=a1;vsd(i,2)=a2;vsd(i,3)=a3;vsd(i,4)=a4;
    vsd(i,5)=1.D0/6.D0*det;
    vsd(i,6)=6.0*vsd(i,5)/(a1+a2+a3+a4);

     enddo
     
end subroutine compute_v_s_d

subroutine neighbor()
!compute the element connection matrices of the meshes
 use meshinfo
 implicit none
    integer :: i, j, k,iv,iv1, iv2,iv3,iv4, mut,vnel(nv) 
    integer, allocatable :: vtoe(:,:)!e element ; v node; f face

    vnel = 0
    do i=1,nel
        do j = 1, 4
            vnel(etov(i,j)) = vnel(etov(i,j)) + 1
        end do
    end do

    allocate(vtoe(nv,maxval(vnel)))
    vtoe = -1
    vnel = 0

    do i = 1, nel
        do j = 1, 4
            iv = etov(i,j)
            vnel(iv) = vnel(iv) + 1
            vtoe(iv,vnel(iv)) = i
        end do
    end do

    do i = 1, nel
        !do j = 1, 4  !j is node
            iv1=etov(i,1);iv2=etov(i,3);iv3=etov(i,2);  !face 1
            mut = 0
            do k = 1, vnel(iv1)
                if(vtoe(iv1,k) == i) cycle
                if(any(vtoe(iv1,k)==vtoe(iv2,:)))then   !then see the third vertex
                    if(any(vtoe(iv1,k)==vtoe(iv3,:))) then
                        mut = vtoe(iv1,k);exit
                    endif
                endif
            end do
            if(mut == 0)then
                ftoe(i,1) = i
                ftof(i,1) = 1
                vtoh(i,1) = 1
            else
                ftoe(i,1) = mut
                !find ftof
                call find1(etov(mut,:),iv1,iv3,iv2,ftof(i,1),vtoh(i,1))
            endif
            
            iv1=etov(i,1);iv2=etov(i,2);iv3=etov(i,4);  !face 2
            mut = 0
            do k = 1, vnel(iv1)
                if(vtoe(iv1,k) == i) cycle
                if(any(vtoe(iv1,k)==vtoe(iv2,:)))then   
                    if(any(vtoe(iv1,k)==vtoe(iv3,:))) then
                        mut = vtoe(iv1,k);exit
                    endif
                endif
            end do
            if(mut == 0)then
                ftoe(i,2) = i
                ftof(i,2) = 2
                vtoh(i,2) = 1
            else
                ftoe(i,2) = mut
                !find ftof
                call find1(etov(mut,:),iv1,iv3,iv2,ftof(i,2),vtoh(i,2))
            endif
            
            iv1=etov(i,1);iv2=etov(i,4);iv3=etov(i,3);  !face 3
            mut = 0
            do k = 1, vnel(iv1)
                if(vtoe(iv1,k) == i) cycle
                if(any(vtoe(iv1,k)==vtoe(iv2,:)))then   
                    if(any(vtoe(iv1,k)==vtoe(iv3,:))) then
                        mut = vtoe(iv1,k);exit
                    endif
                endif
            end do
            if(mut == 0)then
                ftoe(i,3) = i
                ftof(i,3) = 3
                vtoh(i,3) = 1
            else
                ftoe(i,3) = mut
                !find ftof
                call find1(etov(mut,:),iv1,iv3,iv2,ftof(i,3),vtoh(i,3))
            endif
            
            iv1=etov(i,2);iv2=etov(i,3);iv3=etov(i,4);  !face 4
            mut = 0
            do k = 1, vnel(iv1)
                if(vtoe(iv1,k) == i) cycle
                if(any(vtoe(iv1,k)==vtoe(iv2,:)))then   
                    if(any(vtoe(iv1,k)==vtoe(iv3,:))) then
                        mut = vtoe(iv1,k);exit
                    endif
                endif
            end do
            if(mut == 0)then
                ftoe(i,4) = i
                ftof(i,4) = 4
                vtoh(i,4) = 1
            else
                ftoe(i,4) = mut
                !find ftof
                call find1(etov(mut,:),iv1,iv3,iv2,ftof(i,4),vtoh(i,4))
            endif           
    enddo
    
contains
    subroutine find1(a,a1,a2,a3,fx,fy)
    implicit none
    integer a(4),a1,a2,a3,fx,fy
    
    if(a(1)==a1) then
            if(a2==a(2)) then 
                fx=2;fy=1; 
            elseif(a2==a(3)) then
                fx=1;fy=1;
            else
                fx=3;fy=1;
            endif
    elseif(a(1)==a2) then
            if(a3==a(2)) then 
                fx=2;fy=3 
            elseif(a3==a(3)) then
                fx=1;fy=3
            else
                fx=3;fy=3;
            endif  
    elseif(a(1)==a3) then
            if(a1==a(2)) then 
                fx=2;fy=2 
            elseif(a1==a(3)) then
                fx=1;fy=2 
            else
                fx=3;fy=2 
            endif  
    else
        fx=4;
            if(a1==a(2)) then 
                fy=1 
            elseif(a1==a(3)) then
                fy=2 
            else
                fy=3 
            endif 
        
    endif
    end subroutine find1
end subroutine neighbor 

!============================================main============================================    
program pre_tetdriver
use meshinfo
implicit none
integer i,iel,i1,ip
double precision start,finish,cputime
integer n1(npr),n7(npr,npr)
!******************read in mesh infomation begin*****
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
OPEN(file='epart',UNIT=12)  
        DO I=1,nel
		READ(12,*) epart(i)
        end do
close(12)
!******************read in mesh infomation  end*****
!******************initalize*************
call neighbor()
call compute_v_s_d()
call cpu_time(start)
do i=1,nel
    cc(i)=3.D0
enddo
!=====================compute n1 and n7=================
n1=0;n7=0;
do i=1,nel
    n1(epart(i)+1)=n1(epart(i)+1)+1;
enddo
do ip=1,npr
    do i=1,nel
        if(any(epart(FToE(i,:))==ip-1) .and. epart(i)/=ip-1) then
            n7(ip,epart(i)+1)=n7(ip,epart(i)+1)+1;
        end if
    enddo
enddo
!==============output files===========================
OPEN(file='ftoe1.txt',UNIT=10)  
        DO I=1,Nel
		 WRITE(10,*) ftoe(i,:)
        end do         
close(10)
OPEN(file='ftof1.txt',UNIT=10)  
        DO I=1,Nel
		 WRITE(10,*) ftof(i,:)
        end do         
close(10)
OPEN(file='vtoh1.txt',UNIT=10)  
        DO I=1,Nel
		 WRITE(10,*) vtoh(i,:)
        end do         
close(10)
OPEN(file='vsd1.txt',UNIT=10)  
        DO I=1,Nel
         WRITE(10,"(1X,6E15.7)") vsd(i,:)
         !WRITE(10,"(1X,E15.7)") vsd(i,6)
        end do         
close(10)
OPEN(file='velo1.txt',UNIT=10)  
        DO I=1,Nel
         WRITE(10,"(1X,E15.7)") cc(i)
        end do         
close(10)
OPEN(file='n1.txt',UNIT=10)  
        DO I=1,npr
         WRITE(10,*) n1(i)
        end do         
close(10)
OPEN(file='n7.txt',UNIT=10)  
        DO I=1,npr
		 WRITE(10,*) n7(i,:)
        end do         
close(10)


call cpu_time(finish)

end program pre_tetdriver
 


  

 

  
