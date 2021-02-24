!This code is to use Metis to divide the meshes for 50 processors
!The meaning of some variables:
!Nloc: number of basis functions for each element. Here, we use P2, so nloc is 10;
!Nface: number of faces for tetrahedron 
!nel: number of tetrahedrons used in this example
!nv: number of vertices 
!input: etov3D.txt, matrix for element to vertices 
!output: epart, the processor to which each tetrahedron belongs
program test
implicit none
integer,parameter :: Nloc=10,Nface=4 
integer,parameter :: nel=3072000,nv=531441
integer EToV(nel,Nface),epart(nel),i
!******************read in mesh infomation begin*****
OPEN(file='etov3D.txt',UNIT=10)  
        DO I=1,Nel
		READ(10,*) etov(I,:)
        end do         
close(10)

call Update_Parts(nel, nv, 3, etov, 50, epart)  
! 50 is the number of processors.

    open(unit=22, file='epart', status='unknown')
    do i = 1, nel
        write(22,*) epart(i)
    end do
    close(22)
end program test
    
subroutine Update_Parts(ne, nn, ncom, etov, nparts, epart)
    implicit none
    integer, intent(in) :: nn, ne, ncom, nparts, etov(ne,4)
    integer :: eptr(0:ne), eind(0:4*ne-1), objval, epart(0:ne-1), npart(0:nn-1)
    integer, pointer :: vwgt(:)=>null(),vsize(:)=>null(),options(:)=>null()
    integer, pointer :: tpwgts(:)=>null()
    integer :: i

    
    do i = 0, ne-1
        eind(4*i:4*i+3) = etov(i+1,:) - 1
        eptr(i) = 4*i
    end do


    eptr(ne) = 4*ne

    call METIS_PartMeshDual(ne,nn,eptr,eind,vwgt,vsize,ncom,nparts,tpwgts,options,&
                            & objval,epart,npart)
end subroutine Update_Parts