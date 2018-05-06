module matrix_utils

contains

!=======================================================================

! Matrix manipulation utility subroutines

! PJM 6/2002

! Convention - follows standard convention in that 
! matrices are stored by columns:

!  matrix(i,j) =   (   1,1          1,ny )
!                  (                     )
!                  (         i,j         )
!                  (                     )
!                  (  nx,1         nx,ny )
!
!
! The subscripts _ij refer to the ith column of the jth row
!
! If you pass a matrix to a subroutine which is expecting a vector,
! this is the vector that will be read:
! 
!  vector(k) = (1,1; 2,1; 3,1; ...,     i;j,     ..., nx-1,ny; nx,ny)
!                 1    2    3       (j-1)*nx + i                    etc
!=======================================================================

      subroutine MatrixCopy(nx,ny,A,B)

	integer nx,ny,i,j
	double precision A(nx,ny),B(nx,ny)

! B_ij = A_ij  

      do j=1,ny
        do i=1,nx
          B(i,j) = A(i,j)
	  enddo
	enddo
	
      return	
	end subroutine MatrixCopy
	
!=======================================================================

      subroutine Transpose(nx,ny,M,T)

	integer nx,ny,i,j
	double precision M(nx,ny),T(ny,nx)

! T_ij = M_ji  

      do j=1,ny
        do i=1,nx
          T(j,i) = M(i,j)
	  enddo
	enddo
	
      return	
	end subroutine Transpose
	
!=======================================================================

      subroutine LInvert(n,L,Linv)

	integer n,i,j
	double precision Id(n,n),L(n,n),Linv(n,n)

      do j=1,n
        do i=1,n
          Id(i,j) = 0.0
	  enddo
        Id(j,j) = 1.0
	enddo
	
	do j=1,n
	  call LFwdSub(n,L,Linv(1,j),Id(1,j))
	enddo  
	
      return	
	end subroutine LInvert
	
!=======================================================================

      subroutine LFwdSub(n,L,a,b)

	integer n,i,j
	double precision L(n,n),a(n),b(n),sum

! L_ij * a_j = b_i  ---->  a_i = L^-1_ij*b_j
	
	a(1) = b(1)/L(1,1)
	
	do i=1,n
	  
	  sum = 0.0
	  do j=1,i-1
	    sum = sum + L(i,j)*a(j)
	  enddo
	  
	  a(i) = (b(i) - sum) / L(i,i)
	  
	enddo
		
      return	
	end subroutine LFwdSub
	
!=======================================================================

      subroutine MatrixSum(iflag,nx,ny,M1,M2,M3)

	integer iflag,nx,ny,i,j
	double precision M1(nx,ny),M2(nx,ny),M3(nx,ny)

!  Addition:
      if (iflag.ge.0) then
      
        do j=1,ny
	    do i=1,nx
            M3(i,j) = M1(i,j) + M2(i,j)
	    enddo
	  enddo
	
!  Subtraction:
	elseif (iflag.lt.0) then
	
        do j=1,ny
	    do i=1,nx
            M3(i,j) = M1(i,j) - M2(i,j)
	    enddo
	  enddo
	
	endif
	
      return	
	end subroutine MatrixSum
	
!=======================================================================

      subroutine MatrixProduct(nx1,ny1,M1,nx2,ny2,M2,M3)

	integer nx1,ny1,nx2,ny2,i,j,k
	double precision M1(nx1,ny1),M2(nx2,ny2),M3(nx1,ny2),sum

	if (ny1.ne.nx2) stop 'Error in MatrixProduct: ny1 != nx2'

! M3_ij = M1_ik * M2_kj     

      do j=1,ny2
        do i=1,nx1
          
	    sum = 0.0
	    do k=1,ny1
	      if (M1(i,k).ne.0.0.and.M2(k,j).ne.0.0) then
		  sum = sum + M1(i,k)*M2(k,j)
		endif  
	    enddo
          M3(i,j) = sum

	  enddo
	enddo
	
      return	
	end subroutine MatrixProduct
	
!=======================================================================

      subroutine LMatrixProduct(n,L,M1,M2)

	integer n,i,j,k
	double precision L(n,n),M1(n,n),M2(n,n),sum

! M2_ij = L_ik * M1_kj     

      do j=1,n
        do i=1,n
          
	    sum = 0.0
	    do k=1,i
	      if (M1(k,j).ne.0.0) then
		  sum = sum + L(i,k)*M1(k,j)
		endif  
	    enddo
          M2(i,j) = sum

	  enddo
	enddo
	
      return	
	end subroutine LMatrixProduct
	
!=======================================================================

      subroutine VecMatProd(nx,ny,a,M,b)

	integer nx,ny,i,j
	double precision M(nx,ny),a(nx),b(ny),sum

! b_j = a_i * M_ij     

      do j=1,ny
	    
	  sum = 0.0
	  do i=1,nx
	    sum = sum + a(i)*M(i,j)
	  enddo
        b(j) = sum

	enddo
	
      return	
	end subroutine VecMatProd
	
!=======================================================================

      subroutine MatVecProd(nx,ny,M,a,b)

	integer nx,ny,i,j
	double precision M(nx,ny),a(nx),b(ny),sum

! b_i = M_ij * a_j     

      do i=1,nx
	    
	  sum = 0.0
	  do j=1,ny
	    sum = sum + M(i,j)*a(j)
	  enddo
        b(i) = sum

	enddo
	
      return	
	end subroutine MatVecProd
	
!=======================================================================

      subroutine LMatVecProd(n,L,a,b)

	integer n,i,j
	double precision L(n,n),a(n),b(n),sum

! b_i = L_ij * a_j     

      do i=1,n
	    
	  sum = 0.0
	  do j=1,i
	    sum = sum + L(i,j)*a(j)
	  enddo
        b(i) = sum

	enddo
	
      return	
	end subroutine LMatVecProd
	
!=======================================================================

	subroutine Matrix2Vector(nx,ny,matrix,vector)
	
      implicit none

	integer nx,ny,i,j,k
	double precision matrix(nx,ny),vector(nx*ny)

	k=0
	do j=1,ny
        do i=1,nx
	    k = k+1
	    vector(k) = matrix(i,j)
	  enddo
	enddo  

      return 
	end subroutine Matrix2Vector

!=======================================================================

	subroutine Vector2Matrix(nx,ny,vector,matrix)
	
      implicit none

	integer nx,ny,i,j,k
	double precision matrix(nx,ny),vector(nx*ny)

	k=0
	do j=1,ny
        do i=1,nx
	    k = k+1
	    matrix(i,j) = vector(k)
	  enddo
	enddo  

      return 
	end subroutine Vector2Matrix

!=======================================================================
! Given the vector index k, what are the matrix indices i and j?
	
	subroutine Vec2Mat_ID(i,j,k,nx,ny)
	
      implicit none

	integer nx,ny,i,j,k

      i = mod(k,nx)
	if (i.eq.0) i = nx
	
	j = (k - i)/nx + 1
	if (j.eq.ny+1) j = ny

      return 
	end subroutine Vec2Mat_ID

!=======================================================================
! Given the matrix indices i and j, what is the vector index k?

	subroutine Mat2Vec_ID(i,j,k,nx,ny)
	
      implicit none

	integer nx,ny,i,j,k

      k = (j-1)*nx + i

      return 
	end subroutine Mat2Vec_ID

!=======================================================================
! Given the vector index k, what are the matrix indices i and j?
	
	subroutine Vec2Mat_ID0(i,j,k,nx,ny)
	
      implicit none

	integer nx,ny,i,j,k
      
	i = mod(k+1,nx)
	if (i.eq.0) i = nx
	
	j = (k+1 - i)/nx + 1
	if (j.eq.ny+1) j = ny

	i = i-1
	j = j-1
	
      return 
	end subroutine Vec2Mat_ID0

!=======================================================================
! Given the matrix indices i and j, what is the vector index k?

	subroutine Mat2Vec_ID0(i,j,k,nx,ny)
	
      implicit none

	integer nx,ny,i,j,k

      k = j*nx + i

      return 
	end subroutine Mat2Vec_ID0

!=======================================================================

	subroutine Diag(n,v,A)
	
      implicit none

	integer i,j,n
	double precision v(n),A(n,n)

	do j=1,n
        do i=1,n
	    if (i.eq.j) then
	      A(i,j) = v(i)
	    else
	      A(i,j) = 0.0
	    endif
        enddo
      enddo

      return 
	end subroutine Diag

!=======================================================================

	subroutine QuadForm(n,v1,A,v2,q)
	
      implicit none

	integer i,n
	double precision v1(n),A(n,n),v2(n),q
	double precision v3(n)

      call MatVecProd(n,n,A,v2,v3)
	q = 0.0
	do i=1,n
	  q = q + v1(i)*v3(i)
	enddo

      return 
	end subroutine QuadForm

!=======================================================================

	function KDelta(i,j)
	
      implicit none

	integer i,j,KDelta

      if (i.eq.j) then
	  KDelta = 1
	else   
	  KDelta = 0
      endif

      return 
	end function KDelta

!=======================================================================

	subroutine Flush(n,A)
	
      implicit none

	integer i,n
	double precision A(n)

      do i=1,n
	  A(i) = 0.0
      enddo

      return 
	end subroutine Flush

!=======================================================================

	subroutine Scale(factor,n,A)
	
      implicit none

	integer i,n
	double precision factor,A(n)

      do i=1,n
	  A(i) = A(i)*factor
      enddo

      return 
	end subroutine Scale

!=======================================================================

	function VecInnerProd(n,v1,v2)
	
      implicit none

	integer i,n
	double precision v1(n),v2(n),VecInnerProd

	VecInnerProd = 0.0
	do i=1,n
	  VecInnerProd = VecInnerProd + v1(i)*v2(i)
      enddo

      return 
	end function VecInnerProd

!=======================================================================

	subroutine VecOuterProd(n,v1,v2,A)
	
      implicit none

	integer i,j,n
	double precision v1(n),v2(n),A(n,n)

      do i=1,n
        do j=1,n
	    A(i,j) = v1(i)*v2(j)
        enddo
      enddo

      return 
	end subroutine VecOuterProd

!=======================================================================

	subroutine MatZip(m,A,n,iB,B)
	
      implicit none

	integer i,m,n,iB(m)
	double precision A(m),B(m)

      n = 0
      do i=1,m
	  if (A(i).ne.0.0) then
	    n = n + 1
	    iB(n) = i
	    B(n) = A(i)
        endif
      enddo

      return 
	end subroutine MatZip

!=======================================================================

	subroutine MatUnZip(n,iB,B,m,A)
	
      implicit none

	integer i,m,n,iB(m)
	double precision A(m),B(m)

      do i=1,m
        A(i) = 0.0
	enddo
      do i=1,n
        A(iB(i)) = B(i)
	enddo

      return 
	end subroutine MatUnZip

!=======================================================================
end module matrix_utils
