!=======================================================================

module aux

   !use kind_def

   implicit none
   
   double precision waldramkernel(46,2)
      
   data waldramkernel(1,1), waldramkernel(1,2) / -1.5000000000000000, 0.1254756225359338 /
   data waldramkernel(2,1), waldramkernel(2,2) / -1.3999999999999999, 1.7732218327721565 /
   data waldramkernel(3,1), waldramkernel(3,2) / -1.3000000000000000, 9.7911003963584733 /
   data waldramkernel(4,1), waldramkernel(4,2) / -1.2000000000000000, 23.0862102646245653 /
   data waldramkernel(5,1), waldramkernel(5,2) / -1.1000000000000001, 26.2969194116565461 /
   data waldramkernel(6,1), waldramkernel(6,2) / -1.0000000000000000, 16.1597604224756282 /
   data waldramkernel(7,1), waldramkernel(7,2) / -0.8999999999999999, 11.5002120716553495 /
   data waldramkernel(8,1), waldramkernel(8,2) / -0.7999999999999999, 16.1395950462539410 /
   data waldramkernel(9,1), waldramkernel(9,2) / -0.7000000000000000, 14.2294111234126301 /
   data waldramkernel(10,1), waldramkernel(10,2) / -0.6000000000000000, 21.3191166925691604 /
   data waldramkernel(11,1), waldramkernel(11,2) / -0.5000000000000000, 47.9381539355568194 /
   data waldramkernel(12,1), waldramkernel(12,2) / -0.3999999999999999, 66.9984462069081133 /
   data waldramkernel(13,1), waldramkernel(13,2) / -0.2999999999999998, 102.7599058376543582 /
   data waldramkernel(14,1), waldramkernel(14,2) / -0.2000000000000000, 142.1346559684091631 /
   data waldramkernel(15,1), waldramkernel(15,2) / -0.0999999999999999, 158.8036672582770166 /
   data waldramkernel(16,1), waldramkernel(16,2) / 0.0000000000000000, 196.7335418893877090 /
   data waldramkernel(17,1), waldramkernel(17,2) / 0.1000000000000001, 220.2885297254324257 /
   data waldramkernel(18,1), waldramkernel(18,2) / 0.2000000000000002, 230.7847691041149290 /
   data waldramkernel(19,1), waldramkernel(19,2) / 0.3000000000000000, 254.6761866181214202 /
   data waldramkernel(20,1), waldramkernel(20,2) / 0.4000000000000001, 285.5909183516143912 /
   data waldramkernel(21,1), waldramkernel(21,2) / 0.5000000000000000, 306.5397682353604978 /
   data waldramkernel(22,1), waldramkernel(22,2) / 0.6000000000000001, 263.5430514469585432 /
   data waldramkernel(23,1), waldramkernel(23,2) / 0.7000000000000002, 215.8620644876568804 /
   data waldramkernel(24,1), waldramkernel(24,2) / 0.8000000000000003, 214.1539635652230800 /
   data waldramkernel(25,1), waldramkernel(25,2) / 0.9000000000000004, 223.2812424605899082 /
   data waldramkernel(26,1), waldramkernel(26,2) / 1.0000000000000000, 197.4033683429077541 /
   data waldramkernel(27,1), waldramkernel(27,2) / 1.1000000000000001, 176.1115164212970399 /
   data waldramkernel(28,1), waldramkernel(28,2) / 1.2000000000000002, 192.7560062015037090 /
   data waldramkernel(29,1), waldramkernel(29,2) / 1.3000000000000003, 170.2673658182914096 /
   data waldramkernel(30,1), waldramkernel(30,2) / 1.4000000000000004, 112.7633077284649374 /
   data waldramkernel(31,1), waldramkernel(31,2) / 1.5000000000000000, 100.1245442371498200 /
   data waldramkernel(32,1), waldramkernel(32,2) / 1.6000000000000001, 101.4844369036505611 /
   data waldramkernel(33,1), waldramkernel(33,2) / 1.7000000000000002, 74.8303376393097039 /
   data waldramkernel(34,1), waldramkernel(34,2) / 1.8000000000000003, 45.2039365870605394 /
   data waldramkernel(35,1), waldramkernel(35,2) / 1.9000000000000004, 25.1798931285096224 /
   data waldramkernel(36,1), waldramkernel(36,2) / 2.0000000000000000, 19.4339463329444646 /
   data waldramkernel(37,1), waldramkernel(37,2) / 2.1000000000000001, 19.9591335112653304 /
   data waldramkernel(38,1), waldramkernel(38,2) / 2.2000000000000002, 19.4896523186443997 /
   data waldramkernel(39,1), waldramkernel(39,2) / 2.3000000000000003, 20.7650617282226442 /
   data waldramkernel(40,1), waldramkernel(40,2) / 2.4000000000000004, 21.7149930122883781 /
   data waldramkernel(41,1), waldramkernel(41,2) / 2.5000000000000000, 27.9962788848808906 /
   data waldramkernel(42,1), waldramkernel(42,2) / 2.6000000000000005, 30.3564711426984957 /
   data waldramkernel(43,1), waldramkernel(43,2) / 2.7000000000000002, 17.8909338951846131 /
   data waldramkernel(44,1), waldramkernel(44,2) / 2.7999999999999998, 4.6192991108049730 /
   data waldramkernel(45,1), waldramkernel(45,2) / 2.9000000000000004, 0.4815605920096815 /
   data waldramkernel(46,1), waldramkernel(46,2) / 3.0000000000000000, 0.0193687353722131 /

   contains

!-----*-----------------------------------------------------------------

      subroutine read_file(filename,ncols,nrows,catalogue,status)

! Based on subroutine read_bc03 from JZ optical work
! Reads a data file that has a header (skip header on #s)

      implicit none

      integer ifile, status
      character*64 filename
      character*82 line
      integer n,nrows,ncols,i
      double precision, allocatable :: catalogue(:,:)
      double precision, allocatable :: cat(:,:)

      integer, parameter :: MAXCAT = 200

      allocate(cat(MAXCAT,ncols))
      n = 0

      ifile = 10
      open (ifile, file=filename, status='OLD', iostat=status)
      if (status.ne.0) then
         !l = chr_lenb(filename)
         write(*,*) 'read error'
         !write(iout,*) '*** problem opening '//filename(1:l)
         status = 0 ! equals: call io_clrerr(status)
         return
      else
         do while (status.eq.0)
            read(ifile,'(A)',iostat=status) line
!            write(*,*) line
            if (line(1:1).ne.'#') then
               n=n+1
               read(line(3:5),*) cat(n,1)
               read(line(24:30),*) cat(n,2)
               read(line(36:41),*) cat(n,3)
               read(line(46:52),*) cat(n,4)
               read(line(58:63),*) cat(n,5)
               read(line(67:81),*) cat(n,6)
            endif
         enddo
         close(ifile)
         status = 0
      endif

      nrows = n - 1
      write(*,*) nrows, 'rows read from', filename

      allocate(catalogue(nrows,ncols))

      do i = 1, nrows
         catalogue(i,:) = cat(i,:)
      end do

      deallocate(cat)

   end subroutine read_file

!-----*-----------------------------------------------------------------

   subroutine read_kernel(filename,ncols,nrows,catalogue,nhead,head,status)

! Based on subroutine read_bc03 from JZ optical work
! Reads a data file that has a header (skip header on #s)

      implicit none

      integer ifile, status
      character*64 filename,origfile
      character*60 line
      integer n,nrows,ncols,i,nhead,key,temp
      double precision, allocatable :: catalogue(:,:)
      double precision, allocatable :: cat(:,:)
      double precision, allocatable :: head(:)

      integer, parameter :: MAXCAT = 200

      allocate(cat(MAXCAT,ncols))
      n = 0

      ifile = 10
      open (ifile, file=filename, status='OLD', iostat=status)
      if (status.ne.0) then
         !l = chr_lenb(filename)
         write(*,*) 'read error'
         !write(iout,*) '*** problem opening '//filename(1:l)
         status = 0 ! equals: call io_clrerr(status)
         return
      else
         do while (status.eq.0)
            read(ifile,'(A)',iostat=status) line
            if (line(1:1).ne.'#') then
! ... Skip header comments denoted by '#'
               if (line(1:3).eq.'key') then
                  read(line(10:11),'(1I10)') key
                  allocate(head(key-1))
                  head=0.0d0
                  !if(allocated(head)) write(*,*) 'head allocated OK'
               elseif (line(1:5).eq.'ngrid') then
                  read(line(10:13),'(3I10)') temp
                  head(1)=dble(temp)
               elseif (line(1:8).eq.'filename') then
                  read(line(10:),'(A)') origfile
               elseif (line(1:5).eq.'dgrid') then
                  read(line(10:20),'(F10.6)') head(2)
               elseif (line(1:8).eq.'grid_min') then
                  read(line(10:20),'(F10.6)') head(3)
               elseif (line(1:8).eq.'grid_max') then
                  read(line(10:20),'(F10.6)') head(4)
               elseif (line(1:5).eq.'sigma') then
                  read(line(10:20),'(F10.6)') head(5)
               else
                  n=n+1
                  read(line,'(F24.16,F24.16)') (cat(n,i),i=1,2)
               endif
            endif
         enddo
         close(ifile)
         status = 0
      endif

      nrows = n - 1
      write(*,*) nrows, 'rows read from', filename
      write(*,*) '  based on data from', origfile

      if(nrows.ne.head(1)) then
         write(*,*) 'Number of rows does not match header value'
         stop
      endif

      allocate(catalogue(nrows,ncols))

      do i = 1, nrows
         catalogue(i,:) = cat(i,:)
      end do

      deallocate(cat)

   end subroutine read_kernel

!-----*-----------------------------------------------------------------

   subroutine setWaldramKernal(ncols,nrows,catalogue,nhead,head,status)

      implicit none

      integer status
      integer n,nrows,ncols,i,nhead,key
      double precision, allocatable :: catalogue(:,:)
      double precision, allocatable :: head(:)

      integer, parameter :: MAXCAT = 200

      n = 0
	key = 6
      allocate(head(key-1))
	head = 0.0d0
      head(1) = 46.d0
	head(2)=0.1d0
	head(3)=-1.5d0
	head(4)=3.d0
 	head(5)=0.1d0
         
      n = 47
      status = 0
      nrows = n - 1

      if(nrows.ne.head(1)) then
         write(*,*) 'Number of rows does not match header value'
         stop
      endif

      allocate(catalogue(nrows,ncols))
      do i = 1, nrows
         catalogue(i,:) = waldramkernel(i,:)
      end do

   end subroutine setWaldramKernal

!-----------------------------------------------------------------------

      subroutine setWaldramGaussianKernel(ncols,nrows,catalogue,nhead,head,status,mu,sigma)

! ...... Calculate the probability distribution of Waldram 07 spectral index prior multiplied by
! ...... by a Gaussian of given mean & sigma

	implicit none
      integer status
      integer n,nrows,ncols,i,nhead,key
      double precision, allocatable :: catalogue(:,:)
      double precision, allocatable :: cat(:,:)
      double precision, allocatable :: head(:)
	INTEGER  :: LEN
	double precision            :: sigma,mu
	double precision , DIMENSION(:) , ALLOCATABLE ::P_alpha_Gauss, P_alpha_tot
      integer, parameter :: MAXCAT = 200

      allocate(cat(MAXCAT,ncols))
      n = 0
	key = 6
      allocate(head(key-1))
	head = 0.0d0
      head(1) = 46.d0
	head(2)=0.1d0
	head(3)=-1.5d0
	head(4)=3.d0
 	head(5)=0.1d0
      
      LEN = 46
      ALLOCATE(P_alpha_Gauss(LEN), P_alpha_tot(LEN))
      
      allocate(catalogue(nrows,ncols))
      do i = 1, LEN
		P_alpha_Gauss(i)=Gaussfun(waldramkernel(i,1),sigma,mu)
            P_alpha_tot(i)=waldramkernel(i,2) *P_alpha_Gauss(i)
            cat(i,1) = waldramkernel(i,1)
            cat(i,2) = P_alpha_tot(i)
         	catalogue(i,:) = cat(i,:)
	enddo
         	
	DEALLOCATE(P_alpha_Gauss, P_alpha_tot, cat)
	status = 0

      end subroutine setWaldramGaussianKernel

!-----------------------------------------------------------------------

      SUBROUTINE locate(xx,n,x,j)
         ! How to search an ordered table
         ! See Numerical Recipes 77 - 3.4, page 111
         ! Bisection method goes as log2(n) tries
         integer j,n
         double precision x,xx(n)
         ! Given an array xx(1:n), and given a value x, returns
         ! a value j such that x is between xx(j) and xx(j+1).
         ! xx(1:n) must be monotonic, either increasing or decreasing.
         ! j=0 or j=n is returned to indicate that x is out of range.
         integer jl,jm,ju
         jl=0                        ! Initialize lower
         ju=n+1                      ! and upper limits.
10       if(ju-jl.gt.1)then          ! If we are not yet done,
            jm=(ju+jl)/2             ! compute a midpoint,
            if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
               jl=jm                 ! and replace either the lower limit
            else                     
               ju=jm                 ! or the upper limit, as appropriate.
            endif
            goto 10                  ! Repeat until
            endif                    ! test condition 10 is satisfied.
            if(x.eq.xx(1))then       ! Then set the output
               j=1
            else if(x.eq.xx(n))then
               j=n-1
            else
               j=jl
            endif
            return                   ! and return.
      END SUBROUTINE locate

!-----*-----------------------------------------------------------------

      subroutine nearest(a,b,d)

      implicit none
      
      double precision a
      double precision b(:)
      double precision, allocatable :: db(:)
      integer i,c(1),d
      
      allocate(db(size(b)))
      
      do i = 1, size(b)
         db(i) = abs(a-b(i))
      end do
      
      c = minloc(db(:))
      
      !write(*,*) c(1), b(c(1))
      
      deallocate(db)
      
      d = c(1)

      end subroutine nearest

!-----*-----------------------------------------------------------------

      subroutine lint(xx,yy,zz,n,x,y,z)
! ... Linearly interpolates cf Numerical Recipes 77 3.3 - 3.3.1 - page 107

      implicit none

      integer n
      double precision x,xx(n),yy(n),zz(n)
      double precision y,z

      integer j
      double precision A, B

      call locate(xx,n,x,j)

      if ((j > 0) .AND. (j < n)) then
         A = (xx(j+1)-x) / (xx(j+1)-xx(j))
         B = 1 - A

         y = A*yy(j) + B*yy(j+1)
         z = A*zz(j) + B*zz(j+1)

      elseif (j == n .OR. j == 0) then
         y = -99
         z = -99
      endif

      end subroutine lint

!-----*-----------------------------------------------------------------

	FUNCTION Gaussfun(alphaa,sigmaa,muu)
	double precision, INTENT(IN) :: alphaa,sigmaa,muu
	double precision , PARAMETER :: pi=3.14159d0
      double precision :: Gaussfun
	Gaussfun = (1.0 / SQRT(2.0*pi))*EXP(((-1.0)*(alphaa-muu)*(alphaa-muu))/(2.0 *sigmaa*sigmaa))
 
	END FUNCTION Gaussfun

!-----------------------------------------------------------------------

end module aux

!=======================================================================
