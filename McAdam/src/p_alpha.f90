!=======================================================================

module p_alpha

! A module to normalize a kernel, and return a spectral index, given
! prior limits and a cumulative probability for the spectral index.
! The kernel file 'kernel.dat' is generated using the gen_p_alpha program.
! 
! v0.1 JZ 070928

!=======================================================================

!On regenerating p_alpha, need to:
!0. Call this init_p_alpha.
!1. Read in kernel.
!2. Acquire prior limits.
!3. Decide whether you need probability or probability density [ANS:PD].
!4. Normalize kernel over these limits.

![On evaluating p_alpha, need to:]
!0. Call this p_alpha.
!1. Acquire alpha.
!2. Find nearest grid points (binary chop?).
!3. Interpolate between two nearest points.
!4. That's your p_alpha.

!=======================================================================

   use aux

   implicit none

   contains

!-----------------------------------------------------------------------

      subroutine init_p_alpha(pmin,pmax,mgrid,kernel,mu,sigma)

! ...... Reads kernel from file and normalizes according to prior limits

         use aux

         implicit none

         logical, parameter :: verbose=.false. ! TURN ON FOR DEBUGGING

         integer status
         integer nrows,ncols,nhead
         double precision, allocatable :: catalogue(:,:)
         double precision, allocatable :: head(:)
         double precision pmin,pmax
         integer ngrid,igrid,imin,imax,mgrid,i
         double precision dgrid,grid_min,grid_max,sigma1
         double precision kernel(:,:)
         double precision norm
         double precision, parameter :: epsilon = 1.0d-12
         double precision mu, sigma

! ...... Read kernel from file
         !filename='/home/cosmos-med/tws29/McAdam_v3/src/kernel.dat'
         !filename='/home/cosmos-med/tws29/McAdam_v3/src/spectral_kernel.dat'
         ncols=2
         !call read_kernel(filename,ncols,nrows,catalogue,nhead,head,status)
         if( mu == 0d0 ) then
         	call setWaldramKernal(ncols,nrows,catalogue,nhead,head,status)
	   else
         	call setWaldramGaussianKernel(ncols,nrows,catalogue,nhead,head,status,mu,sigma)
         endif
! ...... Assign kernel parameters and report
         ngrid=int(head(1))
         dgrid=head(2)
         grid_min=head(3)
         grid_max=head(4)
         sigma1=head(5)
         
         if (verbose) then
            !write(*,*) 'From ', filename,':'
            write(*,*) 'ngrid=', ngrid
            write(*,*) 'dgrid=', dgrid
            write(*,*) 'grid_min=', grid_min
            write(*,*) 'grid_max=', grid_max
            write(*,*) 'sigma=', sigma1
            write(*,*)
         endif

! ...... Check that prior is within limits of kernel; stop if not
         if(pmin<grid_min) then
            write(*,*) 'Prior range is below limit of kernel'
            write(*,*) 'pmin=',pmin
            write(*,*) 'grid_min=', grid_min
            write(*,*) 'Aborting...'
            stop
         elseif(pmax>grid_max) then
            write(*,*) 'Prior range is above limit of kernel'
            write(*,*) 'pmax=', pmax
            write(*,*) 'grid_max=', grid_max
            write(*,*) 'Aborting...'
            stop
         endif

! ...... Identify indices of prior limits in kernel by binary chop

         call locate(catalogue(:,1),ngrid,pmin,imin)
         call locate(catalogue(:,1),ngrid,pmax,imax)

         if (imin==0.or.imin==ngrid) then
            write(*,*) 'Prior is beyond limit of kernel'
            write(*,*) 'Aborting...'
            stop
         endif
         if (imax==0.or.imax==ngrid) then
            write(*,*) 'Prior is beyond limit of kernel'
            write(*,*) 'Aborting...'
            stop
         endif

         mgrid=imax-imin+1
         if(mgrid>100) then
         	write(*,*)"problem in init_p_alpha, mgrid greater than 100"
            stop
         endif
         if (verbose) then
            write(*,*) 'pmin is at imin=', imin
            write(*,*) 'pmax is at imax=', imax
            write(*,*) 'mgrid=', mgrid
         endif

! ...... Allocate kernel memory and transfer kernel
         !allocate(kernel(mgrid,3))
         kernel=0.0d0

         norm=0.0d0
         do igrid=1,mgrid
            kernel(igrid,1)=catalogue(imin-1+igrid,1)
            kernel(igrid,2)=catalogue(imin-1+igrid,2)
            !norm=norm+kernel(igrid,2)
            ! Assume bin 'values' are the left hand edge - need to start the cumulative probability at 0
            if (igrid.gt.1) then
              norm=norm+kernel(igrid-1,2)
            endif
            kernel(igrid,3)=norm
            !write(*,*) (kernel(igrid,i),i=1,3)
         enddo

! ...... Divide probs_i alone by dgrid element
         kernel(:,2) = kernel(:,2)/dgrid
         kernel(:,3) = kernel(:,3)/dgrid
         norm = norm/dgrid

         if (verbose) write(*,*) 'norm=', norm
         !if (norm.ne.sum(kernel(:,2))) then
! ALERT: ! CAN'T GET RID OF THIS, BUT IT'S ONLY 2.10^-16!
         !   write(*,*) 'norm=', norm
         !   write(*,*) 'Sum of probs =', sum(kernel(:,2))
         !   write(*,*) 'Normalization discrepancy'
         !   write(*,*) 'Aborting...'
         !   stop
         !endif

! ...... Normalize probs_i:
         kernel(:,2)=kernel(:,2)/norm
         kernel(:,3)=kernel(:,3)/norm
         norm=norm/norm
         if (verbose) then
            if (sum(kernel(:,2)) - norm > epsilon) then
         !if (norm.ne.sum(kernel(:,2))) then
! ALERT: I CAN'T GET RID OF THIS, BUT IT'S ONLY 2.10^-16!
               write(*,*) 'Normalization discrepancy:'
               write(*,*) (sum(kernel(:,2)) - norm)
               write(*,*) 'Sum of probs =', sum(kernel(:,2))
               write(*,*) 'norm=', norm
               write(*,*) 'epsilon=', epsilon
            else
               write(*,*) 'kernel is now normalized over prior range OK'
            endif
         endif

         if (verbose) then
            do igrid=1,mgrid
               write(*,*) (kernel(igrid,i),i=1,3)
            enddo
         endif
         deallocate(catalogue)

      end subroutine init_p_alpha

!-----------------------------------------------------------------------

      subroutine query_p_alpha(pcumalpha,nkernel,kernel,alpha,palpha)

! ...... Given p_cumul(alpha), returns alpha and p(alpha), from kernel

         use aux

         implicit none

         integer nkernel
         double precision kernel(:,:)
         double precision alpha,palpha,pcumalpha

! ...... Check input probability is 0<cumul(alpha)<1

         if (pcumalpha<0.0d0.or.pcumalpha>1.0d0) then
            write(*,*) 'pcum(alpha) out of range: pcum(alpha)=',pcumalpha
            write(*,*) 'Aborting...'
            stop
         endif

! ...... Linearly interpolate p(alpha) in kernel, given alpha

         !call locate(kernel(:,1),nkernel,alpha,ia)

         !write(*,*) 'alpha is at ialpha=', ia
         !write(*,*) 'alpha=', alpha
         !write(*,*) 'lower=', kernel(ia,1)
         !write(*,*) 'upper=', kernel(ia+1,1)

         call lint(kernel(:,3),kernel(:,1),kernel(:,2),nkernel,pcumalpha,alpha,palpha)

! ...... Check that p(alpha) and pcum(alpha) are 0=<P(a)=<1.

         if (palpha<0.0d0.or.palpha>1.0d0) then
            !write(*,*) 'p(alpha) out of range: p(alpha)=',palpha
            !write(*,*) 'pcumalpha, alpha=', pcumalpha, alpha
            alpha=-99.d0
         endif

      end subroutine query_p_alpha

!-----------------------------------------------------------------------

      subroutine makeWaldramGaussianKernel(mu,sigma)

! ...... Calculate the probability distribution of Waldram 07 spectral index prior multiplied by
! ...... by a Gaussian of given mean & sigma

	implicit none
	INTEGER  :: i, status
	INTEGER  :: LEN
	double precision            :: sigma,mu
	double precision , DIMENSION(:) , ALLOCATABLE ::alpha,P_alpha_Waldram,P_alpha_Gauss, P_alpha_tot
	character*60 line
      
      
      OPEN(unit=5,file='/home/cosmos-med/tws29/McAdam_v3/src/kernel.dat',status='old',iostat=status)
	OPEN(unit=10,file='/home/cosmos-med/tws29/McAdam_v3/src/spectral_kernel.dat',status='unknown')
      
      if (status.ne.0) then
		write(*,*) 'read error'
         	status = 0 ! equals: call io_clrerr(status)
         	return
      else
      	i=0
		do
            	read(5,'(A)',iostat=status) line
                  if(status.ne.0) exit
            	if (line(1:1).ne.'#') then
! ... Skip header comments denoted by '#'
               		if (line(1:3).eq.'key') then
                  		write(10,'(A)') line
               		elseif (line(1:5).eq.'ngrid') then
                  		read(line(10:13),'(3I10)') LEN
                              ALLOCATE(alpha(LEN), P_alpha_Gauss(LEN), P_alpha_Waldram(LEN), P_alpha_tot(LEN))
                  		write(10,'(A)') line
               		elseif (line(1:8).eq.'filename') then
                  		write(10,'(A)') line
               		elseif (line(1:5).eq.'dgrid') then
                  		write(10,'(A)') line
               		elseif (line(1:8).eq.'grid_min') then
                  		write(10,'(A)') line
               		elseif (line(1:8).eq.'grid_max') then
                  		write(10,'(A)') line
               		elseif (line(1:5).eq.'sigma') then
                  		write(10,'(A)') line
               		else
                  		i=i+1
                  		read(line,'(F24.16,F24.16)') alpha(i),P_alpha_Waldram(i)
                              P_alpha_Gauss(i)=Gaussfun(alpha(i),sigma,mu)
					P_alpha_tot(i)=P_alpha_Waldram(i) *P_alpha_Gauss(i)
					WRITE(10,'(F24.16,F24.16)')alpha(i),P_alpha_tot(i)
               		endif
            	endif
         	enddo
         	close(5)
         	close(10)
         	DEALLOCATE(alpha, P_alpha_Gauss, P_alpha_Waldram, P_alpha_tot)
         	status = 0
      endif

      end subroutine makeWaldramGaussianKernel

!-----------------------------------------------------------------------

      subroutine copySpectralPrior

! ...... make a copy of kernal.dat

	implicit none
	INTEGER  :: iostatus
	character*60 line
      

      OPEN(unit=5,file='/home/cosmos-med/tws29/McAdam_v3/src/kernel.dat',status='old',iostat=iostatus)
	OPEN(unit=10,file='/home/cosmos-med/tws29/McAdam_v3/src/spectral_kernel.dat',status='unknown')
      if (iostatus.ne.0) then
         	iostatus = 0 ! equals: call io_clrerr(status)
         	return
      else
		do
            	read(5,'(A)',iostat=iostatus) line
                  if(iostatus.ne.0) exit
                  write(10,'(A)') line
         	enddo
         	close(5)
         	close(10)
         	iostatus = 0
      endif

      end subroutine copySpectralPrior

!-----------------------------------------------------------------------

end module p_alpha

!=======================================================================
