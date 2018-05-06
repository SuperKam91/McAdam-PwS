!=======================================================================

program gen_p_alpha

! A program to generate an unnormalized spectral index kernel
! from a catalogue of spectral indices.
! The output is sent to 'kernel.dat', which can be queried with
! the p_alpha module.
! 
! v0.1 JZ 070928

! 0. Call this gen_p_alpha.
! 1. Read in data from ascii.
! 2. Decide grid start, stop and increment.
! 3. Calculate kernel value at each point in grid.
! 4. Normalize probability kernel.
! 5. Write to file the array, and keywords for grid properties.

!-----------------------------------------------------------------------

   !use kind_def
   use constants

   use aux

   implicit none

   integer status
   integer ncols,nrows
   character*64 filename,outfile
   double precision, allocatable :: catalogue(:,:)
   integer irow,icol,igrid,ifile
   double precision grid_min,grid_max,dgrid
   integer ngrid
   double precision, allocatable :: grid(:,:)

!-----------------------------------------------------------------------

!..Read in the catalogue of spectral indices

   filename='alpha-1522.asc'
   ncols=6
   call read_file(filename,ncols,nrows,catalogue,status)

   write(*,*) nrows, 'sources'
!   do irow=1,nrows
!      write(*,*) (catalogue(irow,icol),icol=1,ncols)
!   enddo
   write(*,*) 'Mininum value', minval(catalogue(:,6))
   write(*,*) 'Maximum value',maxval(catalogue(:,6))

!..Define the grid for the kernel and allocate memory
   grid_min = -1.50d0
   grid_max = 3.0d0
   dgrid = 0.1d0

   ngrid = 1+int((grid_max-grid_min)/dgrid)

   write(*,*) 'Allocating', ngrid, 'grid points'

   allocate(grid(ngrid,2))
   grid=0.0d0

! Evaluate kernel at each grid point and write to file
   outfile='kernel.dat'
   ifile=10
   write(*,*) 'Writing kernel to', outfile
   open(ifile,file=outfile,status='replace')
   write(ifile,'(A)') '# Spectral index kernel for McADAM'
   write(ifile,'(A)') '# NB probs are deliberately unnormalized'
   write(ifile,'(A)') '#    because the prior range is not established until later'
   write(ifile,'(A)') '# Columns are: alpha, p(alpha)'
   write(ifile,'(A8,1X,A1)') 'key     ', '6'
   write(ifile,'(A8,1X,A)') 'filename', filename
   write(ifile,'(A8,1X,I3)') 'ngrid   ', ngrid
   write(ifile,'(A8,1X,F10.6)') 'dgrid   ', dgrid
   write(ifile,'(A8,1X,F10.6)') 'grid_min', grid_min
   write(ifile,'(A8,1X,F10.6)') 'grid_max', grid_max
   write(ifile,'(A8,1X,F10.6)') 'sigma   ', dgrid

   do igrid=1,ngrid
      grid(igrid,1) = grid_min + (igrid-1)*dgrid
      grid(igrid,2) = 0.0d0
      do irow=1,nrows
         grid(igrid,2)=grid(igrid,2)+ &
              exp((catalogue(irow,6)-grid(igrid,1))* &
              (catalogue(irow,6)-grid(igrid,1))/(-2*dgrid*dgrid))
      enddo         
      grid(igrid,2)=grid(igrid,2)/(2*pi*dgrid*dgrid)
      !write(*,*) grid(igrid,1),grid(igrid,2)
      write(ifile,'(F24.16,1X,F24.16)') grid(igrid,1),grid(igrid,2)
   enddo
   close(ifile)

!-----------------------------------------------------------------------

end program gen_p_alpha

!=======================================================================
