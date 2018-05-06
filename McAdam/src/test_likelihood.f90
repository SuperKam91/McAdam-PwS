!=====================================================================
!
!                            test_likelihood
!
!   Call McAdam likelihood in a randomised fashion or for a given 
!   set of cubes for testing purposes (YCP, 14/9/15)
!
!=======================================================================

      program test_likelihood 	
       
       use params
       use like
       use ReadWrite
       use Nested

 	implicit none

	integer index,i,j,k
        character*100 input_file     ! input file

        integer :: myrank,root

        double precision, allocatable, dimension(:)       :: Cube, Cube2, Cube3           ! Stacked cubes of all parameters
        double precision                            :: lhood1, lhood2, lhood3, lhood4 !likelihood values to output
        integer, dimension(1) :: seed
        integer               :: io_stat1, size, idum
        logical determ


!-----------------------------------------------------------------------

#ifdef MPI
    !MPI initializations
    call MPI_INIT(errcode)
    if (errcode/=MPI_SUCCESS) then
    	write(*,*)'Error starting MPI program. Terminating.'
    	call MPI_ABORT(MPI_COMM_WORLD,errcode)
    end if

    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, errcode)

    ! Find out who is the root
    ! ~~~~~~~~~~~~~~~~~~~~~~~~
    !
    ! 1) Get the rank of each process
    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, errcode ) 
    !
    ! 2) Find the lowest rank using MPI_MIN reduction
    call MPI_ALLREDUCE( myrank, root, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, errcode )
    !
    ! 3) Lowest rank is the root 
    is_root = (myrank==root)

    if (is_root .and. size>1) write(*,*) 'Warning: running on multiple nodes not recommended, your output will be jumbled'

#else
    ! If we're running without MPI, then the single node is the root
    is_root = .true.
#endif

!-----------------------------------------------------------------------

        ! If an include file is provided as an argument, read in parameters (otherwise assume everything's already included when compiling
        if(iargc()==1) then
          ! Put the command line provided file name into input_file
          call getarg(1,input_file) 
          if (is_root) write(*,*) 'Using ', trim(input_file)
          ! Call this subroutine to read the input file and set up the settings and priors
          call read_params(trim(input_file))
        else
          if (is_root) write(*,*) 'Warning: not using an input file'
        endif
	
	call ReadInData
      
	if (is_root) write(*,*)
	
        index = 1
	call Initialise(index)
 
	if( GL == 1 .and. varyzs == 1 ) then
		allocate(zsrcID(Ngals))
		zsrcID(1:Ngals) = int(z_s(1:Ngals))
		if( maxval(zsrcID(1:Ngals)) /= nzsplanes .or. minval(zsrcID(1:Ngals)) < 1 ) then
			call halt_program('source redshift IDs in the galaxy catalogues have incorrect entries')
		endif
	endif
	
        determ=.true.
        if (is_root) write(*,*) 'Running with determ = ', determ

        if (n_rseed .ne. -1) then
          if (is_root) write(*,*) 'n_rseed = ', n_rseed
          seed(1) = n_rseed
          call random_seed(put=seed)
        endif
	
	if( GL == 1 ) then 
		call halt_program("Haven't coded for GL likelihood")
	else
		k = 0
		do i = 1, tot_dim
      			call PClass(i,j)
      			if( j /= 0 ) k=k+1
      		enddo
                if (GasModel==5 .or. GasModel==6) then !adapted for GM=6 kj 01/02/17
                   k = k + 3
                else
                   k = k + 7
                endif

                ! Make a cube
                allocate(Cube(n_totPar))
                allocate(Cube2(n_totPar))
                allocate(Cube3(n_totPar))

                if (.not.determ) then

                   open(10, form='unformatted', file='cubes.dat', status='replace')
                   do idum = 1, 100

                      do i = 1, edim
                         call random_number(Cube(i))
                      enddo
                      write(*,*) 'First Cube', Cube
                      write(10) Cube

                      ! Save unit hypercube parameters in Cube2
                      Cube2 = Cube
                      ! Call the likelihood
                      i=0
                      call FUserbuild(1.d0,i,lhood1,i,i,i,Cube,i)
                      write(*,*) 'First Cube (physical)', Cube
                      write(*,*) 'lhood1 = ', lhood1

                      ! Restore first hypercube
                      Cube(1:edim) = Cube2(1:edim)
                      ! Change the fast parameters only
                      do i = eslow+1, edim
                         call random_number(Cube(i))
                      enddo
                      write(*,*) 'Second Cube', Cube
                      write(10) Cube

                      ! Save new hypercube in Cube3
                      Cube3 = Cube
                      i=0
                      call FUserbuild(1.d0,i,lhood2,i,i,i,Cube,i)
                      write(*,*) 'Second Cube (physical)', Cube
                      write(*,*) 'lhood2 = ', lhood2

                      ! Change all parameters
                      do i = 1, edim
                         call random_number(Cube(i))
                      enddo
                      write(*,*) 'Third Cube', Cube
                      write(10) Cube

                      i=0
                      call FUserbuild(1.d0,i,lhood3,i,i,i,Cube,i)
                      write(*,*) 'Third Cube (physical)', Cube
                      write(*,*) 'lhood3 = ', lhood3

                      ! Go back to 2nd point (change to all parameters)
                      Cube = Cube3
                      write(*,*) 'Fourth Cube', Cube
                      write(10) Cube
                      i=0
                      call FUserbuild(1.d0,i,lhood4,i,i,i,Cube,i)
                      write(*,*) 'Fourth Cube (physical)', Cube
                      write(*,*) 'lhood4 = ', lhood4

                      write(*,*) 'lhood2 - lhood4 = ', lhood2-lhood4
                      write(*,*) '======================='
                   enddo

                   close(10)

                else
                
                   open(10, file='cubes.dat', form='unformatted', action='read', iostat=io_stat1)
                   do while (io_stat1==0)

                      read(10, iostat=io_stat1) Cube
                      if (io_stat1/=0) exit
                      write(*,*) 'First Cube', Cube

                      ! Check the cube is filled correctly
                      do i = 1, edim
                        if (Cube(i)<0d0 .or. Cube(i)>1d0) then
                          call halt_program('These input cubes do not seem appropriate for this include file')
                        endif
                      enddo

                      ! Call the likelihood
                      i=0
                      call FUserbuild(1.d0,i,lhood1,i,i,i,Cube,i)
                      write(*,*) 'First Cube (physical)', Cube
                      write(*,*) 'lhood1 = ', lhood1

                      read(10, iostat=io_stat1) Cube
                      write(*,*) 'Second Cube', Cube

                      i=0
                      call FUserbuild(1.d0,i,lhood2,i,i,i,Cube,i)
                      write(*,*) 'Second Cube (physical)', Cube
                      write(*,*) 'lhood2 = ', lhood2

                      read(10, iostat=io_stat1) Cube
                      write(*,*) 'Third Cube', Cube

                      i=0
                      call FUserbuild(1.d0,i,lhood3,i,i,i,Cube,i)
                      write(*,*) 'Third Cube (physical)', Cube
                      write(*,*) 'lhood3 = ', lhood3

                      read(10, iostat=io_stat1) Cube
                      write(*,*) 'Fourth Cube', Cube

                      i=0
                      call FUserbuild(1.d0,i,lhood4,i,i,i,Cube,i)
                      write(*,*) 'Fourth Cube (physical)', Cube
                      write(*,*) 'lhood4 = ', lhood4

                      write(*,*) 'lhood2 - lhood4 = ', lhood2-lhood4
                      write(*,*) '======================='

                   enddo
                   close(10)

                endif

                deallocate(Cube)
                deallocate(Cube2)
                deallocate(Cube3)
		
	endif
	 
	if (is_root) write(*,*) 
	if (is_root) write(*,*)
	if (is_root) write(*,*)
      
      !deallocate lookup table variables
      if(GL) then
      	if(.not.simulate .and. sn_lim>0.) deallocate(snlook)
            deallocate(SigCritZd,SigCritZs,lookSigCrit,lookD)
      endif
            
      do i = 1, NAtoms
		if( NFWstyle(i) == 2 .and. ( z_PriorType(i) == 8 .and. Mass_PriorType(i,2) == 8 ) ) then
      		deallocate(lookM,lookZ)
                  exit
		endif
	enddo
#ifdef MPI
	call MPI_FINALIZE(errcode)
#endif
      
	end program test_likelihood
