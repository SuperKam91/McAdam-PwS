!=====================================================================
!
!                            McAdam
!
!  Bayesian parameterised analysis of astronomical data sets
!
!  Author:	Phil Marshall
!  Editors:	Nutan Rajguru, Joern Geisbuesch, Marko Velic
!  Uses:	BayeSys3 (Skilling 2003)
! 		Jelly (MacLachlan 2003)
!
!  History: 
!	8/2002	Original version, weak lensing + SZ effect
!	1/2003	Multiple atoms combined into single parameter set
!	3/2003	Multiple data files added
!	4/2003	Group release, CVS access
!	3/2004	FFTW added 
!	6/2004	More SZ beta model styles, convergence GL constraint 
!
!=======================================================================

      program McAdam 	
       
       use params
       use MakeData1
       use like
       use ReadWrite
       use p_alpha
       use massfunction
       use nestwrapper
       use Nested
       use nested_sampling_module,   only: NestedSampling
       use settings_module,          only: program_settings,initialise_settings
       use random_module,            only: initialise_random
       !use mass_stats
       use pws_wrapper

 	implicit none

	integer index,i,j,k,context
        type(program_settings)                    :: settings
        double precision, dimension(4)            :: output_info
        character*100 input_file     ! input file
        logical file_exists

        integer :: myrank,root

!-----------------------------------------------------------------------

#ifdef MPI
	!MPI initializations
	call MPI_INIT(errcode)
	if (errcode/=MPI_SUCCESS) then
     		write(*,*)'Error starting MPI program. Terminating.'
     		call MPI_ABORT(MPI_COMM_WORLD,errcode)
  	end if

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

#else
    ! If we're running without MPI, then the single node is the root
    is_root = .true.
#endif

!-----------------------------------------------------------------------
	index = 1
    if(is_root) call Welcome

        ! If an include file is provided as an argument, read in parameters (otherwise assume everything's already included when compiling
        if(iargc().ge.1) then
          ! Put the command line provided file name into input_file
          call getarg(1,input_file) 
          inquire(file=input_file, exist=file_exists)
          if (.not.file_exists) call halt_program(trim(input_file)//' does not exist')
          if(is_root) write(*,*) 'Using ', trim(input_file)
          ! Call this subroutine to read the input file and set up the settings and priors
          call read_params(trim(input_file))
        else
            if(is_root)  write(*,*) 'Warning: not using an input file'
        endif

	if(simulate) then 
	  call Initialise(index)
	  call MakeData
	  stop
	endif
	
	call ReadInData
      
    if(is_root) write(*,*)
	
	call Initialise(index)
	if(edim>15) n_pWrap=1
 
        ! Setup for PolyChord sampler
	if (which_sampler == 'P') then
          call poly_settings(settings)
          if (n_rseed == -1) then
            call initialise_random()
          else
            call initialise_random(n_rseed)
          endif
        endif

	if( GL == 1 .and. varyzs == 1 ) then
		allocate(zsrcID(Ngals))
		zsrcID(1:Ngals) = int(z_s(1:Ngals))
		if( maxval(zsrcID(1:Ngals)) /= nzsplanes .or. minval(zsrcID(1:Ngals)) < 1 ) then
            if(is_root) write(*,*)'source redshift IDs in the galaxy catalogues have incorrect entries'
#ifdef MPI
     			call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
			stop
		endif
	endif
	
	call srand(n_rseed)
	
        if (is_root) call write_paramnames
        if (is_root) call write_ranges

	if( GL == 1 ) then 
		call nestRun(n_IS,n_mmodal,n_ceff,nest_nlive,n_tol,n_efr,edim,n_totPar,2,n_maxModes,500,-1d99, &
            	n_root,n_rseed,n_pWrap,n_fb,.true.,.true.,.false.,-1d10,0,NGeoPars+2,MassLim,MassMin,MassMax, &
		getloglike,dumper,context)
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
		
                if (which_sampler == 'M') then
                    if(is_root) write(*,*) 'Using MultiNest for sampling'
		    call nestRun(n_IS,n_mmodal,n_ceff,nest_nlive,n_tol,n_efr,edim,n_totPar,2,n_maxModes,500,-1d99, &
            	    n_root,n_rseed,n_pWrap,n_fb,.true.,.true.,.false.,-1d10,0,k,MassLim,MassMin,MassMax, &
		    getloglike,dumper,context)
                elseif (which_sampler == 'P') then
                    if(is_root) write(*,*) 'Using PolyChord for sampling'
                    output_info = NestedSampling(getloglike_poly,settings,MPI_COMM_WORLD)
                else
                    if(is_root) write(*,*) 'Unknown sampling type ', which_sampler
                endif
	endif
	 
    if(is_root) write(*,*) 
    if(is_root) write(*,*)
    if(is_root) write(*,*)
      
      ! Clean up Planck stuff
      if (PL == 1 ) then
        call release_pws_patch()
      endif

      !deallocate lookup table variables
      if(GL==1) then
      	if(.not.simulate .and. sn_lim>0.) deallocate(snlook)
            deallocate(SigCritZd,SigCritZs,lookSigCrit,lookD)
      endif
      
      deallocate(n_pWrap)
      
      do i = 1, NAtoms
		if( NFWstyle(i) == 2 .and. ( z_PriorType(i) == 8 .and. Mass_PriorType(i,2) == 8 ) ) then
      		deallocate(lookM,lookZ)
                  exit
		endif
	enddo
#ifdef MPI
	call MPI_FINALIZE(errcode)
#endif
      
	end program McAdam
