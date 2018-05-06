module nestwrapper

!  Nested sampling includes
   use settings_module,          only: program_settings,initialise_settings
   use params
   use like
   implicit none
   
 contains

! Wrapper around Likelihood Function for MultiNest.
! Cube(1:n_dim) has nonphysical parameters
! scale Cube(1:n_dim) & return the scaled parameters in Cube(1:n_dim) &
! additional parameters in Cube(n_dim+1:nPar)
! return the log-likelihood in lnew
subroutine getLogLike(Cube,n_dim,nPar,lnew,context)

   integer n_dim,nPar,context
   double precision lnew,Cube(nPar)
   integer i
   
   !call your loglike function here   
   !lnew=loglike(likeindx,Cube)
   
   call FUserbuild(1.d0,i,lnew,i,i,i,Cube,i)
   if(lnew<=-1.d10) then
   	lnew=-1d90
   endif

end subroutine getLogLike

!-----*-----------------------------------------------------------------

! Wrapper around likelihood function for PolyChord.
! hyp is the unit hypercube, theta are the corresponding physical parameters, phi are derived parameters

function getLogLike_Poly(hyp,theta,phi)
!        use priors_module
        implicit none
        double precision, intent(in),  dimension(:) :: hyp           !> Unit hypercube parameters
        double precision, intent(out), dimension(:) :: theta         !> Physical parameters
        double precision, intent(out), dimension(:) :: phi           !> Output derived parameters
        double precision, dimension(n_totPar)       :: Cube           ! Stacked cube of all parameters
        double precision                            :: getLogLike_Poly ! loglikelihood value to output
        integer                                     :: i

   Cube(1:edim)=hyp(1:edim)
   call FUserbuild(1.d0,i,getLogLike_Poly,i,i,i,Cube,i)
   if(getLogLike_Poly<=-1.d10) then
   	getLogLike_Poly=-1d90
   endif

   theta(1:edim) = Cube(1:edim)
   if (GasModel==1 .or. GasModel==5 .or. GasModel==6 ) then !adapted for GM=6 kj 01/02/17
     phi(1:n_totPar-edim-3) = Cube(edim+1:n_totPar-3)
   else
     phi(1:n_totPar-edim) = Cube(edim+1:n_totPar)
   endif

 end function getLogLike_Poly

!-----*-----------------------------------------------------------------

! dumper, called after every updInt*10 iterations

subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context)

	implicit none

	integer nSamples				! number of samples in posterior array
	integer nlive					! number of live points
	integer nPar					! number of parameters saved (physical plus derived)
	double precision, pointer :: physLive(:,:)	! array containing the last set of live points
	double precision, pointer :: posterior(:,:)	! array with the posterior distribution
	double precision, pointer :: paramConstr(:)	! array with mean, sigmas, maxlike & MAP parameters
	double precision maxLogLike			! max loglikelihood value
	double precision logZ				! log evidence
	double precision logZerr			! error on log evidence
	integer context					! not required by MultiNest, any additional information user wants to pass
	
end subroutine dumper

!-----*-----------------------------------------------------------------

! Settings for PolyChord

subroutine poly_settings(settings)

    implicit none
    type(program_settings)                    :: settings
    integer                                   :: i
 
    settings%nDims         = edim                 ! Number of dimensions
    settings%nDerived      = aux_dim*NAtoms       ! Number of derived parameters
    !settings%nlive         = settings%nDims*5    ! Number of live points 
    settings%nlive         = nest_nlive           ! Number of live points
    !settings%num_repeats   = settings%nDims*5     ! Number of repeats (~5*nDims)
    if (nest_nrep == 0) then
      settings%num_repeats = eslow*5
    else
      settings%num_repeats   = nest_nrep
    endif
    settings%do_clustering = n_mmodal              ! do clustering?

    !settings%base_dir      = 'chains'             ! directory to put chains files
    i = index(n_root, '/', .true.)
    settings%base_dir      = n_root(1:i-1)
    !settings%file_root     = 'test'               ! root name for all files
    settings%file_root     = n_root(i+1:len_trim(n_root))

    settings%write_resume  = .true.              ! write a resume file?
    settings%read_resume   = .true.              ! resume from resume file?
    settings%write_live    = .true.               ! write the live points file for monitoring
    settings%write_stats   = .true.               ! write a run-time statistics file?

    settings%equals        = .true.               ! create equally weighted posteriors?
    settings%posteriors    = .true.               ! create weighted posteriors?
    settings%cluster_posteriors = n_mmodal        ! record posterior clustering?

    settings%feedback      = 1                    ! degree of feedback
    settings%update_resume = settings%nlive       ! how often to update the resume file
    settings%update_posterior = settings%nlive    ! how often to update the posterior file
                                                  !  (-1 indicates only produce it at the end)
    settings%boost_posterior= 1d0                 ! what factor to increase posteriors numbers by


    ! Multi-speed likelihoods
    !  - many likelihoods exhibit a hierarchy of parameter speeds, with some being
    !    easier to vary than others. PolyChord can exploit this.
    !  - Slow parameters should come first, followed by faster ones
    ! In our case cluster parameters are slow because of the FFT, while 
    ! radio source parameters are fast because the likelihood is analytical
    allocate(settings%grade_frac(2),settings%grade_index(2)) 
    settings%grade_frac=[0.5d0, 0.5d0]
    settings%grade_index=[1, eslow+1]

    ! Initialise the program
    call initialise_settings(settings)

end subroutine poly_settings

end module nestwrapper
