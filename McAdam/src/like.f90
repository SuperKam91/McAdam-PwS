module like
	use params
	use constants
	use fftw_inc
	use CheckPars1
	use ReadWrite
	use MassModels
	use GasModels
	use priors1
	use massfunction
        use utilities
	use RandomNS
        use mkl_dss
	use pws_wrapper
	      
	implicit none

	integer nTotVis
	double precision y2jy_nu0,f_nu

contains

subroutine FUserbuild(Cool,Nsystem,Lhood,Nstore,idummy,ndummy,Cube1,retcode)

    implicit none

#ifdef MPI
    include 'mpif.h'
#endif

    integer Nsystem,Nstore,idummy,ndummy,retcode,i,k,i1,j,flag,list(NAtoms),id,err
    double precision Cool,Lhood,Cube1(:),Cube(tot_dim),Cube2(tot_dim)
    double precision SZLhood,urv,zsrc(nzsplanes), PLLhood

    integer, parameter :: PL_Npars = 6
    !double precision   :: PL_Cube(PL_Npars*NAtoms + 1) !For now I will implement this so that PwS takes a hyperparameter whether this method is being used or not (to make implementation in PwS easier). If it isn't being used, hyperparameter will always be one
    double precision   :: PL_Cube(PL_Npars*NAtoms)

    !-----------------------------------------------------------------------

    id = 0
#ifdef MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD, id, err)
#endif

    !-----------------------------------------------------------------------
    flag=0
    Cube(1:edim)=Cube1(1:edim)
    Cube(edim+1:tot_dim)=0.d0
    GLLhood=0.d0
    SZLhood=0.d0
    PLLhood=0.d0

    !	add the parameters with delta priors/being marginalised outside of MultiNEST

    do i=1,tot_dim
        call PClass(i,j)
        if(j<=0) Cube(i+1:tot_dim)=Cube(i:tot_dim-1)
        ! draw from uniform distribution between zero and one for parameters marginalised outside MultiNEST
        if (j<0) then
            Cube(i)=ranmarNS(id)
        endif
    enddo

    ! 	Additional prior constraints: check parameters and return both an
    ! 	additional contribution to the likelihood and a flag value (for
    ! 	truncated distributions)

    !       Save the old old cluster parameters in case the new cube is rejected
    if (SZ==1) then
        !       Save the old cluster parameters for comparison with the new ones to enable fast/slow likelihood separation
        GeoPars_old = GeoPars
        MassPars_old = MassPars
        GasPars_old = GasPars
        TPars_old = TPars
        z_old = z
    endif

    if (hyperparameters == 1) hyperPars_old = hyperPars

    call CheckPars(Cube,flag)
    if(flag==1) goto 999

    if( GL == 1 ) then
        if( varyzs == 1 ) then
            zsmin=100.0
            zsmax=0.0
            do i=1,nzsplanes
                urv=ranmarNS(id)
                zsrc(i)=Prior(zs_PriorType(i),urv,zs_Prior(i,1),zs_Prior(i,2),zs_Prior(i,3))
                if( zsrc(i) < 0.0 ) then
                    flag = 1
                    goto 999
                endif
                if( zsrc(i) < zsmin ) zsmin = zsrc(i)
                if( zsrc(i) > zsmax ) zsmax = zsrc(i)
            enddo

            do i=1,Ngals
                z_s(i) = zsrc(zsrcID(i))
            enddo
        endif

        do i = 1, NAtoms
            if( z(i) > zsmax ) then
                flag=1
                goto 999
            endif
        enddo
    endif

!   Calculate Planck likelihood (first to adjust sampled parameters to the gridded ones)
    if(PL==1) then
      ! Put parameters into a new cube for PwS likelihood
      ! For simplicity never assume UPP
      ! Order of Planck cube is: Xcoord, Ycoord, Y5r500, thetaS, beta, alpha
      ! Reverse Xcoord direction for consistency with Healpix (?)
      do j = 1, NAtoms
        if (SZ==0 .or. GasModel==5 .or. GasModel==6) then
	  ! Need this to fill in the auxiliary parameters
          call MakeGasDistributions(j,flag)
          if (flag.ne.0) goto 999
        endif
        k = (j-1)*PL_Npars + 1
        PL_Cube(k) = GeoPars(j,1)*-1.d0
        PL_Cube(k+1) = GeoPars(j,2)
        PL_Cube(k+2) = Ytot_Planck !relies on MakeGasDistributions being called a-priori for GasModel==5/6
        PL_Cube(k+3) = thetas_Planck !or CheckGasPars for GasModel==3
        PL_Cube(k+4) = GasPars(j,5)
        PL_Cube(k+5) = GasPars(j,4)
      enddo
      k = PL_Npars * NAtoms
      !whilst not implementing hyperparameter method, don't need to pass hyperparameter to PwS
      !if (hyperparameters == 1) then 
	!if (SZ == 1) then
	!	PL_Cube(k+1) = hyperPars(1, 2) !second hyperparameter if AMI and Planck data is being used
	!else
	!	PL_Cube(k+1) = hyperPars(1, 1) !first hyperparameter if only Planck data is being used
	!endif
      !else
	!PL_Cube(k+1) = 1.d0 !for current implementation, always set hyperparameter for Planck a value
      !endif
      call PlanckLhood(PL_Cube, flag, PLLhood, PL_Cube)
      if (flag.ne.0) goto 999
      ! Substitute the gridded alpha and beta values for the SZLhood calculation for consistency?
!     do j = 1, NAtoms
!        k = (j-1)*PL_Npars + 1
!        GasPars(j,4) = PL_Cube(k+5)
!        GasPars(j,5) = PL_Cube(k+4)
!     enddo
    endif

    ! 	Make predicted data:

    if(GL==1) call PredictShear(Cube,flag)
    if(SZ==1) call PredictVisibilities(Cube,flag)
    if(flag==1) goto 999 

    ! 	Skip L calculation if sampling prior:

    if(SamplePrior) then
      Lhood=0d0
      goto 999
    endif

    ! 	Calculate Likelihood	

    !if(GL==1 .and. ModelClass==1 .and. .not.survey) call ShearLhood(GLLhood)
    !if(GL==1 .and. ModelClass==1) call ShearLhood(GLLhood)
    if(SZ==1) call VisibilitiesLhood(SZLhood)

    Lhood=GLLhood+SZLhood+PLLhood    
	
    ! 	Return

    999	if(flag==0) then
        retcode=1
    else
        retcode=0
        Lhood=B_A_D
        if (verbose) write(*,*) 'Bad set of parameters, discarding'
        !       Restore the original old cluster parameters
        if (SZ==1) then
            GeoPars = GeoPars_old
            MassPars = MassPars_old
            GasPars = GasPars_old
            TPars = TPars_old
            z = z_old
        endif
        if (hyperparameters == 1) hyperPars = hyperPars_old
        return
    endif

    !	return the scaled (physical) parameters for nested sampling
    do i=1,tot_dim
        call Rescale_nest(Cube(i),i,i1)
    enddo

    list = 0
    list(NAtoms) = 1
    do i = 2, NAtoms !only does something if NAtoms >= 2. If not skips loop
        do j = 1, NAtoms
            if( list(j) == 0 ) cycle
            if( Cube( (i - 1) * NPars + 1 ) < Cube( (list(j) - 1) * NPars + 1 ) ) then
                k = j - 1
                exit !terminates inner-most do-loop
            endif

            if( j == NAtoms ) then
                k = NAtoms
                exit
            endif
        enddo !don't believe this needs to be changed for hyperparameter method. Determines k value

        list(1:k - 1) = list(2:k)
        list(k) = i
    enddo

    do i = 1, NAtoms !the following doesn't do anything if Natoms = 1
        do j = 1, NPars
            Cube2( (i - 1) * NPars + j ) = Cube( (list(i) - 1) * NPars + j ) !this and the above is confusing I don't know (care) what it does
        enddo
    enddo
    Cube2(NAtoms * NPars + 1: tot_dim) = Cube(NAtoms * NPars + 1: tot_dim)
    Cube = Cube2

    !	remove the parameters with delta priors/being marginalised outside of MultiNEST

    k=0
    do i=1,tot_dim
        call PClass(i,j)
        if(j>0) then
            k=k+1
            Cube1(k)=Cube(i)
        endif
    enddo

    if((SZ==1).or.(PL==1)) then
        do i=1,NAtoms
            Cube1(k+1:k+aux_dim)=aux(i,1:aux_dim)
            k=k+aux_dim
        enddo
        IF(GasModel==1 .OR. GasModel==5 .OR. GasModel==6 )THEN	!adapted for GM=6 kj 01/02/17	
            Cube1(k+1)=ymap(nx/2,ny/2)
            Cube1(k+2)=y2Jy_nu0*Cube1(k+1)
            Cube1(k+3)=f_nu*Cube1(k+1)*TCMB
        ENDIF		
    endif

    if( GL == 1 .and. varyzs == 1 ) Cube1(n_totPar-nzsplanes+1:n_totPar) = zsrc(1:nzsplanes)

end subroutine FUserbuild

!=======================================================================

subroutine PredictShear(Cube,flag)
	
	implicit none

	integer i,flag
	double precision Cube(tot_dim)
	logical calcLike
	
!-----------------------------------------------------------------------

! 	First zero the working arrays:
	
	if(.not.survey) then
		gamma1=0d0
	  	gamma2=0d0
	  	kappa=0d0
        	g1=0d0
        	g2=0d0
	else
        	gamma1s=0d0
        	gamma2s=0d0
	endif
      
	do i=1,Natoms
	
! 	get the angular diameter distance

		call lookUp1D(Dn,z(i),lookD(:,1),lookD(:,2),D)

!       	Sort out elliptical geometry:
        	if(GeoModel==2) call EllGeometry(i)
	    
!       	Now increase the shear and convergence at each galaxy position:
	  	if(i==Natoms) then
        		calcLike=.true.
        	else
        		calcLike=.false.
        	endif
	  	call LensFields(i,flag,calcLike)
	enddo
      
end subroutine PredictShear

!=======================================================================

	subroutine PredictVisibilities(Cube,flag)
	
	implicit none

	integer i,flag,ipar
      double precision Cube(*)

!-----------------------------------------------------------------------
      
	noCluster = .true.
    !noCluster = any(GeoPars/=GeoPars_old) .and. any(MassPars/=MassPars_old) .and. any(GasPars/=GasPars_old) .and. any(TPars/=TPars_old)  .and. any(z/=z_old) 

        do i=1,NAtoms
           do ipar=1, nGeoPars
             if (.not. GeoPars(i,ipar)==GeoPars_old(i,ipar)) then
                noCluster = .false.
             endif
           enddo
           do ipar=1, nMassPars
             if (.not. MassPars(i,ipar)==MassPars_old(i,ipar)) then
                noCluster = .false.
             endif
           enddo
           do ipar=1, nGasPars
             if (.not. GasPars(i,ipar)==GasPars_old(i,ipar)) then
                noCluster = .false.
             endif
           enddo
           do ipar=1, nTPars
             if (.not. TPars(i,ipar)==TPars_old(i,ipar)) then
                noCluster = .false.
             endif
           enddo
           if (.not. z(i) == z_old(i)) then
             noCluster = .false.
           endif
        enddo

        if (.not. noCluster ) then
          !	initialize working arrays
          pvisr = 0d0
          pvisi = 0d0
          ymap = 0d0
  	  do i=1,NAtoms
		!Sort out elliptical geometry:
		if( GeoModel == 2 ) call EllGeometry(i)
                  
  		call MakeGasDistributions(i,flag)
		
        	if( flag == 1 ) then
			return
		elseif( flag == 0 ) then
       			!Generate Comptonisation parameter map:
          	        call MakeYMap(i)
		elseif( flag == 2 ) then
                        noCluster = .true.
			flag = 0
		endif
	  enddo
        endif

	!Scale, Fourier transform and sample to get visibilities:
        if( .not.noCluster ) call MakeVisibilities
      
	!Add point sources if required:
	if(Nuisance==1.and.SourceSubtract==1) call AddSources
      
 999	return
      end subroutine PredictVisibilities

!=======================================================================

	subroutine ShearLhood(GLLhood)
	
	implicit none

	integer i,j,k,count
      double precision GLLhood,Chisq
      double precision el1,el2,wte1,wte2
      integer eflag
	
      Chisq=0d0
	
      if(survey .or. arrGal) then
      	do i=1,nxpix
            	do j=1,nypix
                  	if(survey) then
                        	count=1
				else
                        	count=ptInPix(i,j,1)
				endif
                        
                        do k=1,count
                        	if(arrGal) then
						el1=e1(ptInPix(i,j,k+1))
						el2=e2(ptInPix(i,j,k+1))
						wte1=wt1(ptInPix(i,j,k+1))
						wte2=wt2(ptInPix(i,j,k+1))
                                    eflag=GLeflag(ptInPix(i,j,k+1))
					elseif(survey) then
						el1=e1s(i,j)
						el2=e2s(i,j)
						wte1=wt1s(i,j)
						wte2=wt2s(i,j)
                                    eflag=GLeflags(i,j)
					endif
                                    
                  		if(eflag==0) then
                        		Chisq=Chisq+(el1)*wte1*(el1)
	    					Chisq=Chisq+(el2)*wte2*(el2)
					endif
				enddo
			enddo
		enddo
	else
		do k=1,Ngals
	  		if(GLeflag(k)==0) then
	    			!Chisq=Chisq+(e1(k)-g1(k))*wt1(k)*(e1(k)-g1(k))
	    			!Chisq=Chisq+(e2(k)-g2(k))*wt2(k)*(e2(k)-g2(k))
	    			Chisq=Chisq+(e1(k)-gamma1(k))*wt1(k)*(e1(k)-gamma1(k))
	    			Chisq=Chisq+(e2(k)-gamma2(k))*wt2(k)*(e2(k)-gamma2(k))
                        if(gamma1(k)/=0. .or. gamma2(k)/=0.) count=count+1
	  		endif	    
		enddo
	endif
      
	GLLhood=GLLhood0-Chisq/2.
		
	return
	end subroutine ShearLhood
	
!=======================================================================
! Note: when move was made to multiple covariance matrices, the memory
! usage grew considerably, leading (I think) to a bug in this part of
! the code where the b vector was being calculated but then replaced
! with INFs and NANs when compiled with g77. Using the -fno-automatic
! seems to have fixed this. 

 subroutine VisibilitiesLhood(SZLhood)

   implicit none

   integer i,j,k,m, error
   double precision SZLhood,Chisq

   ! dense matrix variables and work space      
   double precision a(2*nTotVis),b(2*nTotVis),row(2*nTotVis),sum

   !-----------------------------------------------------------------------
   if(IncludeCMB==1) then

      Chisq=0.d0
      b=0d0 
 
      ! Fill in residuals vector
      j=0
      do m=1,nVisfiles    
         do i=1,Nvis(m)
            j=j+1
            a(j)=visr(m,i)-pvisr(m,i) !real part
            a(j+nTotVis)=visi(m,i)-pvisi(m,i) !imaginary part
         enddo
      enddo

      IF(dsflag.EQ. 0)THEN  

         !               dense matrix format:
         !  		Compute chi-squared as the dot product b^T b, where
         !   		Lb=a and a is the vector of residuals:	

         ! Approximation: CMB=0 for Ryle Telescope dishes. Only apply
         ! full covariance matrix if IncludeCMB(m)=1, 
         ! otherwise just use noise values in the data file.

         k=1
         do i=1,2*nTotVis
            do j=1,i
               row(j)=LCM(k)
               k=k+1
            enddo

            if(i==1) then
               b(1)=a(1)/row(1)
            else
               sum=0.d0
               do j=1,i-1
                  sum=sum+row(j)*b(j)
               enddo
               b(i)=(a(i)-sum)/row(i)
            endif
         enddo

      ELSEIF(dsflag.EQ.1)THEN	

         ! New way using intel sparse solver (YCP 28/9/2015)
         error = dss_solve_real(handle, MKL_DSS_FORWARD_SOLVE, a, 1, b)
         IF (error /= MKL_DSS_SUCCESS) call halt_program('dss_solve_real failed')

      ENDIF

      ! Calculate chi^2 using the solution vector
      do i=1,2*nTotVis
         Chisq=Chisq+b(i)*b(i)
      enddo
      if (hyperparameters == 1 .and. SZ == 1) then
      	SZLhood=SZLhood0 + (2*nTotVis)/2 * log(alpha_L1) - alpha_L1 * Chisq/2.d0
      else      	
	SZLhood=SZLhood0-Chisq/2.d0	
      endif

   else
      Chisq=0.0	
      !	Compute chisquared for diagonal thermal noise only covariance
      !    	matrix, ie the usual independent Gaussian errors routine:

      do m=1,Nvisfiles
10       do k=1,Nvis(m)
            if(SZeflag(m,k)==0) then
               Chisq=Chisq+(visr(m,k)-pvisr(m,k))*viswt(m,k)*(visr(m,k)-pvisr(m,k))
               Chisq=Chisq+(visi(m,k)-pvisi(m,k))*viswt(m,k)*(visi(m,k)-pvisi(m,k))
            endif
         enddo
      enddo
      if (hyperparameters == 1 .and. SZ == 1) then
      	SZLhood=SZLhood0 + (2*nTotVis)/2 * log(alpha_L1) - alpha_L1 * Chisq/2.d0
      else      	
	SZLhood=SZLhood0-Chisq/2.d0	
      endif   
   endif

 end subroutine VisibilitiesLhood
	
!=======================================================================

	subroutine Initialise(index)
	
	implicit none
	
	integer i,j,k,m,index,iend,Wastage
	double precision SZLhood,GLLhood
	double precision d1,d2

!-----------------------------------
!initialization stuff for both lensing & SZ
      
        tot_dim=NDim
        tot_atoms=NAtoms
      
	if( GL==1 ) then
		Mass = 1
	else
		Mass = MMass
	endif
      
	if( SZ ==1 .and. NSrc > 0 ) then
		Nuisance = 1
		SourceSubtract = 1
	else
		Nuisance = 0
		SourceSubtract = 0
	endif
      
!	reduce the dimensionality if there are parameters with delta priors
	edim=0
        eslow=0
	do i=1,tot_dim
		call PClass(i,j)
		if(j>0) then
                  edim=edim+1
                  if (i.le.NAtoms*NPars) then !treat hyperparameters as fast parameters in PolyChord
                    eslow=eslow+1
                  endif
                endif
	end do
      
	!set the total number of parameters
	if((SZ==1).or.(PL==1)) then
		IF(GasModel==1 .OR. GasModel==5 .OR. GasModel==6)THEN	!adapted for GM=6 kj 01/02/17	
			n_totPar=edim+aux_dim*NAtoms+3
		ELSE
			n_totPar=edim+aux_dim*NAtoms
		ENDIF

! choose sparse or dense matrix routines
         if (SZ==1) then
		nTotVis=sum(Nvis)
                IF(nTotVis.LE.1000)THEN
	            dsflag=0      ! dense matrix approach
                ELSE
                    dsflag=1        ! sparse matrix approach
                ENDIF	
	endif
			
	elseif(GL==1) then
		n_totPar=edim
	endif
	if(GL==1 .and. varyzs==1) n_totpar=n_totpar+nzsplanes
	
	allocate(n_pWrap(edim))
	n_pWrap=0
      
	if(.not.simulate) then
                call CheckzPriors
                call CheckMassPriors
		call CheckGeoPriors
                call CheckSrcPriors

		if (PL==1) call CheckPlanckPriors
		!no checks needed on hyperpriors as of yet. KJ 06/06/17
		if( GL == 1 .and. varyzs == 1 ) call CheckzsPriors
             
		!if the background galaxies all have the same redshift
		if( GL==1 .and. .not.n_mmodal .and. zdmax > zsmin ) zdmax = zsmax
            
		if( zdmin == zdmax ) then
			Dn = 1
			SCZdn = 1
		else  
          		Dn = n
			SCZdn = n
		endif
            
        if(is_root) write(*,*)"making Dlookup"
		allocate(lookD(Dn,2))
		call makeDlookup(zdmin,zdmax)
        if(is_root) write(*,*)"Dlookup done"
	
	endif
!-----------------------------------

! Gravitational lensing stuff:
	  
    GLLhood0=0.d0
	if(GL==1 .and. .not.simulate) then
			
	    !allocate memory for sigCrit lookup table
	    !varying source redshift
	    if(vary_zs==1 .or. varyzs==1) then
            	SCZsn=n
	    else
            	SCZsn=1
	    endif
            allocate(SigCritZd(SCZdn),SigCritZs(SCZsn),lookSigCrit(SCZdn,SCZsn))
            
            !make critical density lookup table
            if(is_root) write(*,*)"making Siglookup"
            call makeSigCritlookup(zdmin,zdmax)
            if(is_root) write(*,*)"Siglookup done"
                  
      	    !allocate memory for the S/N lookup table
      	    if(MassModel==1 .and. NFWstyle(1)==2 .and. sn_lim>0. .and. .not.simulate) then
            	if(Mass_PriorType(1,1)==0) then
               		snn=1
		else
                  	snn=n
		endif
            	allocate(snlook(snn,2))
            	!make the sn lookup table
                if(is_root) write(*,*)"making SNlookup"
            	call makeSNlookup
                if(is_root) write(*,*)"SNlookup done"
	    endif
	  
            if(survey) then
            	if(Natoms /= 1) then
                  	call halt_program("Natoms should be set to 1 in survey mode. Change Natoms in .inc file & compile again. Aborting.")
                endif
            
            	gamma1s=0d0
                gamma2s=0d0
	    else
	  	gamma1=0d0
        	gamma2=0d0
        	kappa=0d0
        	g1=0d0
        	g2=0d0
            endif

	    GLLhood0=0d0
	    if(simulate) goto 5
	  
            if(survey) then
            	do i=1,nxpix
                  	do j=1,nypix
                        	if(GLeflags(i,j)==0) then
	      				GLLhood0=GLLhood0-log(e1errs(i,j))
	      				GLLhood0=GLLhood0-log(e2errs(i,j))
	    			endif
			enddo
		enddo
                GLLhood0=GLLhood0-((2.*nxpix*nypix)/2.)*log(TwoPi)
	    elseif(arrGal) then
            	m=0
            	do i=1,nxpix
                  	do j=1,nypix
                        	do k=1,ptInPix(i,j,1)
                        		if(GLeflag(ptInPix(i,j,k+1))==0) then
                                                m=m+1
	      					GLLhood0=GLLhood0-log(e1err(ptInPix(i,j,k+1)))
	      					GLLhood0=GLLhood0-log(e2err(ptInPix(i,j,k+1)))
	    				endif
				enddo
			enddo
		enddo
                GLLhood0=GLLhood0-((2.*m)/2.)*log(TwoPi)
	    else
            	m = 0
	  	do k=1,Ngals
	  	    if(GLeflag(k)==0) then
	      		GLLhood0=GLLhood0-log(e1err(k))
	      		GLLhood0=GLLhood0-log(e2err(k))
                        m = m + 1
	    	    endif
	  	enddo
      		GLLhood0=GLLhood0-((2.*m)/2.)*log(TwoPi)
	    endif

	    if(index.ne.0.and.is_root) write(*,*) '       GLLhood0=',GLLhood0	 
            if(index.ne.0.and.is_root) write(*,*)
	  
	endif
	
!-----------------------------------

! SZ stuff:
	  
5	if( SZ == 1 ) then
      
      	    !initialize the spectral index prior routine
      	    do m=1,NSrc
      		if(abs(Src_PriorType(m,4))==7) then
            		!call copySpectralPrior
                        d1=0d0
                        d2=0d0
                  	call init_p_alpha(prior_min,prior_max,nkernel(m),kernel(m,:,:),d1,d2)
		elseif(abs(Src_PriorType(m,4))==12) then
                  	!call makeWaldramGaussianKernel(Src_Prior(m,4,1),Src_Prior(m,4,2))
                  	call init_p_alpha(prior_min,prior_max,nkernel(m),kernel(m,:,:),Src_Prior(m,4,1),Src_Prior(m,4,2))
            	endif
	    enddo
	
	    do m=1,NVisfiles 
        	if(Nvis(m)>large) then
                if(is_root) write(*,*) 'Error: not enough memory allocated for SZ data;'
                if(is_root) write(*,*) '  Nvis(',m,')=',Nvis(m)
                if(is_root) write(*,*) '        large=',large
                if(is_root) write(*,*) 'Increase value of large in .inc file and '
                if(is_root) write(*,*) '  recompile, but note that total memory'
                if(is_root) write(*,*) '  allocated is ~4*Nvisfiles*large^2*8 bytes.'
            		call halt_program()
          	endif
            enddo  

	    SZLhood0=0.d0
   	    pvisr = 0d0
            pvisi = 0d0
            ymap = 0d0

	    if(simulate) goto 10

!           Loop over all data files:
	    Wastage=0

	    if(IncludeCMB==1) then	
              call ReadCovMat	            
  	    else
!             Use noise values in visibility data file in diagonal likelihood	
              do m=1,Nvisfiles
 7                do k=1,Nvis(m)
                     if(SZeflag(m,k)==0) then
              	        SZLhood0=SZLhood0-log(visrms(m,k))
              	     endif
            	  enddo

            	  j=nint(large*(2*large+1)*8.0/1048576.0)
            	  Wastage=Wastage+j
			
                  SZLhood0= SZLhood0-((2*Nvis(m))/2)*log(TwoPi)
	      enddo
      	      write(*,*) '       SZLhood0=',SZLhood0 
	    endif            
        
! 	    Gridding setup:
	
 10    	    trans(1)=-float(nx/2+1)*cell
       	    trans(2)=cell
       	    trans(3)=0.
       	    trans(4)=-float(ny/2+1)*cell
       	    trans(5)=0.
       	    trans(6)=cell
	
       	    uvcell=1.0/(nx*cell*sec2rad)
       	    uvtrans(1)=-float(nx/2+1)*uvcell
       	    uvtrans(2)=uvcell
       	    uvtrans(3)=0.
       	    uvtrans(4)=-float(ny/2+1)*uvcell
       	    uvtrans(5)=0.
       	    uvtrans(6)=uvcell

!      	    Compute SZ conversion factor:

       	    do m=1,Nvisfiles
       		SZscale(m)=y2Jy(cell,nu(m))
		SZscale_BA(m)=dT2di(cell,nu(m))		
            enddo
       
       	    y2jy_nu0=y2Jy_g(cell,nu0(0),f_nu)

!           Set FFTW to work finding the most efficient transform...

!           if(FFTW) then
            call dfftw_plan_dft_2d(fftwplan,nx,ny,arr,arr,FFTW_BACKWARD,FFTW_MEASURE)
!           endif
	
        endif

!-----------------------------------
      
! Initialise working arrays:

        if( Gas == 1 ) then
          do i = 1, n
      	    logr(i) = log10(rmin) + ( log10(rlimit) - log10(rmin) ) * ( i - 1 ) / dble( n - 1 )
	    r(i) = 10.d0 ** logr(i)
	  enddo
        endif

        IF(GasModel==3) THEN
          DO i = 1, n
            logtheta_Planck(i) = log10(thetamin_Planck) + ( log10(thetalimit_Planck) - log10(thetamin_Planck) ) * ( i - 1 ) / dble( n - 1 )
	    theta_Planck(i) = 10.d0 ** logtheta_Planck(i)
	  ENDDO
	
	  !y0_int_Planck = (1.d0/a_GNFW)*Gammafun((-c_GNFW	+1.d0)/a_GNFW) *&
	  !Gammafun((b_GNFW-1.d0)/a_GNFW)/Gammafun((-c_GNFW + b_GNFW)/a_GNFW)
	
	  !Ytot_int_Planck=(1.d0/a_GNFW)*Gammafun((-c_GNFW	+3.d0)/a_GNFW) *&
	  !Gammafun((b_GNFW-3.d0)/a_GNFW)/Gammafun((-c_GNFW + b_GNFW)/a_GNFW)
	
        ENDIF
	
        if(.not.simulate) then
            if(is_root) write(*,*) '         Lhood0=', SZLhood0+GLLhood0 !This does not include Planck likelihood	 
	  	GLLhood=0.d0
	  	SZLhood=0.d0
        	if( GL == 1 ) call ShearLhood(GLLhood)
	  	if( SZ == 1 ) call VisibilitiesLhood(SZLhood)
	  	NullEv=GLLhood+SZLhood !this will be -inf if hyperparameters are included (as it calculates log(alpha) and all params are 0 at this point).
	  	GLNullEv=GLLhood
        if(is_root) write(*,*) '  Null Evidence=',NullEv !This does not include Planck likelihood
	
        endif
      
        end subroutine Initialise

!=======================================================================

	subroutine Welcome
	
	write(*,'(a)')'****************************************************************'
	write(*,'(a)')
	write(*,'(a)')'                           McAdam'
	write(*,'(a)')
	write(*,'(a)')'     Bayesian parameterised fitting of astronomical datasets.'
	write(*,'(a)')
	write(*,'(a)')'               P.J. Marshall et al (April 2003)'
	write(*,'(a)')
	write(*,'(a)')'****************************************************************'
	write(*,'(a)')
	
      return
      end subroutine Welcome

!=======================================================================


end module like
