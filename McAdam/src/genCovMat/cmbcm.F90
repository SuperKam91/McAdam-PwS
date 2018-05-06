program cmbcm
	
      	use globals
      	use inputs
      	use integrals
      
      	implicit none

#ifdef MPI
  	include 'mpif.h'
  	integer errcode
#endif
      
      	double precision, dimension(:,:), allocatable :: covar
	
	double precision, dimension(:), allocatable ::samat	
	
	
      	integer i,j,p,INFO,intPerBatch,indx
      	character*500 incFile
	logical done
	
#ifdef MPI	
	!MPI initializations
	call MPI_INIT(errcode)
	if (errcode/=MPI_SUCCESS) then
     		write(*,*)'Error starting MPI program. Terminating.'
     		call MPI_ABORT(MPI_COMM_WORLD,errcode)
  	end if
	call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, errcode)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_nthreads, errcode)
#else
	mpi_nthreads=1
	my_rank=0
#endif
      
      	! initialise variables
	PI=4d0*atan(1d0)        ! this really does equal PI
      	tcmb=2.726d0

	if( my_rank == 0 ) then
		! read input parameters
		call getarg(1,incFile)
      		!incFile = '/home/cosmos-med/ff235/McAdam_v3/bin/linux/cmbcm/genCMBCovMat_mpi/input.inp'
      		open(unit=input_unit,file=incFile,status='old')
      
		! read no. of files to be read & the no. of annuli
      		read(input_unit,*) nfmax
      		read(input_unit,*) namax
	endif
	
	
#ifdef MPI
    	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	
	call MPI_BCAST(nfmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(namax,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
#endif
	
	
      	nf=nfmax
      	na=namax
      
	! allocate memory
	allocate(freq(nfmax),sigma(nfmax),antwdth(nfmax),f_ra(nfmax),f_dec(nfmax))
	allocate(r(namax+1))
      	
	if( my_rank == 0 ) then
		allocate(cps1d(namax),confps1d(namax),gps1d(namax))
      		allocate(datafn(nfmax),noisefn(nfmax))
		
		cps1d = 0d0
		confps1d = 0d0
		gps1d = 0d0
		
		! read the input file
		call fread_inps(nfmax,nf,nskip,freq,freq0,sigma,antwdth,f_ra,f_dec,ra_ref,dec_ref,antmax, &
      		scaled,scalen,datafn,noisefn,iaddn)
	endif
	
#ifdef MPI
    	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	call MPI_BCAST(imagFlag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(nskip,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(iaddn,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(freq0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(ra_ref,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(dec_ref,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(antmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(scaled,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(scalen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(conf_cl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	
	do i = 1, nfmax
		call MPI_BCAST(freq(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(sigma(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(antwdth(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(f_ra(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(f_dec(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	enddo
#endif


	! 1d PS are in l^2 C_l/2pi in dT/T units; since l=2pi u this means they
	! are in 2pi u^2 S(u) in dT/T units. Calculate conversion factor to
	! give u^2 S(u) in (Jy/sr)^2 at freq0.

	!...convert to (K)^2 units
      	convf = tcmb**2
	!...convert to (Jy/sr)^2 units
      	convf = convf*(dt_to_di(tcmb,freq0)*1d26)**2
	!...divide by 2pi
      	convf = convf/2d0/PI

	! find limits on uv positions
      	if( my_rank == 0 ) call uvlimits(nf,nskip,datafn,uvmin,uvmax)
	
#ifdef MPI
    	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	
	call MPI_BCAST(uvmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(uvmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
#endif

        ! calculate annuli to use	
      	call fannuli(uvmin,uvmax,antmax,r,namax,na,cps1d,confps1d,conf_cl)
      
      	!...count no. of visibilities
      	if( my_rank == 0 ) call count_ndata(nvmax,nf,na,nskip,r,datafn,noisefn)
	
#ifdef MPI
    	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	
	call MPI_BCAST(nvmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
#endif
	
	nv = nvmax
	
	!...allocate memory
	
      	allocate(freqn(nvmax))
	allocate(u(nvmax),v(nvmax))
	allocate(visr(nvmax),visi(nvmax),rms(nvmax))
	
	if( my_rank == 0 ) then
		! read in input data (assuming units are lambda and Jy), plot uv tracks
		call read_data(nvmax,namax,na,nf,nskip,nv,freqn,u,v,visr,visi,rms,freq,r,scaled,scalen,datafn,noisefn,iaddn)

		! read in parameters controlling likelihood calculation
		call fread_params(diagflag,dsflag,nvmax)
		
		! check if the covariance matrix file exists
		IF(dsflag.EQ. 0)THEN 
			 inquire(file=outputfn, exist=done)
		ELSEIF(dsflag.EQ. 1)THEN 
			 inquire(file=outputfn2, exist=done)
		ENDIF
			      
	endif
	
	
#ifdef MPI
    	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	call MPI_BCAST(done,1,MPI_LOGICAL,0,MPI_COMM_WORLD,errcode)
#endif
	
	if( done ) then
		if( my_rank == 0 ) write(*,*)"File "//trim(outputfn)//" exists. Please delete this file if you want to re-generate the covariance matrix"
#ifdef MPI
		call MPI_BARRIER(MPI_COMM_WORLD,errcode)
		call MPI_FINALIZE(errcode)
#endif
		stop
	endif
	
	if( my_rank == 0 ) then

		!...count no. of non-zero correlations
		call count_integrals(nzmax,nvmax,nf,freqn,u,v,freq,antwdth,diagflag)
      
      		write(*,*)nvmax,nzmax
      
      		close(input_unit)

		! make mask for non-zero elements of covariance matrix

		write(*,*) 'Making look-up table of required integrals ...'
	
	endif
	
	
#ifdef MPI
    	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	
	call MPI_BCAST(diagflag,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(dsflag,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(nzmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	do i = 1, nv
		call MPI_BCAST(freqn(i),1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(u(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(v(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(visr(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(visr(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
		call MPI_BCAST(rms(i),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	enddo
#endif
	

	intPerBatch = 100
	intFileUnit = 18
	
	!...memory allocation
	if( my_rank == 0 ) then
		j = intPerBatch
		
		!open the file to store integrals
		open(unit=intFileUnit, file=integralsfn, form='unformatted', status='unknown')
		
		!get the index of the last integral written
		if( getLastIntegralIndex(na, intFileUnit, indx) ) then
			close(intFileUnit)
			open(unit=intFileUnit, file=integralsfn, form='unformatted', status='unknown')
			call moveIntegralFilePointer(na, intFileUnit, indx)
		endif
	else
		j = intPerBatch
	endif
	
#ifdef MPI
    	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	
	call MPI_BCAST(indx,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
#endif

      	allocate(irn(j),icn(nvmax+1))
	allocate(freqi(j),freqj(j))
	allocate(integral1r(j,namax),integral1i(j,namax),integral2r(j,namax),integral2i(j,namax))
      	call makeintegrals(nv,nf,na,freqn,u,v,freq,antwdth,sigma,f_ra,f_dec,ra_ref,dec_ref,r,diagflag,intPerBatch,intFileUnit,indx+1)
	if( my_rank == 0 ) close(intFileUnit)
      	
	
	if( my_rank == 0 ) then
      		allocate(cm_rr(nzmax),cm_ii(nzmax),cm_ri(nzmax),cm_ir(nzmax))
		
		open(unit=intFileUnit, file=integralsfn, form='unformatted', status='old')
		
		deallocate(irn)
		
		allocate(irn(nzmax))
		
      		call makecm(nv,na,rms,cps1d,gps1d,confps1d,sp,convf,freq0,intFileUnit)
		
		close(intFileUnit)
      
      		allocate(covar(2*nv,2*nv))
      
      		!write out the covariance matrix in lower triangular form
      		covar=0.d0

      		do j=1,nv
        		do p=icn(j),icn(j+1)-1
        			i=irn(p)
            			!real real part
            			covar(i,j)=cm_rr(p)
            			!imaginary imaginary part
            			covar(i+nv,j+nv)=cm_ii(p)
            			!real imaginary part
            			covar(i,j+nv)=cm_ri(p)
            			!imaginary real part
            			covar(i+nv,j)=cm_ir(p)
	  		enddo
		enddo
		
	IF(dsflag.EQ. 0)THEN      
	      
      		!compute the Cholesky decomposition of the covariance matrix 'C'
      		!C = L^T.L
      		!uses the LAPACK routine DPOTRF
		
		write(*,*)"Performing Cholesky decomposition of the covariance matrix"

      		call DPOTRF('L',2*nv,covar,2*nv,INFO)

		if( INFO /= 0 ) then
			write(*,*)"problem in performing Cholesky decomposition of the covariance matrix."
			write(*,*)"Aborting!"
#ifdef MPI
			call MPI_ABORT(MPI_COMM_WORLD, errcode)
#endif
			stop
		endif

		!write(*,*)INFO
      
      		!open the file to store the covariance matrix
      		open(unit=18, file=outputfn, form='unformatted', status='unknown')
      		write(18)2*nv
      
      		do i=1,2*nv
        		write(18)covar(i,1:i)
		enddo
      
      		close(18)
               deallocate(covar)
      ENDIF
      	       
             deallocate(irn,icn)

	     
	IF(dsflag.EQ.1)THEN	
	
	     DO j=1,2*nv
		  DO i=j,2*nv
		        IF(covar(i,j).NE.0.) THEN
			   nz=nz+1
			 ENDIF
		    ENDDO 
	     ENDDO	
		      
              write(*,*)'nz=  ', nz
	      
	      ALLOCATE(samat(nz))
	      ALLOCATE(irn(nz))
	      ALLOCATE(icn(nz))	
			      
              nz=0
	      samat=0d0
	      irn=0
	      icn=0	      
	      
	     DO j=1,2*nv
	          icn(j)=nz+1
		  DO i=j,2*nv
		        IF(covar(i,j).NE.0.) THEN
			   nz=nz+1
			   samat(nz)=covar(i,j)
	                    irn(nz)=i
			ENDIF
		  ENDDO
	     ENDDO
	     icn(2*nv+1)=nz+1
      	 	  open(unit=69, file=outputfn2)
		  
write(69,*)size(irn),size(icn),size(samat) ,nz   
	     
	     DO i=1, nz 
		  write(69,*)irn(i),icn(i),samat(i) 
             ENDDO
	     close(69)
	       deallocate(covar)
	     deallocate(samat)
	ENDIF	  
    

	       
      		deallocate(freq,sigma,antwdth,f_ra,f_dec)
      		deallocate(datafn,noisefn)
		deallocate(r)
		deallocate(cps1d,confps1d,gps1d)
      		deallocate(irn,icn)
		deallocate(freqi,freqj)
      		deallocate(cm_rr,cm_ii,cm_ri,cm_ir)
		deallocate(integral1r,integral1i,integral2r,integral2i)
      		deallocate(freqn)
		deallocate(u,v)
		deallocate(visr,visi,rms)
	endif
	
	
#ifdef MPI
		
	if( my_rank /= 0 ) then
		deallocate(freq,sigma,antwdth,f_ra,f_dec,r)
      		deallocate(freqn,u,v,visr,visi,rms)
		deallocate(irn,icn,freqi,freqj,integral1r,integral2r,integral1i,integral2i)
	endif

	if( my_rank == 0 ) call MPI_ABORT(MPI_COMM_WORLD, errcode)
	
	call MPI_BARRIER(MPI_COMM_WORLD,errcode)

	call MPI_FINALIZE(errcode)

#endif
      
end program cmbcm
