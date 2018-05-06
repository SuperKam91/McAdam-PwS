module ReadWrite

      use params
      use utilities
      use constants
      use telescope1
      use MassModels
#ifdef MPI
      use Nested, only: MPI_COMM_WORLD
#endif
      use mkl_dss
      use fits_handling
      use pws_wrapper

    integer, parameter :: params_unit = 21
    ! Allow strings to be really long to allow for arrays of file names
    character(len=2), parameter :: comment='#!'
    character(len=3), parameter :: equals='=:/'
    integer :: io_stat   ! check to see if we've reached the end of the file
    logical :: is_formatted
    integer, parameter :: bufLen = 20

contains
!=======================================================================
     subroutine ReadCovMat

     implicit none
     integer iostatus, i, j, k, NSZcovmats, nTotVis, Wastage
     double precision row(large*Nvisfiles*2)
     character(len=100)  nTotVis_str, i_str
     integer n_irn, n_icn, n_samat, nz_ichol, n_ichol
     !integer, dimension(:), allocatable     :: irn,icn
     !double precision,  dimension(:), allocatable     :: samat
     external MKL_CVT_TO_NULL_TERMINATED_STR
     double precision, allocatable :: statOut(:)
     character(len=15) statIn
     integer :: error, perm(1), buff(bufLen)

     nTotVis = sum(Nvis)
     NSZcovmats=0
     Wastage=0
     if(dsflag.eq. 0)then      
			
         if(is_root) write(*,*) 'dense matrix approach'
!	Read in Cholesky decomposition of covariance matrices, one row at
!   	a time, and compute the determinant by the product of the
!   	diagonal values. Also store the matrix for future use:  
	open(unit=9,form='unformatted',file=covmatfile,status='old',iostat=iostatus)
	if( iostatus /= 0 ) then
	    call halt_program("ERROR: could not open the file "//trim(covmatfile)//". Aborting")
	endif
	read(9) i
	if(i /= 2*nTotVis) then
                write(nTotVis_str, '(I0)') nTotVis
                write(i_str, '(I0)') i
                if(is_root) write(*,*)  'Array size mismatch-expecting '//trim(nTotVis_str)
                if(is_root) write(*,*)  '                      received  '//trim(i_str)
                if(is_root) write(*,*)  'Possible causes: '
                if(is_root) write(*,*)  '  1) vis.fits/.LCM file mismatch'
                if(is_root) write(*,*)  '  2) byte-swapped .LCM file'
            	call halt_program()
	endif
	k=1
	do i=1,2*nTotVis
	    read(9)row(1:i)
	    SZLhood0=SZLhood0+log(row(i))
	    do j=1,i
                LCM(k)=row(j)
              	k=k+1
            enddo
	enddo
	close(9)
            
        i=nint(nTotVis*(2*nTotVis+1)*8.0/1048576.0)
        write(i_str, '(I0)') i
        !write(*,'(a,i3,a)') ' Read in covariance matrix (',i,' Mb)'
        if(is_root) write(*,*) 'Read in covariance matrix ('//trim(i_str)//' Mb)'
        if(is_root) write(*,*) 'from file '//trim(covmatfile)

        NSZcovmats=NSZcovmats+1
        j=nint(large*Nvisfiles*(2*large*Nvisfiles+1)*8.0/1048576.0)
        Wastage=Wastage+(j-i)

        SZLhood0=-((2.d0*NtotVis)/2.d0)*log(TwoPi)-SZLhood0
        if(is_root) write(*,*) '       SZLhood0=',SZLhood0

     elseif(dsflag.eq.1)then
	
        if(is_root) write(*,*) 'sparse matrix approach'
        open(unit=9,file=covmatfile,status='old',iostat=iostatus)
        read(9,*)n_irn , n_icn , n_samat ,nz_ichol
 	allocate(irn(nz_ichol )) 
	allocate(icn(nz_ichol )) 	
	allocate(samat(nz_ichol  ))	   
	
	do i=1, nz_ichol	     
	    read(9,*)irn(i),icn(i),samat(i)
	enddo
        close(9)	  
	  	 
	n_ichol=2*nTotVis

       ! New way using intel sparse matrix routines (YCP, 28/9/2015)
       ! Initialize the solver.
        error = dss_create( handle, MKL_DSS_DEFAULTS )
        if (error /= MKL_DSS_SUCCESS) call halt_program('dss_create failed')
        ! Define the non-zero structure of the matrix.  Note that I've swapped the column and row vectors because this routine expects SLAP row format while the old one expects SLAP column format, but it's symmetric so this is fine
        error = dss_define_structure( handle, MKL_DSS_SYMMETRIC, icn(1:n_ichol+1), n_ichol, &
      & n_ichol, irn, nz_ichol )
        if (error /= MKL_DSS_SUCCESS) call halt_program('dss_define_structure failed')
        ! Reorder the matrix, find a permutation matrix to minimize fill-in during solve.  Note perm is not used unless you change the options from defaults.
        perm=0
        error = dss_reorder( handle, MKL_DSS_DEFAULTS, perm )
        if (error /= MKL_DSS_SUCCESS) call halt_program('dss_reorder failed')
        ! Factor the matrix.
	!write(*,*) 'reordered structure'
        error = dss_factor_real( handle, MKL_DSS_POSITIVE_DEFINITE, samat )
        if (error /= MKL_DSS_SUCCESS) call halt_program('dss_factor_real failed')
	!write(*,*) 'factored'
        ! Get the determinant of the matrix
        allocate(statOut( 2 ) )
        statIn = 'determinant'
        call mkl_cvt_to_null_terminated_str(buff,bufLen,statIn);
	!write(*,*) 'mkl cvt'
        error = dss_statistics(handle, MKL_DSS_DEFAULTS, buff, statOut )
	!write(*,*) 'got stats'
        if (error /= MKL_DSS_SUCCESS) call halt_program('dss_statistics failed')
        !WRITE(*,"('pow of determinant is '(5F10.3))") ( statOut(1) )
        !WRITE(*,"('base of determinant is '(5F10.3))") ( statOut(2) )
        SZLhood0=statOut(1)+dlog10(statOut(2))
        SZLhood0=SZLhood0/dlog10(exp(1d0))
	SZLhood0=-((2.d0*NtotVis)/2.d0)*log(TwoPi)-SZLhood0/2d0
        if(is_root) write(*,*) '       SZLhood0=',SZLhood0

     endif  
	   			
     if(NSZcovmats>0) then
         i=Nvisfiles*nint(large*Nvisfiles*(2*large*Nvisfiles+1)*8.0/1048576.0)
         write(i_str, '(I0)') i
         if(is_root) write(*,*) 'Memory allocated='//trim(i_str)//' Mb'
         write(i_str, '(I0)') Wastage
         if(is_root) write(*,*) '  ('//trim(i_str)//' Mb left unused).'
         if(is_root) write(*,*)
    endif

end subroutine ReadCovMat

!=======================================================================
	subroutine ReadInData
	
	implicit none
	
	integer i,j,m,iend
	double precision ra,dec,ra2,dec2
	double complex vis_buffer
      double precision weight_buffer,uv_buffer(2),rms_buffer
      integer baseline,nerr
      character(len=100) string
      double precision OBSRA,OBSDEC,CRVAL4,CRVAL5,CRVAL6,PSCAL1,PSCAL2 !,PSCAL4,PZERO4,PSCAL5,PZERO5
      character(len=100) TELESCOP !,OBJECT,INSTRUME

!-----------------------------------------------------------------------
      
! Lensing data:

	if(GL==1) then
		if(survey) then
      		call readcat_survey(nxpix,nypix,e1s,e2s,e1errs,e2errs,z_ss,nerr,GLdatafile)
                  i=nxpix*nypix
		elseif(arrGal) then
            	call readcat_im(m,x,y,e1,e2,e1err,e2err,z_s,ptInPix,nerr,GLdatafile)
                  i=m
            else
      		call readcat(Ngals,x,y,e1,e2,e1err,e2err,z_s,nerr,GLdatafile)
                  i=ngals
		endif
            
            string = GLdatafile
      	do iend=100,1,-1
        		if(string(iend:iend).ne.' ') goto 10
      	enddo 
10      if(is_root) write(*,*) 'Read in ',i,' complex ellipticities '
        if(is_root) write(*,*) 'from file ',string(1:iend),'; '
        if(is_root) write(*,*) nerr,' were given zero weight.'
        if(is_root) write(*,*)
	endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Planck data:
        if ( PL == 1 ) then
           ! initialise PwS parameters from file
           write(*,*) 'Initialising PwS parameters from ', trim(PLdatafile)
           call initializePwS(PLdatafile)
           write(*,*) 'Getting Planck data patch'
           call get_pws_patch(PLPatch)
           write(*,*) 'Finished getting Planck data patch'
        endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! SZ data:

	if( SZ == 1 ) then
	
! Loop over all files;       
if(is_root) write(*,*) 'Reading SZ data'
	do j=1,Nvisfiles
        if(is_root) write(*,*) 'opening '//trim(visdatafile(j))
        call phitsread_open(visdatafile(j), CRVAL4, CRVAL5, CRVAL6, OBSRA, OBSDEC, TELESCOP, PSCAL1, PSCAL2)
        !write(*,*) trim(visdatafile(j)), CRVAL4, CRVAL5, CRVAL6, OBSRA, OBSDEC, trim(TELESCOP), PSCAL1, PSCAL2
        nu(j) = CRVAL4
        !calculate the sigma of the primary beam in arcsec
        !pb_sig(j) = 1.2*180.*clight*3600./(nu(j)*3.7*Pi*2.35482)
        pb_x(j) = CRVAL5
        pb_y(j) = CRVAL6
        ra = map_x*deg2rad
        dec = map_y*deg2rad
        ra2 = pb_x(j)*deg2rad
        dec2 = pb_y(j)*deg2rad
        call calc_ra_dec_off(ra,dec,ra2,dec2,pb_x_offset(j),pb_y_offset(j)) 	  
! NB. RA-Dec is left handed, need -ve sign here. 
        pb_x_offset(j) = -pb_x_offset(j)*rad2sec
        pb_y_offset(j) = pb_y_offset(j)*rad2sec

        ph_x(j) = OBSRA
        ph_y(j) = OBSDEC
        ra = map_x*deg2rad
        dec = map_y*deg2rad
        ra2 = ph_x(j)*deg2rad
        dec2 = ph_y(j)*deg2rad
        call calc_ra_dec_off(ra,dec,ra2,dec2,ph_x_offset(j),ph_y_offset(j)) 	  
! NB. RA-Dec is left handed, need -ve sign here. 
        ph_x_offset(j) = -ph_x_offset(j)*rad2sec
        ph_y_offset(j) = ph_y_offset(j)*rad2sec
        telescope(j) = TELESCOP
        
! Set up beam
        pb_sig(j) = 1.2*180.*clight*3600./(nu(j)*3.7*Pi*2.35482)
        ! uv-limits not used at the moment
        lmin(j) = 0.d0
        lmax(j) = 0.d0
        if(telescope(j)(1:2)=='AM'.and.telescope(j)(1:6)/='AMI-LA') then
           pb_sig(j)=pb_ami(1)/(nu(j)*1.d-9)+pb_ami(2)
           pb_sig(j)=pb_sig(j)*60. ! arcsec
           if(is_root.and.verbose) write(*,*) trim(telescope(j))//' primary beam =',pb_sig(j)
        elseif(telescope(j)(1:2)=='LA'.or.telescope(j)(1:6)=='AMI-LA') then
           pb_sig(j)=pb_la(1)/(nu(j)*1.d-9)+pb_la(2)
           pb_sig(j)=pb_sig(j)*60. ! arcsec
           if(is_root.and.verbose) write(*,*) trim(telescope(j))//' primary beam =',pb_sig(j)
        endif

        if(GCOUNT.ne.Nvis(j)) then
            if(is_root) write(*,*) 'Incompatible data sizes:'
            if(is_root) write(*,*) '  In data file Nvis = ',GCOUNT
            if(is_root) write(*,*) '  In McAdam.inc Nvis = ',Nvis(j)
            if(is_root) write(*,*)
          call halt_program
        endif

        m = 0
        do i=1,Nvis(j)

          call phitsread_one(i,vis_buffer,uv_buffer,weight_buffer,rms_buffer,baseline, CRVAL4, PSCAL1, PSCAL2)

          u(j,i) = 1.d0*uv_buffer(1)
          v(j,i) = 1.d0*uv_buffer(2)
          visr(j,i) = dble(vis_buffer)
          visi(j,i) = aimag(vis_buffer)
          visrms(j,i) = 1.d0/sqrt(weight_buffer)
          viswt(j,i) = 1.d0*weight_buffer

!       if(i.lt.11) write(*,*) i,u(j,i),v(j,i),visrms(j,i),viswt(j,i)

          if(weight_buffer.ne.0.0) then
            SZeflag(j,i) = 0
          else
            SZeflag(j,i) = 1
            m = m + 1
          endif

        enddo

        call phitsread_close

        do iend=100,1,-1
          if(visdatafile(j)(iend:iend).ne.' ') exit
        enddo

	enddo
      
      endif

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      return
      end subroutine ReadInData

!=======================================================================

	subroutine writecat(m,x,y,e1,e2,e1err,e2err,z_s,filename)

	implicit none
	
	integer m
	double precision x(m),y(m),e1(m),e2(m),e1err(m),e2err(m),z_s(m)
	
	integer i
	character(len=100) filename
	
!       write(*,*) 'In writecat:'
!       write(*,*) ' filename = ',filename
!       write(*,*) ' m = ',m
!       write(*,*) ' x(1) = ',x(1)
!       write(*,*) ' y(1) = ',y(1)
!       write(*,*) ' e1(1) = ',e1(1)
!       write(*,*) ' e1err(1) = ',e1err(1)
!       write(*,*) ' e2(1) = ',e2(1)
!       write(*,*) ' e2err(1) = ',e2err(1)
!       write(*,*) ' zs(1) = ',zs(1)
      
	open(unit=23,file=filename,form='formatted',status='unknown')
	write(23,'(a)')'Mock galaxy catalogue from McAdam'
	write(23,*) m
	write(23,*) 
	write(23,'(a)')'         x        y       e1    e1err       e2    e2err     sigc'
	write(23,*)
	do i=1,m
	  write(23,5) x(i),y(i),e1(i),e1err(i),e2(i),e2err(i),z_s(i)
	enddo
 5    format(4x,f7.1,3x,f7.1,3x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5,3x,f7.1)	
	close(23)  
			
	return
	end subroutine writecat

!=======================================================================

	subroutine writecat_survey(nx,ny,gm1s,gm2s,err,zsrc,filename)

	implicit none
	
	integer nx,ny
	double precision gm1s(nx,ny),gm2s(nx,ny),err,zsrc
	
	integer i,j
	character(len=100) filename
      
	open(unit=23,file=filename,form='formatted',status='unknown')
	write(23,'(a)')'Mock shear map from McAdam'
	write(23,*) nx,ny
	write(23,*) 
	write(23,'(a)')'	i	j	gamma1	e1err		gamma2	e2err    zs'
	write(23,*)
	do i=1,nx
      	do j=1,ny
	  		write(23,5) (2.*i-1.-nxpix)*pix_sizex/2.,(2.*j-1.-nypix)*pix_sizey/2.,gm1s(i,j),err,gm2s(i,j),err,zsrc
		enddo
	enddo
 !5    format(4x,i6,3x,i6,3x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5)	
 5    format(2x,f8.3,3x,f8.3,3x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5)	
	close(23)  
			
	return
	end subroutine writecat_survey

!=======================================================================

	subroutine writekappa(m,x,y,kap,filename)

	implicit none
	
	integer m
	double precision x(m),y(m),kap(m)
	
	integer i
	character(len=100) filename
      
	open(unit=23,file=filename,form='formatted',status='unknown')
!	write(23,'(a)')'True kappa'
!	write(23,*) nx,ny
!	write(23,*) 
!	write(23,'(a)')'	x	y	kappa'
!	write(23,*)
	do i=1,m
		write(23,5) x(i),y(i),kap(i)
	enddo
 5    format(4x,f7.1,3x,f7.1,3x,f12.9)	
	close(23)  
			
	return
	end subroutine writekappa

!=======================================================================

	subroutine writekappa_survey(nx,ny,kap,filename)

	implicit none
	
	integer nx,ny
	double precision kap(nx,ny)
	
	integer i,j
	character(len=100) filename
      
	open(unit=23,file=filename,form='formatted',status='unknown')
!	write(23,'(a)')'True kappa'
!	write(23,*) nx,ny
!	write(23,*) 
!	write(23,'(a)')'	i	j	kappa'
!	write(23,*)
	do i=1,nx
      	do j=1,ny
	  		write(23,5) i,j,kap(i,j)
		enddo
	enddo
 5    format(4x,i6,3x,i6,3x,f12.9)	
	close(23)  
			
	return
	end subroutine writekappa_survey

!=======================================================================

	subroutine readcat(m,x,y,e1,e2,e1err,e2err,z_s,nerr,filename)

	implicit none
	
	integer m,nerr
	double precision x(m),y(m),e1(m),e2(m),e1err(m),e2err(m),z_s(m)
	
	integer i
	character(len=100) filename
	
      open(unit=23,file=filename,form='formatted',status='unknown')
      read(23,*)
      read(23,*) i
      read(23,*) 
      read(23,*)
      read(23,*)
      if(i.ne.m) then
          if(is_root) write(*,*) 'Incompatible data sizes in readcat:'
          if(is_root) write(*,*) '  In data file Ngals = ',i
          if(is_root) write(*,*) '  In McAdam.inc Ngals = ',m
        call halt_program
      endif   
	
      zsmax=0.
      zsmin=10000.
      nerr=0
	do i=1,m
	  read(23,*) x(i),y(i),e1(i),e1err(i),e2(i),e2err(i),z_s(i)
        
        !x(i)=(2*x(i)-1-nxpix)*pix_sizex/2.
        !y(i)=(2*y(i)-1-nypix)*pix_sizey/2.
        
        !x(i)=x(i)+pix_sizex/2.
        !y(i)=y(i)+pix_sizey/2.
        
        if(e1err(i)/=0. .and. e2err(i)/=0.) then
!		if(e1err(i)<=sigmaobs .or. e2err(i)<=sigmaobs) then
!	      	write(*,*) 'Obs. error overestimated in McAdam.inc...'
!                  write(*,*)e1err(i),e2err(i),sigmaobs
!	      	stop
!	    	endif
	    	GLeflag(i) = 0
	    	wt1(i) = 1.0/(e1err(i)*e1err(i))
	    	wt2(i) = 1.0/(e2err(i)*e2err(i))
	  else
	    	GLeflag(i) = 1
	    	wt1(i) = 0.0
	    	wt2(i) = 0.0
	    	nerr = nerr + 1
	  endif  
        
        if(varyzs==1) then
        	!z_s(i)=0.
	  else
        	if(z_s(i)<zsmin) zsmin=z_s(i)
        	if(z_s(i)>zsmax) zsmax=z_s(i)
	  endif
        
	enddo
      
      if(varyzs==1 .and. nzsplanes==1) then
      	vary_zs=0
      elseif(zsmin/=zsmax) then
      	vary_zs=1
	endif
      
	close(23)  
      
      return
      end subroutine readcat
      
!=======================================================================

	subroutine readcat_survey(nx,ny,e1s,e2s,e1errs,e2errs,z_ss,nerr,filename)

	implicit none
	
	integer nx,ny !no. of pixels in x & y
      integer nerr
	double precision e1s(nx,ny),e2s(nx,ny),e1errs(nx,ny),e2errs(nx,ny),z_ss(nx,ny)
	
	integer i,j,iostatus
	character(len=100) filename
      double precision d1,d2,d3,d4,d5
	
      open(unit=23,file=filename,form='formatted',status='unknown')
      read(23,*)
      read(23,*) i,j
      read(23,*) 
      read(23,*)
      read(23,*)
      
      if (i/=nx .or. j/=ny) then
          if(is_root) write(*,*) 'Incompatible data sizes in readcat_survey:'
          if(is_root) write(*,*) '  In data file Ngals = ',i,'X',j
          if(is_root) write(*,*) '  In McAdam.inc Ngals = ',nx,'X',ny
        call halt_program
      endif
      
      GLeflags=1
	
      zsmax=0.
      zsmin=10000.
      
      nerr=nx*ny
	do
	  read(23,*,IOSTAT=iostatus) i,j,d1,d2,d3,d4,d5
        
        !i=i+1
        !j=j+1
        
        !end of file?
        if(iostatus<0) exit
        
        if(i<1 .or. i>nx .or. j<1 .or. j>nx) then
            if(is_root) write(*,*)"error in readcat_survey, pixel value out of bounds"
            if(is_root) write(*,*)i,j,nx,ny
            call halt_program
	  endif
        
        e1s(i,j)=d1
        e2s(i,j)=d3
        e1errs(i,j)=d2
        e2errs(i,j)=d4
        z_ss(i,j)=d5
        
        if (e1errs(i,j)/=0. .and. e2errs(i,j)/=0.) then
		if (e1errs(i,j)<=sigmaobs .or. e2errs(i,j)<=sigmaobs) then
            if(is_root) write(*,*) '   Obs. error overestimated in McAdam.inc...'
	      	call halt_program
	    	endif
	    	GLeflags(i,j) = 0
	    	wt1s(i,j) = 1.0/(e1errs(i,j)*e1errs(i,j))
	    	wt2s(i,j) = 1.0/(e2errs(i,j)*e2errs(i,j))
            nerr=nerr-1
	  else
	    	wt1s(i,j) = 0.0
	    	wt2s(i,j) = 0.0
	  endif
        
        if(varyzs==1) then
        	!z_ss(i,j)=0.
	  else
        	if(z_ss(i,j)<zsmin) then
        		zsmin=z_ss(i,j)
        	elseif(z_ss(i,j)>zsmax) then
        		zsmax=z_ss(i,j)
	  	endif
	  endif
        
	enddo
      
      if(varyzs==1 .and. nzsplanes==1) then
      	vary_zs=0
      elseif(zsmin/=zsmax) then
      	vary_zs=1
	endif
      
	close(23)  
      
      return
      end subroutine readcat_survey
      
!=======================================================================

	subroutine readcat_im(m,x,y,e1,e2,e1err,e2err,z_s,ptInPix,nerr,filename)

	implicit none
	
	integer m,ptInPix(:,:,:),nerr
	double precision x(ngals),y(ngals),e1(ngals),e2(ngals),e1err(ngals),e2err(ngals),z_s(ngals)
	
	integer i,j,k,iostatus
	character(len=100) filename
	
      open(unit=23,file=filename,form='formatted',status='unknown')
	
      zsmax=0.
      zsmin=10000.
      ptInPix=0
      i=0
      nerr=0
	do
        i=i+1
	  read(23,*,IOSTAT=iostatus) x(i),y(i),e1(i),e2(i),e1err(i),z_s(i)
        e2err(i)=e1err(i)
        
        if(iostatus<0) then
        	i=i-1
            m=i
            exit
	  endif
        
        j=int((survey_xl/pix_sizex)*x(i))+1-nxmin+1
        k=int((survey_yl/pix_sizey)*y(i))+1-nymin+1
        
        if(j>nxpix .or. k>nypix .or. j<1 .or. k<1) then
        	i=i-1
            cycle
	  endif
        
        ptInPix(j,k,1)=ptInPix(j,k,1)+1
        ptInPix(j,k,ptInPix(j,k,1)+1)=i
        
        if(ptInPix(j,k,1)>50) then
        	call halt_program("more than 50 galaxies in a pixel")
	  endif
        
        x(i)=-survey_xl/2.+survey_xl*x(i)
        y(i)=-survey_yl/2.+survey_yl*y(i)
        
        if (e1err(i)/=0. .and. e2err(i)/=0.) then
		if (e1err(i)<=sigmaobs .or. e2err(i)<=sigmaobs) then
            if(is_root) write(*,*) 'Obs. error overestimated in McAdam.inc...'
            if(is_root) write(*,*)e1err(i),e2err(i),sigmaobs
                call halt_program
	    	endif
	    	GLeflag(i)=0
	    	wt1(i)=1.0/(e1err(i)*e1err(i))
	    	wt2(i)=1.0/(e2err(i)*e2err(i))
	  else
	    	GLeflag(i)=1
	    	wt1(i)=0.0
	    	wt2(i)=0.0
	    	nerr=nerr+1
	  endif
        
        if(varyzs==1) then
        	!z_s(i)=0.
	  else
        	if(z_s(i)<zsmin) then
        		zsmin=z_s(i)
        	elseif(z_s(i)>zsmax) then
        		zsmax=z_s(i)
	  	endif
	  endif
        
	enddo
      
      if(varyzs==1 .and. nzsplanes==1) then
      	vary_zs=0
      elseif(zsmin/=zsmax) then
      	vary_zs=1
	endif
      
	close(23)  
      
      return
      end subroutine readcat_im
      
!=======================================================================

    subroutine read_params(file_name)
        implicit none
        
        character(len=*)   :: file_name !> The name of the file
        integer            :: ifile, iatom, ipar, isrc, i, j
        character(len=3)   :: ifiles, iatoms, ipars, isrcs
        double precision   :: prior_lims(2)
        character(len=100) :: out_file

        i = index(file_name, '.', .true.)

        if (file_name(i+1:i+3)=='dat') then
            if(is_root) write(*,*) 'Assuming this is an unformatted file'
          open(unit=params_unit, form='unformatted', file=file_name, action='read',iostat=io_stat)
          is_formatted=.false.
        else
            if(is_root) write(*,*) 'Assuming this is a formatted file'
          open(unit=params_unit, form='formatted', file=file_name, action='read',iostat=io_stat)
          is_formatted=.true.
          out_file=file_name(1:i)//'dat'
        endif

        verbose = get_logical( 'verbose', .false.)

        if (verbose.and.is_root) write(*,*) 'Parameters read from file are: '
        ModelClass = get_integer( 'ModelClass')
        if (ModelClass .ne. 1) call halt_program('Only ModelClass=1 (clusters) is available')

        Atoms = get_integer( 'Atoms')
        if (Atoms .lt. 1) call halt_program('No objects to fit')

        NAtoms = get_integer( 'NAtoms')
        if (NAtoms == 0 ) call halt_program('No objects to fit')

        allocate(z(NAtoms))
        allocate(z_old(NAtoms))
        allocate(z_PriorType(NAtoms))
        allocate(z_Prior(NAtoms,3))
        allocate(BetaStyle(NAtoms))
        allocate(TStyle(NAtoms))
        allocate(NFWstyle(NAtoms))
        allocate(SISstyle(NAtoms))
        allocate(CPLstyle(NAtoms))
        z=0d0

! Cluster geometry
        GeoModel = get_integer( 'GeoModel')
        NGeoPars = 2*GeoModel
        allocate(Geo_PriorType(NAtoms,NGeoPars))
        allocate(Geo_Prior(NAtoms,NGeoPars,3))
        allocate(Geo_Tri_Prior(NAtoms,2,3))
        allocate(GeoPars(NAtoms,NGeoPars))
        allocate(GeoPars_old(NAtoms,NGeoPars))
        do iatom = 1, NAtoms
          write(iatoms,'(I0)') iatom
          do ipar = 1, NGeoPars
            write(ipars,'(I0)') ipar
            Geo_PriorType(iatom,ipar) = get_integer( 'Geo_PriorType('//trim(iatoms)//','//trim(ipars)//')')
            call get_doubles( 'Geo_Prior('//trim(iatoms)//','//trim(ipars), 3, Geo_Prior(iatom, ipar,:), .false.)
          enddo
        enddo
        GeoPars=0d0

! Cluster mass
        MassModel = get_integer( 'MassModel')
        NMassPars = get_integer( 'NMassPars')
        allocate(Mass_PriorType(NAtoms,NMassPars))
        allocate(Mass_Prior(NAtoms,NMassPars,3))
        allocate(MassPars(NAtoms,NMassPars))
        allocate(MassPars_old(NAtoms,NMassPars))
        do iatom = 1, NAtoms
          write(iatoms,'(I0)') iatom
          NFWstyle(iatom) = get_integer( 'NFWstyle('//trim(iatoms)//')')
          do ipar = 1, NMassPars
            write(ipars,'(I0)') ipar
            Mass_PriorType(iatom,ipar) = get_integer( 'Mass_PriorType('//trim(iatoms)//','//trim(ipars)//')')
            call get_doubles( 'Mass_Prior('//trim(iatoms)//','//trim(ipars), 3, Mass_Prior(iatom, ipar,:), .false.)
          enddo
        enddo
        MassPars=0d0

! Cluster gas
        GasModel = get_integer( 'GasModel')
        NGasPars = get_integer( 'NGasPars')
        allocate(Gas_PriorType(NAtoms,NGasPars))
        allocate(Gas_Prior(NAtoms,NGasPars,3))
        allocate(GasPars(NAtoms,NGasPars))
        allocate(GasPars_old(NAtoms,NGasPars))
        do iatom = 1, NAtoms
          write(iatoms,'(I0)') iatom
          BetaStyle(iatom) = get_integer( 'BetaStyle('//trim(iatoms)//')', 0)
          TStyle(iatom) = get_integer( 'TStyle('//trim(iatoms)//')', 0)
          do ipar = 1, NGasPars
            write(ipars,'(I0)') ipar
            Gas_PriorType(iatom,ipar) = get_integer( 'Gas_PriorType('//trim(iatoms)//','//trim(ipars)//')')
            call get_doubles( 'Gas_Prior('//trim(iatoms)//','//trim(ipars), 3, Gas_Prior(iatom, ipar,:), .false.)
          enddo
        enddo
        GasPars=0d0

! Cluster temperature
        TModel = get_integer( 'TModel')
        NTPars = get_integer( 'NTPars')
        if (NTPars > 0) then
          allocate(T_PriorType(NAtoms,NTPars))
          allocate(T_Prior(NAtoms,NTPars,3))
          allocate(TPars(NAtoms,NTPars))
          allocate(TPars_old(NAtoms,NTPars))
          do iatom = 1, NAtoms
            write(iatoms,'(I0)') iatom
            do ipar = 1, NTPars
              write(ipars,'(I0)') ipar
              T_PriorType(iatom,ipar) = get_integer( 'T_PriorType('//trim(iatoms)//','//trim(ipars)//')')
              call get_doubles( 'T_Prior('//trim(iatoms)//','//trim(ipars), 3, T_Prior(iatom, ipar,:), .false.)
            enddo
          enddo
        else
          allocate(T_PriorType(NAtoms,1))
          allocate(T_Prior(NAtoms,1,3))
          allocate(TPars(NAtoms,1))
          allocate(TPars_old(NAtoms,1))
        endif
        TPars=0d0

! Cluster redshift
        do iatom = 1, NAtoms
          write(iatoms,'(I0)') iatom
          z_PriorType(iatom) = get_integer( 'z_PriorType('//trim(iatoms)//')')
          call get_doubles( 'z_Prior('//trim(iatoms), 3, z_Prior(iatom, :), .false.)
        enddo

! Derived parameters
        aux_dim = get_integer( 'aux_dim')
        allocate(aux(NAtoms, aux_dim))
        do iatom = 1, NAtoms
          do j = 1, aux_dim
            aux(iatom,j)=0d0
          enddo
        enddo

! SZ setup
        SZ = get_integer( 'SZ', 1)
        if (SZ .ne. 0 ) then
          Gas = get_integer( 'Gas', 1)
          Temperature = get_integer( 'Temperature', 1)
          MMass = get_integer( 'MMass')
          Nvisfiles = get_integer( 'Nvisfiles')
          allocate(visdatafile(Nvisfiles))
          ! Try two methods for specifying visdatafiles
          do ifile = 1, Nvisfiles
            write(ifiles, '(I0)') ifile
            visdatafile(ifile) = get_string( 'visdatafile('//trim(ifiles)//')', .true., .true., '')
            if (visdatafile(ifile) == '') exit
          enddo
          if (visdatafile(1) == '') then
            call get_strings( 'visdatafile', Nvisfiles, visdatafile)
          endif
          allocate(Nvis(Nvisfiles))
          call get_integers( 'Nvis', Nvisfiles, Nvis)
          large = get_integer( 'large')
          allocate(u(Nvisfiles,large))
          allocate(v(Nvisfiles,large))
          allocate(visr(Nvisfiles,large))
          allocate(visi(Nvisfiles,large))
          allocate(visrms(Nvisfiles,large))
          allocate(viswt(Nvisfiles,large))
          allocate(pvisr(Nvisfiles,large))
          allocate(pvisi(Nvisfiles,large))
          allocate(pvisr1(Nvisfiles,large))
          allocate(pvisi1(Nvisfiles,large))
          allocate(LCM(large*Nvisfiles*(2*large*Nvisfiles+1)))
          allocate(SZscale(Nvisfiles))
          allocate(SZscale_BA(Nvisfiles))
          allocate(pb_x(Nvisfiles))
          allocate(pb_y(Nvisfiles))
          allocate(pb_x_offset(Nvisfiles))
          allocate(pb_y_offset(Nvisfiles))
          allocate(ph_x(Nvisfiles))
          allocate(ph_y(Nvisfiles))
          allocate(ph_x_offset(Nvisfiles))
          allocate(ph_y_offset(Nvisfiles))
          allocate(telescope(Nvisfiles))
          allocate(SZeflag(Nvisfiles,large))
          allocate(nu(Nvisfiles))
          allocate(pb_sig(Nvisfiles))
          allocate(lmin(Nvisfiles))
          allocate(lmax(Nvisfiles))
          cell = get_double( 'cell') ! arcsec
          thetamin_Planck = cell/2d0/60d0 ! arcmin
          IncludeCMB = get_integer( 'IncludeCMB')
          covmatfile = get_string( 'covmatfile', .true.)
          map_x = get_double( 'map_x')
          map_y = get_double( 'map_y')
	endif

! Planck setup
        PL = get_integer('PL', 0)
        if (PL.ne.0) then
          PLdatafile = get_string('PLdatafile')
          PLPatch = get_integer('PLPatch')
        endif

! hyperparameters
        hyperparameters = get_integer('hyperparameters', 0)
        ! No. of hyperparameters for multi-dataset analysis:
	if (hyperparameters == 1) then
		if (SZ + PL == 2) then 
			Nhyper = 2
		else
			Nhyper = 1 !can still use hyperparameter even for one dataset, it will act as a reciprocal scaling of the covariance matrix.
		endif
	else
		Nhyper = 0
	endif
        allocate(hyper_PriorType(1,Nhyper))
        allocate(hyper_Prior(1,Nhyper,3)) !could use one-less dimension for these arrays since they aren't dependent on NAtoms, but
        allocate(hyperPars(1,Nhyper))	  !kept like this for consistency with other parameters
        allocate(hyperPars_old(1,Nhyper))  

	write(1,'(I0)') 1 !this almost definitely isn't necessary
        do ipar = 1, Nhyper
          write(ipars,'(I0)') ipar
          hyper_PriorType(1,ipar) = get_integer( 'hyper_PriorType(1,'//trim(ipars)//')')
          call get_doubles( 'hyper_Prior(1,'//trim(ipars), 3, hyper_Prior(1, ipar,:), .false.)
        enddo
        hyperPars=0d0

! Common parameters
        if ((PL.ne.0).or.(SZ.ne.0)) then
          Gas = get_integer( 'Gas', 1)
          Temperature = get_integer( 'Temperature', 1)
          MMass = get_integer( 'MMass')
          MassLim = get_double( 'MassLim')
          MassMin = get_double( 'MassMin')
          MassMax = get_double( 'MassMax')
          znow = get_logical( 'znow')
          mass_function = get_integer( 'mass_function')
        endif

! Weak lensing setup
        GL = get_integer( 'GL', 0)
        Varyzs = get_integer( 'Varyzs', 0)
        if (GL .ne. 0) then
          GLdatafile = get_string( 'GLdatafile', .true.)
          Ngals = get_integer( 'Ngals', 0)
          survey = get_logical( 'survey', .false.)
          arrGal = get_logical( 'arrGal', .false.)
          nxpix = get_integer( 'nxpix', 0)
          nypix = get_integer( 'nypix', 0)
          nxmin = get_integer( 'nxmin', 0)
          nymin = get_integer( 'nymin', 0)
          pix_sizex = get_double( 'pix_sizex')
          pix_sizey = get_double( 'pix_sizey')
          survey_xl = get_double( 'survey_xl', 1.d0)
          survey_yl = get_double( 'survey_yl', 1.d0)
          if (.not. survey) then
            allocate(x(Ngals))
            allocate(y(Ngals))
            allocate(e1(Ngals))
            allocate(e2(Ngals))
            allocate(e1err(Ngals))
            allocate(e2err(Ngals))
            allocate(wt1(Ngals))
            allocate(wt2(Ngals))
            allocate(GLeflag(Ngals))
            allocate(gamma1(Ngals))
            allocate(gamma2(Ngals))
            allocate(kappa(Ngals))
            allocate(g1(Ngals))
            allocate(g2(Ngals))
            allocate(z_s(Ngals))
          else
            allocate(e1s(nxpix,nypix))
            allocate(e2s(nxpix,nypix))
            allocate(g1s(nxpix,nypix))
            allocate(g2s(nxpix,nypix))
            allocate(kappas(nxpix,nypix))
            allocate(gamma1s(nxpix,nypix))
            allocate(gamma2s(nxpix,nypix))
            allocate(e1errs(nxpix,nypix))
            allocate(e2errs(nxpix,nypix))
            allocate(wt1s(nxpix,nypix))
            allocate(wt2s(nxpix,nypix))
            allocate(z_ss(nxpix,nypix))
            allocate(GLeflags(nxpix,nypix))
            allocate(ptInPix(nxpix,nypix,50))
          endif
          if (Varyzs == 1) then
            nzsplanes = get_integer( 'nzsplanes', 0)
            allocate(zs_PriorType(nzsplanes))
            allocate(zs_Prior(nzsplanes,3))
          endif
        endif

! Point source setup
        NSrc = get_integer( 'NSrc')
	if (SZ==0) NSrc=0
        allocate(nu0(0:NSrc))
        nu0(0) = get_double( 'nu0(0)')
        if (NSrc > 0) then
          allocate(Src_PriorType(NSrc,4))
          allocate(Src_Prior(NSrc,4,3))
          allocate(SrcPars(NSrc,4))
          allocate(flux0(NSrc))
          allocate(kernel(NSrc,100,3))
          allocate(nkernel(NSrc))
          ! Try two different methods for getting prior_min and max
          ! If specified individually
          prior_min = get_double( 'prior_min', -10.d0, .true.)
          prior_max = get_double( 'prior_max', -10.d0, .true.)
          ! If specified as a pair
          if (prior_min == -10.d0) then
            call get_doubles( 'prior_min,prior_max', 2, prior_lims)
            prior_min = prior_lims(1)
            prior_max = prior_lims(2)
          endif
          do isrc = 1, NSrc
            write(isrcs,'(I0)') isrc
            nu0(isrc) = get_double( 'nu0('//trim(isrcs)//')')
            do ipar = 1, 4
              write(ipars,'(I0)') ipar
              Src_PriorType(isrc,ipar) = get_integer( 'Src_PriorType('//trim(isrcs)//','//trim(ipars)//')')
              call get_doubles( 'Src_Prior('//trim(isrcs)//','//trim(ipars), 3, Src_Prior(isrc, ipar,:), .false.)
            enddo
          enddo
        endif

! Whether to simulate or analyse real data
        simulate = get_logical( 'simulate')!, .false.)

! No. of parameters per object:
	NPars=NGeoPars+NMassPars*max(GL,MMass)+max(SZ,PL)*(NGasPars*Gas+NTPars*Temperature)+1

! No. of nuisance parameters:
	NNuisance = SZ*4*NSrc + Varyzs

! Total no. of dimensions of parameter space:
	NDim = NAtoms*NPars + NNuisance + Nhyper !only one hyperparameter per dataset, regardless of number of objects

! Whether to sample the prior or analyse data
        SamplePrior = get_logical( 'SamplePrior', .false.)

! Nested sampling parameters
        which_sampler = get_string( 'which_sampler', .true., .true., 'M')
        n_IS = get_logical( 'n_IS', .false.)
        n_mmodal = get_logical( 'n_mmodal', .true.)
        n_ceff = get_logical( 'n_ceff', .false.)
        if (which_sampler == 'P') then
          nest_nlive = get_integer( 'nest_nlive', 100)
          nest_nrep = get_integer( 'nest_nrep', 0)
        else
          nest_nlive = get_integer( 'nest_nlive', 1000)
        endif
        
        n_efr = get_double( 'n_efr', 0.8d0)
        n_rseed = get_integer( 'n_rseed', -1)
        n_tol = get_double( 'n_tol', 0.5d0)
        n_fb = get_logical( 'n_fb', .true.)
        n_maxModes = get_integer( 'n_maxModes', 20)
        n_root = get_string( 'n_root', .true.)

        close(params_unit)

        if (is_formatted.and.is_root) call write_unformatted(out_file)

    end subroutine read_params
!
!=======================================================================
    function get_string(key_word,verb1,exact1,dflt,ith)

        implicit none
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        logical, intent(in), optional :: verb1, exact1 !> verbose output, exact matching
        character(len=*),intent(in),optional  :: dflt  !> keyword to search for
        integer,intent(in), optional :: ith       !> Get the ith instance of this string

        character(len=STR_LENGTH) :: get_string  ! string following keyword

        logical :: verb, exact

        verb = .false.
        if (present(verb1)) verb = verb1
        verb = (verb .and. verbose)

        if (is_formatted) then

          exact = .true.
          if (present(exact1)) exact = exact1

          ! Look for the keyword after the current read-point in the file
          if (present(ith)) then
            get_string = search_text(key_word,exact,ith)
          else
            get_string = search_text(key_word,exact)
          endif

          ! If not found, check from the beginning
          if(trim(get_string)=='') then
            rewind(params_unit)
            if (present(ith)) then
              get_string = search_text(key_word,exact,ith)
            else
              get_string = search_text(key_word,exact)
            endif
          endif

          if(trim(get_string)==''.and. present(dflt)) get_string=trim(dflt)

        else
          read(params_unit) get_string
        endif

        get_string = adjustl(get_string)

        if (verb .and. is_root) write(*,*) key_word//'='//trim(get_string)

    end function get_string
!=======================================================================
    function search_text(key_word,exact,ith)

        implicit none
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        logical, intent(in) :: exact !> exact matching
        integer,intent(in), optional :: ith       !> Get the ith instance of this string

        character(len=STR_LENGTH) :: search_text ! string following keyword

        character(len=STR_LENGTH) :: keyword    !> keyword to search for

        character(len=STR_LENGTH) :: line_buffer, temp_line     ! Line buffer
        integer :: i_equals, i_equals2 ! placement of equals signs
        integer :: end_str  ! account for possible comment after value
        integer :: counter
        logical :: continued, found

        write(keyword,'(A)') key_word

        search_text = ''

        counter=1
        continued=.false.

        io_stat = 0
        do while(io_stat==0) 
            ! Read in the next line
            read(params_unit,'(A)',iostat=io_stat) temp_line

            ! Skip any comment lines
            if( scan(trim(temp_line),comment) == 1 ) cycle

            ! Remove any tabs, quotation marks
            call remove_bad_chars(temp_line)

            ! Remove any in-line comments
            end_str = scan(temp_line,comment)
            if (end_str .ne. 0) then
               temp_line = temp_line(:end_str-1)
            endif

            ! Remove parameter( or data in case of using old-style include files
            temp_line = adjustl(temp_line)
            if (temp_line(1:9)=='parameter') then
                end_str = scan(temp_line, '(')
                temp_line = temp_line(end_str+1:)
            elseif (temp_line(1:4) == 'data') then
                temp_line = adjustl(temp_line(5:))
            endif

            ! Check for breaks over lines
            end_str = scan(temp_line, '&')
            if (end_str .ne. 0) temp_line = temp_line(:end_str-1)
            if (continued) then
              line_buffer = trim(line_buffer) // trim(temp_line)
            else
              line_buffer = temp_line
            endif
            if (end_str .ne. 0) then
              continued = .true.
              cycle
            else
              continued = .false.
            endif

            ! Search for equals signs
            i_equals = scan(line_buffer,equals)
            if(i_equals==0) cycle

            ! check to see if this matches our keyword, two chances to account for, eg map_x,map_y
            !if( trim(adjustl(line_buffer(:i_equals-1))) == trim(adjustl(keyword)) ) then
            if ( index(line_buffer(:i_equals-1), trim(adjustl(keyword))) == 0) then
               ! check for a second '=' match
               i_equals2 = scan(line_buffer(i_equals+1:), equals)
               if (i_equals2.ne.0) then
                 end_str = scan(line_buffer, ',')
                 line_buffer = line_buffer(end_str+1:)
                 i_equals = scan(line_buffer, equals)
               endif
            endif

            ! Now recheck the adjusted line or check the original line
            if (exact) then
               found = ( trim(adjustl(line_buffer(:i_equals-1))) == trim(adjustl(keyword)) )
            else
               found = ( index(line_buffer(:i_equals-1), trim(adjustl(keyword))) .ne. 0)
            endif
            if (found) then
                if(present(ith)) then
                    
                    if(counter==ith) then
                        search_text = adjustl(line_buffer(i_equals+1:))
                        exit
                    else
                        counter=counter+1
                    end if
                else
                    search_text = adjustl(line_buffer(i_equals+1:))
                    exit
                end if

            end if

        end do

        ! Get rid of a trailing slash or bracket if present
        if (trim(search_text)/='') then
            end_str = len_trim(search_text)
            if (search_text(end_str:end_str) == '/') search_text = search_text(:end_str-1)
            if (search_text(end_str:end_str) == ')') search_text = search_text(:end_str-1)
        endif

        search_text = adjustl(search_text)

    end function search_text
!=======================================================================
    subroutine get_strings(key_word,num,strings)
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        integer, intent(in)          :: num       !> number of integers
        character(len=STR_LENGTH), intent(out), dimension(num) :: strings

        character(len=STR_LENGTH) :: string1, string  ! string following keyword
        integer :: i, end_str

        ! Get the string
        if (verbose .and. is_root) write(*,*) key_word//':'
        if (is_formatted) then
          string = get_string(key_word)
          string1 = string
          if (trim(string)=='') call halt_program('ini error: no keyword '//trim(key_word))

          i=1
          do while( trim(string)/='' )
            if (i .gt. num) then
                if(is_root) write(*,*) 'Warning: more than', num, 'elements found in', trim(string1)
              exit
            endif
            end_str = scan(string, ',')
            if (end_str == 0) end_str = len_trim(string)+1
            strings(i)=string(1:end_str-1)
            if (verbose .and. is_root) write(*,*) i, trim(strings(i))
            call next_element(string,',')
            i=i+1
          end do
        else
          do i = 1, num
            read(params_unit) strings(i)
            if (verbose .and. is_root) write(*,*) i, trim(strings(i))
          enddo
        endif

    end subroutine get_strings
!=======================================================================
    subroutine remove_bad_chars(line_buffer)
        implicit none
        character(len=STR_LENGTH),intent(inout) :: line_buffer
        character(len=STR_LENGTH)               :: temp_line
        integer                                 :: i
        character(len=3)                        :: bad_chars

        bad_chars = achar(9)//'''"'

        temp_line=''
        do i = 1, len(line_buffer)
          if (scan(line_buffer(i:i), bad_chars) == 0) then
            temp_line=trim(temp_line)//line_buffer(i:i)
          endif
        enddo
        line_buffer = temp_line
    end subroutine remove_bad_chars
!=======================================================================
    subroutine next_element(line_buffer,delimiter) 
        implicit none
        character(len=STR_LENGTH),intent(inout)  :: line_buffer ! Line buffer
        character :: delimiter
        integer   :: del_pos

        del_pos = scan(line_buffer,delimiter)
        if (del_pos > 0) then
          line_buffer = trim(line_buffer(del_pos+1:)) ! Find the next element
        else
          line_buffer = ''
        endif

    end subroutine next_element
!=======================================================================
    function get_integer(key_word,dflt,exact1)

        character(len=*),intent(in)  :: key_word  !> keyword to search for
        integer,intent(in),optional :: dflt
        logical,intent(in),optional :: exact1
        
        character(len=STR_LENGTH) :: string  ! string following keyword
        integer :: get_integer  ! integer following keyword
        logical :: exact

        if(present(dflt)) get_integer=dflt

        if (is_formatted) then
          ! True by default
          exact = .true.
          if (present(exact1)) exact = exact1

          string = get_string(key_word,.false.,exact)
          if(trim(string)/='') then
            read(string,*) get_integer
          else if(present(dflt)) then
            get_integer = dflt
          else
            call halt_program('ini error: no keyword '//trim(key_word))
          end if
        else
          read(params_unit) get_integer
        endif
        if (verbose.and.is_root) write(*,*) key_word, '=', get_integer

    end function get_integer
!=======================================================================
    subroutine get_integers(key_word,num,ints,exact1)
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        integer, intent(in)          :: num       !> number of integers
        integer, intent(out), dimension(num) :: ints
        logical, intent(in), optional :: exact1 !> look for exact match if true

        character(len=STR_LENGTH) :: string1, string  ! string following keyword

        integer :: i
        logical :: exact

        if (verbose.and.is_root) write(*,*) key_word//':'
        if (is_formatted) then
          exact = .true.
          if (present(exact1)) exact = exact1

          ! Get the string
          string = get_string(key_word,.false.,exact)
          string1 = string
          if (trim(string)=='') call halt_program('ini error: no keyword '//trim(key_word))

          ! Initialise array
          ints = 0

          i=1
          do while( trim(string)/='' )
            if (i .gt. num) then
              write(*,*) 'Warning: more than', num, 'elements found in', trim(string1)
              exit
            endif
            read(string,*) ints(i)
            if (verbose.and.is_root) write(*,*) i, ints(i)
            call next_element(string,',')
            i=i+1
          end do
        else
          do i = 1, num
            read(params_unit) ints(i)
            if (verbose.and.is_root) write(*,*) i, ints(i)
          enddo
        endif

    end subroutine get_integers
!=======================================================================
    function get_double(key_word,dflt,exact1)
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        double precision,intent(in),optional :: dflt
        logical, intent(in), optional :: exact1 !> look for exact match if true

        character(len=STR_LENGTH) :: string  ! string following keyword
        double precision :: get_double  ! double following keyword
        logical :: exact

        if (is_formatted) then
          exact = .true.
          if (present(exact1)) exact = exact1

          string = get_string(key_word,.false.,exact)
          if(trim(string)/='') then
            read(string,*) get_double
          else if(present(dflt)) then
            get_double = dflt
          else
            call halt_program('ini error: no keyword '//trim(key_word))
          end if
        else
          read(params_unit) get_double
        endif
        if (verbose.and.is_root) write(*,*) key_word, '=', get_double

    end function get_double
!=======================================================================
    subroutine get_doubles(key_word,num,doubles,exact1)
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        integer, intent(in)          :: num       !> number of doubles
        double precision, intent(out), dimension(num) :: doubles
        logical, intent(in), optional :: exact1 !> look for exact match if true

        character(len=STR_LENGTH) :: string1, string  ! string following keyword

        integer :: i
        logical :: exact

        if (verbose.and.is_root) write(*,*) key_word//':'
        ! Initialise doubles
        doubles = 0.d0
        if (is_formatted) then
          exact = .true.
          if (present(exact1)) exact = exact1

          ! Get the string
          string = get_string(key_word,.false.,exact)
          string1 = string
          if (trim(string)=='') call halt_program('ini error: no keyword '//trim(key_word))

          i=1
          do while( trim(string)/='' )
            if (i .gt. num) then
                if(is_root) write(*,*) 'Warning: more than ', num, ' elements found in', trim(string1)
              exit
            endif
            read(string,*) doubles(i)
            if (verbose.and.is_root) write(*,*) i, doubles(i)
            call next_element(string,',')
            i=i+1
          end do
        else
          do i = 1, num
            read(params_unit) doubles(i)
            if (verbose.and.is_root) write(*,*) i, doubles(i)
          enddo
        endif

    end subroutine get_doubles
!=======================================================================
    function get_logical(key_word,dflt)
        
        character(len=*),intent(in)  :: key_word  !> keyword to search for
        logical,intent(in),optional :: dflt

        character(len=STR_LENGTH) :: string  ! string following keyword
        logical :: get_logical  ! logical following keyword

        if (is_formatted) then
          string = get_string(key_word)
          if(trim(string)/='') then
            read(string,*) get_logical
          else if(present(dflt)) then
            get_logical = dflt
          else
            call halt_program('ini error: no keyword '//trim(key_word))
          end if
        else
          read(params_unit) get_logical
        endif
        if (verbose.and.is_root) write(*,*) key_word, '=', get_logical

    end function get_logical
!=======================================================================
    subroutine write_unformatted(out_file)

     implicit none
     character(len=*),intent(in) :: out_file
     integer iatom, ipar, iprior, isrc
     character(len=STR_LENGTH)   :: temp_str

     ! Write to unformatted file
     open(unit=params_unit, form='unformatted', file=out_file, status='replace')

     ! Note that all strings need to be the same length, otherwise it crashes on reading
     write(params_unit) verbose
     write(params_unit) ModelClass
     write(params_unit) Atoms
     write(params_unit) NAtoms
! Cluster geometry
     write(params_unit) GeoModel
     do iatom = 1, NAtoms
        do ipar = 1, NGeoPars
           write(params_unit) Geo_PriorType(iatom,ipar)
           do iprior = 1, 3
             write(params_unit) Geo_Prior(iatom, ipar, iprior)
           enddo
        enddo
     enddo
! Cluster mass
     write(params_unit) MassModel
     write(params_unit) NMassPars
     do iatom = 1, NAtoms
       write(params_unit) NFWstyle(iatom)
       do ipar = 1, NMassPars
          write(params_unit) Mass_PriorType(iatom,ipar)
          do iprior = 1, 3
            write(params_unit) Mass_Prior(iatom, ipar, iprior)
          enddo
       enddo
     enddo
! Cluster gas
     write(params_unit) GasModel
     write(params_unit) NGasPars
     do iatom = 1, NAtoms
        write(params_unit) BetaStyle(iatom)
        write(params_unit) TStyle(iatom)
        do ipar = 1, NGasPars
           write(params_unit) Gas_PriorType(iatom,ipar)
           do iprior = 1, 3
              write(params_unit) Gas_Prior(iatom, ipar, iprior)
           enddo
        enddo
     enddo
! Cluster temperature
     write(params_unit) TModel
     write(params_unit) NTPars
     if (NTPars > 0) then
        do iatom = 1, NAtoms
           do ipar = 1, NTPars
              write(params_unit) T_PriorType(iatom,ipar)
              do iprior = 1, 3
                write(params_unit) T_Prior(iatom, ipar, iprior)
              enddo
           enddo
        enddo
     endif
! Cluster redshift
     do iatom = 1, NAtoms
        write(params_unit) z_PriorType(iatom)
        do iprior = 1, 3
          write(params_unit) z_Prior(iatom, iprior)
        enddo
     enddo
! Derived parameters
     write(params_unit) aux_dim 
! SZ setup
     write(params_unit) SZ
     if ((SZ == 1) .or. (PL == 1)) then
       write(params_unit) Gas
       write(params_unit) Temperature
       write(params_unit) MMass
       write(params_unit) Nvisfiles
       do ifile = 1, Nvisfiles
         write(params_unit) visdatafile(ifile)
       enddo
       do ifile = 1, Nvisfiles
         write(params_unit) Nvis(ifile)
       enddo
       write(params_unit) large
       write(params_unit) cell
       write(params_unit) IncludeCMB
       temp_str = covmatfile
       write(params_unit) temp_str
       write(params_unit) map_x
       write(params_unit) map_y
       write(params_unit) MassLim
       write(params_unit) MassMin
       write(params_unit) MassMax
       write(params_unit) znow
       write(params_unit) mass_function
     endif

! Planck setup
     write(params_unit) PL
     if (PL.ne.0) then
       temp_str = PLdatafile
       write(params_unit) temp_str
       write(params_unit) PLPatch
     endif

! hyperparameter setup
	write(params_unit) hyperparameters
	if (hyperparameters == 1) then
     		write(params_unit) Nhyper
        	do ipar = 1, Nhyper
           		write(params_unit) hyper_PriorType(1,ipar)
           		do iprior = 1, 3
              			write(params_unit) hyper_Prior(1, ipar, iprior)
           		enddo
        	enddo
	endif

! Weak lensing setup
     write(params_unit) GL
     write(params_unit) Varyzs
     if (GL .ne. 0) then
          write(params_unit) GLdatafile
          write(params_unit) Ngals
          write(params_unit) survey
          write(params_unit) arrGal
          write(params_unit) nxpix
          write(params_unit) nypix
          write(params_unit) nxmin
          write(params_unit) nymin
          write(params_unit) pix_sizex
          write(params_unit) pix_sizey
          write(params_unit) survey_xl
          write(params_unit) survey_yl
          if (Varyzs == 1) then
            write(params_unit) nzsplanes
          endif
     endif

! Point source setup
     write(params_unit) NSrc
     write(params_unit) nu0(0)
     if (NSrc > 0) then
          write(params_unit) prior_min
          write(params_unit) prior_max
          do isrc = 1, NSrc
            write(params_unit) nu0(isrc)
            do ipar = 1, 4
              write(params_unit) Src_PriorType(isrc,ipar)
              do iprior = 1, 3
                 write(params_unit) Src_Prior(isrc, ipar, iprior)
              enddo
            enddo
          enddo
     endif

! Whether to simulate or analyse real data
     write(params_unit) simulate
! Whether to sample the prior or analyse data
     write(params_unit) SamplePrior 
! Nested sampling parameters
     temp_str = which_sampler
     write(params_unit) temp_str
     write(params_unit) n_IS
     write(params_unit) n_mmodal
     write(params_unit) n_ceff
     write(params_unit) nest_nlive
     if (which_sampler == 'P') then
         write(params_unit) nest_nrep
     endif
     write(params_unit) n_efr
     write(params_unit) n_rseed
     write(params_unit) n_tol
     write(params_unit) n_fb
     write(params_unit) n_maxModes
     temp_str = n_root
     write(params_unit) temp_str
     close(unit=params_unit)

    end subroutine write_unformatted
!=======================================================================
    subroutine write_paramnames
      implicit none

      integer ival, i, j, ii, ival1
      logical nuis
      character(len=3) :: js, iis
      character(len=4) :: ext
      character(len=50) :: parname,fmt,parunit !made longer as it was truncating Y_tot name. KJ 05/06/17

      if (GL == 1) then
        write(*,*) 'Do not know how to assign parameter names for gravitational lensing parameters, not writing', trim(n_root)//'.paramnames'
        return
      endif

      open(unit=params_unit, form='formatted', file=trim(n_root)//'.paramnames', status='replace')
      ival1 = 0
      do ival = 1, Ndim
         nuis = .false.
         parunit=''
         if (ival>NAtoms*NPars + Nhyper) nuis=.true. !don't treat hyperparameters as nuissance parameters for now. KJ 06/06/17
         if (.not. nuis) then
	    if (ival<= NAtoms*NPars) then
            	j=(ival-1)/NPars+1
            	write(js,'(I0)') j
            	i=ival-(j-1)*NPars
	    else
		i=ival
	    endif !parameters are ordered Atom1par1,...,Atom1parn,Atomnpar1,...,Atomnparn,hyperparameters,nuissance
            if(i <= NGeoPars) then
               ii = i
               if (Geo_PriorType(j,ii)<=0) cycle
               if (ii == 1) then
                  parname='x_0'
                  parunit='\mathrm{arcsec}'
               elseif (ii == 2) then
                  parname='y_0'
                  parunit='\mathrm{arcsec}'
               elseif (ii == 3) then
                  parname='\theta'
                  parunit='\mathrm{deg}'
               elseif (ii == 4) then
                  parname='f'
               endif
            elseif (i <= NGeoPars+1 ) then
               if (z_PriorType(j)<=0) cycle
               parname='z'
            else if( Mass == 1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass ) ) then
               ii = i - NGeoPars - 1
               if (Mass_PriorType(j,ii)<=0) cycle
               ! Not used for newer models?  Not sure what these are
               write(iis,'(I0)') ii
               parname='MassPars('//trim(js)//','//trim(iis)//')'
            else if( ((SZ == 1) .or. (PL == 1)) .and. Gas == 1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass + NGasPars * Gas ) ) then
               ii = i - NGeoPars - 1 - NMassPars * Mass
               if (Gas_PriorType(j,ii)<=0) cycle
               if (GasModel==3 .or. GasModel==5 .or. GasModel==6) then !adapted for GM=6 kj 01/02/17
                   ! GNFW-Planck model
                   if (GasModel==3) then
                      if (ii==1) then
                         parname='\theta_{s}'
                         parunit='\mathrm{arcmin}'
                      elseif (ii==2) then
                         parname='Y_{\mathrm{tot}}'
                         parunit='\mathrm{arcmin}^2'
                      endif
                   ! DM-GNFW model & Ein_DM-GNFW model
                   else if(GasModel==5 .or. GasModel==6) then !adapted for GM=6 kj 01/02/17
                      if (ii==1) then
                         parname='M_{T,200}'
                         parunit='M_{\odot}'
                      elseif (ii==2) then
                         parname='f_{\mathrm{gas},200}'
                      endif
		      if(GasModel==6) then
                         if (ii==7) then
                            parname='\alpha_{Ein}'
                         endif
                      endif
                   endif
                   ! Common parameters for these three
                   if (ii==3) then
                       parname='\gamma'
                   elseif (ii==4) then
                       parname='\alpha'
                   elseif (ii==5) then
                       parname='\beta'
                   elseif (ii==6) then
                       parname='c_{500}'
                   endif
               elseif (GasModel==1 .and. BetaStyle(1)==4) then
                   if (ii==1) then
                     parname='r_{c}'
                     parunit='\mathrm{kpc}'
                   elseif (ii==2) then
                     parname='\beta'
                   elseif (ii==3) then
                     parname='M_{T,200}'
                     parunit='M_{\odot}'
                   elseif (ii==4) then
                     parname='f_{\mathrm{gas},200}'
                   endif
               else
                   ! These models aren't used much anymore so haven't bothered spelling out the parameters
                   write(iis,'(I0)') ii
                   parname='GasPars('//trim(js)//','//trim(iis)//')'
               endif
            else if( ((SZ == 1) .or. (PL == 1)) .and. Gas == 1 .and. Temperature==1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass + ( NGasPars + NTPars ) * Gas ) ) then
                ii = i - NGeoPars - 1 - NMassPars * Mass - NGasPars * Gas
                if (T_PriorType(j,ii)<=0) cycle
                write(iis,'(I0)') ii
                parname='TPars('//trim(js)//','//trim(iis)//')'
	   !hyperparameters
            else if ((i > NPars * NAtoms) .and. hyperparameters == 1) then
	       j = 1
		write(js,'(I0)') j
               ii = i - NPars * NAtoms
               if (hyper_PriorType(j,ii)<=0) cycle
               if (ii == 1) then
                  parname='\alpha_{\mathrm{Lhood 1}}'
                  parunit=''
	       endif
	       if (ii == 2) then
                  parname='\alpha_{\mathrm{Lhood 2}}'
                  parunit=''
	       endif
	    endif
            if (NAtoms > 1) then
              parname=parname//','//trim(js)
            endif
         ! Nuisance parameters
         else
            i=ival-(NAtoms*NPars + Nhyper)
            j=(i-1)/4+1
            if(SZ==1 .and. SourceSubtract==1 .and. i.le.4*NSrc) then	
                ii=i-(j-1)*4
                if (Src_PriorType(j,ii)<=0) cycle
                write(js,'(I0)') j
                if (ii==1) then
                  parname='RA_{'//trim(js)//'}'
                  parunit='\mathrm{deg}'
                elseif (ii==2) then
                  parname='\delta_{'//trim(js)//'}'
                  parunit='\mathrm{deg}'
                elseif (ii==3) then
                  parname='S_{'//trim(js)//'}'
                  parunit='\mathrm{Jy}'
                elseif (ii==4) then
                  parname='\alpha_{'//trim(js)//'}'
                endif
            endif
         endif
         ival1 = ival1 + 1
         fmt='(2X,A1,I3.3,6X,A)'
         if (len(trim(parunit)).gt.0) then
           parname=trim(parname)//' / '//trim(parunit)
         endif
         write(params_unit,fmt) 'p', ival1, trim(parname)
      enddo

      ! Auxiliary parameters - name ending in '*' means getdist.py interprets them as derived
      ext=''
      do j = 1, NAtoms
        write(js,'(I0)') j
        if (NAtoms > 1) then
          ext=','//trim(js)
        endif
        fmt='(2X,A1,I3.3,A1,5X,A)'
        if(GasModel==3)then
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+1, '*', 'y_{\mathrm{coeff}}'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+2, '*', 'y_{0}'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+3, '*', '5.0\times\theta_{500}'//trim(ext)//' / \mathrm{arcmin}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+4, '*', 'Y_{\mathrm{tot}}'//trim(ext)//' / \mathrm{arcmin}^2'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+5, '*', '\theta_{500}'//trim(ext)//' / \mathrm{arcmin}'		
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+6, '*', 'Y_{500}'//trim(ext)//' / \mathrm{arcmin}^2'
        elseif(GasModel==5 .or. GasModel==6)then	
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+1, '*', 'D'//trim(ext)//' / \mathrm{Mpc}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+2, '*', 'c_{200}'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+3, '*', 'GasPars('//trim(js)//',1)'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+4, '*', 'r_{200}'//trim(ext)//' / \mathrm{Mpc}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+5, '*', 'M_{g,200}'//trim(ext)//' / M_{\odot}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+6, '*', 'T_{g,200}'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+7, '*', 'M_{T,500}'//trim(ext)//' / M_{\odot}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+8, '*', 'f_{g,500}'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+9, '*', 'r_{500}'//trim(ext)//' / \mathrm{Mpc}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+10, '*', 'M_{g,500}'//trim(ext)//' / M_{\odot}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+11, '*', 'T_{g,500}'//trim(ext)	 	 	 		
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+12, '*', 'y_{0}'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+13, '*', 'Y_{\mathrm{sph},200}'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+14, '*', 'Y_{\mathrm{sph},500}'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+15, '*', '\theta_{200}'//trim(ext)//' / \mathrm{arcmin}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+16, '*', 'Y_{\mathrm{sph},200} / \mathrm{arcmin}^2'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+17, '*', '\theta_{500}'//trim(ext)//' / \mathrm{arcmin}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+18, '*', 'Y_{\mathrm{sph},500} / \mathrm{arcmin}^2'//trim(ext)
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+19, '*', 'P_{ei}'//trim(ext)
	   if (GasModel==5) then
              write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+20, '*', '\rho_{s}'//trim(ext)
              write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+21, '*', 'r_{s}'//trim(ext)//' / \mathrm{Mpc}'
	   else
	      write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+20, '*', '\rho_{-2}'//trim(ext)
              write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+21, '*', 'r_{-2}'//trim(ext)//' / \mathrm{Mpc}'
           endif
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+22, '*', 'r_{p}'//trim(ext)//' / \mathrm{Mpc}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+23, '*', '\theta_{p}'//trim(ext)//' / \mathrm{arcmin}'
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+24, '*', '\rho_{\mathrm{crit},z}'//trim(ext)	
           write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+25, '*', 'Y_{\mathrm{cyl}} / \mathrm{arcmin}^2'//trim(ext)
	   write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+26, '*', 'Y_{\mathrm{tot}} / \mathrm{arcmin}^2'//trim(ext)
           ! Extra parameters added to the end of the cube at the end of FUserbuild, MultiNEST only
           if (which_sampler=='M') then
             write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+27, '*', 'y_{\mathrm{cent,map}} / \mathrm{arcmin}^2'//trim(ext)
             write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+28, '*', '\mathrm{central\, cluster\, flux} / \mathrm{Jy}'//trim(ext)
             write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+29, '*', '\mathrm{central\, cluster\, temp} / \mathrm{K}'//trim(ext)
           endif
        else
           do ival = 1, aux_dim
              write(iis,'(I0)') ival
              write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+ival, '*', 'aux('//trim(js)//','//trim(iis)//')'
           enddo
           ! Extra parameters added to the end of the cube at the end of FUserbuild, MultiNEST only
           if (which_sampler=='M' .and. GasModel==1) then
             write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+ival, '*', 'y_{\mathrm{cent,map}} / \mathrm{arcmin}^2'//trim(ext)
             write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+ival+1, '*', '\mathrm{central\, cluster\, flux} / \mathrm{Jy}'//trim(ext)
             write(params_unit,fmt) 'p',edim+(j-1)*aux_dim+ival+2, '*', '\mathrm{central\, cluster\, temp} / \mathrm{K}'//trim(ext)
           endif
        endif
      enddo  

      close(params_unit)
    end subroutine write_paramnames
!=======================================================================
    subroutine write_ranges
      implicit none

      integer ival, i, j, ii, ival1, prior_type
      double precision  :: pmin, pmax
      character(len=30) :: pmin_s, pmax_s
      logical nuis, theta
      character(len=30) :: fmt

      if (GL == 1) then
        write(*,*) 'Do not know how to assign parameter names for gravitational lensing parameters, not writing', trim(n_root)//'.ranges'
        return
      endif

      open(unit=params_unit, form='formatted', file=trim(n_root)//'.ranges', status='replace')
      ival1 = 0
      do ival = 1, Ndim

         nuis = .false.
         pmin_s = 'N'
         pmax_s = 'N'
         theta=.false.
         if (ival>NAtoms*NPars + Nhyper) nuis=.true.
         if (.not. nuis) then
             if (ival<= NAtoms*NPars) then
            	j=(ival-1)/NPars+1
            	i=ival-(j-1)*NPars
	    else
		i=ival
	    endif
            if(i <= NGeoPars) then
               ii = i
               if (Geo_PriorType(j,ii)<=0) cycle
               prior_type=Geo_PriorType(j,ii)
               if (prior_type.eq.9) then
                 pmin=Geo_Tri_Prior(j,ii,1)
                 pmax=Geo_Tri_Prior(j,ii,2)
               elseif (prior_type.eq.12) then
                 pmin=Geo_Prior(j,ii,1)
                 pmax=Geo_Prior(j,ii,3)
               else
                 pmin=Geo_Prior(j,ii,1)
                 pmax=Geo_Prior(j,ii,2)
               endif
            elseif (i <= NGeoPars+1 ) then
               if (z_PriorType(j)<=0) cycle
               prior_type=z_PriorType(j)
               if (prior_type.eq.12) then
                 pmin=z_Prior(j,1)
                 pmax=z_Prior(j,3)
               else
                 pmin=z_Prior(j,1)
                 pmax=z_Prior(j,2)
               endif
            else if( Mass == 1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass ) ) then
               ii = i - NGeoPars - 1
               if (Mass_PriorType(j,ii)<=0) cycle
               ! Not used for newer models?  Not sure what these are
               prior_type=Mass_PriorType(j,ii)
               if (prior_type.eq.12) then
                 pmin=Mass_Prior(j,ii,1)
                 pmax=Mass_Prior(j,ii,3)
               else
                 pmin=Mass_Prior(j,ii,1)
                 pmax=Mass_Prior(j,ii,2)
               endif
            else if( ((SZ == 1) .or. (PL == 1)) .and. Gas == 1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass + NGasPars * Gas ) ) then
               ii = i - NGeoPars - 1 - NMassPars * Mass
               if (Gas_PriorType(j,ii)<=0) cycle
               prior_type=Gas_PriorType(j,ii)
               if (prior_type.eq.12) then
                 pmin=Gas_Prior(j,ii,1)
                 pmax=Gas_Prior(j,ii,3)
               else
                 pmin=Gas_Prior(j,ii,1)
                 pmax=Gas_Prior(j,ii,2)
               endif
               if (GasModel==3.and.ii==1) then
                 ! GNFW-Planck model, hard cut on theta_s
                 theta=.true.
               endif
            else if( ((SZ == 1) .or. (PL == 1)) .and. Gas == 1 .and. Temperature==1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass + ( NGasPars + NTPars ) * Gas ) ) then
                ii = i - NGeoPars - 1 - NMassPars * Mass - NGasPars * Gas
                if (T_PriorType(j,ii)<=0) cycle
                prior_type=T_PriorType(j,ii)
                if (prior_type.eq.12) then
                 pmin=T_Prior(j,ii,1)
                 pmax=T_Prior(j,ii,3)
                else
                 pmin=T_Prior(j,ii,1)
                 pmax=T_Prior(j,ii,2)
                endif
            else if ((i > NPars * NAtoms) .and. hyperparameters == 1) then !hyperparameters
	       j = 1
               ii = i - NPars * NAtoms
               if (hyper_PriorType(j,ii)<=0) cycle
               prior_type=hyper_PriorType(j,ii)
               if (prior_type.eq.12) then !probably not necessary but here for consistency
                 pmin=hyper_Prior(j,ii,1)
                 pmax=hyper_Prior(j,ii,3)
               else
                 pmin=hyper_Prior(j,ii,1)
                 pmax=hyper_Prior(j,ii,2)
               endif
	    endif
         ! Nuisance parameters
         else
            i=ival-(NAtoms*NPars + Nhyper)
            j=(i-1)/4+1
            if(SZ==1 .and. SourceSubtract==1 .and. i.le.4*NSrc) then	
                ii=i-(j-1)*4
                if (Src_PriorType(j,ii)<=0) cycle
                ! For spectral indices only, prior type 12 is 10C prior convolved with Gaussian so no limits
                if ((ii.eq.4).and.(Src_PriorType(j,ii).eq.12)) then
                  ival1 = ival1 + 1
                  cycle
                endif
                prior_type=Src_PriorType(j,ii)
                if (prior_type.eq.12) then
                 pmin=Src_Prior(j,ii,1)
                 pmax=Src_Prior(j,ii,3)
                else
                 pmin=Src_Prior(j,ii,1)
                 pmax=Src_Prior(j,ii,2)
                endif
            endif
         endif
         ival1 = ival1 + 1
         select case (prior_type)
           case (3, 4, 6, 13, 14)
             ! No hard limits in the prior, no need to write anything out
             cycle
           case (7)
             ! Spectral index prior
             pmin=prior_min
             pmax=prior_max
             write(pmin_s,'(ES)') pmin
             write(pmax_s,'(ES)') pmax
           case (12)
             ! Truncated Gaussian, third parameter is either upper or lower limit depending on relation to mean
             if (pmax.gt.pmin) then
               write(pmax_s,'(ES)') pmax
             else
               write(pmin_s,'(ES)') pmax
               if ((theta).and.(pmax.lt.1.3)) then
                 write(pmin_s,'(ES)') 1.3d0
               endif
             endif
           case default
             ! All other priors have 1st parameter = min, 2nd = max
             if ((theta).and.(pmin.lt.1.3)) then
               pmin=1.3d0
             endif
             write(pmin_s,'(ES)') pmin
             write(pmax_s,'(ES)') pmax
         end select
         fmt='(2X,A1,I3.3,4X,A,2X,A)'
         write(params_unit,fmt) 'p', ival1, trim(pmin_s), trim(pmax_s)
      enddo

      close(params_unit)
    end subroutine write_ranges

!=======================================================================
end module ReadWrite
