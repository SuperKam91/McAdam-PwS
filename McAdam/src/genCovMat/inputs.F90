module inputs

	use globals
	use confusion_noise
	use utilities
        use telescope1
        use params
        use fits_handling

contains

!-----*-----------------------------------------------------------------

subroutine read_params(diagflag,dsflag)
	implicit none
	integer   diagflag,dsflag
	character*500 fileRoot

	! option to write likelihood function to screen and a file

	write(*,'(a,$)') ' Root for the files: '
	read(*,'(a)') fileRoot
	
	integralsfn = trim(fileRoot)//'Integrals.dat'
	outputfn = trim(fileRoot)//'Choles.lcm'

	! choose diagonal or full covariance matrix

	write(*,'(a,$)')' Use full (1) or diagonal (0) covariance matrix: '
	read(*,*) diagflag

	! choose sparse or dense matrix routines

	write(*,'(a,$)') ' Use sparse (1) or dense (0) matrix routines: '
	read(*,*) dsflag
	
end subroutine read_params

!-----*-----------------------------------------------------------------

subroutine fread_params(diagflag,dsflag, nvmax)
	implicit none
	integer   diagflag,dsflag, nvmax
	character*500 fileRoot
   
	! choose diagonal(diagflag=0) or full (diagflag=1) covariance matrix

	diagflag=1

	! choose sparse (dsflag=1)  or dense (dsflag=0) matrix routines
	
        IF(nvmax.LE.1000)THEN
	
	dsflag=0      ! dense matrix approach
	
       ELSE	

        dsflag=1      ! sparse matrix approach
	
      ENDIF	   
   
	
	write(*,'(a,$)') ' Root for the output files: '
	read(input_unit,'(a)') fileRoot
	write(*,*) trim(fileRoot)
	
	integralsfn = trim(fileRoot)//'Integrals.dat'
	write(*,'(a,$)') ' Output filename for the integrals: '
	write(*,*) trim(integralsfn)
		
     IF(dsflag.EQ. 0)THEN 
     
         WRITE(*,*)'dense matrix approach' 	
	outputfn = trim(fileRoot)//'Choles.lcm'
	write(*,'(a,$)') ' Output filename for the covariance matrix: '
	write(*,*) trim(outputfn)
		
     ELSEIF(dsflag.EQ.1)THEN
     
         WRITE(*,*)'sparse matrix approach' 	
        outputfn2=trim(fileRoot)//'sacov.txt'
	write(*,'(a,$)') ' Output filename  for the incomplete Cholesky decomposition preconditioner: '	
	write(*,*) trim(outputfn2)
	
    ENDIF		
    
end subroutine fread_params

!-----*-----------------------------------------------------------------

subroutine read_lim(cpsmin,cpsmax,cpsstep,ncps,gpsmin,gpsmax,gpsstep,ngps,spmin,spmax,nsp,spstep,galyn)
	implicit none
	integer  ncps,ngps,nsp
	double precision   cpsmin,cpsmax,cpsstep
	double precision   gpsmin,gpsmax,gpsstep
	double precision   spmin,spmax,spstep
	character galyn*1

	! ... input cmb spectrum limits in each annulus and no. of steps
	write(*,*)
	write(*,'(a,$)') ' Input range of l^2 C_l/2pi for CMB: '
	read(*,*) cpsmin,cpsmax
	write(*,'(a,$)') ' Input number of CMB PS steps: '
	read(*,*) ncps
	cpsstep=(cpsmax-cpsmin)/float(ncps-1)

	! ... input galaxy spectrum limits in each annulus and no. of steps
	write(*,'(a,$)') ' Estimate PS of Galaxy? (y/n) [n]: '
	read(*,'(a1)') galyn
	if (galyn.eq.'y' .or. galyn.eq.'Y') then
		write(*,'(a,$)') ' Input range of l^2 C_l/2pi for Galaxy: '
	  	read(*,*) gpsmin,gpsmax
	  	write(*,'(a,$)') ' Input number of GAL PS steps: '
	  	read(*,*) ngps
	  	gpsstep=(gpsmax-gpsmin)/float(ngps-1)
	  	write(*,'(a,$)') ' Input limits on GAL flux spectral index: '
	  	read(*,*) spmin,spmax
	  	write(*,'(a,$)') ' Input number of steps for spectral index: '
	  	read(*,*) nsp
	  	spstep=(spmax-spmin)/float(nsp-1)        
	else
	  	nsp=1  
	  	ngps=1
	  	gpsmin=0.
	  	gpsstep=0.
	  	spmin=0.
	  	spstep=0.
	endif
	
end subroutine read_lim

!-----*-----------------------------------------------------------------

subroutine read_inps(nfmax,nf,nskip,freq,freq0,sigma,antwdth,antmax,scaled,scalen,datafn,noisefn,iaddn)
	implicit none
	integer  nfmax,nf,nskip,iaddn
	double precision   freq0,antmax,scaled,scalen
	double precision   freq(nfmax),sigma(nfmax),antwdth(nfmax)
	character*500 datafn(nfmax),noisefn(nfmax)

	integer  i
	double precision   PI,dtor
	
! initialise variables

	PI=4d0*datan(1d0)        ! this really does equal PI
	dtor=PI/180d0            ! degrees/radians conversion

! read inputs

	write(*,*)
	write(*,*) ' Telescope primary beam modelled as', &
		     ' Gaussian at all observing frequencies'
	write(*,*)

	write(*,'(a,$)')' Input number of observing frequencies (4 max.): '
	read(*,*) nf

	write(*,'(a,$)')' Input observing frequencies (GHz): '
	read(*,*) (freq(i),i=1,nf)

	write(*,'(a,$)')' Input FWHM of primary beam at each frequency (degrees): '
	read(*,*) (sigma(i),i=1,nf) 

	write(*,'(a,$)')' Input antenna diameter at each frequency (wavelengths): '
	read(*,*) (antwdth(i),i=1,nf)

	write(*,'(a,$)')' Input frequency for calculation (GHz): '
	read(*,*) freq0

! calculate sigma of primary and 1/e diameter of antenna function

	write(*,*)
	antmax=0d0
	do i=1,nf
	  	sigma(i) = sigma(i)*dtor
	  	sigma(i) = sigma(i)/(2*sqrt(2d0*dlog(2d0)))
		!antwdth(i) = 2d0/(PI*sigma(i)*dsqrt(2d0))
		!write(*,'(a,a,f6.2,a,f6.2,a)') '  1/e diameter of aperture ','function at ',freq(i),' GHz = ',antwdth(i),' wavelengths'
	  	if (antwdth(i).gt.antmax) antmax=antwdth(i)
	enddo
 
! ask for data file names and open

	write(*,*)
	write(*,*) 'Input visibilities data file for each frequency: '
	do i=1,nf
	  	read(*,150) datafn(i) 
	enddo
	write(*,*) 'Input visibilities noise file for each frequency: '
	do i=1,nf
	  	read(*,150) noisefn(i) 
	enddo
 150  	format(a)

	write(*,'(a,$)') ' Input factor for scaling data : '
	read(*,*) scaled

	write(*,'(a,$)') ' Input factor for scaling noise: '
	read(*,*) scalen

	write(*,'(a,$)') ' Add noise to data (1) or not (0): '
	read(*,*) iaddn

	write(*,'(a,$)') ' Input data points skip (0 = no skip): '
	read(*,*) nskip

	return
end subroutine read_inps

!-----*-----------------------------------------------------------------

subroutine fread_inps(nfmax,nf,nskip,freq,freq0,sigma,antwdth, &
		f_ra,f_dec,ra_ref,dec_ref,antmax,scaled,scalen,datafn,noisefn,iaddn)
	implicit none
	integer  nfmax,nf,nskip,iaddn
	double precision   freq0,antmax,scaled,scalen,ra_ref,dec_ref
	double precision   freq(nfmax),sigma(nfmax),antwdth(nfmax),f_ra(nfmax),f_dec(nfmax)
	character*500 datafn(nfmax),noisefn(nfmax)
        double precision OBSRA,OBSDEC,CRVAL4,CRVAL5,CRVAL6,PSCAL1,PSCAL2
        character(len=100) TELESCOP

	integer  i,j
	double precision   PI,dtor,slim,d1,d2,d3
	integer confmodel
        logical new_fmt
	
	! initialise variables

	PI=4d0*datan(1d0)        ! this really does equal PI
	dtor=PI/180d0            ! degrees/radians conversion

	! read inputs

	write(*,*)
	write(*,*) ' Telescope primary beam modelled as', &
		     ' Gaussian at all observing frequencies'
	write(*,*)
	
	ra_ref=0.d0
	dec_ref=0.d0
	imagFlag=.false.

        ! Check whether input file is new or old format
        read(input_unit,150) datafn(1)
        backspace(input_unit)
        inquire(file=datafn(1), exist=new_fmt)
        if (new_fmt) then
          i = index(datafn(1), '.', .true.)
          if (datafn(1)(i+1:i+4)/='fits') then
            write(*,*) 'Warning: assuming '//trim(datafn(1))//' is a FITS file'
          endif
          do i=1,nf
            read(input_unit,150) datafn(i)
            call phitsread_open(datafn(i), CRVAL4, CRVAL5, CRVAL6, OBSRA, OBSDEC, TELESCOP, PSCAL1, PSCAL2)
            !write(*,*) trim(datafn(i)), CRVAL4, CRVAL5, CRVAL6, OBSRA, OBSDEC, trim(TELESCOP), PSCAL1, PSCAL2
            freq(i) = CRVAL4*1d-9 ! GHz
            f_ra(i) = CRVAL5*dtor ! radians
            f_dec(i) = CRVAL6*dtor
	    ! Calculate reference RA and DEC and check for mosaicing
	    ra_ref=ra_ref+f_ra(i)	      
	    if(.not.imagFlag .and. ra_ref /= f_ra(i)*i) imagFlag=.true.
	    dec_ref=dec_ref+f_dec(i)
	    if((.not.imagFlag) .and. (  (  dec_ref - (f_dec(i)*i)  ) .GE. 9.0d-16  ) ) imagFlag=.true.
            ! Work out primary beam sigma and antenna width 
            if(TELESCOP(1:2)=='AM') then
               sigma(i)=pb_ami(1)/freq(i)+pb_ami(2)
               antwdth(i)=3.7*CRVAL4/clight
            elseif(TELESCOP(1:2)=='LA') then
               sigma(i)=pb_la(1)/freq(i)+pb_la(2)
               antwdth(i)=12.8*CRVAL4/clight
            endif
            sigma(i)=sigma(i)/60. ! deg
            call phitsread_close
          enddo
        else
	  do i=1,nf
		read(input_unit,*) freq(i) 	!frequency in FHz
	      
		read(input_unit,*) sigma(i)	!FWHM in degrees
	      
		read(input_unit,*) antwdth(i)	!antenna width
	      
		read(input_unit,*) d1,d2,d3	!RA(h,m,s)
	      	
		!convert to radians
	      	call sla_ctf2r(d1, d2, d3, f_ra(i), j)
	      	if(j /= 0) then
	      		write(*,*)"Incorrect RA"
	            	stop
		endif
	      
	      	!reference RA
	      	ra_ref=ra_ref+f_ra(i)
	      
	      	!mosaicking
	      	if(.not.imagFlag .and. ra_ref /= f_ra(i)*i) imagFlag=.true.
		
	      	read(input_unit,*) d1,d2,d3	!DEC(d,m,s)
	      	!convert to radians
	      	call sla_daf2r(d1, d2, d3, f_dec(i), j)
	      	if(j /= 0) then
	      		write(*,*)"Incorrect DEC"
	            	stop
		endif
	      
	      	!reference DEC
	      	dec_ref=dec_ref+f_dec(i)
	      
	      	!mosaicking
	     ! 	if(.not.imagFlag .and. dec_ref /= f_dec(i)*i) imagFlag=.true.
		if((.not.imagFlag) .and. (  (  dec_ref - (f_dec(i)*i)  ) .GE. 9.0d-16  ) ) imagFlag=.true.	     
	     
	  enddo
        endif
	
	ra_ref=ra_ref/nf
	dec_ref=dec_ref/nf
	
	read(input_unit,*) freq0		!reference frequency in GHz
	read(input_unit,*) slim 		!limiting flux for source confusion in mJy
	read(input_unit,*) confmodel 		!confusion noise model
	
	if( confmodel < 1 .or. confmodel > 3 ) then
		write(*,*)"Wrong model index for confusion noise."
		write(*,*)"Aborting!"
		stop
	endif
	
	!calculate the conf_cl
	if(slim>0.d0) then
		call calc_c_conf(confmodel,freq0,dble(slim),conf_cl)
	else
		conf_cl=0.d0
	endif

! calculate sigma of primary and 1/e diameter of antenna function

	write(*,*)
	antmax=0d0
	do i=1,nf
	  	sigma(i) = sigma(i)*dtor
	  	if (.not.new_fmt) sigma(i) = sigma(i)/(2*sqrt(2d0*dlog(2d0)))
		!antwdth(i) = 2d0/(PI*sigma(i)*dsqrt(2d0))
		!write(*,'(a,a,f6.2,a,f6.2,a)') '  1/e diameter of aperture ','function at ',freq(i),' GHz = ',antwdth(i),' wavelengths'
	  	if (antwdth(i) > antmax) antmax=antwdth(i)
	enddo
 
! ask for data file names and open

	do i=1,nf
	  	if (.not.new_fmt) read(input_unit,150) datafn(i)
	  	noisefn(i)=datafn(i)
	enddo
 150  	format(a)
	
	scaled=1
	scalen=1
	iaddn=0
	nskip=0
	
	write(*,*)"number of frequencies=",nf
	write(*,*)"frequencies in GHz: ",freq(1:nf)
	write(*,*)"sigma in degrees :",sigma(1:nf)
	write(*,*)"antenna width in wavelengths: ",antwdth(1:nf)
	write(*,*)"reference frequency in GHz: ",freq0
	write(*,*)"Reference centre: ", ra_ref / MIN2RAD, dec_ref / MIN2RAD, " arcmin ", ra_ref, dec_ref, "radians"
	write(*,*)"limiting flux for source confusion (in mJy): ",slim
	write(*,*)"data file named: "
	do i=1,nf
		write(*,*)trim(datafn(i))
	enddo
	
	return
end subroutine fread_inps

!-----*-----------------------------------------------------------------

subroutine uvlimits(nf,nskip,datafn,uvmin,uvmax)
	implicit none
	integer  nf,nskip
	double precision   uvmin,uvmax
	character*500 datafn(nf)

	integer  i,j,idum,nunit,ind
	double precision   uvdist,upos,vpos,vr,vi,rms
	double complex vis_buffer
        double precision uv_buffer(2)
        double precision OBSRA,OBSDEC,CRVAL4,CRVAL5,CRVAL6,PSCAL1,PSCAL2
        character(len=100) TELESCOP

! read data files to find min and max uv distances

	uvmin=1d10
	uvmax=0d0
	do i=1,nf
          ind = index(datafn(i), '.', .true.)
          if (datafn(1)(ind+1:ind+4)/='fits') then
	  	nunit=12
	  	open(unit=nunit, file=datafn(i), status='OLD') 
 350    	continue
	  	do j=1,nskip
	    		read(nunit,*,end=360) idum,upos,vpos,vr,vi,rms
	  	enddo
	  	read(nunit,*,end=360) idum,upos,vpos,vr,vi,rms
	  	uvdist=sqrt(upos**2+vpos**2)
	  	if (uvdist.lt.uvmin) uvmin=uvdist
	  	if (uvdist.gt.uvmax) uvmax=uvdist
	  	goto 350
 360    	close(nunit)
          else
                call phitsread_open(datafn(i), CRVAL4, CRVAL5, CRVAL6, OBSRA, OBSDEC, TELESCOP, PSCAL1, PSCAL2)
                do j=1,GCOUNT
                    call phitsread_one(j,vis_buffer,uv_buffer,rms,rms,idum, CRVAL4, PSCAL1, PSCAL2)
                    uvdist=sqrt(uv_buffer(1)**2+uv_buffer(2)**2)
	  	    if (uvdist.lt.uvmin) uvmin=uvdist
	  	    if (uvdist.gt.uvmax) uvmax=uvdist
                enddo
                call phitsread_close
          endif
	enddo
	
end subroutine uvlimits

!-----*-----------------------------------------------------------------

subroutine annuli(uvmin,uvmax,antmax,r,namax,na,cps1d)
	implicit none
	integer  namax,na
	double precision   uvmin,uvmax,antmax
	double precision   r(namax),cps1d(namax)

	integer  i
	double precision   umin,umax,ustep
	character autoyn*1

	! write out automatic u range and number of annuli
	
	write(*,*)
	write(*,*) ' CMB Power spectrum required in several bins'
	write(*,*)

	umin=uvmin-antmax/2d0
	if (umin.lt.0d0) umin=0d0
	umax=uvmax+antmax/2d0
	if(umax*2.*PI>8000d0) umax=8000d0/(2d0*PI)
	if(umin*2.*PI<2d0) umin=2d0/(2d0*PI)
	
	na=int((umax-umin)/antmax)
	write(*,'(a,f7.2,a,f7.2,a)') '  uv-range measured from ',uvmin,' to ',uvmax,' wavelengths'
	write(*,'(a,f7.2,a,f7.2,a)') '  power spectrum required from ',umin,' to ',umax,' wavelengths'
	write(*,*) ' number of bins = ',na
	if (na.gt.namax) then
	      	write(*,*)'ERROR - too many annuli, using ',namax,'annuli'
		na=namax
	endif
	ustep=(log(umax)-log(umin))/float(na)
	write(*,'(a,f6.2,a)') '  width  of bins = ',ustep,' wavelengths'
	do i=1,na+1
	  	r(i)=umin+ustep*(i-1)
	enddo

	! use different number of annuli if required

	write(*,*)
	write(*,'(a,$)')' Use automatic bin width above? (y/n) [y] '
	read(*,'(a1)') autoyn
	if (autoyn.eq.'n' .or. autoyn.eq.'N') then
	  	write(*,'(a,i5,a)')' Input number of bins in this range (',namax,'max.): '
	  	read(*,*) na
	  	if (na.gt.namax) then
	      		write(*,*)'ERROR - too many annuli, using ',namax,'annuli'
	  		na=namax
	  	endif
	  
	  	ustep=(log(umax)-log(umin))/float(na)
	  	write(*,*)
	  	write(*,'(a,f6.2,a,f6.2,a)') '  u-range = ',umin,' to ',umax,' wavelengths'
	  	write(*,*) ' number of bins = ',na
	  	write(*,'(a,f6.2,a)') '  width of bin = ',ustep,' wavelengths'
	  	write(*,*)
	  	do i=1,na+1
	    	r(i)=exp(log(umin)+(i-1.)*ustep)
	  	enddo
	endif
	
	! input power spectrum

	! 1d PS are in l^2 C_l/2pi in dT/T units; since l=2pi u this means they
	! are in 2pi u^2 S(u) in dT/T units. Calculate conversion factor to
	! give u^2 S(u) in (Jy/sr)^2 at freq0.
	
	write(*,*)
	write(*,*)'Input the power spectrum in l^2 C_l/2pi in dT/T units &
		in the following bins'
	do i=1,na
		write(*,*)'from ',r(i),'to ',r(i)+ustep
	      	read(*,*)cps1d(i)
	enddo
	
end subroutine annuli

!-----*-----------------------------------------------------------------

subroutine fannuli(uvmin,uvmax,antmax,r,namax,na,cps1d,confps1d,conf_cl)
	implicit none

#ifdef MPI
  	include 'mpif.h'
  	integer errcode
#endif

	integer  namax,na
	double precision   uvmin,uvmax,antmax
	double precision   r(namax+1),cps1d(namax),confps1d(namax),conf_cl

	integer  i, lmax
	double precision   umin,umax,ustep
	character psfile*500
	
	if( my_rank == 0 ) then
		! write out automatic u range and number of annuli
	
		write(*,*)
		write(*,*) ' CMB Power spectrum required in several bins'
		write(*,*)
		
		write(*,*)"input the name of the file with CMB PS in l^2 C_l/2pi in dT/T units"
		read(input_unit,'(a)') psfile
		
		lmax = find_lmax(psfile)
		if( lmax < 2 ) lmax = 8000
	endif
	
#ifdef MPI
    	call MPI_BARRIER(MPI_COMM_WORLD,errcode)
	
	call MPI_BCAST(lmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
#endif

	umin=uvmin-antmax/2d0
	if (umin.lt.0d0) umin=0d0
	umax=uvmax+antmax/2d0
	if(umax*2.*PI>dble(lmax)) umax=dble(lmax)/(2d0*PI)
	if(umin*2.*PI<2d0) umin=2d0/(2d0*PI)
	
	ustep=(log(umax)-log(umin))/float(na)
	
	if( my_rank == 0 ) then
		write(*,*)
		write(*,'(a,E20.12,a,E20.12,a)') '  u-range = ',umin,' to ',umax,' wavelengths'
		write(*,*) ' number of bins = ',na
		write(*,'(a,f6.2,a)') '  width of bin = ',ustep,' wavelengths'
		write(*,*)
	endif
	
	do i=1,na+1
		r(i)=exp(log(umin)+(i-1.)*ustep)
	enddo

	if( my_rank == 0 ) then	
		! input power spectrum

		! 1d PS are in l^2 C_l/2pi in dT/T units; since l=2pi u this means they
		! are in 2pi u^2 S(u) in dT/T units. Calculate conversion factor to
		! give u^2 S(u) in (Jy/sr)^2 at freq0.

		call calc_binned_dl(na,psfile,umin,umax,ustep,cps1d,confps1d,conf_cl)
	endif
	
end subroutine fannuli

!-----*-----------------------------------------------------------------

! note: 7/5/97 - rms now refers to noise on real or imaginary part

subroutine read_data(nvmax,namax,na,nf,nskip,nv,freqn,u,v, &
		visr,visi,rms,freq,r,scaled,scalen,datafn,noisefn,iaddn)
	implicit none
	integer  nvmax,namax,na,nf,nv,nskip,iaddn
	integer  freqn(nvmax)
	double precision   scaled,scalen
	double precision   u(nvmax),v(nvmax)
	double precision   visr(nvmax),visi(nvmax)
	double precision   rms(nvmax)
	double precision   freq(nf),r(namax+1)
	character*500 datafn(nf),noisefn(nf)

	integer  i,j,idum,count,ind
	double precision   upos,vpos,vr,vi,vrn,vin,rmserr,uvdist,rdum
	double precision   uplot,vplot
	double complex vis_buffer
        double precision uv_buffer(2)
        double precision OBSRA,OBSDEC,CRVAL4,CRVAL5,CRVAL6,PSCAL1,PSCAL2
        character(len=100) TELESCOP
	

	! initialise variables

	nv=0

	!...begin loop over observed frequencies
	do i=1,nf
	  count=0
          ind = index(datafn(i), '.', .true.)
          if (datafn(1)(ind+1:ind+4)/='fits') then
		!...open data and noise files
	  	open(unit=12, file=datafn(i), status='OLD') 
	  	open(unit=13, file=noisefn(i), status='OLD') 
	  	!open(unit=14, file=trim(datafn(i))//'.cov', status='REPLACE')
		!...read the data and noise
 150    	continue         
	  	do j=1,nskip
	    		read(12,*,end=160) idum,upos,vpos,vr,vi,rdum
	    		read(13,*,end=160) idum,rdum,rdum,vrn,vin,rmserr
	  	enddo
	  	read(12,*,end=160) idum,upos,vpos,vr,vi,rdum
	  	read(13,*,end=160) idum,rdum,rdum,vrn,vin,rmserr
		!...ignore data point if outside region of power spectrum estimation
	  	uvdist=dsqrt(upos**2+vpos**2)
	  	if ((uvdist >= r(na+1)).or.(uvdist < r(1))) goto 150
		!...update counters
	  	count=count+1
	  	nv=nv+1
	    	!write(14,'(i6,5E20.12)') count,upos,vpos,vr,vi,rmserr
	  	if (nv>nvmax) then
			write (*,*) ' ERROR: more than ',nvmax,' visibilities'
	      		stop
	  	endif
		!...scale visibility and noise
	  	vr=vr*scaled
	  	vi=vi*scaled
	  	vrn=vrn*scalen
	  	vin=vin*scalen
	  	rmserr=rmserr*scalen
		!...plot measured uv-position
	  	uplot=dble(upos)
	  	vplot=dble(vpos)
	  	!call pgpoint(1,uplot,vplot,1)  ! plot data point
	      
	  	u(nv)=upos
	  	v(nv)=vpos
	  	if (iaddn.eq.1) then
	  		visr(nv)=vr+vrn
	  		visi(nv)=vi+vin
	  	else
	  		visr(nv)=vr
	  		visi(nv)=vi
	  	endif
	  	rms(nv)=rmserr   ! input rms on Re/Im part
		!rmsa(nv)=rmserr/dsqrt(2.0) ! input rms on modulus
	  	freqn(nv)=i
	  	!call pgsci(mod(1,7)+1)
	  	!call pgpoint(1,uplot,vplot,1)         ! plot data point
	  	goto 150
		!...close data and noise files
 160    	close(12)
	  	close(13)
	  	!close(14)
          else
                call phitsread_open(datafn(i), CRVAL4, CRVAL5, CRVAL6, OBSRA, OBSDEC, TELESCOP, PSCAL1, PSCAL2)
                do j=1,GCOUNT
                    call phitsread_one(j,vis_buffer,uv_buffer,rmserr,rdum,idum, CRVAL4, PSCAL1, PSCAL2)
                    upos=uv_buffer(1)
                    vpos=uv_buffer(2)
                    vr=real(vis_buffer)
                    vi=aimag(vis_buffer)
                    uvdist=sqrt(upos**2+vpos**2)
  	  	    if((uvdist >= r(na+1)) .and. (uvdist < r(1))) cycle
                    rmserr=1d0/sqrt(rmserr)
		    !...update counters
	  	    count=count+1
	  	    nv=nv+1
	    	    !write(14,'(i6,5E20.12)') count,upos,vpos,vr,vi,rmserr
	  	    if (nv>nvmax) then
		        write (*,*) ' ERROR: more than ',nvmax,' visibilities'
	      		stop
	  	    endif
		    !...scale visibility and noise
	  	    vr=vr*scaled
	  	    vi=vi*scaled
                    vrn=vr
                    vin=vi
	  	    rmserr=rmserr*scalen
		    !...plot measured uv-position
	  	    uplot=dble(upos)
	  	    vplot=dble(vpos)
	  	    !call pgpoint(1,uplot,vplot,1)  ! plot data point
	      
	  	    u(nv)=upos
	  	    v(nv)=vpos
	  	    if (iaddn.eq.1) then
	  		visr(nv)=vr+vrn
	  		visi(nv)=vi+vin
	  	    else
	  		visr(nv)=vr
	  		visi(nv)=vi
	  	    endif
	  	    rms(nv)=rmserr   ! input rms on Re/Im part
		    !rmsa(nv)=rmserr/dsqrt(2.0) ! input rms on modulus
	  	    freqn(nv)=i
	  	    !call pgsci(mod(1,7)+1)
	  	    !call pgpoint(1,uplot,vplot,1)         ! plot data point

                enddo
                call phitsread_close
          endif
	  !...write out information 
	  write(*,'(a,f6.2,a,i5,a)') '  At ',freq(i),' GHz read ',count,' visibilities'
	!...close loop over observed frequencies
	enddo

	! write out total number of visibilities

	if (nf.gt.1) then
	  	write(*,*)
	  	write(*,'(a,i6)') '  Total number of visibilities = ',nv
	endif
	
!      datav(1,1:nv)=u(1:nv)
!      datav(2,1:nv)=v(1:nv)
!      datav(3,1:nv)=visr(1:nv)
!      datav(4,1:nv)=visi(1:nv)
!      datav(5,1:nv)=rms(1:nv)
!      do i=1,nv
!      	do j=i+1,nv
!            	if(u(j)**2+v(j)**2<datav(1,i)**2+datav(2,i)**2) then
!                  	temp(1:5)=datav(1:5,i)
!                  	datav(1:5,i)=datav(1:5,j)
!                  	datav(1:5,j)=temp(1:5)
!			endif
!		enddo
!	enddo
!      
!      open(unit=15, file='data.vbn', status='replace')
!      do i=1,nv
!      	write(15,'(i6,2F10.3,3E15.6)'),i,datav(1:5,i)
!	enddo
!      close(15)

	return
end subroutine read_data

!-----*-----------------------------------------------------------------

! note: 7/5/97 - rms now refers to noise on real or imaginary part

subroutine count_ndata(nvmax,nf,na,nskip,r,datafn,noisefn)
	implicit none
	integer  na,nf,nvmax,nskip,ind
	double precision   r(na+1)
	character*500 datafn(nf),noisefn(nf)
	double complex vis_buffer
        double precision uv_buffer(2)
        double precision OBSRA,OBSDEC,CRVAL4,CRVAL5,CRVAL6,PSCAL1,PSCAL2
        character(len=100) TELESCOP

	integer  i,j,idum
	double precision   upos,vpos,vr,vi,vrn,vin,rmserr,uvdist,rdum

	! initialise variables

	nvmax=0

	!...begin loop over observed frequencies
	do i=1,nf
	!...open data and noise files
          ind = index(datafn(i), '.', .true.)
          if (datafn(1)(ind+1:ind+4)/='fits') then
	  	open(unit=12, file=datafn(i), status='OLD') 
	  	open(unit=13, file=noisefn(i), status='OLD') 
		!...read the data and noise
 150    	continue         
	  	do j=1,nskip
	    		read(12,*,end=160) idum,upos,vpos,vr,vi,rdum
	    		read(13,*,end=160) idum,rdum,rdum,vrn,vin,rmserr
	  	enddo
	  	read(12,*,end=160) idum,upos,vpos,vr,vi,rdum
	  	read(13,*,end=160) idum,rdum,rdum,vrn,vin,rmserr
		!...ignore data point if outside region of power spectrum estimation
	  	uvdist=dsqrt(upos**2+vpos**2)
	  	if((uvdist >= r(na+1)) .or. (uvdist < r(1))) goto 150
		!...update counters
	  	nvmax=nvmax+1
	  	goto 150
		!...close data and noise files
 160    	close(12)
	  	close(13)
          else
                call phitsread_open(datafn(i), CRVAL4, CRVAL5, CRVAL6, OBSRA, OBSDEC, TELESCOP, PSCAL1, PSCAL2)
                do j=1,GCOUNT
                    call phitsread_one(j,vis_buffer,uv_buffer,rdum,rdum,idum, CRVAL4, PSCAL1, PSCAL2)
                    uvdist=sqrt(uv_buffer(1)**2+uv_buffer(2)**2)
  	  	    if((uvdist < r(na+1)) .and. (uvdist >= r(1))) then
		      !...update counters
	  	      nvmax=nvmax+1
                    endif
                enddo
                call phitsread_close
          endif
	!...close loop over observed frequencies
	enddo
	
end subroutine count_ndata

!-----*-----------------------------------------------------------------
	!calculate the binned PS values
subroutine calc_binned_dl(na,psfile,umin,umax,ustep,cps1d,confps1d,conf_cl)
	implicit none
	
	integer na
	double precision umin,umax,ustep,cps1d(na),confps1d(na),conf_cl
	character*500 psfile
	
	integer i,j,iostatus
	integer lmax,lmin
	double precision dl(nint(umax*2*3.1416)),PI,sig,ulo,uhi
	
	PI=4d0*atan(1d0)
	lmax=nint(umax*2*PI)
	
	open(unit=33, file=psfile, status='OLD')
	do i=1,lmax
		read(33,*,IOSTAT=iostatus)j,dl(i)
	      	if(iostatus<0) then
	      		write(*,*)"CMB PS values in file:",psfile," do not go up to the required lmax of ",lmax
	            	stop
		endif
	enddo
	close(33)
	
	confps1d=0.d0
	uhi=umin
	lmax=nint(umin*2*PI)-1
	
	write(*,*)"total CMB PS bins: ",na
	write(*,*)"CMB & CONF PS in l^2 C_l/2pi in dT/T units in each bin"

	do i=1,na
		lmin=lmax+1
	      	ulo=uhi
	      	uhi=exp(log(umin)+i*ustep)
	      	lmax=max(nint(uhi*2*PI),lmin)
	      
	      	!CMB PS
		cps1d(i)=sum(dl(lmin:lmax))/dble(lmax-lmin+1)
	      	sig=sum(dl(lmin:lmax)**2)/dble(lmax-lmin+1)
	      	sig=sqrt(sig-cps1d(i)**2)
	      
	      	!Confusion noise PS
	      	do j=lmin,lmax
	      		confps1d(i)=confps1d(i)+conf_cl*dble(j*(j+1))
		enddo
	      	confps1d(i)=confps1d(i)/dble(lmax-lmin+1)
	      	!convert it in K^2/2pi from uK^2
	      	!confps1d(i)=(confps1d(i)*(1d-6)**2)/(2.d0*PI)
	
		write(*,'(3i7,2E20.12)')i,lmin,lmax,cps1d(i),confps1d(i)
	      	
	enddo
end subroutine calc_binned_dl

!-----*-----------------------------------------------------------------
	!calculate the max l value
integer function find_lmax(psfile)
	implicit none
	
	character*500 psfile
	integer i,iostatus,lmax
	double precision cl
	
	lmax = 0
	open(unit=33, file=psfile, status='OLD')
	do
		read(33,*,IOSTAT=iostatus)i, cl
	      	if(iostatus<0) exit
		if( i > lmax ) lmax = i
	enddo
	close(33)
	
	find_lmax = lmax
	
end function find_lmax
	            
!-----*-----------------------------------------------------------------
end module inputs
