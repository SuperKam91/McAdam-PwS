module MassModels
	use params
      use constants
      use matrix_utils
      use utilities
      use cosmology

contains

!=======================================================================

	subroutine LensFields(i,flag)

	implicit none
	
      integer i,i1,j,flag,k
	double precision x0,y0,thetaE
	double precision rs,ps,vdsq,k0,gamma,cc
	double precision rc,Sigma0,rproj
	parameter(rproj=1.0)
	double precision SigmaCrit
		
	double precision dx(2),rr,cosTwophi,sinTwophi,sn,d1
	double precision Twophi,rhocritz

	double precision drdx(2),d2rdxdy(2,2),drdxdrdy(2,2)
	double precision temp1(2,2),temp2(2,2)
	double precision psiprime,psi2prime,d2psidxdy(2,2)
      double precision NFWpsiprime,NFWpsi2prime,sigerr,totSN

!-----------------------------------------------------------------------
      if (Varyzs==1) then
        	!if the source plane is at a redshift less than the cluster, no shear
        	if(zs<z) return
	endif
      
      rhocritz=rhocritofz(z)
      
! 	Set up atom parameters for the ith atom:
	
	x0=GeoPars(i,1)
	y0=GeoPars(i,2)

      if (MassModel==1) then

!     NFW:

		if (NFWstyle(i)==1) then 
	    		rs=MassPars(i,1) 
	    		ps=MassPars(i,2) 
	  	elseif (NFWstyle(i)==2) then
	    		cc=MassPars(i,1)
	    		M200=MassPars(i,2)
	    		if (cc<=0. .or. M200<=0.) then
	      		flag=1
            		airballs=airballs+1
	      		goto 999
	    		endif
	    		r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
	    		rs=r200/cc
	    		ps=M200/(4.0*pi*rs*rs*rs*(dlog(1.+cc)-cc/(cc+1.)))
          
	  	elseif (NFWstyle(i)==3) then

! 			M200 now contains M100, so that Komatsu and Seljak's relation can be used!!
	    		M200=MassPars(i,2)
		    	cc=6.0*(M200/1.e14)**(-0.2)
		    	r200=(M200/(100.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
		    	rs=r200/cc
		    	ps=M200/(4.0*Pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
		elseif (NFWstyle(i)==4) then
	    		thetaE=MassPars(i,1)
	    		M200=MassPars(i,2)
	    		if (thetaE.le.0.0.or.M200.le.0.0) then
	      		flag=1
            		airballs=airballs+1
	      		goto 999
	    		endif
	    	
            	rE=ThetaE*sec2rad*D
	    		scE=sigcritE
	    		r200=(M200/(200.0*(4.0*Pi/3.0)*rhocritz))**(1.0/3.0)
	    		rs=FINDR(rsfunc,1.0d-4,1.0d4,1.0d-4)
	    		if (rs==0.0) then
	      		flag=1
            		airballs=airballs+1
	      		goto 999
	    		else  
	      		cc=r200/rs
	      		ps=M200/(4.0*Pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
	    		endif
        	else
          		STOP 'Unrecognised NFWstyle.'
	  	endif  

      elseif (MassModel==2) then

!     	SIS:

        	if (SISstyle(i)==1) then 
! 	    		vdsq=MassPars(i,1)*SigmaSq2G
	    		vdsq=200.0*Pi*rhocritz*MassPars(i,1)*MassPars(i,1)/48.0
	    		vdsq=vdsq**(0.3333)
        	elseif (SISstyle(i)==2) then 
	    		thetaE=MassPars(i,1)
	    		vdsq=thetaE*sec2rad*D*sigcritE/2.d0
        	else
          		STOP 'Unrecognised SISstyle.'
	  	endif
	  
        
      elseif (MassModel==3) then

!     	Power law profile-negative numbers indicate CIS:

        	if (CPLstyle(i)==1) then 
		    	rc=MassPars(i,1)/1000.d0
	    		Mproj=MassPars(i,2)
		    	slope=MassPars(i,3)
		    	cc=(rproj/rc)*(rproj/rc)
	    		Sigma0=Mproj/(Pi*rc*rc*(cc/(1.0+cc))*(1+cc)**slope)
		elseif (CPLstyle(i)==2) then 
	    		thetaE=MassPars(i,1)
			Mproj=MassPars(i,2)
	    		slope=MassPars(i,3)
	    		rE=ThetaE*sec2rad*D
	    		rp=rproj
	    		scE=sigcritE
	    		rc=FINDR(rcfunc,1.0d-4,1.0d4,1.0d-4)
	    		if (rc==0.0) then
	      		flag=1
            		airballs=airballs+1
	      		goto 999
	    		else  
	      		cc=(rproj/rc)*(rproj/rc)
	      		Sigma0=Mproj/(Pi*rc*rc*(cc/(1.0+cc))*(1+cc)**slope)
          		endif
        	elseif (CPLstyle(i)==-1) then 
	    		rc=MassPars(i,1)/1000.d0
	    		M200=MassPars(i,2)
	    		slope=MassPars(i,3)
          		if (slope.ne.0.5d0) then
            		write(*,*) ' Fix slope at 0.5 to use CIS model.'
            		stop
          		endif  
	    		r200=(M200/(200.0*(4.0*Pi/3.0)*rhocritz))**(1.0/3.0)
	    		cc=(r200/rc)*(r200/rc)
	    		Sigma0=200.0*rhocritz*(2.0*pi/3.0)*rc*(1.0+cc)
        	else
          		STOP 'Unrecognised CPLstyle.'
	  	endif  
	endif 

!-----------------------------------------------------------------------

	if(.not.simulate .and. sn_lim>0. .and. GeoModel==1) then
      	!look-up the S/N cut-off
      	call lookUp1D(n,cc,snlook(:,1),snlook(:,2),sn)
	endif
      
      if (Varyzs==1) then
      	!lookup SigmaCrit which is the same for all the galaxies
	  	call lookUp2D(SCZdn,SCZsn,z,zs,SigCritZd,SigCritZs,lookSigCrit,SigmaCrit)
        	SigmaCrit=SigmaCrit*1.d12
	elseif(vary_zs==0) then
      	!lookup SigmaCrit which is the same for all the galaxies
	  	call lookUp2D(SCZdn,SCZsn,z,z_s(1),SigCritZd,SigCritZs,lookSigCrit,SigmaCrit)
        	SigmaCrit=SigmaCrit*1.d12
	endif
      
      !convert the S/N cut-off to arcsec
	sn_asec=sn*r200/(sec2rad*D)
            
	!find the no. of pixels this cut-off radius correspond to
	snx_npix=sn_asec/pix_sizex
	sny_npix=sn_asec/pix_sizey
      
      !find out the pixel no. of the cluster center
      x0_npix=x0/nxpix+nxpix/2
      y0_npix=y0/nypix+nypix/2
      
      !loop over the interesting region
      do i1=1,nxpix
      	!x distance of the pixel under consideration in arcsec
      	xk=x0+(i1-1.)*pix_sizex
            !x pixel no.
            xk_npix=x0_npix+i1-1
            
      	do j1=1,nypix
      		!y distance of the pixel under consideration in arcsec
      		yk=y0+(j1-1.)*pix_size
            	!y pixel no.
            	yk_npix=y0_npix+j1-1

			!Find radial distance from centre of cluster in Mpc:
                  dx(1)=(xk-x0)*sec2rad*D
			dx(2)=(yk-y0)*sec2rad*D 
			rr=sqrt(dx(1)*dx(1)+dx(2)*dx(2))
      		
                  !S/N cut-off
      		if(.not.simulate .and. sn_lim>0.) then
				if(sn*r200<rr) exit
			end if
                  
			Twophi=2.0*atan2(1.0*dx(2),1.0*dx(1))
			cosTwophi=cos(Twophi)
			sinTwophi=sin(Twophi)

      		if (rr==0.) cycle
                  
                  kp=kappa(xk_npix,yk_npix)+NFWconvergence(rr,ps,rs,SigmaCrit)
                  gamma=NFWshear(rr,ps,rs,SigmaCrit)
                  gm1=cosTwophi*gamma
                  gm2=sinTwophi*gamma
                  
                  !exploit spherical symmetry
        		if(xk_npix>0 .and. xk_npix<=nxpix .and. yk_npix>0 .and. yk_npix<=nypix) then
                  	kappak(xk_npix,yk_npix)=kp
	  			gamma1k(xk_npix,yk_npix)=gm1
	  			gamma2k(xk_npix,yk_npix)=gm2
			endif
                  
                  if(xk_npix>0 .and. xk_npix<=nxpix .and. -yk_npix+2*y0_npix>0 .and. -yk_npix+2*y0_npix<=nypix) then
                  	kappak(xk_npix,-yk_npix+2*y0_npix)=kp
	  			gamma1k(xk_npix,-yk_npix+2*y0_npix)=gm1
	  			gamma2k(xk_npix,-yk_npix+2*y0_npix)=gm2
			endif
                  
                  if(-xk_npix+2*x0_npix>0 .and. -xk_npix+2*x0_npix<=nxpix .and. yk_npix>0 .and. yk_npix<=nypix) then
                  	kappak(-xk_npix+2*x0_npix,yk_npix)=kp
	  			gamma1k(-xk_npix+2*x0_npix,yk_npix)=gm1
	  			gamma2k(-xk_npix+2*x0_npix,yk_npix)=gm2
			endif
                  
                  if(-xk_npix+2*x0_npix>0 .and. -xk_npix+2*x0_npix<=nxpix .and. -yk_npix+2*y0_npix>0 .and. -yk_npix+2*y0_npix<=nypix) then
                  	kappak(-xk_npix+2*x0_npix,-yk_npix+2*y0_npix)=kp
	  			gamma1k(-xk_npix+2*x0_npix,-yk_npix+2*y0_npix)=gm1
	  			gamma2k(-xk_npix+2*x0_npix,-yk_npix+2*y0_npix)=gm2
			endif
		enddo
	enddo
	
 999  return
	end subroutine LensFields

!=======================================================================

	subroutine makeSNlookup

	implicit none
	
      integer i,j,flag,k,i1,j1
	double precision rs,ps,cc,rr,sn,da,sigerr,totSN,m1,l1,r1,rhocritz
      double precision cmin,cmax,cinc,M,conc

!-----------------------------------------------------------------------
	if(NFWstyle(1)==2) then
            if(Mass_PriorType(1,1)==0) then
               	j=1
		else
                  j=n
		endif
	endif
      
      cmin=Mass_Prior(1,1,1)
      cmax=Mass_Prior(1,1,2)
      
      !set up the step size
      cinc=(cmax-cmin)/j
      
	do j1=1,j
            !S/N Cut-off is independent of M200 if r200 is expressed in terms of r_s.
            !pick any mass
      	M200=1.
            !find c in the grid
            if(j1==1) then
            	conc=cmin
		elseif(j1==j) then
            	conc=cmax
		else
            	conc=cmin+cinc/2.+cinc*(j1-1.)
		endif
            
            !calculate rhocrit(z) & the angular diameter distance
            rhocritz=1.
		!D=1.
            !calculate r200 & rs
		r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
		rs=r200/conc
            !calculate total S/N (at 30*r200)
      	totSN=NFWSN(30.*r200,rs,0.5d0,1.d0,1.d0,40.d0,1.d0)
                  
            l1=0.1
	      r1=30.
      	do
      		m1=(l1+r1)/2.
            	rr=m1*r200
            
	      	sn=NFWSN(rr,rs,0.5d0,1.d0,1.d0,40.d0,1.d0)/totSN
      		if(sn<sn_lim*(1+.1/100.) .and. sn>sn_lim*(1-.1/100.)) then
            		snlook(j1,1)=conc
            		snlook(j1,2)=m1
                        exit
			elseif(sn<sn_lim) then
      	      	l1=m1
			elseif(sn>sn_lim) then
	            	r1=m1
			endif
		enddo
 	enddo
      	
	end subroutine makeSNlookup
	
!=======================================================================

	subroutine EllGeometry(i)
	
	implicit none

	integer i	
	double precision theta,f,RR(2,2),SS(2,2),temp1(2,2),temp2(2,2)

!-----------------------------------------------------------------------
	
! Construct the transformation matrix Q for elliptical geometry fo
! the ith atom:

!   theta=orientation angle of major axis anticlockwise from x axis
!   f=axis ratio b/a < 1
	
	theta=GeoPars(i,3)*Pi/180.d0
	f=GeoPars(i,4)
	
	RR(1,1)=cos(theta)	  
	RR(1,2)=sin(theta)	  
	RR(2,1)=-sin(theta)	  
	RR(2,2)=cos(theta)	  
	SS(1,1)=sqrt(f)	  
	SS(1,2)=0.0	  
	SS(2,1)=0.0	  
	SS(2,2)=1.0/sqrt(f)	  
	  
	call MatrixProduct(2,2,SS,2,2,RR,temp1)
	call Transpose(2,2,temp1,temp2)
	call MatrixProduct(2,2,temp2,2,2,temp1,Q)
	
	return
	end subroutine EllGeometry
	
!=======================================================================

	function NFWconvergence(r,ps,rs,sigcrit)
	
	implicit none
	
	double precision r,ps,rs,sigcrit
	double precision x,k0,kappa
	double precision t1,t2,bracket
	double precision NFWconvergence,atanh
	
	k0=rs*ps/sigcrit
	x=r/rs
	
	if (x==1.) then
	  
	  kappa=2.d0*k0/3.d0
	
	elseif (x<1.) then
	
	  t1=1.0/(1.0-x*x)
	  t2=dsqrt((1.0-x)/(1.0+x))
	  bracket=1.d0-2*dsqrt(t1)*atanh(t2)
	  kappa=-2.d0*k0*bracket*t1

 	elseif (x>1.) then
	
	  t1=1.d0/(x*x-1.d0)
	  t2=dsqrt((x-1.d0)/(1.d0+x))
	  bracket=1.d0-2*dsqrt(t1)*atan(1.0*t2)
	  kappa=2.d0*k0*bracket*t1
	  
	endif  
	
	NFWconvergence=kappa
	
	return
	end function NFWconvergence

!-----------------------------------------------------------------------

	function NFWshear(r,ps,rs,sigcrit)
	
	implicit none
	
	double precision r,ps,rs,sigcrit
	double precision NFWshear,atanh
	double precision x,k0,gamma
	double precision t1,t2,bracket
	
	k0=rs*ps/sigcrit
	x=r/rs
	
	if (x==1.d0) then
	 
	  gamma=k0*(10.d0/3.d0+4.d0*dlog(0.5d0))
	  
	elseif (x<1.d0) then
	
	  t1=1.0/(1.0-x*x)
	  t2=dsqrt((1.d0-x)/(1.d0+x))
        bracket=atanh(t2)*dsqrt(t1)
	  bracket=bracket*(8.d0/(x*x)-4.d0*t1)
	  bracket=bracket+4.d0*dlog(x/2.d0)/(x*x)
        bracket=bracket+2.d0*t1
        gamma=k0*bracket
	  
	elseif (x>1.d0) then
	
	  t1=1.d0/(x*x-1.d0)
	  t2=dsqrt((x-1.d0)/(1.0+x))
	  bracket=atan(1.0*t2)*dsqrt(t1)
        bracket=bracket*(8.d0/(x*x)+4.d0*t1)
        bracket=bracket+4*dlog(x/2.d0)/(x*x)
        bracket=bracket-2.d0*t1
        gamma=k0*bracket

 	endif  
	
	NFWshear=gamma
	
	return
	end function NFWshear

!-----------------------------------------------------------------------
	!returns the function r*g(r)^2 in the calculation of shear at radius r (in kpc)
	double precision function gs(r)
	
	implicit none
	
	double precision r
      
      gs=r*((NFWshear(r,1.d0,1.d0,1.d0))**2.)
      
	end function gs

!-----------------------------------------------------------------------
	!returns the signal-to-noise of a NFW cluster at 
	double precision function NFWSN(rr,rs,rwidth,rhos,sigcrit,ndensity,sigerr)
	
	implicit none
	
      !input parameters
	double precision rr !distance (in Mpc) at which the S/N needs to be calculated
      double precision rs !scale radius (in Mpc) of NFW profile, =r200/c
      double precision rwidth !the width of the annulus (in arcmin) in which the S/N is to be measured
      double precision rhos !scale density of NFW profile (in M_sun/Mpc^3)
      double precision sigcrit !critical surface mass density of the lens (in M_sun/Mpc^2)
      double precision ndensity !no. density of galaxies (in arcmin^-2)
      double precision sigerr !dispersion of the ellipticities
      
      !work parameters
      double precision S,rin,rout,nd,rw
      double precision eps
      parameter(eps=1d-4)
      
      !convert number density from arcmin^-2 to Mpc^-2
      !nd=ndensity/((60.*sec2rad*D)**2.)
      
      !convert the rwidth from arcmin to Mpc
      !rw=rwidth*60.*sec2rad*D
      
      !calculate the inner & outer radii of the annulus
      !rin=(rr-rw/2.)/rs
      !rout=(rr+rw/2.)/rs
      
      rin=0.00001
      rout=rr/rs
      !rout=rr*60.*sec2rad*D/rs
      
      !performs the integral over x.g(x) from rin to rout
      call qtrap(gs,rin,rout,eps,S)
      
      !calculate the S/N
      !NFWSN=2.*sqrt(Pi*nd)*rhos*(rs**2.)*sqrt(S)/(sigcrit*sigerr*D)
      NFWSN=sqrt(S)
      
	end function NFWSN

!-----------------------------------------------------------------------

	function NFWpsiprime(x,k0)
	
	double precision x,k0
	double precision k1
	double precision t1,t2,bracket
	double precision NFWpsiprime,atanh
	
	k1=4.0*k0
	
	if (x==1.0) then
	  
	  NFWpsiprime=k1*(1.d0-log(2.0))
	
	elseif (x<1.0) then
	
	  t1=2.d0*dsqrt(1.d0/(1.d0-x*x))
	  t2=dsqrt(abs((1.d0-x)/(1.d0+x)))
	  bracket=dlog(0.5d0*x)+t1*atanh(t2)
	  NFWpsiprime=k1*bracket/x

 	elseif (x>1.0) then
	
	  t1=2.d0*dsqrt(1.d0/(x*x-1.d0))
	  t2=dsqrt((x-1.d0)/(1.d0+x))
	  bracket=dlog(0.5d0*x)+t1*atan(t2)
	  NFWpsiprime=k1*bracket/x
	  
	endif  
		
	return
	end function NFWpsiprime

!-----------------------------------------------------------------------

	function NFWpsi2prime(x,k0)
	
	double precision x,k0
	double precision k1
	double precision t1,t2,t3,bracket
	double precision NFWpsi2prime,atanh
	
	k1=4.d0*k0
	
	if (x==1.0) then
	  
	  NFWpsi2prime=k1*(log(2.0)-2.0/3.0)
	
	elseif (x<1.0) then
	
	  t1=2.d0*(1.d0-2*x*x)/(1.d0-x*x)
	  t1=t1*dsqrt(1.d0/(1.d0-x*x))
	  t2=dsqrt((1.d0-x)/(1.d0+x))
	  t3=(1.d0/(1.d0-x*x))-1.d0
	  bracket=-dlog(0.5d0*x)-t1*atanh(t2)-t3
	  NFWpsi2prime=k1*bracket/(x*x)
	  
 	elseif (x>1.0) then
	
	  t1=2.d0*(1.d0-2*x*x)/(1.d0-x*x)
	  t1=t1*dsqrt(1.d0/(x*x-1.d0))
	  t2=dsqrt((x-1.d0)/(1.d0+x))
	  t3=(1.d0/(1.d0-x*x))-1.d0
	  bracket=-dlog(0.5d0*x)-t1*atan(t2)-t3
	  NFWpsi2prime=k1*bracket/(x*x)
	  
	endif  
	  		
	return
	end function NFWpsi2prime

!=====================================================================

      function PowerLawconvergence(r,rc,k0,a1)
	
	implicit none
	
	double precision r,rc,k0,a1
	double precision x,kappa
	double precision PowerLawconvergence
	
	x=r/rc
	kappa=k0*((1.d0+a1*x*x)/((1.d0+x*x)**(2.d0-a1)))
	
	PowerLawconvergence=kappa
			
	return
	end function PowerLawconvergence

!---------------------------------------------------------------------

      function PowerLawshear(r,rc,k0,a1)
	
	implicit none
	
	double precision r,rc,k0,a1
	double precision x,gamma
	double precision PowerLawshear
	
	x=r/rc
! 	gamma=2.d0*k0*(1.d0-a1)*x*x/((1.d0+x*x)**(2.d0-a1))
	gamma=k0*(1.d0-a1)*x*x/((1.d0+x*x)**(2.d0-a1))
			
	PowerLawshear=gamma
! 	write(*,*) 'In PowerLawShear, r,rc,k0,a1,gamma: ',r,rc,k0,a1,gamma
	
	return
	end function PowerLawshear

!=====================================================================
	
	function rcfunc(rc)
	
	implicit none
	
	double precision rcfunc,rc
	double precision sE,sR,a1,term1,term2	
	
	sE=rE/rc
	sR=rp/rc
	a1=slope
	
	term1=1.0/(1.0+sE*sE)**(a1-1.d0)

! 	write(*,*) 'rE=',rE
! 	write(*,*) 'xE=',xE
! 	write(*,*) 'r200=',r200
! 	write(*,*) '!=',!
! 	write(*,*) 'sigcrit4=',sigcrit4
! 	write(*,*) 'M4pi=',M4pi
! 	write(*,*) 'g=',g
	
	
	term2=(Mproj/(Pi*scE))*(1.0+sR*sR)/(rc*rc*sR*sR*(1.d0+sR*sR)**a1)
! 	write(*,*) 'term1,term2=',term1,term2
	
	rcfunc=term1-term2
	
	return
	end function rcfunc
     
!=======================================================================
	subroutine makeSigCritlookup
      
      	implicit none
            
            integer i,j,k
            double precision zdmin,zdmax
            
            !determine the min & max values that lens redshift can take
            zdmin=z_Prior(1)
      	zdmax=z_Prior(2)
      	
            !iterate over lens redshift
            do i=1,SCZdn
                  if(i==SCZdn) then
                  	SigCritZd(i)=zdmax
			else
                  	SigCritZd(i)=10.**(log10(zdmin)+ &
                        (log10(zdmax)-log10(zdmin))*(i-1.)/(SCZdn-1.))
			endif
                  
                  !iterate over source redshift
                  do j=1,SCZsn
                  	if(i==1) then
                  		if(j==SCZsn) then
                  			SigCritZs(j)=zsmax
					else
                  			SigCritZs(j)=10.**(log10(zsmin)+ &
                              	(log10(zsmax)-log10(zsmin))*(j-1.)/(SCZsn-1.))
					endif
				endif
                        
                        !can't calculate sigCrit if zd<zs
                        if(SigCritZd(i)==SigCritZs(j)) then
                        	lookSigCrit(i,j)=0.
                              cycle
				!calculate sigCrit
                        else
                        	lookSigCrit(i,j)=criticaldensity(SigCritZd(i),SigCritZs(j))
				endif
			enddo
		enddo
                        
	end subroutine makeSigCritlookup

!-----------------------------------------------------------------------
	subroutine makeDlookup
      
      	implicit none
            
            integer i,j,k
            double precision zdmin,zdmax
            
            !determine the min & max values that lens redshift can take
            zdmin=z_Prior(1)
      	zdmax=z_Prior(2)
      	
            !iterate over lens redshift
            do i=1,Dn
                  if(i==Dn) then
                  	lookD(i,1)=zdmax
			else
                  	lookD(i,1)=10.**(log10(zdmin)+ &
                        (log10(zdmax)-log10(zdmin))*(i-1.)/(Dn-1.))
			endif
                  lookD(i,2)=angdist(lookD(i,1))
		enddo
                        
	end subroutine makeDlookup

!-----------------------------------------------------------------------
end module MassModels
