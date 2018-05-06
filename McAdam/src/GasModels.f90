module GasModels
	use params
      	use constants
      	use telescope1
	use utilities
      	use matrix_utils
      	use fftw_inc
      	!use rfft
      	use cosmology
	use CheckPars1
      
      	implicit none
      	double precision lmin1,lmax1,gd1,gd2

contains

!=======================================================================

	subroutine MakeGasDistributions(k,flag)

! Make 1-dimensional distributions needed for modelling cluster gas in
! the kth atom.
! Almost certainly needs tidying up, including further modularisation.
! Models to be included should be a HSE profile based on beta model for
! gas.
! Note parts returning non-zero flags moved to CheckGasPars module to facilitate fast/slow parameter separation in PolyChord (YCP 31/7/15)

	implicit none
			
        integer i,j,k,flag
	
        double precision A1,Gm,yc,index_g
        parameter(A1=1.424d-19,Gm=3.12d-14)
	
        double precision T0,logT0,rr,rlimit1
        double precision prefactor,prefactor2
	double precision eps, r1, r2
	parameter(eps=1d-4)

	double precision :: rsdivr200, frsdivr200, f500, p, xf500	!for analytical r_200 -> r_500 mapping for NFW DM kj 01/02/17
	double precision, parameter :: a_1=0.5116e0, a_2=-0.4283e0, a_3=-3.13e-3, a_4=-3.52e-5
	double precision :: X0, X !for Newton raphson r_200 -> r_500 mapping for Einasto DM kj 12/02/17
	
!-----------------------------------------------------------------------
	
      
        Phi = 0.d0
	Rhogas = 0.d0
	T = 0.d0
	logPhi = 0.d0
	logRhogas = 0.d0
	logT = 0.d0
	
!-----------------------------------
! Isothermal beta models:
!-----------------------------------

! Hobson and McLachlan beta model, normalised to central y parameter:
! Parameters of kth atom are 
!    GasPars(k,1)=apparent core radius (arcmin)
!    GasPars(k,2)=central y parameter\
! Beta is fixed at 2/3, and profile is "truncated" at rt=3*rc
! Note no need for D or z to take sensible values-the r array is
! interpreted as being in arcminutes not Mpc. This is ok even at z=0.05,
! where rmax=100 corresponds to 4 Mpc or so. Also note that this profile
! falls to 0.1% of its maximum at about 18 core radii.

      call CheckGasPars(k,flag)
!  Should never return a non-zero flag here because it's already been checked earlier, but just in case...
      if (.not. flag == 0) return


      if(GasModel==0) then
		rc=GasPars(k,1)
	  	yc=GasPars(k,2)
	   
        	do i=1,n
	    		Yarray(i)=BetaModelHM(r(i),rc,yc)
	    		if(Yarray(i)<0.0) then
	      		write(*,*) 'about to take log of negative number: '
	      		write(*,*) '  i,Yarray(i)=',i,Yarray(i)
	      		write(*,*) '  yc,rc=',yc,rc
	      		write(*,*) '  GasPars=',GasPars
	    		endif
	    		logYarray(i)=phlog10(Yarray(i))
	  	enddo
	  
	  	flag=0
	  	return
    ELSEIF(GasModel==4 )THEN
! Starting Beta_Atomic model

        thetac_BA= GasPars(k,1)
        Beta_BA=GasPars(k,2)
	dT0_BA=GasPars(k,3)*(1.d-6)
	aux(k,1)=dT0_BA
	  	flag=0
	  	return   	    	
!-----------------------------------

! Isothermal beta model, normalised to gas mass within 3D radius rmass if 
! Mass=0, or within 3D radius r200 if Mass=1; projection is numerical.
! Note tying of beta model core radii together by using CPLstyle-1:
         
      elseif (GasModel==1) then		
         
		!null run
		if( Mgas200 == 0d0 ) then
			flag = 2
			return
		endif

!       	Integrate profile to ensure model has correct gas mass
	  	do i=1,n
	    		Rhogas(i)=BetaModel3D(r(i),rc,beta,1.d0)
	    		logRhogas(i)=phlog10(Rhogas(i))
	  	enddo
        	
            	rr=r200
        	Rhogas0=Mgas200/GasMass(rr,rc,beta)
        	Rhogas_central=Rhogas0
	  	logRhogas0=phlog10(Rhogas0)
	  
!       	Tabulate T profile and compute pressure:        
	  	do i=1,n
		     	if(Tmodel==1.and. TStyle(k)==0) then
            			T(i)=TPars(k,1)
          		elseif(Tmodel==2.and. TStyle(k)==0) then
	      			index_g=TPars(k,2)
            			T(i)=TPars(k,1)*(Rhogas(i))**(index_g-1)
			elseif(Tmodel==0 .and. TStyle(k)==1) then
				T(i)=(8.2d0)* ( (M200*1.d-15*h)**(2.d0/3.d0) )*(((Om*(1.d0+z(k))**3.d0) +OL)** (1.d0/3.d0))       !scaling relation see Voit, 2004 eq. (59)
			elseif(Tmodel==0 .and. TStyle(k)==2) then           !HSE       		
				T(i)=((4.d0*pi*Gm*200.d0*rhocritz)/(9.d0*beta))*(r200*r200 +rc*rc)
			else
	      			stop 'Unknown temperature model.'
          		endif
          		
                  prefactor=A1*T(i)
	    		Rhogas(i)=Rhogas0*Rhogas(i)
	    		logRhogas(i)=logRhogas0+logRhogas(i)
	    		Pgas(i)=prefactor*Rhogas(i)
	    		logPgas(i)=phlog10(Pgas(i))
	  	enddo
	       	T200=T(n)
!-----------------------------------
! Hydrostatic equilibrium models, numerical projection:
!-----------------------------------

! 	Isothermal temperature profile:

	elseif(GasModel==2.and.TModel==1) then
	
!     		First calculate gas density and temperature:	
	
        	T(1:n)=TPars(k,1)
	    
        	if(MassModel==1) then
           	   
	    		prefactor=-4.d0*Pi*Gmu*ps*rs*rs/TPars(k,1)
	    
	    		do i=1,n
	      			rr=r(i)/rs
	      			Rhogas(i)=1.0*dexp(prefactor*((rr-log(1.0+rr))/rr))
	      			logRhogas(i)=phlog10(Rhogas(i))
	    		enddo

			!Now normalise to Mgas200, and calculate pressure:
	    		rr=r200
          		Rhogas0=GasMass(rr,rc,beta)
          		Rhogas_central=Rhogas0

 	    		if(Rhogas0 > 1d-20) then
	      			Rhogas0=GasPars(k,1)/GasMass(rr,rc,beta)
            			Rhogas_central=Rhogas0
	      			logRhogas0=phlog10(Rhogas0)
	      			prefactor=A1*TPars(k,1)
	      			do i=1,n
	        			Rhogas(i)=Rhogas0*Rhogas(i)
	        			logRhogas(i)=logRhogas(i)+logRhogas0
	        			Pgas(i)=1.0*prefactor*Rhogas(i)
	        			logPgas(i)=phlog10(Pgas(i))
	      			enddo
	    		else
	      			Rhogas0=0d0
            			Rhogas_central=Rhogas0
	      			logRhogas0=-20d0
	      			Rhogas(1:n)=0d0
	        		logRhogas(1:n)=-45d0
	        		Pgas(1:n)=0d0
	        		logPgas(1:n)=-45d0
	        		Yarray(1:n)=0d0
	        		logYarray(1:n)=-45d0
				return
	    		endif
	  	endif
	   	   	
!-----------------------------------

! 	Polytropic temperature profile:
         
	elseif(GasModel==2.and.TModel==2) then
		  
!     		First calculate potential, then the temperature, then the density:	
	
        		if(MassModel==1) then
        
	    			if(TPars(k,2)==1.0d0) then
	    				prefactor=-4.d0*Pi*Gmu*ps*rs*rs/TPars(k,1)
					T0=TPars(k,1)
					do i=1,n
	      					T(i)=TPars(k,1)
						rr=r(i)/rs
	      					Rhogas(i)=1.0d0*exp(prefactor*((rr-log(1.0d0+rr))/rr))
	      					logRhogas(i)=phlog10(Rhogas(i))
	    				enddo
	    			else
	    				prefactor=4.d0*Pi*Gmu*ps*rs*rs*(TPars(k,2)-1.0d0)/TPars(k,2)
	    				T0=1.d0*TPars(k,1)
	    				logT0=phlog10(T0)
	    				prefactor2=(1.d0/(TPars(k,2)-1.0d0))
	    				do i=1,n
	      					rr=r(i)/rs
	      					T(i)=TPars(k,1)-1.d0*prefactor*(((rr-log(1.d0+rr))/rr))
	      					logT(i)=phlog10(T(i))
	      					logRhogas(i)=prefactor2*(logT(i)-logT0)
	      					Rhogas(i)=10.d0**(logRhogas(i))
	    				enddo
	    			endif

!         			Now normalise to Mgas200, and calculate pressure:

	    			rr=1.0*r200
          			Rhogas0=GasMass(rr,rc,beta)
          			Rhogas_central=Rhogas0
 	    			if(Rhogas0 > 1d-20) then
	      				Rhogas0=GasPars(k,1)/GasMass(rr,rc,beta)
            				Rhogas_central=Rhogas0
	      				logRhogas0=phlog10(Rhogas0)
	      				prefactor=A1
	      				do i=1,n
	        				Rhogas(i)=Rhogas0*Rhogas(i)
	        				logRhogas(i)=logRhogas(i)+logRhogas0
	        				Pgas(i)=1.0*prefactor*T(i)*Rhogas(i)
	        				logPgas(i)=phlog10(Pgas(i))
	      				enddo
	    			else
	      				Rhogas0=0.0
            				Rhogas_central=Rhogas0
	      				logRhogas0=-20d0
	      				Rhogas(1:n)=0d0
	        			logRhogas(1:n)=-45d0
	        			Pgas(1:n)=0d0
	        			logPgas(1:n)=-45d0
	        			Yarray(1:n)=0d0
	        			logYarray(1:n)=-45d0
					return
	    			endif
	  		endif
	   	   	
!-----------------------------------------------------------
!GasModel=3
!
! GNFW-Planck model:

        ELSEIF(GasModel==3 )THEN
	    
!null run
		if( Ytot_Planck == 0.d0 ) then
			flag = 2
			return
		endif
	     
!	     WRITE(*,*)'ycoeff_Planck=  ', ycoeff_Planck
     
	     y0_Planck=2.d0 * ycoeff_Planck * thetas_Planck * y0_int_Planck    

             IF(a_GNFW==1.0620d0.AND. b_GNFW==5.4807d0 .AND. c_GNFW==0.3292d0 .AND. c500_GNFW== 1.156)THEN 

	        Y500_Planck=Ytot_Planck/1.814d0
		
		!write(*,*)'magicY500=', Y500_Planck
	     
	     ELSE
	     
	         Y500_Planck=(4.d0*pi)* ycoeff_Planck*(thetas_Planck*thetas_Planck)*&
	                                   CalcGNFW_PlancksphVolY(1.0d-2 ,theta500_Planck )
		
		!write(*,*)'nonmagicY500=', Y500_Planck			   
	    
	    ENDIF
	                
  
!             DO i=1,n
!	     yintegrand_Planck(i)=ycoeff_Planck*((theta_Planck(i)/thetas_Planck)**(-c_GNFW))*&	     
!	 ( (1.d0 + ((theta_Planck(i)/thetas_Planck)**(a_GNFW)))**((c_GNFW -b_GNFW)/a_GNFW))
	 
!	     logyintegrand_Planck(i)=phlog10(yintegrand_Planck(i))	

!	     WRITE(*,*)'i= ',i , 'theta_Planck(i)=', theta_Planck(i), 'yintegrand_Planck(i)= ' , yintegrand_Planck(i)	     
!	     ENDDO

!--------------------------------------------------------------------
! GasModel=5
!
!DM-GNFW Model:
! This model uses the GNFW profile for the electron pressure and dark matter density profile .
!---------------------------------------------------------------------

        ELSEIF(GasModel==5 )THEN

	      c200_DM=5.26d0*(((MT200_DM*h)/1.d14)**(-0.1d0))*(1.d0/(1.d0 +z(k))) !Neto et~al. 2007 for relaxed clusters			
	      Mg200_DM=MT200_DM*fg200_DM     !M_sun
!null run
		if( Mg200_DM == 0.d0 ) then
			flag = 2
			return
		endif
			      
	      r200_DM =((3.d0*MT200_DM)/(4.d0*pi*200.d0*rhocritz))**(1.d0/3.d0)   !Mpc
	      	      
	      rs_DM=r200_DM/c200_DM                   !Mpc
	      
	      rhos_DM=(200.d0/3.d0)*((r200_DM/rs_DM)**3.d0)*&
		   ( rhocritz/(DLOG(1.d0 + r200_DM/rs_DM) - ( 1.d0/(1.d0 + rs_DM/r200_DM) )))   !M_sunMpc-3
		   
	      !r500_DM=r200_DM/1.5d0                 !Mpc. The coeff 1.5 was derived from testing the model over a wide
	                                            ! range of cluster masses see paper for the derivation
	      !Correct calculation for r_500 following Hu Kravtsov implemented
	      !adapted for GM=5 kj 01/02/17
	      rsdivr200 = rs_DM/r200_DM
	      frsdivr200 = (rsdivr200)**(3.0)*(log(1 + (rsdivr200)**(-1.0)) - (1 + rsdivr200)**(-1.0))
	      f500 = 5.0*frsdivr200/2.0
	      p = a_2 + a_3*log(f500) + a_4 * (log(f500))**2.0
	      xf500 = (a_1 * (f500)**(2.0*p) + (9.0/16.0))**(-0.5) + 2*f500
	      r500_DM = rs_DM/xf500
		   
              rp_GNFW=r500_DM/c500_GNFW    !Mpc   
	     		       
               !DO i=1,n
               !    Pe_GNFW(i)=GNFWmodel3D(r(i),rp_GNFW,a_GNFW,b_GNFW,c_GNFW,1.0d0)
		  
	       !    logPe_GNFW(i)=phlog10(Pe_GNFW(i))
		   	       	       
               !ENDDO
	       
	      Pei_GNFW=	((mu_m/mu_e)*(G*rhos_DM*rs_DM*rs_DM*rs_DM)*Mg200_DM)/ &
	                 DM_GNFWgasVol(r200_DM,rs_DM,rp_GNFW,a_GNFW,b_GNFW,c_GNFW)	           
		      
	      MT500_DM=(4.d0*pi/3.d0)*(500.d0*rhocritz)*(r500_DM*r500_DM*r500_DM)        !M_sun
	                
	     
	      Mg500_DM= (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/rhos_DM)*(1.0d0/(rs_DM*rs_DM*rs_DM))* &
	                  DM_GNFWgasVol(r500_DM,rs_DM,rp_GNFW,a_GNFW,b_GNFW,c_GNFW)         !M_sun
			       
	      fg500_DM=Mg500_DM/MT500_DM
	     
	      Tg200_DM=(4.d0*pi*0.6*m_p*G*rhos_DM)*(rs_DM*rs_DM*rs_DM)*&
	    ((DLOG(1.0 +(r200_DM/rs_DM)) - (1.0/(1.0 +(rs_DM/r200_DM))))/r200_DM)*&           
        (1.0 +((r200_DM/rp_GNFW)**(a_GNFW)))*&
        (((b_GNFW*((r200_DM/rp_GNFW)**(a_GNFW))) +c_GNFW)**(-1.0))
	
	      Tg200_DM=Tg200_DM *(m_sun *Mpc2m*Mpc2m)*(J2keV)      !keV
	    
	      Tg500_DM=(4.d0*pi*0.6*m_p*G*rhos_DM)*(rs_DM*rs_DM*rs_DM)*&
	    ((DLOG(1.0 +(r500_DM/rs_DM)) - (1.0/(1.0 +(rs_DM/r500_DM))))/r500_DM)*&           
        (1.0 +((r500_DM/rp_GNFW)**(a_GNFW)))*&
        (((b_GNFW*((r500_DM/rp_GNFW)**(a_GNFW))) +c_GNFW)**(-1.0))
	
	      Tg500_DM=Tg500_DM *(m_sun *Mpc2m*Mpc2m)*(J2keV)       !keV 
	 
	      Pei_GNFW= Pei_GNFW*(m_sun/Mpc2m) *(J2keV)       !keVm-3
		   IF(Pei_GNFW.LE. 0d0)THEN
		      write(*,*)c200_DM ,MT200_DM, fg200_DM
		      write(*,*)Mg200_DM, r200_DM,rs_DM
		      write(*,*)rhos_DM, r500_DM, rp_GNFW
		      write(*,*)MT500_DM,Mg500_DM,fg500_DM 
		      write(*,*)Tg200_DM,Tg500_DM 
		      write(*,*)'error Pei_GNFW= ',Pei_GNFW
		      stop
		    ENDIF	 
	 
	
	 
	     y0_GNFW=((2.d0*sigma_T*Pei_GNFW*rp_GNFW)/(mec2*m2Mpc))*(1/a_GNFW)*&
	           (Gammafun((1.d0-c_GNFW)/a_GNFW) *&
		   Gammafun((b_GNFW-1.d0)/a_GNFW)/Gammafun((b_GNFW-c_GNFW)/a_GNFW))
		   	 
	     ycoeff_GNFW=(sigma_T*Pei_GNFW)/(mec2*m2Mpc)
	
	     Ysph500_GNFW=4.d0*pi*ycoeff_GNFW*CalcGNFWsphVolY(1.0d-3 ,r500_DM )

             Ysph200_GNFW=4.d0*pi*ycoeff_GNFW*CalcGNFWsphVolY(1.0d-3 ,r200_DM )	 
	     !calculated for input to PwS. Using same naming convention as gasmodel == 3
	     Ytot_Planck = 4.d0 * pi * ycoeff_GNFW * rp_GNFW**(3.d0) *&
                        gammaFun((3.d0 - c_GNFW) / a_GNFW) * gammaFun((b_GNFW - 3.d0) / a_GNFW)&
                        / (a_GNFW * gammaFun((b_GNFW - c_GNFW) / a_GNFW)) 
	     Ytot_Planck = Ytot_Planck /(D*D) *(rad2min*rad2min) !arcmin^2 !
             !n.b. theta_s in Planck notation is theta_p in ours!
	     thetas_Planck = rp_GNFW / D * rad2min

	     
		!uncomment for debugging purposes
	        !write(*,*)"NFW"
		!write(*,*)"MT200 = ",MT200_DM
		!write(*,*)"fg200 = ",fg200_DM
		!write(*,*)"a_GNFW = ", a_GNFW
		!write(*,*)"b_GNFW = ", b_GNFW 
		!write(*,*)"c_GNFW = ", c_GNFW
		!write(*,*)"c500_GNFW = ", c500_GNFW 		
		!write(*,*)"c200 = ",c200_DM
		!write(*,*)"Mg200 = ", Mg200_DM
		!write(*,*)"r200 = ", r200_DM
		!write(*,*)"rs = ", rs_DM
		!write(*,*)"rhos = ", rhos_DM
		!write(*,*)"r500 = ", r500_DM
		!write(*,*)'Pei_GNFW= ',Pei_GNFW
		!write(*,*)"Mg500 = ", Mg500_DM
		!write(*,*)"Tg200 = ", Tg200_DM
		!write(*,*) "Tg500 = ", Tg500_DM 
		!write(*,*)"y_0 = ", y0_GNFW
		!write(*,*)"Y_200 = ",Ysph200_GNFW
		!write(*,*)"Y_500 = ",Ysph500_GNFW	  
	  	!DO i=1,n
		!yintegrand_GNFW(i)=ycoeff_GNFW *Pe_GNFW(i) 
	        	  
	   	!logyintegrand_GNFW(i)= phlog10(yintegrand_GNFW(i))
		!ENDDO	       			
	

	ELSEIF(GasModel==6 )THEN !adapted for GM=6 kj 01/02/17

		c200_DM = 10.d0 ** (0.459d0 + 0.518d0 * exp(-0.49d0 * z(k)**(1.303d0)) &
                        + (-0.13d0 + 0.029d0 * z(k)) * log10(MT200_DM * h / 1.d12))!Dutton et al. 2014 simulations fitted with Einasto profile			
	        Mg200_DM=MT200_DM*fg200_DM     !M_sun
	
		r200_DM =((3.d0*MT200_DM)/(4.d0*pi*200.d0*rhocritz))**(1.d0/3.d0)   !Mpc

		rm2_DM = r200_DM/c200_DM !MPc, analgous to r200_DM from GM 5

		rhom2_DM = (200.d0 / 3.d0) * (r200_DM / rm2_DM) ** (3.d0) * & 
			   rhocritz / (1.d0 / aEin_DM * (aEin_DM / 2.d0)**(3.d0/aEin_DM) * &
		           exp(2.d0 / aEin_DM) & 
			   * gammaFun(3.d0 / aEin_DM) * gammds(2.d0 / aEin_DM * (r200_DM / rm2_DM)**aEin_DM, 3.d0 / aEin_DM, flag))

		if(flag == 1 .or. flag == 2) then
			write(*,*) "Gamma calculation for rhom2_DM has screwed up"
			flag = 1
			return
		endif

		X0 = r200_DM / (1.5d0 * rm2_DM) !Initial guess for r500 for Newton Raphson
		call find_root( X_func, X0, 1.0d-4, 50, X, flag )
		if(flag == 1) then
			write(*,*) "Newton Raphson to determine r500 failed to converge."
			return
		endif
		r500_DM = X * rm2_DM
		rp_GNFW = r500_DM / c500_GNFW
	  
		!do i=1,n
		!           Pe_GNFW(i) = GNFWmodel3D(r(i),rp_GNFW,a_GNFW,b_GNFW,c_GNFW,1.0d0) 
		!	  
		!	   logPe_GNFW(i) = phlog10(Pe_GNFW(i))
			   	       	       
		!enddo 

		
		Pei_GNFW = ((mu_m/mu_e)*(G*rhom2_DM*rm2_DM*rm2_DM*rm2_DM)*Mg200_DM)* &
			   1.d0 / aEin_DM * (aEin_DM / 2.d0) ** (3.d0 / aEin_DM) * exp(2.d0 / aEin_DM) / &
			   Ein_GNFWgasVol(a_GNFW, b_GNFW, c_GNFW, aEin_DM, r200_DM, rm2_DM, rhom2_DM, rp_GNFW)
		MT500_DM=(4.d0*pi/3.d0)*(500.d0*rhocritz)*(r500_DM*r500_DM*r500_DM)        !M_sun
			        
		     
		Mg500_DM= (mu_e/mu_m)*(1.d0/G)*(Pei_GNFW/rhom2_DM)*(1.0d0 / &
			  (1.d0 / aEin_DM * (aEin_DM / 2.d0) ** (3.d0 / aEin_DM) * exp(2.d0 / aEin_DM) * rm2_DM*rm2_DM*rm2_DM))* &
			  Ein_GNFWgasVol(a_GNFW, b_GNFW, c_GNFW, aEin_DM, r500_DM, rm2_DM, rhom2_DM, rp_GNFW)        !M_sun
				       
		fg500_DM=Mg500_DM/MT500_DM

		Tg200_DM= 4.d0 * pi * 0.6d0 * m_p * G * rhom2_DM * rm2_DM ** (3.d0) * 1.d0 / aEin_DM * &
			  (aEin_DM / 2.d0) ** (3.d0 / aEin_DM) * exp(2.d0 / aEin_DM) *&
			  gammaFun(3.d0 / aEin_DM) * gammds(2.d0 / aEin_DM * (r200_DM / rm2_DM) ** (aEin_DM), 3.d0 / aEin_DM, flag) / &
			  r200_DM * (1.d0 + (r200_DM / rp_GNFW) ** a_GNFW) * &
			  (b_GNFW * (r200_DM / rp_GNFW) ** a_GNFW + c_GNFW) **(-1.d0)

		if(flag == 1 .or. flag == 2) then
			flag = 1
			write(*,*) "Gamma calculation for Tg200_DM has screwed up"
			return
		endif

		Tg200_DM=Tg200_DM *(m_sun *Mpc2m*Mpc2m)*(J2keV)      !keV

		Tg500_DM= 4.d0 * pi * 0.6d0 * m_p * G * rhom2_DM * rm2_DM ** (3.d0) * 1.d0 / aEin_DM * &
			  (aEin_DM / 2.d0) ** (3.d0 / aEin_DM) * exp(2.d0 / aEin_DM) *&
			  gammaFun(3.d0 / aEin_DM) * gammds(2.d0 / aEin_DM * (r500_DM / rm2_DM) ** (aEin_DM), 3.d0 / aEin_DM, flag) / &
			  r500_DM * (1.d0 + (r500_DM / rp_GNFW) ** a_GNFW) * &
			  (b_GNFW * (r500_DM / rp_GNFW) ** a_GNFW + c_GNFW) **(-1.d0)

		if(flag == 1 .or. flag == 2) then
			flag = 1
			write(*,*) "Gamma calculation for Tg500_DM has screwed up"
			return
		endif

		Tg500_DM=Tg500_DM *(m_sun *Mpc2m*Mpc2m)*(J2keV)      !keV
		Pei_GNFW= Pei_GNFW*(m_sun/Mpc2m) *(J2keV)       !keVm-3
			   if(Pei_GNFW <= 0.d0) then
			      write(*,*)c200_DM ,MT200_DM, fg200_DM
			      write(*,*)Mg200_DM, r200_DM
			      write(*,*)rm2_DM, r500_DM, rp_GNFW
			      write(*,*)MT500_DM,Mg500_DM,fg500_DM 
			      write(*,*)Tg200_DM,Tg500_DM 
			      write(*,*)'error Pei_GNFW= ',Pei_GNFW
			      stop
			    endif

		y0_GNFW=((2.d0*sigma_T*Pei_GNFW*rp_GNFW)/(mec2*m2Mpc))*(1/a_GNFW)*&
	           (Gammafun((1.d0-c_GNFW)/a_GNFW) *&
		   Gammafun((b_GNFW-1.d0)/a_GNFW)/Gammafun((b_GNFW-c_GNFW)/a_GNFW))
		   	 
		ycoeff_GNFW=(sigma_T*Pei_GNFW)/(mec2*m2Mpc)

		Ysph500_GNFW=4.d0*pi*ycoeff_GNFW*CalcGNFWsphVolY(1.0d-3 ,r500_DM )
		  
		Ysph200_GNFW=4.d0*pi*ycoeff_GNFW*CalcGNFWsphVolY(1.0d-3 ,r200_DM )
             !calculated for input to PwS. Using same naming convention as gasmodel == 3
	     Ytot_Planck = 4.d0 * pi * ycoeff_GNFW * rp_GNFW**(3.d0) *&
                        gammaFun((3.d0 - c_GNFW) / a_GNFW) * gammaFun((b_GNFW - 3.d0) / a_GNFW)&
                        / (a_GNFW * gammaFun((b_GNFW - c_GNFW) / a_GNFW)) 
             Ytot_Planck = Ytot_Planck /(D*D) *(rad2min*rad2min) !arcmin^2
             !n.b. theta_s in Planck notation is theta_p in ours!
	     thetas_Planck = rp_GNFW / D  * rad2min	   
		!below can be uncommented for debugging purposes
		!write(*,*)"Einasto"
		!write(*,*)"MT200 = ",MT200_DM
		!write(*,*)"fg200 = ",fg200_DM
		!write(*,*)"a_GNFW = ", a_GNFW
		!write(*,*)"b_GNFW = ", b_GNFW 
		!write(*,*)"c_GNFW = ", c_GNFW
		!write(*,*)"c500_GNFW = ", c500_GNFW 
		!write(*,*)"aEin_DM = ", aEin_DM		
		!write(*,*)"c200 = ",c200_DM
		!write(*,*)"Mg200 = ", Mg200_DM
		!write(*,*)"r200 = ", r200_DM
		!write(*,*)"rm2 = ", rm2_DM
		!write(*,*)"rho_m2 = ", rhom2_DM
		!write(*,*)"r500 = ", r500_DM
		!write(*,*)'Pei_GNFW= ',Pei_GNFW
		!write(*,*)"Mg500 = ", Mg500_DM
		!write(*,*)"Tg200 = ", Tg200_DM
		!write(*,*) "Tg500 = ", Tg500_DM 
		!write(*,*)"y_0 = ", y0_GNFW
		!write(*,*)"Y_200 = ",Ysph200_GNFW
		!write(*,*)"Y_500 = ",Ysph500_GNFW
		!stop	   

		!do i=1,n
		!	yintegrand_GNFW(i)=ycoeff_GNFW *Pe_GNFW(i) 
		!        	  
		!   	logyintegrand_GNFW(i)= phlog10(yintegrand_GNFW(i))
		!enddo	      			
	endif	
     
!---------------------------------------------------------	
! 	Finally project pressure to generate y arrays:

	if(GasModel == 1 .OR. GasModel==2) then
		do i=1,n
	    		uu=r(i)
          		rlimit1 = sqrt( max(rlimit * rlimit - uu * uu, 0d0) )
          
	    		if( rlimit1 >0d0 ) call qtrap(PgasIntegrand,-rlimit1,rlimit1,eps,Yarray(i))
	    		if(Yarray(i) > 1d10) then
	      			write(*,*) 'i,r(i),rlimit1,Yarray(i)=',i,r(i),rlimit1,Yarray(i)
	      			write(*,*) 'QTRAP error. Printing pressure/density arrays:'
	      			do j=1,n
	          			write(*,*)j,r(j),Pgas(j),logPgas(j),logRhogas(j),Rhogas(j)
		  		enddo
	      			write(*,*) 'rhogas0=',rhogas0
	      			write(*,*) 'prefactor=',prefactor
	      			write(*,*) 'r200,rc,rlimit=',r200,rc,rlimit
	      			write(*,*) 'GasPars(1)=',GasPars(k,1)
	      			write(*,*) 'TPars(1)=',TPars(k,1)
	      			write(*,*) 'MassPars(1)=',MassPars(k,1)
	      			write(*,*) 'MassPars(2)=',MassPars(k,2)
	    		endif
	    		logYarray(i)=phlog10(Yarray(i))
	  	enddo 
	ELSEIF(GasModel==3) THEN
	   DO i=1,n
	      uu=theta_Planck(i)
	      rlimit1=sqrt( max(thetalimit_Planck * thetalimit_Planck - uu * uu, 0d0) )
              r1=sqrt( max(thetamin_Planck*thetamin_Planck-uu*uu, 0d0) )
              r2=sqrt( max(thetamax_Planck*thetamax_Planck-uu*uu, 0d0) )
	      if( rlimit1 >0d0 ) then
                yarray_Planck(i)=qtrap2(yintegral_Planck, r1, r2)
                if (r1.gt.0d0) yarray_Planck(i)=yarray_Planck(i)+GNFWmodel3D(r1,thetas_Planck,a_GNFW,b_GNFW,c_GNFW,ycoeff_Planck)*r1
              endif
	      if(yarray_Planck(i) > 1d10 .or. yarray_Planck(i)/=yarray_Planck(i)) then
	          write(*,*) 'i,theta_Planck(i),rlimit1,yarray_Planck(i)=',i,theta_Planck(i),rlimit1,yarray_Planck(i)		   
	      	  write(*,*) 'QTRAP2 error.'			
	      endif
!    Take advantage of symmetry, quicker?
              yarray_Planck(i)=2.d0*yarray_Planck(i)
              logyarray_Planck(i)=phlog10(yarray_Planck(i))
	   ENDDO
	   				
	ELSEIF(GasModel==5 .or. GasModel==6) THEN !adapted for GM=6 kj 01/02/17
	   DO i=1,n	
	      uu=r(i)
              ! Don't need to integrate any closer than the cell size to fill up the ymap
              if (uu/(sec2rad*D).gt.cell/4) then
          	rlimit1 = sqrt( max(rlimit * rlimit - uu * uu, 0d0) )
                r1=sqrt( max(rmin*rmin-uu*uu, 0d0) )
                r2=sqrt( max(rmax*rmax-uu*uu, 0d0) )
	    	if( rlimit1 >0d0 ) then
                    yarray_GNFW(i)=qtrap2(GNFWyintegrand, r1, r2)
                    if (r1.gt.0d0) yarray_GNFW(i)=yarray_GNFW(i)+GNFWmodel3D(r1,rp_GNFW,a_GNFW,b_GNFW,c_GNFW,ycoeff_GNFW)*r1
                    yarray_GNFW(i)=yarray_GNFW(i)*2d0
                endif
	        if(yarray_GNFW(i) > 1d10 .or. yarray_GNFW(i)/=yarray_GNFW(i)) then
	          write(*,*) 'i,r(i),rlimit1,yarray_GNFW(i)=',i,r(i),rlimit1,yarray_GNFW(i)		   
	      	  write(*,*) 'QTRAP2 error.'			
	        endif
                logyarray_GNFW(i)=phlog10(yarray_GNFW(i))
              endif
	  ENDDO
	
	endif
	flag=0
      
      	!store derived parameters
	IF(GasModel ==1)THEN	
      	aux(k,1)=rc/(sec2rad*D) !r_{core} in arcsec
      	aux(k,2)=rhogas_central !central gas density in [h Msun Mpc^{-3}]
      	aux(k,3)=rhogas_central/(1.14d0*m_p*(Mpc2m**3)) !central electron number density in [m^{-3}]
      	aux(k,4)=totMass(1500d0,r200,rc,rhocritz,aux(k,22))
      	aux(k,5)=totMass(1000d0,r200,rc,rhocritz,aux(k,23))
      	aux(k,6)=totMass(500d0,r200,rc,rhocritz,aux(k,24))
      	aux(k,7)=totMass(200d0,r200,rc,rhocritz,aux(k,25))
      	aux(k,8)=totMass(178d0,r200,rc,rhocritz,aux(k,26))
      	aux(k,9)=totMass(150d0,r200,rc,rhocritz,aux(k,27))
      	call getGassMass(rc,beta,rhogas_central,r200,aux(k,10:14))
      	aux(k,15)=aux(k,14)
      	aux(k,14)=aux(k,13)
      	aux(k,13)=Mgas200
      	do i=16,21
      		if(aux(k,i-6)==0d0 .and. aux(k,i-12)==0d0) then
            		aux(k,i)=0d0
		else
            		aux(k,i)=aux(k,i-6)/aux(k,i-12)
		endif
      		if(.not.(aux(k,i)/h>=0d0 .and. aux(k,i)/h<=1d0)) then
            		flag=1
			return
		endif
	enddo
      	aux(k,28)=rhocritz
        aux(k,29)=T200 	
	ELSEIF(GasModel==3)THEN



	aux(k,1)=ycoeff_Planck	
	aux(k,2)= y0_Planck  
	aux(k,3)=5.0*theta500_Planck           !arcmin	    
	aux(k,4)=Ytot_Planck                   !arcmin2
	aux(k,5)=theta500_Planck               !arcmin		
	aux(k,6)=Y500_Planck                   !arcmin2
	
	ELSEIF(GasModel==5 .or. GasModel==6)THEN	!adapted for GM=6 kj 01/02/17
	
		!aux(k,1)=alpha_L1 !testing                                      
		!aux(k,2)=alpha_L2 !purposes
		aux(k,1)=D !angular diameter distance Mpc
		aux(k,2)=c200_DM
		aux(k,3)= GasPars(k,1)
		aux(k,4) =r200_DM	
		aux(k,5) =Mg200_DM
		aux(k,6) =Tg200_DM
	      	aux(k,7) =MT500_DM
		aux(k,8) =fg500_DM
	      	aux(k,9) =r500_DM
	   	aux(k,10) =Mg500_DM
	   	aux(k,11) =Tg500_DM	 	 	 		
	      	aux(k,12)=y0_GNFW
	      	aux(k,13)=Ysph200_GNFW
		aux(k,14)=Ysph500_GNFW
		aux(k,15)=(r200_DM/D)*rad2min
		aux(k,16)=(Ysph200_GNFW/(D*D))*(rad2min*rad2min)
		aux(k,17)=(r500_DM/D)*rad2min
		aux(k,18)=(Ysph500_GNFW/(D*D))*(rad2min*rad2min)
		aux(k,19)=Pei_GNFW
		if(GasModel==5) then !adapted for GM=6 kj 01/02/17
			aux(k,20)=rhos_DM
			aux(k,21)=rs_DM
		else
			aux(k,20)=rhom2_DM
			aux(k,21)=rm2_DM
		endif
		aux(k,22)=rp_GNFW 
		aux(k,23)=(rp_GNFW/D)*rad2min !thetas_Planck/ thetap_GNFW
		aux(k,24)=rhocritz
		!25th aux param is Y_cycl which is added below
		aux(k,26)=Ytot_Planck                   !arcmin2	
	
	ENDIF
                
 999	return 
	end subroutine MakeGasDistributions

!=======================================================================

	subroutine MakeYMap(k)

! Generate a map of comptonisation parameter y for a given set of
! cluster parameters corresponding to the kth atom

	implicit none
	
	integer i,j,k	
	double precision x0,y0
	double precision xx,yy,dx(2),rr,angfactor
	double precision A1,rp1,Gm
	parameter(A1=1.424d-19,rp1=1.0)
        parameter(Gm=3.12d-14)
!-----------------------------------------------------------------------

	x0=GeoPars(k,1)
	y0=GeoPars(k,2)

!------------------------------------------------------------------------
      
	if(GasModel==0 .or.GasModel==3 ) then
		angfactor=1.0/60.0
	elseif(GasModel==1 .or. GasModel==2 .or. GasModel==5 .or. GasModel==6 )THEN !adapted for GM=6 kj 01/02/17
        	angfactor=sec2rad*D
     	elseif(GasModel==4)THEN
                angfactor=1.0		
      	endif
	IF(GasModel==4)THEN
	   	do i=1,nx
			do j=1,ny
			
		 		xx=trans(1)+i*trans(2)
		 		yy=trans(4)+j*trans(6)
		 		dx(1)=(xx-x0)
		 		dx(2)=(yy-y0)
				
				IF(GeoModel==1)THEN
		 			rr=(dx(1)*dx(1)+dx(2)*dx(2))*angfactor*angfactor
				ELSEIF( GeoModel==2)THEN
					call QuadForm(2,dx,Q,dx,rr)
		 			rr=rr*angfactor*angfactor
				ENDIF 
				
				ymap(i,j)=dT0_BA*((1.d0+rr/( thetac_BA*thetac_BA))**((1.d0 -3.d0*Beta_BA)/2.d0))  !temperature decrement map
				 
                        enddo
		enddo	 
        	RETURN
   
   	ELSE             
!-----------------------------------
	map_sum_GNFW=0d0
		do i=1,nx
			do j=1,ny
		 		xx=trans(1)+i*trans(2)
		 		yy=trans(4)+j*trans(6)
		 		dx(1)=(xx-x0)
		 		dx(2)=(yy-y0)
				if(GeoModel==1) then ! Spherical geometry:
		 			rr=sqrt(dx(1)*dx(1)+dx(2)*dx(2))*angfactor
				else if(GeoModel==2) then ! Elliptical geometry:
					call QuadForm(2,dx,Q,dx,rr)
		 			rr=sqrt(rr)*angfactor
				endif
			     
			     	IF(GasModel==1 .OR. GasModel==2)THEN				
		 			ymap(i,j)=ymap(i,j)+yfunc(rr)
			     	ELSEIF(GasModel==3)THEN
                                        !!!!! THIS IS A BUG, MAP NOT CENTRED ON CLUSTER
					!IF( i==nx/2+1 .AND.j==ny/2+1)THEN
                                        IF (abs(dx(1)).lt.cell/2.0 .and. abs(dx(2)).lt.cell/2.0) THEN
				 		ymap(i,j)=y0_Planck
					ELSE
		 				ymap(i,j)=ymap(i,j)+yfunc_Planck(rr)
					ENDIF
					
			      ELSEIF(GasModel==5 .or. GasModel==6)THEN !adapted for GM=6 kj 01/02/17
		                     !IF( i==nx/2+1 .AND.j==ny/2+1)THEN
                                     IF (abs(dx(1)).lt.cell/2.0 .and. abs(dx(2)).lt.cell/2.0) THEN
		                         ymap(i,j)=y0_GNFW
		                     ELSE
		                         ymap(i,j)=ymap(i,j)+yGNFWfunc2(rr)
		                    ENDIF
					
				ENDIF
              map_sum_GNFW = map_sum_GNFW + ymap(i,j)							     				
	      		enddo
	   	enddo

		! Calculating volume integrated y parameter both in cylindrical and sphericalgeometry at r500 and r200.
     		IF(GasModel==1)THEN     
            		IF(TStyle(k)==0 .OR. TStyle(k)==2) THEN      !Assuming HSE at r500 and r200 and const. temperature throughout the cluster
	       			r500=aux(k,24)
	       		ELSEIF(TStyle(k) == 1) THEN
	        		r500=max((9.d0/(4.d0*pi*Gm*500.d0*rhocritz))*beta*T200-rc*rc,0.d0)
             			r500=sqrt(r500)  !Assuming HSE at r500 and virial at r200 and const. temperature throughout the cluster(T500=T200)
	    		ENDIF
			
			!cylindrical	    
	    
	    		IF(beta.NE.1.) THEN
               
               			Ycyl500=(pi*ymap(nx/2,ny/2)*rc*rc) *(1.0/(1.5-1.5*beta)) * &
                        	((1.0+(r500*r500)/(rc*rc))**(1.5 -1.5*beta) -1.0)

               			Ycyl200=(pi*ymap(nx/2,ny/2)*rc*rc) *(1.0/(1.5-1.5*beta)) * &
                        	((1.0+(r200*r200) / (rc*rc))**(1.5 -1.5*beta) -1.0)

            		ELSEIF(beta.EQ.1.) THEN
              
				Ycyl500=(pi*ymap(nx/2,ny/2)*rc*rc)*DLOG(1.0 +(r500*r500) /(rc*rc))
              			Ycyl200=(pi*ymap(nx/2,ny/2)*rc*rc)*DLOG(1.0 +(r200*r200) /(rc*rc))

			ENDIF

			Ysph500= A1*T200*Rhogas0* CalcBETAsphVolY(1.0d-4 ,r500 )
			Ysph200= A1*T200*Rhogas0* CalcBETAsphVolY(1.0d-4 ,r200 )
         	     
			aux(k,30)=Ycyl500
			aux(k,31)=Ycyl200
      			aux(k,32)=Ysph500
			aux(k,33)=Ysph200
		ELSEIF(GasModel==5 .or. GasModel==6) THEN !adapted for GM=6 kj 01/02/17
		Ycyl_GNFW = D*D* map_sum_GNFW * cell * cell * sec2rad * sec2rad
                aux(k,25)=  Ycyl_GNFW		
 		ENDIF

	ENDIF      
	
	end subroutine MakeYMap
		
!=======================================================================

	subroutine MakeVisibilities

! Convert a map of comptonisation parameter y to microwave sky
! brightness, fourier transform and sample at the relevant uv points.
! Might well be worth investigating more efficient fourier transform
! routines, if dataset is small.

	implicit none
	
	integer i,j,m
	double precision szsky(nx,ny),re1(nx,ny),im(nx,ny)
	double precision xx,yy,r2,beam,taper,fu(large),fv(large),rr,ii
	double precision ll,mm,theta      
      integer k
      logical verbose1
!-----------------------------------------------------------------------

      ! Set to true to write *lots* of output...
      verbose1 = .false.	

! Make sure map has zero mean, and convert to Janskys:
      
!      sum=0.0
!      do i=1,nx
!         do j=1,ny
!            sum=sum+ymap(i,j)
!         end do
!      end do
!      sum=sum/float(nx*ny)	
	do i=1,nx
         do j=1,ny
            !szsky(i,j)=(ymap(i,j)-sum)*SZscale
            !!FF Hack
            !szsky(i,j)=(ymap(i,j)-sum)
            szsky(i,j)=(ymap(i,j))
         end do
      end do

! Loop over data files:

      do m=1,Nvisfiles
	  if(verbose1) write(*,*) 'for pointing ',m,'...'

!  Apply primary beam 
          beam=pb_sig(m)
          lmax1=lmax(m)
          lmin1=lmin(m)
        if(verbose1) then
        write(*,*) ' telescope=',trim(telescope(m))
        write(*,*) ' SZscale=',SZscale(m)
        write(*,*) ' beam size=',beam
        write(*,*) ' map centre (deg)=',map_x,map_y
        write(*,*) ' beam centre (deg)=',pb_x(m),pb_y(m)
        write(*,*) ' phase centre (deg)=',ph_x(m),ph_y(m)
        write(*,*) ' beam offset (sec)=',pb_x_offset(m),pb_y_offset(m)
        write(*,*) ' phase offset (sec)= ',ph_x_offset(m),ph_y_offset(m)
        write(*,*) ' model offset (sec)=',GeoPars(1,1),GeoPars(1,2)
        endif
        beam=2.0*beam*beam
!write(file,'(i,a)') m,'.ymap'
!open(2,file=file,form='formatted',status='replace')

     
        do j=1,ny
           do i=1,nx
           
!write(2,'(2i6,E20.6)')i,j,szsky(i,j)*SZscale(m)/ &
!	((cell*sec2rad)*(cell*sec2rad))* &
!	((clight/nu(1))**2.)*(10.**-26)/(kboltzmann*2.)

             xx=trans(1)+i*trans(2)-pb_x_offset(m)
             yy=trans(4)+j*trans(6)-pb_y_offset(m)
             r2=xx*xx+yy*yy
             taper=exp(-r2/beam)
IF(GasModel==4)THEN
              re1(i,j)=taper*szsky(i,j)*SZscale_BA(m)
              im(i,j)=0.0
ELSE	      
             re1(i,j)=taper*szsky(i,j)*SZscale(m)
             im(i,j)=0.0

             if(verbose1) then
             if(i==nx/2.and.j==ny/2) then
               write(*,*) 'szsky central pixel:',szsky(nx/2,ny/2)*SZscale(m)
               write(*,*) 'taper at this point=',taper
               write(*,*) 're map central pixel:',re1(nx/2,ny/2)
             endif
             endif
ENDIF
           enddo
        enddo
! close(2)
        
! FFT into aperture plane:
      
! Old code (profile) used NR FFT, painfully slow.
!
!        call xzapit(re,nx,ny)
!         nn(1)=nx
!         nn(2)=ny
!         k=1
!         do j=1,ny
!            do i=1,nx
!                work(k)=re(i,j)
!                k=k+1
!                work(k)=im(i,j)
!                k=k+1
!            enddo
!         enddo
!         call fourn2(work,nn,2,1)
!         k=1
!         do j=1,ny
!            do i=1,nx
!               re(i,j)=work(k)
!               k=k+1
!               im(i,j)=work(k)
!               k=k+1
!            end do
!         end do
!         call xzapit(re,nx,ny)
!         call xzapit(im,nx,ny)
 	
! New FFT code: Hobson's rfft.f code used, but only if the FFTW library
! is not available.

!      if(FFTW) then

        do j=1,ny
          do i=1,nx
            arr(i,j)=dcmplx(re1(i,j),im(i,j))
          end do
        end do
        call cheqboard(arr,nx,ny)
        call dfftw_execute(fftwplan)
        call cheqboard(arr,nx,ny)
        do j=1,ny
          do i=1,nx
            re1(i,j)=dble(arr(i,j))
            im(i,j)=aimag(arr(i,j))
          end do
        end do
!
!      else
!
!        nd=2
!        nn(1)=nx
!        nn(2)=ny
!        k=0
!        do j=1,ny
!          do i=1,nx
!            k=k+1
!            wkspce(k)=dcmplx(re1(i,j),im(i,j))
!          end do
!        end do
!        job=+1                   ! forward transform (-ve exponential)
!        iform=0                  ! data are real
!        call cheqboard(wkspce,nx,ny)
!        call fourt(wkspce,nn,nd,job,iform,work)
!        call cheqboard(wkspce,nx,ny)
!        k=0
!        do j=1,ny
!          do i=1,nx
!            k=k+1
!            re1(i,j)=dble(wkspce(k))
!            im(i,j)=aimag(wkspce(k))
!          end do
!        end do

!      endif
      
!  Now sample aperture plane at given u-v points, and rotate to phase
!  centre:
	
        ll=sin(-ph_x_offset(m)*sec2rad)*cos(-ph_y_offset(m)*sec2rad)
        mm=sin(-ph_y_offset(m)*sec2rad)

	
!      write(file,'(i,a)') m,'.vis'
!	open(2,file=file,form='formatted',status='replace')

        do k=1,Nvis(m)

          if(SZeflag(m,k)==0) then

            fu(k)=-1.0*u(m,k)
            fv(k)=1.0*v(m,k)
            
            call extract_visibility(re1,im,nx,ny,fu(k),fv(k),cell,rr,ii)
		
            theta=1.0*TwoPi*(-u(m,k)*ll+v(m,k)*mm)
            
!           Phase rotation by theta:
!           (vr+i*vi)*(cos+i*sin)=(vr*cos-vi*sin)+i*(vr*sin+vi*cos)

            pvisr(m,k)=dble(rr)*dble(cos(theta))-dble(ii)*dble(sin(theta))
            pvisi(m,k)=dble(rr)*dble(sin(theta))+dble(ii)*dble(cos(theta))
            
!            write(2,'(i4,2F10.3,2E14.4,F4.1)')k,u(m,k),v(m,k),pvisr(m,k),pvisi(m,k),0.0

          else
            pvisr(m,k)=0.d0
            pvisi(m,k)=0.d0
          endif

        enddo
!      close(2)
	
      enddo
!      stop

      ! Store cluster-only visibilities so that you can add different sets of point sources to them when implementing fast/slow parameters in PolyChord
      pvisr1=pvisr
      pvisi1=pvisi
      
	return
	end subroutine MakeVisibilities
	
!=======================================================================

	subroutine AddSources

! Add a number of point sources to the predicted visibilities.

	implicit none
	
	integer i,j,m
	double precision ll,mm,ff,r2,beam,taper,ra,dec,sra,sdec,ra_off,dec_off,flux
	double precision theta,lambda
	
!-----------------------------------------------------------------------

! Initialise visibilities to cluster-only visibilities
       pvisr = pvisr1
       pvisi = pvisi1

! Loop over pointings

      do m=1,Nvisfiles

       beam=pb_sig(m)
       beam=2.d0*beam*beam

       !lambda=clight/CRVAL4
       lambda=clight/nu(m)

! Loop over sources, and visibilities, calculating amplitude and phase
! factors at each uv point and applying these to the predicted data:
! Cannot use small angle approximation for VSA fields!

       do i=1,NSrc

!       Apply primary beam to source flux:

         ra=pb_x(m)*deg2rad
         dec=pb_y(m)*deg2rad
         sra=SrcPars(i,1)*deg2rad
         sdec=SrcPars(i,2)*deg2rad

         call calc_sepn(ra,dec,sra,sdec,r2)
         r2=r2*rad2sec
         r2=r2*r2
         taper=exp(-r2/beam)
         !ff=SrcPars(i,3)*taper
         !calculate the source flux at the channel frequency using spectral index
         flux=flux0(i)*((nu(m)/nu0(i))**(-SrcPars(i,4)))
         ff=flux*taper

!        Find direction cosines of source position vector
!        For small angles, ll \approx x/radians etc:
!        ll=SrcPars(i,1)*sec2rad
!        mm=SrcPars(i,2)*sec2rad
!        Otherwise, do offsets properly and use full formula:
!        (Note that it is relative to the phase centre we calculate)

         ra=ph_x(m)*deg2rad
         dec=ph_y(m)*deg2rad
         call calc_ra_dec_off(ra,dec,sra,sdec,ra_off,dec_off)
         ra_off=-1.0*ra_off

         ll=sin(ra_off)*cos(dec_off)
         mm=sin(dec_off)
         
!        Finally add the point source flux directly to all visibilities

         do j=1,NVis(m)
           theta=1.0*TwoPi*(-u(m,j)*ll+v(m,j)*mm)
           pvisr(m,j)=pvisr(m,j)+ff*cos(theta)
           pvisi(m,j)=pvisi(m,j)+ff*sin(theta)
         enddo

       enddo

      enddo 
       
	return
	end subroutine AddSources
	
!-----------------------------------------------------------------------

      subroutine extract_visibility(re_map,im_map,nx,ny,u,v,cell,re1,im)

! De-grids a visibility at specified u,v from an unconvolved gridded 
! aperture

      implicit none

	integer nx,ny
      double precision re_map(nx,ny),im_map(nx,ny)
      double precision u,v,cell,re1,im,uvdist,l
      double precision sec2rad
      parameter(sec2rad=4.8481368d-6)
      double precision centre,x,y,uvcell

      uvcell=1.0/(nx*cell*sec2rad)
      centre=(nx/2)+1
      x=(u/uvcell)+centre
      y=(v/uvcell)+centre

      if((x<0).or.(y<0).or.(x.gt.nx).or.(y.gt.ny)) then
         re1=0.0
         im=0.0
         write(*,*) 'extract_visibility: aperture plane overshoot.'
         write(*,*) x,y,nx,ny
      else
         uvdist=sqrt(u*u+v*v)
         l=2.0*pi*uvdist
         if((l>lmax1 .or. l<lmin1) .and. .false.) then
         	re1=0.
            im=0.
	   else
         	call interp2d(re_map,nx,ny,x,y,re1)
         	call interp2d(im_map,nx,ny,x,y,im)
  	   endif
      end if

      return
	end subroutine extract_visibility

!======================================================================


	function BetaModel2D(r,rc,beta,y0)
	
	implicit none
	
	double precision r,rc,beta,y0
	double precision BetaModel2D
	
	BetaModel2D=y0/((1.0+(r/rc)*(r/rc))**(3.0*beta/2.0-0.5))

	return
	end function BetaModel2D
      
!=======================================================================

	function BetaModelHM(r,rc,y0)
	
	implicit none
	
	double precision r,rc,y0,rt,bracket
	double precision BetaModelHM
	
	rt=3.0*rc
      bracket=1.0/sqrt(rc*rc+r*r)	
      bracket=bracket-1.0/sqrt(rt*rt+r*r)	
      BetaModelHM=bracket*y0*rc*rt/(rt-rc)
       
	return
	end function BetaModelHM
      
!=======================================================================

	function BetaModel3D(r,rc,beta,rhogas0)
	
	implicit none
	
	double precision r,rc,beta,rhogas0
	double precision BetaModel3D
	
	BetaModel3D=rhogas0/((1.0+(r/rc)*(r/rc))**(3.0*beta/2.0))
 
	return
	end function BetaModel3D
	
!=======================================================================

      function BetaMass(rc,beta,y0,rlimit1)
      
      implicit none
      
      double precision rc,beta,y0,rlimit1,result,BetaMass
      double precision eps
      parameter(eps=1d-4)
      
      a=rc
      b=beta
      c=y0
      call qtrap(Mgas_Beta_Integrand,0.d0,rlimit1,eps,result)
      
      BetaMass=result
      
      return
      end function BetaMass

!=======================================================================

      function GasMass(rlimit1,rc,beta)
      
      implicit none
      
      double precision rlimit1,GasMass,result,rc,beta
	double precision eps
	parameter(eps=1d-4)
	
      gd1=rc
      gd2=beta

      !call qtrap(Mgas_Integrand,0.d0,rlimit1,eps,result)
      call qtrap(betaProfileInt,0.d0,rlimit1,eps,result)
      GasMass=result
      
      return
      end function GasMass

!=======================================================================

      function Mgas_Beta_Integrand(r)
      
      implicit none
      
      double precision r,Mgas_Beta_Integrand
      double precision TwoPi
      parameter(TwoPi=6.283185307)
      
      Mgas_Beta_Integrand=TwoPi*r*c/((1.0+(r/a)*(r/a))**(3.0*b/2.0-0.5))
      
      return
      end function Mgas_Beta_Integrand

!=======================================================================

      function Mgas_Integrand(rr)
      
      implicit none
      
      double precision rr,Mgas_Integrand
      double precision t3
      
	if(rr<rmin) then 
	  t3=Rhogasfunc(rmin)
        Mgas_Integrand=4.0*Pi*rmin*rmin*t3
	else
	  t3=Rhogasfunc(rr)
        Mgas_Integrand=4.0*Pi*rr*rr*t3
	endif
	
      return
      end function Mgas_Integrand

!======================================================================

	function yfunc(rr)
	
	implicit none
	
	double precision rr,yfunc,result
	
	if(rr<rmin) then 
	  call interp1d_even(logYarray,logr,n,phlog10(rmin),result)
      else if(rr>rlimit) then
        yfunc=0.
        return
	else
	  call interp1d_even(logYarray,logr,n,phlog10(rr),result)
	endif
	yfunc=10.0**result
	
	return
	end function yfunc

!======================================================================

function Rhogasfunc(rr)
	
	implicit none
	
	double precision rr,Rhogasfunc,result
	
	if(rr<rmin) then 
	  	call interp1d_even(logRhogas,logr,n,phlog10(rmin),result)
	else if(rr>rlimit) then
        	Rhogasfunc=0.
        	return
	else
	  	call interp1d_even(logRhogas,logr,n,phlog10(rr),result)
	endif
	Rhogasfunc=10.d0**result
	
end function Rhogasfunc

!======================================================================

	function Pgasfunc(rr)
	
	implicit none
	
	double precision rr,Pgasfunc,result
	
	if(rr<rmin) then 
	  	call interp1d_even(logPgas,logr,n,phlog10(rmin),result)
	else if(rr>rlimit) then
        	Pgasfunc=0.
        	return
      	else
	  	call interp1d_even(logPgas,logr,n,phlog10(rr),result)
	endif
	
	Pgasfunc=10.0**result
	
	end function Pgasfunc

!======================================================================

function PgasIntegrand(zz)
	
	implicit none
	
	double precision PgasIntegrand,zz,rr

	rr = sqrt( uu * uu + zz * zz )
	if( rr < rmin ) then 
	  	rr = rmin
	  	PgasIntegrand = Pgasfunc(rr)
	else if( rr > rlimit ) then
       	 	PgasIntegrand = 0.d0
	else
	  	PgasIntegrand = Pgasfunc(rr) 
	endif
	
end function PgasIntegrand

!======================================================================

	double precision function totMass(x,r200,rc,rhocrit,r)
	
	implicit none
	
	double precision x,r200,rc,rhocrit
      double precision r

	r=rx(x,r200,rc)
      totMass=(4d0*pi/3d0)*(r**3)*x*rhocrit
      
	end function totMass

!======================================================================

	double precision function rx(x,r200,rc)
	
	implicit none
	
	double precision x,r200,rc

	rx=sqrt(max(200d0*(r200**2+rc**2)/x-rc**2,0.d0))
      
	end function rx

!======================================================================

	subroutine getGassMass(rc,beta,rho0,r200,Mgas)
	
	implicit none
	
	double precision rc,beta,r200,rho0
      double precision Mgas(5) !Mgas value at r1500, r1000, r150, r178 & r500
      double precision rhi,rlo
	double precision eps
	parameter(eps=1d-4)
	
      gd1=rc
      gd2=beta
      
      !r1500
      rlo=0.d0
      rhi=rx(1500d0,r200,rc)
	call qtrap(betaProfileInt,rlo,rhi,eps,Mgas(1))
      Mgas(1)=Mgas(1)*rho0
      
      !r1000
      rlo=rhi
      rhi=rx(1000d0,r200,rc)
	call qtrap(betaProfileInt,rlo,rhi,eps,Mgas(2))
      Mgas(2)=Mgas(2)*rho0+Mgas(1)
      
      !r500
      rlo=rhi
      rhi=rx(500d0,r200,rc)
	call qtrap(betaProfileInt,rlo,rhi,eps,Mgas(3))
      Mgas(3)=Mgas(3)*rho0+Mgas(2)
      
      !r178
      rlo=rhi
      rhi=rx(178d0,r200,rc)
	call qtrap(betaProfileInt,rlo,rhi,eps,Mgas(4))
      Mgas(4)=Mgas(4)*rho0+Mgas(3)
      
      !r150
      rlo=rhi
      rhi=rx(150d0,r200,rc)
	call qtrap(betaProfileInt,rlo,rhi,eps,Mgas(5))
      Mgas(5)=Mgas(5)*rho0+Mgas(4)
      
	end subroutine getGassMass

!======================================================================

	double precision function betaProfileInt(r)
	
	implicit none
	
	double precision r
	
      if(r<rmin) then
		betaProfileInt=4d0*pi*(rmin**2)/((1d0+((rmin/gd1)**2))**(1.5d0*gd2))
	else
		betaProfileInt=4d0*pi*(r**2)/((1d0+((r/gd1)**2))**(1.5d0*gd2))
      endif
      
	end function betaProfileInt

!======================================================================
!======================================================================
          FUNCTION  CalcBETAsphVolY(r_delta1 ,r_delta2 )

          IMPLICIT NONE

          double precision , PARAMETER   :: eps=1.0d-4
          double precision               :: lim1 , lim2 , S
          double precision               :: r_delta1 ,r_delta2
          double precision               :: CalcBETAsphVolY

          lim1=r_delta1
          lim2=r_delta2

          CALL qtrap(BYsphVolInt,lim1,lim2,eps,S)
          CalcBETAsphVolY=S

          END FUNCTION CalcBETAsphVolY
!===============================================================================================
         FUNCTION BYsphVolInt(r3D)

          IMPLICIT NONE
          double precision               :: r3D , BYsphVolInt 

          BYsphVolInt=(4.0d0 *pi*r3D*r3D)/&
            ((1.d0+((r3D*r3D)/(rc*rc)))** &
                                   (3.d0*beta/2.d0))
				   
          END FUNCTION BYsphVolInt

!================================================================================================
      FUNCTION yintegral_Planck(zz)
	
	IMPLICIT NONE
	
        double precision, dimension(:), intent(in) :: zz
        double precision, dimension(size(zz)) :: yintegral_Planck
        double precision, dimension(size(zz)) :: rr
        integer :: imax, nz, i

        rr=SQRT(uu*uu +zz*zz)

        where (rr<thetamin_Planck) rr=thetamin_Planck
        imax = ifirstloc(rr>thetalimit_Planck)-1
        nz = size(zz)

        yintegral_Planck(1:imax)=ycoeff_Planck*((rr(1:imax)/thetas_Planck)**(-c_GNFW))*&	     
	 ( (1.d0 + ((rr(1:imax)/thetas_Planck)**(a_GNFW)))**((c_GNFW -b_GNFW)/a_GNFW))

        if (imax<nz) yintegral_Planck(imax+1:nz)=0d0
	
     END FUNCTION yintegral_Planck

!========================================================

	FUNCTION yfunc_Planck(rr)
	
	IMPLICIT NONE
	
	double precision     :: rr,yfunc_Planck,result
	
	if(rr<thetamin_Planck) then 
	  call interp1d_even(logyarray_Planck,logtheta_Planck,n,phlog10(thetamin_Planck),result)
      else if(rr>thetalimit_Planck) then
        yfunc_Planck=0.
        return
	else
	  call interp1d_even(logyarray_Planck,logtheta_Planck,n,phlog10(rr),result)
	endif
	yfunc_Planck=10.0**result
	
	return
	END FUNCTION yfunc_Planck

!=========================================================================

          FUNCTION  CalcGNFW_PlancksphVolY(theta1 ,theta2 )

          IMPLICIT NONE

          double precision               :: lim1 , lim2 , S
          double precision               :: theta1 ,theta2
          double precision               :: CalcGNFW_PlancksphVolY

          lim1=theta1
          lim2=theta2

          CalcGNFW_PlancksphVolY=qtrap2(GNFW_PlanckYsphVolInt,lim1,lim2)

          END FUNCTION CalcGNFW_PlancksphVolY
!===============================================================================================

         FUNCTION GNFW_PlanckYsphVolInt(theta3D)
          IMPLICIT NONE
          double precision, dimension(:), intent(in)               :: theta3D
          double precision, dimension(size(theta3D))               :: GNFW_PlanckYsphVolInt

          GNFW_PlanckYsphVolInt=((theta3D/thetas_Planck)**(2.0-c_GNFW))*&
	                         ((1.d0+((theta3D/thetas_Planck)**(a_GNFW)))**((c_GNFW - b_GNFW)/(a_GNFW)))
				   
          END FUNCTION GNFW_PlanckYsphVolInt

!================================================================================================

	FUNCTION DM_GNFWgasVol(radius,rs_DM,rp_GNFW,a_GNFW,b_GNFW,c_GNFW) !adapted to use qtrap2, which takes an array returning function as an argument instead of a function which relies on interpolation. kj 02/02/17
	
        IMPLICIT NONE

        double precision , PARAMETER   :: eps=1.0d-4
        !double precision               :: result
        double precision               :: radius,rs_DM,rp_GNFW,a_GNFW,b_GNFW,c_GNFW	
        double precision               ::DM_GNFWgasVol
        !CALL qtrap(DM_GNFWsphVolInt,1.0d-4,radius,eps,result)
        !DM_GNFWgasVol=result
	DM_GNFWgasVol=qtrap2(DM_GNFWsphVolInt,1.0d-4,radius)	
		
	
       END FUNCTION DM_GNFWgasVol	
!=========================================================================

	FUNCTION DM_GNFWsphVolInt(zz) !adapted for use in qtrap2 function. Now takes array or radii as input argument, and returns array of values. kj 02/02/17.
	
        IMPLICIT NONE	
	!double precision               :: r         
        !double precision               ::DM_GNFWsphVolInt
        double precision, dimension(:), intent(in)             :: zz !need intent(in) array to fit interface of qtrap2
	double precision, dimension(size(zz))       ::rr !array that is actually used, as it may need to be changed in function at limits         
        double precision, dimension(size(zz))       ::DM_GNFWsphVolInt
	integer ::	 imax, nr
		
	       
        !IF(r<rmin)THEN
	!DM_GNFWsphVolInt=((rmin*rmin*rmin)/((DLOG(1.0 + (rmin/rs_DM))) - (1.0/(1.0 +(rs_DM/rmin))))) *&
        !                 ((rmin/rp_GNFW)**(-1.0*c_GNFW))*&
        !                ( (1.0 + (rmin/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW +b_GNFW -c_GNFW)/a_GNFW))*&
        !                    ((b_GNFW*((rmin/rp_GNFW)**(a_GNFW))) + c_GNFW)
	
	!ELSEIF(r>rmax) THEN
	
	! DM_GNFWsphVolInt=0.d0  
	   
        !ELSE	
	
	!DM_GNFWsphVolInt=((r*r*r)/((DLOG(1.0 + (r/rs_DM))) - (1.0/(1.0 +(rs_DM/r))))) *&
        !                 ((r/rp_GNFW)**(-1.0*c_GNFW))*&
        !                ( (1.0 + (r/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW +b_GNFW -c_GNFW)/a_GNFW))*&
        !                    ((b_GNFW*((r/rp_GNFW)**(a_GNFW))) + c_GNFW)
	                  
	
	
	!ENDIF

	rr=zz !use values passed to function
	where(rr < rmin) rr = rmin
	imax = ifirstloc(rr>rmax)-1
        nr = size(rr)
	DM_GNFWsphVolInt(1:imax) = ((rr(1:imax)*rr(1:imax)*rr(1:imax))/((DLOG(1.0 + (rr(1:imax)/rs_DM))) - &
			 (1.0/(1.0 +(rs_DM/rr(1:imax)))))) *&
                         ((rr(1:imax)/rp_GNFW)**(-1.0*c_GNFW))*&
                        ( (1.0 + (rr(1:imax)/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW +b_GNFW -c_GNFW)/a_GNFW))*&
                            ((b_GNFW*((rr(1:imax)/rp_GNFW)**(a_GNFW))) + c_GNFW)
	if (imax<nr) DM_GNFWsphVolInt(imax+1:nr)=0d0
	
	END FUNCTION DM_GNFWsphVolInt	       
!================================================================================================
          FUNCTION GNFWmodel3D(r,rp_GNFW,a_GNFW,b_GNFW,c_GNFW,Pei)

        IMPLICIT NONE
        double precision            :: r,rp_GNFW,a_GNFW,b_GNFW,c_GNFW, Pei
        double precision            :: GNFWmodel3D

        GNFWmodel3D = Pei/(( (r/rp_GNFW)**(c_GNFW))* &
                               ((1.d0+((r/rp_GNFW)**(a_GNFW)))** &
                               ((b_GNFW - c_GNFW)/a_GNFW))) 
      END FUNCTION  GNFWmodel3D 	  
!=======================================================================================       
          FUNCTION CalcGNFWsphVolY(r_delta1 , r_delta2)

          IMPLICIT NONE

          double precision               :: lim1 , lim2 , S
          double precision               :: r_delta1 ,r_delta2
          double precision               :: CalcGNFWsphVolY

          lim1=r_delta1
          lim2=r_delta2

          CalcGNFWsphVolY=qtrap2(GNFWYsphVolInt,lim1,lim2)

          END FUNCTION CalcGNFWsphVolY
!================================================================================================
          FUNCTION GNFWYsphVolInt(r3D)

          IMPLICIT NONE
          double precision, dimension(:), intent(in)              :: r3D
          double precision, dimension(size(r3D))                  :: GNFWYsphVolInt

          GNFWYsphVolInt=(r3D*r3D)/   &
             (((r3D/rp_GNFW)**(c_GNFW))*  &
             ( ( 1.d0 +((r3D/rp_GNFW)**(a_GNFW)) ) &
	     **((b_GNFW -c_GNFW)/a_GNFW) )) 
          
          END FUNCTION GNFWYsphVolInt
!=================================================================================================
          FUNCTION GNFWyintegrand(zz)

          IMPLICIT NONE

          double precision, dimension(:), intent(in) :: zz
          double precision, dimension(size(zz)) :: GNFWyintegrand, GNFWmodel3D
          double precision, dimension(size(zz)) :: rr
          integer :: i, imax, nz

          rr=SQRT(uu*uu +zz*zz)

          where (rr<rmin) rr=rmin
          imax = ifirstloc(rr>rmax)-1
          nz = size(zz)

          GNFWmodel3D(1:imax) = 1d0/(( (rr(1:imax)/rp_GNFW)**(c_GNFW))* &
                               ((1.d0+((rr(1:imax)/rp_GNFW)**(a_GNFW)))** &
                               ((b_GNFW - c_GNFW)/a_GNFW))) 
          GNFWyintegrand(1:imax)=ycoeff_GNFW*GNFWmodel3D(1:imax)

          if (imax<nz) GNFWyintegrand(imax+1:nz)=0d0

         END  FUNCTION GNFWyintegrand

!================================================================================================	  

          FUNCTION yGNFWfunc2(rr)
	
          IMPLICIT NONE
	
          double precision              :: rr
          double precision              :: yGNFWfunc2
          double precision              :: result
	
           IF(rr<rmin) THEN 
              CALL Interp1d_even(logyarray_GNFW,logr,n,phlog10(rmin),result)
           ELSEIF(rr>rmax) then
             yGNFWfunc2=0.
             RETURN
	   ELSE
	      CALL Interp1d_even(logyarray_GNFW,logr,n,phlog10(rr),result)
           ENDIF
	
           yGNFWfunc2=10.0**result
	
	
        END FUNCTION yGNFWfunc2
	  
!adapted for GM=6 kj 01/02/17
!================================================================================================

	function Ein_GNFWgasVol(a_GNFW, b_GNFW, c_GNFW, aEin_DM, radius, rm2_DM, rhom2_DM, rp_GNFW)
	!uses qtrap2
        implicit none

        double precision , parameter   :: eps=1.0d-4
        !double precision               :: result
        double precision               :: radius, rm2_DM,rp_GNFW
	double precision 	       :: a_GNFW,b_GNFW,c_GNFW,aEin_DM, rhom2_DM	
        double precision               ::Ein_GNFWgasVol
        !call qtrap(Ein_GNFWsphVolInt,1.0d-4,radius,eps,result)
        !Ein_GNFWgasVol=result	
	Ein_GNFWgasVol=qtrap2(Ein_GNFWsphVolInt,1.0d-4,radius)		
	
       end function Ein_GNFWgasVol	
			
!=========================================================================

	function Ein_GNFWsphVolInt(zz)
	!commented out parts are from when qtrap was used instead of qtrap2
        implicit none	
        !double precision               :: r         
        !double precision               :: Ein_GNFWsphVolInt
	integer 		       :: flag
	double precision, dimension(:), intent(in)             :: zz !need intent(in) array to fit interface of qtrap2
	double precision, dimension(size(zz))       ::rr !array that is actually used, as it may need to be changed in function at limits
	double precision, dimension(size(zz)) :: Ein_GNFWsphVolInt
	double precision, dimension(size(zz)) :: gamma_arr !used to store gammds values as it doesn't take arrays as arguments
	integer			       :: imax, nr, i=1
		
	       
        !if(r<rmin)then
	!Ein_GNFWsphVolInt=rmin**3 * (b_GNFW * (rmin / rp_GNFW) ** a_GNFW + c_GNFW) / &
	!		 (gammaFun(3.d0 / aEin_DM) * &
	!		 gammds(2.d0 / aEin_DM * (rmin / rm2_DM) ** aEin_DM, 3.d0 / aEin_DM, flag) * &
	!		  (rmin / rp_GNFW) ** (c_GNFW) * (1.d0 + (rmin / rp_GNFW)**a_GNFW)** &
	!		  ((a_GNFW + b_GNFW - c_GNFW) / a_GNFW))
	
	!elseif(r>rmax) then
	
	! Ein_GNFWsphVolInt=0.d0  
	   
        !else	
	
	!Ein_GNFWsphVolInt=r**3 * (b_GNFW * (r / rp_GNFW) ** a_GNFW + c_GNFW) / &
	!		 (gammaFun(3.d0 / aEin_DM) * &
	!		 gammds(2.d0 / aEin_DM * (r / rm2_DM) ** aEin_DM, 3.d0 / aEin_DM, flag) * &
	!		  (r / rp_GNFW) ** (c_GNFW) * (1.d0 + (r / rp_GNFW)**a_GNFW)** & 
	!		  ((a_GNFW + b_GNFW - c_GNFW) / a_GNFW))
	                  
	
	
	!endif
	
	!if(flag == 1 .or. flag == 2) then
	!	write(*,*) "Gamma calculation for Ein_GNFWsphVolInt has screwed up"
	!	flag = 1
	!	return
	!endif
	rr = zz !use values passed into function
	where(rr < rmin) rr = rmin
	imax = ifirstloc(rr>rmax)-1
        nr = size(rr)
	do i=1,imax
		gamma_arr(i) = gammds(2.d0 / aEin_DM * (rr(i) / rm2_DM) ** aEin_DM, 3.d0 / aEin_DM, flag)
	enddo 

	Ein_GNFWsphVolInt(1:imax) = rr(1:imax)**3 * (b_GNFW * (rr(1:imax) / rp_GNFW) ** a_GNFW + c_GNFW) / &
			 (gammaFun(3.d0 / aEin_DM) * gamma_arr(1:imax) *&
			  (rr(1:imax) / rp_GNFW) ** (c_GNFW) * (1.d0 + (rr(1:imax) / rp_GNFW)**a_GNFW)** & 
			  ((a_GNFW + b_GNFW - c_GNFW) / a_GNFW))
	if(flag == 1 .or. flag == 2) then
		write(*,*) "Gamma calculation for Ein_GNFWsphVolInt has screwed up"
		flag = 1
		return
	endif
	if (imax<nr) Ein_GNFWsphVolInt(imax+1:nr)=0d0

	end function Ein_GNFWsphVolInt	          
!================================================================================================
!=======================================================================
!Function used to determine r500 using Newton Raphson to invert function. X = r500 / rm2

	function X_func(X)
	double precision :: X_func
	double precision :: X
	integer		 :: flag

	X_func = X**3.d0 / (gammaFun(3.d0 / aEin_DM) * &
		 gammds(2.d0 / aEin_DM * (X) ** (aEin_DM), 3.d0 / aEin_DM, flag)) &
		 - (3d0*rhom2_DM / 500d0 * (aEin_DM/2d0)**(3d0/aEin_DM) &
		 * exp(2d0/aEin_DM) / (rhocritz*aEin_DM))
	if(flag == 1 .or. flag == 2) then
		flag = 1
		write(*,*) "Gamma calculation for r500 Newton Raphson has screwed up"
		return
	endif

	end function X_func
	
!========================================================================


end module GasModels
