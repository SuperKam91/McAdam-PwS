module CheckPars1
	use params
	use constants
	use utilities
	use matrix_utils
	use MassModels
	use priors1
	use p_alpha
	use massfunction
        use ReadWrite
      
	implicit none
	
	logical crapAlpha
	logical :: noCluster=.false.

contains
!=======================================================================

subroutine CheckPars(Cube,flag)

	implicit none
	
	double precision Cube(Ndim)
	integer flag,j,k

! Rescale parameters:

	do j = 1, Ndim
		crapAlpha = .false.
	  	
		call Rescale(Cube(j),j)
		
		if( crapAlpha ) then
			flag = 1
            		return
	  	endif
	enddo
      
!   If we are fitting mass distribution, check various mass parameter:

	if( Mass == 1 ) call CheckMassPars(flag)
      
!   If we are fitting gas distribution, check various gas parameters:

	if( Gas == 1 ) then
          do k = 1, NAtoms
            call CheckGasPars(k,flag)
          enddo
        endif

        call CheckzPars(flag)

	if (hyperparameters == 1) call Checkhypers(flag)
      
end subroutine CheckPars
	
!=======================================================================
! Check redshift

subroutine CheckzPars(flag)

	implicit none
	
	integer flag,k 

        do k = 1, NAtoms
          if (z(k).le.0) then
            flag = 1
          endif
        enddo

end subroutine CheckzPars
      
!=======================================================================
!! Check various mass parameters:

subroutine CheckMassPars(flag)

	implicit none
	
	integer flag,k 

!-----------------------------------------------------------------------

! Loop over all atoms:

	do k=1,Natoms
      
      		!incompatible priors for M & Z?
		if( ( MassModel /= 1 .and. ( z_PriorType(k) == 8 .or. Mass_PriorType(k,2) == 8 ) ) .or. &
      		( MassModel == 1 .and. NFWstyle(k) /= 2 .and. ( z_PriorType(k) == 8 .or. Mass_PriorType(k,2) == 8 ) ) &
            	.or. ( MassModel == 1 .and. NFWstyle(k) == 2 .and. ( ( z_PriorType(k) == 8 .and. &
            	Mass_PriorType(k,2) /= 8 ) .or. ( z_PriorType(k) /= 8 .and. Mass_PriorType(k,2) == 8 ) ) ) ) then
			call halt_program("incompatible priors used for M200 & z")
		endif
     	enddo 
      		
end subroutine CheckMassPars
	
!=======================================================================
! Check various gas parameters:

subroutine CheckGasPars(k,flag)

	implicit none
	
	integer flag, k  
	double precision index_g, rr
        double precision cc, thetaE

	if( znow ) then
      		rhocritz=rhocrit
	else
      		rhocritz=rhocritofz(z(k))
	endif
	
      call lookUp1D(Dn,z(k),lookD(:,1),lookD(:,2),D)

        if(GasModel==1) then

 	  	if(Mass==1.and.MassModel==3.and.CPLstyle(k)==-1) then 
	    		rc=MassPars(k,1)/1000.d0
        	else
          		rc=GasPars(k,1)/1000.d0
	  	endif
        
        	beta=GasPars(k,2) 
		  	    
        	if(Mass==1) then
	   		if(MassModel==1) then
				!NFW profile for potential:
          			if(NFWstyle(k)==1) then 
	      				call halt_program('Cannot calculate r200 in this NFW style.')
	    			elseif(NFWstyle(k)==2) then
	      				M200=MassPars(k,2)
	      				r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
	    			elseif(NFWstyle(k)==3) then
	      				M200=MassPars(k,1)
	      				r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
	    			endif
         		elseif(MassModel==2) then
				!SIS profile for potential:
          			if(SISstyle(k)==1) then 
	      				M200=MassPars(k,1)
	      				r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
          			elseif(SISstyle(k)==2) then 
	      				r200=sqrt(MassPars(k,1)*sec2rad*D*SigcritE/2.0)
	      				r200=r200*sqrt(3.0/(200.0*rhocritz*pi))
	    			endif
         		elseif(MassModel==3.and.CPLstyle(k)==-1) then
 	     			M200=MassPars(k,2)
           			r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
         		else  
	    			call halt_program('Can only do joint fit with NFW, SIS and CIS potentials!')
	   		endif

	  	else
          		if(Betastyle(k)==1) then 
				!Unknown potential, normalise to fixed radius:
	      			r200=rmass
	    		elseif(Betastyle(k)==2.or.Betastyle(k)==3) then      
				!HSE potential, compute r200:
		 		r200=max((9.d0/(4.d0*pi*Gm*200.d0*rhocritz))*beta*TPars(k,1)-rc*rc,0.d0)
             			r200=sqrt(r200)
            			if(Tmodel==2) then 
              				index_g=TPars(k,2)
              				rr=2.0/(3.0*beta*(index_g-1.0)+2.0)
              				r200=(r200*TPars(k,2))**rr
				endif
          			M200=4.0*pi*200.0*rhocritz*r200*r200*r200/3.0
			elseif(Betastyle(k)==4) then
				M200=GasPars(k,3)	
				r200=((3.d0*M200)/(4.d0*pi*200.d0*rhocritz))**(1.d0/3.d0)				
          		else 
	      			call halt_program('Unrecognised Beta style.')
	    		endif
        	endif
		
		if( M200 <= 0d0 .or. r200 <= 0d0 .or. beta <= 0d0 .or. rc <= 0d0 ) then
            		flag = 1
                        if (verbose) write(*,*) 'M200, r200, beta, rc = ', M200, r200, beta, rc
                  	return
		endif	

		!Gas density profile normalisation: by mass within fixed radius rmass 
		!(Betastyle=1), by mass within r200 from HSE (Betastyle=2), or by gas 
		!fraction within within r200 from HSE:    
        	if(Betastyle(k)==1) then
          		Mgas200=GasPars(k,3)
        	elseif(Betastyle(k)==2) then
          		Mgas200=GasPars(k,3)
        	elseif(Betastyle(k)==3) then
          		Mgas200=GasPars(k,3)*M200
         	elseif(Betastyle(k)==4) then
			fgas200=GasPars(k,4)
			IF(fgas200 .LT. 0.0)THEN
			 	flag=1
                                if (verbose) write(*,*) 'fgas200=', fgas200
			 	RETURN
			ENDIF	
			Mgas200=M200*fgas200
        	endif
				
		if( Mgas200 < 0d0 ) then
            		flag = 1
                        if (verbose) write(*,*) Mgas200
                  	return
		endif	

! 	Isothermal temperature profile:
	elseif(GasModel==2.and.TModel==1) then
	    
        	if(MassModel==1) then

!       		NFW profile for potential:
          		if(NFWstyle(k)==1) then 
	      			rs=MassPars(k,1) 
	      			ps=MassPars(k,2) 
	    		elseif(NFWstyle(k)==2) then
	      			cc=MassPars(k,1)
	      			M200=MassPars(k,2)
	      			r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
	      			rs=r200/cc
	      			ps=M200/(4.0*pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
	    		elseif(NFWstyle(k)==3) then
				!M200 now contains M100, so that Komatsu and Seljak's relation can be used!!
            			M200=MassPars(k,2)
            			cc=6.0*(M200/1d14)**(-0.2)
            			r200=(M200/(100.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
            			rs=r200/cc
            			ps=M200/(4.0*Pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
          		elseif(NFWstyle(k)==4) then
            			thetaE=MassPars(k,1)
            			M200=MassPars(k,2)
            			rE=ThetaE*sec2rad*D
            			scE=sigcritE
            			r200=(M200/(200.0*(4.0*Pi/3.0)*rhocritz))**(1.0/3.0)
            			rs=FINDR(rsfunc,1.0d-4,1.0d4,1.0d-4)
            			cc=r200/rs
            			ps=M200/(4.0*Pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
          		else
            			call halt_program('Unrecognised NFWstyle.')
          		endif
                  
                  	!sanity check
            		if( rs <= 0d0 .or. ps <= 0d0 ) then
                                if (verbose) write(*,*) 'rs, ps=', rs, ps
            			flag = 1
                  		return
			endif
		else
	    		call halt_program('Cannot do non-NFW potentials!')
                endif

! 	Polytropic temperature profile:
	elseif(GasModel==2.and.TModel==2) then
		  
!     		First calculate potential, then the temperature, then the density:	
	
        		if(MassModel==1) then

!       			NFW profile for potential:
          			if(NFWstyle(k)==1) then 
	      				rs=MassPars(k,1) 
	      				ps=MassPars(k,2) 
	    			elseif(NFWstyle(k)==2) then
	      				cc=MassPars(k,1)
	      				M200=MassPars(k,2)
	      				r200=(M200/(200.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
	      				rs=r200/cc
	      				ps=M200/(4.0*pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
	    			elseif(NFWstyle(k)==3) then
! 					M200 now contains M100, so that Komatsu and Seljak's relation can be used!!
            				M200=MassPars(k,2)
            				cc=6.0*(M200/1d14)**(-0.2)
            				r200=(M200/(100.0*(4.0*pi/3.0)*rhocritz))**(1.0/3.0)
            				rs=r200/cc
            				ps=M200/(4.0*Pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
          			elseif(NFWstyle(k)==4) then
            				thetaE=MassPars(k,1)
            				M200=MassPars(k,2)
            				rE=ThetaE*sec2rad*D
            				scE=sigcritE
            				r200=(M200/(200.0*(4.0*Pi/3.0)*rhocritz))**(1.0/3.0)
            				rs=FINDR(rsfunc,1.0d-4,1.0d4,1.0d-4)
            				cc=r200/rs
            				ps=M200/(4.0*Pi*rs*rs*rs*(log(1.0+cc)-cc/(cc+1.0)))
          			else
            				call halt_program('Unrecognised NFWstyle.')
          			endif
         
                  		!sanity check
            			if( rs <= 0d0 .or. ps <= 0d0 ) then
                                        if (verbose) write(*,*) 'rs, ps = ', rs, ps
            				flag = 1
                  			return
				endif
	  		else
	    			call halt_program('Cannot do non-NFW potentials!')
                        endif
! GNFW-Planck model:

        ELSEIF(GasModel==3 )THEN
            thetas_Planck = GasPars(k,1)          !arcmin  
            Ytot_Planck= GasPars(k,2)             !arcmin2
	    c_GNFW=GasPars(k,3)
	    a_GNFW=GasPars(k,4)
	    b_GNFW=GasPars(k,5)
	    c500_GNFW=GasPars(k,6)
	    
! Some checks to avoid the profile blowing up
            if ( a_GNFW < 0.1d0 .or. b_GNFW < 0. .or. c_GNFW < 0. ) then
                if (verbose) write(*,*) 'a_GNFW, b_GNFW, c_GNFW = ', a_GNFW, b_GNFW, c_GNFW
                flag = 1
                return
            endif

            if ( abs(b_GNFW-3.0d0)/a_GNFW<0.01 .or. abs(b_GNFW-1.0d0)/a_GNFW<0.01 ) then
                if (verbose) write(*,*) 'abs(b_GNFW-3.0d0)/a_GNFW, abs(b_GNFW-1.0d0)/a_GNFW = ', abs(b_GNFW-3.0d0)/a_GNFW, abs(b_GNFW-1.0d0)/a_GNFW
                flag = 1
                return
            endif

            if ( abs(c_GNFW-1.0d0)/a_GNFW<0.01 .or. abs(b_GNFW-c_GNFW)/a_GNFW<0.01 ) then
                if (verbose) write(*,*) 'abs(c_GNFW-1.0d0)/a_GNFW, abs(b_GNFW-c_GNFW)/a_GNFW = ', abs(c_GNFW-1.0d0)/a_GNFW, abs(b_GNFW-c_GNFW)/a_GNFW
                flag = 1
                return
            endif

! Hard-coded limit on thetas, still necessary?
	    if( thetas_Planck <= 1.30d0 ) then
                if (verbose) write(*,*) 'thetas_Planck', thetas_Planck
	    	flag = 1
		return
	    endif

	    if( Ytot_Planck < 0.d0 ) then
                if(verbose) write(*,*) 'Ytot_Planck', Ytot_Planck
		flag = 1
		return
	    endif

            theta500_Planck=thetas_Planck *c500_GNFW
	    y0_int_Planck = (1.d0/a_GNFW)*Gammafun((-c_GNFW	+1.d0)/a_GNFW) *&
	Gammafun((b_GNFW-1.d0)/a_GNFW)/Gammafun((-c_GNFW + b_GNFW)/a_GNFW)
	
            Ytot_int_Planck=(1.d0/a_GNFW)*Gammafun((-c_GNFW	+3.d0)/a_GNFW) *&
	Gammafun((b_GNFW-3.d0)/a_GNFW)/Gammafun((-c_GNFW + b_GNFW)/a_GNFW)

	    ycoeff_Planck=(1.d0/(4.d0*pi)) *(Ytot_Planck/(thetas_Planck*thetas_Planck*thetas_Planck))*(1.d0/Ytot_int_Planck)
            if ( ycoeff_Planck < 0.d0 ) then
               if (verbose) write(*,*) 'ycoeff_Planck = ', ycoeff_Planck
               flag = 1
               return
            endif
	    	    
!DM-GNFW Model:
!adapted for GM=6 kj 01/02/17
        ELSEIF(GasModel==5 .or. GasModel==6 )THEN
	      MT200_DM = GasPars(k,1)   !M_sun
	      fg200_DM = GasPars(k,2)
  	      c_GNFW=GasPars(k,3)
	      a_GNFW=GasPars(k,4)
	      b_GNFW=GasPars(k,5)
	      c500_GNFW=GasPars(k,6)

! Some checks to avoid the profile blowing up
              if ( a_GNFW < 0.1d0 .or. b_GNFW < 0. .or. c_GNFW < 0. .or. c500_GNFW < 0. ) then
                if (verbose) write(*,*) 'a_GNFW, b_GNFW, c_GNFW, c500_GNFW = ', a_GNFW, b_GNFW, c_GNFW, c500_GNFW
                flag = 1
                return
              endif

              if ( abs(b_GNFW-3.0d0)/a_GNFW<0.01 .or. abs(b_GNFW-1.0d0)/a_GNFW<0.01 ) then
                if (verbose) write(*,*) 'abs(b_GNFW-3.0d0)/a_GNFW, abs(b_GNFW-1.0d0)/a_GNFW = ', abs(b_GNFW-3.0d0)/a_GNFW, abs(b_GNFW-1.0d0)/a_GNFW
                flag = 1
                return
              endif

              if ( abs(c_GNFW-1.0d0)/a_GNFW<0.01 .or. abs(b_GNFW-c_GNFW)/a_GNFW<0.01 ) then
                if (verbose) write(*,*) 'abs(c_GNFW-1.0d0)/a_GNFW, abs(b_GNFW-c_GNFW)/a_GNFW = ', abs(c_GNFW-1.0d0)/a_GNFW, abs(b_GNFW-c_GNFW)/a_GNFW 
                flag = 1
                return
              endif
	      
              IF(fg200_DM .LT. 0.0 .OR. MT200_DM .LT. 0.0)THEN
                 if (verbose) write(*,*) 'fg200_DM, MT200_DM = ', fg200_DM, MT200_DM
		 flag=1
		 RETURN
	      ENDIF	
	
	      !Ein_DM-GNFW Model:
	      IF(GasModel==6 )THEN
		 aEin_DM = GasPars(k,7)
		 if (aEin_DM < 0.05d0 .or. aEin_DM > 10.d0) then
	            if (verbose) write(*,*) 'aEin_DM = ', aEin_DM
		    flag = 1
		    return
	         endif
	      ENDIF
	ENDIF
	
end subroutine CheckGasPars

!=======================================================================	

!=======================================================================
! Check hyperparameters

subroutine Checkhypers(flag) !should probably pass dummy arguments for alpha_L1 and alpha_L2, but will keep consistent with rest

	implicit none
	integer flag 
	alpha_L1 = hyperPars(1,1)
	if (alpha_L1 < -1.d3 .or. alpha_L1 > 1.d3) flag = 1
	if (Nhyper == 2) then 
		alpha_L2 = hyperPars(1,2)
		if (alpha_L2 < -1.d3 .or. alpha_L2 > 1.d3) flag = 1
	endif
	return

end subroutine Checkhypers
      
!=======================================================================


subroutine Rescale(value,ival)

    implicit none

    integer ival,i,ii,j,k,tot_pars,tmp_prior
    logical nuis
    double precision value,dummy
    double precision a(2), b(2), c(2), p(2)
    double precision cos2, sin2, sig_x2, sig_y2, B1, C1

    tot_pars=NPars

    if( ival > NAtoms * tot_pars + Nhyper ) then
        nuis=.true.
    else
        nuis=.false.
    endif

    if( .not. nuis ) then

        !     	Rescale the parameters of the atom, labelled j=1..NAtoms; i runs from 1 to NPars:

        if (ival<= NAtoms*NPars) then
            	j=(ival-1)/NPars+1
            	i=ival-(j-1)*NPars
	else !hyperparameter
		i=ival
	endif
        if(i <= NGeoPars) then	
            ii = i
            tmp_prior=abs(Geo_PriorType(j,ii))
            if( tmp_prior == 9 .and. ii <= 2 ) then
                GeoPars(j,ii)=Prior(1,value,Geo_Tri_Prior(j,ii,1),Geo_Tri_Prior(j,ii,2),Geo_Tri_Prior(j,ii,3))

                !cluster position in triangle constraint
                if( ii == 2 ) then
                    a(1:2) = Geo_Prior(j,1:2,1)
                    b(1:2) = Geo_Prior(j,1:2,2)
                    c(1:2) = Geo_Prior(j,1:2,3)
                    p(1:2) = GeoPars(j,1:2)
                    if( .not.inTriangle(a, b, c, p) ) crapAlpha = .true.
                endif
            else
                GeoPars(j,ii)=Prior(tmp_prior,value,Geo_Prior(j,ii,1),Geo_Prior(j,ii,2),Geo_Prior(j,ii,3))
            endif

        else if( i <= NGeoPars + 1 ) then
            tmp_prior=abs(z_PriorType(j))
            if(z_PriorType(j)==8) then
                !mass function prior then call the lookup table
                z(j) = value
                !call lookUpZ( value, MassPars(j,2), z(j) )
            else
                z(j) = Prior(tmp_prior, value, z_Prior(j,1), z_Prior(j,2),z_Prior(j,3))
                if( z(j) < 0.01d0 .or. z(j) > zdmax ) crapAlpha = .true.
            endif


        else if( Mass == 1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass ) ) then

            ii = i - NGeoPars - 1
            tmp_prior=abs(Mass_PriorType(j,ii))

            if(ii==2 .and. tmp_prior==8) then
                !mass function prior for M200
                !call the lookup table
                call lookUpM(value, MassPars(j,2))
                if(abs(z_PriorType(j))==8) then
                    call lookUpZ(z(j), MassPars(j,2), z(j))
                    if( z(j) < 0.01d0 .or. z(j) > zdmax ) crapAlpha = .true.
                endif
            else
                MassPars(j,ii)=Prior(tmp_prior,value,Mass_Prior(j,ii,1),Mass_Prior(j,ii,2),Mass_Prior(j,ii,3))
            endif

        else if( ((SZ == 1) .or. (PL == 1)) .and. Gas == 1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass + NGasPars * Gas ) ) then

            ii = i - NGeoPars - 1 - NMassPars * Mass
            tmp_prior=abs(Gas_PriorType(j,ii))

            if((BetaStyle(j)==4 .and. ii==3 .and. tmp_prior==8).or. (GasModel==5 .and. ii==1.and. tmp_prior==8) .or. (GasModel==6 .and. ii==1.and. tmp_prior==8)) then !adapted for GM=6 kj 01/02/17
                !mass function prior for M200
                !call the lookup table
                call lookUpM(value, GasPars(j,ii))
                if(abs(z_PriorType(j))==8) then
                    call lookUpZ(z(j), GasPars(j,ii), z(j))
                    if( z(j) < 0.01d0 .or. z(j) > zdmax ) crapAlpha = .true.
                endif
            else if (abs(Gas_PriorType(j,1))==14 .and. (ii==1 .or. ii==2)) then
                ! joint 2D prior on Ytot and theta_s
                if (ii==1) then
                    ! store the unit hypercube value temporarily
                    GasPars(j,1) = value
                else if (ii==2) then
                    sin2=sin(Gas_Prior(j,1,3))*sin(Gas_Prior(j,1,3))
                    cos2=cos(Gas_Prior(j,1,3))*cos(Gas_Prior(j,1,3))
                    sig_y2=Gas_Prior(j,2,2)*Gas_Prior(j,2,2)
                    sig_x2=Gas_Prior(j,1,2)*Gas_Prior(j,1,2)
                    ! marginalised distribution of log10(theta)
                    GasPars(j,1) = GaussianPrior(GasPars(j,1), Gas_Prior(j,1,1), sqrt(2*(sin2*sig_y2+cos2*sig_x2)))
                    B1 = sin2/2./sig_x2 + cos2/2./sig_y2
                    C1 = -sin(2*Gas_Prior(j,1,3))/4./sig_x2+sin(2*Gas_Prior(j,1,3))/4./sig_y2
                    ! conditional distribution of log10(Ytot)
                    GasPars(j,2) = GaussianPrior(value, Gas_Prior(j,2,1)-C1*(GasPars(j,1)-Gas_Prior(j,1,1))/B1, 1./sqrt(B1))
                    GasPars(j,1) = 10.d0**GasPars(j,1)
                    GasPars(j,2) = 10.d0**GasPars(j,2)
                endif
            else	 
                GasPars(j,ii)=Prior(tmp_prior,value,Gas_Prior(j,ii,1),Gas_Prior(j,ii,2),Gas_Prior(j,ii,3))
            endif

        else if( ((SZ == 1) .or. (PL == 1))  .and. Gas == 1 .and. Temperature==1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass + ( NGasPars + NTPars ) * Gas ) ) then

            ii = i - NGeoPars - 1 - NMassPars * Mass - NGasPars * Gas	 
            tmp_prior = abs(T_PriorType(j,ii))
            TPars(j,ii) = Prior(tmp_prior,value,T_Prior(j,ii,1),T_Prior(j,ii,2),T_Prior(j,ii,3))
            !ensure T is always positive
            if( ii == 1 .and. TPars(j,ii) <= 0d0 ) crapAlpha = .true.

	!hyperparameters
        else if ((i > NPars * NAtoms) .and. hyperparameters == 1) then	
		j = 1
                ii = i - NPars * NAtoms
		tmp_prior=abs(hyper_PriorType(j,ii)) !don't see why this is necessary but will use
		hyperPars(j,ii)=Prior(tmp_prior,value,hyper_Prior(j,ii,1),hyper_Prior(j,ii,2),hyper_Prior(j,ii,3))
        endif

    else if(nuis) then

        !       	Rescale the nuisance parameters:

        i=ival-(NAtoms*NPars + Nhyper)
        j=(i-1)/4+1

        if(SZ==1 .and. SourceSubtract==1 .and. i.le.4*NSrc) then	

            ii=i-(j-1)*4
            tmp_prior=abs(Src_PriorType(j,ii))

            !sanity checks for ratio prior
            if(tmp_prior == 13) then
                if(ii /= 3) then
                    call halt_program("ratio prior can only be used for source fluxes.")
                else
                    if(nint(Src_Prior(j,ii,1)) > NSrc .or. nint(Src_Prior(j,ii,1)) < 1 .or. nint(Src_Prior(j,ii,1)) == j) then
                        call halt_program("Incorrect source index in the first argument of ratio prior.")
                    else if(Src_Prior(j,ii,2) <= 0d0) then
                        call halt_program("ratio in ratio prior must be positive.")
                    endif
                endif
            endif

            !source spectral index
            if(ii==4) then
                if(tmp_prior == 0) then
                    SrcPars(j,ii)=Src_Prior(j,ii,1)
                else if(tmp_prior==7 .or. tmp_prior==12) then
                    call query_p_alpha(value,nkernel(j),kernel(j,:,:),SrcPars(j,ii),dummy)
                    crapAlpha=.false.
                    if(SrcPars(j,ii)==-99.d0) then
                        crapAlpha=.true.
                        if (verbose) write(*,*) 'hypercube, j, ii, SrcPars(j,ii) = ', value, j, ii, SrcPars(j,ii)
                        return
                    endif
                else
                    SrcPars(j,ii)=Prior(tmp_prior,value,Src_Prior(j,ii,1),Src_Prior(j,ii,2),Src_Prior(j,ii,3))
                endif
            else
                if(tmp_prior/=13) SrcPars(j,ii)=Prior(tmp_prior,value,Src_Prior(j,ii,1),Src_Prior(j,ii,2),Src_Prior(j,ii,3))
                if(ii==3 .and. SrcPars(j,ii)<0d0) then
                    if (verbose) write(*,*) 'hypercube, j, ii, SrcPars(j,ii) = ', value, j, ii, SrcPars(j,ii)
                    crapAlpha=.true.
                    return
                endif
                if(ii==3 .and. tmp_prior/=13) flux0(j)=SrcPars(j,3)

                ! check if the ratio prior has been imposed on any source fluxes linked with the flux of this source
                if(ii==3) then
                    do k=1,NSrc
                        if(abs(Src_PriorType(k,3))==13 .and. nint(Src_Prior(k,3,1))==j) then
                            SrcPars(k,3) = Src_Prior(k,3,2)*SrcPars(j,3)
                            if(SrcPars(k,3)<0d0) then
                                if (verbose) write(*,*) 'k, SrcPars(k,3) = ', k, SrcPars(k,3)
                                crapAlpha=.true.
                                return
                            endif
                            flux0(k)=SrcPars(k,3)
                        endif
                    enddo
                endif
            endif
        endif

    endif

end subroutine Rescale
	
!=======================================================================

subroutine Rescale_nest(value,ival1,cont)

    implicit none

    integer ival1,ival,i,ii,j,tot_pars,cont
    logical nuis
    double precision value

    ival=ival1
    tot_pars=NPars

    !----------------------------------------------------------------------

    ! 	First determine whether 1D coordinate ival refers to an atomic (or hyperparameter) or 
    ! 	a nuisance parameter:	

    if(ival>NAtoms*tot_pars + Nhyper) then
        nuis=.true.
    else
        nuis=.false.
    endif

    !-----------------------------------------------------------------------


    if( .not. nuis ) then

        !     	Rescale the parameters of the atom, labelled j=1..NAtoms; i runs from 1 to NPars:

        if (ival<= NAtoms*NPars) then
            	j=(ival-1)/NPars+1
            	i=ival-(j-1)*NPars
	else !hyperparameter
		i=ival
	endif

        if(i <= NGeoPars) then

            ii = i
            value = GeoPars(j,ii)

        else if( i <= NGeoPars + 1 ) then

            value = z(j)

        else if( Mass == 1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass ) ) then

            ii = i - NGeoPars - 1
            value = MassPars(j,ii)

        else if( ((SZ == 1) .or. (PL ==1)) .and. Gas == 1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass + NGasPars * Gas ) ) then

            ii = i - NGeoPars - 1 - NMassPars * Mass
            value = GasPars(j,ii)

        else if( ((SZ == 1) .or. (PL == 1)) .and. Gas == 1 .and. Temperature==1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass + ( NGasPars + NTPars ) * Gas ) ) then

            ii = i - NGeoPars - 1 - NMassPars * Mass - NGasPars * Gas
            value = TPars(j,ii)

	!hyperparameters
        else if ((i > NPars * NAtoms) .and. hyperparameters == 1) then	
		j = 1
                ii = i - NPars * NAtoms
		value = hyperPars(j,ii)
        endif

    else if(nuis) then

        !       	Rescale the nuisance parameters:
        i=ival-tot_atoms*NPars + Nhyper
        j=(i-1)/4+1

        if(SZ==1 .and. SourceSubtract==1 .and. i.le.4*NSrc) then	

            ii=i-(j-1)*4
            value = SrcPars(j,ii)

        endif
    endif

end subroutine Rescale_nest

!=======================================================================

subroutine PClass(ival,Ptype)
	
	implicit none
	
	integer ival,i,ii,j,tot_pars
	integer Ptype
	logical nuis

	tot_pars=NPars
	
      	if(ival>tot_atoms*tot_pars + Nhyper) then
        	nuis=.true.
      	else
        	nuis=.false.
      	endif

	if( .not. nuis ) then
	   	if (ival<= NAtoms*NPars) then
            		j=(ival-1)/NPars+1
            		i=ival-(j-1)*NPars
		else !hyperparameter
			i=ival
		endif
         	
            	if( i <= NGeoPars ) then	
         		ii = i	 
			Ptype = Geo_PriorType(j,ii)
            	elseif( i <= NGeoPars + 1 ) then
            		Ptype = z_PriorType(j)
         	elseif( Mass == 1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass ) ) then
            		ii = i - NGeoPars - 1
            		Ptype = Mass_PriorType(j,ii)
         	elseif( ((SZ == 1).or.(PL ==1)) .and. Gas == 1 .and. i <= ( NGeoPars + 1 + NMassPars * Mass + NGasPars * Gas ) ) then
			ii = i - NGeoPars - 1 - NMassPars * Mass	 
			Ptype = Gas_PriorType(j,ii)
      		elseif( ((SZ == 1).or.(PL ==1)) .and. Gas == 1 .and. Temperature ==1 .and.i <= ( NGeoPars + 1 + NMassPars * Mass + ( NGasPars + NTPars ) * Gas ) ) then
			ii = i - NGeoPars - 1 - NMassPars * Mass - NGasPars * Gas	 
			Ptype = T_PriorType(j,ii)
	   	!hyperparameters
       		else if ((i > NPars * NAtoms) .and. hyperparameters == 1) then	
			j = 1
        	        ii = i - NPars * NAtoms
			Ptype = hyper_PriorType(j,ii)
       		 endif

      	elseif(nuis) then
         	i=ival-(NAtoms*NPars + Nhyper)
	   	j=(i-1)/4+1
        
         	if(SZ==1 .and. SourceSubtract==1 .and. i.le.4*NSrc) then	
	    		ii=i-(j-1)*4
            		Ptype=Src_PriorType(j,ii)
			if(Ptype==13) Ptype=0
		else
			Ptype = 0
		endif
	endif

end subroutine PClass
	
!=======================================================================
! Check redshift priors

subroutine CheckzPriors

        implicit none
	integer i
      
	zdmin=1.d10
	zdmax=0.d0
	do i = 1, NAtoms
	    if(z_priorType(i)==0) then
		if( z_prior(i,1) < zdmin ) zdmin = z_prior(i,1)
		if( z_prior(i,1) > zdmax ) zdmax = z_prior(i,1)
	    elseif(z_priorType(i) == 1 .or. z_priorType(i) == 2 .or. z_priorType(i) == 8 ) then
	        if( z_prior(i,1) < zdmin ) zdmin = z_prior(i,1)
	        if( z_prior(i,2) > zdmax ) zdmax = z_prior(i,2)
	    elseif( z_priorType(i) == 3 .or. z_priorType(i) == 4 ) then
		if( max( 0.01d0, z_prior(i,1) - 5d0 * z_prior(i,2) ) < zdmin ) zdmin = max( 0.01d0, z_prior(i,1) - 5d0 * z_prior(i,2) )
		if( z_prior(i,1) + 5d0 * z_prior(i,2) > zdmax ) zdmax = z_prior(i,1) + 5d0 * z_prior(i,2)
	    else
                call halt_program("Can not use a prior other than delta, uniform, log uniform, Gaussian, lognormal or mass function for the cluster redshift. Aborting.")
	    endif
	enddo
      
end subroutine CheckzPriors
	
!=======================================================================
! Check redshift priors for GL when varyzs==1

subroutine CheckzsPriors
        implicit none
        integer j
	
        zsmin=1.d10
	zsmax=0.d0
	do j=1,nzsplanes
	    if(zs_priorType(j)==0) then
		if( zs_prior(j,1) < zsmin ) zsmin = zs_prior(j,1)
		if( zs_prior(j,1) > zsmax ) zsmax = zs_prior(j,1)
	    elseif(zs_priorType(j) == 1 .or. zs_priorType(j) == 2 ) then
		if( zs_prior(j,1) < zsmin ) zsmin = zs_prior(j,1)
		if( zs_prior(j,2) > zsmax ) zsmax = zs_prior(j,2)
	    elseif( zs_priorType(j) == 3 .or. zs_priorType(j) == 4 ) then
		if( max( 0.01d0, zs_prior(j,1) - 5d0 * zs_prior(j,2) ) < zsmin ) zsmin = max( 0.01d0, zs_prior(j,1) - 5d0 * zs_prior(j,2) )
		if( zs_prior(j,1) + 5d0 * zs_prior(j,2) > zsmax ) zsmax = zs_prior(j,1) + 5d0 * zs_prior(j,2)
	    else
	        call halt_program("Can not use a prior other than delta, uniform, log uniform, Gaussian, &
					or lognormal for the source redshift. Aborting.")
	    endif
	enddo

end subroutine CheckzsPriors

!=======================================================================
!! Check mass priors

subroutine CheckMassPriors

	implicit none
	
	integer k, M200_par
        double precision Mmin, Mmax

!-----------------------------------------------------------------------


! Loop over all atoms:

        if (Mass == 1) then
	  do k=1,Natoms
      		!incompatible priors for M & Z?
		if( ( MassModel /= 1 .and. ( z_PriorType(k) == 8 .or. Mass_PriorType(k,2) == 8 ) ) .or. &
      		( MassModel == 1 .and. NFWstyle(k) /= 2 .and. ( z_PriorType(k) == 8 .or. Mass_PriorType(k,2) == 8 ) ) &
            	.or. ( MassModel == 1 .and. NFWstyle(k) == 2 .and. ( ( z_PriorType(k) == 8 .and. &
            	Mass_PriorType(k,2) /= 8 ) .or. ( z_PriorType(k) /= 8 .and. Mass_PriorType(k,2) == 8 ) ) ) ) then
			call halt_program("incompatible priors used for M200 & z")
		endif
     	  enddo 
        endif

	!allocate memory if mass function priors used for both M & z
	if( ( MassModel == 1 ) .or. ( ((SZ == 1) .or. (PL == 1)) .and. Gas == 1 ) ) then
	    Mmin = 1.d90
	    Mmax = 0d0
	    do k = 1, NAtoms
		if( (Mass==1) .and. (MassModel==1) ) then
		    if( NFWstyle(k) == 2 .and. z_PriorType(k) == 8 .and. Mass_PriorType(k,2) == 8 ) then
			if( Mass_Prior(k,2,1) < 0d0 ) then
			    call halt_program("Incorrect lower limit on M200 prior. Aborting.")
			endif
			if( Mass_Prior(k,2,1) < Mmin ) Mmin = Mass_Prior(k,2,1)
			if( Mass_Prior(k,2,2) > Mmax ) Mmax = Mass_Prior(k,2,2)
		    endif
		endif
		if( (BetaStyle(k) == 4 .or. GasModel==5 .or. GasModel==6) .and. z_PriorType(k) == 8) then !adapted for GM=6 kj 01/02/17
                    if (BetaStyle(k)==4) M200_par=3
                    if (GasModel==5 .or. GasModel==6) M200_par=1
                    if (Gas_PriorType(k,M200_par) == 8 ) then
		      if( Gas_Prior(k,M200_par,1) < 0d0 ) then
			call halt_program("Incorrect lower limit on M200 prior. Aborting.")
                      endif
		    endif
		    if( Gas_Prior(k,M200_par,1) < Mmin ) Mmin = Gas_Prior(k,M200_par,1)
		    if( Gas_Prior(k,M200_par,2) > Mmax ) Mmax = Gas_Prior(k,M200_par,2)
		endif

	    enddo
                  
	    if( Mmax > 0d0 ) then
		allocate(lookM(n,2),lookZ(n,n,2))
		call makeMZlookup(Mmin, Mmax, zdmin, zdmax)
		write(*,*)"MZlookup done"
	    endif
	endif
      		
end subroutine CheckMassPriors
!=======================================================================
! Sanity check for triangle priors, and set them up if present and sensible

subroutine CheckGeoPriors

    implicit none
    integer  i, j, k, idx(1)

    do i=1,NAtoms
	! uniform sampling in triangle check
	if( ( Geo_PriorType(i,1) == 9 .and. Geo_PriorType(i,2) /= 9 ) .or. &
	( Geo_PriorType(i,1) /= 9 .and. Geo_PriorType(i,2) == 9 ) ) then
		call halt_program("Geo_PriorType should be set to 9 for both x and y positions if &
				uniform sampling in a triangle is desired")
	else if( Geo_PriorType(i,1) == 9 .and. Geo_PriorType(i,2) == 9 ) then
	    do j = 1, 2
		do k = j + 1, 3
		    if( Geo_Prior(i,j,1) == Geo_Prior(i,k,1) .and. Geo_Prior(i,j,2) == Geo_Prior(i,k,2) ) then
			call halt_program("two vertices of the triangle can not be same")
		    endif
		enddo
	    enddo
	    do j = 1, 2
		idx = minloc(Geo_Prior(i,j,1:3))
		Geo_Tri_Prior(i,j,1) = Geo_Prior(i,j,idx(1))
		idx = maxloc(Geo_Prior(i,j,1:3))
		Geo_Tri_Prior(i,j,2) = Geo_Prior(i,j,idx(1))
  	   enddo
	endif
    enddo

end subroutine CheckGeoPriors
!=======================================================================

subroutine CheckSrcPriors

    implicit none
    integer  i

    if (NSrc > 0) then
      ! Switch Gaussian priors on source fluxes to truncated Gaussian priors
      do i = 1, NSrc
        if (Src_PriorType(i,3) == 3) then
          Src_PriorType(i,3) = 12
          Src_Prior(i,3,3) = 0.
          if (verbose.and.is_root) write(*,*) 'Src_PriorType(i,3) is now 12 for source ', i
        endif
      enddo
    endif

end subroutine CheckSrcPriors

!=======================================================================
! Some sanity checks for running Planck data
subroutine CheckPlanckPriors

     implicit none

     if (GasModel.ne.3 .and. GasModel.ne.5 .and. GasModel.ne.6) then
       write(*,*) 'Can only currently handle GasModel = 3/5/6 with Planck data'
       stop
     endif
     if (NAtoms.gt.1) then
       write(*,*) 'Can only currently handle NAtoms = 1 with Planck data'
       stop
     endif
     if (NGeoPars.gt.2) then
       write(*,*) 'Can only currently model spherically symmetric clusters with Planck data'
       stop
     endif

end subroutine CheckPlanckPriors

!=======================================================================

end module CheckPars1
