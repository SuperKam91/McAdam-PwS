module cosmology

	use constants
      use utilities
      
contains

!=======================================================================

      function Hubble(zloc)

      	implicit none
      	
      	double precision Hubble,zloc,Q

      	if (wa /= 0.) then
        	Q = (1.+zloc)**(3.*(w0+wa))*exp(-3.*wa*zloc/(1.+zloc))
      	else
        	Q = (1.+zloc)**(3.*w0)
      	endif
	
! 	General expression:
      
      	Hubble = (1.+zloc)*(1.+zloc)*(1.+zloc)*(Om + OL*Q) + Ok*(1.+zloc)*(1.+zloc)
      	Hubble = h*100.*sqrt(Hubble)

	end function Hubble
     
!=======================================================================
! Function returns rhocrit(zz)

      function rhocritofz(zz)

      	implicit none
      
      	double precision rhocritofz,zz
      
      	if (zz == 0.) then
        	rhocritofz = rhocrit
      	else
        	rhocritofz = rhocrit*(Hubble(zz)/(h*100.))**2.
	endif
            
      end function rhocritofz

!=======================================================================
! Function returns rhobar(zz)

      function rhobarofz(zz)

      	implicit none
      
      	double precision rhobarofz,zz
      
      	rhobarofz = Omofz(Om,zz)*rhocritofz(zz)
            
      end function rhobarofz

!=======================================================================
! Function returns Omega_m(zz)

      function Omofz(Omega_m,zz)

      	implicit none
      
      	double precision Omofz,Omega_m,zz
      
      	if (zz == 0.) then
        	Omofz = Omega_m
      	else  
        	Omofz = Omega_m*(1.+zz)**3./(1.-Omega_m+Omega_m*(1.+zz)**3.)
	endif
            
      end function Omofz

!=======================================================================
! Function the (angular diameter distance) * h as a function of redshift z

      function calcAngDiamDis(z1,z2)

      	implicit none
	
      	double precision calcAngDiamDis,z1,z2,S,lim1,lim2
      	double precision eps
      	parameter(eps=1d-3)
      
      	lim1=1./(1.+z2)
      	lim2=1./(1.+z1)
      	call qtrap(AngDiamInt,lim1,lim2,eps,S)
      	calcAngDiamDis=clight*S/(h*100000.*(1.+z2))

	end function calcAngDiamDis
	
!=======================================================================
      
      function AngDiamInt(r)
      
      	implicit none
      
      	double precision r,AngDiamInt
      
      	AngDiamInt=1./sqrt(Om*r+OL*r**(1.-3.*w0))
      
      end function AngDiamInt
	
!=======================================================================
! Comoving distance in h^-1 Mpc:

	function rcomoving(zz)
      
      	implicit none
      
      	double precision rcomoving,zz
      	double precision eps
      	parameter(eps=1d-4)
      
      	call qtrap(conH,0.d0,zz,eps,rcomoving)
            
      end function rcomoving

!=======================================================================

	function conH(zz)
      
      	implicit none
      
      	double precision conH,zz
      
      	if (zz == 0.0) then
        	conH = 2.998d5/100.0
      	else
        	conH = 2.998d5/Hubble(zz)
		endif
            
      end function conH

!=======================================================================
! Return lookback time in Gyr:
	function lookbacktime(zz)
      
     	 	implicit none
      
      	double precision lookbacktime,zz
      	double precision eps
      	parameter(eps=1d-3)
      
      	call qtrap(timeintegrand,0.d0,zz,eps,lookbacktime)
      
      end function lookbacktime

!=======================================================================

	function timeintegrand(zz)
      
      	implicit none
      
      	double precision timeintegrand,zz
      
      	timeintegrand = 100./(0.1021*Hubble(zz)*(1.+zz))
      
      end function timeintegrand

!=======================================================================

	function Scurvature(r)
      
      	implicit none
      
      	double precision Scurvature,r
      	integer k
      
		if (Ok == 0.) then
	  		k = 0
		elseif (Ok > 0.) then
	  		k = -1
		elseif (Ok < 0.) then
	  		k = 1
		endif 
      
		if (k == -1) then
	  		Scurvature = sinh(r)
		elseif (k == 0) then
	  		Scurvature = r
		elseif (k == 1) then
	  		Scurvature = sin(r)
		endif  
      
      end function Scurvature

!=======================================================================

	function angdist(zz)
      
      	implicit none
      
      	double precision angdist,zz
      	double precision r,R0
      
	      if (Ok == 0.) then
      	  	R0 = conH(0.d0)
      	else
        		R0 = conH(0.d0)/sqrt(abs(Ok))
      	endif
      	r = rcomoving(zz)
      	angdist = R0*Scurvature(r/R0)/(1.+zz)
      
      end function angdist

!=======================================================================

	function angdistdiff(z1,z2)
      
      	implicit none
      
	      double precision angdistdiff,z1,z2
      	double precision r1,r2,R0
      
	      if (Ok == 0.0) then
        		R0 = conH(0.d0)
      	else
        		R0 = conH(0.d0)/sqrt(abs(Ok))
      	endif
      	r1 = rcomoving(z1)
      	r2 = rcomoving(z2)
      	angdistdiff = R0*Scurvature((r2-r1)/R0)/(1.+z2)
      
      end function angdistdiff

!=======================================================================

	function criticaldensity(z1,z2)
      
      	implicit none
            
      	double precision criticaldensity,z1,z2
      
      	criticaldensity = angdist(z2)/(angdist(z1)*angdistdiff(z1,z2))
      	criticaldensity = criticaldensity*1.6620d6
      
      end function criticaldensity

!=======================================================================

	function lumdist(zz)
      
      	implicit none
      
      	double precision lumdist,zz
      	double precision R0
      
      	if (Ok == 0.) then
        		R0 = conH(0.d0)
      	else
        		R0 = conH(0.d0)/sqrt(abs(Ok))
      	endif
      	lumdist = R0*Scurvature(rcomoving(zz)/R0)*(1.+zz)
      
      end function lumdist

!=======================================================================
! 	Differential comoving volume for unit solid angle:

	function dVdzcomoving(zz)
      
      	implicit none
      
      	double precision dVdzcomoving,zz
      	double precision S,r,R0
      
      	if (Ok == 0.0) then
        		R0 = conH(0.d0)
      	else
        		R0 = conH(0.d0)/sqrt(abs(Ok))
      	endif
      	r = rcomoving(zz)
      	S = R0*Scurvature(r/R0)
      	dVdzcomoving = conH(zz)*S*S
      
      end function dVdzcomoving

!=======================================================================
! 	Comoving volume between 2 redshifts for unit solid angle:

	subroutine Vcomoving(z1,z2,Volume)
      
      	implicit none
      
      	double precision z1,z2,Volume
      	double precision eps
      	parameter(eps=1d-3)
      
      	call qtrap(dVdzcomoving,z1,z2,eps,Volume)
      
      end subroutine Vcomoving

!=======================================================================

end module cosmology
