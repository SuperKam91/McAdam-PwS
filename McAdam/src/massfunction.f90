module massfunction

	use cosmology
	use constants
	use params
      
      	double precision zglob,MM,Mglob,Minfglob
      
	contains

!=======================================================================

function PSfunc(M,zloc)

	!The Evrard approximation to the halo mass function, as described in
	!Allen et al 2002. Function returns dn/dMdz in h^4 M_sun^-1 Mpc^-3. Flat
	!cosmology only at present.

      	implicit none
      
      	double precision M,zloc,PSfunc
	
	double precision Omz,g0,gz,sigma8z,sigmaz,rhobar,rhobarz,R
	double precision BigGamma,LittleGamma,x,A,B,c,eps,f,dndM
	
!-----------------------------------------------------------------------
	
	! First find density parameter at redshift zloc:

	Omz = Omofz(Om,zloc)
      
	! Now find growth factors at redshift 0 and zloc:

	g0 = growth(0.d0)
	gz = growth(zloc)

	! Now find sigma8 at redshift zloc:

	sigma8z = sigma8*(gz/g0)/(1.+zloc)

	! Now need the radius which corresponds to M...

 	rhobar = Om*rhocrit
	R = (3.*M/(4.*pi*rhobar))**(1./3.)

	! ...in order to calculate shape parameters:

      	BigGamma = Om*h*(2.7/2.726)*(2.7/2.726)*exp(-Ob-sqrt(h/0.5)*Ob/Om)
	LittleGamma = (0.3*BigGamma+0.2)*(2.92+log10(R*h/8.))

	! So now can calculate pieces to go in Evrard formula:

	rhobarz = rhobar*(1.+zloc)**3
	!Better to have comoving number density:
	!rhobarz = rhobar
	sigmaz = sigma8z*(R*h/8.)**(-LittleGamma)

	
	! Interpolations for A B and eps:
	
	if( mass_function == 1 .or. mass_function == 2 ) then
		! Evrard approximation to Press-Schechter
		if( mass_function == 1 ) then
			x = (1.-Omz)/0.7
	      		A = (1.-x)*0.27 + x*0.22
      			B = (1.-x)*0.65 + x*0.73
      			eps = (1.-x)*3.77 + x*3.86
		! Jenkins mass function
		elseif( mass_function == 2 ) then
			A = 0.315
			B = 0.61
			eps = 3.8
		endif
		f = A*exp(-1.*abs(log(1./sigmaz)+B)**eps)
	!Tinker mass function
	elseif( mass_function == 3 ) then
		A = 0.1d0 * log10(200d0) - 0.05d0
		b = 1d0 + ( log10(200d0) - 1.6d0 )**(-1.5d0)
		c = 1.2d0 - ( -log10(200d0) + 2.35d0 )**1.6d0
		eps = 1.43d0 + ( log10(200d0) - 2.3d0 )**1.5d0
		f = A * ( ( ( sigmaz / b )**(-eps) ) + 1d0 ) * exp( -c / ( sigmaz**2d0 ) )
	else
		write(*,*)"Incorrect value set for variable mass_function in the include file"
		stop
	endif

	dndM = LittleGamma*rhobarz*f/(3.*M*M)

	! Done!

	PSfunc = dndM	
	
end function PSfunc

!=======================================================================
!Function returns linear growth factor g(zloc), from CPT92

function growth(zloc)

      	implicit none
      
      	double precision growth,zloc,Omz
      
      	Omz = Omofz(Om,zloc)
     	growth = (5./2.)*Omz/(1./70. + 209.*Omz/140. - Omz**2/140. + Omz**(4./7.))
      
end function growth

!=======================================================================
!Function returns dn(logM,zloc)/dlogM in h^3 Mpc^-3

function PSlogfunc(logM,zloc)

      	implicit none
      
      	double precision PSlogfunc,logM,zloc
      	double precision M
      
      	M = exp(logM)
      	PSlogfunc = M*PSfunc(M,zloc)
            
end function PSlogfunc

!=======================================================================
!Integrate over masses to get n(>M,zloc) in h^3 Mpc^-3

subroutine ngtMofz(M,Minfinity,zloc,Result)
      
      	implicit none
      
      	double precision M,zloc,Result
      	double precision Minfinity,check1,check2
      	!parameter(Minfinity=1.d17)
      	double precision eps
      	parameter(eps=1d-3)
      
      	zglob = zloc
      
      	check1 = ngtMofzIntegrand(log(M))
      	check2 = ngtMofzIntegrand(log(Minfinity))
      	if (check1 < 1d-12 .and. check2 < 1d-12) then
        	Result = 0.
      	else
        	call qtrap(ngtMofzIntegrand,log(M),log(Minfinity),eps,Result)
      	endif
      
end subroutine ngtMofz

!-----------------------------------------------------------------------

function ngtMofzIntegrand(logM)
      
      	implicit none
      
      	double precision ngtMofzIntegrand,logM
      
      	ngtMofzIntegrand = PSlogfunc(logM,zglob)
      
end function ngtMofzIntegrand

!=======================================================================
!Differential  dN(M,<zloc)/dlogM in Msun^-1 sr^-1

subroutine dNltzdM(M,zloc,Result)
      
      	implicit none
      
      	double precision M,zloc,Result
      	double precision eps
      	parameter(eps=1d-3)
      
      	MM = log(M)
      	call qtrap(dNltzdMIntegrand,0.d0,zloc,eps,Result)
      
end subroutine dNltzdM

!-----------------------------------------------------------------------

function dNltzdMIntegrand(zloc)
      
      	implicit none
      
      	double precision zloc,dNltzdMIntegrand 
      
      	dNltzdMIntegrand = PSlogfunc(MM,zloc)
      	dNltzdMIntegrand = dNltzdMIntegrand*dVdzcomoving(zloc)
		
end function dNltzdMIntegrand

!=======================================================================
!Differential  dN(>M,zloc)/dz in sr^-1

subroutine dNgtMdz(M,Minfinity,zloc,Result)
      
      	implicit none
      
	double precision M,Minfinity,zloc,Result
      
      	call ngtMofz(M,Minfinity,zloc,Result)
      
      	Result = Result*dVdzcomoving(zloc)
            
end subroutine dNgtMdz

!=======================================================================
!Integral  N(>M,<zloc) in sr^-1

subroutine NgtMltz(M,Minfinity,zloc,Result)
      
      	implicit none
      
      	double precision M,Minfinity,zloc,Result
      	double precision eps
      	parameter(eps=1d-3)
      
      	MM=M
        Minfglob=Minfinity
            
      	call qtrap(NgtMltzIntegrand,0.d0,zloc,eps,Result)
            
end subroutine NgtMltz

!-----------------------------------------------------------------------

function NgtMltzIntegrand(zloc)

	implicit none

	double precision NgtMltzIntegrand,zloc
	double precision Minfinity

	Minfinity=Minfglob

	call dNgtMdz(MM,Minfinity,zloc,NgtMltzIntegrand)
      
end function NgtMltzIntegrand

!-----------------------------------------------------------------------
subroutine makeMZlookup(Mmin, Mmax, zmin, zmax)
      
      	implicit none
            
	integer i,j
	double precision Mmin,Mmax,Minc,zmin,zmax,zinc
	double precision normM,normz
            
!	Mmin=Mass_Prior(1,2,1) !smallest mass (prior)
!      	Mmax=Mass_Prior(1,2,2) !largest mass (prior)
!      	zmin=zdmin !smallest redshift (prior)
!      	zmax=zdmax !largest redshift (prior)
      
	!set up the step size
	!n is the no. of steps, generally set to 512
	Minc=(Mmax-Mmin)/n
      	zinc=(zmax-zmin)/n
      
	do i=1,n
            	!look-up table for M & its cumulative probability
		!lookM(i,1)=Mmin+Minc/2.+Minc*(i-1.)
		if(i==n) then
                  	lookM(i,1)=Mmax
		else
                  	lookM(i,1)=10.**(log10(Mmin)+(log10(Mmax)-log10(Mmin))*(i-1.)/(n-1.))
		endif
                  
		!calculate the cumulative Pr(M)
		if(i==1) then
			lookM(i,2)=cumPrM(log(Mmin),log(lookM(i,1)))
		else
                  	lookM(i,2)=lookM(i-1,2)+cumPrM(log(lookM(i-1,1)),log(lookM(i,1)))
		endif
                  
		!now calculate the cumulative probability of z given M
		!set Mglob so that M is fixed
		Mglob=log(lookM(i,1))
                  
		!iterate over the z values from zmin to zmax
            	do j=1,n
            		!look-up table for M & its cumulative probability
                  	!lookZ(i,j,1)=zmin+zinc/2.+zinc*(j-1.)
                        if(j==n) then
                  		lookZ(i,j,1)=zmax
			else
                              	lookZ(i,j,1)=10.**(log10(zmin)+(log10(zmax)-log10(zmin))*(j-1.)/(n-1.))
			endif
                        
                  	!calculate the cumulative Pr(z)
                        if(j==1) then
				if( zmin == zmax ) then
					lookZ(i,j,2) = 0d0
				else
					lookZ(i,j,2)=cumPrz(zmin,lookZ(i,j,1))
				endif
			else
				if( zmin == zmax ) then
					lookZ(i,j,2) = 1d0
				else
                        		lookZ(i,j,2)=lookZ(i,j-1,2)+cumPrz(lookZ(i,j-1,1),lookZ(i,j,1))
				endif
			end if
		end do
                  
            	!normalizing constants for the probabilities
            	normZ=lookZ(i,n,2)+cumPrz(lookZ(i,n,1),zmax)
                  
		!normalize the cumulative probabilities of z given M                  
		lookZ(i,1:n,2)=lookZ(i,1:n,2)/normz
	end do
		
	!normalizing constants for the probabilities
	normM=lookM(n,2)+cumPrM(log(lookM(n,1)),log(Mmax))
            
	!normalize the cumulative probabilities of M
	lookM(1:n,2)=lookM(1:n,2)/normM
                        
end subroutine makeMZlookup

!-----------------------------------------------------------------------
!dNdMdz returns the dN/dlogMdz
function dNdMdz(logM,zz)
      
      	implicit none
            
	double precision dNdMdz,logM,zz
            
	dNdMdz=PSlogfunc(logM,zz)*dVdzcomoving(zz)
	
end function dNdMdz

!-----------------------------------------------------------------------
!dNdMdz returns the dN/dlogMdz with the M value taken as the global one & z passed
function dNdMdz1(zz)
      
      	implicit none
            
	double precision dNdMdz1,zz,logM
            
	logM = Mglob
	dNdMdz1=dNdMdz(logM,zz)
	
end function dNdMdz1

!-----------------------------------------------------------------------
!Prz returns the un-normalized probability of z given M which is equal to dN/DMdz
!z is passed as an argument & M is taken to be the global variable Mglob
function Prz(zz)
      
      	implicit none
            
	double precision Prz,zz,logM
            
	logM = Mglob
	Prz=dNdMdz(logM,zz)
	
end function Prz

!-----------------------------------------------------------------------
!PrM returns the un-normalized cumulative probability of z given M
function cumPrz(zmin,zz)
      
      	implicit none
            
	double precision zmin,zz,cumPrz
	double precision eps
	parameter(eps=1d-3)
            
	call qtrap(Prz,zmin,zz,eps,cumPrz)
            
end function cumPrz

!-----------------------------------------------------------------------
!PrM returns the un-normalized probability of M by marginalizing dN/DMdz over z
function PrM(logM)
      
      	implicit none
            
	double precision logM,zmin,zmax,PrM
	double precision eps
	parameter(eps=1d-3)
            
	!get the min & max redshift values from the prior
	zmin=zdmin
	zmax=zdmax
            
	!set the global M value, so that qtrap calls a function with takes only 1 variable
	Mglob=logM
        
	if( zmin == zmax ) then
		PrM = dNdMdz1(zmin)
	else
		call qtrap(dNdMdz1,zmin,zmax,eps,PrM)
	endif
            
end function PrM

!-----------------------------------------------------------------------
!PrM returns the un-normalized cumulative probability of M
function cumPrM(logMmin,logM)
      
      	implicit none
            
	double precision logMmin,logM,cumPrM
	double precision eps
	parameter(eps=1d-3)
            
	call qtrap(PrM,logMmin,logM,eps,cumPrM)
            
end function cumPrM

!-----------------------------------------------------------------------

subroutine lookUpMZ(cPM,cPZ,M,zz)
      
      	implicit none
      
      	double precision cPM !cumulative probability of M
	double precision cPZ !cumulative probability of z given M
	double precision M !M corresponding to cPM
	double precision zz !z corresponding to cPZ
	integer i,j,k
	double precision zz1,zz2
      	
	!Lookup the Mass table
	i=binSearch(n,lookM(:,2),cPM)
            
	!Lookup the z table for M>M_req
	j=binSearch(n,lookZ(i,:,2),cPZ)
	
	!Lookup the z table for M<M_req
	k=binSearch(n,lookZ(i-1,:,2),cPZ)
            
	!sanity check
	if(i==0 .or. j==0 .or. k==0 .or. i>n .or. j>n .or. k>n) then
            	write(*,*)"problem in lookUpMZ, value to look for is not within the table bounds"
	endif
            
	!interpolation to get the z values from z tables
	zz1=lookZ(i,j-1,1)+(lookZ(i,j,1)-lookZ(i,j-1,1))*(cPZ-lookZ(i,j-1,2))/(lookZ(i,j,2)-lookZ(i,j-1,2))
	zz2=lookZ(i+1,k-1,1)+(lookZ(i+1,k,1)-lookZ(i+1,k-1,1))*(cPZ-lookZ(i+1,k-1,2))/(lookZ(i+1,k,2)-lookZ(i+1,k-1,2))
            
	!interpolation between the z values for M<M_req & M>M_req
	zz=zz1+(zz2-zz1)*(cPM-lookM(i-1,2))/(lookM(i,2)-lookM(i-1,2))
	M=lookM(i-1,1)+(lookM(i,1)-lookM(i-1,1))*(cPM-lookM(i-1,2))/(lookM(i,2)-lookM(i-1,2))
      
end subroutine lookUpMZ

!-----------------------------------------------------------------------

subroutine lookUpM(cPM,M)
      
      	implicit none
      
      	double precision cPM !cumulative probability of M
	double precision M !M corresponding to cPM
	integer i
      	
	!Lookup the Mass table
	i=binSearch(n,lookM(:,2),cPM)
            
	!sanity check
	if(i==0 .or. i>n) then
            	write(*,*)"problem in lookUpM, value to look for is not within the table bounds"
	endif
            
	M=lookM(i-1,1)+(lookM(i,1)-lookM(i-1,1))*(cPM-lookM(i-1,2))/(lookM(i,2)-lookM(i-1,2))
      
end subroutine lookUpM

!-----------------------------------------------------------------------

subroutine lookUpZ(cPZ,M,zz)
      
      	implicit none
      
	double precision cPZ !cumulative probability of z given M
	double precision M !M corresponding to cPM
	double precision zz !z corresponding to cPZ
	integer i,j,k
	double precision zz1,zz2
      	
	!Lookup the Mass table
	i=binSearch(n,lookM(:,1),M)
            
	!Lookup the z table for M>M_req
	j=binSearch(n,lookZ(i,:,2),cPZ)
	
	!Lookup the z table for M<M_req
	k=binSearch(n,lookZ(i-1,:,2),cPZ)
            
	!sanity check
	if(i==0 .or. j==0 .or. k==0 .or. i>n .or. j>n .or. k>n) then
            	write(*,*)"problem in lookUpZ, value to look for is not within the table bounds"
	endif
            
	!interpolation to get the z values from z tables
	zz1=lookZ(i,j-1,1)+(lookZ(i,j,1)-lookZ(i,j-1,1))*(cPZ-lookZ(i,j-1,2))/(lookZ(i,j,2)-lookZ(i,j-1,2))
	zz2=lookZ(i+1,k-1,1)+(lookZ(i+1,k,1)-lookZ(i+1,k-1,1))*(cPZ-lookZ(i+1,k-1,2))/(lookZ(i+1,k,2)-lookZ(i+1,k-1,2))
            
	!interpolation between the z values for M<M_req & M>M_req
	zz=zz1+(zz2-zz1)*(M-lookM(i-1,1))/(lookM(i,1)-lookM(i-1,1))
      
end subroutine lookUpZ

!-----------------------------------------------------------------------

!=======================================================================

end module massfunction
