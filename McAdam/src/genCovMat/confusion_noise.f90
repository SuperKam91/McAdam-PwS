module kind_def

	integer, parameter :: dp = kind(1.0D0)
  	integer, parameter :: rk = SELECTED_REAL_KIND(12,200)
  	!integer :: num_params  = 2

end module kind_def

module consts
	
	use kind_def

   	implicit none

   	double precision, parameter :: clight = 299792458.d0
   	double precision, parameter :: kboltzmann = 1.380658e-23

end module consts

module conf_globals

   	use kind_def

   	implicit none

   	double precision a,b

end module conf_globals

module aux

   	!use kind_def
   	use consts

   	implicit none

contains

!=======================================================================       
! Numerical recipes integration routine-the supplied real function
! func is integrated between limits A and B, and the result returned as
! the integral S. Progressively finer discretisations are used until 
! either
!
!     | S_{j}-S_{j-1} | 
!     ------------------- < eps
!          | S_{j-1} |
!
! or the range has to be split into more than 2^jmax sections. For
! smooth functions the latter shouldnt happen...

subroutine qtrap(func,lim1,lim2,S)
      
	implicit none
	
	double precision S,lim1,lim2
	double precision a1,b1,eps,olds,sign,t1,t2
	integer j,jmax
	parameter (eps=1.e-3,jmax=25)
      
      	INTERFACE
		function func(zz)
            		double precision zz,func
            	end function func
      	end INTERFACE

!-----------------------------------

	a1=lim1
	b1=lim2
	
	! First catch some stupidities:

	if (a1.eq.b1) then
	  	S=0.0
	  	goto 30
	elseif (a1.gt.b1) then        
	  	olds=b1
	  	b1=a1
	  	a1=olds
	  	sign=-1.0
	else 
	  	sign=1.0
	endif  
	 
	olds=-1.e30
	s=0.0
      	do j=1,jmax
        	call trapzd(func,a1,b1,S,j)
	  	if (j.eq.1) then
	    		t1=func(a1)
	    		t2=func(b1)
	  	endif  
		!write(*,*) '   QTRAP: j,S,a1,f(a1),b1,f(b1)=',j,S,a1,t1,b1,t2
		!write(*,*) '   QTRAP: j,S=',j,S
        	if (abs(s-olds).lt.eps*abs(olds).or.abs(s-olds).eq.0.0) goto 20
        	if (j.eq.jmax) goto 10
        	if (abs(s).lt.1d-45) then
	    		s=0.0
	    		goto 30
        	endif
	  	olds=s
      	enddo
 10   	write(*,*) 'QTRAP error: too many steps...'
      	write(*,*) '   S=',S
      	write(*,*) '   oldS=',oldS
      	write(*,*) '   % difference=',100.0*(S-oldS)/oldS
      	write(*,*) '   limits were ',a1,b1
      	write(*,*) '   function values at limits were ',t1,t2
	t1=func(0.5*(a1+b1))
      	write(*,*) '   function value at midpoint was ',t1
	S=2e10
	goto 30
!	pause 
      
 20   	S=S*sign
 
 30   	return
end subroutine qtrap

!-----------------------------------------------------------------------
      
subroutine trapzd(func,a1,b1,S,n)
	
	implicit none
	
	double precision a1,b1,S
	integer n
	integer it,j
	double precision tnm,del,x,sum
      
      	INTERFACE
            	function func(zz)
            		double precision zz,func
            	end function func
      	end INTERFACE
	
      	if (n.eq.1) then
        	S=0.5*(b1-a1)*(func(a1)+func(b1))
        	it=1
     	else
	  	it=2**(n-2)
        	tnm=it*1.0
        	del=(b1-a1)/tnm
        	x=a1+0.5*del
        	sum=0.
        	do 11 j=1,it
          		sum=sum+func(x)
          		x=x+del
11      	continue
        	s=0.5*(s+(b1-a1)*sum/tnm)
      	endif
      
	return
end subroutine trapzd
	
!=======================================================================

function integrate(y,x,n,xmin,xmax)
!
!  Numerical integration of function y(x): y and x are arrays of
!  length n, and y is integrated between xmin and ymax using the
!  trapezium rule.
!
      	implicit none

      	integer n,i,j,k
      	double precision y(n),x(n)
      	double precision sum,integrate
      	double precision xmin,xmax,ymin,ymax
!
!-----------------------------------------------------------------------
!
      	sum=0.0
!
      	do i=1,n
      		if (xmin.lt.x(i)) goto 10
      	enddo
!
 10   	j=i
      	ymin=y(j)+(y(j-1)-y(j))*(x(j)-xmin)/(x(j)-x(j-1))
      	sum=(ymin+y(j))*(x(j)-xmin)/2.0
!
      	do i=j,n
      		if (xmax.le.x(i)) goto 20
      	enddo
!
 20   	k=i-1
      	ymax=y(k)-(y(k)-y(k+1))*(xmax-x(k))/(x(k+1)-x(k))
      	sum=sum+(ymax+y(k))*(xmax-x(k))/2.0   
!       
      	do i=j,k-1
      		sum=sum+(y(i+1)+y(i))*(x(i+1)-x(i))/2.0
      	enddo
!
      	integrate=sum
!
      	return
end function integrate
!
!=======================================================================

end module aux

!=======================================================================

module my_types

	use kind_def

   	implicit none

   	type telescope_beam
      		integer id             ! Tag
      		character*8 telescope  ! Telescope name
      		double precision synth_beam_sq_arcsec ! Synthesized beam in sq. arcsec
      		double precision synth_beam_sq_arcmin ! Synthesized beam in sq. arcmin
      		double precision synth_beam_sr ! Synthesized beam in sr
   	end type telescope_beam

   	type power_law
      		integer id        ! Tag
      		double precision a        ! Normalisation
      		double precision b        ! Index
      		double precision slow     ! Lower flux limit (Jy) below which not to be used
      		character*16 name ! Reference
   	end type power_law

end module my_types

!=======================================================================

module count_models

   	use kind_def
   	use my_types

   	implicit none

   	type(power_law) :: count_model(3)

contains

subroutine init_models()

	count_model(1) = power_law(1,80.0d0,-2.0d0,10.0d-6,'Taylor01')
      	count_model(2) = power_law(2,51.0d0,-2.15d0,10.0d-6,'Waldram03')
      	count_model(3) = power_law(3,376.0d0,-1.8d0,10.0d-6,'Cam10C')

end subroutine init_models

end module count_models

!=======================================================================


module telescopes

   	use kind_def
   	use my_types

   	implicit none

   	type(telescope_beam) tel_synth(4)

contains

subroutine init_telescopes()

    	tel_synth(1) = telescope_beam(1,'SA',240.0d0**2,4.0d0**2,1.4d-6)
    	tel_synth(2) = telescope_beam(2,'LA',30.0d0**2,0.5d0**2,2.1d-8)
    	tel_synth(3) = telescope_beam(3,'RT',25.0d0*19.0d0,0.417d0*0.317d0,1.1d-8)
    	tel_synth(4) = telescope_beam(4,'VSA-SEX',324.0d0*402.0d0,5.40d0*6.70d0,3.1d-6)

end subroutine init_telescopes

end module telescopes

!=======================================================================


module confusion_noise

! A module to generate the confusion noise given a telescope,
! a source count model and a limiting flux.
! v0.1 JZ 070903
! v0.2 JZ 070919 updated to output C_Conf as per MPH's email 070904

use kind_def
use aux

   	implicit none

contains

!-----------------------------------------------------------------------

function integrand(flux)

      	use conf_globals, only : a,b

      	implicit none

      	double precision integrand,flux
	
      	!integrand = diff_count*(flux**2)
      	integrand = a*flux**(2.0d0+b)

end function integrand

!-----------------------------------------------------------------------

function db_dt(freq)

      	use consts

      	implicit none

      	double precision db_dt
      	double precision freq
      	double precision one_over_lambda

      	one_over_lambda = freq/clight

	! ... db_dt is in units of Jy K^{-1} sr^{-1}

	! ... NB T is antenna temperature not thermodynamic temperature!!

      	db_dt = 2*kboltzmann*one_over_lambda*one_over_lambda*1.0d26
	!      write(*,*) db_dt

end function db_dt

!-----------------------------------------------------------------------

function conf_var(flux_low,flux_lim,diff_count_type)

      	use conf_globals, only : a,b
      	implicit none

      	double precision conf_var
      	integer diff_count_type

      	double precision flux_low,flux_lim

	! NB Lower limit has been fixed at 1.0e-20 because blows up if 0.0e0!

      	!call qtrap(integrand,flux_low,flux_lim,integral)

      	!conf_var = dble(integral)
	!write(*,*) conf_var
	
	
	conf_var = dble( a * ( flux_lim**( 3.0 + b ) ) / ( 3.0 + b ) )

end function conf_var

!-----*-----------------------------------------------------------------

!calculates conversion factor between equivalent CMBR thermodynamic
!temperature fluctuation and antenna temperature fluctuation at freq

function dtt_to_dta(tcmb,freq)

      	implicit none

      	double precision dtt_to_dta
      	double precision tcmb,freq

      	double precision x
      	double precision, parameter :: honk=0.0479927 ! this is h/k for frequency in GHz

      	x=honk*freq/tcmb
      	dtt_to_dta=(x**2)*exp(x)/((exp(x)-1)**2)

end function dtt_to_dta


!-----------------------------------------------------------------------

subroutine calc_c_conf(diff_count_type,freq,slim,c_conf)

      	use my_types
      	use count_models
      	use conf_globals, only : a,b
	use globals, only : tcmb, PI
      	implicit none

      	double precision slim, slow
      	double precision freq,freqloc
      	integer diff_count_type
      	double precision c_conf, c_conf_ta, factor
      	logical, parameter :: verbose = .false.

	! ... Retrieve source count models
      	call init_models()
      	a=count_model(diff_count_type)%a
      	b=count_model(diff_count_type)%b

	! ... Convert flux mJy -> Jy
      	slim=slim/1000
      	slow=count_model(diff_count_type)%slow

	! ... Calculate T_A -> T_Th conversion factor of freq/GHz
      	factor=1.0d0/dtt_to_dta(tcmb,freq)

	! ... Convert frequency GHz -> Hz
      	freqloc=freq*1.0d9

      	c_conf=dble(conf_var(dble(slow),slim,diff_count_type)/ &
           (db_dt(freqloc)*db_dt(freqloc)))
      	
	c_conf_ta=c_conf

	! ... Convert c_conf K^2 -> \mu K^2
      	!c_conf=c_conf*1.0d12
      	!c_conf_ta=c_conf

	! ... Convert c_conf \mu K^2 : antenna T to thermodynamic T
      	c_conf=c_conf*factor*factor

	! ... make it dimensionless & divide by 2\pi
      	c_conf=c_conf/(2d0*PI*tcmb**2)

      	if(verbose) then
         	write(*,*)
         	!write(*,*) 'SUMMARY (Verbose =',verbose,'):'
         	write(*,*) 'frequency/Hz= ', freqloc
         	write(*,*) 'flux_low/Jy= ', slow
         	write(*,*) 'flux_lim/Jy= ', slim
         	write(*,*) 'Model ', diff_count_type
         	write(*,*) 'a= ', a
         	write(*,*) 'b= ', b
         	write(*,*) 'c_conf/uK^2 (T_antenna)= ', c_conf_ta
         	write(*,*) 'c_conf/uK^2 (T_thermo)= ', c_conf
         	write(*,*)
      	endif
	
	
!	open(unit=99, file='conf.dat', status='replace')
!	open(unit=98, file='cl8000.dat', status='old')
!	do i = 1, 8000
!		read(98,*)j, c_cmb
!		write(99,*)i*(i+1)*c_conf*1.0d12*(2.726**2), c_cmb*1.0d12*(2.726**2)
!		!write(99,*)i*(i+1)*c_conf, c_cmb
!	enddo
!	close(99)
!	close(98)
!	stop
	

end subroutine calc_c_conf

end module confusion_noise

!=======================================================================
