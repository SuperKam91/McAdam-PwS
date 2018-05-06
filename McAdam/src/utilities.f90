module utilities
	use params
      use constants
      use telescope1
	
      
contains

      subroutine interp2d(map,nx,ny,x,y,value)
      
! Subroutine to interpolate a value from a map using bilinear 
! interpolation

      implicit none

      integer nx,ny
	double precision map(nx,ny)
      double precision x, y, value

      double precision t, u, t1, u1
      integer i, j

      i=int(x)
      j=int(y)
      t=x-i
      u=y-j
      t1=1-t
      u1=1-u
      
      value=t1*u1*map(i,j)+t*u1*map(i+1,j)+t*u*map(i+1,j+1)+t1*u*map(i,j+1)

      end subroutine interp2d

!======================================================================

subroutine interp1d(f,x,n,x0,f0)
      
! Subroutine to interpolate a value f0 from a function f using linear 
! interpolation

      	implicit none

	integer n
	double precision f(n),x(n)
      	double precision x0,f0
	integer i

	! First find two x values between which x0 lies:

	if( x(1) >= x0 ) then
	  	f0 = f(1)
		return
	elseif( x(n) <= x0 ) then
	  	f0 = f(n)
	  	return
 	else
	  	call locate(x,n,x0,i)
	endif
	
	! Now do the interpolation:

	f0=f(i)+(f(i+1)-f(i))*(x0-x(i))/(x(i+1)-x(i))
 
end subroutine interp1d

!-----------------------------------------------------------------------

subroutine interp1d_even(f,x,n,x0,f0)

! Subroutine to interpolate a value f0 from a function f using linear
! interpolation, on an evenly spaced grid

      	implicit none

	integer n
	double precision f(n),x(n)
      	double precision x0,f0
        double precision xmin,xmax,dx,wx
	integer i1

        xmin=x(1)
        xmax=x(n)
        dx=x(2)-x(1)

        if(x0.ge.xmin.and.x0.le.xmax) then
               i1= int((x0-xmin)/dx)+1
               wx=(x0-x(i1))/dx
               f0=(1-wx)*f(i1)+wx*f(i1+1)
        else
               write(*,*) 'x = ', x0, '  is out of table range' 
               stop
        endif
        return

end subroutine interp1d_even

!-----------------------------------------------------------------------

      subroutine locate(xx,n,x,j)
      
	implicit none
	
	integer n,j
	double precision xx(n),x
      integer jl,ju,jm
	
	jl=0
      ju=n+1
	
11    if ((ju-jl).gt.1) then
        jm=(ju+jl)/2
        if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
          jl=jm
        else
          ju=jm
        endif
        go to 11
      endif
      
	j=jl
      
	return
      end subroutine locate

!=======================================================================


      function polar(r,theta,Amp)

        double precision var_i, var_j, c1
        parameter(var_i=20., var_j=20., c1=1.)
        double precision r, theta, Amp
        double precision polar

        polar=c1*Amp*r*cos(theta)*exp( -(r**2*(cos(theta))**2 &
        /(2*var_i) )-(r**2*(sin(theta))**2/(2*var_j)) )

        return
        end function polar

!=======================================================================

	function phlog10(x)
	
	implicit none
		
	double precision x,phlog10

	if (x.lt.0.0) then
	  write(*,*) 'phlog10: error, tried to take log of ',x
	  stop
	elseif (x.lt.1e-45_8) then
	  phlog10=-45.0
	else
	  phlog10=log10(x)
	endif
	   
	return
	end function phlog10

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

subroutine qtrap(func,lim1,lim2,eps,S)
      
	implicit none
	
	double precision S,lim1,lim2
	double precision a1,b1,eps,olds,sign,t1,t2
	integer j,jmax
	parameter (jmax=25)
      
      	INTERFACE
		function func(zz)
            		double precision zz,func
            	end function func
      	end INTERFACE

!-----------------------------------

	a1 = lim1
	b1 = lim2
	
	! First catch some stupidities:

	if( a1 == b1 ) then
	  	S = 0.d0
		goto 30
	elseif( a1 > b1 ) then        
	  	olds = b1
	  	b1 = a1
	  	a1 = olds
	  	sign = -1.d0
	else 
	  	sign = 1.d0
	endif  
	 
	olds = -1.d30
	s = 0.d0
      
      	do j = 1,jmax
        	call trapzd(func,a1,b1,S,j)
		if ( j == 1 ) then
			t1 = func(a1)
			t2 = func(b1)
		endif
		 
	        if( abs( s - olds ) < eps * abs(olds) .or. abs( s - olds ) == 0.d0 ) goto 20
	        if( j == jmax ) goto 10
	  
		olds = s
      	enddo
      
 10   	write(*,*) 'QTRAP error: too many steps...'
      	write(*,*) '   S = ',S
      	write(*,*) '   oldS = ',oldS
      	write(*,*) '   % difference = ',100.0*(S-oldS)/oldS
      	write(*,*) '   limits were ',a1,b1
	write(*,*) '   function values at limits were ',t1,t2
      	t1 = func(0.5*(a1+b1))
      	write(*,*) '   function value at midpoint was ',t1
	goto 30
      
 20   	S = S * sign
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
11      continue
        s=0.5*(s+(b1-a1)*sum/tnm)
      endif
      
	return
      end subroutine trapzd
	
!=======================================================================
! Newer Numerical Recipes routine using array-based functions
      FUNCTION qtrap2(func,a,b)
        IMPLICIT NONE
        double precision, INTENT(IN) :: a,b
        double precision :: qtrap2
        INTERFACE
           FUNCTION func(x)
             double precision, DIMENSION(:), INTENT(IN) :: x
             double precision, DIMENSION(size(x)) :: func
           END FUNCTION func
        END INTERFACE
        INTEGER, PARAMETER :: JMAX=25
        double precision, PARAMETER :: EPS=1.0d-4
        double precision :: olds
        INTEGER :: j
        olds=0.0
        do j=1,JMAX
           call trapzd2(func,a,b,qtrap2,j)
           if (j > 5) then
              if (abs(qtrap2-olds) < EPS*abs(olds) .or. &
                   (qtrap2 == 0.0 .and. olds == 0.0)) RETURN
           end if
           olds=qtrap2
        end do

      	write(*,*) 'QTRAP2 error: too many steps...'
      	write(*,*) '   S = ', qtrap2
      	write(*,*) '   oldS = ', olds
      	write(*,*) '   % difference = ',100.0*(qtrap2-olds)/olds
      	write(*,*) '   limits were ',a,b
	write(*,*) '   function values at limits were ',func( (/ a,b /) )
      	write(*,*) '   function value at midpoint was ',func( (/0.5*(a+b)/) )
        stop

      END FUNCTION qtrap2
!=======================================================================

      SUBROUTINE trapzd2(func,a,b,s,n)
        IMPLICIT NONE
        double precision, INTENT(IN) :: a,b
        double precision, INTENT(INOUT) :: s
        INTEGER, INTENT(IN) :: n
        INTERFACE
           FUNCTION func(x)
             implicit none
             double precision, DIMENSION(:), INTENT(IN) :: x
             double precision, DIMENSION(size(x)) :: func
           END FUNCTION func
        END INTERFACE
        double precision :: del,fsum
        INTEGER :: it

       if (n == 1) then
           s=0.5d0*(b-a)*sum(func((/a,b/)))
        else
           it=2**(n-2)
           del=(b-a)/it
           fsum=sum(func(arth(a+0.5d0*del,del,it)))
           s=0.5d0*(s+del*fsum)
        end if
      END SUBROUTINE trapzd2

!=======================================================================

      FUNCTION arth(first,increment,n)
        double precision, INTENT(IN) :: first,increment
        INTEGER, INTENT(IN) :: n
        double precision, DIMENSION(n) :: arth
        INTEGER :: k
        if (n > 0) arth(1)=first
        do k=2,n
           arth(k)=arth(k-1)+increment
        end do
      END FUNCTION arth

!=======================================================================
      FUNCTION ifirstloc(mask)
        implicit none
        LOGICAL, DIMENSION(:), INTENT(IN) :: mask
        INTEGER :: ifirstloc
        INTEGER, DIMENSION(1) :: loc
        loc=maxloc(merge(1,0,mask))
        ifirstloc=loc(1)
        if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
        return
      END FUNCTION ifirstloc
!=======================================================================

      function integrate(y,x,n,xmin,xmax)
!
!  Numerical integration of function y(x): y and x are arrays of
!  length n, and y is integrated between xmin and ymax using the
!  trapezium rule.
!
      implicit none
!       
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
 10   j=i
      ymin=y(j)+(y(j-1)-y(j))*(x(j)-xmin)/(x(j)-x(j-1))
      sum=(ymin+y(j))*(x(j)-xmin)/2.0
!
      do i=j,n
      if (xmax.le.x(i)) goto 20
      enddo
!
 20   k=i-1
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
       subroutine CALC_SEPN(ra1,dec1,ra2,dec2,sepn)
!
!  VSA Project-Reduce program. Finds the separation between two points on
!  the sky-answer returned in radians
!
!  History: 22/8/00-original version [KG]
!           5/12/00-It used to return NaN, when it should 0 (occasionally)
!                      Added line 7 [AS]

      double precision ra1,ra2,dec1,dec2,sepn,cos_sepn
            
      cos_sepn=sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(ra1-ra2)
      if (cos_sepn.gt.1.00000000) cos_sepn=1.0 
      sepn=acos(cos_sepn)

	return
      end subroutine CALC_SEPN

!-----------------------------------------------------------------------

      subroutine CALC_RA_DEC_OFF (ra,dec,ra2,dec2,ra_off,dec_off)
!
! VSA PROJECT
!
! Calculates ra_off,dec_off for an observation centered at ra,dec with
! the offset source at ra2,dec2
!
!
! History: 23/3/01-First version [as, after being crossed with fringe rotation]
!

      double precision ra,dec,ra2,dec2,ra_off,dec_off
      double precision sdo,cdo,srao

! sin of dec out
      sdo=-1.0*sin(dec)*cos(dec2)*(cos(ra)*cos(ra2)+sin(ra)*sin(ra2))+sin(dec2)*cos(dec)
! cos of dec out
      cdo=sqrt(1.0-sdo**2)
! sin of ra out (cf. 'srao' in slavic languages)
      srao=cos(dec2)*(-1.0*sin(ra)*cos(ra2)+cos(ra)*sin(ra2))/cdo

      dec_off=asin(sdo)
      ra_off=asin(srao)
            
	return
	end subroutine CALC_RA_DEC_OFF

!=======================================================================

	function y2Jy(cell,freq)
	
	implicit none
	
	double precision cell,freq
	double precision omega,ICMB,x,g,y2Jy
	
! Solid angle subtended by one pixel in steradians:
 
      omega=1.d0*(cell*sec2rad)*(cell*sec2rad)

! Blackbody formula for CMB radiation 
! (specific intensity in Jy sr^-1):

	!x=1.d0*CRVAL4*Planckfactor
	!ICMB=2.d0*hplanck*CRVAL4*CRVAL4*CRVAL4/(clight*clight)
      x=1.d0*freq*Planckfactor
	ICMB=2.d0*hplanck*freq*freq*freq/(clight*clight)
	ICMB=ICMB/(dexp(x)-1.d0)
	ICMB=ICMB*1.0d26

! Frequency dependent SZ factor (-2 in RJ region):    

	g=x*dexp(x)*(x*(dexp(x)+1.d0)/(dexp(x)-1.d0)-4.d0)
	g=g/(dexp(x)-1.d0)

	y2Jy=g*ICMB*omega 

      return
	end function y2Jy
	
!======================================================================

	double precision function y2Jy_g(cell,freq,g)
	
	implicit none
	
	double precision cell,freq
	double precision omega,ICMB,x,g
	
! Solid angle subtended by one pixel in steradians:
 
      omega=1.d0*(cell*sec2rad)*(cell*sec2rad)

! Blackbody formula for CMB radiation 
! (specific intensity in Jy sr^-1):

	!x=1.d0*CRVAL4*Planckfactor
	!ICMB=2.d0*hplanck*CRVAL4*CRVAL4*CRVAL4/(clight*clight)
      x=1.d0*freq*Planckfactor
	ICMB=2.d0*hplanck*freq*freq*freq/(clight*clight)
	ICMB=ICMB/(dexp(x)-1.d0)
	ICMB=ICMB*1.0d26

! Frequency dependent SZ factor (-2 in RJ region):    

	g=x*dexp(x)*(x*(dexp(x)+1.d0)/(dexp(x)-1.d0)-4.d0)
	g=g/(dexp(x)-1.d0)

	y2Jy_g=g*ICMB*omega 

      return
	end function y2Jy_g
	
!======================================================================
           FUNCTION dT2di(cell,freq)
             IMPLICIT NONE
	     double precision       ::dT2di , cell , freq
	     double precision       ::omega,x,dICMB2dT
 ! Solid angle subtended by one pixel in steradians:
 
      omega=1.d0*(cell*sec2rad)*(cell*sec2rad)
      x=1.d0*freq*Planckfactor
      dICMB2dT=  &
((2.d0*kboltzmann*kboltzmann*kboltzmann*TCMB*TCMB)/(clight*clight*hplanck*hplanck) )* &
      (x*x*x*x*dexp(x))/( (dexp(x) -1.d0)*(dexp(x) -1.d0))
      

     dICMB2dT=dICMB2dT*1.0d26
      
  dT2di= dICMB2dT * omega
  
         END FUNCTION dT2di	
!======================================================================
  
      subroutine cheqboard(array,nx,ny)
      
      implicit none
  
      integer nx,ny
      double complex array(nx,ny)
  
      integer i,j
  
      do i=1,nx-1,2
        do j=1,ny-1,2
          array(i+1,j)=-array(i+1,j)
          array(i,j+1)=-array(i,j+1)
        end do
      end do
                                                                                
      end subroutine cheqboard
  
!======================================================================

      function nintegrate(y,x,n,xmin,xmax)
!
!  Numerical integration of function y(x): y and x are arrays of
!  length n, and y is integrated between xmin and ymax using the
!  trapezium rule.
!
      implicit none
!	
      integer n,i,j,k
      double precision y(n),x(n)
      double precision sum,nintegrate
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
 10   j=i
      ymin=y(j)+(y(j-1)-y(j))*(x(j)-xmin)/(x(j)-x(j-1))
      sum=(ymin+y(j))*(x(j)-xmin)/2.0
!
      do i=j,n
        if (xmax.le.x(i)) goto 20
      enddo
!
 20   k=i-1
      ymax=y(k)-(y(k)-y(k+1))*(xmax-x(k))/(x(k+1)-x(k))
      sum=sum+(ymax+y(k))*(xmax-x(k))/2.0	
!	
      do i=j,k-1
        sum=sum+(y(i+1)+y(i))*(x(i+1)-x(i))/2.0
      enddo
!
      nintegrate=sum
!
      return
      end function nintegrate
!
!======================================================================

	function atanh(x)
	
	double precision x
	double precision atanh
	
	if (x*x.ge.1.d0) then
	  write(*,*) 'This should not happen!, x=',x
	  atanh=0.0
	else
	  atanh=0.5d0*dlog((1.d0+x)/(1.d0-x))
	endif
	
	return
	end function atanh
	
!=======================================================================

      FUNCTION FINDR(FUNC,X1,X2,XACC)
	
	IMPLICIT NONE
	
	double precision FINDR,X1,X2,XACC
	double precision FMID,F,DX,XMID
	INTEGER J,JMAX	
	PARAMETER (JMAX=40)
      
      INTERFACE
            function func(zz)
            	double precision zz,func
            end function func
      end INTERFACE
      
	FMID=FUNC(X2)
      F=FUNC(X1)
      IF ((F*FMID).GT.0.) THEN
! 	  write(*,*) 'f(1)=',F
! 	  write(*,*) 'f(2)=',FMID
! 	  write(*,*) 'Root must be bracketed for bisection.'
	  FINDR=0.0
	  RETURN
	ENDIF  
	
      IF(F.LT.0.)THEN
        FINDR=X1
        DX=X2-X1
      ELSE
        FINDR=X2
        DX=X1-X2
      ENDIF
      DO 11 J=1,JMAX
        DX=DX*.5
        XMID=FINDR+DX
        FMID=FUNC(XMID)
        IF(FMID.LE.0.)FINDR=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.) RETURN
11    CONTINUE
      write(*,*) 'too many bisections'
      END FUNCTION FINDR

!=======================================================================
	
	function rsfunc(rs)
	
	implicit none

	double precision rsfunc,rs
	double precision c1,xE,g,term1,term2,atanh	
	
	xE=rE/rs
	!=r200/rs
	
	if (xE.lt.0.01d0) then
	  g=-0.25*xE*xE*(1.0+dlog(0.25*xE*xE))
	elseif (xE.lt.1.d0) then
	  g=dlog(xE/2.d0)
	  g=g+(2.d0/dsqrt(1.d0-xE*xE))*atanh(dsqrt((1.d0-xE)/(1.d0+xE)))
	elseif (xE.gt.1000.0d0) then
	  g=dlog(xE/2.d0)+1.570796327/xE
	elseif (xE.gt.1.d0) then
	  g=dlog(xE/2.d0)
	  g=g+(2.d0/dsqrt(xE*xE-1.d0))*atan(1.0*dsqrt((xE-1.d0)/(xE+1.d0)))
	else
	  g=dlog(xE/2.d0)
	endif
	term1=xE*xE*scE/(4.0*rs*g)
	
	
	term2=M200/(4.0*pi*rs*rs*rs*(dlog(1.d0+c1)-c1/(c1+1.d0)))
	
	rsfunc=term1-term2
	
	return
	end function rsfunc

!=======================================================================
      !Real Incomplete Beta Function
      !Numerical Recipes
      
      double precision FUNCTION betai(a,b,x)
      
	IMPLICIT NONE
	double precision a,b,x
	double precision bt
      
      if(x<0. .or. x>1.) write(*,*) 'bad argument x in betai'
      
	if (x==0. .or. x==1.) then
		bt=0.0
	else
		bt=(x**a)*((1.-x)**b)
	end if
	if (x<(a+1.)/(a+b+2.)) then
		betai=bt*betacf(a,b,x)/a
	else
		betai=1.-bt*betacf(b,a,1.-x)/b
	end if
	END FUNCTION betai
	
!=======================================================================
      !Used by betai: Evaluates continued fraction for incomplete beta function by modified
	!Lentz's method
      !Numerical Recipes
      
      double precision FUNCTION betacf(a,b,x)
	
      INTEGER MAXIT
	double precision a,b,x,EPS,FPMIN
	PARAMETER (MAXIT=100,EPS=3.d-7,FPMIN=1.d-30)

	INTEGER m,m2
	double precision aa,c,d,del,h,qab,qam,qap
	
      qab=a+b 	!These q's will be used in factors that occur in the coefficients
	qap=a+1.
	qam=a-1.
	c=1. 		!First step of Lentz\u2019s method.
	d=1.-qab*x/qap
	if(abs(d)<FPMIN) d=FPMIN
	d=1./d
	h=d
	do m=1,MAXIT
		m2=2*m
		aa=m*(b-m)*x/((qam+m2)*(a+m2))
		d=1.+aa*d !One step (the even one) of the recurrence.
		if(abs(d)<FPMIN) d=FPMIN
		c=1.+aa/c
		if(abs(c)<FPMIN) c=FPMIN
		d=1./d
		h=h*d*c
		aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
		d=1.+aa*d !Next step of the recurrence (the odd one).
		if(abs(d)<FPMIN) d=FPMIN
		c=1.+aa/c
		if(abs(c)<FPMIN)c=FPMIN
		d=1./d
		del=d*c
		h=h*del
		if(abs(del-1.)<EPS) then
            	betacf=h
                  return
		end if
	enddo
	write(*,*) 'a or b too big, or MAXIT too small in betacf'
 	betacf=h
	return
	END FUNCTION betacf
	
!=======================================================================
	!root funding using the Brent method
      !Numerical Recipes

      FUNCTION zbrent(funct,x1,x2,tol)
      INTEGER ITMAX
      double precision zbrent,tol,x1,x2,EPS
      PARAMETER (ITMAX=100,EPS=3.d-8)
      INTEGER iter
      double precision a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      
      INTERFACE
            function funct(zz)
            	double precision zz,funct
            end function funct
      end INTERFACE
      
      a=x1
      b=x2
      fa=funct(a)
      fb=funct(b)
      if((fa>0. .and. fb>0.) .or. (fa<0. .and. fb<0.))write(*,*) &
      'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb>0. .and. fc>0.) .or. (fb<0. .and. fc<0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc)<abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm)<=tol1 .or. fb==0.)then
          zbrent=b
          return
        endif
        if(abs(e)>=tol1 .and. abs(fa)>abs(fb)) then
          s=fb/fa
          if(a==c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p>0.) q=-q
          p=abs(p)
          if(2.*p<min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d)>tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=funct(b)
11    continue
      write(*,*) 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END function zbrent
	
!=======================================================================
	!binary search
      function binSearch(n,list,value)
      
      implicit none
      
      double precision binSearch !largest i for which which has list(i)>=value
      integer n !num of elements in the list
      double precision list(n) !arrray that needs to be searched
      double precision value !value ot be searched
      integer m1,l1,r1

      !sanity check
      if(list(1)>value) then
      	binSearch=0 !1, surely should be 0?
        return
      elseif(list(n)<value) then
      	binSearch=n+1
        return
      elseif(n==1) then
        if (list(n)==value) then
          binSearch=1
          return
        else
          write(*,*) 'Cannot interpolate value', value, 'in list containing one element = ', list(1)
          stop
        endif
      endif
      
      m1=n/2
      l1=1
      r1=n
      do
      	if(value==list(m1)) then
            	binSearch=m1
                  return
		elseif((l1-r1)**2==1) then
            	binSearch=r1
                  return
      	elseif(value<list(m1)) then
            	r1=m1
		elseif(value>list(m1)) then
            	l1=m1
		endif
            m1=(l1+r1)/2
	enddo
      
	end function binSearch
	
!=======================================================================
!	lookup a 2D table, output 2 values
	subroutine lookUp2D(num1,num2,val1,val2,tab1,tab2,tab3,v)
      
      	implicit none
      	
            integer num1 !total no. of entries in dim 1
            integer num2 !total no. of entries in dim 2
      	double precision val1 !1st value to look for
            double precision val2 !2nd value to look for
            double precision v !output
            double precision tab1(num1) !1st dim
            double precision tab2(num2) !2nd dim
            double precision tab3(num1,num2) !lookup table
            
            integer i,j
            double precision x,y
      	
            !Lookup 1st dim
            i=binSearch(num1,tab1,val1)
      	
            !Lookup 2nd dim
            j=binSearch(num2,tab2,val2)
            
            if(tab1(i)==val1 .and. tab2(j)==val2) then
            	v=tab3(i,j)
		elseif(tab1(i)==val1 .and. .not.(j==0 .or. j==num2)) then
            	v=tab3(i,j-1)+(tab3(i,j)-tab3(i,j-1))*(val2-tab2(j-1))/(tab2(j)-tab2(j-1))
		elseif(tab2(j)==val2 .and. .not.(i==0 .or. i==num1)) then
            	v=tab3(i-1,j)+(tab3(i,j)-tab3(i-1,j))*(val1-tab1(i-1))/(tab1(i)-tab1(i-1))
            !sanity check
            elseif(i==0 .or. j==0 .or. i>num1 .or. j>num2) then
            	write(*,*)"problem in lookup2D, value to look for isn't &
                  	within the table bounds"
		else
            	!bilinear interpolation
			x=(i-1.)+(val1-tab1(i-1))/(tab1(i)-tab1(i-1))
			y=(j-1.)+(val2-tab2(j-1))/(tab2(j)-tab2(j-1))
	            call interp2d(tab3,num1,num2,x,y,v)
		endif
      
      end subroutine lookUp2D
	
!=======================================================================
!	lookup a 2D table, output 2 values
	subroutine lookUp1D(num1,val1,tab1,tab3,v)
      
      	implicit none
      	
            integer num1 !total no. of entries in dim 1
      	double precision val1 !1st value to look for
            double precision v !output
            double precision tab1(num1) !1st dim
            double precision tab3(num1) !lookup table
            
            integer i
      	
            !Lookup 1st dim
            i=binSearch(num1,tab1,val1)
            if(tab1(i)==val1) then
			v=tab3(i)
		!sanity check
            elseif(i==0 .or. i>num1) then
            	write(*,*)"problem in lookup1D, value to look for isn't &
                  	within the table bounds"
		else
            	!bilinear interpolation
			v=tab3(i-1)+(tab3(i)-tab3(i-1))*(val1-tab1(i-1))/(tab1(i)-tab1(i-1))
		endif
      
      end subroutine lookUp1D
	
!=======================================================================
!	check if the given point p is inside the triangle with vertices a, b & c
	logical function inTriangle(a, b, c, p)
      
	implicit none
	
	! input variables
	double precision a(2), b(2), c(2) ! the 3 vertices of the triangle
	double precision p(2) ! the given point
	
	! work variables
	double precision v0(2), v1(2), v2(2), dot00, dot01, dot02, dot11, dot12, invDenom, u ,v
	
	v0 = c - a
	v1 = b - a
	v2 = p - a

	! compute dot products
	dot00 = dot_product(v0, v0)
	dot01 = dot_product(v0, v1)
	dot02 = dot_product(v0, v2)
	dot11 = dot_product(v1, v1)
	dot12 = dot_product(v1, v2)

	! Compute barycentric coordinates
	invDenom = 1d0 / (dot00 * dot11 - dot01 * dot01)
	u = (dot11 * dot02 - dot01 * dot12) * invDenom
	v = (dot00 * dot12 - dot01 * dot02) * invDenom

	! Check if point is in triangle
	if( ( u > 0d0 ) .and. ( v > 0d0 ) .and. ( u + v < 1d0 ) ) then
		inTriangle = .true.
	else
		inTriangle = .false.
	endif

      
	end function inTriangle
	
!=======================================================================
   FUNCTION Gammafun(x)
    
    IMPLICIT NONE
     double precision   Gammafun , x
     INTEGER    k,nn
     double precision  w, y , fun
     double precision   p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13
     PARAMETER(p0 = 0.999999999999999990d+00, p1 = -0.422784335098466784d+00)
     PARAMETER(p2 = -0.233093736421782878d+00,p3 = 0.191091101387638410d+00)
     PARAMETER(p4 = -0.024552490005641278d+00,p5 = -0.017645244547851414d+00)
     PARAMETER(p6 = 0.008023273027855346d+00 ,p7 = -0.000804329819255744d+00)
     PARAMETER(p8 = -0.000360837876648255d+00,p9 = 0.000145596568617526d+00)
     PARAMETER(p10 = -0.000017545539395205d+00,p11 = -0.000002591225267689d+00)
     PARAMETER(p12 = 0.000001337767384067d+00,p13 = -0.000000199542863674d+00)     
     
      nn = nint(x - 2)
      w = x - (nn + 2)
      y = ((((((((((((p13 * w + p12) * w + p11) * w + p10) * &
         w + p9) * w + p8) * w + p7) * w + p6) * w + p5) * &
         w + p4) * w + p3) * w + p2) * w + p1) * w + p0      
      if (nn .gt. 0) then
          w = x - 1
          do k = 2, nn
              w = w * (x - k)
          end do
      else
          w = 1
          do k = 0, -nn - 1
              y = y * (x + k)
          end do
      end if
!      write(*,*)w,y
        fun= w / y
	Gammafun =fun
      end function Gammafun     

!=======================================================================
    subroutine halt_program(message)
        implicit none

#ifdef MPI
        include 'mpif.h'
#endif
        character(LEN=*), intent(in), optional :: message

#ifdef MPI
        integer :: errorcode=1
        integer :: mpierror
#endif

        if (present(message)) then
            write(*,'( 20("=") )')
            write(*,'(A)') trim(adjustl(message))
            write(*,'( 20("=") )')
        end if

#ifdef MPI
        call MPI_ABORT(MPI_COMM_WORLD,errorcode,mpierror)
#else
        stop 1
#endif

        end subroutine halt_program

!======================================

!adapted for GM=6 kj 01/02/17
function alngam ( xvalue, ifault )

!*****************************************************************************80
!
!! ALNGAM computes the logarithm of the gamma function.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Allan Macleod.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, double precision XVALUE, the argument of the Gamma function.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, double precision ALNGAM, the logarithm of the gamma function of X.
!
  implicit none

  double precision alngam
  double precision, parameter :: alr2pi = 0.918938533204673D+00
  integer ifault
  double precision, dimension ( 9 ) :: r1 = (/ &
    -2.66685511495D+00, &
    -24.4387534237D+00, &
    -21.9698958928D+00, &
     11.1667541262D+00, &
     3.13060547623D+00, &
     0.607771387771D+00, &
     11.9400905721D+00, &
     31.4690115749D+00, &
     15.2346874070D+00 /)
  double precision, dimension ( 9 ) :: r2 = (/ &
    -78.3359299449D+00, &
    -142.046296688D+00, &
     137.519416416D+00, &
     78.6994924154D+00, &
     4.16438922228D+00, &
     47.0668766060D+00, &
     313.399215894D+00, &
     263.505074721D+00, &
     43.3400022514D+00 /)
  double precision, dimension ( 9 ) :: r3 = (/ &
    -2.12159572323D+05, &
     2.30661510616D+05, &
     2.74647644705D+04, &
    -4.02621119975D+04, &
    -2.29660729780D+03, &
    -1.16328495004D+05, &
    -1.46025937511D+05, &
    -2.42357409629D+04, &
    -5.70691009324D+02 /)
  double precision, dimension ( 5 ) :: r4 = (/ &
     0.279195317918525D+00, &
     0.4917317610505968D+00, &
     0.0692910599291889D+00, &
     3.350343815022304D+00, &
     6.012459259764103D+00 /)
  double precision x
  double precision x1
  double precision x2
  double precision, parameter :: xlge = 5.10D+05
  double precision, parameter :: xlgst = 1.0D+30
  double precision xvalue
  double precision y

  x = xvalue
  alngam = 0.0D+00
!
!  Check the input.
!
  if ( xlgst <= x ) then
    ifault = 2
    return
  end if

  if ( x <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  ifault = 0
!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
  if ( x < 1.5D+00 ) then

    if ( x < 0.5D+00 ) then

      alngam = - log ( x )
      y = x + 1.0D+00
!
!  Test whether X < machine epsilon.
!
      if ( y == 1.0D+00 ) then
        return
      end if

    else

      alngam = 0.0D+00
      y = x
      x = ( x - 0.5D+00 ) - 0.5D+00

    end if

    alngam = alngam + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
                  y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )

    return

  end if
!
!  Calculation for 1.5 <= X < 4.0.
!
  if ( x < 4.0D+00 ) then

    y = ( x - 1.0D+00 ) - 1.0D+00

    alngam = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
                  x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
!
!  Calculation for 4.0 <= X < 12.0.
!
  else if ( x < 12.0D+00 ) then

    alngam = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
                  x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
!
!  Calculation for 12.0 <= X.
!
  else

    y = log ( x )
    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

    if ( x <= xlge ) then

      x1 = 1.0D+00 / x
      x2 = x1 * x1

      alngam = alngam + x1 * ( ( &
             r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )

    end if

  end if

  return
end function alngam

!=======================================================================

function gammds ( x, p, ifault )

!*****************************************************************************80
!
!! GAMMDS computes the incomplete Gamma integral.
!
!  Discussion:
!
!    The parameters must be positive.  An infinite series is used.
!
!  Auxiliary function:
!
!    ALNGAM = CACM algorithm 291
!
!  Modified:
!
!    22 January 2008
!
!  Author:
!
!    Chi Leung Lau
!    Modifications by John Burkardt
!
!  Reference:
!
!    Chi Leung Lau,
!    Algorithm AS 147:
!    A Simple Series for the Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 29, Number 1, 1980, pages 113-114.
!
!  Parameters:
!
!    Input, double precision X, P, the arguments of the incomplete
!    Gamma integral.  X and P must be greater than 0.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no errors.
!    1, X <= 0 or P <= 0.
!    2, underflow during the computation.
!
!    Output, double precision GAMMDS, the value of the incomplete
!    Gamma integral.
!
  implicit none

  double precision a
  !double precision alngam
  !interface
  !      function alngam(xvalue, ifault)
  !	double precision :: xvalue
  !	integer :: ifault
  !	end function alngam
  !  end interface
  double precision arg
  double precision c
  double precision, parameter :: e = 1.0D-09
  double precision f
  double precision gammds
  integer ifault
  integer ifault2
  double precision p
  double precision, parameter :: uflo = 1.0D-37
  double precision x
!
!  Check the input.
!
  if ( x <= 0.0D+00 ) then
    ifault = 1
    gammds = 0.0D+00
    return
  end if

  if ( p <= 0.0D+00 ) then
    ifault = 1
    gammds = 0.0D+00
    return
  end if
!
!  ALNGAM is the natural logarithm of the gamma function.
!
  ifault2 = 0
  arg = p * log ( x ) - alngam ( p + 1.0D+00, ifault2 ) - x

  if ( arg < log ( uflo ) ) then
    gammds = 0.0D+00
    ifault = 2
    return
  end if

  f = exp ( arg )

  if ( f == 0.0D+00 ) then
    gammds = 0.0D+00
    ifault = 2
    return
  end if

  ifault = 0
!
!  Series begins.
!
  c = 1.0D+00
  gammds = 1.0D+00
  a = p

  do

    a = a + 1.0D+00
    c = c * x / a
    gammds = gammds + c

    if ( c <= e * gammds ) then
      exit
    end if

  end do

  gammds = gammds * f

  return
end function gammds

!=======================================================================

! newton.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Straightforward implementation of the Newton-Raphson method
!     for finding roots of an equation
!
subroutine find_root( f, xinit, tol, maxiter, result, flag )

    interface
        function f(x)
	double precision :: f
	double precision :: x
	end function f
    end interface

    double precision, intent(in)     :: xinit
    double precision, intent(in)     :: tol
    integer, intent(in)  :: maxiter
    double precision, intent(out)    :: result
    integer, intent(out) :: flag

    double precision                 :: eps = 1.0e-4
    double precision                 :: fx1
    double precision                 :: fx2
    double precision                 :: fprime
    double precision                 :: x
    double precision                 :: xnew
    integer              :: i

    result  = 0.0
    flag = 1

    x = xinit
    do i = 1,max(1,maxiter)
        fx1    = f(x)
        fx2    = f(x+eps)
        !write(*,*) i, fx1, fx2, eps
        fprime = (fx2 - fx1) / eps

        xnew   = x - fx1 / fprime

        if ( abs(xnew-x) <= tol ) then
            flag = 0
            result  = xnew
	    !write(*,*) "X from Newton Raphson is ",result
	    !write(*,*) "number of iterations taken: ", i
            exit
        endif

        x = xnew
        !write(*,*) i, x
     enddo

end subroutine find_root

end module utilities
