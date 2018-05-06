module numerics
	
contains
!-----*-----------------------------------------------------------------
	
subroutine qcrude(func,a,b,s,n)
	implicit none
	integer  n
	double precision   a,b,s
	integer  j
	
	INTERFACE
	function func(zz)
	double precision zz,func
	end function func
	end INTERFACE
	
	do j=1,n+1
		call trapzd(func,a,b,s,j)
	enddo
	
end subroutine qcrude
	
!-----*-----------------------------------------------------------------
	
SUBROUTINE qtrap(func,a,b,s,eps,ifail,j)
	INTEGER ifail,JMAX
	double precision a,b,s,EPS
	PARAMETER (JMAX=20)
	INTEGER j
	double precision olds
	
	INTERFACE
		function func(zz)
		double precision zz,func
		end function func
	end INTERFACE
	
	ifail=0
	olds=-1.d30
	
	do j=1,JMAX
		call trapzd(func,a,b,s,j)
		if (dabs(s-olds) <= EPS*dabs(olds)) return
		olds=s
	enddo
	
	ifail=1
	
END SUBROUTINE qtrap
!  (!) Copr. 1986-92 Numerical Recipes Software #!{_W(Y"1i4129.
	
!-----*-----------------------------------------------------------------
	
	SUBROUTINE cqtrap(func,a,b,s_r,s_i,eps,ifail,j)
	INTEGER ifail,JMAX
	double precision a,b,s_r,s_i,EPS
	PARAMETER (JMAX=25)
	INTEGER j
	double precision os_r,os_i
	logical done_r,done_i
	
	INTERFACE
		double complex function func(zz)
		double precision zz
		end function func
	end INTERFACE
	
	ifail=0
	done_r=.false.
	done_i=.false.
	os_r= -1.d30
	os_i= -1.d30
	
	do j=1,JMAX
		call ctrapzd(func,a,b,s_r,s_i,j,.not.done_r,.not.done_i)
	
		if(.not.done_r) then
			if(abs(s_r-os_r) <= EPS*abs(os_r)) then
				done_r=.true.
			else
				os_r=s_r
			endif
	  	endif
	
		if(.not.done_i) then
			if(abs(s_i-os_i) <= EPS*abs(os_i)) then
				done_i=.true.
			else
				os_i=s_i
			endif
	  	endif
	
		if(done_r .and. done_i) return
	enddo
	
	ifail=1
	write(*,*)"CQTRAP did not converge"
END SUBROUTINE cqtrap
!  (!) Copr. 1986-92 Numerical Recipes Software #!{_W(Y"1i4129.
	
!-----*-----------------------------------------------------------------
	
SUBROUTINE qsimp(func,a,b,s,eps,ifail,j)
	INTEGER ifail,JMAX
	double precision a,b,s,EPS
	PARAMETER (JMAX=20)
	INTEGER j
	double precision os,ost,st
	
	INTERFACE
		function func(zz)
		double precision zz,func
		end function func
	end INTERFACE
	
	ifail=0
	
	!sanity check
	if(a==b) then
		s=0.d0
		return
	endif
	
	ost=-1.d30
	os= -1.d30
	
	do j=1,JMAX
		call trapzd(func,a,b,st,j)
		s=(4d0*st-ost)/3d0
		if (dabs(s-os) <= EPS*dabs(os)) return
		os=s
		ost=st
	enddo
	
	ifail=1
	
END SUBROUTINE qsimp
!  (!) Copr. 1986-92 Numerical Recipes Software #!{_W(Y"1i4129.
	
!-----*-----------------------------------------------------------------
	
!numerical integration of a complex function
SUBROUTINE cqsimp(func,a,b,s_r,s_i,eps,ifail,j)
	INTEGER ifail,JMAX
	double precision a,b,s_r,s_i,EPS
	PARAMETER (JMAX=20)
	INTEGER j
	double precision os_r,os_i,ost_r,ost_i,st_r,st_i
	logical done_r,done_i
	
	INTERFACE
		double complex function func(zz)
		double precision zz
		end function func
	end INTERFACE
	
	ifail=0
	done_r=.false.
	done_i=.false.
	ost_r=-1.d30
	ost_i=-1.d30
	os_r= -1.d30
	os_i= -1.d30
	
	do j=1,JMAX
		call ctrapzd(func,a,b,st_r,st_i,j,.not.done_r,.not.done_i)
		
		if(.not.done_r) then
			s_r=(4d0*st_r-ost_r)/3d0
			if(abs(s_r-os_r) <= EPS*abs(os_r)) done_r=.true.
	  	endif
	
		if(.not.done_i) then
			s_i=(4d0*st_i-ost_i)/3d0
			if(abs(s_i-os_i) <= EPS*abs(os_i)) done_i=.true.
	  	endif
	
		if(done_r .and. done_i) return
	
		if(.not.done_r) then
			if(j==JMAX) write(*,*)"CQSIMP:",s_r,os_r,st_r,ost_r
			os_r=s_r
			ost_r=st_r
	  	endif
	
		if(.not.done_i) then
			if(j==JMAX) write(*,*)"CQSIMP:",s_i,os_i,st_i,ost_i
			os_i=s_i
			ost_i=st_i
	  	endif
	enddo
	
	ifail=1
	write(*,*)"CQSIMP did not converge"
	
END SUBROUTINE cqsimp
!  (!) Copr. 1986-92 Numerical Recipes Software #!{_W(Y"1i4129.
	
!-----*-----------------------------------------------------------------
	
! Romberg integration routine - prone to failure owing to poor error
! estimate in POLINT routine. 
SUBROUTINE qromb(func,a,b,ss,eps,ifail,j)
	INTEGER ifail,JMAX,JMAXP,K,KM
	double precision a,b,ss,EPS
	PARAMETER (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
	INTEGER j
	double precision dss,h(JMAXP),s(JMAXP)
	
	INTERFACE
		function func(zz)
		double precision zz,func
		end function func
	end INTERFACE
	
	ifail=0
	h(1)=1d0
	
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		
		if (j.ge.K) then
			call polint(h(j-KM),s(j-KM),K,0d0,ss,dss)
			if (dabs(dss).le.EPS*dabs(ss)) return
		endif
	
		s(j+1)=s(j)
		h(j+1)=0.25d0*h(j)
	enddo
	
	ifail=1
	
END SUBROUTINE qromb
!  (!) Copr. 1986-92 Numerical Recipes Software #!{_W(Y"1i4129.
	
!-----*-----------------------------------------------------------------
	
! this version of qromb uses the simple convergence test in qtrap
! rather than relying in the error estimate in polint.
! Should reproduce QSIMP results of K=3.
	
SUBROUTINE qrombnew(func,a,b,ss,eps,ifail,j)
	INTEGER ifail,JMAX,JMAXP,K,KM
	double precision a,b,ss,EPS
	PARAMETER (JMAX=20, JMAXP=JMAX+1, K=3, KM=K-1)
	INTEGER j
	double precision oldss,dss,h(JMAXP),s(JMAXP)
	
	INTERFACE
	function func(zz)
		double precision zz,func
		end function func
	end INTERFACE
	
	ifail=0
	oldss=-1.d30
	h(1)=1d0
	
	do j=1,JMAX
		call trapzd(func,a,b,s(j),j)
		
		if (j.ge.K) then
			call polint(h(j-KM),s(j-KM),K,0d0,ss,dss)
			if (dabs(ss-oldss) <= EPS*dabs(oldss)) return
			oldss=ss
		endif
	
		s(j+1)=s(j)
		h(j+1)=0.25d0*h(j)
	enddo
	
	ifail=1
	
END SUBROUTINE qrombnew
	
!-----*-----------------------------------------------------------------
	
SUBROUTINE polint(xa,ya,n,x,y,dy)
	INTEGER n,NMAX
	double precision dy,x,y,xa(n),ya(n)
	PARAMETER (NMAX=10)
	INTEGER i,m,ns
	double precision den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
	
	ns=1
	dif=dabs(x-xa(1))
	
	do i=1,n
		dift=dabs(x-xa(i))
		if (dift.lt.dif) then
			ns=i
			dif=dift
		endif
	
		!(i)=ya(i)
		d(i)=ya(i)
	enddo
	
	y=ya(ns)
	ns=ns-1
	
	do m=1,n-1
		do i=1,n-m
			ho=xa(i)-x
			hp=xa(i+m)-x
			w=c(i+1)-d(i)
			den=ho-hp
			if(den.eq.0d0) write(*,*) 'failure in polint'
			den=w/den
			d(i)=hp*den
			c(i)=ho*den
		enddo
	
		if (2*ns.lt.n-m) then
			dy=c(ns+1)
		else
			dy=d(ns)
			ns=ns-1
		endif
	
		y=y+dy
	enddo
	
	return
	
END SUBROUTINE polint
!  (!) Copr. 1986-92 Numerical Recipes Software #!{_W(Y"1i4129.
	
!-----*-----------------------------------------------------------------
	
SUBROUTINE trapzd(func,a,b,s,n)
	INTEGER n
	double precision a,b,s
	INTEGER it,j
	double precision del,sum,tnm,x
	
	INTERFACE
		function func(zz)
		double precision zz,func
		end function func
	end INTERFACE
	
	if (n==1) then
		s=0.5d0*(b-a)*(func(a)+func(b))
	else
		it=2**(n-2)
		tnm=it
		del=(b-a)/tnm
		x=a+0.5d0*del
		sum=0d0
		
		do j=1,it
			sum=sum+func(x)
			x=x+del
		enddo
	
		s=0.5d0*(s+(b-a)*sum/tnm)
	endif
	
END SUBROUTINE trapzd
!  (!) Copr. 1986-92 Numerical Recipes Software #!{_W(Y"1i4129.
	
!-----*-----------------------------------------------------------------
	
!trapezoidal integration of a complex function
SUBROUTINE ctrapzd(func,a,b,s_r,s_i,n,cal_r,cal_i)
	INTEGER n
	logical cal_r,cal_i !calculate the real & imag parts respectively?
	double precision a,b,s_r,s_i
	INTEGER it,j
	double precision del,sum_r,sum_i,tnm,x
	double complex cmplx_a,cmplx_b
	
	INTERFACE
		double complex function func(zz)
		double precision zz
		end function func
	end INTERFACE
	
	if (n==1) then
		cmplx_a=func(a)
		cmplx_b=func(b)
		if(cal_r) s_r=0.5d0*(b-a)*(dble(cmplx_a)+dble(cmplx_a))
		if(cal_i) s_i=0.5d0*(b-a)*(dimag(cmplx_a)+dimag(cmplx_a))
	else
		it=2**(n-2)
		tnm=it
		del=(b-a)/tnm
		x=a+0.5d0*del
		sum_r=0d0
		sum_i=0d0
		do j=1,it
			cmplx_a=func(x)
			if(cal_r) sum_r=sum_r+dble(cmplx_a)
			if(cal_i) sum_i=sum_i+dimag(cmplx_a)
			x=x+del
	  	enddo
		if(cal_r) s_r=0.5d0*(s_r+(b-a)*sum_r/tnm)
		if(cal_i) s_i=0.5d0*(s_i+(b-a)*sum_i/tnm)
	endif
	
END SUBROUTINE ctrapzd
	
!-----*-----------------------------------------------------------------
	
!compute the modulus of a complex number z = x+ iy 
! minimising round-off errors - compare NR.  
! PARAMS 
!     z [I]   complex number z
! RETURNS
!     modulus of z
!
double precision function modulus(z)
	implicit none
	
	double complex z
	double precision x, y, temp
	
	x = abs(dble(z))
	y = abs(dimag(z))
	
	if(x == 0.d0) then
		modulus = y
	elseif(y == 0.d0) then
		modulus = x
	elseif(x > y) then
		temp = y / x
		modulus = x * sqrt(1.d0 + temp*temp)
	else
		temp = x / y
		modulus = y * sqrt(1.d0 + temp*temp)
	endif
	
end function modulus
	
!-----*-----------------------------------------------------------------
	
end module numerics
	
