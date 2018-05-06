module priors1

contains

!=======================================================================
! Prior distribution functions: r is a uniform deviate from the unit
! Cube, Prior returns the corresponding deviate drawn from the desired
! distribution (labelled k).
!
! PJM 9/2002
!
      function Prior(k,r,x1,x2,x3)

      implicit none

      integer k
      double precision r,x1,x2,x3,Prior

      if (k==0) then
       Prior=DeltaFunctionPrior(r,x1,x2)
      elseif (k==1) then
       Prior=UniformPrior(r,x1,x2)
      elseif (k==2) then
       Prior=LogPrior(r,x1,x2)
      elseif (k==3) then
       Prior=GaussianPrior(r,x1,x2)
      elseif (k==4) then
       Prior=LogNormalPrior(r,x1,x2)
      elseif (k==5) then
       Prior=SinPrior(r,x1,x2)
      elseif (k==6) then
       Prior=CauchyPrior(r,x1,x2)
      elseif (k==7) then
       write(*,*)"can't use the spectral index prior for any other parameter"
       stop
      elseif (k==8) then
       write(*,*)"can't use the mass function prior for any parameter other than M200 & z"
       stop
      elseif (k==9) then
       write(*,*)"can't use the triangle prior for any parameter other than x0 and y0"
       stop
      elseif (k==10) then
       Prior=ExpPrior(r,x1,x2,x3)
      elseif (k==11) then
       Prior=PowerPrior(r,x1,x2,x3)
      elseif (k==12) then
       Prior=TruncatedGaussianPrior(r,x1,x2,x3)
      elseif (k==13) then
       write(*,*)"can't use the ratio prior for any parameter other than point source flux"
       stop
      elseif (k==14) then
       write(*,*)"can't use the 2D elliptical Gaussian prior for any parameter other than Ytot and theta_s"
       stop
      endif

      return
      end function Prior

!=======================================================================
! Uniform[0:1]  ->  Delta[x1]

      function DeltaFunctionPrior(r,x1,x2)

      implicit none

      double precision r,x1,x2,DeltaFunctionPrior

      DeltaFunctionPrior=x1

      return
      end function DeltaFunctionPrior

!=======================================================================
! Uniform[0:1]  ->  Uniform[x1:x2]

      function UniformPrior(r,x1,x2)

      implicit none

      double precision r,x1,x2,UniformPrior

      UniformPrior=x1+r*(x2-x1)

      return
      end function UniformPrior

!=======================================================================
! Uniform[0:1]  ->  LogUniform[x1:x2]

      function LogPrior(r,x1,x2)

      implicit none

      double precision r,x1,x2,LogPrior
      double precision lx1,lx2

      if (r < 0.0d0) then
       LogPrior=-1.0d32
      else
       lx1=dlog10(x1)
       lx2=dlog10(x2)
       LogPrior=10.d0**(lx1+r*(lx2-lx1))
      endif

      return
      end function LogPrior

!=======================================================================
! Uniform[0:1]  ->  Sin[x1:x2]  (angles in degrees):

      function SinPrior(r,x1,x2)

      implicit none

      double precision r,x1,x2,SinPrior
      double precision cx1,cx2,deg2rad
      parameter(deg2rad=0.017453292)

      cx1=cos(x1*deg2rad)
      cx2=cos(x2*deg2rad)
      SinPrior=1.d0*acos(cx1+r*(cx2-cx1))

      return
      end function SinPrior

!=======================================================================
! Uniform[0:1]  ->  Gaussian[mean=mu,variance=sigma**2]

      function GaussianPrior(r,mu,sigma)

      implicit none

      double precision r,mu,sigma,GaussianPrior
      double precision SqrtTwo
      parameter(SqrtTwo=1.414213562d0)

      if (r < 1.0d-16.or.(1.0d0-r) < 1.0d-16) then
       GaussianPrior=-1.0d32
      else
       GaussianPrior=mu+sigma*SqrtTwo*dierfc(2.d0*(1.d0-r))
      endif

      return
      end function GaussianPrior

!=======================================================================
! x = Uniform[0:1]  ->  Gaussian[mean=mu,variance=sigma**2] with (x > limit if limit < mu) 
! or (x < limit if limit > mu)

      function TruncatedGaussianPrior(r,mu,sigma,limit)

      implicit none

      double precision r,mu,sigma,limit,TruncatedGaussianPrior
      double precision SqrtTwo, a, b, fact
      parameter(SqrtTwo=1.414213562d0)

      a = limit - mu
      fact = 1.d0
      if (a < 0) then
        fact = -1.d0
        a = fact*a
        r = 1 - r
      endif
       
      b = fact*SqrtTwo*sigma*dierfc(2-(1+erf(a/SqrtTwo/sigma))*r)
      TruncatedGaussianPrior = mu + b

      return
      end function TruncatedGaussianPrior

!=======================================================================
! Uniform[0:1]  ->  Cauchy[mean=x0,FWHM=2*gamma]

      function CauchyPrior(r,x0,gamma)

      implicit none

      double precision r,x0,gamma,CauchyPrior
      double precision Pi
      parameter(Pi=3.141592654)

      CauchyPrior=x0+gamma*tan(Pi*(r-0.5))

      return
      end function CauchyPrior

!=======================================================================
! Uniform[0:1]  ->  LogNormal[mode=a,width parameter=sigma]

      function LogNormalPrior(r,a,sigma)

      implicit none

      double precision r,a,sigma,LogNormalPrior
      double precision SqrtTwo,bracket
      parameter(SqrtTwo=1.414213562d0)

      bracket=sigma*sigma+sigma*SqrtTwo*dierfc(2.d0*r)
      LogNormalPrior=a*dexp(bracket)

      return
      end function LogNormalPrior

!=======================================================================
! Uniform[0:1]  ->  lambda * exp(-lambda * x)[min=a, max=b]

      function ExpPrior(r,a,b,lambda)

      implicit none

      double precision ExpPrior,r,a,b,lambda,temp
      
      if( a > b ) then
      	temp = a
      	a = b
      	b = temp
      elseif( a == b ) then
      	ExpPrior = a
      	return
      endif

      if( a < 0d0 ) then
	write(*,*)'One of the limits for exponential prior is less than 0. Aborting'
      	stop
      endif

      ExpPrior = a - log( 1d0 - r + r * exp( -lambda * (b - a) ) ) / lambda

      end function ExpPrior

!=======================================================================
! Uniform[0:1]  ->  x^{-alpha}[min=a, max=b]

      function PowerPrior(r,a,b,alpha)

      implicit none

      double precision PowerPrior,r,a,b,alpha,p,temp
      
      if( a > b ) then
      	temp = a
      	a = b
      	b = temp
      elseif( a == b ) then
      	PowerPrior = a
      	return
      endif

      if( a < 0d0 ) then
	write(*,*)'One of the limits for power prior is less than 0. Aborting'
      	stop
      endif

      p = 1d0 - alpha
      PowerPrior = ( ( 1d0 - r ) * ( a ** p ) + r * ( b ** p ) ) ** ( 1d0 / p )

      end function PowerPrior

!=======================================================================       
! Inverse of complimentary error function in double precision
!
      double precision function dierfc(y)
      double precision y
      double precision infinity
      double precision x,z,w,u,s,t
      double precision qa,qb,qc,qd,q0,q1,q2,q3,q4
      double precision pa,pb,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22
      parameter (infinity=5.0d0)
      parameter (qa=9.16461398268964d-01, &
      qb=2.31729200323405d-01, &
      qc=4.88826640273108d-01, &
      qd=1.24610454613712d-01, &
      q0=4.99999303439796d-01, &
      q1=1.16065025341614d-01, &
      q2=1.50689047360223d-01, &
      q3=2.69999308670029d-01, &
      q4=-7.28846765585675d-02)
      parameter (pa=3.97886080735226000d+00, &
      pb=1.20782237635245222d-01, &
      p0=2.44044510593190935d-01, &
      p1=4.34397492331430115d-01, &
      p2=6.86265948274097816d-01, &
      p3=9.56464974744799006d-01, &
      p4=1.16374581931560831d+00, &
      p5=1.21448730779995237d+00, &
      p6=1.05375024970847138d+00, &
      p7=7.13657635868730364d-01, &
      p8=3.16847638520135944d-01, &
      p9=1.47297938331485121d-02, &
      p10=-1.05872177941595488d-01, &
      p11=-7.43424357241784861d-02)
      parameter (p12=2.20995927012179067d-03, &
      p13=3.46494207789099922d-02, &
      p14=1.42961988697898018d-02, &
      p15=-1.18598117047771104d-02, &
      p16=-1.12749169332504870d-02, &
      p17=3.39721910367775861d-03, &
      p18=6.85649426074558612d-03, &
      p19=-7.71708358954120939d-04, &
      p20=-3.51287146129100025d-03, &
      p21=1.05739299623423047d-04, &
      p22=1.12648096188977922d-03)
      if (y==0.0) then
        x=infinity
        goto 999
      endif  
      z=y
      if (y > 1) z=2-y
      w=qa-log(z)
      u=sqrt(w)
      s=(qc+log(u))/w
      t=1/(u+qb)
      x=u*(1-s*(0.5d0+s*qd))-((((q4*t+q3)*t+q2)*t+q1)*t+q0)*t
      t=pa/(pa+x)
      u=t-0.5d0
      s=(((((((((p22*u+p21)*u+p20)*u+p19)*u+p18)*u+p17)*u+p16)*u+p15)*u+p14)*u+p13)*u+p12
      s=((((((((((((s*u+p11)*u+p10)*u+p9)*u+p8)*u+p7)*u+p6)*u+p5)*u+p4)*u+p3)*u+p2) &
      *u+p1)*u+p0)*t-z*exp(x*x-pb)
      x=x+s*(1+x*s)
      if (y > 1) x=-x
 999  dierfc=x
      return
      end function dierfc
!=======================================================================
! From Numerical Recipes, eg http://www.aip.de/groups/soe/local/numres/bookfpdf/f6-2.pdf
! Returns the complementary error function erfc(x) with fractional error everywhere less than 1.2 × 10−7.
      FUNCTION erf(x) 

      implicit none

      double precision erfcc,x,erf
      double precision t,z
      
      z=abs(x)
      t=1./(1.+0.5*z) 
      erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
      t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
      t*(1.48851587+t*(-.82215223+t*.17087277))))))))) 
      if (x.lt.0.) erfcc=2.-erfcc
      erf = 1 - erfcc
      return
      END function erf
!=======================================================================
end module priors1
