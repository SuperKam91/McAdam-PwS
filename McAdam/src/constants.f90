module constants

	use params
	
	implicit none
	
! Global variables
	integer, allocatable :: zsrcID(:)

! Physical constants:

	double precision TCMB
	parameter(TCMB=2.726)
	double precision hplanck
	parameter(hplanck=6.6260755d-34)
	double precision kboltzmann
	parameter(kboltzmann=1.380658d-23)
	double precision clight
	parameter(clight=299792458.)
	double precision planckfactor
	parameter(planckfactor=hplanck/(kboltzmann*TCMB))
        double precision sigma_T	
	parameter(sigma_T=6.6524d-29)
        double precision mec2	
	parameter(mec2=511.d0)
        double precision mu_e	
	parameter(mu_e=1.14*1.67262158d-27)
        double precision mu_m	
	parameter(mu_m=0.6*1.67262158d-27)
        double precision m_sun	
	parameter(m_sun=1.98892d+30)			
! Numerical constants:

      double precision Pi
      parameter(Pi = 3.1415926535d0)
      double precision SqrtPi
      parameter (SqrtPi = 1.772453851d0)
      double precision TwoPi
      parameter(TwoPi = 6.283185307d0)
      double precision SqrtTwoPi
      parameter (SqrtTwoPi = 2.506628275d0)

! Conversion factors

      double precision sec2rad,rad2sec
      parameter(sec2rad=4.8481368d-6,rad2sec=1.0/sec2rad)
      double precision min2rad,rad2min
      parameter(min2rad=sec2rad*60.0,rad2min=1.0/min2rad)
      double precision deg2rad,rad2deg
      parameter(deg2rad=sec2rad*3600.0,rad2deg=1.0/deg2rad)
      double precision Mpc2m !mega parsec to metres
      parameter(Mpc2m=3.08568025d22)
      double precision m2Mpc !metres to mega parsec
      parameter(m2Mpc=3.24077649d-23)      
      double precision J2keV     
      parameter(J2keV=6.24251d+15)	
	double precision SigmaSq2G
      parameter (SigmaSq2G = 1.154d8)
	double precision Gmu
      parameter (Gmu=3.12d-14)
      double precision c_mpc !speed of light in Mpc/s
      parameter (c_mpc=9.7156d-15)
      double precision G !Gravitational constant in Mpc^3 M_sun^-1 s^-2
      parameter (G=4.518d-48)
      double precision m_p !proton mass in solar masses
      parameter(m_p=8.40969762d-58)

! Cosmology

	double precision h,Om,OL,Ok,Ob,rhocrit,w0,wa,sigma8
	!parameter(Om=0.20,OL=0.70,rhocrit=2.7752e11,w0=-1.)
      parameter(h=0.7,Om=0.30,OL=0.70,Ok=1.-Om-OL,Ob=0.041)
      parameter(w0=-1.,wa=0.,sigma8=0.8)
      parameter(rhocrit=3.*(h*100./3.0857e19)**2./(8.*Pi*G))
	

end module constants
