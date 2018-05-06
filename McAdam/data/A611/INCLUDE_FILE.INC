! Date: 2014-02-12 18:42:01
!This is a SZ run 
! Object identified as A611
! Weights for P0001_ch8    9157515.0  
! Weights for P0001_ch5    20294234.0  
! Weights for P0001_ch4    21449182.0  
! Weights for P0001_ch7    20460228.0  
! Weights for P0001_ch6    21860902.0  
! Weights for P0001_ch3    6497759.0  
! All sources outside of 0.1 of the primary beam are put as delta priors 
! All sources below 0.400561540269mJy on the SA map (before PB correction) are set as delta priors 
!Cluster temp is 9.7000KeV 
!Cluster redshift is 0.2880 
! Additional comment : 
! 24 sources with delta priors on flux and spectral index 
! 25 total number of sources put into McAdam 
 !======================================================================= 
module params

!=======================================================================
!
!                             McAdam.inc
!
! Include file for McAdam, PJM 7/2003
! Now chose which class of model also, NR 2/2005
!
! Available Classes of Models:
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     1 = Clusters
!
! Available (Cluster) Models:
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^
!   Geometry
!     1 = Circular symmetry, gas and mass have same centroids
!     2 = Elliptical symmetry, gas and mass have same centroids
!
!   Mass
!     1 = NFW profile for mass: parameters c, M_200
!     2 = SIS profile for mass: parameter M_200
!     3 = CPL profile for mass: parameters rcore (kpc) ,alpha,M_proj
!
!   Gas
!     0 = H+M02 survey beta model: parameters thetacore,y0
!     1 = Beta model: parameters rcore,beta, +
!        (parameter 3 = 3d gas mass within r200 or rmass, or gas fraction,
!         or total mass within r200 depending on style
!          Also depending on  style there is a parameter 4 = gas mass fraction
!          within r200)
!     2 = Hydrostatic equilibrium: parameter Mgas_200
!
!     3 = GNFW_Planck model: GNFW pressure profile, parameters are: scaling angular
!         radius(thetas) in arcmin and total volume integrated y parameter(Ysph) in arcmin^2
!     4 = Beta_Atomic/Blobology model: parameters are:thetac(arcsec), beta and central temperature decrement (K)
!
!  Temperature models:
!
!     1 = Isothermal: parameter T / keV
!     2 = Polytropic: parameter T0 / keV & gamma
!
!   Sources
!     x,y,flux all to be fitted
!       -x,y in degrees on the sky
!       -flux in Jy
!
! Available Priors:
! ^^^^^^^^^^^^^^^^
!     0 = Delta function prior (parameter fixed) - (fixedpt,*)
!     1 = Uniform prior - (xmin,xmax)
!     2 = Uniform prior in log_10 - (xmin,xmax)
!     3 = Gaussian prior - (xmean,xstdev)
!     4 = LogNormal prior - (xmean,xwidth)
!     5 = Sinusoidal prior, uniform in cos - (thetamin,thetamax)
!     6 = Cauchy prior - (xmean,xwidth)
!     7 = Prior for spectral indices from Waldram 07
!     8 = joint prior for M200 & z using the mass function (with mass_function = 1, 2 & 3 for Evrard, Jenkins & Tinker mass functions respectively)
!     9 = joint prior on x & y with (x,y) lying uniformly in a triangle defined by the 3 vertices
!     10 = Exponential prior = lambda * exp(-lambda * x) (min, max, lambda)
!     11 = Power law prior = x^{-alpha} (min, max, alpha)
!     12 = Truncated Gaussian prior - (xmean,xstdev,xlimit)
!     or, for spectral indices only:
!     12 = Waldram 07 prior multiplied by Gaussian - (xmean, xstdev)
!     13 = Ratio prior on point source flux
!     14 = 2D elliptical Gaussian prior in log space on theta_s and Y_{tot} (for Planck GNFW model), (centre, width, angle) for both parameters
!
!=======================================================================

! Data:

	integer ModelClass
	integer GL,SZ
	logical simulate,survey,arrGal

! Data filenames:
      integer Nvisfiles
	character*100 GLdatafile,visdatafile,covmatfile

! GL data arrays:
	integer Ngals
	double precision x,y
	double precision e1,e2,e1err,e2err,wt1,wt2
	integer GLeflag
	double precision gamma1,gamma2,kappa,g1,g2
	double precision sigmaobs
	parameter(sigmaobs=0.20)
	logical errcorrect
	parameter(errcorrect=.false.)
	double precision sigcritE
      double precision z_s
      double precision zsmin,zsmax
      !lookup tables
      double precision, dimension(:), allocatable :: SigCritZd,SigCritZs
      double precision, dimension(:,:), allocatable :: lookM,snlook,lookSigCrit,lookD
      double precision, dimension(:,:,:), allocatable :: lookZ
      integer Dn,SCZdn,SCZsn,snn
      !survey data
      integer nxpix,nypix
      integer nxmin,nymin
      double precision pix_sizex,pix_sizey,survey_xl,survey_yl
      double precision e1s,e2s,g1s,g2s,kappas,gamma1s,gamma2s,e1errs,e2errs,wt1s,wt2s,z_ss
      integer GLeflags,ptInPix

! SZ data arrays/maps:
	integer Nvis,large
	double precision u,v,visr,visi,visrms,viswt,SZscale,SZscale_BA
	integer SZeflag
	integer nx,ny
	parameter(nx=256,ny=256)
	double precision ymap(nx,ny),cell,trans,uvcell,uvtrans
	double precision pvisr,pvisi,pvisr1,pvisi1
	double precision LCM
	integer IncludeCMB
	double precision map_x,map_y
      double precision pb_x,pb_y,pb_x_offset,pb_y_offset
      double precision ph_x,ph_y,ph_x_offset,ph_y_offset
	character*100 telescope

! Cluster data and working arrays
	double precision Q
	integer n
	parameter(n=256)
	integer aux_dim               !number of derived parameters
!	parameter (aux_dim=33)       ! for isothermal beta model
!	parameter (aux_dim=6)      ! for GNFW_PLanck model
!	parameter (aux_dim=1)        ! for Beta_Atomic/blobology model
        parameter (aux_dim=25)      ! for DM_GNFW model
	double precision r,rmin,rmax,rlimit,Phi,T,Rhogas,Pgas,Yarray,Rhogas_central,sn_lim,aux
	double precision logr,logPhi,logT,logRhogas,logPgas,logYarray
	double precision uu,loguu , T200 , Mgas200 , fgas200
	parameter(rmin=0.0001,rmax=20.0,rlimit=20.0)
      parameter(sn_lim=0.98)!lensing S/N limit
      double precision r200,M200,r500,rs,ps
      double precision rc,beta,A1,Gm,rhocritz, Rhogas0,logRhogas0
      double precision Ysph500,Ysph200,Ycyl500,Ycyl200
      double precision rE,scE
      double precision rp,Mproj,slope
      double precision a,b,c
      integer    blocks,bsize
      parameter (blocks=1,bsize=2880)
      integer mbuff
      parameter (mbuff=blocks*bsize/4)
      double precision     buffer(mbuff)
      integer    ibuff,nblock
      integer mass_function
      double precision MassLim, MassMin, MassMax
      logical znow
! Cluster data and working arrays for GNFW_Planck model
       REAL *8                  ::  a_GNFW                 ! The constants in GNFW pressure profile: alpha, beta, gamma
       REAL *8                  ::  b_GNFW
       REAL *8                  ::  c_GNFW
       REAL *8                  ::  c500_GNFW                      ! Gas concentration parameter
       REAL *8   , PARAMETER    ::  thetamin_PLanck=0.5       !arcmin
       REAL *8   , PARAMETER    ::  thetamax_Planck=20.0      !arcmin
       REAL *8   , PARAMETER    ::  thetalimit_Planck=20.0     !arcmin


       REAL *8                  ::  thetas_Planck
       REAL *8                  ::  theta500_Planck
       REAL *8                  ::  Y500_Planck
       REAL *8                  ::  Ytot_Planck
       REAL *8                  ::  Ytot_int_Planck
       REAL *8                  ::  y0_Planck , y0_int_Planck
       REAL *8                  ::  ycoeff_Planck
       REAL *8                  ::  map_sum_Planck

! arrays for GNFW_Planck model


       REAL *8                 ::   yintegrand_Planck , logyintegrand_Planck
       REAL *8                 ::   theta_Planck , logtheta_Planck
       REAL *8                 ::   yarray_Planck , logyarray_Planck

! Cluster data and working arrays for Beta_Atomic model
       REAL *8                 ::  thetac_BA , Beta_BA , dT0_BA
       REAL *8                 ::  x0_BA , y0_BA

! cluster data and working arrays for DM_GNFW model

       REAL *8                  ::  rs_DM , rhos_DM, Pei_GNFW, rp_GNFW
       REAL *8                  ::  c200_DM
       REAL *8                  ::  MT200_DM,fg200_DM
       REAL *8                  ::  r200_DM, Mg200_DM, Tg200_DM
       REAL *8                  ::  rhogasi_DM, nei_DM,Tei_DM
       REAL *8                  ::  r500_DM,MT500_DM,Mg500_DM,fg500_DM,Tg500_DM
       REAL *8                  ::  y0_GNFW,ycoeff_GNFW, Ysph500_GNFW,Ysph200_GNFW
       REAL *8                  :: Ycyl_GNFW
       REAL *8                  ::  map_sum_GNFW

! arrays for DM-GNFW model

       REAL *8                  ::  Pe_GNFW,logPe_GNFW,yintegrand_GNFW,logyintegrand_GNFW
       REAL *8                  ::  yarray_GNFW,logyarray_GNFW




!  File control
      integer ifile, filep

! McAdam workspace:
	double complex arr
      integer*8 fftwplan

! FITS header definitions:
	logical SIMPLE,EXTEND
      integer BITPIX,BKSIZE,NAXIS,NAXES(7),PCOUNT,GCOUNT
      character*100 OBJECT,TELESCOP,INSTRUME
      real OBSRA,OBSDEC,CRVAL4,CRVAL5,CRVAL6,PSCAL1,PSCAL2,PSCAL4,PZERO4,PSCAL5,PZERO5
      double precision nu,nu0,flux0,pb_sig

! Ouput sample filenames:
	character*100 root

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Object parameters:

! ??? Are there objects to be fitted?
	integer Atoms
      parameter(Atoms=1)

! ??? How many?
	integer NAtoms,tot_atoms
	parameter(NAtoms=1)

!   Integer flags for cluster distributions that are to be inferred:
	integer Mass,MMass
	integer Gas
	integer Temperature

!   Integer flags for the object models:
	integer GeoModel
	integer MassModel
	integer GasModel
	integer TModel

!   Number of parameters per object (depends on class):
	integer NPars
	integer NPars2

!   Number of parameters per object of each type, and the arrays to
!   store these parameters:
	integer NGeoPars
	double precision GeoPars, GeoPars_old, GeoPars_old1
	integer NMassPars
	double precision MassPars, MassPars_old, MassPars_old1
	integer NGasPars
	double precision GasPars, GasPars_old, GasPars_old1
	integer NTPars
	double precision TPars, TPars_old, TPars_old1
	double precision z(NAtoms),z_old(NAtoms),z_old1(NAtoms),D

! Available Model Classes:
! ^^^^^^^^^^^^^^^^^^^^^^^
!     1 = Clusters

! Cluster Models
! ^^^^^^^^^^^^^^

! Available geometry models:
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
!     1 = Circular symmetry (2 pars, x/arcsec & y/arcsec)
!     2 = Elliptical symmetry (4 pars, x,y plus theta/degrees & axisratio)
! ??? Which geometry model?
	parameter(GeoModel=1)
	parameter(NGeoPars=2*GeoModel)

! Available mass models:
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
!     1 = NFW (2 pars, M200/h-1 Msun & c)
!     2 = SIS (1 par, M200/h-1 Msun)
!     3 = Cored Power Law (3 pars, Mprojected/h-1 Msun, r_c & alpha)
! ??? Which mass model?
	parameter(MassModel=1)
	parameter(NMassPars=2)

! Available gas models:
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
!     0 = Hobson and McLachlan survey Beta model (2 pars, y0 and theta_c)
!     1 = Beta model (3 pars, r_c, beta and Mgas(r200)/h-2 Msun) or (4 pars, r_c, beta , Mtot(r200)/h-1Msun , fgas(r200)/h-1
!     2 = Hydrostatic equilibrium with (NFW) potential (1 par, Mgas/h-2 Msun)
!     3 = GNFW_Planck model (2 pars , thetas/arcmin , Ytot/arcmin2)
!     4 = Beta_Atomic/Blobology model (3 pars: thetac/arcsec , beta , central temperature decrement DeltaT0 /K)
!     5= DM_GNFW model (3 pars: Mtot(r200), fgas(r200)

! ??? Which gas density profile model?
!
!	parameter(GasModel=1)
!	parameter(NGasPars=3)     ! for betastayle(1)=2
!	parameter(NGasPars=4)     ! for betastyle(1)=4

!	parameter(GasModel=3)     ! for GNFW_Planck model
!	parameter(NGasPars=6) 	  ! six gas parameters for GNFW_Planck model: thetas(arcmin) and Ytot(arcmin^2), gamma, alpha, beta, c500 (GNFW shape parameters)

!        parameter(GasModel=4)   !  Beta_Atomic/blobology Model
!        parameter(NGasPars=3)   !  parameters are angular core radius, beta and central temperature decrement

          parameter(GasModel=5)   !  for DM_GNFW model
          parameter(NGasPars=6)   !  parameters are M_tot(r200), fgas(r200), gamma, alpha, beta, c500 (GNFW shape parameters)



! Gas profile styles:
!   different atoms may require different styles:
!    I. Style for Beta model- BetaStyle
!    II. Style for Temperature - TStyle
!   Beta model, GasPars(3) is gas mass within rmass (BetaStyle=1,TStyle=0 ) or r200 in HSE (BetaStyle=2,TStyle=0 )
!   or gas fraction within r200 (BetaStyle=3,TStyle=0)
!   or total mass within r200 and gas mass fraction within r200 in Beta_Virial (BetaStyle=4 ,TStyle=1 )
!   or total mass within r200 and gas mass fraction within r200 in Beta_HSE (BetaStyle=4 ,TStyle=2 )
	integer BetaStyle(NAtoms) , TStyle(NAtoms)
      double precision rmass
      parameter(rmass=1.0)


! Available Temperature models:
! ^^^^^^^^^^^^^^^^^^^^^^^^^^
!     0 = No temperature parameterisation (but need NTpars > 0)
!     1 = Isothermal (1 par, T/keV)
!     2 = Polytropic  (2 pars, T/keV & index gamma)
! ??? Which temperature model?
!	parameter(TModel=1)          ! when sampling from Temperature:(Betastyle=1,2,3)
	parameter(TModel=0)          ! when not sampling from Temperature:(Betastyle=4)
	parameter(NTPars=1)

! Mass profile styles:
!   Weak lensing only use 2(NFW), 1(SIS), 1(CPL)
!   CPLstyle=-1 is needed for CIS profile with M200 calculated.
!   Weak lensing + Einstein radius/arcsec use 4,2,2
!   Obviously different atoms require different styles
	integer NFWstyle(50),SISstyle(NAtoms),CPLstyle(NAtoms)


! Nuisance parameters:
!^^^^^^^^^^^^^^^^^^^^^

! ??? Are there nuisance parameters to be fitted?
	integer Nuisance

!   Number of nuisance parameters - set automatically:
	integer NNuisance

!   Integer flag for presence of point sources in the radio data?
	integer SourceSubtract
! ??? How many sources are there?
	integer NSrc
	parameter(NSrc=25)
	double precision SrcPars

!   Integer flag for varying the redshift of the lensed galaxies?
	integer Varyzs,Vary_zs
	double precision zs,zdmin,zdmax

! ??? How many source planes are there? Relevant only if GL=1 & Varyzs=1
	integer nzsplanes
	parameter(nzsplanes=0)


! Priors and plotting ranges, one for each parameter:
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!  XXX_Prior(i,j,k) -
!       i = atom number
!       j = parameter number, eg beta model, r!->1, beta->2
!       k = 1,2 - prior definition
!
! Cluster arrays:
!^^^^^^^^^^^^^^^^

	integer Geo_PriorType(NAtoms,NGeoPars)
	double precision Geo_Prior(NAtoms,NGeoPars,3)

	integer Mass_PriorType(NAtoms,NMassPars)
	double precision Mass_Prior(NAtoms,NMassPars,3)

	integer Gas_PriorType(NAtoms,NGasPars)
	double precision Gas_Prior(NAtoms,NGasPars,3)

	integer T_PriorType(NAtoms,NTPars)
	double precision T_Prior(NAtoms,NTPars,3)

	integer zs_PriorType(nzsplanes)
	double precision zs_Prior(nzsplanes,3)

	integer z_PriorType(NAtoms)
	double precision z_Prior(NAtoms,3)

! Source arrays:
!^^^^^^^^^^^^^^^

	integer Src_PriorType(NSrc,4)
	double precision Src_Prior(NSrc,4,3)

      !spectral index prior array
      double precision kernel(NSrc,100,3)
      integer nkernel(NSrc)
      double precision prior_min, prior_max

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Example: simulated AMI+VLT cluster

! ??? Make simulated data? or infer from real data?

!     Ensure the simlate= occurs only once in this file for
!     McCompile to operate correctly...
	data simulate/.false./

! Output sample filenames root_?_?_par.dat,
! plus statistics root_?_lhood.dat, root_?_evidence.dat
! ??? Root of output file name?
      parameter(root='data/090807gl')

! Which class of models? 1 = clusters
      parameter(ModelClass=1)

! ??? Infer gas and/or mass and/or temperature distributions?

! ??? Have weak lensing data?
      parameter(GL=0)
! ??? Where is the weak lensing data?
      parameter(GLdatafile='data/GL/sim/090707/2.cat')
! ??? How many galaxies in catalogue?
      parameter(Ngals=0)
! ??? Allow redshift of sources to vary (for prediction)?
	parameter(Varyzs=0)

! ??? Survey mode in GL?
      data survey /.false./
! ??? arrange galaxies intelligently
	parameter(arrGal=.false.)
!no. of pixels in x & y direction
	parameter(nxpix=0,nypix=0)
      parameter(nxmin=0,nymin=0)
!x & y pixel sizes in arcsec
	parameter(pix_sizex=10.546875,pix_sizey=10.546875)
!survey length & breadth
	parameter(survey_xl=1.,survey_yl=1.)

! ??? Have SZ data?
      parameter(SZ=1)
      parameter(Gas=1)
!     parameter(Temperature=1)          ! when sampling from Temperature
      parameter(Temperature=0)          ! when not sampling from Temperature
! ??? Model Mass
      parameter(MMass=0)
! ??? How many SZ data files?
      parameter(Nvisfiles=6)
! ??? Where are the data?
      dimension visdatafile(Nvisfiles)
      data visdatafile/ &
         'DATADIR/P0001_ch3b.fits', &
         'DATADIR/P0001_ch4b.fits', &
         'DATADIR/P0001_ch5b.fits', &
         'DATADIR/P0001_ch6b.fits', &
         'DATADIR/P0001_ch7b.fits', &
        'DATADIR/P0001_ch8b.fits'  &  
      /



! ??? How many uv points in each file?
      dimension Nvis(Nvisfiles)
      data Nvis /873, 962, 1023, 1099, 1147, 979/

! For minimal memory usage, set large >= the largest of the Nvis values
! This can be very wasteful but is needed to allow easy array
! allocation...
! ??? Maximum vis data array size?
      parameter(large=1147)
! 128x128 map is used in FFT - need it big enough to cover x,y priors
! ??? Comptonisation map pixel scale / arcsec?
	parameter(cell=30.0)

! Need to specify which files have an associated covariance matrix,
! setting IncludeCMB=0 neglects CMB and just uses 6th column rms
! reciever noise values in diagonal chisquared computation:
! ??? Include CMB in likelihood via a non-diagonal covariance matrix?
      parameter(IncludeCMB=1)

! ??? Matrix files location? - watch out for byte-swapping...
        data covmatfile/'DATADIR/Covsacov.txt' /


! ??? Centre of coordinate system (RA, Dec in degrees)?
! x,y are measured in arcsec from this, in right handed sense
      parameter(map_x=120.2475,map_y=36.05027778)

! Clusters example....

! set the MassLim (mass in h^{-1}MSun above which we expect to detect all the clusters)
	parameter(MassLim = 2d14)
! set the MassMin (mass in h^{-1}MSun the minimum mass a cluster can have)
	parameter(MassMin = 1d14)
! set the MassMax (mass in h^{-1}MSun the maximum mass a cluster can have)
	parameter(MassMax = 5d15)
! if znow = T, \rho_{crit} is defined to be the present critical density
! if znow = F, \rho_{crit} is defined to be the critical density at the cluster redshift
	parameter(znow=.false.)

! set what mass function to use (relevant only if using joint prior for M200 & z using the mass function)
! mass_function = 1, Evrard approximation to Press-Schechter
! mass_function = 2, Jenkins
! mass_function = 3, Tinker
	parameter(mass_function = 3)

! Example mock data = NFW mass profile + gas in hydrostatic equilibrium
! Model = Isothermal beta model in cored isothermal sphere potential.

! Clump 1:

! ??? Geo
!      data Geo_PriorType(1,1) / 9 /
!      data Geo_Prior(1,1,1),Geo_Prior(1,1,2),Geo_Prior(1,1,3) / -30d0, 30d0, 0d0 /
!      data Geo_PriorType(1,2) / 9 /
!      data Geo_Prior(1,2,1),Geo_Prior(1,2,2),Geo_Prior(1,2,3) / -30d0, -30d0, 30d0 /
!      data Geo_PriorType(1,1) / 3 /
!      data Geo_Prior(1,1,1),Geo_Prior(1,1,2) / 0.d0, 60.d0 /
!      data Geo_PriorType(1,2) / 3 /
!      data Geo_Prior(1,2,1),Geo_Prior(1,2,2) / 0.d0, 60.d0 /
!      data Geo_PriorType(2,1) / 3 /
!      data Geo_Prior(2,1,1),Geo_Prior(2,1,2) / 0.d0, 60.d0 /
!      data Geo_PriorType(2,2) / 3 /
!      data Geo_Prior(2,2,1),Geo_Prior(2,2,2) / 0.d0, 60.d0 /

      data Geo_PriorType(1,1) / 3 /
      data Geo_Prior(1,1,1),Geo_Prior(1,1,2) / 0.d0, 60.d0 /
      data Geo_PriorType(1,2) / 3 /
      data Geo_Prior(1,2,1),Geo_Prior(1,2,2) / 0.d0, 60.d0 /
!      data Geo_PriorType(1,3) / 1 / ! orientation angle phi
!      data Geo_Prior(1,3,1),Geo_Prior(1,3,2) / 0.d0, 180.d0 /

!      data Geo_PriorType(1,4) / 1 / ! ratio of semi-minor to semi-major axes f
!      data Geo_Prior(1,4,1),Geo_Prior(1,4,2) / 0.5d0, 1d0 /



! ??? NFW profile
      data NFWstyle(1) / 2 /
      data Mass_PriorType(1,1) / 1 /
      data Mass_Prior(1,1,1),Mass_Prior(1,1,2) / 0.1, 15. /
      data Mass_PriorType(1,2) / 3 /
      data Mass_Prior(1,2,1),Mass_Prior(1,2,2) / 5.d14,2.5d14 /
!      data NFWstyle(2) / 2 /
!      data Mass_PriorType(2,1) / 1 /
!      data Mass_Prior(2,1,1),Mass_Prior(2,1,2) / 0.1, 15. /
!      data Mass_PriorType(2,2) / 3 /
!      data Mass_Prior(2,2,1),Mass_Prior(2,2,2) / 5.d14,2.5d14 /

! Beta profile
! ( betastyle 2 TStyle 0 ):
! Parameters rcore (kpc), beta, gas mass(r200)(h-2Msun)
!      data BetaStyle(1) / 2 /
!      data TStyle(1) / 0 /     ! sampling from temperature
!      data Gas_PriorType(1,1) / 1 /
!      data Gas_Prior(1,1,1),Gas_Prior(1,1,2) / 10.0d0, 1000.0d0 /
!      data Gas_PriorType(1,2) / 1 /
!      data Gas_Prior(1,2,1),Gas_Prior(1,2,2) / 0.3d0, 1.5d0 /
!      data Gas_PriorType(1,3) / 2 /
!      data Gas_Prior(1,3,1),Gas_Prior(1,3,2) / 1.d12, 5.0d14 /
!----- In case of multi atoms---------------------------------
!      data BetaStyle(2) / 2 /
!       data TStyle(1) / 0 /
!      data Gas_PriorType(2,1) / 1 /
!      data Gas_Prior(2,1,1),Gas_Prior(2,1,2) / 10.0d0, 1000.0d0 /
!      data Gas_PriorType(2,2) / 1 /
!      data Gas_Prior(2,2,1),Gas_Prior(2,2,2) / 0.3d0, 1.5d0 /
!      data Gas_PriorType(2,3) / 2 /
!      data Gas_Prior(2,3,1),Gas_Prior(2,3,2) / 1.d12, 5.0d14 /
!-------------------------------------------------------------
! Beta profile and Virial theorem
! ( betastyle 4 , TStyle 1):
! or
! Beta profile and HSE
! ( betastyle 4 , TStyle 2):
! 4 Parameters rcore (h-1kpc), beta, tot mass(r200)(h-1Msun), fgas(r200) (h-1)

!      	data BetaStyle(1) / 4 /
!        data TStyle(1) / 1 /         ! uses mass - temperature scaling relation to calculate temperature
!       data TStyle(1) / 2 /         ! uses HSE to calculate temperature
!      	data Gas_PriorType(1,1) / 1 /
!      	data Gas_Prior(1,1,1),Gas_Prior(1,1,2) / 10.0d0, 1000.d0 /     !rc(h-1Kpc)
!      	data Gas_PriorType(1,2) / 1 /
!      	data Gas_Prior(1,2,1),Gas_Prior(1,2,2) / 0.3d0, 2.5d0 /        !beta
!      	data Gas_PriorType(1,3) / 2 /
!      	data Gas_Prior(1,3,1),Gas_Prior(1,3,2) / 3.16d13, 5.0d15 /      !M200(h-1Msun)
!      	data Gas_PriorType(1,4) / 3 /
!      	data Gas_Prior(1,4,1),Gas_Prior(1,4,2) / 0.123, 0.0308  /          !fgas200(h-1) updated from 0.0864 to 0.0864/0.7=0.123 also error updated from 0.0216 to 0.0308 on 08/04/2011
!=========================================================================================
! GNFW_Planck model:six parameters: scaling angular radius(thetas /arcmin), total volume
! integrated Y parameter (Ytot/arcmin^2), gamma, alpha, beta (GNFW shape parameters),
! c500 (concentration parameter)
! Joint 2D elliptical Gaussian prior in log-space
!       data Gas_PriorType(1,1) / 14 /
!       data Gas_Prior(1,1,1),Gas_Prior(1,1,2) ,Gas_Prior(1,1,3) / 0.6171, 0.1153, 0.7011 /        !thetas_Planck(arcmin)
!       data Gas_PriorType(1,2) / 14 /
!       data Gas_Prior(1,2,1),Gas_Prior(1,2,2) ,Gas_Prior(1,2,3)/ -2.743, 0.2856, 0.7011 /         !Ytot_Planck (arcmin2)
!       data Gas_PriorType(1,3) / 0 /
!       data Gas_Prior(1,3,1),Gas_Prior(1,3,2) / 0.3081d0, 0.0 /  !c_gnfw
!       data Gas_PriorType(1,4) / 0 /
!       data Gas_Prior(1,4,1),Gas_Prior(1,4,2) / 1.0510d0, 0.0 /  !a_gnfw
!       data Gas_PriorType(1,5) / 0 /
!       data Gas_Prior(1,5,1),Gas_Prior(1,5,2) / 5.4905d0, 0.0 /  !b_gnfw
!       data Gas_PriorType(1,6) / 0 /
!       data Gas_Prior(1,6,1),Gas_Prior(1,6,2) / 1.177d0, 0.0 /  !c500_gnfw
!===========================================================================================
!Beta_Atomic/blobology Model:
!The priors for Beta_Atomic Model are: thetac(arcsec), beta and central temperature decrement (K)
!     	data Gas_PriorType(1,1) / 1/
!      	data Gas_Prior(1,1,1),Gas_Prior(1,1,2) / 20.0d0, 500.d0 /     !thetac('')
!      	data Gas_PriorType(1,2) / 1 /
!     	data Gas_Prior(1,2,1),Gas_Prior(1,2,2) / 0.4d0, 2.5d0 /        !beta
!      	data Gas_PriorType(1,3) / 1 /
!        data Gas_Prior(1,3,1),Gas_Prior(1,3,2) /-3000.d0 , -10.d0 /    !deltaT0(microK)
!========================================================================================
! DM_GNFW Model
! The priors for DM_GNFW Model are:  Mtot(r200)and fgas(r200), gamma, alpha, beta (GNFW shape parameters),
! c500 (concentration parameter)
        data Gas_PriorType(1,1) / 2 /
        data Gas_Prior(1,1,1),Gas_Prior(1,1,2) /1.0d14, 6.0d15 /    !M200(Msun)
        data Gas_PriorType(1,2) / 3 /
         data Gas_Prior(1,2,1),Gas_Prior(1,2,2) /0.13 ,0.02 /    !fgas200
        data Gas_PriorType(1,3) / 0 /
        data Gas_Prior(1,3,1),Gas_Prior(1,3,2) / 0.3081d0, 0.0 /  !c_gnfw
        data Gas_PriorType(1,4) / 0 /
        data Gas_Prior(1,4,1),Gas_Prior(1,4,2) / 1.0510d0, 0.0 /  !a_gnfw
        data Gas_PriorType(1,5) / 0 /
        data Gas_Prior(1,5,1),Gas_Prior(1,5,2) / 5.4905d0, 0.0 /  !b_gnfw
        data Gas_PriorType(1,6) / 0 /
        data Gas_Prior(1,6,1),Gas_Prior(1,6,2) / 1.177d0, 0.0 /  !c500_gnfw
!======================================================================================

! Temperature Model (isothermal)
        data T_PriorType(1,1) / 3 /
        data T_Prior(1,1,1),T_Prior(1,1,2) / 9.7000d0, 1.0d0 /

! redshift
  	data z_PriorType(1) / 0 /
        data z_Prior(1,1),z_Prior(1,2) / 0.2880d0  ,  0.2880d0 /

! central frequency
	data nu0(0) / 15.7d9 /

! Point sources (RA, dec, flux):
	!prior ranges for the spectral index prior
	data prior_min, prior_max / -1.d0, 2.5d0 /

! Point Source

! Point Source 1
        data Src_PriorType(1,1) / 0 /
        data Src_Prior(1,1,1),Src_Prior(1,1,2) / 1.200329d+02, 0.001389 /
        data Src_PriorType(1,2) / 0 /
        data Src_Prior(1,2,1),Src_Prior(1,2,2) / 3.606900d+01, 0.001389 /
        data Src_PriorType(1,3) / 3 /
        data Src_Prior(1,3,1),Src_Prior(1,3,2) / 6.406854d-03 , 2.562741d-03 /
        data nu0(1) / 1.552000d+10 /
        data Src_PriorType(1,4) / 3 /
        data Src_Prior(1,4,1),Src_Prior(1,4,2) / -1.192123d-02 , 1.869074d-01 /

! Point Source 2
        data Src_PriorType(2,1) / 0 /
        data Src_Prior(2,1,1),Src_Prior(2,1,2) / 1.201693d+02, 0.001389 /
        data Src_PriorType(2,2) / 0 /
        data Src_Prior(2,2,1),Src_Prior(2,2,2) / 3.624016d+01, 0.001389 /
        data Src_PriorType(2,3) / 0 /
        data Src_Prior(2,3,1),Src_Prior(2,3,2) / 2.916799d-03 , 1.166720d-03 /
        data nu0(2) / 1.557000d+10 /
        data Src_PriorType(2,4) / 0 /
        data Src_Prior(2,4,1),Src_Prior(2,4,2) / 1.419320d+00 , 5.384796d-01 /

! Point Source 3
        data Src_PriorType(3,1) / 0 /
        data Src_Prior(3,1,1),Src_Prior(3,1,2) / 1.199632d+02, 0.001389 /
        data Src_PriorType(3,2) / 0 /
        data Src_Prior(3,2,1),Src_Prior(3,2,2) / 3.618514d+01, 0.001389 /
        data Src_PriorType(3,3) / 0 /
        data Src_Prior(3,3,1),Src_Prior(3,3,2) / 3.019666d-03 , 1.207866d-03 /
        data nu0(3) / 1.538000d+10 /
        data Src_PriorType(3,4) / 0 /
        data Src_Prior(3,4,1),Src_Prior(3,4,2) / -1.386926d+00 , 1.009091d+00 /

! Point Source 4
        data Src_PriorType(4,1) / 0 /
        data Src_Prior(4,1,1),Src_Prior(4,1,2) / 1.205026d+02, 0.001389 /
        data Src_PriorType(4,2) / 0 /
        data Src_Prior(4,2,1),Src_Prior(4,2,2) / 3.614918d+01, 0.001389 /
        data Src_PriorType(4,3) / 0 /
        data Src_Prior(4,3,1),Src_Prior(4,3,2) / 2.857467d-03 , 1.142987d-03 /
        data nu0(4) / 1.557000d+10 /
        data Src_PriorType(4,4) / 0 /
        data Src_Prior(4,4,1),Src_Prior(4,4,2) / 7.258344d-01 , 5.513006d-01 /

! Point Source 5
        data Src_PriorType(5,1) / 0 /
        data Src_Prior(5,1,1),Src_Prior(5,1,2) / 1.201793d+02, 0.001389 /
        data Src_PriorType(5,2) / 0 /
        data Src_Prior(5,2,1),Src_Prior(5,2,2) / 3.623451d+01, 0.001389 /
        data Src_PriorType(5,3) / 0 /
        data Src_Prior(5,3,1),Src_Prior(5,3,2) / 2.272726d-03 , 9.090905d-04 /
        data nu0(5) / 1.558000d+10 /
        data Src_PriorType(5,4) / 0 /
        data Src_Prior(5,4,1),Src_Prior(5,4,2) / 1.031309d+00 , 6.423221d-01 /

! Point Source 6
        data Src_PriorType(6,1) / 0 /
        data Src_Prior(6,1,1),Src_Prior(6,1,2) / 1.199835d+02, 0.001389 /
        data Src_PriorType(6,2) / 0 /
        data Src_Prior(6,2,1),Src_Prior(6,2,2) / 3.597604d+01, 0.001389 /
        data Src_PriorType(6,3) / 0 /
        data Src_Prior(6,3,1),Src_Prior(6,3,2) / 2.397663d-03 , 9.590654d-04 /
        data nu0(6) / 1.562000d+10 /
        data Src_PriorType(6,4) / 0 /
        data Src_Prior(6,4,1),Src_Prior(6,4,2) / 1.998252d-01 , 6.507359d-01 /

! Point Source 7
        data Src_PriorType(7,1) / 0 /
        data Src_Prior(7,1,1),Src_Prior(7,1,2) / 1.199516d+02, 0.001389 /
        data Src_PriorType(7,2) / 0 /
        data Src_Prior(7,2,1),Src_Prior(7,2,2) / 3.611259d+01, 0.001389 /
        data Src_PriorType(7,3) / 0 /
        data Src_Prior(7,3,1),Src_Prior(7,3,2) / 1.870947d-03 , 7.483788d-04 /
        data nu0(7) / 1.562000d+10 /
        data Src_PriorType(7,4) / 0 /
        data Src_Prior(7,4,1),Src_Prior(7,4,2) / 2.058270d+00 , 9.097151d-01 /

! Point Source 8
        data Src_PriorType(8,1) / 0 /
        data Src_Prior(8,1,1),Src_Prior(8,1,2) / 1.200486d+02, 0.001389 /
        data Src_PriorType(8,2) / 0 /
        data Src_Prior(8,2,1),Src_Prior(8,2,2) / 3.583793d+01, 0.001389 /
        data Src_PriorType(8,3) / 0 /
        data Src_Prior(8,3,1),Src_Prior(8,3,2) / 1.800681d-03 , 7.202723d-04 /
        data nu0(8) / 1.566000d+10 /
        data Src_PriorType(8,4) / 0 /
        data Src_Prior(8,4,1),Src_Prior(8,4,2) / -3.814897d-01 , 8.547186d-01 /

! Point Source 9
        data Src_PriorType(9,1) / 0 /
        data Src_Prior(9,1,1),Src_Prior(9,1,2) / 1.205515d+02, 0.001389 /
        data Src_PriorType(9,2) / 0 /
        data Src_Prior(9,2,1),Src_Prior(9,2,2) / 3.606450d+01, 0.001389 /
        data Src_PriorType(9,3) / 0 /
        data Src_Prior(9,3,1),Src_Prior(9,3,2) / 1.701458d-03 , 6.805831d-04 /
        data nu0(9) / 1.559000d+10 /
        data Src_PriorType(9,4) / 0 /
        data Src_Prior(9,4,1),Src_Prior(9,4,2) / 1.108764d+00 , 9.230058d-01 /

! Point Source 10
        data Src_PriorType(10,1) / 0 /
        data Src_Prior(10,1,1),Src_Prior(10,1,2) / 1.202289d+02, 0.001389 /
        data Src_PriorType(10,2) / 0 /
        data Src_Prior(10,2,1),Src_Prior(10,2,2) / 3.629587d+01, 0.001389 /
        data Src_PriorType(10,3) / 0 /
        data Src_Prior(10,3,1),Src_Prior(10,3,2) / 1.422945d-03 , 5.691779d-04 /
        data nu0(10) / 1.584000d+10 /
        data Src_PriorType(10,4) / 0 /
        data Src_Prior(10,4,1),Src_Prior(10,4,2) / -5.650678d-01 , 1.280430d+00 /

! Point Source 11
        data Src_PriorType(11,1) / 0 /
        data Src_Prior(11,1,1),Src_Prior(11,1,2) / 1.202461d+02, 0.001389 /
        data Src_PriorType(11,2) / 0 /
        data Src_Prior(11,2,1),Src_Prior(11,2,2) / 3.593108d+01, 0.001389 /
        data Src_PriorType(11,3) / 0 /
        data Src_Prior(11,3,1),Src_Prior(11,3,2) / 9.141387d-04 , 3.656555d-04 /
        data nu0(11) / 1.563000d+10 /
        data Src_PriorType(11,4) / 0 /
        data Src_Prior(11,4,1),Src_Prior(11,4,2) / 1.053630d+00 , 7.750558d-01 /

! Point Source 12
        data Src_PriorType(12,1) / 0 /
        data Src_Prior(12,1,1),Src_Prior(12,1,2) / 1.201848d+02, 0.001389 /
        data Src_PriorType(12,2) / 0 /
        data Src_Prior(12,2,1),Src_Prior(12,2,2) / 3.628001d+01, 0.001389 /
        data Src_PriorType(12,3) / 0 /
        data Src_Prior(12,3,1),Src_Prior(12,3,2) / 1.474000d-03 , 5.896000d-04 /
        data nu0(12) / 1.581162d+10 /
        data Src_PriorType(12,4) / 0 /
        data Src_Prior(12,4,1),Src_Prior(12,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 13
        data Src_PriorType(13,1) / 0 /
        data Src_Prior(13,1,1),Src_Prior(13,1,2) / 1.203537d+02, 0.001389 /
        data Src_PriorType(13,2) / 0 /
        data Src_Prior(13,2,1),Src_Prior(13,2,2) / 3.609409d+01, 0.001389 /
        data Src_PriorType(13,3) / 0 /
        data Src_Prior(13,3,1),Src_Prior(13,3,2) / 5.750000d-04 , 2.300000d-04 /
        data nu0(13) / 1.581162d+10 /
        data Src_PriorType(13,4) / 0 /
        data Src_Prior(13,4,1),Src_Prior(13,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 14
        data Src_PriorType(14,1) / 0 /
        data Src_Prior(14,1,1),Src_Prior(14,1,2) / 1.202145d+02, 0.001389 /
        data Src_PriorType(14,2) / 0 /
        data Src_Prior(14,2,1),Src_Prior(14,2,2) / 3.586446d+01, 0.001389 /
        data Src_PriorType(14,3) / 0 /
        data Src_Prior(14,3,1),Src_Prior(14,3,2) / 4.990000d-04 , 1.996000d-04 /
        data nu0(14) / 1.581162d+10 /
        data Src_PriorType(14,4) / 0 /
        data Src_Prior(14,4,1),Src_Prior(14,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 15
        data Src_PriorType(15,1) / 0 /
        data Src_Prior(15,1,1),Src_Prior(15,1,2) / 1.203207d+02, 0.001389 /
        data Src_PriorType(15,2) / 0 /
        data Src_Prior(15,2,1),Src_Prior(15,2,2) / 3.607354d+01, 0.001389 /
        data Src_PriorType(15,3) / 0 /
        data Src_Prior(15,3,1),Src_Prior(15,3,2) / 5.136129d-04 , 2.054452d-04 /
        data nu0(15) / 1.561000d+10 /
        data Src_PriorType(15,4) / 0 /
        data Src_Prior(15,4,1),Src_Prior(15,4,2) / -3.092238d-01 , 1.462213d+00 /

! Point Source 16
        data Src_PriorType(16,1) / 0 /
        data Src_Prior(16,1,1),Src_Prior(16,1,2) / 1.201604d+02, 0.001389 /
        data Src_PriorType(16,2) / 0 /
        data Src_Prior(16,2,1),Src_Prior(16,2,2) / 3.618312d+01, 0.001389 /
        data Src_PriorType(16,3) / 0 /
        data Src_Prior(16,3,1),Src_Prior(16,3,2) / 4.280000d-04 , 1.712000d-04 /
        data nu0(16) / 1.581162d+10 /
        data Src_PriorType(16,4) / 0 /
        data Src_Prior(16,4,1),Src_Prior(16,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 17
        data Src_PriorType(17,1) / 0 /
        data Src_Prior(17,1,1),Src_Prior(17,1,2) / 1.201252d+02, 0.001389 /
        data Src_PriorType(17,2) / 0 /
        data Src_Prior(17,2,1),Src_Prior(17,2,2) / 3.601024d+01, 0.001389 /
        data Src_PriorType(17,3) / 0 /
        data Src_Prior(17,3,1),Src_Prior(17,3,2) / 4.090000d-04 , 1.636000d-04 /
        data nu0(17) / 1.581162d+10 /
        data Src_PriorType(17,4) / 0 /
        data Src_Prior(17,4,1),Src_Prior(17,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 18
        data Src_PriorType(18,1) / 0 /
        data Src_Prior(18,1,1),Src_Prior(18,1,2) / 1.202994d+02, 0.001389 /
        data Src_PriorType(18,2) / 0 /
        data Src_Prior(18,2,1),Src_Prior(18,2,2) / 3.617277d+01, 0.001389 /
        data Src_PriorType(18,3) / 0 /
        data Src_Prior(18,3,1),Src_Prior(18,3,2) / 3.730000d-04 , 1.492000d-04 /
        data nu0(18) / 1.581162d+10 /
        data Src_PriorType(18,4) / 0 /
        data Src_Prior(18,4,1),Src_Prior(18,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 19
        data Src_PriorType(19,1) / 0 /
        data Src_Prior(19,1,1),Src_Prior(19,1,2) / 1.202868d+02, 0.001389 /
        data Src_PriorType(19,2) / 0 /
        data Src_Prior(19,2,1),Src_Prior(19,2,2) / 3.596247d+01, 0.001389 /
        data Src_PriorType(19,3) / 0 /
        data Src_Prior(19,3,1),Src_Prior(19,3,2) / 3.470000d-04 , 1.388000d-04 /
        data nu0(19) / 1.581162d+10 /
        data Src_PriorType(19,4) / 0 /
        data Src_Prior(19,4,1),Src_Prior(19,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 20
        data Src_PriorType(20,1) / 0 /
        data Src_Prior(20,1,1),Src_Prior(20,1,2) / 1.201749d+02, 0.001389 /
        data Src_PriorType(20,2) / 0 /
        data Src_Prior(20,2,1),Src_Prior(20,2,2) / 3.589287d+01, 0.001389 /
        data Src_PriorType(20,3) / 0 /
        data Src_Prior(20,3,1),Src_Prior(20,3,2) / 3.470000d-04 , 1.388000d-04 /
        data nu0(20) / 1.581162d+10 /
        data Src_PriorType(20,4) / 0 /
        data Src_Prior(20,4,1),Src_Prior(20,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 21
        data Src_PriorType(21,1) / 0 /
        data Src_Prior(21,1,1),Src_Prior(21,1,2) / 1.202193d+02, 0.001389 /
        data Src_PriorType(21,2) / 0 /
        data Src_Prior(21,2,1),Src_Prior(21,2,2) / 3.610514d+01, 0.001389 /
        data Src_PriorType(21,3) / 0 /
        data Src_Prior(21,3,1),Src_Prior(21,3,2) / 3.170000d-04 , 1.268000d-04 /
        data nu0(21) / 1.581162d+10 /
        data Src_PriorType(21,4) / 0 /
        data Src_Prior(21,4,1),Src_Prior(21,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 22
        data Src_PriorType(22,1) / 0 /
        data Src_Prior(22,1,1),Src_Prior(22,1,2) / 1.201642d+02, 0.001389 /
        data Src_PriorType(22,2) / 0 /
        data Src_Prior(22,2,1),Src_Prior(22,2,2) / 3.599703d+01, 0.001389 /
        data Src_PriorType(22,3) / 0 /
        data Src_Prior(22,3,1),Src_Prior(22,3,2) / 2.650000d-04 , 1.060000d-04 /
        data nu0(22) / 1.581162d+10 /
        data Src_PriorType(22,4) / 0 /
        data Src_Prior(22,4,1),Src_Prior(22,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 23
        data Src_PriorType(23,1) / 0 /
        data Src_Prior(23,1,1),Src_Prior(23,1,2) / 1.203329d+02, 0.001389 /
        data Src_PriorType(23,2) / 0 /
        data Src_Prior(23,2,1),Src_Prior(23,2,2) / 3.608933d+01, 0.001389 /
        data Src_PriorType(23,3) / 0 /
        data Src_Prior(23,3,1),Src_Prior(23,3,2) / 2.430000d-04 , 9.720000d-05 /
        data nu0(23) / 1.581162d+10 /
        data Src_PriorType(23,4) / 0 /
        data Src_Prior(23,4,1),Src_Prior(23,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 24
        data Src_PriorType(24,1) / 0 /
        data Src_Prior(24,1,1),Src_Prior(24,1,2) / 1.201747d+02, 0.001389 /
        data Src_PriorType(24,2) / 0 /
        data Src_Prior(24,2,1),Src_Prior(24,2,2) / 3.597491d+01, 0.001389 /
        data Src_PriorType(24,3) / 0 /
        data Src_Prior(24,3,1),Src_Prior(24,3,2) / 2.420000d-04 , 9.680000d-05 /
        data nu0(24) / 1.581162d+10 /
        data Src_PriorType(24,4) / 0 /
        data Src_Prior(24,4,1),Src_Prior(24,4,2) / 5.000000d-01 , 0.000000d+00 /

! Point Source 25
        data Src_PriorType(25,1) / 0 /
        data Src_Prior(25,1,1),Src_Prior(25,1,2) / 1.202822d+02, 0.001389 /
        data Src_PriorType(25,2) / 0 /
        data Src_Prior(25,2,1),Src_Prior(25,2,2) / 3.603959d+01, 0.001389 /
        data Src_PriorType(25,3) / 0 /
        data Src_Prior(25,3,1),Src_Prior(25,3,2) / 2.390000d-04 , 9.560000d-05 /
        data nu0(25) / 1.581162d+10 /
        data Src_PriorType(25,4) / 0 /
        data Src_Prior(25,4,1),Src_Prior(25,4,2) / 5.000000d-01 , 0.000000d+00 /


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! No. of parameters per object:

	parameter(NPars=NGeoPars+NMassPars*max(GL,MMass)+SZ*(NGasPars*Gas+NTPars*Temperature)+1)

! No. of nuisance parameters:

	parameter(NNuisance = SZ*4*NSrc + Varyzs)

! Total no. of dimensions of parameter space:

	integer NDim,tot_dim,edim,efast,eslow
	parameter(NDim = NAtoms*NPars + NNuisance)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Sampling parameters:
!^^^^^^^^^^^^^^^^^^^^^

! Useful statistics:
	double precision GLLhood,GLLhood0,SZLhood0,GLNullEv,NullEv,Lbar

! ??? Visualise degeneracies in prior by sampling with no data?
	logical SamplePrior
	parameter (SamplePrior=.false.)

! ??? Debugging (may not be very helpful...):
	logical verbose
      parameter(verbose=.false.)




! ??? Use Nested Sampling?

      ! M for MultiNest, P for PolyChord
      character*1 which_sampler
      parameter(which_sampler='P')

      !whether to Importance Nested Sampling
      logical n_IS
 	parameter(n_IS=.false.)

      !whether to do multimodal sampling
      logical n_mmodal
 	parameter(n_mmodal=.true.)

      !constant efficiency sampling
      logical n_ceff
 	parameter(n_ceff=.false.)

      !max no. of live points
      integer nest_nlive
	parameter(nest_nlive=100)

      ! number of repeats (for PolyChord) - set to 0 or leave out for default based on number of slow parameters
      integer nest_nrep
        parameter(nest_nrep=0)

      !required efficiency
      !relevant only for the ellipsoidal sampling
      double precision n_efr
	parameter(n_efr=0.8)

      !the following parameters should best be left to their default values in
      !most cases

      !seed for nested sampler, -ve means take it from sys clock
	integer n_rseed, n_totPar
	parameter(n_rseed=-1)

      !min no. of points per cluster
      integer n_minp

      !evidence tolerance factor
      double precision n_tol
      parameter(n_tol=0.5)

      !root for saving posterior files
	character*100 n_root

  parameter(n_root='OUTDIRECTORY')


      integer, allocatable :: n_pWrap(:)

      !feedback on the sampling progress?
      logical n_fb
      parameter(n_fb=.true.)

      !max modes expected, for memory allocation
      integer n_maxModes
      parameter(n_maxModes=20)


      ! .....Sparse  matrix  and Incomplete Cholesky decomposition  variables and working arrays


      	INTEGER    	                               :: diagflag,dsflag

	INTEGER                                        :: n_irn, n_icn, n_samat

	INTEGER                                         :: n_ichol, nz_ichol


      	INTEGER, DIMENSION(:),   ALLOCATABLE :: irn,icn
	REAL *8 ,  DIMENSION(:),   ALLOCATABLE ::samat(:)

	REAL*8                                   ::  lndet,  lndet_ichol
        INTEGER                                 :: nmax_ichol, nzmax_ichol
        INTEGER                                 :: lenw
	INTEGER                                 :: leniw
	INTEGER                                 :: isym,iter,ierr,iunit,itmax
	INTEGER                                 :: itol
	INTEGER                                 :: locdin,locdz,locel,lociel,lociw,locjel,locp,locr
	INTEGER                                 :: locw,locz,locrb, locib,nl
	REAL*8                                   :: err, tol
      	INTEGER, DIMENSION(:),   ALLOCATABLE ::iwork
        REAL*8,  DIMENSION(:),   ALLOCATABLE :: rwork


!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Common blocks:

	common /Work/ Q(2,2),r(n),Phi(n),T(n),Rhogas(n),Pgas(n),Yarray(n)
	common /Work/ logr(n),logPhi(n),logT(n),logRhogas(n),logPgas(n),logYarray(n)
	common /Projection/ uu,loguu , T200 , Mgas200 ,fgas200
      common /dim/ tot_dim,edim,tot_atoms,efast,eslow

	common /GravLens/ x(Ngals),y(Ngals), &
      e1(Ngals),e2(Ngals),e1err(Ngals),e2err(Ngals), &
      gamma1(Ngals),gamma2(Ngals),kappa(Ngals),z_s(Ngals), &
      g1(Ngals),g2(Ngals),wt1(Ngals),wt2(Ngals),GLeflag(Ngals)

	common /GravLensSurvey/ e1s(nxpix,nypix),e2s(nxpix,nypix), &
      gamma1s(nxpix,nypix),gamma2s(nxpix,nypix),g1s(nxpix,nypix),g2s(nxpix,nypix), &
      kappas(nxpix,nypix),e1errs(nxpix,nypix),e2errs(nxpix,nypix),wt1s(nxpix,nypix), &
     	wt2s(nxpix,nypix),z_ss(nxpix,nypix),GLeflags(nxpix,nypix),ptInPix(nxpix,nypix,50)

	common /SZ/  u(Nvisfiles,large),v(Nvisfiles,large), &
      visr(Nvisfiles,large),visi(Nvisfiles,large), &
      visrms(Nvisfiles,large),viswt(Nvisfiles,large), &
      pvisr(Nvisfiles,large),pvisi(Nvisfiles,large), &
      pvisr1(Nvisfiles,large),pvisi1(Nvisfiles,large), &
     LCM(large*Nvisfiles*(2*large*Nvisfiles+1)),SZscale(Nvisfiles),SZscale_BA(Nvisfiles), &
      ymap,trans(6),uvcell,uvtrans(6), &
      pb_x(Nvisfiles),pb_y(Nvisfiles), &
      pb_x_offset(Nvisfiles),pb_y_offset(Nvisfiles), &
      ph_x(Nvisfiles),ph_y(Nvisfiles), &
      ph_x_offset(Nvisfiles),ph_y_offset(Nvisfiles), &
      telescope(Nvisfiles), &
      SZeflag(Nvisfiles,large),aux(NAtoms,aux_dim)

	common /Object/ GeoPars(NAtoms,NGeoPars),GeoPars_old(NAtoms,NGeoPars), GeoPars_old1(NAtoms,NGeoPars), &
      MassPars(NAtoms,NMassPars), MassPars_old(NAtoms,NMassPars), MassPars_old1(NAtoms,NMassPars), &
      GasPars(NAtoms,NGasPars), GasPars_old(NAtoms,NGasPars), GasPars_old1(NAtoms,NGasPars),&
      TPars(NAtoms,NTPars), TPars_old(NAtoms,NTPars), TPars_old1(NAtoms,NTPars)

	common /Nuisance/ SrcPars(NSrc,4),zs,z,z_old,z_old1,D

	common /Statistics/ GLLhood,GLLhood0,SZLhood0,OldLhood,NullEv,Lbar, &
      Shots,Bricks,Airballs,Hoops, &
      Ns_afterburn,burnin,Evidfile

	common /Plotting/ Geo_Ranges,Mass_Ranges,Gas_Ranges,T_Ranges, &
      Src_Ranges,zs_Ranges,dev

      common /solvers/ r200,M200,r500,rs,ps
      common/solvers2/rc,beta,A1,Gm,rhocritz,Rhogas0 ,logRhogas0

      common/derived/Ysph500,Ysph200,Ycyl500,Ycyl200

      common /solve/ rE,scE

      common /solverc/ rp,Mproj,slope

      common /betamodel/ a,b,c

      common /fits_buff/ buffer,ibuff,nblock

      common /fits_file/ ifile, filep

      common /fitsvars/ SIMPLE,BITPIX,NAXIS,NAXES,BKSIZE, &
                       PCOUNT,GCOUNT,EXTEND,OBJECT,TELESCOP,INSTRUME, &
                       OBngalSRA,OBSDEC,CRVAL4,CRVAL5,CRVAL6, &
                       PSCAL1,PSCAL2,PSCAL4,PZERO4,PSCAL5,PZERO5, &
                       nu(Nvisfiles),pb_sig(Nvisfiles),nu0(0:NSrc),flux0(NSrc)

      common arr(nx,ny)
      common fftwplan
! GNFW_Planck common blocks

      COMMON /PlanckPars/  thetas_Planck ,&
                          c_GNFW, a_GNFW, b_GNFW, c500_GNFW, &
                           Y500_Planck  , y0_Planck, &
			  ycoeff_Planck  , map_sum_Planck , lim1 , lim2, &
			  theta500_Planck,   y0_int_Planck ,&
			  Ytot_Planck ,Ytot_int_Planck

      COMMON /Planckarrs/ yintegrand_Planck(n), logyintegrand_Planck(n), &
      theta_Planck(n),logtheta_Planck(n),yarray_Planck(n),logyarray_Planck(n)

!Beta_Atomic/Blobology common blocks

     COMMON /BetaAtomic/ thetac_BA , Beta_BA , dT0_BA ,x0_BA,y0_BA


 ! DM-GNFW model common blocks

      COMMON /DM_GNFWPars/rs_DM , rhos_DM, Pei_GNFW, rp_GNFW , c200_DM ,&
                           MT200_DM,fg200_DM, r200_DM, Mg200_DM, Tg200_DM, &
			   rhogasi_DM, nei_DM,Tei_DM, r500_DM,MT500_DM,Mg500_DM,fg500_DM,Tg500_DM, &
			   y0_GNFW,ycoeff_GNFW, Ysph500_GNFW,Ysph200_GNFW, &
			   Ysph500_GNFW_ang,Ysph200_GNFW_ang ,Ycyl_GNFW, Ycyl_GNFW_ang, &
			   map_sum_GNFW

      COMMON /DM_GNFWarrs/Pe_GNFW(n),logPe_GNFW(n),yintegrand_GNFW(n),logyintegrand_GNFW(n), &
                           yarray_GNFW(n),logyarray_GNFW(n)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	double precision B_A_D
!       parameter(B_A_D=-1.7976931348623157E+308)
      parameter(B_A_D=-1.0d64)

!=======================================================================

end module params
