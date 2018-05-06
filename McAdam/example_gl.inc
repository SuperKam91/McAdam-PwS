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
	double precision pvisr,pvisi
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
	parameter (aux_dim=1)        ! for Beta_Atomic/blobology model		
	double precision r,rmin,rmax,rlimit,Phi,T,Rhogas,Pgas,Yarray,Rhogas_central,sn_lim,aux
	double precision logr,logPhi,logT,logRhogas,logPgas,logYarray
	double precision uu,loguu , T200 , Mgas200 , fgas200
	parameter(rmin=0.001,rmax=100.0,rlimit=20.0)
      parameter(sn_lim=0.98)!lensing S/N limit
      double precision r200,M200,r500
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
       REAL *8   , PARAMETER    ::  a_GNFW=1.0620d0        ! The constants in GNFW pressure profile: alpha, beta, gamma, Arnaud et al. 2009
       REAL *8   , PARAMETER    ::  b_GNFW=5.4807d0
       REAL *8   , PARAMETER    ::  c_GNFW=0.3292d0
       REAL *8   , PARAMETER    ::  c500_GNFW= 1.156          ! Gas concentration parameter, Arnaud et al. 2009            
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
	double precision GeoPars
	integer NMassPars
	double precision MassPars
	integer NGasPars
	double precision GasPars	
	integer NTPars
	double precision TPars
	double precision z(NAtoms),D

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
! ??? Which gas density profile model?
!
!	parameter(GasModel=1)
!	parameter(NGasPars=3)     ! for betastayle(1)=2
!	parameter(NGasPars=4)     ! for betastyle(1)=4

!	parameter(GasModel=3)     ! for GNFW_Planck model
!	parameter(NGasPars=2)	  !two gas parameters for GNFW_Planck model: thetas(arcmin) and Ytot(arcmin^2)

        parameter(GasModel=4)   !  Beta_Atomic/blobology Model
        parameter(NGasPars=3)   !  parameters are angular core radius, beta and central temperature decrement
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
	parameter(NSrc=3)
	double precision SrcPars
				
!   Integer flag for varying the redshift of the lensed galaxies?
	integer Varyzs,Vary_zs
	double precision zs,zdmin,zdmax

! ??? How many source planes are there? Relevant only if GL=1 & Varyzs=1
	integer nzsplanes
	parameter(nzsplanes=2)
	

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
      parameter(GL=1)
! ??? Where is the weak lensing data?
      parameter(GLdatafile='data/GL/sim/090707/2.cat')
! ??? How many galaxies in catalogue?
      parameter(Ngals=5760)
! ??? Allow redshift of sources to vary (for prediction)?
      parameter(Varyzs=1)
      
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
      parameter(SZ=0)
      parameter(Gas=0)
!     parameter(Temperature=0)          ! when sampling from Temperature
      parameter(Temperature=0)          ! when not sampling from Temperature
! ??? Model Mass
      parameter(MMass=0)
! ??? How many SZ data files?      
      parameter(Nvisfiles=0)
! ??? Where are the data?
      dimension visdatafile(Nvisfiles)
 


! ??? How many uv points in each file?
      dimension Nvis(Nvisfiles)

! For minimal memory usage, set large >= the largest of the Nvis values
! This can be very wasteful but is needed to allow easy array
! allocation...
! ??? Maximum vis data array size?
      parameter(large=0)
! 128x128 map is used in FFT - need it big enough to cover x,y priors
! ??? Comptonisation map pixel scale / arcsec?
	parameter(cell=30.0)
      
! Need to specify which files have an associated covariance matrix,
! setting IncludeCMB=0 neglects CMB and just uses 6th column rms
! reciever noise values in diagonal chisquared computation:
! ??? Include CMB in likelihood via a non-diagonal covariance matrix? 
      parameter(IncludeCMB=0)
      
! ??? Matrix files location? - watch out for byte-swapping...


! ??? Centre of coordinate system (RA, Dec in degrees)? 
! x,y are measured in arcsec from this, in right handed sense
      parameter(map_x=60.d0,map_y=40.d0)

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
	parameter(mass_function = 2)

! Example mock data = NFW mass profile + gas in hydrostatic equilibrium
! Model = Isothermal beta model in cored isothermal sphere potential.

! Clump 1:
 
! ??? Geo		          
      data Geo_PriorType(1,1) / 3 /
      data Geo_Prior(1,1,1),Geo_Prior(1,1,2) / 0.d0, 60.d0 /
      data Geo_PriorType(1,2) / 3 /
      data Geo_Prior(1,2,1),Geo_Prior(1,2,2) / 0.d0, 60.d0 /
 
! ??? NFW profile
      data NFWstyle(1) / 2 /
      data Mass_PriorType(1,1) / 1 /
      data Mass_Prior(1,1,1),Mass_Prior(1,1,2) / 0.1, 15. /
      data Mass_PriorType(1,2) / 3 /
      data Mass_Prior(1,2,1),Mass_Prior(1,2,2) / 5.d14,2.5d14 /
      
! cluster redshift
      data z_PriorType(1) / 1 /
      data z_Prior(1,1),z_Prior(1,2) / 0.01d0, 2d0 /
      
! source redshift, relevant only if varyzs = 1
      data zs_PriorType(1) / 1 /
      data zs_Prior(1,1),zs_Prior(1,2) / 0.01d0, 2d0 /
      data zs_PriorType(2) / 1 /
      data zs_Prior(2,1),zs_Prior(2,2) / 0.01d0, 2d0 /

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! No. of parameters per object:

	parameter(NPars=NGeoPars+NMassPars*max(GL,MMass)+SZ*(NGasPars*Gas+NTPars*Temperature)+1)

! No. of nuisance parameters:

	parameter(NNuisance = SZ*4*NSrc)
      
! Total no. of dimensions of parameter space:

	integer NDim,tot_dim,edim
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
      
      !whether to Importance Nested Sampling
      logical n_IS
 	parameter(n_IS=.false.)
      
      !whether to do mode separation
      logical n_mmodal 
 	parameter(n_mmodal=.false.)
      
      !constant efficiency sampling
      logical n_ceff
 	parameter(n_ceff=.false.)
	
      !max no. of live points
      integer nest_nlive
	parameter(nest_nlive=1000)
      
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

  parameter(n_root='chains/example_gl_varyzs')

      
      integer, allocatable :: n_pWrap(:)
      
      !feedback on the sampling progress?
      logical n_fb 
      parameter(n_fb=.true.)
      
      !max modes expected, for memory allocation
      integer n_maxModes 
      parameter(n_maxModes=1)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Common blocks:
	
	common /Work/ Q(2,2),r(n),Phi(n),T(n),Rhogas(n),Pgas(n),Yarray(n)
	common /Work/ logr(n),logPhi(n),logT(n),logRhogas(n),logPgas(n),logYarray(n)
	common /Projection/ uu,loguu , T200 , Mgas200 ,fgas200
      common /dim/ tot_dim,edim,tot_atoms

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
      LCM(large*Nvisfiles*(2*large*Nvisfiles+1)),SZscale(Nvisfiles),SZscale_BA(Nvisfiles), &
      ymap,trans(6),uvcell,uvtrans(6), &
      pb_x(Nvisfiles),pb_y(Nvisfiles), &
      pb_x_offset(Nvisfiles),pb_y_offset(Nvisfiles), &
      ph_x(Nvisfiles),ph_y(Nvisfiles), &
      ph_x_offset(Nvisfiles),ph_y_offset(Nvisfiles), &
      telescope(Nvisfiles), &
      SZeflag(Nvisfiles,large),aux(NAtoms,aux_dim)
          
	common /Object/ GeoPars(NAtoms,NGeoPars), &
      MassPars(NAtoms,NMassPars), &
      GasPars(NAtoms,NGasPars), &
      TPars(NAtoms,NTPars)

	common /Nuisance/ SrcPars(NSrc,4),zs,z,D
     
	common /Statistics/ GLLhood,GLLhood0,SZLhood0,OldLhood,NullEv,Lbar, &
      Shots,Bricks,Airballs,Hoops, &
      Ns_afterburn,burnin,Evidfile
	
	common /Plotting/ Geo_Ranges,Mass_Ranges,Gas_Ranges,T_Ranges, &
      Src_Ranges,zs_Ranges,dev
      
      common /solvers/ r200,M200,r500
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
