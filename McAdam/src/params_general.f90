INCLUDE "mkl_dss.f90" 
module params

  use mkl_dss

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
!     5 = DM_GNFW model: NFW profile for DM, GNFW model for gas pressure.  Parameters are: MT_200, fgas_200, gamma, alpha, beta, c_500
!
!  Temperature models:
!
!     1 = Isothermal: parameter T / keV
!     2 = Polytropic: parameter T0 / keV & gamma
!
!   Sources
!     x,y,flux,spectral index all to be fitted
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
!     12 = Prior for spectral indices from Waldram 07 multipled by a Gaussian
!     13 = Ratio prior on point source flux
!     14 = 2D elliptical Gaussian prior in log space on theta_s and Y_{tot} (for Planck GNFW model), (centre, width, angle) for both parameters
!
!=======================================================================

! MPI Info
    ! flag to show if node is root or not
    logical is_root

! Data:

	integer ModelClass
	integer GL,SZ
	logical simulate,survey,arrGal

! Data filenames:
      integer Nvisfiles
	character*100 GLdatafile,covmatfile
    integer, parameter :: STR_LENGTH = 1000
    character(len=STR_LENGTH), dimension(:), allocatable :: visdatafile

! GL data arrays:
	integer Ngals
	double precision, dimension(:), allocatable ::  x,y
	double precision, dimension(:), allocatable ::  e1,e2,e1err,e2err,wt1,wt2
	integer, dimension(:), allocatable :: GLeflag
	double precision, dimension(:), allocatable ::  gamma1,gamma2,kappa,g1,g2
	double precision sigmaobs
	parameter(sigmaobs=0.20d0)
	logical errcorrect
	parameter(errcorrect=.false.)
	double precision sigcritE
      double precision, dimension(:), allocatable :: z_s
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
      double precision, dimension(:,:), allocatable ::  e1s,e2s,g1s,g2s,kappas,gamma1s,gamma2s,e1errs,e2errs,wt1s,wt2s,z_ss
      integer, dimension(:,:), allocatable :: GLeflags
      integer, dimension(:,:,:), allocatable :: ptInPix

! SZ data arrays/maps:
	integer, allocatable :: Nvis(:)
        integer large
	double precision, dimension(:,:), allocatable :: u,v,visr,visi,visrms,viswt
        double precision, dimension(:), allocatable :: SZscale,SZscale_BA
	integer, dimension(:,:), allocatable :: SZeflag
	integer nx,ny
	parameter(nx=256,ny=256)
	double precision ymap(nx,ny),cell,trans(6),uvcell,uvtrans(6)
	double precision, dimension(:,:), allocatable :: pvisr,pvisi,pvisr1,pvisi1
	double precision map_x,map_y
        double precision, dimension(:), allocatable :: pb_x,pb_y,pb_x_offset,pb_y_offset
        double precision, dimension(:), allocatable :: ph_x,ph_y,ph_x_offset,ph_y_offset
	character(len=100), dimension(:), allocatable :: telescope

! Planck data information
      integer PLPatch
      integer PL
      character*100 PLdatafile

! hyperparameter variables/ arrays

	integer hyperparameters
	integer Nhyper
	integer, dimension(:,:), allocatable :: hyper_PriorType
        double precision, dimension(:,:,:), allocatable :: hyper_Prior
        double precision, dimension(:,:), allocatable :: hyperPars
        double precision, dimension(:,:), allocatable :: hyperPars_old
	double precision :: alpha_L1, alpha_L2

! Covariance matrix working arrays/variables
        TYPE(MKL_DSS_HANDLE) :: handle ! Allocate storage for the sparse matrix solver handle.
	double precision, dimension(:), allocatable :: LCM
        integer, dimension(:), allocatable     :: irn,icn
        double precision,  dimension(:), allocatable     :: samat
	integer IncludeCMB, dsflag

! Cluster data and working arrays
	double precision Q(2,2)
	integer n
	parameter(n=256)
	integer aux_dim               !number of derived parameters

	double precision r(n),rmin,rmax,rlimit,Phi(n),T(n),Rhogas(n),Pgas(n),Yarray(n),Rhogas_central,sn_lim
        double precision, dimension(:,:), allocatable ::  aux
	double precision logr(n),logPhi(n),logT(n),logRhogas(n),logPgas(n),logYarray(n)
	double precision uu,loguu , T200 , Mgas200 , fgas200
	parameter(rmin=0.0001d0,rmax=20.0d0,rlimit=20.0d0)
      parameter(sn_lim=0.98d0)!lensing S/N limit
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
       double precision                  ::  a_GNFW                 ! The constants in GNFW pressure profile: alpha, beta, gamma
       double precision                  ::  b_GNFW
       double precision                  ::  c_GNFW
       double precision                  ::  c500_GNFW                      ! Gas concentration parameter
       double precision                  ::  thetamin_PLanck             !arcmin
       double precision   , PARAMETER    ::  thetamax_Planck=20.0d0      !arcmin
       double precision   , PARAMETER    ::  thetalimit_Planck=20.0d0     !arcmin


       double precision                  ::  thetas_Planck
       double precision                  ::  theta500_Planck
       double precision                  ::  Y500_Planck
       double precision                  ::  Ytot_Planck
       double precision                  ::  Ytot_int_Planck
       double precision                  ::  y0_Planck , y0_int_Planck
       double precision                  ::  ycoeff_Planck
       double precision                  ::  map_sum_Planck

! arrays for GNFW_Planck model
       double precision                 ::   yintegrand_Planck(n) , logyintegrand_Planck(n)
       double precision                 ::   theta_Planck(n) , logtheta_Planck(n)
       double precision                 ::   yarray_Planck(n) , logyarray_Planck(n)

! Cluster data and working arrays for Beta_Atomic model
       double precision                 ::  thetac_BA , Beta_BA , dT0_BA
       double precision                 ::  x0_BA , y0_BA

! cluster data and working arrays for DM_GNFW model
       double precision                  ::  rs_DM , rhos_DM, Pei_GNFW, rp_GNFW
       double precision                  ::  c200_DM
       double precision                  ::  MT200_DM,fg200_DM
       double precision                  ::  r200_DM, Mg200_DM, Tg200_DM
       double precision                  ::  rhogasi_DM, nei_DM,Tei_DM
       double precision                  ::  r500_DM,MT500_DM,Mg500_DM,fg500_DM,Tg500_DM
       double precision                  ::  y0_GNFW,ycoeff_GNFW, Ysph500_GNFW,Ysph200_GNFW
       double precision                  ::  Ycyl_GNFW
       double precision                  ::  map_sum_GNFW

! arrays for DM-GNFW model
       double precision                  ::  Pe_GNFW(n),logPe_GNFW(n),yintegrand_GNFW(n),logyintegrand_GNFW(n)
       double precision                  ::  yarray_GNFW(n),logyarray_GNFW(n)

! additional parameters for Ein_DM model
       double precision			 ::  aEin_DM, rm2_DM, rhom2_DM !adapted for GM=6 kj 01/02/17


!  File control
      integer ifile, filep

! McAdam workspace:
	double complex arr(nx,ny)
      integer*8 fftwplan

! FITS header definitions:
	logical SIMPLE,EXTEND
      integer BITPIX,BKSIZE,NAXIS,NAXES(7),PCOUNT,GCOUNT
      !character*100 OBJECT,TELESCOP,INSTRUME
      !double precision OBSRA,OBSDEC,CRVAL4,CRVAL5,CRVAL6,PSCAL1,PSCAL2,PSCAL4,PZERO4,PSCAL5,PZERO5
      double precision, dimension(:), allocatable :: nu,nu0,flux0,pb_sig,lmin,lmax

! Ouput sample filenames: Not actually used?
!	character*100 root

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Object parameters:

! ??? Are there objects to be fitted?
	integer :: Atoms

! ??? How many?
	integer :: NAtoms,tot_atoms

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
	double precision, dimension(:,:), allocatable :: GeoPars, GeoPars_old
	integer NMassPars
	double precision, dimension(:,:), allocatable :: MassPars, MassPars_old
	integer NGasPars
	double precision, dimension(:,:), allocatable :: GasPars, GasPars_old
	integer NTPars
	double precision, dimension(:,:), allocatable :: TPars, TPars_old
	double precision, allocatable :: z(:),z_old(:)
        double precision D

! Available Model Classes:
! ^^^^^^^^^^^^^^^^^^^^^^^
!     1 = Clusters

! Cluster Models
! ^^^^^^^^^^^^^^

! Gas profile styles:
!   different atoms may require different styles:
!    I. Style for Beta model- BetaStyle
!    II. Style for Temperature - TStyle
!   Beta model, GasPars(3) is gas mass within rmass (BetaStyle=1,TStyle=0 ) or r200 in HSE (BetaStyle=2,TStyle=0 )
!   or gas fraction within r200 (BetaStyle=3,TStyle=0)
!   or total mass within r200 and gas mass fraction within r200 in Beta_Virial (BetaStyle=4 ,TStyle=1 )
!   or total mass within r200 and gas mass fraction within r200 in Beta_HSE (BetaStyle=4 ,TStyle=2 )
	integer, allocatable :: BetaStyle(:) , TStyle(:)
      double precision rmass
      parameter(rmass=1.0d0)

! Mass profile styles:
!   Weak lensing only use 2(NFW), 1(SIS), 1(CPL)
!   CPLstyle=-1 is needed for CIS profile with M200 calculated.
!   Weak lensing + Einstein radius/arcsec use 4,2,2
!   Obviously different atoms require different styles
	integer, allocatable :: NFWstyle(:)
        integer, allocatable :: SISstyle(:),CPLstyle(:)


! Nuisance parameters:
!^^^^^^^^^^^^^^^^^^^^^

! ??? Are there nuisance parameters to be fitted? Set automatically
	integer Nuisance

!   Number of nuisance parameters - set automatically:
	integer NNuisance

!   Integer flag for presence of point sources in the radio data?
	integer SourceSubtract
! ??? How many sources are there?
	integer NSrc
	double precision, dimension(:,:), allocatable :: SrcPars

!   Integer flag for varying the redshift of the lensed galaxies?
	integer Varyzs,Vary_zs
	double precision zs,zdmin,zdmax

! ??? How many source planes are there? Relevant only if GL=1 & Varyzs=1
	integer nzsplanes


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

	integer, allocatable :: Geo_PriorType(:,:)
	double precision, dimension(:,:,:), allocatable :: Geo_Prior, Geo_Tri_Prior

	integer, allocatable ::  Mass_PriorType(:,:)
	double precision, allocatable ::  Mass_Prior(:,:,:)

	integer, allocatable :: Gas_PriorType(:,:)
	double precision, allocatable :: Gas_Prior(:,:,:)

	integer, allocatable :: T_PriorType(:,:)
	double precision, allocatable :: T_Prior(:,:,:)

	integer, allocatable :: zs_PriorType(:)
	double precision, allocatable :: zs_Prior(:,:)

	integer, allocatable :: z_PriorType(:)
	double precision, allocatable :: z_Prior(:,:)

! Source arrays:
!^^^^^^^^^^^^^^^

	integer, allocatable :: Src_PriorType(:,:)
	double precision, allocatable :: Src_Prior(:,:,:)

      !spectral index prior array
      double precision, allocatable :: kernel(:,:,:)
      integer, allocatable :: nkernel(:)
      double precision prior_min, prior_max

      ! Numbers of parameters
	integer NDim,tot_dim,edim,eslow

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Sampling parameters:
!^^^^^^^^^^^^^^^^^^^^^

! Useful statistics:
	double precision GLLhood,GLLhood0,SZLhood0,GLNullEv,NullEv,Airballs,Lbar

! ??? Visualise degeneracies in prior by sampling with no data?
	logical SamplePrior

! ??? Debugging (may not be very helpful...):
	logical verbose


! ??? Use Nested Sampling?

      ! M for MultiNest, P for PolyChord
      character*1 which_sampler

      !whether to Importance Nested Sampling
      logical n_IS

      !whether to do multimodal sampling
      logical n_mmodal

      !constant efficiency sampling
      logical n_ceff

      !max no. of live points
      integer nest_nlive

      ! number of repeats (for PolyChord)
      integer nest_nrep

      !required efficiency
      !relevant only for the ellipsoidal sampling
      double precision n_efr

      !the following parameters should best be left to their default values in
      !most cases

      !seed for nested sampler, -ve means take it from sys clock
	integer n_rseed, n_totPar

      !min no. of points per cluster
      integer n_minp

      !evidence tolerance factor
      double precision n_tol

      !root for saving posterior files
	character*100 n_root

      integer, allocatable :: n_pWrap(:)

      !feedback on the sampling progress?
      logical n_fb

      !max modes expected, for memory allocation
      integer n_maxModes

	double precision B_A_D
!       parameter(B_A_D=-1.7976931348623157E+308)
      parameter(B_A_D=-1.0d64)

!=======================================================================

end module params
