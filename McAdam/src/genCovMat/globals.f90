module globals

	!makeintegrals
      	logical 	imagFlag
	integer  	nzmax
      	integer, dimension(:),   allocatable :: irn,icn
      	double precision,  dimension(:),   allocatable :: freqi,freqj
	integer		intFileUnit !file unit for the file with integrals
      
      	!makecm
      	double precision,  dimension(:),   allocatable :: cm_rr,cm_ii,cm_ri,cm_ir
      
      	!calcints
      	double precision     	ar1,ar2,ar3,ar3_real,ar3_img,b,fmax
     	integer    	nvmax,nfmax,namax
      	double precision,  dimension(:,:), allocatable :: integral1r,integral1i,integral2r,integral2i
      	integer, dimension(:),   allocatable :: freqn
      	double precision,  dimension(:),   allocatable :: u,v
      	double precision,  dimension(:),   allocatable :: visr,visi,rms
      	double precision,  dimension(:),   allocatable :: r
      	double precision,  dimension(:),   allocatable :: cps1d,confps1d,gps1d
      	double precision 	sp,conf_cl
      	double precision,  dimension(:),   allocatable :: freq,sigma,antwdth,f_ra,f_dec

	! ... other variables
      	integer    	itsp,itc,itg,ncps,ngps,nsp,iaddn
      	integer    	nv,nf,na,nz,nskip,diagflag,dsflag
      	double precision     	PI,antmax,uvmin,uvmax,freq0,ra_ref,dec_ref,leff,scaled,scalen
      	double precision     	convf,tcmb,chisqr,lndetr,chisqi,lndeti,ll,normr,normi
      	double precision     	cpsmin,cpsmax,cpsstep
      	double precision     	gpsmin,gpsmax,gpsstep
      	double precision     	spmin,spmax,spstep
      	character  	galyn*1,outyn*1,wryn*1,restart*1
      	character*500, dimension(:),   allocatable :: datafn,noisefn
      	character  	outputfn*500,intfile*500,integralsfn*500
	character*100  outputfn2
      	integer	input_unit
      	parameter	(input_unit=44)
      
      	!constants
      	double precision 	MIN2RAD
      	parameter	(MIN2RAD=3.14159265358979323844d0/180.d0/60.d0)
  
  	!MPI
	integer my_rank
  	integer mpi_nthreads !total no. of mpi processors
      
end module globals
