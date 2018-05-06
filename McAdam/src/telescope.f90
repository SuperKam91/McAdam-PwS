module telescope1	
	
      implicit none

! Include file for SZ data analysis / FITS file manipulation
      
! Primary beam width (Gaussian sigma) in arcsec:
! CBI pb is at 31 GHz
! VSA pbs are at 34 GHz
 
	double precision pb_vsa_ext
      !parameter(pb_vsa_ext=3229.0)
      parameter(pb_vsa_ext=3134.0)
	double precision pb_vsa_sea
      !parameter(pb_vsa_sea=1682.0)
      parameter(pb_vsa_sea=1834.54) 
	!double precision pb_ami(6)
	!data pb_ami/549.087,527.429,507.810,489.719,473.412,458.124/
      !data pb_ami/545.265,517.236,512.141,501.949,476.470,430.607/
      double precision pb_ami(2)
      data pb_ami /101.1,1.89/
        !double precision pb_la(6)
        !data pb_la/155.136,149.609,144.622,140.098,135.977,132.207/
      double precision pb_la(2)
      data pb_la /24.905,0.79/
	double precision pb_ryl
      parameter(pb_ryl=153.0)
	double precision pb_cbi
      parameter(pb_cbi=1152.0)
      !PB cut-ff in meters
      double precision lmin_ami(6)
      data lmin_ami/542.93,586.,586.,586.,586.,586./
      double precision lmax_ami(6)
      data lmax_ami/6667.78,7496.,7496.,7496.,7496.,7496./

! Visibility data limits:	
	integer maxvis
      parameter (maxvis = 200000)
      integer max_buff_size
      parameter (max_buff_size = 2101248)

	logical BLOCKED 
	parameter(BLOCKED=.TRUE.)
	logical GROUPS 
	parameter(GROUPS=.TRUE.)

	character*100 BUNIT
	parameter(BUNIT='        ')
	double precision BZERO
	parameter(BZERO=0.0)
	double precision BSCALE
	parameter(BSCALE=1.0)
	integer DECIMALS
	parameter(DECIMALS=10)

	double precision EPOCH
	parameter(EPOCH=1950.0)
	                      
	character*100 CTYPE2                                               
	parameter(CTYPE2 ='COMPLEX ')
  	double precision CRVAL2
	parameter(CRVAL2=1.0)                                               
	double precision CRDEL2
	parameter(CRDEL2=1.0)                                               
	double precision CRPIX2
	parameter(CRPIX2=1.0)                                               
	double precision CROTA2
	parameter(CROTA2=0.0)
	                                               
	character*100 CTYPE3
	parameter(CTYPE3 ='STOKES  ')                
  	double precision CRVAL3
	parameter(CRVAL3=1.0)                                               
	double precision CRDEL3
	parameter(CRDEL3=1.0)                                               
	double precision CRPIX3
	parameter(CRPIX3=1.0)                                               
	double precision CROTA3
	parameter(CROTA3=0.0)
	                                               
	character*100 CTYPE4
	parameter(CTYPE4 ='FREQ    ')                
	double precision CRDEL4
	parameter(CRDEL4=1.0)                                               
	double precision CRPIX4
	parameter(CRPIX4=1.0)                                               
	double precision CROTA4
	parameter(CROTA4=0.0)
	                                               
	character*100 CTYPE5
	parameter(CTYPE5='RA      ')                
	double precision CRDEL5
	parameter(CRDEL5=1.0)                                               
	double precision CRPIX5
	parameter(CRPIX5=1.0)                                               
	double precision CROTA5
	parameter(CROTA5=0.0)
	                                               
	character*100 CTYPE6
	parameter(CTYPE6='DEC     ')                
	double precision CRDEL6
	parameter(CRDEL6=1.0)                                               
	double precision CRPIX6
	parameter(CRPIX6=1.0)                                               
	double precision CROTA6
	parameter(CROTA6=0.0)
	                                               
	character*100 PTYPE1
	parameter(PTYPE1='UU      ')                
	double precision PZERO1
	parameter(PZERO1=0.0)                                               
	character*100 PTYPE2
	parameter(PTYPE2='VV      ')                
	double precision PZERO2
	parameter(PZERO2=0.0)                                               
	character*100 PTYPE3
	parameter(PTYPE3='WW      ')                
	double precision PSCAL3
	parameter(PSCAL3=1.000000000E-10)                                               
	double precision PZERO3
	parameter(PZERO3=0.0)
	character*100 PTYPE4
	parameter(PTYPE4='DATE    ')                
	character*100 PTYPE5
	parameter(PTYPE5='DATE    ')                
	character*100 PTYPE6
	parameter(PTYPE6='BASELINE')                
	double precision PSCAL6
	parameter(PSCAL6=1.0)                                               
	double precision PZERO6
	parameter(PZERO6=0.0)



end module telescope1
