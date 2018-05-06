program makeCat
      
      use RandomNS
      
      integer i,j,k,iostat
      double precision e1hat,e2hat,e1,e2,x,y,z,kappa,g1,g2,mode,noise1,noise2
      double precision sig_obs,sig_int,error
      
      sig_obs=0.2
      sig_int=0.25
      
      call initRandomNS(1)
      
      open(unit=23,file='gal.011.01.01.cat',status='unknown')
      open(unit=24,file='white_thousand.cat',status='unknown')
      open(unit=25,file='white_thousand_true.cat',status='unknown')
      
      do
      	read(23,*,IOSTAT=iostatus)x,y,z,kappa,g1,g2
            if(iostatus<0) exit
            do
      		e1hat=Gaussian1NS(0)*sig_int
      		e2hat=Gaussian1NS(0)*sig_int
            	mode=sqrt(e1hat*e1hat+e2hat*e2hat)
      		!deal with bad intrinsic galaxies (hack)
			if(ehat<1.) exit
		enddo
      	call lensdble(e1hat,e2hat,g1,g2,e1,e2)
            
            !alter data inside strong lensed region
		mode=e1*e1+e2*e2
		if(mode>=1.0) then
	   		e1=e1/mode
	   		e2=e2/mode
	 	endif
            
            write(25,'(4x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5)')x,y,kappa,e1,e2,z
            
	 	!add experimental noise:
		noise1=Gaussian1NS(0)*sig_obs
		noise2=Gaussian1NS(0)*sig_obs
	  	e1=e1+noise1
		e2=e2+noise2
            
            error=sqrt(sig_obs**2+sig_int**2)
            
           	write(24,'(4x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5,3x,f8.5)')x,y,e1,e2,error,z
      enddo
      
      close(23)
      close(24)
      close(25)
      call killRandomNS
      
end program makeCat
        
!=======================================================================

        subroutine lensdble(eps_01,eps_02,g1,g2,eps1,eps2)
        
        implicit none
        
        double precision eps1,eps2,g1,g2,eps_01,eps_02
        double precision re,im
        double complex eps,g,eps_0
        
        re = 1.0*eps_01
        im = 1.0*eps_02
        eps_0 = dcmplx(re,im)
        
        re = 1.0*g1
        im = 1.0*g2
        g = dcmplx(re,im)
        
        eps = (eps_0 + g) / (1.0 + conjg(g) * eps_0)
        
        eps1 = dble(eps)
        eps2 = aimag(eps)
        
        return  
        end subroutine lensdble
        
!=======================================================================
