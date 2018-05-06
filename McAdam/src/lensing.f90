module lensing
	use utilities
      use constants

contains
!=======================================================================

        subroutine lens(eps_01,eps_02,g1,g2,eps1,eps2)
        
        implicit none
        
        double precision eps1,eps2,g1,g2,eps_01,eps_02
        double complex eps,g,eps_0
        
        eps_0 = dcmplx(eps_01,eps_02)
        g = dcmplx(g1,g2)
        
        eps = (eps_0 + g) / (1.0 + conjg(g) * eps_0)
        
        eps1 = dble(eps)
        eps2 = aimag(eps)
        
        return  
        end subroutine lens
        
!=======================================================================

        subroutine unlens(eps1,eps2,g1,g2,eps_01,eps_02)
        
        implicit none
        
        double precision eps1,eps2,g1,g2,eps_01,eps_02
        double complex eps,g,eps_0
        
        eps = dcmplx(eps1,eps2)
        g = dcmplx(g1,g2)
        
        eps_0 = (eps - g) / (1.0 - conjg(g) * eps)
        
        eps_01 = dble(eps_0)
        eps_02 = aimag(eps_0)
        
        return  
        end subroutine unlens
        
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

        subroutine unlensdble(eps1,eps2,g1,g2,eps_01,eps_02)
        
        implicit none
        
        double precision eps1,eps2,g1,g2,eps_01,eps_02
        double precision re,im
        double complex eps,g,eps_0
        
        re = 1.0*eps1
        im = 1.0*eps2
        eps = dcmplx(re,im)
        re = 1.0*g1
        im = 1.0*g2
        g = dcmplx(re,im)
        
        eps_0 = (eps - g) / (1.0 - conjg(g) * eps)
        
        eps_01 = 1.d0*dble(eps_0)
        eps_02 = 1.d0*aimag(eps_0)
        
        return  
        end subroutine unlensdble
        
!=======================================================================
!
      subroutine frwgeom(z1,z2,D1,D2,D12,sigcrit)
!
!  Work out sigma critical (see lab book) and angular diameter distances.
!  f is the integrand of the chi integral, without the factors of c
!  and Rzero. Chi and angdist all take their usual meanings - some
!  interpolation is performed to get the best estimate of the proper 
!  volume at each z. 
!
      implicit none
!	
      double precision z2,z1,D1,D2,D12
	double precision alpha,sigcrit
	integer k,i,j,n1	
	parameter(n1=1000)      
      double precision z(2),g,f(n1),zz(n1),t,chi(2),sofchi(2)
      double precision test,zmin,zmax
!
!-----------------------------------------------------------------------
!
! Initialise variables:
!
	z(1) = z1
	z(2) = z2
!
! Find alpha:
!
	test = Om + Ol - 1.0
	if (test.eq.0.0) then
	  k = 0
	  alpha = 1.0
	elseif (test.lt.0.0) then
	  k = -1
	  alpha = 1.0/sqrt(-test)
	elseif (test.gt.0.0) then
	  k = 1
	  alpha = 1.0/sqrt(test)
	endif 
!
! Now work out chi etc for each redshift:
!
	zmin = 0.0
!
	do j=1,2
!
        zmax = z(j)
!	  
	  do i=1,n1
          zz(i) = zmax*float(i-1)/float(n1-1)
          t = zz(i)
          f(i) = Om*t*(1+t)*(1+t) - Ol*t*(2+t) + (1+t)*(1+t)
	    f(i) = 1.0/sqrt(f(i))
        enddo
!
        chi(j) = integrate(f,zz,n1,zmin,zmax)
	  chi(j) = chi(j)/alpha	 
!
	  if (k.eq.-1) then
	    sofchi(j) = sinh(chi(j))
	  elseif (k.eq.0) then
	    sofchi(j) = chi(j)
	  elseif (k.eq.1) then
	    sofchi(j) = sin(chi(j))
	  endif  
!
	enddo
!
! First calculate distances; c/H = 2.9979e3 Mpc if H = 100 km/s/Mpc:
!
	D1 = alpha*(clight/(h*100.*1000.))*sofchi(1)/(1.0+z1)
	D2 = alpha*(clight/(h*100.*1000.))*sofchi(2)/(1.0+z2)

! Now combine to form g:
!
	g = (1+z1)*sofchi(2)/sofchi(1)
!
	t = chi(2) - chi(1)
	if (k.eq.-1) then
	  t = sinh(t)
	elseif (k.eq.0) then
	  t = t
	elseif (k.eq.1) then
	  t = sin(t)
	endif  
!
	g = g/t
!	
! Use t to return distance between 2 redshifts:
!
	D12 = alpha*(clight/(h*100.*1000.))*t/(1.0+z2)
!	
! Finally combine to return sigma crit in units of  h Msun pc^-2:
!
	sigcrit = 554.4*g/alpha
!
      return
      end subroutine frwgeom
!	
!=======================================================================
	
	subroutine HDFnz(maglim,z1,n1,meanz)
	
! Read in catalogue from http://bat.phys.unsw.edu.au/~fsoto/hdfcat.html
! (published in ApJ 1999, 513, 34) and work out some statistics

	implicit none
      
	double precision maglim,z1,n1,meanz
	
	integer ngals
	parameter(ngals=1067)
	double precision z(ngals),ab(ngals),xpos(ngals),ypos(ngals)
      integer i,j,flag(ngals)

	character*100 filename
	integer idum1,idum2
	character cdum1*41,cdum2*7

	double precision ablim	
	double precision count,sum
      double precision median(ngals)
	
	
!-----------------------------------------------------------------------
	
	filename = 'photz_t4.dat'
			
	open(unit=23,file=filename,form='formatted',status='old')
	do i=1,19
	  read(23,*)
	enddo
	do i=1,ngals
	  read(23,5) idum1,xpos(i),ypos(i),cdum1,ab(i),cdum2,z(i),idum2
	enddo
	close(23) 
 5    format(i4,1x,f6.1,1x,f6.1,a41,f5.2,a7,f5.3,1x,i1)	

! Discard galaxies from outer regions - should be 121:

	do i=1,ngals
	  if (xpos(i).lt.300.0.or.xpos(i).gt.3800.0.or. &
        ypos(i).lt.300.0.or.ypos(i).gt.3800.0) then
          z(i) = -1.0
	  elseif (xpos(i).gt.1900.0.and.ypos(i).gt.2000.0) then
          z(i) = -1.0
	  endif
	enddo
		
      do j=1,ngals
	  flag(j) = 0
	  if (z(i).lt.0.0) flag(i) = 1
	enddo  

! Separate out galaxies with ab less than ablim
! AB_lambda = -2.5 log (S_v / Jy) + 8.90
! lambda = 8140 Angstroms
!
! Peacock gives Mould I band mag = 8.45 - 2.5 log (S_v / Jy)
! with effective wavelength 8300 A
!
! Assume S_v is constant, then mag_I = AB - 0.45

	ablim = maglim + 0.45
			
	do i=1,ngals
	  if (ab(i).gt.ablim) then
	    flag(i) = 1
	  endif  
	  if (z(i).le.z1) then
	    flag(i) = 1
	  endif  
	enddo

	count = 0.0
	sum = 0.0
	do i=1,ngals
	  if (flag(i).ne.1) then
	    count = count + 1.0
	    sum = sum + z(i)
	    median(i) = dble(z(i))
	  else
	    median(i) = 0.0
	  endif  
	enddo

      call sort(median,ngals)
	
	meanz = median(nint(count/2.0))
	
!	if (count.gt.0.0) then
!	  meanz = sum/count
!      else
!	  meanz = 0.0
!	endif  	
	
	n1 = count/3.92
	
	return
	end subroutine HDFnz
	
!=======================================================================
! 
!       subroutine sort(A,n)
! 
! !  Sort array into descending order.
! 
!       implicit none
! 
!       integer n
!       real A(n)
! 
!       real*8 AA(n)
!       integer rank(n)
!       integer ifail,i
!       external  m01daf,m01eaf
! 
! !-----------------------------------------------------------------------
! 
!       do i=1,n        
!         rank(i) = 0
!         AA(i) = dble(A(i))
!       enddo  
! 
! !  Rank the data in order of descending S:
! 
!       ifail = 0
!       call m01daf(AA,1,n,'d',rank,ifail)
! 
! !  Now rearrange data into sorted list:
! 
!       call m01eaf(AA,1,n,rank,ifail)
! 
!       do i=1,n        
!         A(i) = float(AA(i))
!       enddo  
! 
!       return
!       end
!                
! !=======================================================================
! 
!       subroutine sort2(A,B,n)
! 
! !  Sort 2 arrays into descending order.
! 
!       implicit none
! 
!       integer n
!       real A(n),B(n)
! 
!       real*8 AA(n),BB(n)
!       integer rank(n)
!       integer ifail,i
!       external  m01daf,m01eaf
! 
! !-----------------------------------------------------------------------
! 
!       do i=1,n        
!         rank(i) = 0
!         AA(i) = dble(A(i))
!         BB(i) = dble(B(i))
!       enddo  
! 
! !  Rank the data in order of descending S:
! 
!       ifail = 0
!       call m01daf(AA,1,n,'d',rank,ifail)
! 
! !  Now rearrange data into sorted list:
! 
!       call m01eaf(AA,1,n,rank,ifail)
!       call m01eaf(BB,1,n,rank,ifail)
! 
!       do i=1,n        
!         A(i) = float(AA(i))
!         B(i) = float(BB(i))
!       enddo  
! 
!       return
!       end
!                
!=======================================================================
      SUBROUTINE SORT(RA,n1)
      INTEGER n1
      DOUBLE PRECISION, DIMENSION(n1) :: RA
      INTEGER L, IR, I, J
      DOUBLE PRECISION RRA
      L=n1/2+1
      IR=n1
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            goto 30
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
 30   call flip(ra,n1)	
      return
	END SUBROUTINE SORT
      
      
      SUBROUTINE SORT2(n1,RA,RB)
      INTEGER n1
      DOUBLE PRECISION, DIMENSION(n1) :: RA,RB
      INTEGER L, IR, I, J
      DOUBLE PRECISION RRA, RRB
      L=n1/2+1
      IR=n1
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
 30   call flip(ra,n1)	
      call flip(rb,n1)
      return	
      END SUBROUTINE SORT2
      
      
      subroutine flip(x,n1)
      integer i,j,n1
      double precision x(n1),temp
      do i=1,n1/2
        j = n1-(i-1)
        temp = x(i)
        x(i) = x(j)
        x(j) = temp
      enddo
      return
      end subroutine flip
!=======================================================================
end module lensing
