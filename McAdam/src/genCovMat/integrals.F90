module integrals
	
      use globals
      use numerics
      use utilities

contains

!-----*-----------------------------------------------------------------

subroutine count_integrals(nzmax,nv,nf,freqn,u,v,freq,antwdth,diagflag)
	implicit none
      	integer  nzmax,nv,nf,diagflag
      	integer  freqn(nv)
      	double precision   u(nv),v(nv)
      	double precision   freq(nf),antwdth(nf)

	! ... for independent bins (saves space)
	!      integer  maxbin1(nzmax),maxbin2(nzmax)
	!      real*8   integral1(nzmax,3),integral2(nzmax,3),temp(3)

	! ... for general bins

      	integer  i,j
      	integer  ifr,jfr
      	double precision   ui,uj,vi,vj,anti,antj,amax,sep1,sep2

	nzmax=0
      
	!...loop over data points
      	do j=1,nv
        	jfr=freqn(j)
        	antj=antwdth(jfr)
        	uj=u(j)
        	vj=v(j)
        	do i=j,nv  ! only calculate elements in lower triangle
          		if((diagflag == 0) .and. (i /= j)) cycle
          		ifr=freqn(i)
          		anti=antwdth(ifr)
          		ui=u(i)
          		vi=v(i)
          		amax=max(anti,antj)
			!...find separation of ith and jth points and skip if too large
          		sep1=sqrt((ui-uj)**2+(vi-vj)**2)
          		sep2=sqrt((ui+uj)**2+(vi+vj)**2)
          
          		if((sep1>amax).and.(sep2>amax)) cycle
          
			!...store useful numbers
          		nzmax=nzmax+1

		!...close loops over data points
		enddo
      	enddo
      	write(*,*) 'Total number of non-zero matrix elements = ',nzmax
      
end subroutine count_integrals

!-----*-----------------------------------------------------------------

subroutine makeintegrals(nv,nf,na,freqn,u,v,freq,antwdth,sigma,f_ra,f_dec,ra_ref,dec_ref,ra,diagflag,intPerBatch,intFileUnit,indx)
	implicit none

#ifdef MPI
  	include 'mpif.h'
  	integer mpi_status(MPI_STATUS_SIZE), errcode
#endif
      	integer  nv,nf,na,diagflag,intPerBatch,intFileUnit,indx
      	integer  freqn(nv)
      	double precision   u(nv),v(nv)
      	double precision   freq(nf),antwdth(nf),sigma(nf),f_ra(nf),f_dec(nf)
      	double precision   ra_ref,dec_ref
      	double precision   ra(na+1)

	! ... for independent bins (saves space)
	!      integer  maxbin1(nzmax),maxbin2(nzmax)
	!      real*8   integral1(nzmax,3),integral2(nzmax,3),temp(3)

	! ... for general bins
      	double precision   tempr(namax),tempi(namax)

      	integer  i,j,k,m,n,p,kmax
      	integer  ifr,jfr
      	double precision   ui,uj,vi,vj,si,sj,xi,xj,yi,yj,fi,fj,anti,antj,amax,sep1,sep2
      	double precision   xpos(nf),ypos(nf)
	logical done

	! calculate mask
	!open(unit=23,file='covar.out',status='replace')
      
      	!determine relative coordinates of pointing centres
  	do i=1,nf
    		call sky_coords(xpos(i), ypos(i), f_ra(i), f_dec(i), ra_ref, dec_ref)
	enddo
	
	done = .false.
      	nz = 0
	n = 0
	k = 0
	!...loop over data points
      	do j=1,nv
        	icn(j)=nz+1
        	jfr=freqn(j)
        	fj=freq(jfr)
        	antj=antwdth(jfr)
        	sj=sigma(jfr)
        	xj=xpos(jfr)
        	yj=ypos(jfr)
        	uj=u(j)
        	vj=v(j)	
        	do i=j,nv  ! only calculate elements in lower triangle
          		if((diagflag == 0).and.(i /= j)) cycle
          		ifr=freqn(i)
          		fi=freq(ifr)
          		anti=antwdth(ifr)
          		si=sigma(ifr)
          		xi=xpos(ifr)
          		yi=ypos(ifr)
          		ui=u(i)
          		vi=v(i)
          		amax=max(anti,antj)
			!...find separation of ith and jth points and skip if too large
          		sep1=dsqrt((ui-uj)**2+(vi-vj)**2)
          		sep2=dsqrt((ui+uj)**2+(vi+vj)**2)
          
          		if((sep1>amax).and.(sep2>amax)) cycle
			!...store useful numbers
          		nz=nz+1

			if( nz > ( n * mpi_nthreads + my_rank ) * intPerBatch .and. nz <= ( n * mpi_nthreads + my_rank + 1 ) * intPerBatch ) then
          			if(nz>nzmax) then
            				write(*,*)'ERROR: insufficient sparse storage at i,j = ',i,j
            				write(*,*)nz
#ifdef MPI
					call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
					stop
          			endif
				
				k = k + 1
				
				!...check if these integrals were calculated previously
				if( nz < indx .and. n /= ( indx - 1 ) / ( intPerBatch * mpi_nthreads ) ) then
					if( ( k == intPerBatch .or. nz == nzmax ) .and. .not.done ) then
						k = 0
						n = n + 1
					
						if( ( n * mpi_nthreads + my_rank ) * intPerBatch >= nzmax ) then
							if( my_rank > 0 ) then
								exit
							else
								done = .true.
							endif
						endif
					endif
					
					cycle
				endif
				
          			irn(k)=i
          			freqi(k)=fi
          			freqj(k)=fj
				
				!...calculate integrals

				!...for independent bins
				!          if(sep1>2d0*amax) then
				!            do m=1,3
				!              temp(m)=0d0
				!            enddo
				!          else
				!            call calcints(ui,vi,si,uj,vj,sj,na,ra,temp,kmax)
				!            maxbin1(nz)=kmax
				!          endif                
				!          do m=1,3
				!              integral1(nz,m)=temp(m)
				!          enddo
				!          if(sep2>2d0*amax) then
				!            do m=1,3
				!              temp(m)=0d0
				!            enddo
				!          else
				!            call calcints(ui,vi,si,-uj,-vj,sj,na,ra,temp,kmax)
				!            maxbin2(nz)=kmax
				!          endif                
				!          do m=1,3
				!              integral2(nz,m)=temp(m)
				!          enddo
	
				!...for general bins
	          		!if(sep1>2d0*amax) then
	            		!do m=1,na
	              			!temp(m)=0d0
	            		!enddo
	          		!else
				if( nz >= indx ) then
		          		if(.not.imagFlag) then
		          			!no mosaicking
		          			call calcints(ui,vi,si,uj,vj,sj,na,ra,tempr,kmax)                    	

		            			do m=1,na
							integral1r(k,m)=tempr(m)
		          			enddo
			   	 	else
		          			!mosaicking
		          			call calcints6(ui,vi,si,xi,yi,uj,vj,sj,xj,yj,na,ra,tempr,tempi,kmax)
						
		            			do m=1,na
							integral1r(k,m)=tempr(m)
							integral1i(k,m)=tempi(m)
		          			enddo
			    		endif
		          
		          		if(.not.imagFlag) then
		          			!no mosaicking
		            			call calcints(ui,vi,si,-uj,-vj,sj,na,ra,tempr,kmax)
		          
		          			do m=1,na
		              				integral2r(k,m)=tempr(m)
		          			enddo
			   	 	else
		          			!mosaicking						
		          			call calcints6(ui,vi,si,xi,yi,-uj,-vj,sj,xj,yj,na,ra,tempr,tempi,kmax)           
		          			
		            			do m=1,na
							integral2r(k,m)=tempr(m)
							integral2i(k,m)=tempi(m)
						
						enddo
			    		endif
				endif
				!calculated the batch then send the integrals to the root node
				if( ( k == intPerBatch .or. nz == nzmax ) .and. .not.done ) then
					if( my_rank /= 0 ) then
#ifdef MPI
						!send the integrals to the root node
						call MPI_SEND(irn(1:intPerBatch),intPerBatch,MPI_INTEGER,0,my_rank,MPI_COMM_WORLD,errcode)
						call MPI_SEND(freqi(1:intPerBatch),intPerBatch,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
						call MPI_SEND(freqj(1:intPerBatch),intPerBatch,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
						call MPI_SEND(integral1r(1:intPerBatch,1:na),intPerBatch*na,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
						call MPI_SEND(integral2r(1:intPerBatch,1:na),intPerBatch*na,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)

						if(imagFlag) then
							call MPI_SEND(integral1i(1:intPerBatch,1:na),intPerBatch*na,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
							call MPI_SEND(integral2i(1:intPerBatch,1:na),intPerBatch*na,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
						endif
#endif
					else
						do p = 0, mpi_nthreads - 1
							!receive the integrals
							if( p > 0 ) then
#ifdef MPI
								call MPI_RECV(irn(1:intPerBatch),intPerBatch,MPI_INTEGER,p,p,MPI_COMM_WORLD,mpi_status,errcode)
								call MPI_RECV(freqi(1:intPerBatch),intPerBatch,MPI_DOUBLE_PRECISION,p,p,MPI_COMM_WORLD,mpi_status,errcode)
								call MPI_RECV(freqj(1:intPerBatch),intPerBatch,MPI_DOUBLE_PRECISION,p,p,MPI_COMM_WORLD,mpi_status,errcode)
								call MPI_RECV(integral1r(1:intPerBatch,1:na),intPerBatch*na,MPI_DOUBLE_PRECISION,p,p,MPI_COMM_WORLD,mpi_status,errcode)
								call MPI_RECV(integral2r(1:intPerBatch,1:na),intPerBatch*na,MPI_DOUBLE_PRECISION,p,p,MPI_COMM_WORLD,mpi_status,errcode)
			
								if(imagFlag) then
									call MPI_RECV(integral1i(1:intPerBatch,1:na),intPerBatch*na,MPI_DOUBLE_PRECISION,p,p,MPI_COMM_WORLD,mpi_status,errcode)
									call MPI_RECV(integral2i(1:intPerBatch,1:na),intPerBatch*na,MPI_DOUBLE_PRECISION,p,p,MPI_COMM_WORLD,mpi_status,errcode)
								endif
#endif
							endif
							
							!write integrals
							do m = 1, intPerBatch
								if( ( n * mpi_nthreads + p ) * intPerBatch + m > nzmax ) exit
								
								if( ( n * mpi_nthreads + p ) * intPerBatch + m >= indx ) then
									write(intFileUnit) ( n * mpi_nthreads + p ) * intPerBatch + m
									write(intFileUnit) irn(m), freqi(m), freqj(m)
							
									write(intFileUnit) integral1r(m,1:na)
									write(intFileUnit) integral2r(m,1:na)
									if(imagFlag) then
										write(intFileUnit) integral1i(m,1:na)
										write(intFileUnit) integral2i(m,1:na)
									endif
								endif
							enddo
							if( ( n * mpi_nthreads + p + 1 ) * intPerBatch  > nzmax ) exit
						enddo
					endif
					k = 0
					n = n + 1
#ifdef MPI
					call MPI_BARRIER(MPI_COMM_WORLD,errcode)
#endif
					
					if( ( n * mpi_nthreads + my_rank ) * intPerBatch >= nzmax ) then
						if( my_rank > 0 ) then
							exit
						else
							done = .true.
						endif
					endif
				endif
				
			endif

		!...close loops over data points
 	    	enddo
		if (nv < 100) then
		       write(*,*) 'Fraction of table completed = ',float(j)/float(nv)
		else
			if( my_rank == 0 .and. mod( j, nv / 100 ) == 0 ) write(*,*) 'Fraction of table completed = ',float(j)/float(nv)
		endif
      	enddo
	
      	icn(nv+1)=nz+1
      	if( my_rank == 0 ) write(*,*) 'Total number of non-zero matrix elements = ',nz
      
end subroutine makeintegrals

!-----*-----------------------------------------------------------------

subroutine combineIntegrals(na,intPerNode)
	implicit none

#ifdef MPI
  	include 'mpif.h'
  	integer mpi_status(MPI_STATUS_SIZE), errcode
#endif
      	integer  na, intPerNode
	integer i,j
	
			
	if( my_rank /= 0 ) then
		!send the integrals to the root node
#ifdef MPI
		call MPI_SEND(irn(1:intPerNode),intPerNode,MPI_INTEGER,0,my_rank,MPI_COMM_WORLD,errcode)
		call MPI_SEND(freqi(1:intPerNode),intPerNode,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
		call MPI_SEND(freqj(1:intPerNode),intPerNode,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
#endif		
		if(.not.imagFlag) then
#ifdef MPI
	          	call MPI_SEND(integral1r(1:intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
			call MPI_SEND(integral2r(1:intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
#endif
		else
#ifdef MPI
	          	call MPI_SEND(integral1r(1:intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
	          	call MPI_SEND(integral1i(1:intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
	          	call MPI_SEND(integral2r(1:intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
	          	call MPI_SEND(integral2i(1:intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,errcode)
#endif
		endif
	else
		do i = 1, mpi_nthreads - 1
			j = i * intPerNode
			
			!receive the integrals
#ifdef MPI
			call MPI_RECV(irn(j+1:j+intPerNode),intPerNode,MPI_INTEGER,i,i,MPI_COMM_WORLD,mpi_status,errcode)
			call MPI_RECV(freqi(j+1:j+intPerNode),intPerNode,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,mpi_status,errcode)
			call MPI_RECV(freqj(j+1:j+intPerNode),intPerNode,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,mpi_status,errcode)
#endif
			
			if(.not.imagFlag) then
#ifdef MPI
				call MPI_RECV(integral1r(j+1:j+intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,mpi_status,errcode)
				call MPI_RECV(integral2r(j+1:j+intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,mpi_status,errcode)
#endif
			else
#ifdef MPI
				call MPI_RECV(integral1r(j+1:j+intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,mpi_status,errcode)
				call MPI_RECV(integral2r(j+1:j+intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,mpi_status,errcode)
				call MPI_RECV(integral1i(j+1:j+intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,mpi_status,errcode)
				call MPI_RECV(integral2i(j+1:j+intPerNode,1:na),intPerNode*na,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,mpi_status,errcode)
#endif
			endif
		enddo
	endif
      
end subroutine combineIntegrals

!-----*-----------------------------------------------------------------

subroutine makecm(nv,na,rms,cps1d,confps1d,gps1d,sp,convf,freq0,intFileUnit)
      	implicit none

#ifdef MPI
  	include 'mpif.h'
  	integer errcode
#endif
      	integer  nv,na,intFileUnit
      	double precision   rms(nv)
      	double precision   cps1d(na),confps1d(na),gps1d(na)
      	double precision   sp,convf,freq0

	! ... for independent bins (saves space)
	!      integer  maxbin1(nzmax),maxbin2(nzmax)
	!      real*8   integral1(nzmax,3),integral2(nzmax,3)

	! ... for general bins

      	integer  i,j,m,p,k
      	double precision   cn,gn,fi,fj,value1,value2,pstot,tcmb
	double precision int1r(na),int2r(na),int1i(na),int2i(na)

	! ... for independent bins
	!      common /maxbincom/   maxbin1,maxbin2

	! initialise variables

      	tcmb=2.726d0

	! calculate covariance matrix
      	do j=1,nv
        	do p=icn(j),icn(j+1)-1
		
			!read relevant information from the integrals file
			read(intFileUnit) k
			if( k /= p ) then
				write(*,*)"Problem in reading the integrals file. Aborting"
#ifdef MPI
				call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
				stop
			endif
			read(intFileUnit) i, fi, fj
			read(intFileUnit) int1r(1:na)
			read(intFileUnit) int2r(1:na)
			if( imagFlag ) then
				read(intFileUnit) int1i(1:na)
				read(intFileUnit) int2i(1:na)
			endif
			irn(p) = i
		
			!...access relevant information from arrays
          		!i=irn(p)
          		!fi=freqi(p)
          		!fj=freqj(p)
			
			!...for independent bins..
			!          kmax1=maxbin1(p)
			!          kmax2=maxbin2(p)
			!.........................
          		cn=(cmbfcn(tcmb,fi)/cmbfcn(tcmb,freq0))*(cmbfcn(tcmb,fj)/cmbfcn(tcmb,freq0))
          		gn=((fi/freq0)**sp)*((fj/freq0)**sp)
			!...calculate the corresponding element of covariance matrix

			!...for independent bins
			!          value1=0d0
			!          do m=1,3
			!            nbin=kmax1-2+m
			!            if((nbin.ge.1).and.(nbin.le.na)) then
			!              pstot=(cn*cps1d(nbin)+gn*gps1d(nbin))*convf
			!              value1=value1+pstot*integral1(p,m) 
			!            endif
			!          enddo
			!          value2=0d0
			!          do m=1,3
			!            nbin=kmax2-2+m
			!            if((nbin.ge.1).and.(nbin.le.na)) then
			!              pstot=(cn*cps1d(nbin)+gn*gps1d(nbin))*convf
			!              value2=value2+pstot*integral2(p,m) 
			!            endif
			!          enddo

			!...for general bins
          		value1=0d0
          		do m=1,na
            			pstot=(cn*cps1d(m)+gn*gps1d(m)+confps1d(m))*convf
            			value1=value1+pstot*int1r(m)
          		enddo
          		value2=0d0
          		do m=1,na
            			pstot=(cn*cps1d(m)+gn*gps1d(m)+confps1d(m))*convf
            			value2=value2+pstot*int2r(m)
          		enddo

          		cm_rr(p)=0.5*(value1+value2)
          		cm_ii(p)=0.5*(value1-value2)
          
          		if(imagFlag) then
          			value1=0d0
          			do m=1,na
            				pstot=(cn*cps1d(m)+gn*gps1d(m)+confps1d(m))*convf
            				value1=value1+pstot*int1i(m)
          			enddo
          			value2=0d0
          			do m=1,na
            				pstot=(cn*cps1d(m)+gn*gps1d(m)+confps1d(m))*convf
            				value2=value2+pstot*int2i(m)
          			enddo

          			cm_ir(p)=0.5*(value2+value1)
          			cm_ri(p)=0.5*(value2-value1)
			else
          			cm_ir(p)=0d0
          			cm_ri(p)=0d0
          		endif
          
			!...add (diagonal) instrumental noise variance
          		if(i==j) then
            			cm_rr(p)=cm_rr(p)+rms(i)**2
            			cm_ii(p)=cm_ii(p)+rms(i)**2
          		endif
        	enddo
      	enddo
	
end subroutine makecm

!-----*-----------------------------------------------------------------

logical function getLastIntegralIndex(na, intFileUnit, indx)
      	implicit none
      	integer  na, intFileUnit, indx
      	integer  i, j, k, iostatus
      	double precision   d1
	double precision int(na)
	
	j = 0
	k = 0
	do
		getLastIntegralIndex = .true.
		
        	!read relevant information from the integrals file
		read(intFileUnit, IOSTAT = iostatus) i
		if( iostatus /= 0 ) exit
		k = i
		if( k /= j + 1 ) then
			backspace(intFileUnit)
			k = j
			getLastIntegralIndex = .false.
			exit
		endif
		j = k
		read(intFileUnit, IOSTAT = iostatus) i, d1, d1
		if( iostatus /= 0 ) exit
		read(intFileUnit, IOSTAT = iostatus) int(1:na)
		if( iostatus /= 0 ) exit
		read(intFileUnit, IOSTAT = iostatus) int(1:na)
		if( iostatus /= 0 ) exit
		if( imagFlag ) then
			read(intFileUnit, IOSTAT = iostatus) int(1:na)
			if( iostatus /= 0 ) exit
			read(intFileUnit, IOSTAT = iostatus) int(1:na)
			if( iostatus /= 0 ) exit
		endif
		
		getLastIntegralIndex = .false.
	enddo
	
	if( getLastIntegralIndex ) k = k - 1
	indx = k
	
end function getLastIntegralIndex

!-----*-----------------------------------------------------------------

subroutine moveIntegralFilePointer(na, intFileUnit, indx)
      	implicit none
      	integer  na, intFileUnit, indx
      	integer  i, j
      	double precision   d1
	double precision int(na)
	
	do i = 1, indx
        	!read relevant information from the integrals file
		read(intFileUnit) j
		read(intFileUnit) j, d1, d1
		read(intFileUnit) int(1:na)
		read(intFileUnit) int(1:na)
		if( imagFlag ) then
			read(intFileUnit) int(1:na)
			read(intFileUnit) int(1:na)
		endif
	enddo
	
end subroutine moveIntegralFilePointer

!-----*-----------------------------------------------------------------

subroutine spnorm(n,nz,irn,icn,sa,norm,samin,samax)
      	implicit none
      	integer  n,nz,irn(nz),icn(n+1)
      	double precision   sa(nz),norm,samin,samax

      	integer  i,j,p
      	double precision   trace

      	samin=sa(1)
      	samax=sa(1)
      	trace=0d0
      	do j=1,n
        	do p=icn(j),icn(j+1)-1
          		i=irn(p)
          		if(sa(p)<samin) samin=sa(p)
          		if(sa(p)>samax) samax=sa(p)
          		if(i.eq.j) trace=trace+sa(p)
        	enddo
      	enddo
      	norm=trace/dble(n)
      	if(norm.eq.0d0) return
      	do p=1,nz
        	sa(p)=sa(p)/norm
      	enddo
      
end subroutine spnorm

!-----*-----------------------------------------------------------------

! intensity frequency dependence due to CMBR
double precision function cmbfcn(tcmb,freq)
      	implicit none
      	double precision tcmb,freq

      	double precision hk,x
      	parameter (hk=0.0479927)        ! this is h/k for frequency in GHz

      	x=hk*freq/tcmb
      	cmbfcn=(x**4)*exp(x)/((exp(x)-1)**2)
      
end function cmbfcn

!-----*-----------------------------------------------------------------

! calculates conversion factor between equivalent CMBR thermodynamic
! temperature fluctuation in K to intensity fluctuations in W/m^2/Hz/sr
double precision function dt_to_di(tcmb,freq)
      	implicit none
      	double precision tcmb,freq

      	double precision norm
      	parameter (norm=6.675d-20)      ! k**3/(h**2 c**2) in SI units
      
      	dt_to_di=2d0*(tcmb**2)*norm*cmbfcn(tcmb,freq)
      
end function dt_to_di

!-----*-----------------------------------------------------------------

subroutine calcints(ui,vi,si,uj,vj,sj,na,ra,temp,kmax)
      	implicit none

#ifdef MPI
  	include 'mpif.h'
  	integer errcode
#endif
      	integer  na,kmax
      	double precision   ui,vi,si,uj,vj,sj
      	double precision   ra(na+1)
	!...for independent bins
	!      real*8   temp(3)
	!...for general bins
      	double precision   temp(namax)

      	integer  k,ifail,jiter
      	double precision   PI,rlo,rhi,value,eps
      	double precision   ai,aj

      	PI=4d0*atan(1d0)        ! this really does equal PI
      	eps=1.d-2                ! fractional accuracy in integration

      	ai=2d0*(PI**2)*(si**2)
      	aj=2d0*(PI**2)*(sj**2)

      	ar1=ai*(ui**2+vi**2)+aj*(uj**2+vj**2)
      	ar2=ai+aj
      	ar3=sqrt((ai*ui+aj*uj)**2+(ai*vi+aj*vj)**2)
      	b=(2d0*PI*si*sj)**2    
      
      	call fcnmax(na,ra,kmax,fmax)

	!... for general bins
      	do k=1,na
        	rlo=ra(k)
        	rhi=ra(k+1)
        	call qtrap(normfcn,rlo,rhi,value,eps,ifail,jiter)
        	if(ifail==1) then
        		write(*,*)"qtrap failed in calcints"
#ifdef MPI
			call MPI_ABORT(MPI_COMM_WORLD,errcode)
#endif
			stop
	  	endif
        	temp(k)=fmax*value
      	enddo

end subroutine calcints

!-----*-----------------------------------------------------------------

subroutine calcints6(ui,vi,si,xi,yi,uj,vj,sj,xj,yj,na,ra,tempr,tempi,kmax)
      	implicit none
      	integer  na,kmax
      	double precision   ui,vi,si,xi,yi,uj,vj,sj,xj,yj
      	double precision   ra(na+1)
	!...for independent bins
	!      real*8   temp(3)
	!...for general bins
      	double precision   tempr(namax),tempi(namax)

      	integer  k,ifail,jiter
      	double precision   PI,rlo,rhi,eps
      	double precision   ai,aj,si2,sj2,PI2,x1,x2,w1,w2,c,d
      	double complex c1,c2

      	PI=4d0*atan(1d0)        ! this really does equal PI
      	eps=1.d-2                ! fractional accuracy in integration
	
      	si2=si**2
      	sj2=sj**2
      	PI2=PI**2
      
      	ai=2d0*PI2*si2
      	aj=2d0*PI2*sj2

      	ar1=ai*(ui**2+vi**2)+aj*(uj**2+vj**2)
      	ar2=ai+aj
      
  	x1=xj-xi
  	x2=yj-yi
      	w1=si2*ui+sj2*uj
  	w2=si2*vi+sj2*vj 

  	c=4d0*PI2*(w1*w1+w2*w2)-(x1*x1+x2*x2)
  	d=4d0*PI*(w1*x1+w2*x2)
      
      	c1=dcmplx(c,d)
      	c2=sqrt(c1)
      
      	ar3_real=dble(c2)*2.*PI
      	ar3_img=dimag(c2)*2.*PI
      
      	b=4d0*PI2*si2*sj2
      
      	call cfcnmax(na,ra,kmax,fmax)

	!... for general bins
      	do k=1,na
        	rlo=ra(k)
        	rhi=ra(k+1)
        
        	!numerical radial integration
        	call cqtrap(normcfcn,rlo,rhi,tempr(k),tempi(k),eps,ifail,jiter)
        	!call qsimp(normcfcn_real,rlo,rhi,tempr(k),eps,ifail,jiter)
		!if(ifail==1) then
		!	write(*,*)"cqsimp failed in calcints6",rlo,rhi,tempr(k),tempi(k),eps,ifail,jiter
		!	stop
		!endif
        	tempr(k)=tempr(k)*fmax
        	tempi(k)=tempi(k)*fmax
      	enddo

end subroutine calcints6

!-----*-----------------------------------------------------------------

! calculates the maximium value of fcn and annulus in which it occurs
subroutine fcnmax(na,ra,maxbin,fmax)
      	implicit none
      	integer  na,maxbin
      	double precision   ra(na+1),fmax
      	integer  npts,i,k
      	parameter (npts=10)
      	double precision   x,xstep,xlo,xhi,f

      	fmax=fcn(ra(1))

      	do k=1,na
        	xlo=ra(k)
        	xhi=ra(k+1)
        	xstep=(xhi-xlo)/dble(npts-1)
        	do i=1,npts
          		x=xlo+dble(i-1)*xstep
          		f=fcn(x)

          		if(f>fmax) then
            			fmax=f
            			maxbin=k
          		endif
        	enddo
      	enddo
      
end subroutine fcnmax

!-----*-----------------------------------------------------------------

! calculates the maximium value of the module of complex function cfcn and annulus in which it occurs
subroutine cfcnmax(na,ra,maxbin,fmax)
      	implicit none
      	integer  na,maxbin
      	double precision   ra(na+1),fmax
      	integer  npts,i,k
      	parameter (npts=10)
      	double precision   x,xstep,xlo,xhi,f

      	fmax=cdabs(cfcn(ra(1)))

	do k=1,na
        	xlo=ra(k)
        	xhi=ra(k+1)
        	xstep=(xhi-xlo)/dble(npts-1)
        	do i=1,npts
          		x=xlo+dble(i-1)*xstep
          		f=cdabs(cfcn(x))

          		if(f>fmax) then
            			fmax=f
            			maxbin=k
          		endif
        	enddo
      	enddo
      
end subroutine cfcnmax

!-----*-----------------------------------------------------------------

! calculates fcn

double precision function fcn(x)
      	implicit none
      	double precision x
      	double precision arg1,arg2,arg3,PI,logf

      	PI=4d0*datan(1d0)        ! this really does equal PI
      	arg1=ar1
      	arg2=ar2*x**2
      	arg3=2*ar3*x
      	logf=dlog(2d0*PI*b)-arg1-arg2+logbessi0(arg3)-dlog(x)
      	fcn=dexp(logf)
      
end function fcn

!-----*-----------------------------------------------------------------

! compute the real and imaginary parts of the angular integral as a 
! a function of radius 
double complex function cfcn(rho)
	implicit none
      
      	double precision rho
	double precision rintp, iintp
  	double precision PI, c, d, arg
      
      	PI=4d0*datan(1d0)        ! this really does equal PI

  	arg = -ar1 - ar2 * rho * rho
  	c = ar3_real * rho
      	d = ar3_img * rho
  	call CDEFACBESSI0(arg, c, d, rintp, iintp)

  	cfcn = dcmplx(b * 2. * PI * rintp / rho, b * 2. * PI * iintp / rho)

end function cfcn

!-----*-----------------------------------------------------------------

! compute the normalized real and imaginary parts of the angular integral as a 
! a function of radius 
double complex function normcfcn(rho)
	implicit none
      
      	double precision rho
	double precision rintp, iintp
  	double precision PI, c, d, arg
      
      	PI=4d0*atan(1d0)        ! this really does equal PI

  	arg = -ar1 - ar2 * rho * rho - log(fmax)
  	c = ar3_real * rho
      	d = ar3_img * rho
  	call CDEFACBESSI0(arg, c, d, rintp, iintp)
      	rintp = b * 2. * PI * rintp / rho
      	iintp = b * 2. * PI * iintp / rho

  	normcfcn = dcmplx(rintp, iintp)

end function normcfcn

!-----*-----------------------------------------------------------------

! compute the real part of the angular integral as a 
! a function of radius
double precision function normcfcn_real(rho)
	implicit none
      
      	double precision rho
	double precision rintp, iintp
  	double precision PI, c, d, arg
      
      	PI=4d0*atan(1d0)        ! this really does equal PI

  	arg = -ar1 - ar2 * rho * rho - log(fmax)
  	c = ar3_real * rho
      	d = ar3_img * rho
  	call CDEFACBESSI0(arg, c, d, rintp, iintp)

  	normcfcn_real = b * 2. * PI * rintp / rho

end function normcfcn_real

!-----*-----------------------------------------------------------------

! compute the imaginary part of the angular integral as a 
! a function of radius 
double precision function normcfcn_img(rho)
	implicit none
      
      	double precision rho
	double precision rintp, iintp
  	double precision PI, c, d, arg
      
      	PI=4d0*datan(1d0)        ! this really does equal PI

  	arg = -ar1 - ar2 * rho * rho - log(fmax)
  	c = ar3_real * rho
      	d = ar3_img * rho
  	call CDEFACBESSI0(arg, c, d, rintp, iintp)

  	normcfcn_img = b * 2. * PI * iintp / rho

end function normcfcn_img

!-----*-----------------------------------------------------------------

! calculates fcn divided by fmax
double precision function normfcn(x)
      	implicit none
      	double precision x
	double precision arg1,arg2,arg3,PI,logf

      	PI=4d0*datan(1d0)        ! this really does equal PI
      	arg1=ar1
      	arg2=ar2*x**2
      	arg3=2*ar3*x
      	logf=dlog(2d0*PI*b)-arg1-arg2+logbessi0(arg3)-dlog(x)-dlog(fmax)
      	normfcn=dexp(logf)
      
end function normfcn

!-----*-----------------------------------------------------------------

double precision FUNCTION logbessi0(x)
      	double precision x,ax,bessi0
      	DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      	SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      	DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0, &
      	1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      	DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1, &
      	0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1, &
      	-0.1647633d-1,0.392377d-2/
      	if(dabs(x)<3.75d0) then
       	 	y=(x/3.75d0)**2
        	bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
        	logbessi0=dlog(bessi0)
      	else
        	ax=dabs(x)
        	y=3.75d0/ax
        	logbessi0=ax-0.5d0*dlog(ax)+dlog(q1+y*(q2+y*(q3+y*(q4+y*(q5+y* &
        	(q6+y*(q7+y*(q8+y*q9))))))))
      	endif
      
END FUNCTION logbessi0

!-----*-----------------------------------------------------------------

double precision function dmin(x,y)
      	implicit none
      	double precision x,y

      	if(x<y) then
        	dmin=x
	elseif(y.le.x) then
        	dmin=y
      	endif
      
end function dmin

!-----*-----------------------------------------------------------------

double precision function dmax(x,y)
      	implicit none
      	double precision x,y

      	if(x>y) then
        	dmax=x
      	elseif(y.ge.x) then
        	dmax=y
      	endif
      
end function dmax

!-----*-----------------------------------------------------------------

end module integrals
