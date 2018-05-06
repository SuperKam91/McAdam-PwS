program writechol

	integer nv
      	parameter(nv=4283)
      	integer i,j,iostatus
      	double precision covar(2*nv,2*nv)
      
      	open(unit=12,file='covar',status='old')
     	open(unit=13,file='../../../data/AMI/fake/jon221107/AJ/AJ3.lcm', &
      	form='unformatted',status='unknown')
      	write(13)2*nv
      	covar=0.
      	do
		read(12,*,IOSTAT=iostatus)i,j,covar(i,j)
            	!end of file?
        	if(iostatus<0) exit
	enddo
      	close(12)
      
      	do i=1,2*nv
        	write(13)covar(i,1:i)
	enddo
      	close(13)
      
end program writechol
