module fits_handling

        use params
        use utilities

contains
!=======================================================================
! Open fits uv file and read header

        subroutine phitsread_open(filename, CRVAL4, CRVAL5, CRVAL6, OBSRA, OBSDEC, TELESCOP, PSCAL1, PSCAL2)

        implicit none

        character(len=100) filename
        integer iunit,status
        integer rwstat
        character(len=80) ftcom
        double precision OBSRA,OBSDEC,CRVAL4,CRVAL5,CRVAL6,PSCAL1,PSCAL2 !,PSCAL4,PZERO4,PSCAL5,PZERO5
        character(len=100) TELESCOP !,OBJECT,INSTRUME

! Open the existing FITS file with readonly access:

        iunit=15
        rwstat=0
        status=0
        call ftopen(iunit,filename,rwstat,bksize,status)
!         write(*,*) iunit,filename,rwstat,bksize
! 	  write(*,*) 'status = ',status
	if( status /= 0 ) then
        if(is_root) write(*,*) "ERROR: could not open the file "//trim(filename)
        if(is_root) write(*,*) "Aborting"
		call halt_program
	endif

!       read the required primary array keywords
        call ftgprh(iunit,simple,bitpix,naxis,naxes,pcount,gcount,extend,status)
!        write(*,*)iunit,simple,bitpix,naxis,naxes,pcount,gcount,extend,status

! Bunch seemed to have difficulty reading the naxes array - here it is
! hard wired for safety:

        naxes(1) = 0                   ! No image: uv data                              
        naxes(2) = 3                   ! Complex: cos, sin, weight                      
        naxes(3) = 1                   ! Number of polarisations                      
        naxes(4) = 1                   ! Number of frequencies                      
        naxes(5) = 1                   ! RA                      
        naxes(6) = 1                   ! Dec                       
        
	  if(status.eq.104) call halt_program('File does not exist')
        
	  !call ftgkys(iunit,'OBJECT',OBJECT,ftcom,status)
!         write(*,*) OBJECT,status
	  status = 0
	  call ftgkys(iunit,'TELESCOP',TELESCOP,ftcom,status)
!         write(*,*) TELESCOP,status
        !status = 0
 	  !call ftgkys(iunit,'INSTRUME',INSTRUME,ftcom,status)
!         write(*,*) INSTRUME,status
        status = 0
        
	  call ftgkyd(iunit,'OBSRA',OBSRA,ftcom,status)
!         write(*,*) OBSRA,status
        status = 0
        call ftgkyd(iunit,'OBSDEC',OBSDEC,ftcom,status)
!         write(*,*) OBSDEC,status
        status = 0
        call ftgkyd(iunit,'CRVAL4',CRVAL4,ftcom,status)
!         write(*,*) CRVAL4,status
        status = 0
        call ftgkyd(iunit,'CRVAL5',CRVAL5,ftcom,status)
!         write(*,*) CRVAL5,status
        status = 0
        call ftgkyd(iunit,'CRVAL6',CRVAL6,ftcom,status)
!         write(*,*) CRVAL6,status
        status = 0
        call ftgkyd(iunit,'PSCAL1',PSCAL1,ftcom,status)
!         write(*,*) PSCAL1,status
        status = 0
        call ftgkyd(iunit,'PSCAL2',PSCAL2,ftcom,status)
!         write(*,*) PSCAL2,status
        !status = 0
        !call ftgkyd(iunit,'PSCAL5',PSCAL5,ftcom,status)
!         write(*,*) PSCAL5,status
        !status = 0
        !call ftgkyd(iunit,'PZERO5',PZERO5,ftcom,status)
!         write(*,*) PZERO5,status
        !status = 0
	  
!        print*,'Telescope: ',telescop(1:3)
!        print*,'Frequency: ',crval4
!        print*,'Pointing centre: ',crval5,crval6
!        print*,'Phase centre: ',obsra,obsdec
!	  pause
      
        return
	  end subroutine phitsread_open
	  
!------------------------------------------------------------------ 
! read one group of specified uvfile, merge over frequency and return
! visibilities, weights, and uv coordinates 

      subroutine phitsread_one(group,vis_buffer,uv_buffer,weight_buffer, &
      rms_buffer,baseline, CRVAL4, PSCAL1, PSCAL2)

      implicit none
        
      double complex vis_buffer
      double precision weight_buffer, uv_buffer(2)
      double precision rms_buffer
      integer baseline
      double precision freq,uscale,vscale
      integer iunit,status
      integer group,i,j
      integer fpixels(7),lpixels(7),inc(7)
      integer axes(6),npoints
      double precision ivalue(105),parms(6)
      double precision weight
      logical anyflg
      double precision CRVAL4, PSCAL1, PSCAL2

      iunit = 15
	status = 0

      fpixels(1) = 0
      lpixels(1) = 0
      do i = 2,naxis
        fpixels(i-1) = 1
        axes(i-1) = naxes(i)
        lpixels(i-1) = naxes(i)
        inc(i-1) = 1
      end do
      npoints = 1
      do i = 2,naxis
        npoints = npoints * naxes(i)
      end do

	freq = CRVAL4
	uscale = PSCAL1
	vscale = PSCAL2
	
      call ftggpd(iunit,group,1,naxis,parms,status)
      call ftgsvd(iunit,group,naxis-1,axes,fpixels,lpixels,inc,0,ivalue,anyflg,status)
      
      uv_buffer(1) = parms(1)*uscale*freq
      uv_buffer(2) = parms(2)*vscale*freq
      baseline = int(parms(6))
      vis_buffer = (0.0,0.0)
      rms_buffer = 0.0
      weight_buffer = 0.0
      
	do j = 1,npoints, 3
        weight = ivalue(j+2)
        if(weight.lt.0.0) weight = 0.0
!       Calculate rms of *complex* number
        rms_buffer=rms_buffer+(ivalue(j)**2 + ivalue(j+1)**2+ &
        (2*ivalue(j)*ivalue(j+1)))*weight
        weight_buffer = weight_buffer + weight
        vis_buffer = vis_buffer+dcmplx(ivalue(j),ivalue(j+1))*weight
      end do
	
      if(weight_buffer.gt.0) then
        vis_buffer=vis_buffer/weight_buffer
        rms_buffer=rms_buffer/weight_buffer
        rms_buffer=sqrt(rms_buffer/(float(npoints)))
      else
        vis_buffer = (0.0,0.0)
        rms_buffer = 0.0
        weight_buffer = 0.0
      end if
        
      return
	end subroutine phitsread_one

!=======================================================================
! Close fits uv file and read header

      subroutine phitsread_close

      implicit none

      integer iunit,status

! Close the existing FITS file:

 	iunit = 15
	status = 0
	call ftclos(iunit,status)
      if(status.gt.0) then
          if(is_root) write(*,*) 'Error closing vis data file - status = ',status
	  call halt_program('    ...stopping')
	endif
        
      return
	end subroutine phitsread_close

!=======================================================================
end module fits_handling
