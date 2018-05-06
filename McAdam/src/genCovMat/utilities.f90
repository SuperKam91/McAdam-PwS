module utilities
	
contains
	
!-----*-----------------------------------------------------------------
	
!convert RA in (h,m,s) to radians
subroutine sla_ctf2r(hour, min, sec, rad, status)
	implicit none
	
	!input variables
	double precision hour, min, sec
	integer status
	!output variables
	double precision rad !RA in radians
	!work variables
	integer ihour, imin
	double precision d2s
	parameter(d2s = 86400.d0)
	double precision PI
	
	PI = 4.d0 * atan(1.d0)
	
	ihour = int(hour)
	imin = int(min)
	
	status = 0;
	if((sec < 0.d0) .or. (sec >= 60.d0)) status = 3
	if((imin < 0) .or. (imin > 59)) status = 2
	if((ihour < 0) .or. (ihour > 23)) status = 1
	rad = 2.d0 * PI * (60.d0 * (60.d0 * ihour + imin) + sec) / d2s
	
end subroutine sla_ctf2r
	
!-----*-----------------------------------------------------------------
	
!convert DEC in (d,m,s) to radians
subroutine sla_daf2r(deg, amin, asec, rad, status)
	implicit none
	
	!input variables
	double precision deg, amin, asec
	integer status
	!output variables
	double precision rad !DEC in radians
	!work variables
	integer ideg, iamin
	double precision as2r
	parameter(as2r = 4.84813681109535994e-6)
	
	ideg  = int(deg)
	if( ideg < 0d0 ) ideg = 360d0+ideg
	iamin = int(amin)
	
	status = 0
	if((asec < 0.d0) .or. (asec >= 60.d0)) status = 3
	if((iamin < 0) .or. (iamin > 59)) status = 2
	if((ideg < 0) .or. (deg > 360)) status = 1
	rad = as2r * (60.d0 * (60.d0 * ideg + iamin) + asec)
	
end subroutine sla_daf2r
	
!-----*-----------------------------------------------------------------
	
! return the projection of the point given by ra, dec onto the plane 
! tangent to (ra_ref, dec_ref).
subroutine sky_coords(x, y, ra, dec, ra_ref, dec_ref)
	implicit none
	
	double precision x,y,ra,dec,ra_ref,dec_ref
	double precision alpha
	
	alpha = ra - ra_ref
	x = sin(alpha) * cos(dec)
	y = cos(dec_ref) * sin(dec) - sin(dec_ref) * cos(dec) * cos(alpha)
	
end subroutine sky_coords
	
!-----*-----------------------------------------------------------------
subroutine halt_program(message)
        implicit none

#ifdef MPI
        include 'mpif.h'
#endif
        character(LEN=*), intent(in), optional :: message

#ifdef MPI
        integer :: errorcode=1
        integer :: mpierror
#endif

        if (present(message)) then
            write(*,'( 20("=") )')
            write(*,'(A)') trim(adjustl(message))
            write(*,'( 20("=") )')
        end if

#ifdef MPI
        call MPI_ABORT(MPI_COMM_WORLD,errorcode,mpierror)
#else
        stop 1
#endif

end subroutine halt_program

!======================================
	
end module utilities
	
