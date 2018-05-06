module pws_wrapper


contains

!=======================================================================


    subroutine initializePwS(parfile)

    use, intrinsic :: iso_c_binding
    implicit none
    
    interface
      integer(c_int) &
      function AMI_initializePwSlike( string ) bind (C, name = "AMI_initializePwSlike" )
        use iso_c_binding
        character (kind = C_CHAR ) :: string(*)
      end function AMI_initializePwSlike
    end interface

    character*100 parfile
    integer retval

    retval = AMI_initializePwSlike( trim(parfile) // C_NULL_CHAR )
    if (retval.ne.0) then
      write(*,*) 'Error: AMI_initializePwSlike returned', retval
      stop
    endif

    end subroutine initializePwS

!=======================================================================

    subroutine get_pws_patch(patchnum)

    use, intrinsic :: iso_c_binding
    implicit none

    interface
      integer(c_int) &
      function AMI_fetchPatch( patchnum ) bind (C, name = "AMI_fetchPatch" )
        use iso_c_binding
        integer (kind = c_int) :: patchnum
      end function AMI_fetchPatch
    end interface

    integer, intent(in) :: patchnum
    integer retval

    retval = AMI_fetchPatch(patchnum)
    if (retval.ne.0) then
      write(*,*) 'Error fetching PwS patch number', patchnum
      call release_pws_patch()
      stop
    endif

    end subroutine get_pws_patch

!=======================================================================

    subroutine release_pws_patch()

    use, intrinsic :: iso_c_binding
    implicit none

    interface
      subroutine AMI_ReleasePatch() bind (C, name = "AMI_ReleasePatch" )
        use iso_c_binding
      end subroutine AMI_ReleasePatch
    end interface

    call AMI_ReleasePatch()

    end subroutine release_pws_patch

!=======================================================================

    subroutine PlanckLhood(Cube, flag, Lhood, Cube_out)

    use, intrinsic :: iso_c_binding
    implicit none

    interface
      integer(c_int) &
      function AMI_PwSLikelihood( LikeValue, NParams, In_params, Out_params, NoDuocimation ) bind (C, name = "AMI_PwSLikelihood" )
        use iso_c_binding
        integer(kind = c_int), intent(in) :: NParams
        !double precision, intent(in)      :: In_params(NParams + 1) !for current implementation, always set hyperparameter for Planck a value
	double precision, intent(in)      :: In_params(NParams)
        !double precision, intent(out)     :: Out_params(Nparams + 1) !for current implementation, always set hyperparameter for Planck a value
	double precision, intent(out)     :: Out_params(Nparams)
        double precision, intent(out)     :: LikeValue
        integer(kind = c_int), intent(in) :: NoDuocimation
      end function AMI_PwSLikelihood
    end interface

    integer, parameter            :: NParams = 6 !for current implementation, don't count hyperparam as one of parameters in PwS
    integer, intent(out)          :: flag
    integer                       :: i
    double precision, intent(in)  :: Cube(NParams) !think these two need to be changed for hyperparameter implementation
    double precision, intent(out) :: Cube_out(NParams) !think this was what caused previous bug
    double precision, intent(out) :: Lhood

    integer, parameter                :: NoDuocimation = 1

    flag = AMI_PwSLikelihood( Lhood, NParams, Cube, Cube_out, NoDuocimation) !current implementation of PwS doesn't need NParams to be incremented inside it to consider hyperparameter
    if (flag.ne.0) then
      write(*,*) 'Error: AMI_PwSLikelihood returned', flag
      write(*,*) 'Inputs/output parameters were:'
      !do i = 1, NParams + 1 !for current implementation, always set hyperparameter for Planck a value
      do i = 1, NParams  
	write(*,*) i, Cube(i), Cube_out(i)
      enddo
      stop
    endif
    
    end subroutine

!=======================================================================

end module pws_wrapper
