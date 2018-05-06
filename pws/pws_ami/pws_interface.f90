module pws_interface

  use iso_c_binding

  implicit none

  interface
    subroutine pws_init_c(filename, fpatch, lpatch) bind(c, name='pws_init')
      use iso_c_binding
      character(kind=c_char) :: filename(:)
	  integer :: fpatch, lpatch
    end subroutine pws_init_c
  end interface

contains

  subroutine pws_init(filename, fpatch, lpatch)
    character(len=*), intent(in) :: filename
	integer, intent(in) :: fpatch, lpatch

    integer :: i, n
    character(kind=c_char), allocatable :: tmp

    n = len(filename)
    allocate(tmp(n+1))
    do i = 1, n
      tmp(i) = filename(i:i)
    end do
    tmp(n+1) = c_null_char

    call pws_init_c(tmp, fpatch, lpatch)

    deallocate(tmp)

  end subroutine pws_init

end module pws_interface
