!=======================================================================

program test_p_alpha

   ! A wrapper for testing subroutines
   !   init_p_alpha(prior_min,prior_max,nkernel,kernel) and
   !   query_p_alpha(alpha,nkernel,kernel,palpha,pcumalpha)

   !use kind_def

   use p_alpha

   implicit none

   double precision prior_min,prior_max
   integer nkernel,ikernel,i
   double precision, allocatable :: kernel(:,:)
   double precision alpha,palpha,pcumalpha

!..Retrieve prior limits from user:

   prior_min=-1.0d0
   prior_max=2.5d0

   write(*,*) 'Cambridge convention: alpha>0 falling spectrum'
   write(*,*) 'What prior range for alpha: prior_min,prior_pmax (e.g. -1.0d0,2.5d0)'
   read(*,*) prior_min,prior_max

!..Initialize kernel, given the prior limits, and return it

   write(*,*) 'Normalizing kernel to this prior range'
   call init_p_alpha(prior_min,prior_max,nkernel,kernel)

!   do ikernel=1,nkernel
!      write(*,*) (kernel(ikernel,i),i=1,3)
!   enddo

!..Query kernel for p(alpha)

   write(*,*) 'Give a cumulative probability of alpha, 0<cumul(alpha)<1'
   read(*,*) pcumalpha

   !alpha=1.43d0
   call query_p_alpha(pcumalpha,nkernel,kernel,alpha,palpha)

   write(*,*) 'pcumalpha=', pcumalpha
   write(*,*) 'palpha=', palpha
   write(*,*) 'alpha=', alpha

end program test_p_alpha

!=======================================================================
