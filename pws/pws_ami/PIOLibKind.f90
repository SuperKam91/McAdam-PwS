
module piolibkind
  integer, parameter :: PIOFLAG = 1
  integer, parameter :: PIOBYTE = selected_int_kind(2)
  integer, parameter :: PIOSHORT = selected_int_kind(4)
  integer, parameter :: PIOINT = selected_int_kind(9)
  integer, parameter :: PIOLONG = selected_int_kind(18)
  integer, parameter :: PIOFLOAT = kind(0.0)
  integer, parameter :: PIODOUBLE = kind(0.0d0)
  integer, parameter :: DMCPIOSTRINGMAXLEN = 256



  !
  ! ---- > Pointor on array of double <-----
  !
  ! If you define a array of PIOPtrOnArrayXXX is like 2D array of  XXX 
  ! This structure is necessary for wrapper F90 from C without copy memory
  !                                                    ------------------- 
  ! 
  ! Jm. Colley, Fev. 2006
  !

  type PIOPtrOnArrayFlag
     logical(KIND=PIOFLAG), dimension(:), pointer :: IdxSple ! Index sample
  end type PIOPtrOnArrayFlag

  type PIOPtrOnArrayByte
     integer(KIND=PIOBYTE), dimension(:), pointer :: IdxSple ! Index sample
  end type PIOPtrOnArrayByte

  type PIOPtrOnArrayShort
     integer(KIND=PIOSHORT), dimension(:), pointer :: IdxSple ! Index sample
  end type PIOPtrOnArrayShort

  type PIOPtrOnArrayInt
     integer(KIND=PIOINT), dimension(:), pointer :: IdxSple ! Index sample
  end type PIOPtrOnArrayInt

  type PIOPtrOnArrayLong
     integer(KIND=PIOLONG), dimension(:), pointer :: IdxSple ! Index sample
  end type PIOPtrOnArrayLong

  type PIOPtrOnArrayFloat
     real(KIND=PIOFLOAT), dimension(:), pointer :: IdxSple ! Index sample
  end type PIOPtrOnArrayFloat

  type PIOPtrOnArrayDble
     real(KIND=PIODOUBLE), dimension(:), pointer :: IdxSple ! Index sample
  end type PIOPtrOnArrayDble
  

  TYPE PIODBLogger
     character(len=DMCPIOSTRINGMAXLEN)::prefix="NEED_INIT"
     integer::level=0
     integer(kind=PIOLONG)::id
     integer(kind=PIOLONG)::mpi_rank
     integer(kind=PIOLONG)::mpi_size
  END TYPE PIODBLogger

end module piolibkind
