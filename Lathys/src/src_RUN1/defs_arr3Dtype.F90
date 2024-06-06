!!=============================================================
!!=============================================================
!!module: defs_arr3Dtype
!! NAME
!!  defs_arr3Dtype (MMancini)
!!
!! FUNCTION
!!  Contains the type of 3-dimensional vector
!!  and relied functions
!! NOTE
!!  Use standar f2003. arr3Dtype is defined by allocatable array
module defs_arr3Dtype

 use defs_basis

 implicit none
 private


 !!=============================================================
 !!--Type for three-dimensional vector
 type,public :: arr3Dtype
  real(dp),allocatable,dimension(:,:,:) :: x,y,z
 end type arr3Dtype

 public ::              &
      alloc_arr3D,      &!--Allocate an arr3Dtype
      dealloc_arr3D      !--Deallocate an arr3Dtype

contains
 !!#####################################################################

 !!=============================================================
 !!routine: defs_arr3Dtype/alloc_arr3D
 !!
 !! FUNCTION
 !!  allocate arr3Dtype
 !! IN 
 !!  ncm(3)=size of the grid in the 3-dimensions
 !! OUT
 !!   
 !! SIDE EFFECT
 !!  arr3D_in=allocated array
 !!  
 subroutine alloc_arr3D(arr3D_in,ncm)

  integer,intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: arr3D_in
  
  if(.not.(allocated(arr3D_in%x).or.allocated(arr3D_in%y).or.allocated(arr3D_in%z))) then
   allocate(arr3D_in%x(ncm(1),ncm(2),ncm(3)),&
        &   arr3D_in%y(ncm(1),ncm(2),ncm(3)),&
        &   arr3D_in%z(ncm(1),ncm(2),ncm(3)))
  else
   stop " Trying to allocate an arr3Dtype yet allocated"
  endif
  arr3D_in%x = zero
  arr3D_in%y = zero
  arr3D_in%z = zero

 end subroutine alloc_arr3D

 !!=============================================================
 !!routine: defs_arr3Dtype/dealloc_arr3D
 !!
 !! FUNCTION
 !!  deallocate arr3Dtype
 !! IN 
 !! OUT
 !!   
 !! SIDE EFFECT
 !!  arr3D_in=deallocated array
 subroutine dealloc_arr3D(arr3D_in)

  type(arr3Dtype),intent(inout) :: arr3D_in

  if((allocated(arr3D_in%x).and.allocated(arr3D_in%y).and.allocated(arr3D_in%z))) then
   deallocate(arr3D_in%x,arr3D_in%y,arr3D_in%z)
  else
   stop " Trying to deallocate an arr3Dtype not allocated"
  endif
 end subroutine dealloc_arr3D


end module defs_arr3Dtype
