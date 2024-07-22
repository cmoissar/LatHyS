!!=============================================================
!!=============================================================
!!module: defs_counts_types
!! NAME
!!  defs_counts_types (MMancini)
!!
!! FUNCTION
!! Defines type for diagnostic,
!! Initialisation and Clean
module defs_counts_types

 use defs_basis
 
 implicit none
 private

 !!=============================================================
 !!--Type for counting particles
 !!  contains information concerning the ouput particles or pickups
 !!  from a proc region at any time
 type,public :: count_single_type
  integer :: out_dir(nb_voisins) !--Exiting objects in any direction
  integer :: in_dir(nb_voisins)  !--Entering objects in any direction
  integer :: out                 !--Total exiting objects
  integer :: in                  !--Total entering objects
 end type count_single_type

 !!=============================================================
 !!--Type for counting particles
 !!  contains information concerning the ouput particles AND pickups
 !!  from a proc region at any time
 type,public :: count_particle_type
  type(count_single_type) :: np !--Counter for normal particles
  type(count_single_type) :: pp !--Counter for pickups particles
  integer :: out_tot !--Total exiting particles
  integer :: in_tot  !--Total entering particles
  integer :: out_xp  !--Exiting on x-direction (positive)
  integer :: out_xm  !--Exiting on x-direction (negative)
  integer :: kabs    !--Absorbed by the planet
 end type count_particle_type
 
 
 public::                       &
      init_count_particle_type   

contains
 !!#####################################################################

 !!=============================================================
 !!routine: defs_counts_types/init_count_particle_type
 !!
 !! FUNCTION
 !!  initialize count_particle_type
 !! IN 
 !! OUT
 !!   count=initialized count_particle_type
 !! SIDE EFFECT
 subroutine init_count_particle_type(count)
  
  type(count_particle_type),intent(out) :: count
  !--Set to zero particles counters
  count%np%out_dir = 0
  count%np%in_dir  = 0
  count%np%out     = 0
  count%np%in      = 0
  !--Set to zero pickups counters
  count%pp%out_dir = 0
  count%pp%in_dir  = 0
  count%pp%out     = 0
  count%pp%in      = 0
  !--Set to zero total counters
  count%out_tot  = 0
  count%in_tot   = 0
  count%out_xp   = 0
  count%out_xm   = 0
  count%kabs     = 0
 end subroutine init_count_particle_type

end module defs_counts_types
