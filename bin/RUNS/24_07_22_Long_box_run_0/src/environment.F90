!!=============================================================
!!=============================================================
!!module: environment
!! NAME
!!  environment (MMancini, SHess)
!!
!! FUNCTION
!!  Contains planetary environments list and its selector 
!!
!! NOTE
module environment

 use defs_basis
 use defs_species
 use defs_particletype
 use defs_atmospheretype
 use defs_parametre
 use atm_charge_exchange
 use atm_photoproduction
 use env_mars
 use env_mars_try
 use env_venus
 use env_moon
 use env_ganymede
 use env_titan
 use env_mercure
 use env_earth
 use time_variation
 use m_writeout

#include "q-p_common.h"

 implicit none
 private

 type(atmosphere_type),SAVE,public :: atmosphere 
 
 !--Pointers to environment functions
 ! each of these pointer must point toward a function to enable the simulation of the corresponding physical process
 procedure(init_species_mars),pointer           :: planet_species => Null()
 procedure(dealloc_mars),pointer                :: dealloc_planet => Null()
 procedure(exosphere_mars),pointer              :: exosphere => Null()
 procedure(photoproduction_generic),pointer     :: photoproduction => Null()
 procedure(charge_exchange_generic),pointer     :: charge_exchange => Null()
 procedure(create_ionosphere_ganymede),pointer  :: ionosphere => Null()
 procedure(add_b_dipole_ganymede),pointer       :: b_dipole => Null()
 procedure(feed_ionosphere_mars),pointer        :: feed_production => Null()
 procedure(temporal_change_CME),pointer         :: temporal_change => Null()
 procedure(init_change_CME),pointer             :: init_temporal_change => Null()

 !--Density array on which environment density array will point
 real(dp),public,allocatable :: prod_pp(:,:,:,:)
 real(dp),public,allocatable :: density_exo(:,:,:,:)

 public ::                   &
      select_environment,    &
      nullify_environment,   &
      init_species,          &
      add_exosphere,         &
      add_ionosphere,        &
      add_B_dipole,          &
      calc_photoproduction,  &    
      calc_charge_exchange,  &    
      feed_ionosphere,       &
      load_temporal_change,  &
      set_time_variation

contains
 !!#####################################################################

 !!=============================================================
 !!routine: environment/select_environment
 !!
 !! FUNCTION
 !!  Select environment if it exists
 !!
 !!  Replace Obsolete modules photoproduction, charge exchange, exosphere and ionosphere,
 !!  as well as the obsolete function select_init_species 
 !!
 !! IN 
 !!  planetname=the name of the planet given by command line or
 !!  the default defined in module "parametre"
 !! OUT
 !!   
 !! SIDE EFFECT
 !!  ns=number of species
 !!  atmosphere related pointer are assigned
 !!
 subroutine select_environment(planetname,ns)
 use defs_grid,only : ncm
 
  integer,intent(inout) :: ns
  character(len=*),intent(in) :: planetname

  character(len=200) :: msg

  select case(trim(planetname))
  case("mars")
   ns = 2
   planet_species  => init_species_mars
   dealloc_planet  => dealloc_mars
   exosphere       => exosphere_mars
   photoproduction => photoproduction_mars
   
 !  if (trim(ionospherename) /="") then
 !  ionosphere      => load_ionosphere_mars
 !  else
   ionosphere      => create_ionosphere_mars
 !  endif
   charge_exchange => charge_exchange_generic
   feed_production => feed_ionosphere_mars
   b_dipole        => mars_magnetic_field
   call alloc_mars(density_exo,prod_pp,atmosphere,ncm)

  case("venus")   
   ns = 2
   planet_species  => init_species_venus
   dealloc_planet  => dealloc_venus
   exosphere       => exosphere_venus
   photoproduction => photoproduction_venus
   ionosphere      => create_ionosphere_venus
   charge_exchange => charge_exchange_generic
   feed_production => feed_ionosphere_venus
   !temporal_change => temporal_change_IMF_planet
   !init_temporal_change => init_change_IMF_planet
   call alloc_venus(density_exo,prod_pp,atmosphere,ncm)
   
  case("moon")
   ns = 2
   planet_species  => init_species_moon
  case("mercure")
   ns = 2
   planet_species  => init_species_mercure
   !ionosphere      => create_ionosphere_mercure
   b_dipole        => add_b_int_mercure
   !exosphere       => exosphere_mercure
   !photoproduction => photoproduction_mercure
   !feed_production => feed_ionosphere_mercure
   !charge_exchange => charge_exchange_generic
   call alloc_mercure(density_exo,prod_pp,atmosphere,ncm)
  case("ganymede")
   ns = 2
   planet_species  => init_species_ganymede
   dealloc_planet  => dealloc_ganymede
   ionosphere      => create_ionosphere_ganymede
   b_dipole        => add_b_dipole_ganymede
   exosphere       => exosphere_ganymede
   feed_production => feed_ionosphere_ganymede
   photoproduction => photoproduction_ganymede
   charge_exchange => charge_exchange_generic
  call alloc_ganymede(density_exo,prod_pp,atmosphere,ncm)
  case("mars3try")
   ns = 3
   planet_species  => init_species_mars_try
   dealloc_planet  => dealloc_mars
   exosphere       => exosphere_mars
   photoproduction => photoproduction_mars
   ionosphere      => create_ionosphere_mars
   charge_exchange => charge_exchange_mars
   call alloc_mars(density_exo,prod_pp,atmosphere,ncm)!associate pointers
  case("titan")
   ns = 2
   planet_species  => init_species_titan
   dealloc_planet  => dealloc_titan
 !  photoproduction => photoproduction_titan
   exosphere       => exosphere_titan
   !call alloc_titan(density_exo,prod_pp,atmosphere,ncm)!associate pointers
  case("earth")
    ns = 1
    planet_species  => init_species_earth
    b_dipole        => add_b_dipole_earth
    !b_dipole        => Null()
    temporal_change => temporal_change_CME
    init_temporal_change => init_change_CME
  case default
   write(msg,'(a,3x,a,4x,3a)')ch10,&
        "ERROR: Selected Environment:",&
        "Planet '",trim(planetname),"' does Not Exist"
   call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
   stop
  end select

  if (no_env_mag.eq.1) b_dipole  => Null()

 write(msg,'(a,1x,2a,3x,a,a15,a,3x,a,i15)')ch10,&
       "___________________ Selected Environment ___________________",ch10,&
       "Planet:  ",trim(planetname),ch10,&
       "Species: ",ns
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

end subroutine select_environment

!!=============================================================
!  Following routines are used to call the routine assigned to the pointers 
subroutine nullify_environment
   if (associated(dealloc_planet)) call dealloc_planet(0)
   dealloc_planet => Null()
   planet_species => Null()
   exosphere => Null()
   photoproduction => Null()
   charge_exchange => Null()
   ionosphere => Null()
   b_dipole => Null()
   init_temporal_change => Null()
   temporal_change => Null()
end subroutine nullify_environment


!!=============================================================
!! add_exosphere: if the exosphere pointer is associated, 
!! it starts the routine computing the exosphere
!!=============================================================
subroutine add_exosphere(Spe,ncm,gstep,s_min_loc,resistivity)
 use m_timing,only     : time_get
  integer,intent(in) :: ncm(3)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  real(dp),intent(inout) :: resistivity(:,:,:)
!!!!!!!  Change Leblanc F. 01/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(species_type), intent(inout) :: Spe
!  type(species_type), intent(in) :: Spe
!!!!!!!  Change Leblanc F. 01/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  __NO_PLANET_CYCLE
  __GETTIME(72,1)!--Timer start
  if(associated(exosphere)) call exosphere(Spe,ncm,gstep,s_min_loc,resistivity)
  __GETTIME(72,2)!--Timer stop
end subroutine add_exosphere

!!=============================================================
!! feed_ionosphere: if the feed_production pointer is associated, 
!! it starts the routine adding at each time steps the
!! required number of particles (photo_prod, electron impact, ionosphere)
!!=============================================================
subroutine feed_ionosphere(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
 use m_timing,only     : time_get
  real(dp),intent(in) :: gstep(3),s_min_loc(3),s_max_loc(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(inout) ::atmosphere
  __NO_PLANET_CYCLE
  __GETTIME(78,1)!--Timer start
  if(associated(feed_production)) &
&         call feed_production(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
  __GETTIME(78,2)!--Timer stop
 end subroutine feed_ionosphere

!!=============================================================
!! add_ionosphere: if the ionosphere pointer is associated, 
!! it starts the routine computing the ionosphere and adding the
!! particles
!!=============================================================
subroutine add_ionosphere(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot)
 use m_timing,only     : time_get
  real(dp),intent(in) :: gstep(3),s_min_loc(3),s_max_loc(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  __NO_PLANET_CYCLE
  __GETTIME(76,1)!--Timer start
  if(associated(ionosphere)) &
&         call ionosphere(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
  __GETTIME(76,2)!--Timer stop
 end subroutine add_ionosphere

!!=============================================================
!! add_b_dipole: if the b_dipole pointer is associated, 
!! it starts the routine computing the magnetic field and add it
!! (even for multipolar fields)
!!=============================================================
subroutine add_b_dipole(Bfield,ncm,Spe,gstep,s_min_loc)
 use defs_arr3Dtype
  integer, intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type),intent(inout) :: Spe
  __NO_PLANET_CYCLE
  if(associated(b_dipole)) call b_dipole(Bfield,ncm,Spe,gstep,s_min_loc)
end subroutine add_b_dipole

!!=============================================================
!! calc_photoproduction: if the photoproduction pointer is associated, 
!! it starts the routine computing the photoproduction
!!=============================================================
subroutine calc_photoproduction(Spe,ncm,gstep,s_min_loc)
 use m_timing,only     : time_get
  integer,intent(in) :: ncm(3)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type), intent(in) :: Spe
  __NO_PLANET_CYCLE
 __GETTIME(73,1)!--Timer start
  if(associated(photoproduction)) &
        call photoproduction(Spe,ncm,gstep,s_min_loc,atmosphere)
 __GETTIME(73,2)!--Timer stop
end subroutine calc_photoproduction

!!=============================================================
!! calc_charge_exchange: if the charge_exchange pointer is associated, 
!! it starts the routine implementing the particle charge exchange
!!=============================================================
subroutine calc_charge_exchange(nn,kpickup,qsm,irand,ijk,&
        &                    v_p,w,Spe,particule,atmosphere)
  use defs_variable,only :t
  integer,intent(in) :: nn,ijk(3)
  integer,intent(inout) :: irand,kpickup
  real(dp),intent(in) :: qsm,v_p(3),w(8)
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(in) ::atmosphere
  __NO_PLANET_CYCLE
  if((associated(charge_exchange)).and.(t.gt.1)) &
         call  charge_exchange(nn,kpickup,&
                 qsm,irand,ijk,v_p,w,Spe,particule,atmosphere) !load environment dependent data
end subroutine calc_charge_exchange

!!=============================================================
!! init_species: if the planet_species pointer is associated, 
!! it starts the routine defining the planet environment
!!=============================================================
subroutine init_species(Spe,s_centr)
 real(dp),intent(in) :: s_centr(3)
  type(species_type),intent(inout) :: Spe
  if(associated(planet_species)) call planet_species(Spe,s_centr)
  Spe%P%ssl=planet_ssl
  Spe%P%sslat=planet_sslat
  if (atmosphere%n_exc.gt.0) &
atmosphere%exc_reactions(:)%cross_section=atmosphere%exc_reactions(:)%cross_section*dt*Spe%ref%Alfvenspeed*Spe%ref%inv_gyro*1.e9
end subroutine init_species
 !!=============================================================
 
 !!=============================================================
 !! load_temporal_change: if the temporal variation of the incident plasma
 !! is activated then it starts the routine intializing the information
 !!=============================================================
 subroutine load_temporal_change(Spe)
 type(species_type),intent(inout) :: Spe
   if(associated(init_temporal_change)) call init_temporal_change(Spe)
 end subroutine load_temporal_change
 !!=============================================================
 
  !!=============================================================
  !! set_time_variation: if the temporal variation of the incident plasma
  !! is activated then it starts the routine implementing at each time step
  !! the loaded temporal varying infromation
  !!=============================================================
  subroutine set_time_variation(by0,bz0,by1,bz1,vxmean,species,iter)
  real(dp)          , intent(inout) :: by0,bz0,by1,bz1,vxmean
  type(species_type), intent(inout) :: species
  integer           , intent(inout) :: iter
    if(associated(temporal_change)) call temporal_change(by0,bz0,by1,bz1,vxmean,species,iter)
  end subroutine set_time_variation
 !!=============================================================

end module environment
