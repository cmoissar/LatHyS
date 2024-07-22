!!=============================================================
!!=============================================================
!!module: defs_atmospheretype
!! NAME
!!  defs_defs_atmospheretype (SHess)
!!
!! FUNCTION
!!  Contains the type of atmosphere
!!  and relied functions
!! NOTE
!!  Use standar f2003. 
module defs_atmospheretype
 use defs_basis

 implicit none
 private

 !!=============================================================
 type :: atm_species
  character(len=10) :: name
  character(len=500)::description
  logical  :: iono,opaque
  real(dp) :: mass, charge
  real(dp),dimension(:,:,:),pointer :: density=> Null(),prod=> Null()
  real(dp),dimension(:),allocatable :: ion_abs
 end type atm_species 

 type :: atm_photoreaction
  character(len=500)::description
  type(atm_species),pointer :: mother=> Null(),daughter=> Null()
  real(dp),dimension(:),allocatable :: ion_react
  real(dp) :: frequency
 end type atm_photoreaction

 type :: atm_exc_reaction
  character(len=500)::description
  type(atm_species),pointer :: ion=> Null(),neutral=> Null()
  real(dp) :: cross_section,qsm
 end type atm_exc_reaction

 type :: atm_ei_reaction
  character(len=500)::description
  type(atm_species),pointer :: ion=> Null(),neutral=> Null()
  real(dp) :: coeff(10)  ! coefficients of the polynomial sigma=sum(c_n*t^n)
  integer :: n ! number of coefficients (degree +1)
 end type atm_ei_reaction

 type,public :: atmosphere_type
   integer :: n_species,n_pp,nb_lo,n_exc,n_spe_pp,n_pp_freq_fixed,n_ei
   real(dp)                                             :: F107,F107_Avg
   real(dp),dimension(:),allocatable                    :: EUVFLX
   type(atm_species),dimension(:),allocatable           :: species
   type(atm_photoreaction),dimension(:),allocatable     :: photo_reactions
   type(atm_exc_reaction),dimension(:),allocatable      :: exc_reactions
   type(atm_ei_reaction),dimension(:),allocatable       :: ei_reactions
 end type atmosphere_type

public ::  &
    allocate_atmosphere
 !!=============================================================
contains

subroutine allocate_atmosphere(ncm,atmosphere,density,prod_pp)
integer,intent(in) ::ncm(3)
real(dp),dimension(:,:,:,:),allocatable,intent(inout),target ::density
real(dp),dimension(:,:,:,:),allocatable,intent(inout),target ::prod_pp
type(atmosphere_type),intent(inout),target :: atmosphere 
integer ::i
        allocate(density(ncm(1),ncm(2),ncm(3),atmosphere%n_species)); density(:,:,:,:)=zero

    allocate(atmosphere%species(1:atmosphere%n_species))
       ! commented RM problem when we have more ionosphereic species than planetary species
       ! if(atmosphere%n_spe_pp.gt.0) then 
       !   allocate(prod_pp(ncm(1),ncm(2),ncm(3),atmosphere%n_spe_pp)); prod_pp(:,:,:,:)=zero
       ! endif
        allocate(prod_pp(ncm(1),ncm(2),ncm(3),atmosphere%n_pp)); prod_pp(:,:,:,:)=zero
       
        if(atmosphere%n_pp.gt.0) allocate(atmosphere%photo_reactions(1:atmosphere%n_pp))
        if(atmosphere%n_exc.gt.0) allocate(atmosphere%exc_reactions(1:atmosphere%n_exc))
        if(atmosphere%n_ei.gt.0) allocate(atmosphere%ei_reactions(1:atmosphere%n_ei))

        do i=1,atmosphere%n_species
           atmosphere%species(i)%density=>density(:,:,:,i)
           if(atmosphere%nb_lo.gt.0) allocate(atmosphere%species(i)%ion_abs(atmosphere%nb_lo))
         enddo
        if(atmosphere%n_pp.gt.0) then
        allocate(atmosphere%EUVFLX(atmosphere%nb_lo))
        do i=1,atmosphere%n_pp
         allocate(atmosphere%photo_reactions(i)%ion_react(atmosphere%nb_lo))
        enddo
        endif
        do i=1,atmosphere%n_species
                write(atmosphere%species(i)%description,'(a)') ' '
        enddo
        do i=1,atmosphere%n_pp
                write(atmosphere%photo_reactions(i)%description,'(a)') ' '
        enddo
        do i=1,atmosphere%n_exc
                write(atmosphere%exc_reactions(i)%description,'(a)') ' '
        enddo
        do i=1,atmosphere%n_ei
                write(atmosphere%ei_reactions(i)%description,'(a)') ' '
        enddo
end subroutine allocate_atmosphere

end module defs_atmospheretype
