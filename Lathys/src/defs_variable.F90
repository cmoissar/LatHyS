!!=============================================================
!!=============================================================
module defs_variable

 ! Ce module comprend toutes les définitions de variable.
 use defs_basis
 use defs_arr3Dtype
 use defs_particletype
 use defs_species,only             : species_type
 use defs_diag_type,only        : diag_type
 use defs_counts_types,only     : count_particle_type
 !use defs_parametre
 use defs_tregister,only           : treg_type

 implicit none
 public

 type(particletype),dimension(:),allocatable,save :: particule

 type(arr3Dtype),save :: Efield      !--Electric field
 type(arr3Dtype),save :: Bfield      !--Magnetic field
 type(arr3Dtype),save :: Bfield_0      !--Magnetic field
 type(arr3Dtype),save :: Bfield_h    !--Magnetic field tmp
 type(arr3Dtype),save :: vel         !--Speed array

 !--Tmp arrays
 type(arr3Dtype),save :: vela,velf  ! Tmp speed array

 real(dp),allocatable,save :: dn(:,:,:) !--Tableaux de densités 
 
! ---Tableaux de densites electroniques planetaire et incidente
 real(dp),allocatable,save :: dn_e_pl(:,:,:) 
 real(dp),allocatable,save :: dn_e_incdt(:,:,:) 

 !--mais aussi tableaux pour avancer les courants
 !  (cf utilisation du code)
 real(dp),allocatable,save :: dna(:,:,:) !--Tableaux temporaire
 real(dp),allocatable,save :: dnb(:,:,:) !--Tableaux temporaire
 real(dp),allocatable,save :: dnf(:,:,:) !--Tableaux temporaire
 real(dp),allocatable,save :: dnsa(:,:,:)!--Tableaux temporaire
 real(dp),allocatable,save :: dnsb(:,:,:)!--Tableaux temporaire
 real(dp),allocatable,save :: pe(:,:,:)  !--Tableau de pression

 !--Tableaux de résistivité et de taux d'ionisation de l'oxygène
 real(dp),allocatable,save :: resistivity(:,:,:)
 real(dp),allocatable,save :: viscosity(:,:,:)

 !--Determination de la position du choc
 integer,allocatable,save :: pos_choc(:,:) 

 type(species_type),save :: Spe !--Variable for species

 type(treg_type),save :: r_tm  !--Time results variables

 integer,save :: iter !--Iteration du pas de temps
 real(dp),save :: t  !--Temps

 integer,allocatable :: iwr(:) !--Tableau de lissage (nwrm)

 !--Counters for particles
 integer,save :: nptot  !--Nombre de particule totale dans la boite
 integer,save,allocatable :: index_exit(:) !--Indices of out-going particles
 type(count_particle_type) :: Preg         !--Counter of particle outgoing.

 real(dp),save,dimension(3) :: &
      &                         s_min,    &!--Global min in any direction
      &                         s_max,    &!--Global max in any direction
      &                         s_r,      &!--Length in any direction (global)  
      &                         s_min_loc,&!--Min in any direction for any proc
      &                         s_max_loc,&!--Max in any direction for any proc
      &                         e_conv,   &
      &                         e1_conv

 !--Differents parametres plasma-planet
 real(dp),save :: &
      &          vac,rmu0,eps0,rmu0i, &
      &          rho0,cd0,cs,betai,te, &
      &          b0,b1,bx0,by0,bz0,by1,bz1, &
      &          vx0,vy0,vz0, &
      &          theta, vxmean
  
 !--Differents parametres Planet
 real(dp),save :: tau_ionis

 integer,save,allocatable :: diag(:),diag_pickup(:)

! temporal varying parameter
real(dp),dimension(:),allocatable :: By_CME, Bz_CME, By1_CME, Bz1_CME
real(dp),dimension(:),allocatable :: V_CME, Vy_alfven, Vz_alfven
real(dp),dimension(:),allocatable :: Ng_CME
real(dp),dimension(:),allocatable :: vth1_CME, vth2_CME


!-- Diag for distribution function
real(dp),dimension(:,:,:,:,:,:),allocatable :: distrib_func
integer :: n_Ene_distrib = 32,n_theta_distrib = 9,n_phi_distrib = 18,iter_start_distrib
real(dp),dimension(:),allocatable :: Energy_distrib,Theta_distrib,Phi_distrib

 type(diag_type),public,save :: diag_info !--Info for diagnostic
    

! ------ Variables relative to impacting particles flux
integer,save :: n_inte=64,n_inta=30,nphi=80,ntheta=40
integer,save :: nspecies,nionprod_start
real(dp),save :: alt_part_imp=250.d0 ! altitude in km
real(dp),save :: ener_min=2.d0,ener_max
real(dp),save :: dE,dE_E=0.17
real(dp),allocatable,dimension(:),save :: esurf_grid,asurf_grid,phi_grid,&
            theta_grid
real(dp),allocatable,dimension(:,:,:),save  :: flux_surf
real(dp),allocatable,dimension(:,:,:,:),save :: fdde_surf,fdda_surf


end module defs_variable

! contains
!  !!#####################################################################

!  !!=============================================================
!  !!subroutine: variable/test_file_ex
!  !! FUNCTION 
!  !!  Test if a file exist
!  !!
!  !! INPUT
!  !!  namfilrest=file name to test
!  !! OUTPUT
!  !!  file_e=.T. if the file exists, .F. otherwise
!  subroutine test_file(namfilrest,file_e)
!   character(len=*),intent(in) :: namfilrest
!   logical,intent(out) :: file_e

!   !--Test the file existence
!   inquire( file=trim(namfilrest), exist=file_e )
!   if( file_e ) then
!    write(*,*) "dir exists!"
!   else
!    write(*,*) "does not exists!"
!   end if
!  end subroutine test_file



