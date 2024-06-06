!!=============================================================
!!=============================================================
!!module: m_species
!! NAME
!!  m_species (MMancini)
!!
!! FUNCTION
!!  Contains definition of type for species and some function  
!!  working on this type
!!
!! NOTE
module defs_species

 use defs_basis
 use defs_parametre
#include "q-p_common.h"

 implicit none
 private

 !!=============================================================
 !!--Type for species
 !!  contains information concerning a particular species
 type :: partic_species_type
  real(dp) :: betas,     &!--Beta pour chaque espece
       &       rvth,     &!--Rapport des vitesse thermique par espece
       &       vth1,     &!--Vit therm. para par espece
       &       vth2,     &!--Vit therm. perp par espece
       &       sq,       &!--Charge par espece
       &       sm,       &!--masse par espece
       &       rmds,     &!--Rapport des masses
       &       qms,      &!--Rapport charge sur masse par espece
       &       vxs,      &!--Vitesse dirigée suivant X
       &       vys,      &!--Vitesse dirigée suivant Y
       &       vzs,      &!--Vitesse dirigée suivant Z 
       &       flux_tot   !--Rapport des vitesse thermique par espece
  
  !--Ratios
  real(dp) :: percent, &!--% of other species with respect to the 1st
       &      rcharge, &!--ratio of species charge with respect to the 1st 
       &      rmass,   &!--ratio of species mass with respect to the 1st  
       &      rtemp,   &!--ratio of species temperature with respect to the 1st(THesTH)
       &      rspeed,  &!--vtHesvtH
       &      prob      !--Probabilty to extract a species during creation (pro_He)

  integer :: np,  &!--Nombre de particule  a t=0
       &     n1,  &!--Numéro de début des particules par espece
       &     n2,  &!--Numéro de fin des particules par espece
       &     ng    !--Nombre de particule par cellule a t=0

  integer :: ip_inj   !--Injected particles

  character(len=3)  :: name !--Name of the species
 
 end type partic_species_type

 !!=============================================================
 !!--Type for planet
 !!  contains information concerning the planet
 !!  the length are in units c_omegapi
 type :: planet_type
  
  real(dp) :: dipole(3) !-- position of dipole center
  real(dp) :: centr(3)  !--Position of the planet
  real(dp) :: radius,&  !--Planet Radius
       &      r_iono,&  !
       &      r_exo,&   !--Exosphere radius
       &      r_lim,&   !--Planet radius+altitude
       &      speed,&   !--Speed of planet x-direction
       &      ssl, &
       &      sslat
 end type planet_type

 !!=============================================================
 !!--Type for Physical parameters
 !!  contains information concerning the physical parameters
 !!  of dominant species used as reference
 type :: physical_param
  real(dp) :: density !-Input physical density (m^3)
  real(dp) :: mag     !--Input Magnetic Field (input)(T)
  real(dp) :: alfvenspeed !--alfven speed(Km/s)
  real(dp) :: c_omegapi   !--Ions inertial length(Km)
  real(dp) :: inv_gyro    !--Inverse of gyrofrequency(s)
  real(dp) :: maxabs      !--Max absorption distance (Km)
  real(dp) :: conv_t      !--converts temperature in phys. units
 end type physical_param

 !!=============================================================
 !!--Type for species
 !!  contains information concerning the species
 type :: species_type

  real(dp),allocatable :: vplus(:,:) !--Distribution of velocities for entering 
  type(partic_species_type), allocatable :: S(:)
  type(planet_type) :: P           !--Planet infos
  type(physical_param) :: ref      !--Reference physical parameters
  real(dp) :: betae                !--Beta electronique
  real(dp) :: tempe_ratio          !--Electron temperature ratio (Iono/SW)
  real(dp) :: viscosity            !--vicosity of the plasma (dv=dv-v_ion*density_neutral*viscosity)
  integer :: ns                    !--Number of species
  character(len=10) :: planetname
  character(len=50) :: exospherename
 
 end type species_type


 public ::                 &
      species_type,        &!--Type for species
      planet_type,         &!--Type for planet
      alloc_species,       &!--Allocate an species_type
      dealloc_species,     &!--Deallocate an species_type
      print_species         !--Print species_type

#ifdef HAVE_NETCDF
 public ::                 &   
      species_def_dim_cdf, &!--Define dimensions for species_type
      species_def_var_cdf, &!--Define variables for species_type 
      species_put_var_cdf, &!--Write variables species_type
      species_get_var_cdf   !--Read variables species_type
#else
 public ::                 &
      species_put_var_bin, &!--Binary write of species_type
      species_get_var_bin   !--Binary read of species_type
#endif

contains
 !!#####################################################################

 !!=============================================================
 !!routine: m_species/alloc_species
 !!
 !! FUNCTION
 !!  Initialize species type
 !! IN 
 !!  ns=number of species
 !! OUT
 !!  species=variable containg species
 !!   
 !! SIDE EFFECT
 subroutine alloc_species(ns,species)

  integer,intent(in) :: ns !--Number of species
  type(species_type),intent(out) :: species

  species%ns = ns
  species%planetname = ""
  species%betae = 0._dp
  species%tempe_ratio = 1._dp
  species%viscosity =0._dp

  !--Physical_param Initialisation
  species%ref%density     = zero
  species%ref%mag         = zero
  species%ref%alfvenspeed = zero
  species%ref%c_omegapi   = zero
  species%ref%inv_gyro    = zero
  species%ref%maxabs      = zero
  species%ref%conv_t      = zero

  !--Planet Initialisation
  species%P%centr(3) = zero
  species%P%dipole(3) = zero
  species%P%radius   = zero
  species%P%r_iono   = zero
  species%P%r_exo    = zero
  species%P%r_lim    = zero
  species%P%speed    = zero
  species%P%ssl      = zero
  species%P%sslat    = zero

  !--Allocation des tableaux pour les species du vent solaire
  allocate(species%S(ns))
  
  species%S(:)%np = 0 
  species%S(:)%ng = 0
  species%S(:)%n1 = 0 
  species%S(:)%n2 = 0 
  species%S(:)%ip_inj = 0
  species%S(:)%betas = zero
  species%S(:)%vth1 = zero
  species%S(:)%vth2 = zero
  species%S(:)%rvth = zero  
  species%S(:)%rmds = zero  
  species%S(:)%sq   = zero 
  species%S(:)%sm   = zero
  species%S(:)%qms  = zero 
  species%S(:)%vxs  = zero
  species%S(:)%vys  = zero 
  species%S(:)%vzs  = zero
  species%S(:)%flux_tot = zero
  species%S(:)%percent = zero 
  species%S(:)%rcharge = zero
  species%S(:)%rmass   = zero 
  species%S(:)%rtemp   = zero
  species%S(:)%rspeed  = zero
  species%S(:)%prob  = zero

  species%S(:)%name = ""

  !--Allocate Arrays for vplus speed distribution
  allocate(species%vplus(nfl,ns))

 end subroutine alloc_species

 !!=============================================================
 !!routine: m_species/dealloc_species
 !!
 !! FUNCTION
 !!  Deallocate species_type
 !! IN 
 !! OUT
 !! SIDE EFFECT
 !!  Deallocated species_type  
 !!
 subroutine dealloc_species(species)

  type(species_type),intent(inout) :: species

  !--Allocation of arrays for species of solar wind
  deallocate(species%S)   
  deallocate(species%vplus)

 end subroutine dealloc_species


 !!=============================================================
 !!routine: m_species/print_species
 !!
 !! FUNCTION
 !!  Print species_type
 !! IN 
 !! OUT
 !! SIDE EFFECT
 !!
subroutine print_species(species)
 use m_writeout

 type(species_type),intent(in) :: species

 character(len=2) :: str_ns
 character(len=20) :: format1,format2,format3,format4,format5
 character(len=2000) :: msg

 !--Complicate formatting which depends on the number of species
 write(str_ns,'(i2)') species%ns
 format1 = "(a,20x,"//trim(str_ns)//"a14)"
 format2 = "(a,4x,a,9x,"//trim(str_ns)//"i14)"
 format3 = "(4x,a,14x,"//trim(str_ns)//"f14.7)"
 format4 = "(a,4x,a,14x,"//trim(str_ns)//"f14.7)"
 format5 = "(a,4x,a,"//trim(str_ns)//"f14.7)"
 
 write(msg,"(3a,3x,"//trim(format1)//","//trim(format2)//")")ch10,&
      &" __________ Initialisation of the Species ___________________",&
      &ch10,"Species:",species%S(:)%name,&
      &ch10,"number for cell ng ",species%S(:)%ng
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

 write(msg,"("//trim(format3)//",10"//trim(format4)//",2"//trim(format5)//")")&
      &     "betas         ",species%S(:)%betas,&
      &ch10,"mass ratio    ",species%S(:)%rmass,&
      &ch10,"charge ratio  ",species%S(:)%rcharge,&
      &ch10,"qms           ",species%S(:)%qms,&
      &ch10,"percentage    ",species%S(:)%percent,&
      &ch10,"sm            ",species%S(:)%sm,&
      &ch10,"sq            ",species%S(:)%sq,&
      &ch10,"vxs           ",species%S(:)%vxs,&
      &ch10,"rvth          ",species%S(:)%rvth,&
      &ch10,"rmds          ",species%S(:)%rmds,&
      &ch10,"prob. extrac. ",species%S(:)%prob,&
      &ch10,"parallel thermal speed      ",species%S(:)%vth1,&
      &ch10,"perpendicular thermal speed ",species%S(:)%vth2
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

 write(msg,"(2(a,4x,a,14x,f14.7))")&
      &ch10,"sum(psm)      ",sum(species%S(:)%sm*real(species%S(:)%ng,dp)),&
      &ch10,"sum(psq)      ",sum(species%S(:)%sq*real(species%S(:)%ng,dp))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

#ifndef HAVE_NO_PLANET
 write(msg,"(a,3x,a,a40,(a,4x,a,14x,3f8.2),4(a,4x,a,14x,en14.3))")&
      &ch10,"Planet:",species%planetname,&
      &ch10,"position          ",species%P%centr,&
      &ch10,"radius            ",species%P%radius,&
      &ch10,"limit radius      ",species%P%r_lim,&
      &ch10,"exo limit         ",species%P%r_exo,&
      &ch10,"speed             ",species%P%speed
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

 write(msg,"(a,3x,a,a,6(a,4x,a,14x,en14.3,a))")&
      &ch10,"Physical Parameters (dominant): ",species%S(1)%name,&
      &ch10,"density           ",species%ref%density," (m^-3)",&
      &ch10,"magnetic field    ",species%ref%mag," (T)",&
      &ch10,"alfven speed      ",species%ref%alfvenspeed," (Km/s)",&
      &ch10,"Ions inertial len ",species%ref%c_omegapi," (Km)",&
      &ch10,"1/gyrofrequency   ",species%ref%inv_gyro," (s)",&
      &ch10,"dist absorption   ",species%ref%maxabs," (Km)"
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
#else
 write(msg,"(a,3x,2a,4x,a)")&
      &ch10,"Planet:",&
      &ch10,"No Planet Mode: quiet plasma!"
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
#endif

end subroutine print_species

#ifdef HAVE_NETCDF
 !!=============================================================
 !!subroutine: m_species/species_def_dim_cdf
 !! FUNCTION 
 !!
 !! INPUT
 subroutine species_def_dim_cdf(ncid,dimspec)
  use netcdf
  use defs_basic_cdf,only    : test_cdf

  integer,intent(in) :: ncid
  integer,intent(inout),allocatable :: dimspec(:)

  integer :: stId

  allocate(dimspec(6)) 

  stId = nf90_def_dim(ncid, "dim_1", 1, dimspec(1))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "ns" , ns , dimspec(2))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "nfl", nfl, dimspec(3))
  call test_cdf(stId)
  !--Length of species names
  stId = nf90_def_dim(ncid, "dimname",20, dimspec(4))
  call test_cdf(stId)
  !--Length of planet names
  stId = nf90_def_dim(ncid, "dimplanet",20, dimspec(5))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "dim_3", 3, dimspec(6))
  call test_cdf(stId)


 end subroutine species_def_dim_cdf

 !!=============================================================
 !!subroutine: m_species/species_def_var_cdf
 !! FUNCTION 
 !!
 !! INPUT
 subroutine species_def_var_cdf(ncid,varid,dimspec,ii)
  use netcdf
  use defs_basic_cdf,only    : test_cdf

  integer,intent(in)   :: ncid
  integer,intent(inout):: ii
  integer,intent(inout),dimension(:) :: varid,dimspec

  integer :: stId

  stId = nf90_def_var(ncid, "ns",nf90_int,dimspec(1), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "nfl",nf90_int,dimspec(1), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "planetname",nf90_char,dimspec(5), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "speciesname",nf90_char,(/dimspec(4),dimspec(2)/), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "phys_density",QP_NF90_DP,dimspec(1), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "phys_mag",QP_NF90_DP,dimspec(1), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "phys_speed",QP_NF90_DP,dimspec(1), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "phys_length",QP_NF90_DP,dimspec(1), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "phys_time",QP_NF90_DP,dimspec(1), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "phys_temp_e",QP_NF90_DP,dimspec(1), varid(ii))
  call test_cdf(stId); ii = ii+1  
  stId = nf90_def_var(ncid, "absorpt_len",QP_NF90_DP,dimspec(1), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "s_centr",  QP_NF90_DP,dimspec(6),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "r_planet",  QP_NF90_DP,dimspec(1),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "r_lim",  QP_NF90_DP,dimspec(1),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "r_exo",  QP_NF90_DP,dimspec(1),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "r_iono",  QP_NF90_DP,dimspec(1),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "v_planet", QP_NF90_DP,dimspec(6),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "np",nf90_int,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "ng",nf90_int,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "n1",nf90_int,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1           
  stId = nf90_def_var(ncid, "n2",nf90_int,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid,"betae",QP_NF90_DP,dimspec(1), varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid,"betas",QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "rvth",QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "vth1",QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "vth2",QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "rmds",QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "qms", QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "sq" , QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "sm" , QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "vxs", QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1             
  stId = nf90_def_var(ncid, "vys", QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1             
  stId = nf90_def_var(ncid, "vzs", QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "flux_tot", QP_NF90_DP,dimspec(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "percent",         QP_NF90_DP,dimspec(2),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "ratios_charges",  QP_NF90_DP,dimspec(2),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "ratios_masses",   QP_NF90_DP,dimspec(2),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "THesTH",          QP_NF90_DP, dimspec(2),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "vtHesvtH",        QP_NF90_DP, dimspec(2),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "prob_extract",    QP_NF90_DP, dimspec(2),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "vplus",  QP_NF90_DP,(/dimspec(3),dimspec(2)/),varid(ii))
  call test_cdf(stId); ii = ii+1

 end subroutine species_def_var_cdf


 !!=============================================================
 !!subroutine: m_species/species_put_var_cdf
 !! FUNCTION 
 !!
 !! INPUT
 subroutine species_put_var_cdf(species,ncid,varid,dimspec,ii)

  use defs_basic_cdf,only     : test_cdf
  use netcdf


  type(species_type),intent(in) :: species
  integer,intent(in)   :: ncid
  integer,intent(inout):: ii
  integer,intent(in),dimension(:) :: varid
  integer,intent(inout),allocatable:: dimspec(:)

  integer :: stId,is
  character(len=20) :: arrayname(species%ns)

  !--Write the species name in a tmp array
  do is=1,species%ns
   arrayname(is) = species%S(is)%name
  enddo

  stId = nf90_put_var(ncid, varid(ii), species%ns)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), size(species%vplus(:,1)))
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%planetname)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), arrayname)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%ref%density)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%ref%mag)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%ref%alfvenspeed)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%ref%c_omegapi)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%ref%inv_gyro)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), species%ref%conv_t)
  call test_cdf(stId); ii = ii+1   
  stId = nf90_put_var(ncid, varid(ii), species%ref%maxabs)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), species%P%centr)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%P%radius)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%P%r_lim)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%P%r_exo)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%P%r_iono)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%P%speed)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%np)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%ng)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%n1)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%n2)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%betae)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%betas)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%rvth)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%vth1)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%vth2)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%rmds)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%qms)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%sq)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%sm)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%vxs)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%vys)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%vzs)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%flux_tot)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%percent)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%rcharge)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%rmass)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%rtemp)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%rspeed)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), species%S(:)%prob)
  call test_cdf(stId); ii = ii+1
  
  stId = nf90_put_var(ncid, varid(ii), species%vplus)
  call test_cdf(stId); ii = ii+1


  deallocate(dimspec)

 end subroutine species_put_var_cdf


 !!=============================================================
 !!routine: m_species/species_get_var_cdf
 !!
 !! FUNCTION
 !!  Get variable speciestype from a netcdf file
 !!         
 !! IN
 !! OUT
 subroutine species_get_var_cdf(species,ncid)
  use netcdf
  use defs_basic_cdf,only     : test_cdf,get_simple_variable_cdf

  
  type(species_type),intent(inout) :: species
  integer,intent(in) :: ncid

  integer :: nfl,is
  character(len=20) :: arrayname(species%ns)

  call get_simple_variable_cdf(ncid,"ns",  species%ns)
  call get_simple_variable_cdf(ncid,"nfl", nfl)
  call get_simple_variable_cdf(ncid,"planetname",species%planetname)
  call get_simple_variable_cdf(ncid,"phys_density",species%ref%density)
  call get_simple_variable_cdf(ncid,"phys_mag",species%ref%mag)
  call get_simple_variable_cdf(ncid,"phys_speed",species%ref%alfvenspeed)
  call get_simple_variable_cdf(ncid,"phys_length",species%ref%c_omegapi)
  call get_simple_variable_cdf(ncid,"phys_time",species%ref%inv_gyro)
  call get_simple_variable_cdf(ncid,"absorpt_len",species%ref%maxabs)
  call get_simple_variable_cdf(ncid,"s_centr",species%P%centr)
  call get_simple_variable_cdf(ncid,"r_planet",species%P%radius)
  call get_simple_variable_cdf(ncid,"r_lim",species%P%r_lim)
  call get_simple_variable_cdf(ncid,"r_exo",species%P%r_exo)
  call get_simple_variable_cdf(ncid,"r_iono",species%P%r_iono)
  call get_simple_variable_cdf(ncid,"v_planet",species%P%speed)
  call get_simple_variable_cdf(ncid,"np", species%S(:)%np)
  call get_simple_variable_cdf(ncid,"ng", species%S(:)%ng)              
  call get_simple_variable_cdf(ncid,"n1", species%S(:)%n1)
  call get_simple_variable_cdf(ncid,"n2", species%S(:)%n2)
  call get_simple_variable_cdf(ncid,"betae", species%betae)
  call get_simple_variable_cdf(ncid,"betas", species%S(:)%betas)
  call get_simple_variable_cdf(ncid,"rvth",  species%S(:)%rvth)
  call get_simple_variable_cdf(ncid,"vth1",  species%S(:)%vth1)
  call get_simple_variable_cdf(ncid,"vth2",  species%S(:)%vth2)
  call get_simple_variable_cdf(ncid,"rmds",  species%S(:)%rmds)
  call get_simple_variable_cdf(ncid,"qms",   species%S(:)%qms)
  call get_simple_variable_cdf(ncid,"sq",    species%S(:)%sq)
  call get_simple_variable_cdf(ncid,"sm",    species%S(:)%sm)
  call get_simple_variable_cdf(ncid,"vxs",   species%S(:)%vxs)
  call get_simple_variable_cdf(ncid,"vys",   species%S(:)%vys)
  call get_simple_variable_cdf(ncid,"vzs",   species%S(:)%vzs)
  call get_simple_variable_cdf(ncid,"flux_tot", species%S(:)%flux_tot)
  call get_simple_variable_cdf(ncid,"percent", species%S(:)%percent) 
  call get_simple_variable_cdf(ncid,"ratios_charges", species%S(:)%rcharge) 
  call get_simple_variable_cdf(ncid,"ratios_masses",  species%S(:)%rmass) 
  call get_simple_variable_cdf(ncid,"THesTH",   species%S(:)%rtemp)
  call get_simple_variable_cdf(ncid,"vtHesvtH", species%S(:)%rspeed)
  call get_simple_variable_cdf(ncid,"prob_extract", species%S(:)%prob)
  call get_simple_variable_cdf(ncid,"vplus", species%vplus)

  
  !--Get the name of species using a tmp array
  call get_simple_variable_cdf(ncid,"speciesname",arrayname)
  do is=1,species%ns
   species%S(is)%name = trim(arrayname(is))
  enddo

 end subroutine species_get_var_cdf

#else

 !!=============================================================
 !!subroutine: m_species/species_put_var_bin
 !! FUNCTION 
 !!
 !! INPUT
 subroutine species_put_var_bin(species,unit)
  type(species_type),intent(in) :: species
  integer,intent(in) :: unit

  write(unit) species%ns,size(species%vplus(:,1))
  write(unit) species%planetname,&
       &      species%S(:)%name
  write(unit) species%ref%density,species%ref%mag,&
       &      species%ref%alfvenspeed,&
       &      species%ref%c_omegapi,&
       &      species%ref%inv_gyro,&
       &      species%ref%maxabs
  write(unit) species%P%centr,&
       &      species%P%radius,&
       &      species%P%r_lim,&
       &      species%P%r_exo,&
       &      species%P%r_iono,&
       &      species%P%speed
  write(unit) species%S(:)%np,&
       &      species%S(:)%ng,&
       &      species%S(:)%n1,&
       &      species%S(:)%n2
  write(unit) species%betae,species%S(:)%betas, species%S(:)%rvth
  write(unit) species%S(:)%vth1,species%S(:)%vth2
  write(unit) species%S(:)%rmds,species%S(:)%qms
  write(unit) species%S(:)%sq, species%S(:)%sm
  write(unit) species%S(:)%vxs,&
       &      species%S(:)%vys,&
       &      species%S(:)%vzs,&
       &      species%S(:)%flux_tot

  write(unit) species%S(:)%percent,&
       &      species%S(:)%rcharge,&
       &      species%S(:)%rmass  ,&
       &      species%S(:)%rtemp  ,&
       &      species%S(:)%rspeed, &
       &      species%S(:)%prob
  write(unit) species%vplus

 end subroutine species_put_var_bin
 !!=============================================================
 !!subroutine: m_species/species_get_var_bin
 !! FUNCTION 
 !!
 !! INPUT
 subroutine species_get_var_bin(species,unit)
  type(species_type),intent(inout) :: species
  integer,intent(in) :: unit

  integer :: nfl

  read(unit) species%ns,nfl
  read(unit) species%planetname,species%S(:)%name
  read(unit) species%ref%density,species%ref%mag,&
       &     species%ref%alfvenspeed,&
       &     species%ref%c_omegapi,&
       &     species%ref%inv_gyro,&
       &     species%ref%maxabs
  read(unit) species%P%centr,&
       &     species%P%radius,&
       &     species%P%r_lim,&
       &     species%P%r_exo,&
       &     species%P%r_iono,&
       &     species%P%speed
  read(unit) species%S(:)%np,&
       &     species%S(:)%ng,&
       &     species%S(:)%n1,&
       &     species%S(:)%n2
  read(unit) species%betae,species%S(:)%betas, species%S(:)%rvth
  read(unit) species%S(:)%vth1,species%S(:)%vth2
  read(unit) species%S(:)%rmds,species%S(:)%qms
  read(unit) species%S(:)%sq, species%S(:)%sm
  read(unit) species%S(:)%vxs,&
       &     species%S(:)%vys,&
       &     species%S(:)%vzs,&
       &     species%S(:)%flux_tot
  read(unit) species%S(:)%percent,&
       &     species%S(:)%rcharge,&
       &     species%S(:)%rmass  ,&
       &     species%S(:)%rtemp  ,&
       &     species%S(:)%rspeed, &
       &     species%S(:)%prob
  read(unit) species%vplus

 end subroutine species_get_var_bin

#endif

end module defs_species
