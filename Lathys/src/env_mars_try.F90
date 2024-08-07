!!=============================================================
!!=============================================================
!!module: env_mars_try
!! NAME
!!  env_mars_try (MMancini) (name changed by SHess)
!!
!! FUNCTION
!!  Contains definition of type for the environment planet mars_try
!!
!! NOTE
module env_mars_try

 use defs_basis
 use defs_species
 use defs_parametre

#include "q-p_common.h"

 implicit none
 private

 public ::                 &
      init_species_mars_try     !--Intialise the Mars_try

contains
 !!#####################################################################


 !!=============================================================
 !!routine: m_species/init_species_mars_try
 !!
 !! FUNCTION
 !!  Initialize species for mars_try
 !! IN 
 !! OUT
 !! SIDE EFFECT
 !!  
 !!
 subroutine init_species_mars_try(species,s_centr)

  real(dp),intent(in) :: s_centr(3)
  type(species_type),intent(inout) :: species
  
  integer,parameter :: H=1,He=2,H2=3
  integer :: is
  real(dp) :: tot_tmp,vitessethermique

  !--Mars_try contains two species of particles

  !--Allocation of species_type
  call alloc_species(3,species)

  !--Name of the planet
  species%planetname = "mars_try"
  
  !--Intialize Physical parameter used here
  !--Physical density of reference (m-3)
  species%ref%density = 2.5e6  

  !--Ions inertial length
  species%ref%c_omegapi = Sp_Lt*sqrt(epsilon0*pmasse_e2/species%ref%density)/1.e3

  !--Magnetic Field (Tesla)
  species%ref%mag = 3.e-9

  !--Alfven Speed
  species%ref%alfvenspeed = species%ref%mag/&
       &                   (sqrt(mu0*amu_pmass*species%ref%density)*1.e3)

  !--Inverse of gyrofrequency
  species%ref%inv_gyro = amu_pmass/(e_Cb*species%ref%mag)

  !--Assignation of Planet values
  species%P%centr  = s_centr
#ifndef HAVE_NO_PLANET  
  species%P%radius = 3376._dp/species%ref%c_omegapi
  species%P%r_exo  = species%P%radius+200._dp/species%ref%c_omegapi
  species%P%r_iono = species%P%radius   !TO change!!!!!!!!
  species%P%r_lim  = species%P%radius+110._dp/species%ref%c_omegapi
  species%P%speed  = zero 
#endif

  !--Hydrogen and Helium
  species%S(:)%name = (/"H  ","He ","H  "/)
  
  !--Number of particles for cells
  species%S(:)%ng = (/1, 2, 1/) !(/2,2/)

  !--Direct Speeds  (vitesses dirigees?)
  species%S(H)%vxs = 10.0_dp !--  H+
  species%S(He)%vxs = species%S(H)%vxs !-- He++
  species%S(H2)%vxs = 10.0_dp !--  H+
  species%S(:)%vys  = zero
  species%S(:)%vzs  = zero

 !--Percentage of He in the solar wind: n(He++)/(n(He++)+n(H+))
  species%S(He)%percent = .05_dp
  species%S(H)%percent = (one-species%S(He)%percent)/two
  species%S(H2)%percent =  species%S(H)%percent

  !--Ratio of charges and masses and Temperatures
  species%S(:)%rcharge = (/one,two ,one/)
  species%S(:)%rmass   = (/one,four,one/)
  species%S(:)%rtemp   = (/one,four,one/)

  !--Betas
  species%S(H)%betas = 0.44_dp
  species%S(He)%betas = 0.09_dp
  species%S(H2)%betas = 0.5_dp
  species%betae = 2.79_dp

  !--Rapport des vitesses thermque entre parallel et perpendiculair  H+ et He++
  species%S(:)%rspeed = one

  !--Rapport des vitesse thermiques entre H+ et He++
  species%S(:)%rvth = one
  species%S(:)%rmds = one

  !--Rapport charge sur masse
  species%S(:)%qms  = species%S(:)%rcharge/species%S(:)%rmass

  !--Thermal speed (parallel and perpendicular)
  vitessethermique = sqrt( three*species%S(H)%betas/(one+two*species%S(H)%rvth**two))

  species%S(:)%vth1 = species%S(:)%rspeed * vitessethermique

  species%S(:)%vth2 = species%S(:)%rvth*species%S(:)%vth1

  tot_tmp = sum(species%S(:)%rmass*species%S(:)%percent)
 
  !--Macro-particle mass
  species%S(:)%sm =  species%S(:)%rmass*species%S(:)%percent/(tot_tmp*real(species%S(:)%ng,dp))

  !--Macro-particle charge 
  species%S(:)%sq = species%S(:)%qms*species%S(:)%sm

  !--Set probability of extraction
  tot_tmp = sum(real(species%S(:)%ng,dp)*species%S(:)%vxs)
  species%S(:)%prob = (real(species%S(:)%ng,dp)*species%S(:)%vxs)/tot_tmp

  !--Accumulate sum
  species%S(:)%prob = (/(sum(species%S(:is)%prob),is=1,species%ns)/)

 end subroutine init_species_mars_try
end module env_mars_try
