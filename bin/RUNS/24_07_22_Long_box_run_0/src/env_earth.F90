!!=============================================================
!!=============================================================
!!module: env_earth
!! NAME
!!  env_earth (SHess)
!!
!! Gathers all the earth related file into a single one, for more clarity
!!
!! FUNCTION
!!  Contains definition of type for the environment planet earth
!!
!! NOTE
module env_earth

 use defs_basis
 use defs_species
 use defs_atmospheretype
 use defs_parametre
 use m_writeout

#include "q-p_common.h"

 implicit none
 private

 public::                     &

      init_species_earth,      &!--Intialise the earth environment
      add_b_dipole_earth
contains
 !!#####################################################################
 !!=============================================================
 !!routine: env_earth/init_species_earth
 !!
 !! FUNCTION
 !!  Initialize species for earth
 !! IN 
 !! planet center
 !!
 !! OUT
 !! species 
 !!
 subroutine init_species_earth(species,s_centr)

  real(dp),intent(in) :: s_centr(3)
  type(species_type),intent(inout) :: species
  
  integer,parameter :: H=1
  integer :: is
  real(dp) :: tot_tmp,vitessethermique,V_SW

  !--earth contains two species of particles
  !--Allocation of species_type
  call alloc_species(1,species)
  
  !--Name of the planet
  species%planetname = "earth"
  
  !--Intialize Physical parameter used here
  !--Physical density of reference (m-3)
  !-- Densite de l'espece dominante dans le plasma incident
  !-- earth aphelie:30.4e6
  !-- earth perihelie:69.35e6
  species%ref%density = 6e6  !14-12-2006

  !--Ions inertial length (km)
  species%ref%c_omegapi = Sp_Lt*sqrt(epsilon0*pmasse_e2/species%ref%density)/1.e3 !Sp_Lt: speed of light

  !--Magnetic Field (Tesla)
  !-- earth aphelie =21nT
  !-- earth perihelie= 46 nt
  species%ref%mag = 10e-9 !3.75e-9    !14-12-2006

  !--Alfven Speed
  species%ref%alfvenspeed = species%ref%mag/&
       &                   (sqrt(mu0*amu_pmass*species%ref%density)*1.e3)

  !--Inverse of gyrofrequency
  species%ref%inv_gyro = amu_pmass/(e_Cb*species%ref%mag)

  !--Max absorption length (km)
  species%ref%maxabs = 0.0_dp

  !--Assignation of Planet values
  species%P%centr  = s_centr
  !write(*,*)s_centr
#ifndef HAVE_NO_PLANET
  !species%P%radius = 20._dp!24._dp/species%ref%c_omegapi
  species%P%radius = 14._dp!24._dp/species%ref%c_omegapi
  species%P%r_exo  = species%P%radius+0._dp/species%ref%c_omegapi
  species%P%r_iono = species%P%radius   !TO change!!!!!!!!
  species%P%r_lim  = species%P%radius+0._dp/species%ref%c_omegapi
  species%P%speed  = zero 
#endif

  !--Hydrogen and Helium
  species%S(:)%name = (/"H  "/)
  
  !--Number of particles for cells
  species%S(:)%ng = (/10/) !(/6, 6/) !(/2,2/)

  !--Direct Speeds  (vitesses dirigees? en km/s)
  V_SW = 400.0_dp/species%ref%alfvenspeed !-- Solar wind speed, value 14-12-2006

  species%S(H)%vxs = V_SW !--  H+ 
  species%S(:)%vys = zero
  species%S(:)%vzs = zero
 
  !--Percentage of He in the solar wind: n(He++)/(n(He++)+n(H+))
  species%S(H)%percent = one

  !--Ratio of charges and masses and Temperatures
  !-- First column always one
  species%S(:)%rcharge = (/one/) !(/one,two/)
  species%S(:)%rmass   = (/one/) !(/one,four/)
  species%S(:)%rtemp   = (/one/) !(/one,four/)

  !--Betas
  !-- earth Aphelie 0.32 (H) 0.06 (He) 0.47 (e)
  !-- earth perihelie  (H) 
  species%S(H)%betas = 0.5_dp
  species%betae = 0.5_dp

  !--Rapport des vitesse thermiques entre H+ et He++
  species%S(:)%rspeed = one

  !--Rapport des vitesses thermique entre parallel et perpendiculair  H+ et He++
  species%S(:)%rvth = one

  !--Rapport des masses
  species%S(:)%rmds = one

  !--Rapport charge sur masse
  species%S(:)%qms  = species%S(:)%rcharge/species%S(:)%rmass

  !--Thermal speed (parallel and perpendicular) 
  vitessethermique = sqrt( three*species%S(H)%betas/(one+two*species%S(H)%rvth**two))

  species%S(:)%vth1 = species%S(:)%rspeed * vitessethermique

  species%S(:)%vth2 = species%S(:)%rvth * species%S(:)%vth1

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

 end subroutine init_species_earth 
!****************************** END INIT_SPECIES ************************************


 !!=============================================================
 !!subroutine: b_dipole/add_b_dipole_earth
 !! NAME
 !!  add_b_dipole_earth (RModolo,GM Chanteur, E. Richer, S. Hess, MMancini,RAllioux)
 !!
 !! FUNCTION
 !!  Contains Dipole moment calculation for Mercury
 !!
 !! NOTE
 !! The magnetic field input at initialization is derived from Andersson et al, Science,2008, from Mariner 10 and Messenger
 !! update version  : EGU abstract 2010, Alexeev et al, M = 196 nT*R_M^3, offse
!!of 405km Northward., tilt 4° (offset not implemented yet)
 !!
 !!
 subroutine add_b_dipole_earth(Bfield,ncm,Spe,gstep,s_min_loc)
 use defs_arr3Dtype
 use atm_magnetic_fields
  integer, intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type),intent(inout) :: Spe


!   call add_dipole_generic(Bfield,ncm,Spe,gstep,s_min_loc,196.e-9,180.,0.)
!local
  integer :: ii,jj,kk
  real(dp) :: radius,rinclinaison,rphase
  real(dp),dimension(3) :: ss,moment_dip,b_dip,moment_dip_u
  real(dp) :: bme
 __WRT_DEBUG_IN("add_b_dipole_earth")

   bme=150.e-9 ! en Tesla * rayon planete au cube (moment dipolaire de la terre)
  ! bme=196e-9 ! en Tesla * rayon planete au cube (moment dipolaire de mercure)
  !--Initialisation
  !--Dipolar Moment in planetary units (Spe%ref%mag) mu0/4pi*M
   moment_dip_u(1) = 0.
   moment_dip_u(2) = 0.
   moment_dip_u(3) = -1.
   moment_dip = moment_dip_u*bme/Spe%ref%mag*(Spe%P%radius)**three
  !--Main loop
  !--Relative distance from the planet in m
  do kk = 1,ncm(3)-1
   ss(3) = (real((kk-1),dp)*gstep(3) + s_min_loc(3)-Spe%P%centr(3))    
   do jj = 1,ncm(2)-1
    ss(2) = (real((jj-1),dp)*gstep(2) + s_min_loc(2)-Spe%P%centr(2))
    do ii = 1,ncm(1)-1
     ss(1) = (real((ii-1),dp)*gstep(1) + s_min_loc(1)-Spe%P%centr(1))
     !--Distance from the center in radius of planet
     radius = sqrt(dot_product(ss,ss))
     if (radius >= .75*Spe%P%radius) then
      b_dip = (three*ss*dot_product(ss,moment_dip)/(radius*radius)-moment_dip)/radius**three
      Bfield%x(ii,jj,kk) = Bfield%x(ii,jj,kk) + b_dip(1)
      Bfield%y(ii,jj,kk) = Bfield%y(ii,jj,kk) + b_dip(2)
      Bfield%z(ii,jj,kk) = Bfield%z(ii,jj,kk) + b_dip(3)
     else
      Bfield%x(ii,jj,kk) = zero
      Bfield%y(ii,jj,kk) = zero
      Bfield%z(ii,jj,kk) = zero
     endif
    enddo
   enddo
  enddo
  __WRT_DEBUG_OUT("add_b_dipole_earth")
 end subroutine add_b_dipole_earth
end module env_earth
