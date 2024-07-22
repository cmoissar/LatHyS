!!=============================================================
!!=============================================================
!!module: m_titan_env
!! NAME
!!  m_titan_env (MMancini)
!!
!! FUNCTION
!!  Contains definition of type for the environment planet titan
!!  !!!!!!!!A completer
!! NOTE
module env_titan

 use defs_basis
 use defs_species
 use defs_parametre
 use defs_particletype
 use defs_mpitype
 use defs_atmospheretype
 use m_writeout,only       : wrt_debug
 use atm_sections_efficaces
 use mpi

#include "q-p_common.h"

 implicit none
 private
 !--Pointer towards the density_exo allocated in m_exosphere
 real(dp),public,pointer :: density_H2(:,:,:),  & 
          &                 density_CH4(:,:,:), &
          &                 density_N2(:,:,:),  &
          &                 density_Op(:,:,:),  &
          &                 density_Hp(:,:,:), &
          &                 density_N2p(:,:,:), &
          &                 density_CH4p(:,:,:), &
          &                 density_H2p(:,:,:), &
          &                 prod_H2(:,:,:) , &
          &                 prod_CH4(:,:,:), &
          &                 prod_N2(:,:,:) 
 public ::               &
      exosphere_titan,   &!--Compute exosphere densities for Titan
      alloc_titan,       &!--Allocate density arrays for Titan
      dealloc_titan,     & !--Allocate density arrays for Titan
      init_species_titan, &  !--Intialise the Titan
      photoproduction_titan  !-- compute Titan's photoproduction

contains
 !!#####################################################################

 !!=============================================================
 !!routine: m_exo_titan/alloc_exo_titan
 !!
 !! FUNCTION
 !!  associate density arrays for Titan
 !! IN 
 !! OUT
 !!   
 !! SIDE EFFECT
 !! 
 !! NOTE
 subroutine alloc_titan(density,prod_pp,atmosphere,ncm)
 integer,dimension(3),intent(in) :: ncm
  real(dp),intent(inout),allocatable,target :: density(:,:,:,:)
  real(dp),intent(inout),allocatable,target :: prod_pp(:,:,:,:)
 type(atmosphere_type),intent(inout),target :: atmosphere

 atmosphere%n_species=8 ! nombre d'espece neutres et ionisés
 atmosphere%n_spe_pp=3 !  nombre d'especes obtenues par photoproduction
 !atmosphere%n_spe_iono=0 !  nombre d'especes ionospherique
 atmosphere%n_pp = 3 ! nombre de reactions de photoproduction
 atmosphere%nb_lo = 37 ! nombre de longueur d'onde du spectre UV
 atmosphere%n_exc = 0 ! nombre de reaction d'echange de charge

 call allocate_atmosphere(ncm,atmosphere,density,prod_pp)
  density_Op => density(:,:,:,1)
                atmosphere%species(1)%name = "Op    "
                atmosphere%species(1)%mass = 1._dp
                atmosphere%species(1)%charge = 1._dp
                atmosphere%species(1)%opaque = .FALSE.
                atmosphere%species(1)%iono = .TRUE.
                atmosphere%species(1)%prod => prod_pp(:,:,:,1)

  density_Hp => density(:,:,:,2)
                atmosphere%species(2)%name = "Hp    "
                atmosphere%species(2)%mass = 1._dp/16._dp
                atmosphere%species(2)%charge = 1._dp
                atmosphere%species(2)%opaque = .FALSE.
                atmosphere%species(2)%iono = .TRUE.
                atmosphere%species(2)%prod => prod_pp(:,:,:,2)
  
  density_N2p => density(:,:,:,3)
                atmosphere%species(3)%name = "N2p   "
                atmosphere%species(3)%mass = 28._dp/16._dp
                atmosphere%species(3)%charge = 1._dp
                atmosphere%species(3)%opaque = .FALSE.
                atmosphere%species(3)%iono = .TRUE.
                atmosphere%species(3)%prod => prod_pp(:,:,:,3)

   density_CH4p => density(:,:,:,4)
                atmosphere%species(4)%name = "CH4p  "
                atmosphere%species(4)%mass = 1._dp
                atmosphere%species(4)%charge = 1._dp
                atmosphere%species(4)%opaque = .FALSE.
                atmosphere%species(4)%iono = .TRUE.
                atmosphere%species(4)%prod => prod_pp(:,:,:,4)

  density_H2p => density(:,:,:,5)
                atmosphere%species(5)%name = "H2p   "
                atmosphere%species(5)%mass = 2._dp/16._dp
                atmosphere%species(5)%charge = 1._dp
                atmosphere%species(5)%opaque = .FALSE.
                atmosphere%species(5)%iono = .TRUE.
                atmosphere%species(5)%prod => prod_pp(:,:,:,5)

  density_N2 => density(:,:,:,6)
                atmosphere%species(6)%name = "N2    "
                atmosphere%species(6)%mass = 28._dp/16._dp
                atmosphere%species(6)%charge = zero
                atmosphere%species(6)%opaque = .TRUE.
                atmosphere%species(6)%iono = .FALSE.
                !atmosphere%species(6)%prod => prod_pp(:,:,:,6)

  density_CH4 => density(:,:,:,7)
                atmosphere%species(7)%name = "CH4    "
                atmosphere%species(7)%mass = 1._dp
                atmosphere%species(7)%charge = zero
                atmosphere%species(7)%opaque = .TRUE.
                atmosphere%species(7)%iono = .FALSE.
                !atmosphere%species(7)%prod => prod_pp(:,:,:,7)

  density_H2 => density(:,:,:,8)
                atmosphere%species(8)%name = "H2    "
                atmosphere%species(8)%mass = 2._dp/16._dp
                atmosphere%species(8)%charge = 1._dp
                atmosphere%species(8)%opaque = .TRUE.
                atmosphere%species(8)%iono = .FALSE.
                !atmosphere%species(8)%prod => prod_pp(:,:,:,8)


! associate pointes to photoproduction array
prod_N2 => prod_pp(:,:,:,3)
prod_CH4 => prod_pp(:,:,:,4)
prod_H2 => prod_pp(:,:,:,5)

!H2->H2+
atmosphere%photo_reactions(1)%mother => atmosphere%species(8)
atmosphere%photo_reactions(1)%daughter => atmosphere%species(5)

!N2 -> N2+
atmosphere%photo_reactions(2)%mother => atmosphere%species(6)
atmosphere%photo_reactions(2)%daughter => atmosphere%species(3)

! CH4 -> CH4+
atmosphere%photo_reactions(3)%mother => atmosphere%species(7)
atmosphere%photo_reactions(3)%daughter => atmosphere%species(4)



 end subroutine alloc_titan

 !!=============================================================
 !!routine: m_exo_titan/dealloc_exo_titan
 !!
 !! FUNCTION
 !!  deallocate density arrays for Titan
 !! IN 
 !! OUT
 !!   
 !! SIDE EFFECT
 !! 
 !! NOTE
 subroutine dealloc_titan(dummy)
  integer,intent(in) :: dummy
  density_H2 => Null()
  density_CH4 => Null()
  density_N2 => Null()
  prod_H2 => Null()
  prod_CH4 => Null()
  prod_N2 => Null()
 end subroutine dealloc_titan

 !!=============================================================
 !!routine: m_species/init_species_titan
 !!
 !! FUNCTION
 !!  Initialize species for titan
 !! IN 
 !! OUT
 !! SIDE EFFECT
 !!  
 !!
 subroutine init_species_titan(species,s_centr)

  real(dp),intent(in) :: s_centr(3)
  type(species_type),intent(inout) :: species
  
  integer,parameter :: O=1,H=2
  integer :: is
  real(dp) :: tot_tmp,vitessethermique

  !--Titan contains two species of particles

  !--Allocation of species_type
  call alloc_species(2,species)

  !--Name of the planet
  species%planetname = "titan"
  
  !--Intialize Physical parameter used here
  !--Physical density of reference (m-3)
  species%ref%density = 0.3*1e6

  !--Ions inertial length (km)
  species%ref%c_omegapi = Sp_Lt*sqrt(epsilon0*pmasse_e2/species%ref%density)/1.e3

!--Magnetic Field (Tesla)
  species%ref%mag = 6.12*1e-9

  !--Alfven Speed
  species%ref%alfvenspeed = species%ref%mag/&
       &                   (sqrt(mu0*amu_pmass*species%ref%density)*1.e3)

  !--Inverse of gyrofrequency
  species%ref%inv_gyro = amu_pmass/(e_Cb*species%ref%mag)

  !--Assignation of Planet values
  species%P%centr  = s_centr
#ifndef HAVE_NO_PLANET  
  species%P%radius = 2575._dp/species%ref%c_omegapi
  species%P%r_exo  = species%P%radius+1430._dp/species%ref%c_omegapi
  species%P%r_iono = species%P%radius   !TO change!!!!!!!!
  species%P%r_lim  = species%P%radius+700._dp/species%ref%c_omegapi
  species%P%speed  = zero 
#endif

  !--Oxygen and Hydrogen
  species%S(:)%name = (/"O  ","H  "/)
  
  !--Number of particles for cells
  species%S(:)%ng = (/2, 2/) 

  !--Direct Speeds  (vitesses dirigees?)
  species%S(O)%vxs = 110._dp/species%ref%alfvenspeed !-- 
  !--species%S(Hen)%vxs = species%S(H)%vxs !--
  species%S(H)%vxs = Species%S(O)%vxs !--  
  species%S(:)%vys  = zero
  species%S(:)%vzs  = zero

 !--Percentage of Hen in saturn's magnetosphere: 
  species%S(O)%percent = .67_dp
  !--species%S(Hen)%percent = 0.1_dp
  species%S(H)%percent = .33_dp

  !--Ratio of charges and masses and Temperatures
  species%S(:)%rcharge = (/one ,one/)
  species%S(:)%rmass   = (/one,one/16._dp/)
  species%S(:)%rtemp   = (/one,16._dp/)

  !--Betas
  species%S(O)%betas = 3.114
  !species%S(Hen)%betas = 0 
  species%S(H)%betas = 0.111
  species%betae = 0.322

  !--Rapport des vitesses thermque entre parallel et perpendiculair  H+ et Hen++
  species%S(:)%rspeed = one

  !--Rapport des vitesse thermiques entre 
  species%S(:)%rvth = one
  species%S(:)%rmds = one

  !--Rapport charge sur masse
  species%S(:)%qms  = species%S(:)%rcharge/species%S(:)%rmass

  !--Thermal speed (parallel and perpendicular)
  vitessethermique = sqrt( three*species%S(O)%betas/(one+two*species%S(O)%rvth**two))

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

 end subroutine init_species_titan

 !!=============================================================
 !!=============================================================
 !!subroutine: m_exo_titan/exosphere_titan
 !! NAME
 !!  exosphere_titan (RModolo,MMancini)
 !!
 !! FUNCTION
 !!  Contains Exosphere calculatio for Titan
 !!
 !! NOTE
 !!  The densities results are in (cm-3) where the
 !!  during the calculation (for the altitude) the units are Km.
 !!  Written by R. Modolo 11/10/02 inspired from the subroutine of the 2D code

 subroutine exosphere_titan(Spe,ncm,&
      &                    gstep,s_min_loc,  &           
      &                    resistivity&
      &                    )

  integer, intent(in) :: ncm(3)
  !  real(dp),intent(in) :: dt
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  real(dp),intent(inout) :: resistivity(:,:,:)
  type(species_type),intent(inout) :: Spe

  !local
  integer :: ii,jj,kk
  real(dp) :: radius,altitude_km,r_planet_km,inv_r_exo_km
  real(dp) :: ss(3)
!!!!!!!! A CHANGER TOUS LES VALEURS !!!!!!!!!!!!!!!!!!!!!
  real(dp),parameter ::   d0_N2  = 7.2e11
  real(dp),parameter ::   d1_N2  = 2.6e8
  real(dp),parameter ::   d0_CH4 = 7.2e8
  real(dp),parameter ::   d1_CH4 = 7e6 
  real(dp),parameter ::   d0_H2  = 2e9
  real(dp),parameter ::   d1_H2  = 6e5

  real(dp),parameter ::   e0_N2  = one/34.968              
  real(dp),parameter ::   e1_N2  = one/32.769
  real(dp),parameter ::   e0_CH4 = one/41.951
  real(dp),parameter ::   e1_CH4 = one/57.346
  real(dp),parameter ::   e0_H2  = one/338.690
  real(dp),parameter ::   e1_H2  = one/458.771

  real(dp),parameter ::   z0_N2  = 700._dp
  real(dp),parameter ::   z1_N2  = 1265._dp
  real(dp),parameter ::   z0_CH4 = 925._dp
  real(dp),parameter ::   z1_CH4 = 1265._dp
  real(dp),parameter ::   z0_H2  = 700._dp
  real(dp),parameter ::   z1_H2  = 1400._dp
  
  __WRT_DEBUG_IN("exosphere_titan")
  
  !--Initialisation
  ! if (Spe%P%centr(1) >= 30.0_dp*gstep(1)) Spe%P%centr(1) = Spe%P%centr(1) + Spe%P%speed*dt
  ! if (Spe%P%centr(1) <= 30.0_dp*gstep(1)) Spe%P%speed = zero
  
  !--radius of the Planet in Km
  r_planet_km = Spe%P%radius*Spe%ref%c_omegapi
  inv_r_exo_km = one/(Spe%P%r_exo*Spe%ref%c_omegapi)

  !--Main loop
  do kk = 1,ncm(3)-1
   ss(3) = real((kk-1),dp)*gstep(3) + s_min_loc(3) 
   do jj = 1,ncm(2)-1
    ss(2) = real((jj-1),dp)*gstep(2) + s_min_loc(2)
    do ii = 1,ncm(1)-1
     ss(1) = real((ii-1),dp)*gstep(1) + s_min_loc(1)      
     radius = sqrt(dot_product(ss-Spe%P%centr,ss-Spe%P%centr))
     if (radius >= Spe%P%r_lim) then
      !--Altitude in Kilometers
      altitude_km = (radius-Spe%P%radius)*(Spe%ref%c_omegapi)
      
      !--Pour H2
      density_H2(ii,jj,kk) = d0_H2*exp(-(altitude_km-z0_H2)*e0_H2)+&
           &                 d1_H2*exp(-(altitude_km-z1_H2)*e1_H2)
      !--Pour CH4
      density_CH4(ii,jj,kk) = d0_CH4*exp(-(altitude_km-z0_CH4)*e0_CH4)+&
           &                  d1_CH4*exp(-(altitude_km-z1_CH4)*e1_CH4) 
                
      !--For N2
      density_N2(ii,jj,kk) = d0_N2*exp(-(altitude_km-z0_N2)*e0_N2)+&
           &                 d1_N2*exp(-(altitude_km-z1_N2)*e1_N2)

      !if(kk==50 .and. jj==50) print *,"Al_km",Altitude_km,radius,ii,density_O(ii,jj,kk)
      !!print *,((s_cen(1)-nn)*gstep(1)-Spe%P%radius)*Spe%ref%c_omegapi,density_O(nn,int(s_cen(2)*gstep(2))+1,int(s_cen(3)*gstep(3))+1)
     else
      density_H2(ii,jj,kk) = zero
      density_CH4(ii,jj,kk) = zero
      density_N2(ii,jj,kk) = zero
     endif
    enddo
   enddo
  enddo
  print *,' dans exosphere',density_H2(:,:,:) 


  __WRT_DEBUG_OUT("exosphere_titan")
 end subroutine exosphere_titan
 
!****************************** END EXOSPHERE ************************************

 !!=============================================================
 !!routine: env_titan/photoproduction_ganymede
 !!
  subroutine photoproduction_titan(Spe,ncm,gstep,s_min_loc,atmosphere)
  use atm_photoproduction
  integer,intent(in) :: ncm(3)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type), intent(in) :: Spe
  type(atmosphere_type),intent(inout) ::atmosphere
  real(dp) :: F107,F107a,dist_conv
  
  F107  = 74._dp  ! daily F10.7 (e.g. 74)
  F107A = 74._dp  ! 81 day average F10.7 (F10.7A) (e.g. 86)
  dist_conv=9.1_dp**2  ! Flux solaire a l'orbite de l'objet--(1.e4 conversion to photons/(m2*s))

  call flux_solaire_generic(atmosphere%nb_lo,atmosphere,F107,F107a,dist_conv)
  call section_efficace_H2_abs(atmosphere%nb_lo,atmosphere%species(8)%ion_abs)
  call section_efficace_N2_abs(atmosphere%nb_lo,atmosphere%species(6)%ion_abs)
  call section_efficace_CH4_abs(atmosphere%nb_lo,atmosphere%species(7)%ion_abs)
  
  call section_efficace_H2_ion(atmosphere%nb_lo,atmosphere%photo_reactions(1)%ion_react)
  call section_efficace_N2_ion(atmosphere%nb_lo,atmosphere%photo_reactions(2)%ion_react)
  call section_efficace_CH4_ion(atmosphere%nb_lo,atmosphere%photo_reactions(3)%ion_react)
 
                call photoproduction_generic(Spe,ncm,gstep,s_min_loc,atmosphere)! does all the work
 end subroutine photoproduction_titan
 !****************************** END PHOTOPRODUCTION ************************************

end module env_titan
