!!=============================================================
!!=============================================================
!!module: env_mercure
!! NAME
!!  env_mercure (SHess)
!!
!! Gathers all the mercure related file into a single one, for more clarity
!!
!! FUNCTION
!!  Contains definition of type for the environment planet mercure
!!
!! NOTE
module env_mercure

 use defs_basis
 use defs_species
 use defs_atmospheretype
 use defs_parametre
 use m_writeout

#include "q-p_common.h"

 implicit none
 private

 !--Pointer towards the density_exo allocated in m_exosphere
 real(dp),pointer :: density_Hp(:,:,:) , &
                      density_H(:,:,:),  &
                      density_Nap(:,:,:), &
                      density_Na(:,:,:)
 !--Array for Photoproduction Rate
 
 real(dp),pointer :: prod_H(:,:,:), &
                     prod_Na(:,:,:)


 public::                     &
      alloc_mercure,             &!--allocate pp and density arrays for mercure
      dealloc_mercure,           &!--deallocate pp and density arrays for mercure
      init_species_mercure,      &!--Intialise the mercure environment
      create_ionosphere_mercure, &!--Compute ionosphere densities for mercure
          exosphere_mercure,         &!--Compute exosphere densities for Mars
          photoproduction_mercure,   &!--Compute production ratio
          feed_ionosphere_mercure,   &
          !charge_exchange_mercure,   &!--Charge exchange for Mercury
          add_b_int_mercure,      &
          preserve_b_int_mercure
contains
 !!#####################################################################
 !!=============================================================
 !!routine: env_mercure/init_species_mercure
 !!
 !! FUNCTION
 !!  Initialize species for mercure
 !! IN 
 !! planet center
 !!
 !! OUT
 !! species 
 !!
 subroutine init_species_mercure(species,s_centr)

  real(dp),intent(in) :: s_centr(3)
  type(species_type),intent(inout) :: species
  
  integer,parameter :: H=1,He=2
  integer :: is
  real(dp) :: tot_tmp,vitessethermique

  !--mercure contains two species of particles
  !--Allocation of species_type
  call alloc_species(2,species)
  
  !--Name of the planet
  species%planetname = "mercure"
  
  !--Intialize Physical parameter used here
  !--Physical density of reference (m-3)
  !-- Densite de l'espece dominante dans le plasma incident
  !-- Mercure aphelie:30.4e6
  !-- Mercure perihelie:69.35e6
  species%ref%density = 30.4e6  

  !--Ions inertial length (km)
  species%ref%c_omegapi = Sp_Lt*sqrt(epsilon0*pmasse_e2/species%ref%density)/1.e3 !Sp_Lt: speed of light

  !--Magnetic Field (Tesla)
  !-- Mercure aphelie =21nT
  !-- Mercure perihelie= 46 nt
  species%ref%mag = 14.3*1.e-9

  !--Alfven Speed
  species%ref%alfvenspeed = species%ref%mag/&
       &                   (sqrt(mu0*amu_pmass*species%ref%density)*1.e3)

  !--Inverse of gyrofrequency
  species%ref%inv_gyro = amu_pmass/(e_Cb*species%ref%mag)

  !--Max absorption length (km)
  species%ref%maxabs = 00._dp

  !--Assignation of Planet values
  species%P%centr  = s_centr
#ifndef HAVE_NO_PLANET
  species%P%radius = 2440._dp/species%ref%c_omegapi
  species%P%r_exo  = species%P%radius+0._dp/species%ref%c_omegapi
  species%P%r_iono = species%P%radius   !TO change!!!!!!!!
  species%P%r_lim  = species%P%radius+0._dp/species%ref%c_omegapi
  species%P%speed  = zero 
#endif

  !--Hydrogen and Helium
  species%S(:)%name = (/"H  ","He "/)
  
  !--Number of particles for cells
  species%S(:)%ng = (/5, 2/) !(/2,2/)

  !--Direct Speeds  (vitesses dirigees? en km/s)
  species%S(H)%vxs = 430.0_dp/species%ref%alfvenspeed !--  H+
  species%S(He)%vxs = species%S(H)%vxs !-- He++
  species%S(:)%vys  = zero
  species%S(:)%vzs  = zero

  !--Percentage of He in the solar wind: n(He++)/(n(He++)+n(H+))
  species%S(He)%percent = .05_dp
  species%S(H)%percent = one-species%S(He)%percent

  !--Ratio of charges and masses and Temperatures
  !-- First column always one
  species%S(:)%rcharge = (/one,two/)
  species%S(:)%rmass   = (/one,four/)
  species%S(:)%rtemp   = (/one,four/)

  !--Betas
        !-- Mercure Aphelie 0.32 (H) 0.06 (He) 0.47 (e)
        !-- Mercure perihelie  (H) 
  species%S(H)%betas = 0.32_dp
  species%S(He)%betas = 0.06_dp
  species%betae = 0.47_dp

  !--Rapport des vitesse thermiques entre H+ et He++
  species%S(:)%rspeed = one

  !--Rapport des vitesses thermque entre parallel et perpendiculair  H+ et He++
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

 end subroutine init_species_mercure !****************************** END INIT_SPECIES ************************************


 !!=============================================================
 !!routine: env_mercure/alloc_mercure
 !!
 !! FUNCTION
 !!  allocate photoproduction and exosphere arrays for mercure
 !! IN 
 !! density and photoproduction arrays, allocated in environment
 !!   
 subroutine alloc_mercure(density,prod_pp,atmosphere,ncm)
  integer,dimension(3),intent(in) :: ncm
  real(dp),intent(inout),allocatable,target :: density(:,:,:,:)
  real(dp),intent(inout),allocatable,target :: prod_pp(:,:,:,:)
  type(atmosphere_type),intent(inout),target :: atmosphere 

  atmosphere%n_species=4 !nombre d'especes atmospheriques (neutres et ions)
  atmosphere%n_spe_pp=2  !nombre d'especes obtenues par photoproduction
  !atmosphere%n_spe_iono = 0 ! nombre d'especes ionospherique
  atmosphere%n_pp=2      !nombre de reactions de photoproduction
  atmosphere%n_pp_freq_fixed = 1 ! nombre de reactions de photoproduction avec une frequence d'ionisation fixé
  atmosphere%nb_lo=37    !nombre de longueur d'onde du spectre UV
  atmosphere%n_exc=1     !nombre de reactions d'echange de charge
  atmosphere%n_ei=0      !impact electronique
  call allocate_atmosphere(ncm,atmosphere,density,prod_pp)
  !--Associate pointers to exospheric density and ionospheric arrays

  density_Hp   => density(:,:,:,1)
                atmosphere%species(1)%name  ="Hp       "
                atmosphere%species(1)%mass  = 1._dp
                atmosphere%species(1)%charge= one
                atmosphere%species(1)%opaque= .FALSE.
                atmosphere%species(1)%iono  = .TRUE.
                atmosphere%species(1)%prod=> prod_pp(:,:,:,1)
   density_Nap   => density(:,:,:,2)
                atmosphere%species(2)%name  ="Nap      "
                atmosphere%species(2)%mass  = 23._dp
                atmosphere%species(2)%charge= one
                atmosphere%species(2)%opaque= .FALSE.
                atmosphere%species(2)%iono  = .TRUE.
                atmosphere%species(2)%prod=> prod_pp(:,:,:,2)
   density_H   => density(:,:,:,3)
                atmosphere%species(3)%name  ="H         "
                atmosphere%species(3)%mass  = 1._dp
                atmosphere%species(3)%charge= zero
                atmosphere%species(3)%opaque= .FALSE.
                atmosphere%species(3)%iono  = .FALSE.
                
   density_Na   => density(:,:,:,4)
                atmosphere%species(4)%name  ="Na        "
                atmosphere%species(4)%mass  = 23._dp
                atmosphere%species(4)%charge= zero
                atmosphere%species(4)%opaque= .FALSE.
                atmosphere%species(4)%iono  = .FALSE.
 
 !--Associate pointers to photoproduction array 
  prod_H   => prod_pp(:,:,:,1)
  prod_Na   => prod_pp(:,:,:,2)
 
  !--Associate pointers for photoproduction 

  !H->H+
  atmosphere%photo_reactions(1)%mother     =>atmosphere%species(3)
  atmosphere%photo_reactions(1)%daughter   =>atmosphere%species(1)
  
  !Na->Na+
  atmosphere%photo_reactions(2)%mother     =>atmosphere%species(4)
  atmosphere%photo_reactions(2)%daughter   =>atmosphere%species(2)
  atmosphere%photo_reactions(2)%frequency  = 1.62e-5/0.22  
  
  !--Associate pointers for charge exchange
  !H + H+ ->H+ + H
  atmosphere%exc_reactions(1)%qsm        = one
  atmosphere%exc_reactions(1)%ion        =>atmosphere%species(1)
  atmosphere%exc_reactions(1)%neutral    =>atmosphere%species(3)
  atmosphere%exc_reactions(1)%cross_section = 2.E-19

 end subroutine alloc_mercure

 !!=============================================================
 !!routine: env_mercure/dealloc_mercure
 !!
 !! FUNCTION
 !!  deallocate photoproduction and exosphere arrays for mercure
 !!   
 subroutine dealloc_mercure(dummy)
  integer,intent(in) :: dummy
  density_Hp => Null()
  density_H =>  Null()
  prod_H => Null()
 end subroutine dealloc_mercure

 !!=============================================================
 !!routine: env_mercure/create_ionosphere_mercure
 !!
 !! FUNCTION
 !!  generate an ionosphere for mercure
 !!   
subroutine create_ionosphere_mercure(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
  use defs_particletype
  use atm_ionosphere
  use defs_grid

  real(dp),intent(in) :: gstep(3),s_min_loc(3),s_max_loc(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(in) ::atmosphere
  
  integer,dimension(3)  ::min_i,max_i
  integer::i,j,k,cnt
  real(dp),dimension(3) ::s_cen
  real(dp) :: r_lim,vol_unit
  real(dp),parameter :: k4=1.64E-10,k5=9.6E-11,k6=1.1E-9,k7=7.38E-8 !coefficients des reactions
  integer  :: npcell=300 !nombre max de particules par ceellule et par espece
  character(len=100) :: msg,msg2,msg3,msg4
  real(dp) ::r_lim2,rp2,rb,radius
  real(dp),parameter :: density_iono_max = 100.
  real(dp) :: Hscale
     
  Hscale = 1000._dp/(Spe%ref%c_omegapi) !--Scale height of the inosphere=1000km
  r_lim=(Hscale+Spe%P%r_lim)!altitude maximum a remplir
  r_lim2=r_lim**2
  rp2=Spe%P%r_lim**2
  s_cen = Spe%P%centr-s_min_loc
  vol_unit=1E6/Spe%ref%density !normalisation densite

  min_i=max(int((s_cen-r_lim)/gstep)-1,1) !ou commence l'ionosphere, le -1 corrige INT pour des valeurs negatives
  max_i=min(int((s_cen+r_lim)/gstep),ncm-2) !ou finit l'ionosphere
  cnt=0

  if(all(min_i.lt.max_i)) then !y a t-il de l'ionosphere dans la boite?
   do k=min_i(3),max_i(3)+1
          do j=min_i(2),max_i(2)+1
              rb = (float(j)*gstep(2)-s_cen(2))**2+(float(k)*gstep(3)-s_cen(3))**2
                  do i=min_i(1),max_i(1)+1
                        radius = (float(i)*gstep(1)-s_cen(1))**2+rb !squared
                if ((radius < rp2).or.(radius > r_lim2)) CYCLE
        cnt=cnt+1
        density_Hp(i,j,k)= density_iono_max*exp(-(sqrt(radius)-sqrt(rp2))/Hscale)
                enddo
        enddo
   enddo
  npcell=min(npcell,int((size(particule)-nptot)/float(cnt))) !limitation pour ne pas depasser le nombre permis de particules
  density_Hp=density_Hp*vol_unit
  call create_ionosphere_generic(Spe,s_cen,s_min_loc,particule,gstep,min_i,max_i,irand,nptot,npcell,atmosphere)!does all the work
  !call iono_densities_generic(particule,nptot,atmosphere,Spe,gstep,s_min_loc)
  endif
end subroutine create_ionosphere_mercure

!=============================================================
!=============================================================
 !!subroutine: env_mercure/exosphere_mercure
 !! NAME
 !!  exosphere_mercure (RModolo,MMancini,ERicher)
 !!
 !! FUNCTION
 !!  Contains Exosphere calculation for Mercury
 !!
 !! NOTE
 !!  The densities results are in (cm-3) where the
 !!  during the calculation (for the altitude) the units are Km.
 !!  Written by R. Modolo 11/10/02 inspired from the subroutine of the 2D code
 subroutine exosphere_mercure(Spe,ncm,gstep,s_min_loc,resistivity)
  integer, intent(in) :: ncm(3)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  real(dp),intent(inout) :: resistivity(:,:,:)
  type(species_type),intent(in) :: Spe

  integer :: ii,jj,kk
  real(dp) :: radius,altitude_km,r_planet_km,inv_r_exo_km
  real(dp) :: ss(3)
  !real(dp),parameter :: d0_O = 5.85e13,e0_O = one/10.56,&
       !&                d1_O = 7.02e9 ,e1_O = one/33.97,&
       !&                d2_O = 5.23e3 ,e2_O = one/626.2,&
       !&                d3_O = 9.76e2 ,e3_O = one/2790.,&
       !&                d4_O = 3.71e4 ,e4_O = one/88.47,&
       !&                d0_H = 1.5e5  ,e0_H = 25965.,   &
       !&                d1_H = 1.9e4  ,e1_H = 10365,    &
       !&                d0_CO2 = 1.e10,e0_CO2 = one/15.8
           
   real(dp),parameter :: d0_H = 1.e5  ,e0_H = 1292., & !d0_H sans ionisation   T=575K dawn/dusk  !,   &
                         d0_Na = 2.e6, e0_Na = 50., &  ! F. Leblanc, note de planéto 2.e4    
           &                 d1_H = 1.9e4  ,e1_H = 10365.

  __WRT_DEBUG_IN("exosphere_mercure")
  
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
      altitude_km = (radius-Spe%P%r_lim)*(Spe%ref%c_omegapi)
      
      !--Pour l'hydrogene
      !density_H(ii,jj,kk) = d0_H*exp(e0_H*(&
           !&   one/(altitude_km+r_planet_km)-inv_r_exo_km))  !+&
           !&   d1_H*exp(e1_H*(&
           !&   one/(altitude_km+r_planet_km)-inv_r_exo_km))
          density_H(ii,jj,kk) = d0_H*exp(-altitude_km/e0_H)
          
      !--Pour le Sodium
          density_Na(ii,jj,kk) = d0_Na*exp(-altitude_km/e0_Na)

      !--Pour l'Oxygene
      !density_O(ii,jj,kk) = d0_O*exp(-altitude_km*e0_O)+&
      !     &                d1_O*exp(-altitude_km*e1_O)+&
       !    &                d2_O*exp(-altitude_km*e2_O)+&
        !   &                d3_O*exp(-altitude_km*e3_O)+&
       !    &                d4_O*exp(-altitude_km*e4_O)
      !--For CO2
     ! density_CO2(ii,jj,kk) = d0_CO2*exp(-(altitude_km-140._dp)*e0_CO2)
     else
      density_H(ii,jj,kk) = zero
      density_Na(ii,jj,kk) = zero
      !density_O(ii,jj,kk) = zero
      !density_CO2(ii,jj,kk) = zero
     endif
    enddo
   enddo
  enddo
  __WRT_DEBUG_OUT("exosphere_mercure")
 end subroutine exosphere_mercure
 
  !!=============================================================
 !!routine: env_mercure/photoproduction_mercure
 !!
 subroutine photoproduction_mercure(Spe,ncm,gstep,s_min_loc,atmosphere)
  use atm_photoproduction
  use atm_sections_efficaces
  integer,intent(in) :: ncm(3)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type), intent(in) :: Spe
  type(atmosphere_type),intent(inout) ::atmosphere
  real(dp) :: F107,F107a,dist_conv

  F107  = 74._dp  ! daily F10.7 (e.g. 74)
  F107A = 74._dp  ! 81 day average F10.7 (F10.7A) (e.g. 86)
  dist_conv=0.22_dp  ! Flux solaire a l'orbite de l'objet--(1.e4 conversion to photons/(m2*s))
  
  call flux_solaire_generic(atmosphere%nb_lo,atmosphere,F107,F107a,dist_conv)
  call section_efficace_H_ion(atmosphere%nb_lo,atmosphere%photo_reactions(1)%ion_react)  
  call section_efficace_H_abs(atmosphere%nb_lo,atmosphere%species(1)%ion_abs)
  
  !stop
  call photoproduction_generic(Spe,ncm,gstep,s_min_loc,atmosphere)! does all the work
 end subroutine photoproduction_mercure

!!=============================================================

subroutine feed_ionosphere_mercure(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
  use defs_particletype
  use atm_ionosphere
  real(dp),intent(in) :: gstep(3),s_min_loc(3),s_max_loc(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(inout) ::atmosphere

  call Ion_production_generic(Spe,atmosphere,particule,20.,0.01,nptot,irand,gstep,s_min_loc,dummy)

end subroutine feed_ionosphere_mercure

 !!=============================================================
 !!subroutine: b_dipole/add_b_dipole_mercure
 !! NAME
 !!  add_b_dipole_mercure (RModolo,GM Chanteur, E. Richer, S. Hess, MMancini,RAllioux)
 !!
 !! FUNCTION
 !!  Contains Dipole moment calculation for Mercury
 !!
 !! NOTE
 !! The magnetic field input at initialization is derived from Andersson et al, 
!! Science,2008, from Mariner 10 and Messenger
 !! update version  : EGU abstract 2010, Alexeev et al, M = 196 nT*R_M^3, offset 
!!of 405km Northward., tilt 4° (offset not implemented yet)
 !!
 !!
 subroutine add_b_int_mercure(Bfield,ncm,Spe,gstep,s_min_loc)
 use defs_arr3Dtype
 use atm_magnetic_fields
  integer, intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type),intent(inout) :: Spe

  __WRT_DEBUG_IN("add_b_dipole_mercure")

  Spe%P%dipole(1) = Spe%P%centr(1)+0.
  Spe%P%dipole(2) = Spe%P%centr(2)+0.
  Spe%P%dipole(3) = Spe%P%centr(3)+0.
  !Spe%P%dipole(3) = Spe%P%centr(3)+484./Spe%ref%c_omegapi
  !print*,'-------dipole------',Spe%P%dipole

   call add_dipole_generic(Bfield,ncm,Spe,gstep,s_min_loc,196.e-9,180.,0.)
   !call add_dipquad(Bfield,ncm,Spe,gstep,s_min_loc,196.e-9,180.,0.,0.72,0.38)

  __WRT_DEBUG_OUT("add_b_dipole_mercure")
 end subroutine add_b_int_mercure
 
   !!=============================================================
 !!subroutine: env_mercure/preserve_dip_mercure
 !! NAME
 !!  preserve_dip_mercure (E. Richer)
 !!
 !!
 !!
 subroutine preserve_b_int_mercure(Bfield,ncm,Spe,gstep,s_min_loc)
 use defs_arr3Dtype
 use atm_magnetic_fields
 
  integer, intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type),intent(in) :: Spe

  __WRT_DEBUG_IN("preserve_b_int_mercure")

   !call preserve_dip_generic(Bfield,ncm,Spe,gstep,s_min_loc,196.e-9,180.,0.,75.)
   call preserve_dipquad(Bfield,ncm,Spe,gstep,s_min_loc,196.e-9,180.,0.,75.,.72,0.38)

  __WRT_DEBUG_OUT("preserve_b_int_mercure")
 end subroutine preserve_b_int_mercure
 
end module env_mercure
