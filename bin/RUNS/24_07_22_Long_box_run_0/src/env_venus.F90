!!=============================================================
!!=============================================================
!!module: env_venus
!! NAME
!!  env_venus (SAizawa)
!!
!! Gathers all the mars related file into a single one, for more clarity
!!
!! FUNCTION
!!  Contains definition of type for the environment planet venus
!!
!! NOTE
module env_venus

 use defs_basis
 use defs_species
 use defs_atmospheretype
 use atm_external_atmosphere
! use defs_parametre
 use m_writeout
#ifdef HAVE_NETCDF 
    use defs_basic_cdf
    use netcdf
#endif


#include "q-p_common.h"

 implicit none
 private

 integer,parameter :: CO2p=1,Op=2,Hp=3,O2p=4,O=5,CO2=6,Hn=7

 !--Pointer towards the density_exo allocated in m_exosphere
 real(dp),pointer ::        density_H(:,:,:),  &
      &                     density_Hp(:,:,:), &
      &                     density_O(:,:,:),  &
      &                     density_CO2(:,:,:),&
      &                     density_Op(:,:,:), &
      &                     density_O2p(:,:,:),&
      &                     density_CO2p(:,:,:)
 
 !--Array for Photoproduction Rate
 real(dp),pointer :: prod_O(:,:,:),prod_H(:,:,:),prod_CO2(:,:,:)

 public::                     &
      alloc_venus,             &!--allocate pp and density arrays for Venus
      dealloc_venus,           &!--deallocate pp and density arrays for Venus
      init_species_venus,      &!--Intialise the Venus environment
      exosphere_venus,         &!--Compute exosphere densities for Venus
      photoproduction_venus,   &!--Compute production ratio
      charge_exchange_venus,   &!--Charge exchange for Venus
      create_ionosphere_venus, &!--Compute ionosphere densities for Venus
      load_ionosphere_venus,   &!--Load ionsopheric densities for Venus
      feed_ionosphere_venus

contains
 !!#####################################################################
 !!=============================================================
 !!routine: env_venus/init_species_venus
 !!
 !! FUNCTION
 !!  Initialize species for venus
 !! IN 
 !! planet center
 !!
 !! OUT
 !! species 
 !!
 subroutine init_species_venus(species,s_centr)
 use defs_parametre

  real(dp),intent(in) :: s_centr(3)
  type(species_type),intent(inout) :: species
  
  integer,parameter :: H=1,He=2
  integer :: is
  real(dp) :: tot_tmp,vitessethermique

  !--Venus contains two species of particles
  !--Allocation of species_type
  call alloc_species(2,species)
  
  !--Name of the planet
  species%planetname = "venus"
  species%exospherename = "venus"
  
  !--Intialize Physical parameter used here
  !--Physical density of reference (m-3)
  species%ref%density = 1.47e7  ! 4.e6 ! extreme event 20.e6

  !--Ions inertial length (km)
  species%ref%c_omegapi = Sp_Lt*sqrt(epsilon0*pmasse_e2/species%ref%density)/1.e3

  !--Magnetic Field (Tesla)
  species%ref%mag = 7.58e-9 ! 3.e-9 ! extreme event 3.e-9

  !--Alfven Speed
  species%ref%alfvenspeed = species%ref%mag/&
       &                   (sqrt(mu0*amu_pmass*species%ref%density)*1.e3)

  !--Inverse of gyrofrequency
  species%ref%inv_gyro = amu_pmass/(e_Cb*species%ref%mag)

  !--Max absorption length (km)
  species%ref%maxabs = 500._dp

  !--Assignation of Planet values
  species%P%centr  = s_centr
#ifndef HAVE_NO_PLANET
  species%P%radius = 6052._dp/species%ref%c_omegapi
  species%P%r_exo  = species%P%radius+300._dp/species%ref%c_omegapi
  species%P%r_lim  = species%P%radius+200._dp/species%ref%c_omegapi
  if (trim(ionospherename) /="") then
   ! species%P%r_iono = species%P%radius+300._dp/species%ref%c_omegapi
    species%P%r_iono = species%P%radius+400._dp/species%ref%c_omegapi
  else  
    species%P%r_iono = species%P%r_lim+350._dp/species%ref%c_omegapi
  endif  
  species%P%speed  = zero 
#endif

  !--Hydrogen and Helium
  species%S(:)%name = (/"H  ","He "/)
  
  !--Number of particles for cells
  species%S(:)%ng = (/3, 1/) !(/2,2/)

  !--Direct Speeds  (vitesses dirigees?)
  species%S(H)%vxs = 338._dp/species%ref%alfvenspeed !--  H+  ! 12.18
  species%S(He)%vxs = species%S(H)%vxs !-- He++
  species%S(:)%vys  = zero
  species%S(:)%vzs  = zero

  !--Percentage of He in the solar wind: n(He++)/(n(He++)+n(H+))
  species%S(He)%percent = .05_dp
  species%S(H)%percent = one-species%S(He)%percent

  !--Ratio of charges and masses and Temperatures
  species%S(:)%rcharge = (/one,two/)
  species%S(:)%rmass   = (/one,four/)
  species%S(:)%rtemp   = (/one,four/)

  !--Betas
  species%S(H)%betas = 1.1_dp   ! 1.57  ! Extreme event 0.12
  species%S(He)%betas = 0.03_dp  ! 0.03 !               0.024
  species%betae = 0.7_dp        ! 1.08 !               0.084

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
  species%S(:)%vth2 = species%S(:)%rvth * species%S(:)%vth1

  tot_tmp = sum(species%S(:)%rmass*species%S(:)%percent)

  !-- Ionosphere temperature
  species%tempe_ratio = 0.01_dp !0.2eV
  !-- plasma viscosity
  species%viscosity=1.7e-9*species%ref%inv_gyro
  !--Macro-particle mass
  species%S(:)%sm =  species%S(:)%rmass*species%S(:)%percent/(tot_tmp*real(species%S(:)%ng,dp))

  !--Macro-particle charge 
  species%S(:)%sq = species%S(:)%qms*species%S(:)%sm

  !--Set probability of extraction
  tot_tmp = sum(real(species%S(:)%ng,dp)*species%S(:)%vxs)
  species%S(:)%prob = (real(species%S(:)%ng,dp)*species%S(:)%vxs)/tot_tmp

  !--Accumulate sum
  species%S(:)%prob = (/(sum(species%S(:is)%prob),is=1,species%ns)/)

  !-- Normalization factors
   !-- Temperature
  species%ref%conv_t=(species%ref%mag**2)/(species%ref%density*2._dp*mu0*kb_JK)


 end subroutine init_species_venus

 !!=============================================================
 !!routine: env_venus/alloc_venus
 !!
 !! FUNCTION
 !!  allocate photoproduction and exosphere arrays for Venus
 !! IN 
 !! density and photoproduction arrays, allocated in environment
 !!   
 subroutine alloc_venus(density,prod_pp,atmosphere,ncm)
 use defs_parametre

  integer,dimension(3),intent(in) :: ncm
  real(dp),intent(inout),allocatable,target :: density(:,:,:,:)
  real(dp),intent(inout),allocatable,target :: prod_pp(:,:,:,:)
  type(atmosphere_type),intent(inout),target :: atmosphere 

  atmosphere%n_species=7 !nombre d'especes atmospherique (neutres et ions)
  atmosphere%n_spe_pp=3  !nombre d'especes obtenues par photoproduction
  !atmosphere%n_spe_iono = 4 ! nombre d'especes ionospherique
  atmosphere%n_pp=4      !nombre de reactions de photoproduction
  atmosphere%nb_lo=37    !nombre de longueur d'onde du spectre UV
  atmosphere%n_exc=3     !nombre de reactions d'echange de charge
  atmosphere%n_ei=2      !nombre de reactions d'ionisation par impact electronique
  atmosphere%n_pp_freq_fixed=0
 call allocate_atmosphere(ncm,atmosphere,density,prod_pp)

  !--Associate pointers to exospheric density and ionospheric arrays
  density_CO2p => density(:,:,:,CO2p)
                atmosphere%species(CO2p)%name  = "CO2+      "
                atmosphere%species(CO2p)%mass  = 44._dp
                atmosphere%species(CO2p)%charge= one
                atmosphere%species(CO2p)%opaque= .FALSE.
                atmosphere%species(CO2p)%iono  = .TRUE.
                atmosphere%species(CO2p)%prod  => prod_pp(:,:,:,CO2p)
  
  density_Op   => density(:,:,:,Op)
                atmosphere%species(Op)%name  = "O+        "
                atmosphere%species(Op)%mass  = 16._dp
                atmosphere%species(Op)%charge= one
                atmosphere%species(Op)%opaque= .FALSE.
                atmosphere%species(Op)%iono  = .TRUE.
                atmosphere%species(Op)%prod  => prod_pp(:,:,:,Op)
  
  density_Hp    => density(:,:,:,Hp)
                atmosphere%species(Hp)%name  = "H+        "
                atmosphere%species(Hp)%mass  = 1._dp
                atmosphere%species(Hp)%charge= one
                atmosphere%species(Hp)%opaque= .FALSE.
                atmosphere%species(Hp)%iono  = .TRUE.
                atmosphere%species(Hp)%prod  => prod_pp(:,:,:,Hp)
  
  density_O2p  => density(:,:,:,O2p)
                atmosphere%species(O2p)%name  = "O2+       "
                atmosphere%species(O2p)%mass  = 32._dp
                atmosphere%species(O2p)%charge= one
                atmosphere%species(O2p)%opaque= .FALSE.
                atmosphere%species(O2p)%iono  = .TRUE.

  density_O    => density(:,:,:,O)
                atmosphere%species(O)%name  = "O         "
                atmosphere%species(O)%mass  = 16._dp
                atmosphere%species(O)%charge= zero
                atmosphere%species(O)%opaque= .TRUE.
                atmosphere%species(O)%iono  = .FALSE.
                write(atmosphere%species(O)%description,'(a)') &
& '<Profile>Krasnopolsky et al, 2002,Kim et al, 1998</Profile>'//&
& char(10)//'<ModelURL>http://dx.doi.org/10.1029/2001JE001809</ModelURL>'
  density_CO2  => density(:,:,:,CO2)
                atmosphere%species(CO2)%name  ="CO2        "
                atmosphere%species(CO2)%mass  = 44._dp
                atmosphere%species(CO2)%charge= zero
                atmosphere%species(CO2)%opaque= .TRUE.
                atmosphere%species(CO2)%iono  = .FALSE.
        write(atmosphere%species(CO2)%description,'(a)') &
& '<Profile>Ma et al, 2004</Profile>'
  density_H    => density(:,:,:,Hn)
                atmosphere%species(Hn)%name  = "H         "
                atmosphere%species(Hn)%mass  = 1._dp
                atmosphere%species(Hn)%charge= zero
                atmosphere%species(Hn)%opaque= .TRUE.
                atmosphere%species(Hn)%iono  = .FALSE.
                write(atmosphere%species(Hn)%description,'(a)') &
& '<Profile>Krasnopolsky et al, 1993</Profile>'
            
! ionospheric density at initialisation            
!  density_CO2i => density(:,:,:,CO2i)
!                atmosphere%species(CO2i)%name  = "CO2+i     "
!                atmosphere%species(CO2i)%mass  = 44._dp
!                atmosphere%species(CO2i)%charge= one
!                atmosphere%species(CO2i)%opaque= .FALSE.
!                atmosphere%species(CO2i)%iono  = .FALSE.
                  
!  density_Oi   => density(:,:,:,Oi)
!                atmosphere%species(Oi)%name  = "O+i       "
!                atmosphere%species(Oi)%mass  = 16._dp
!                atmosphere%species(Oi)%charge= one
!                atmosphere%species(Oi)%opaque= .FALSE.
!                atmosphere%species(Oi)%iono  = .FALSE.
!                  
!  density_Hi    => density(:,:,:,Hi)
!                atmosphere%species(Hi)%name  = "H+i       "
!                atmosphere%species(Hi)%mass  = 1._dp
!                atmosphere%species(Hi)%charge= one
!                atmosphere%species(Hi)%opaque= .FALSE.
!                atmosphere%species(Hi)%iono  = .FALSE.
                  
!  density_O2i  => density(:,:,:,O2i)
!  		atmosphere%species(O2i)%name  = "O2+i      "
!		atmosphere%species(O2i)%mass  = 32._dp
!		atmosphere%species(O2i)%charge= one
!		atmosphere%species(O2i)%opaque= .FALSE.
!		atmosphere%species(O2i)%iono  = .FALSE.
            
            

  !--Associate pointers to photoproduction array 
  prod_CO2 => prod_pp(:,:,:,CO2p)
  prod_O   => prod_pp(:,:,:,Op)
  prod_H   => prod_pp(:,:,:,Hp)
 
  !--Associate pointers for photoproduction 
  !CO2->CO2+
  atmosphere%photo_reactions(1)%mother     =>atmosphere%species(CO2)
  atmosphere%photo_reactions(1)%daughter   =>atmosphere%species(CO2p)
  !O->O+
  atmosphere%photo_reactions(2)%mother     =>atmosphere%species(O)
  atmosphere%photo_reactions(2)%daughter   =>atmosphere%species(Op)
  !H->H+
  atmosphere%photo_reactions(3)%mother     =>atmosphere%species(Hn)
  atmosphere%photo_reactions(3)%daughter   =>atmosphere%species(Hp)
  !CO2->O+
  atmosphere%photo_reactions(4)%mother     =>atmosphere%species(CO2)
  atmosphere%photo_reactions(4)%daughter   =>atmosphere%species(Op)

  !--Associate pointers for charge exchange
  !H + H+ ->H+ + H
  atmosphere%exc_reactions(1)%qsm        = one
  atmosphere%exc_reactions(1)%ion        =>atmosphere%species(Hp)
  atmosphere%exc_reactions(1)%neutral    =>atmosphere%species(Hn)
  atmosphere%exc_reactions(1)%cross_section = 2.5E-19
  !H+ + O->O+ + H
  atmosphere%exc_reactions(2)%qsm        = one
  atmosphere%exc_reactions(2)%ion        =>atmosphere%species(Hp)
  atmosphere%exc_reactions(2)%neutral    =>atmosphere%species(O)
  atmosphere%exc_reactions(2)%cross_section = 1.E-19!!1.e-15
  !H +O+ ->H+ + O
  atmosphere%exc_reactions(3)%qsm        = one/16._dp
  atmosphere%exc_reactions(3)%ion        =>atmosphere%species(Op)
  atmosphere%exc_reactions(3)%neutral    =>atmosphere%species(Hn)
  atmosphere%exc_reactions(3)%cross_section = 9.E-20

  !associate pointers for electron impact
  !H + e-> H+
  atmosphere%ei_reactions(1)%ion        =>atmosphere%species(Hp)
  atmosphere%ei_reactions(1)%neutral    =>atmosphere%species(Hn)
  atmosphere%ei_reactions(1)%n          = 5
  atmosphere%ei_reactions(1)%coeff(1:5) = &
& (/-1143.33,323.767,-35.0431,1.69073,-0.0306575/)
  write(atmosphere%ei_reactions(1)%description,'(A)') &
& '<ProcessModel>ISSIChallenge</ProcessModel>'

  !O + e-> O+
  atmosphere%ei_reactions(2)%ion        =>atmosphere%species(Op)
  atmosphere%ei_reactions(2)%neutral    =>atmosphere%species(O)
  atmosphere%ei_reactions(2)%n          = 5
  atmosphere%ei_reactions(2)%coeff(1:5) = &
& (/-1233.29,347.764,-37.4128,1.79337,-0.032277/)
  write(atmosphere%ei_reactions(2)%description,'(A)') &
& '<ProcessModel>ISSIChallenge</ProcessModel>'
 

 end subroutine alloc_venus

 !!=============================================================
 !!routine: env_venus/dealloc_venus
 !!
 !! FUNCTION
 !!  deallocate photoproduction and exosphere arrays for Venus
 !!   
 subroutine dealloc_venus(dummy)
  integer,intent(in) :: dummy
  density_O => Null()
  density_H => Null()
  density_Hp => Null()
  density_CO2 => Null()
  prod_O => Null()
  prod_H => Null()
  prod_CO2 => Null()
 end subroutine dealloc_venus

 !!=============================================================
 !!subroutine: env_venus/exosphere_venus
 !! NAME
 !!  exosphere_venus (SAizawa)
 !!
 !! FUNCTION
 !!  Contains Exosphere calculatio for Venus
 !!
 !! NOTE
 !!  The densities results are in (cm-3) where the
 !!  during the calculation (for the altitude) the units are Km.
 !!  Written by S.Aizawa
 subroutine exosphere_venus(Spe,ncm,gstep,s_min_loc,resistivity)
  use defs_variable,only : viscosity
  use defs_parametre,only :dt,exospherename,atmospherename
 
   integer, intent(in) :: ncm(3)
   real(dp),intent(in) :: gstep(3),s_min_loc(3)
   real(dp),intent(inout) :: resistivity(:,:,:)
   type(species_type),intent(inout) :: Spe
 
   integer :: ii,jj,kk,iii,jjj,kkk
   real(dp) :: radius,altitude_km,r_planet_km,inv_r_exo_km,nrm, inv_r_exo_km_O
   real(dp) :: ss(3)
   real(dp),parameter :: d0_O = 7.5e4  ,e0_O = 98389.7710,&
        &                d1_O = 7.02e9 ,e1_O = one/33.97,&
        &                d2_O = 5.23e3 ,e2_O = one/626.2,&
        &                d3_O = 9.76e2 ,e3_O = one/2790.,&
        &                d4_O = 3.71e4 ,e4_O = one/88.47,&
        &                d0_H = 1.32e5 ,e0_H = 138090.9067,   &    ! solar max
       ! &                d0_H = 1.5e5  ,e0_H = 25965.,   &    ! solar max
        &                d1_H = 0.  ,e1_H = 1.,    &   ! solar max
        !&                d0_H = 1.5e5  ,e0_H = 25965,    &    ! solar min
        !&                d1_H = 1.9e4  ,e1_H = 10365,    &    ! solar min
        &                d0_CO2 = 1.0e6,e0_CO2 = 1.0e8  , &
        &                d1_CO2 = 1.67e15, e1_CO2 = one/33.97
  character(len=500) :: msg          

  !!!! variables for SZA dependent exosphere !!!!
  real(dp) :: coef_G, coef_M, coef_mp, coef_kb, pi
  real(dp) :: radius_km, rho, x_km, rho_km, SZA_deg, xx, pwf
  real(dp) :: coef, coef_cold, coef_T_day, coef_beta_day, &
          &   coef_T_night, coef_beta_night, n_0, temp, beta2
 
   __WRT_DEBUG_IN("exosphere_venus")
   
   !--radius of the Planet in Km
   r_planet_km = Spe%P%radius*Spe%ref%c_omegapi

   inv_r_exo_km = one/((Spe%P%r_exo*Spe%ref%c_omegapi)-130.0_dp)
   inv_r_exo_km_O = one/((Spe%P%r_exo*Spe%ref%c_omegapi)-100.0_dp)

   !--initialize density arrays
       density_H(:,:,:) = zero
       density_O(:,:,:) = zero
       density_CO2(:,:,:) = zero

   !--
   nrm=1._dp!/729._dp

   !!!!! coefficient for exospheric profiles !!!!!
   coef_G = 6.67408e-11_dp
   coef_M = 4.867e24_dp
   coef_mp = 1.672e-27_dp
   coef_kb = 1.38e-23_dp
   pi = 4.0_dp*atan(1.0_dp)

   !--Main loop
   do kk = 1,ncm(3)-1
    ss(3) = real((kk-1),dp)*gstep(3) + s_min_loc(3)
    do jj = 1,ncm(2)-1
     ss(2) = real((jj-1),dp)*gstep(2) + s_min_loc(2)
     do ii = 1,ncm(1)-1
       ss(1) = real((ii-1),dp)*gstep(1) + s_min_loc(1)
      
       radius = sqrt(dot_product(ss-Spe%P%centr,ss-Spe%P%centr))
       radius_km = radius*(Spe%ref%c_omegapi)
       xx = ss(1) 
       x_km = (xx - Spe%P%centr(1)) * (Spe%ref%c_omegapi)
       rho_km = sqrt(radius_km**two - x_km**two)
       SZA_deg = 180.0_dp - atan2(rho_km, x_km) * 180.0_dp/pi

      if (radius >= Spe%P%radius) then
       !--Altitude in Kilometers
       altitude_km = (radius-Spe%P%radius)*(Spe%ref%c_omegapi)
       
       !--Pour l'hydrogene 
       if (x_km < 0.0_dp) then
        coef = SZA_deg/90.0_dp
        coef_cold = SZA_deg/180.0_dp
        coef_T_day = 285.0_dp
        coef_beta_day = coef_G*coef_M*coef_mp/(coef_kb*coef_T_day)*1e-3_dp
        coef_T_night = 110.0_dp
        coef_beta_night = coef_G*coef_M*coef_mp/(coef_kb*coef_T_night)*1e-3_dp

        density_H(ii, jj, kk) = density_H(ii, jj, kk) + nrm*(&
         & 1.e-6_dp*exp(-6.2625e-5_dp*(altitude_km+r_planet_km)+15.4817_dp+3.64e4_dp*one/(altitude_km+r_planet_km)) &
         & * (1.0_dp - coef) + coef *&
         & 1.e-6_dp*exp(-8.4607e-5_dp*(altitude_km+r_planet_km)+15.9944_dp +2.9743e4_dp*one/(altitude_km+r_planet_km)) &
         & + (1.0_dp - coef_cold) * &
         & 1.32e5_dp*exp(-coef_beta_day*(inv_r_exo_km-one/(altitude_km+r_planet_km))) &
         & + coef_cold * 2.59e7_dp * exp(-coef_beta_night*(inv_r_exo_km-one/(altitude_km+r_planet_km))))
         
       else
         coef = (180.0_dp - SZA_deg)/90.0_dp
         coef_cold = SZA_deg/180.0_dp
         coef_T_day = 285.0_dp
         coef_beta_day = coef_G*coef_M*coef_mp/(coef_kb*coef_T_day)*1e-3_dp
         coef_T_night = 110.0_dp
         coef_beta_night = coef_G*coef_M*coef_mp/(coef_kb*coef_T_night)*1e-3_dp

        density_H(ii, jj, kk) = density_H(ii, jj, kk) + nrm*(&
        & 1.e-6_dp*exp(-6.2309e-5_dp*(altitude_km+r_planet_km)+15.2723_dp+4.3781e4_dp*one/(altitude_km+r_planet_km)) &
        & * (1.0_dp - coef) + coef *&
        & 1.e-6_dp*exp(-8.4607e-5_dp*(altitude_km+r_planet_km)+15.9944_dp +2.9743e4_dp*one/(altitude_km+r_planet_km)) &
        & + (1.0_dp - coef_cold) * &
        & 1.32e5_dp*exp(-coef_beta_day*(inv_r_exo_km-one/(altitude_km+r_planet_km))) &
        & + coef_cold * 2.59e7_dp * exp(-coef_beta_night*(inv_r_exo_km-one/(altitude_km+r_planet_km))))
   
       endif

       !--Pour l'Oxygene
       coef_T_day = 6400.0_dp
       coef_beta_day = coef_G*coef_M*16.0_dp*coef_mp/(coef_kb*coef_T_day)*1e-3_dp 
       coef_T_night = 4847.0_dp
       coef_beta_night = coef_G*coef_M*16.0_dp*coef_mp/(coef_kb*coef_T_night)*1e-3_dp
       coef = SZA_deg/180.0_dp

       n_0 = 10.0_dp**(14.5_dp+0.46_dp * cos(SZA_deg*pi/180.0_dp) - 6.0_dp)
       if (SZA_deg*pi/180.0_dp < 1.9284_dp) then
          temp = (314.0_dp - 44.1_dp*(SZA_deg*pi/180.0_dp)**two)
       else
           temp = 150
       endif
       beta2 = 1.0e-3_dp*coef_G*coef_M*16.0_dp * coef_mp/(coef_kb*temp)

       if (coef < 0.4) then
         density_O(ii, jj, kk) = density_O(ii, jj, kk) + nrm*(&
             & (1-coef) * 7.5e4_dp * exp(-coef_beta_day*(1.0_dp/(6252.0_dp) - one/(altitude_km+r_planet_km))) &
             & + coef * 2.0e3_dp * exp(-coef_beta_night * (1.0_dp/(6352.0_dp) - one/(altitude_km+r_planet_km))) &
             & + n_0 * exp(-beta2*(1.0_dp/6222.0_dp - one/(altitude_km+r_planet_km))))
       else
         pwf = 5.0_dp * (coef - 0.4_dp)
         density_O(ii, jj, kk) = density_O(ii, jj, kk) + nrm*(&
             & (1.0_dp/10._dp**pwf) *(1.0_dp-coef) * 7.5e4_dp * exp(-coef_beta_day*&
             & (1.0_dp/(6252.0_dp) - one/(altitude_km+r_planet_km))) &
             & + coef * 2.0e3_dp * exp(-coef_beta_night * (1.0_dp/(6352.0_dp) - one/(altitude_km+r_planet_km))) &
             & + n_0 * exp(-beta2*(1.0_dp/6252.0_dp - one/(altitude_km+r_planet_km))))
       endif
        
       !--For CO2
       density_CO2(ii,jj,kk) = density_CO2(ii,jj,kk)+&
       &nrm*(1.0e6_dp*exp(-6.926639e6_dp*(1.0_dp/6252.0_dp-one/(altitude_km+r_planet_km)))&
       + 1.e11*exp(-9.620333e6_dp*(1.0_dp/6192.0_dp - one/(altitude_km+r_planet_km))))

      endif
        ! enddo
        ! enddo
        ! enddo
     enddo
    enddo
   enddo
 
   viscosity(:,:,:)=Spe%viscosity*density_O
   
  __WRT_DEBUG_OUT("exosphere_venus")
  end subroutine exosphere_venus


 !!=============================================================
 !!routine: env_venus/photoproduction_venus
 !!
 !! FUNCTION
 !!  computes photoproduction arrays for Venus
 !!   
 subroutine photoproduction_venus(Spe,ncm,gstep,s_min_loc,atmosphere)
  use atm_photoproduction
  use atm_sections_efficaces
  integer,intent(in) :: ncm(3)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type), intent(in) :: Spe
  type(atmosphere_type),intent(inout) ::atmosphere
  real(dp) :: F107,F107a,dist_conv

  atmosphere%F107     = 70.0_dp  ! daily F10.7 (e.g. 74)
  atmosphere%F107_Avg = 70.0_dp  ! 81 day average F10.7 (F10.7A) (e.g. 86)
  dist_conv=0.49_dp  ! Flux solaire a l'orbite de l'objet--(1.e4 conversion to photons/(m2*s))

  call flux_solaire_generic(atmosphere%nb_lo,atmosphere,atmosphere%F107,atmosphere%F107_Avg,dist_conv)
  call section_efficace_H_abs(atmosphere%nb_lo,atmosphere%species(Hn)%ion_abs)
  call section_efficace_O_abs(atmosphere%nb_lo,atmosphere%species(O)%ion_abs)
  call section_efficace_CO2_abs(atmosphere%nb_lo, atmosphere%species(CO2)%ion_abs)
  call section_efficace_CO2_ion(atmosphere%nb_lo,atmosphere%photo_reactions(1)%ion_react)
  call section_efficace_O_ion(atmosphere%nb_lo, atmosphere%photo_reactions(2)%ion_react)
  call section_efficace_H_ion(atmosphere%nb_lo,atmosphere%photo_reactions(3)%ion_react)
  call section_efficace_O_CO2_ion(atmosphere%nb_lo, atmosphere%photo_reactions(4)%ion_react)
  call photoproduction_generic(Spe,ncm,gstep,s_min_loc,atmosphere)! does all the work
 end subroutine photoproduction_venus

 !!=============================================================
 !!routine: env_venus/chrage_exchange_venus
 !!
 !! FUNCTION
 !!  computes charge exchange for Venus
 !!   
subroutine charge_exchange_venus(nn,kpickup,&
                 qsm,irand,ijk,v_p,ww,Spe,particule,atmosphere)
 use defs_parametre
 use atm_charge_exchange
 use defs_particletype
  integer,intent(in) :: nn,ijk(3)
  integer,intent(inout) :: irand,kpickup
  real(dp),intent(in) :: qsm,v_p(3),ww(8)
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(inout) ::atmosphere
  real(dp) ::Va,conv_fac,vmod
  conv_fac = Spe%ref%c_omegapi*1.e5 ! 1.38e+7
  vmod  = sqrt(sum(v_p*v_p))*dt*conv_fac
   call  charge_exchange_generic(nn,kpickup,qsm,irand,&
      &                     ijk,v_p,ww,Spe,particule,atmosphere)! does all the work
end subroutine

 !!=============================================================
 !!routine: env_venus/create_ionosphere_venus
 !!
 !! FUNCTION
 !!  generate an ionosphere for Venus
 !!   
subroutine create_ionosphere_venus(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
  !use defs_parametre
  use defs_particletype
  use atm_ionosphere
  use defs_grid,only : ncm

  real(dp),intent(in) :: gstep(3),s_min_loc(3),s_max_loc(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(in) ::atmosphere
  
  integer  ::min_i(3),max_i(3),i,j,k,cnt,npcell
  real(dp) ::s_cen(3),r_lim2,rp2,rb,radius,vol_unit,k4,k5,k6,k7,k8,r_lim
  character(len=100) :: msg



!======V environment dependent code from here V===============
 k4=1.64E-10;k5=9.6E-11;k6=9.4E-10;k7=1.6E-7*(300.0_dp/4000._dp)**0.55;k8=1.14E-4*1.0_dp/4000._dp
    !coefficients des reactions (pour k8 cf these de O wittasse)
 npcell=150 !nombre max de particules par cellule et par espece
 r_lim=500._dp !altitude maximum a remplir
!===================^  To here ^================================
  
  r_lim=(r_lim/Spe%ref%c_omegapi+Spe%P%radius)
  r_lim2=(r_lim+sqrt(gstep(1)**2+gstep(2)**2+gstep(3)**2))**2
  rp2=(Spe%P%r_lim-sqrt(gstep(1)**2+gstep(2)**2+gstep(3)**2))**2
  s_cen = Spe%P%centr-s_min_loc
  vol_unit=1E6/Spe%ref%density 
  min_i=max(int((s_cen-r_lim)/gstep)-1,1) !ou commence l'ionosphere, le -1 corrige INT pour des valeurs negatives
  max_i=min(int((s_cen+r_lim)/gstep),ncm-2) !ou finit l'ionosphere
  if(all(min_i.lt.max_i)) then !y a t-il de l'ionosphere dans la boite?
   

!======V environment dependent code from here V===============
  cnt=1; density_CO2p(:,:,:) = 0.;  density_O2p(:,:,:)  = 0.;  density_Op(:,:,:)   = 0.
  do k=min_i(3),max_i(3)+1
          do j=min_i(2),max_i(2)+1
              rb = (float(j-1)*gstep(2)-s_cen(2))**2+(float(k-1)*gstep(3)-s_cen(3))**2
                  do i=min_i(1),max_i(1)+1
                        radius = (float(i-1)*gstep(1)-s_cen(1))**2+rb
                if (((radius < rp2).or.(radius > r_lim2))) CYCLE !!.or.((rb < rp2).and.(float(i)*gstep(1) > s_cen(1)))

        if(density_O(i,j,k).ge.0) density_CO2p(i,j,k)=prod_CO2(i,j,k)/((k4+k5)*density_O(i,j,k))

        if(density_CO2(i,j,k).ne.0) density_Op(i,j,k)=&
        (prod_O(i,j,k)+k5*density_CO2p(i,j,k)*density_O(i,j,k))/&
        & (k6*density_CO2(i,j,k)+k8*density_H(i,j,k))
        if (radius.gt.(300./Spe%ref%c_omegapi+Spe%P%r_lim)**2) &
                & density_Op(i,j,k)=0.

        density_O2p(i,j,k) =(k4*density_CO2p(i,j,k)*density_O(i,j,k)+&
        & k6*density_CO2(i,j,k)*density_Op(i,j,k))/k7 
        density_O2p(i,j,k) =0.5*(sqrt((density_CO2p(i,j,k)+&
        & density_Op(i,j,k))**2+4.0_dp*density_O2p(i,j,k))-&
        & (density_CO2p(i,j,k)+density_Op(i,j,k)))
        
        cnt=cnt+1
                enddo
        enddo
 enddo
! normalisation
  density_CO2p=density_CO2p*vol_unit
  density_O2p=density_O2p*vol_unit
  density_Op=density_Op*vol_unit
!===================^  To here ^================================

 ! npcell=min(npcell,int(((size(particule)-nptot)/3.)/float(cnt)))
  write(msg,'(a,i10)')&
       & " npcell ",npcell
    call wrtout(6,msg,'PERS')
 
  call create_ionosphere_generic(Spe,s_cen,s_min_loc,particule,gstep,min_i,max_i,irand,nptot,npcell,atmosphere)!does all the work
  write(msg,'(a,i10)')&
       & " max weight ",int(maxval(particule(:)%char))
    call wrtout(6,msg,'PERS')
  endif
end subroutine create_ionosphere_venus

!!=============================================================
 !!routine: env_venus/load_ionosphere_venus
 !!
 !! FUNCTION
 !!  load an ionosphere for Venus
 !!   
subroutine load_ionosphere_venus(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
  use defs_parametre,only : ionospherename
  use defs_particletype
  use atm_ionosphere
  use defs_grid,only : ncm

  real(dp),intent(in) :: gstep(3),s_min_loc(3),s_max_loc(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(inout) ::atmosphere
  
  integer  ::min_i(3),max_i(3),npcell
  real(dp) ::s_cen(3),r_lim2,rp2,r_lim
  character(len=500) :: msg
  
  integer :: nAlti, ncid,StId,nTime,nLong,nLati,nHoru,varfieldId,nvar(3),varid
 integer :: nAltiKp1,i,ilat,ialt_min,ialt_max,ilon,iHoru,ir,ip,it,ialt,idone,ilatSub,ilonSub,idone0
 integer,dimension(nf90_max_var_dims) :: DimfieldId

 real(dp),allocatable :: Time(:), Latitude(:), Longitude(:) 
 real(dp),allocatable :: zls(:), dsm(:), dec(:),longsubsol(:)
 real(dp),allocatable :: sza(:,:) 
 real(dp),dimension(:,:,:),allocatable :: Altitude,nHp,nCO2p,nOp,nO2p,Tn, &
                                          Temp_elec,Temp_elec0,Altitude_tmp,nHp_tmp,nCO2p_tmp,&
                                          nO2p_tmp,nOp_tmp
 integer :: ii,jj,kk,tt
 integer :: ilat_GCM,ilon_GCM,ialt_GCM
 real(dp) :: x_cdr,y_cdr,z_cdr,r_cdr,ss(3),Xpc,Ypc,radius
 real(dp),save :: Obliqui = 0.439 ! Obliquitee of the planet
 real(dp) :: clock, szamin,cs0, cs, cmcs, coscmcs,cosclock
 real(dp) :: diff_Lat,diff_Long,diff_Alt,Lon_GEO,Lat_GEO,Alt_GEO_min,Alt_GEO_max 
 real(dp) :: a_CO2,b_CO2,log_ni_av,Alt_av,a_num,a_denom,a_O,b_O,a_O2,b_O2,a_H,b_H

 __WRT_DEBUG_IN("load_ionosphere_venus_LMD")
 
 call wrt_double(6,"Reading file : "//ionospherename,wrtscreen,wrtdisk)
 
 !-- Open NetCDF file
 stId = nf90_open(trim(ionospherename),nf90_nowrite, ncid)
 call test_cdf(stId)
 call get_simple_dimens_cdf(ncid,"Time",nTime)
 call get_simple_dimens_cdf(ncid,"Longitude",nLong)
 call get_simple_dimens_cdf(ncid,"Latitude",nLati)
 
 nHoru = nLong * nLati
 StId = nf90_inq_varId(ncid,"Altitude",varfieldId)
 StId = nf90_Inquire_Variable(ncid,varfieldId,dimids = DimfieldId)
 StId = nf90_Inquire_Dimension(ncid,DimfieldId(1), len=nvar(1))
 StId = nf90_Inquire_Dimension(ncid,DimfieldId(2), len=nvar(2))
 StId = nf90_Inquire_Dimension(ncid,DimfieldId(3), len=nvar(3))
 
 nAlti = nvar(2)
   write(msg,'(2a,4(a,a17,i12))')&
        & ch10," ___________________ Ionosphere file Information  _________________",&
        & ch10, "   nAlti  = ",nAlti,&
        & ch10, "   nLong  = ",nLong,&
        & ch10, "   nLati  = ",nLati,&
        & ch10, "   nTime  = ",nTime
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 allocate(Time(nTime))
 call get_simple_variable_cdf(ncid,"Time",Time(:))
   write(msg,'(2a,2(a,a17,e12.5))')&
        & ch10," ___________________ Ionosphere Time  _________________",&
        & ch10, "   Start  = ",minval(Time(:)),&
        & ch10, "   Stop  = ",maxval(Time(:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)       
 allocate(zls(nTime))
 call get_simple_variable_cdf(ncid,"zls",zls(:))
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ______________ Ionosphere Solar Longitude (Ls)  _____________",&
        & ch10, " Start (degrees)  = ",minval(zls(:)*180.0_dp/pi),&
        & ch10, " Stop             = ",maxval(zls(:)*180.0_dp/pi)
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 allocate(dsm(nTime))
 call get_simple_variable_cdf(ncid,"dsm",dsm(:))
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Ionosphere Sun distance  _________________",&
        & ch10, " Start (AU)  = ",dsm(1),&
        & ch10, " Stop        = ",dsm(nTime)
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 allocate(longsubsol(nTime))
  call get_simple_variable_cdf(ncid,"longsubsol",longsubsol(:))
  write(msg,'(2a,2(a,a20,e12.5))')&
         & ch10," ___________________ Ionosphere Subsolar longitude  _________",&
         & ch10, " Start (AU)  = ",longsubsol(1),&
         & ch10, " Stop        = ",longsubsol(nTime)
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 allocate(dec(nTime))
 call get_simple_variable_cdf(ncid,"dec",dec(:))
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Ionosphere Declination  _________________",&
        & ch10, " Start (AU)  = ",dec(1),&
        & ch10, " Stop        = ",dec(nTime)
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 allocate(sza(nHoru,nTime))
 call get_simple_variable_cdf(ncid,"mu0",sza(:,:))
 sza(1:nHoru,1:nTime) = acos(sza(1:nHoru,1:nTime))
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Ionosphere SZA  _________________",&
        & ch10, " Min (AU)  = ",minval(sza(:,:)),&
        & ch10, " Max       = ",maxval(sza(:,:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 allocate(Longitude(nLong))
 call get_simple_variable_cdf(ncid,"longitude",Longitude(:))
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Ionosphere Longitude  _________________",&
        & ch10, " Min (degrees)  = ",minval(Longitude(:)),&
        & ch10, " Max            = ",maxval(Longitude(:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 allocate(Latitude(nLati))
 call get_simple_variable_cdf(ncid,"Latitude",Latitude(:))
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Ionosphere Latitude  _________________",&
        & ch10, " Min (degrees)  = ",minval(Latitude(:)),&
        & ch10, " Max            = ",maxval(Latitude(:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 ! adding 15 points more in altitude
 allocate(Altitude(nHoru,nAlti+15,nTime))
 stId = nf90_inq_varid(ncid,'Altitude', varid)
 call test_cdf(stId)
 call get_simple_variable_cdf(ncid,"Altitude",Altitude(:,1:nAlti,:))
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Ionosphere Altitude  _________________",&
        & ch10, " Min (km)  = ",minval(Altitude(:,:,:)),&
        & ch10, " Max       = ",maxval(Altitude(:,:,:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
  ! adding 15 points more in altitude
 allocate(nCO2p(nHoru,nAlti+15,nTime))
 call get_simple_variable_cdf(ncid,"co2plus",nCO2p(:,1:nAlti,:))
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Ionosphere nCO2p  _________________",&
        & ch10, " Min (cm-3)  = ",minval(nCO2p(:,:,:)),&
        & ch10, " Max         = ",maxval(nCO2p(:,:,:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
  ! adding 15 points more in altitude
 allocate(nOp(nHoru,nAlti+15,nTime))
 call get_simple_variable_cdf(ncid,"oplus",nOp(:,1:nAlti,:))
 write(msg,'(2a,2(a,a20,e12.5))')&
        & ch10," ___________________ Ionosphere nOp  _________________",&
        & ch10, " Min (cm-3)  = ",minval(nOp(:,:,:)),&
        & ch10, " Max         = ",maxval(nOp(:,:,:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
  ! adding 15 points more in altitude
 allocate(nHp(nHoru,nAlti+15,nTime))
  call get_simple_variable_cdf(ncid,"hplus",nHp(:,1:nAlti,:))
  write(msg,'(2a,2(a,a20,e12.5))')&
         & ch10," ___________________ Ionosphere nHp  _________________",&
         & ch10, " Min (cm-3)  = ",minval(nHp(:,:,:)),&
         & ch10, " Max         = ",maxval(nHp(:,:,:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
  ! adding 15 points more in altitude
  allocate(nO2p(nHoru,nAlti+15,nTime))
  call get_simple_variable_cdf(ncid,"o2plus",nO2p(:,1:nAlti,:))
  write(msg,'(2a,2(a,a20,e12.5))')&
         & ch10," ___________________ Ionosphere nO2p  _________________",&
         & ch10, " Min (cm-3)  = ",minval(nO2p(:,:,:)),&
         & ch10, " Max         = ",maxval(nO2p(:,:,:))
 call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 
 Latitude(1:nLati)  = Latitude(1:nLati)*pi/180.0_dp
 Longitude(1:nLong) = Longitude(1:nLong)*pi/180.0_dp
 
 
 
 ! Extraplating the GCM to higher altitude to be sure that we have GCM input everywhere below 300km
 do ii=1,nHoru
 do kk = nAlti+1,nAlti+15
 Altitude(ii,kk,nTime) = Altitude(ii,nAlti,nTime)+(kk-nAlti)*10.
 enddo
 enddo
 !print *,'Extended Altitude'
 !print *,Altitude(1,nAlti:nAlti+15,nTime)
 
 ! Computing the linear regression and extrapolation of the density
 do ii=1,nHoru
! First part linear regression
   alt_av = 0.
   do kk=nAlti-4,nAlti
     alt_av = alt_av+Altitude(ii,kk,nTime)
   enddo
   alt_av = alt_av/5.
    
! for CO2+      
   log_ni_av = 0.
   a_CO2 = 0.
   b_CO2 = 0.
   a_num = 0.
   a_denom = 0.
   do kk=nAlti-4,nAlti
     log_ni_av = log_ni_av + log(nCO2p(ii,kk,nTime))
   enddo
   log_ni_av = log_ni_av/5.
   do kk=nAlti-4,nAlti
   a_num = a_num + (Altitude(ii,kk,nTime)-alt_av)*(log(nCO2p(ii,kk,nTime))-log_ni_av)
   a_denom = a_denom + (Altitude(ii,kk,nTime)-alt_av)**2
   enddo
   a_CO2 = a_num/a_denom
   b_CO2 = log_ni_av-a_CO2*alt_av
  ! print *, '1/a_CO2,b_CO2',1./a_CO2,b_CO2
  
! for O2+      
   log_ni_av = 0.
   a_O2 = 0.
   b_O2 = 0.
   a_num = 0.
   a_denom = 0.
   do kk=nAlti-4,nAlti
     log_ni_av = log_ni_av + log(nO2p(ii,kk,nTime))
   enddo
   log_ni_av = log_ni_av/5.
   do kk=nAlti-4,nAlti
   a_num = a_num + (Altitude(ii,kk,nTime)-alt_av)*(log(nO2p(ii,kk,nTime))-log_ni_av)
   a_denom = a_denom + (Altitude(ii,kk,nTime)-alt_av)**2
   enddo
   a_O2 = a_num/a_denom
   b_O2 = log_ni_av-a_O2*alt_av  
   
   
! for O+      
   log_ni_av = 0.
   a_O = 0.
   b_O = 0.
   a_num = 0.
   a_denom = 0.
   do kk=nAlti-4,nAlti
     log_ni_av = log_ni_av + log(nOp(ii,kk,nTime))
   enddo
   log_ni_av = log_ni_av/5.
   do kk=nAlti-4,nAlti
   a_num = a_num + (Altitude(ii,kk,nTime)-alt_av)*(log(nOp(ii,kk,nTime))-log_ni_av)
   a_denom = a_denom + (Altitude(ii,kk,nTime)-alt_av)**2
   enddo
   a_O = a_num/a_denom
   b_O = log_ni_av-a_O*alt_av 
   
! for H+      
   log_ni_av = 0.
   a_H = 0.
   b_H = 0.
   a_num = 0.
   a_denom = 0.
   do kk=nAlti-4,nAlti
     log_ni_av = log_ni_av + log(nHp(ii,kk,nTime))
   enddo
   log_ni_av = log_ni_av/5.
   do kk=nAlti-4,nAlti
   a_num = a_num + (Altitude(ii,kk,nTime)-alt_av)*(log(nHp(ii,kk,nTime))-log_ni_av)
   a_denom = a_denom + (Altitude(ii,kk,nTime)-alt_av)**2
   enddo
   a_H = a_num/a_denom
   b_H = log_ni_av-a_H*alt_av   
   
 ! second part extrapolation
   do kk=nAlti+1,nAlti+15
   nCO2p(ii,kk,nTime) = exp(a_CO2*Altitude(ii,kk,nTime)+b_CO2)
   nO2p(ii,kk,nTime) = exp(a_O2*Altitude(ii,kk,nTime)+b_O2)
   nOp(ii,kk,nTime) = exp(a_O*Altitude(ii,kk,nTime)+b_O)
 !  nCO2p(ii,kk,nTime) = nCO2p(ii,nAlti,nTime)
 !  nO2p(ii,kk,nTime) = nO2p(ii,nAlti,nTime)
 !  nOp(ii,kk,nTime) = nOp(ii,nAlti,nTime)

   nHp(ii,kk,nTime) = nHp(ii,nAlti,nTime)
  ! print *,'Add O2 profile ', 1./a_O2,nO2p(ii,nAlti-1:nAlti+10,nTime)
   enddo
 
 enddo
 ! now changing the number of altitude to match the size of the new density grid
 nAlti = nAlti+15

! Philisophy
! for each grid pont (MSO or simulation coord)
! we check if the altitude is lower than 500 km (upper boundary for ionospheric model)
! if no => nothing to do, the density is determined by the exospheric model
! if yes = > from the MSO position we convert it in GEO coord.
!		we find the closest grid point in the GCM input and affect it to the density

! first we determine the location of subsolar latitude and longitude of the GCM
! we also determoine angle cs
! we use only the last GCM time diagnostic (nTime)

cosclock = cos(Obliqui)/cos(dec(nTime))
if (abs(cosclock) <= 1.) then
  clock   = acos(cos(Obliqui)/cos(dec(nTime))) 
else
  if (cosclock > 1) clock = 0.
  if (cosclock < -1) clock = pi
endif



szamin = 3.0_dp*pi
cs0    = 2._dp*pi*(0.5_dp - Time(nTime))
if (cs0.gt. pi) cs0 = cs0 - 2._dp*pi
if (cs0.lt.-pi) cs0 = cs0 + 2._dp*pi
do ilat=1,nLati
  do ilon=1,nLong
    iHoru   = ilon + (ilat - 1)*nLong
    if (sza(iHoru,nTime).lt.szamin) then
      szamin  = sza(iHoru,nTime)
      ilatSub = ilat
      ilonSub = ilon
!   write(6,'(2(1x,i3),3(1x,e12.5))')ilatSub,ilonSub,szamin,Latitude(ilat),Longitude(ilon)
    endif
!  write(6,'(2(1x,i3),3(1x,e12.5))')ilat,ilon,sza(iHoru,iTime)*180.0/Pi,Latitude(ilat)*180.0/Pi,Longitude(ilon)*180.0/Pi
   if (cos(Latitude(ilat)).gt.1.E-06) then
       coscmcs = (cos(sza(iHoru,nTime)) - sin(Latitude(ilat))*sin(dec(nTime))) &
     & / (cos(Latitude(ilat))*cos(dec(nTime)))
!  write(6,'(5(1x,e12.5))')dcos(sza(iHoru,iTime)),dsin(Latitude(ilat)),dsin(dec(iTime)),dcos(Latitude(ilat)),dcos(dec(iTime))
       cmcs    = acos(coscmcs)
       if (Longitude(ilon)-Longitude(ilonSub).gt.pi) cmcs = 2.0*Pi - cmcs
       cs = Longitude(ilon) - cmcs
       if (cs.gt.pi)  cs = pi
       if (cs.lt.-pi) cs = -pi             
    endif 
  enddo
enddo

!print *,'cs,cs0,longsubsol(nTime)',cs,cs0,longsubsol(nTime)
cs0 = longsubsol(nTime)
print *,'clock,cs0,dec(nTime)',clock,cs0,dec(nTime)

 ! reset to zero previously charged neutral density
 ss(:) = zero
 density_Op(:,:,:) = zero
 density_CO2p(:,:,:) = zero
 density_O2p(:,:,:) = zero
 density_Hp(:,:,:) = zero
 
 
  do kk = 1,ncm(3)-1
    do jj = 1,ncm(2)-1
      do ii = 1,ncm(1)-1
            ss(3) = (real((kk-1),dp))*gstep(3) + s_min_loc(3) 
            ss(2) = (real((jj-1),dp))*gstep(2) + s_min_loc(2)
            ss(1) = (real((ii-1),dp))*gstep(1) + s_min_loc(1)      
 
            radius = sqrt(dot_product(ss-Spe%P%centr,ss-Spe%P%centr))
 
        !   if ((radius <= Spe%P%radius+500./Spe%ref%c_omegapi).and.(radius >= maxval(Altitude(:,1,1)))) then
           if ((radius <= Spe%P%radius+500./Spe%ref%c_omegapi).and.(radius > Spe%P%r_lim)) then
        !   if (radius <= Spe%P%radius+500./Spe%ref%c_omegapi) then           
              x_cdr = -(ss(1)-Spe%P%centr(1))
              y_cdr = -(ss(2)-Spe%P%centr(2))
              z_cdr =  (ss(3)-Spe%P%centr(3))
              r_cdr = sqrt(x_cdr**2+y_cdr**2+z_cdr**2)
              x_cdr = x_cdr/r_cdr
              y_cdr = y_cdr/r_cdr
              z_cdr = z_cdr/r_cdr
              ! x_cdr,y_cdr and z_cdr are the grid point position in MSO frame
              ! we determine the latitude and longitude in the GEO (GCM) frame
              Lat_GEO = asin(sin(dec(nTime))*x_cdr+sin(clock)*cos(dec(nTime))*y_cdr + &
                & cos(clock)*cos(dec(nTime))*z_cdr)
              Xpc = cos(dec(nTime))*cos(cs0)*x_cdr + &
                & (-sin(clock)*sin(dec(nTime))*cos(cs0)-cos(clock)*sin(cs0))*y_cdr + &
                & (-cos(clock)*sin(dec(nTime))*cos(cs0)+sin(clock)*sin(cs0))*z_cdr
             Ypc = cos(dec(nTime))*sin(cs0)*x_cdr + &
                & (-sin(clock)*sin(dec(nTime))*sin(cs0)+cos(clock)*cos(cs0))*y_cdr + &
                & (-cos(clock)*sin(dec(nTime))*sin(cs0)-sin(clock)*cos(cs0))*z_cdr
             Lon_GEO = atan(Ypc/Xpc)  
             if (Xpc.lt.0._dp) then
               Lon_GEO = Lon_GEO + pi
             else
               if (Ypc.lt.0._dp) Lon_GEO = Lon_GEO + 2._dp*pi
            endif
            if (Lon_GEO.gt.pi) Lon_GEO = Lon_GEO - 2._dp*pi
            if (abs(Lon_GEO).gt.pi) then
                write(6,'(3(1x,i3),a,4(1x,e12.5))')ii,jj,kk,' Lon_GEO = ',Lon_GEO,Xpc,Ypc,atan(Ypc/Xpc)
                stop
            endif
            ! we determine the index for the corresponding longitude in the GCM
            ilon = 1
            ilat = 1
            diff_Long = 3.*pi
            diff_Lat  = 3.*pi
            do ilon_GCM = 2,nLong
              if ((Longitude(ilon_GCM-1)<Lon_GEO) .and. (Longitude(ilon_GCM) >= Lon_GEO)) ilon = ilon_GCM
              !if (abs(Longitude(ilon_GCM)-Lon_GEO) < diff_Long) then
              !   diff_Long = abs(Longitude(ilon_GCM)-Lon_GEO)
              !   ilon = ilon_GCM
              ! endif
            enddo
            do ilat_GCM = 2,nLati
              
              if ((Latitude(ilat_GCM-1) > Lat_GEO) .and. (Latitude(ilat_GCM) <= Lat_GEO)) ilat = ilat_GCM
              !print *,'Lattitude GCM, LAT_GEO,ilat',Latitude(ilat_GCM),Lat_GEO,ilat
              !if (abs(Latitude(ilat_GCM)-Lat_GEO) < diff_Lat) then
              !	 diff_Lat = abs(Latitude(ilat_GCM)-Lat_GEO)
              !	 ilat = ilat_GCM
              ! endif
            enddo
            iHoru   = ilon + (ilat - 1)*nLong
            
!	    Alt_GEO_min = (radius-Spe%P%radius)*Spe%ref%c_omegapi
 !  	    ialt = 1
 !  	    diff_alt = 4000.
 !  	    do ialt_GCM = 1,nAlti
 !  	      if (abs(Altitude(iHoru,ialt_GCM,nTime)-Alt_GEO_min) < diff_alt) then
 !  	        diff_alt = abs(Altitude(iHoru,ialt_GCM,nTime)-Alt_GEO_min)
 !  	        ialt = ialt_GCM
 !  	      endif  
 !  	    enddo
            
 !  	      density_Op(ii,jj,kk) = nOp(iHoru,ialt,nTime)
 !  	      density_CO2p(ii,jj,kk) = nCO2p(iHoru,ialt,nTime)
 !  	      density_O2p(ii,jj,kk) = nO2p(iHoru,ialt,nTime)   	    
            
            Alt_GEO_min = (radius-Spe%P%radius-gstep(1)/2._dp)*Spe%ref%c_omegapi
            Alt_GEO_max = (radius-Spe%P%radius+gstep(1)/2._dp)*Spe%ref%c_omegapi
            ialt_min = 1
            diff_alt = 4000.
            do ialt_GCM = 1,nAlti
              if (abs(Altitude(iHoru,ialt_GCM,nTime)-Alt_GEO_min) < diff_alt) then
                diff_alt = abs(Altitude(iHoru,ialt_GCM,nTime)-Alt_GEO_min)
                ialt_min = ialt_GCM
              endif  
            enddo
            ialt_max = 1
            diff_alt = 4000.
            do ialt_GCM = 1,nAlti
              if (abs(Altitude(iHoru,ialt_GCM,nTime)-Alt_GEO_max) < diff_alt) then
                diff_alt = abs(Altitude(iHoru,ialt_GCM,nTime)-Alt_GEO_max)
                ialt_max = ialt_GCM
              endif  
            enddo
            
   
        endif ! end on the altitude check     
     
      enddo
    enddo
 enddo
 
 
! where(density_Oi <1.e-5) density_Oi = 0.
! where(density_O2i <1.e-5) density_O2i = 0.
! where(density_CO2i <1.e-5) density_CO2i = 0.
! where(density_Hi <1.e-5) density_Hi = 0.

! dumping the density values as initial density
!density_Op = density_Oi
!density_Hp = density_Hi
!density_O2p = density_O2i
!density_CO2p = density_CO2i

!atmosphere%species(1)%density_ini = density_CO2p
!atmosphere%species(2)%density_ini = density_Op
!atmosphere%species(3)%density_ini = density_Hp
!atmosphere%species(4)%density_ini = density_O2p
!density_ini(:,:,:,2) = density_Op
!density_ini(:,:,:,3) = density_Hp
!density_ini(:,:,:,4) = density_O2p


!!======V environment dependent code from here V===============
! k4=1.64E-10;k5=9.6E-11;k6=1.1E-9;k7=7.38E-8;k8=2.5E-11*sqrt(300.+2000./16.)
!    !coefficients des reactions (pour k8 cf these de O wittasse)
 npcell=150 !nombre max de particules par cellule et par espece
 r_lim=500._dp !altitude maximum a remplir
!!===================^  To here ^================================
!  
  r_lim=(r_lim/Spe%ref%c_omegapi+Spe%P%radius)
  r_lim2=(r_lim+sqrt(gstep(1)**2+gstep(2)**2+gstep(3)**2))**2
  rp2=(Spe%P%r_lim-sqrt(gstep(1)**2+gstep(2)**2+gstep(3)**2))**2
  s_cen = Spe%P%centr-s_min_loc
!  vol_unit=1E6/Spe%ref%density 
  min_i=max(int((s_cen-r_lim)/gstep)-1,1) !ou commence l'ionosphere, le -1 corrige INT pour des valeurs negatives
  max_i=min(int((s_cen+r_lim)/gstep),ncm-2) !ou finit l'ionosphere

  call create_ionosphere_generic(Spe,s_cen,s_min_loc,particule,gstep,min_i,max_i,irand,nptot,npcell,atmosphere)!does all the work
  write(msg,'(a,i10)')&
       & " max weight ",int(maxval(particule(:)%char))
    call wrtout(6,msg,'PERS')
!  endif

deallocate(Time,zls,sza,dsm,dec,Longitude,Latitude,Altitude,nCO2p,nOp,nO2p)
__WRT_DEBUG_OUT("load_ionosphere_venus_LMD")  
end subroutine load_ionosphere_venus

!!================================================================
!!routine : env_venus/Iono_venus
!!
subroutine Iono_venus(i,j,k,atmosphere,l,t3,r2,rb,Spe,s_cen)
  use defs_parametre, only :gstep
  type(species_type),intent(in) :: Spe
  type(atmosphere_type),intent(in) ::atmosphere
  real(dp),intent(in) ::r2,rb
  real(dp),intent(in),dimension(3)::s_cen
  real(dp),intent(inout) ::t3
  integer,intent(in) :: i,j,k,l
  real(dp) ::k4,k5,k6,k7,k8
  logical::no_coll

k4=1.64E-10;k5=9.6E-11;k6=9.4E-10;k7=1.6E-7*(300.0_dp/4000._dp)**0.55;k8=1.14E-4*1.0_dp/4000._dp
!k4=1.64E-10;k5=9.6E-11;k6=1.1E-9;k7=7.38E-8;k8=2.5E-11*sqrt(300.+2000./16.)
 !coefficients des reactions (pour k8 cf these de O wittasse)
        no_coll=((density_O(i,j,k).lt.1E8).and.(density_CO2(i,j,k).lt.1E8).and.(density_H(i,j,k).lt.1E8))
        if (no_coll) then
                if (atmosphere%species(l)%name.eq.("CO2+      ")) &
                & t3=prod_CO2(i,j,k)-(k4+k5)*density_O(i,j,k)*density_CO2p(i,j,k)
                if (atmosphere%species(l)%name.eq.("O+        ")) then
                        t3=prod_O(i,j,k)+k5*density_CO2p(i,j,k)*&
                        & density_O(i,j,k)-(k6*density_CO2(i,j,k)+&
                        & k8*density_H(i,j,k))*density_Op(i,j,k)
                        if (r2.gt.(300./Spe%ref%c_omegapi+Spe%P%radius)**2) & 
                        & t3=0.
                endif
                if (atmosphere%species(l)%name.eq.("O2+       ")) then
                        t3 =(k4*density_CO2p(i,j,k)*density_O(i,j,k)+&
                        & k6*density_CO2(i,j,k)*density_Op(i,j,k))&
                        &-k7*density_O2p(i,j,k)*(density_O2p(i,j,k)+&
                        & density_Op(i,j,k)+density_CO2p(i,j,k))
                        if ((rb.lt.(Spe%P%radius)**2).and.&
                        &((float(i-1)*gstep(1)-s_cen(1)).gt.0)) t3=0.
                endif
        else 
                t3=0.
        endif
end subroutine Iono_venus

!!================================================================
!!routine : env_venus/Iono_loaded_venus
!!
subroutine Iono_loaded_venus(i,j,k,atmosphere,l,t3,r2,rb,Spe,s_cen)
  use defs_parametre,only :gstep,dt
  type(species_type),intent(in) :: Spe
  type(atmosphere_type),intent(in) ::atmosphere
  real(dp),intent(in) ::r2,rb
  real(dp),intent(in),dimension(3)::s_cen
  real(dp),intent(inout) ::t3
  integer,intent(in) :: i,j,k,l
  real(dp) :: prod_unit
  
  prod_unit=1E6/Spe%ref%density*dt*Spe%ref%inv_gyro
    t3 = 0.

end subroutine Iono_loaded_venus

!!
!! --------------- end Iono_Venus_loaded --------------------

subroutine feed_ionosphere_venus(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
  use defs_parametre,only : ionospherename
  use defs_particletype
  use atm_ionosphere
  real(dp),intent(in) :: gstep(3),s_min_loc(3),s_max_loc(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(inout) ::atmosphere
  
  
  if (trim(ionospherename)/="") then
    call Ion_production_generic(Spe,atmosphere,particule,20.,0.01,nptot,irand,gstep,s_min_loc,Iono_loaded_venus)
  else
    call Ion_production_generic(Spe,atmosphere,particule,20.,0.01,nptot,irand,gstep,s_min_loc,Iono_venus)
  endif  

end subroutine feed_ionosphere_venus



end module env_venus

