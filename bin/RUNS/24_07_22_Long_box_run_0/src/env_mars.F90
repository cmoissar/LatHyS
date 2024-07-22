!!=============================================================
!!=============================================================
!!module: env_mars
!! NAME
!!  env_mars (SHess)
!!
!! Gathers all the mars related file into a single one, for more clarity
!!
!! FUNCTION
!!  Contains definition of type for the environment planet mars
!!
!! NOTE
module env_mars

 use defs_basis
 use defs_parametre
 use defs_species
 use defs_atmospheretype
 use atm_external_atmosphere
 use atm_charge_exchange
! use defs_parametre
 use m_writeout
#ifdef HAVE_NETCDF 
    use defs_basic_cdf
    use netcdf
#endif


#include "q-p_common.h"

 implicit none
 private

 integer,parameter :: CO2p=1,Op=2,Hp=3,O2p=4,O=5,CO2=6,Hn=7!,CO2i=8,Oi=9,Hi=10,O2i=11

 !--Pointer towards the density_exo allocated in m_exosphere
 real(dp),pointer ::        density_H(:,:,:),  &
      &                     density_Hp(:,:,:), &
      &                     density_O(:,:,:),  &
      &                     density_CO2(:,:,:),&
      &                     density_Op(:,:,:), &
      &                     density_O2p(:,:,:),&
      &                     density_CO2p(:,:,:)!, &
      !&   		    density_Oi(:,:,:),&
      !&			    density_O2i(:,:,:), &
      !&			    density_CO2i(:,:,:), &
      !&                     density_Hi(:,:,:)
 !--Array for Photoproduction Rate
 real(dp),pointer :: prod_O(:,:,:),prod_H(:,:,:),prod_CO2(:,:,:)

 public::                     &
      alloc_mars,             &!--allocate pp and density arrays for Mars
      dealloc_mars,           &!--deallocate pp and density arrays for Mars
      init_species_mars,      &!--Intialise the Mars environment
      exosphere_mars,         &!--Compute exosphere densities for Mars
      photoproduction_mars,   &!--Compute production ratio
      charge_exchange_mars,   &!--Charge exchange for Mars
      create_ionosphere_mars, &!--Compute ionosphere densities for Mars
      load_ionosphere_mars,   &!--Load ionsopheric densities for Mars
      feed_ionosphere_mars,   &
      mars_magnetic_field
contains
 !!#####################################################################
 !!=============================================================
 !!routine: env_mars/init_species_mars
 !!
 !! FUNCTION
 !!  Initialize species for mars
 !! IN 
 !! planet center
 !!
 !! OUT
 !! species 
 !!
 subroutine init_species_mars(species,s_centr)
 use defs_parametre

  real(dp),intent(in) :: s_centr(3)
  type(species_type),intent(inout) :: species
  
  integer,parameter :: H=1,He=2
  integer :: is
  real(dp) :: tot_tmp,vitessethermique

  !--Mars contains two species of particles
  !--Allocation of species_type
  call alloc_species(2,species)
  
  !--Name of the planet
  species%planetname = "mars"
  species%exospherename = "mars"
  
  !--Intialize Physical parameter used here
  !--Physical density of reference (m-3)
  species%ref%density =1.5e6  ! 4.e6 ! extreme event 20.e6

  !--Ions inertial length (km)
  species%ref%c_omegapi = Sp_Lt*sqrt(epsilon0*pmasse_e2/species%ref%density)/1.e3

  !--Magnetic Field (Tesla)
  species%ref%mag = 2.e-9 ! 3.e-9 ! extreme event 3.e-9

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
  species%P%radius = 3393._dp/species%ref%c_omegapi
  species%P%r_exo  = species%P%radius+200._dp/species%ref%c_omegapi
  species%P%r_lim  = species%P%radius+110._dp/species%ref%c_omegapi
  if (trim(ionospherename) /="") then
   ! species%P%r_iono = species%P%radius+300._dp/species%ref%c_omegapi
    species%P%r_iono = species%P%radius+250._dp/species%ref%c_omegapi
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
  species%S(H)%vxs = 370._dp/species%ref%alfvenspeed !--  H+  ! 12.18
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
  species%S(H)%betas = 1.57_dp   ! 1.57  ! Extreme event 0.12
  species%S(He)%betas = 0.03_dp  ! 0.03 !               0.024
  species%betae = 1.5_dp        ! 1.08 !               0.084

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


 end subroutine init_species_mars

 !!=============================================================
 !!routine: env_mars/alloc_mars
 !!
 !! FUNCTION
 !!  allocate photoproduction and exosphere arrays for Mars
 !! IN 
 !! density and photoproduction arrays, allocated in environment
 !!   
 subroutine alloc_mars(density,prod_pp,atmosphere,ncm)
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
                atmosphere%species(O2p)%prod => prod_pp(:,:,:,O2p)

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
  !O+ + O->O+ + O
!  atmosphere%exc_reactions(4)%qsm        = one/16._dp
!  atmosphere%exc_reactions(4)%ion        =>atmosphere%species(Op)
!  atmosphere%exc_reactions(4)%neutral    =>atmosphere%species(O)
!  atmosphere%exc_reactions(4)%cross_section = zero

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
 

 end subroutine alloc_mars

 !!=============================================================
 !!routine: env_mars/dealloc_mars
 !!
 !! FUNCTION
 !!  deallocate photoproduction and exosphere arrays for Mars
 !!   
 subroutine dealloc_mars(dummy)
  integer,intent(in) :: dummy
  density_O => Null()
  density_H => Null()
  density_Hp => Null()
  density_CO2 => Null()
  !density_Oi => Null()
  !density_O2i => Null()
  !density_CO2i => Null()
  !density_Hi => Null()
  prod_O => Null()
  prod_H => Null()
  prod_CO2 => Null()
 end subroutine dealloc_mars

 !!=============================================================
 !!subroutine: env_mars/exosphere_mars
 !! NAME
 !!  exosphere_mars (RModolo,MMancini)
 !!
 !! FUNCTION
 !!  Contains Exosphere calculatio for Mars
 !!
 !! NOTE
 !!  The densities results are in (cm-3) where the
 !!  during the calculation (for the altitude) the units are Km.
 !!  Written by R. Modolo 11/10/02 inspired from the subroutine of the 2D code
 subroutine exosphere_mars(Spe,ncm,gstep,s_min_loc,resistivity)
!  use defs_parametre
  use defs_variable,only : viscosity
  use defs_parametre,only :dt,exospherename,atmospherename
 
   integer, intent(in) :: ncm(3)
   real(dp),intent(in) :: gstep(3),s_min_loc(3)
   real(dp),intent(inout) :: resistivity(:,:,:)
!!!!!!!! Change by F. Leblanc 01/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   type(species_type),intent(inout) :: Spe
!   type(species_type),intent(in) :: Spe
!!!!!!!! Change by F. Leblanc 01/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   integer :: ii,jj,kk,iii,jjj,kkk
   real(dp) :: radius,altitude_km,r_planet_km,inv_r_exo_km,nrm
   real(dp) :: ss(3)
   real(dp),parameter :: d0_O = 2.33e13,e0_O = one/12.27,&
        &                d1_O = 2.84e9 ,e1_O = one/48.57,&
        &                d2_O = 1.5e4 ,e2_O = one/696.9,&
        &                d3_O = 2.92e3 ,e3_O = one/2891.,&
        &                d4_O = 5.01e4 ,e4_O = one/99.19,&
        &                d0_H = 1.e3  ,e0_H = 9.25e5,   &    ! solar max
        &                d1_H = 3.e4  ,e1_H = 1.48e4,    &   ! solar max
        !&                d0_H = 1.5e5  ,e0_H = 25965,    &    ! solar min
        !&                d1_H = 1.9e4  ,e1_H = 10365,    &    ! solar min
        &                d0_CO2 = 5.88e18,e0_CO2 = one/7, &
        &                d1_CO2 = 3.55e15, e1_CO2 = one/16.67
  character(len=500) :: msg          
 
   __WRT_DEBUG_IN("exosphere_mars")
   
   !--radius of the Planet in Km
   r_planet_km = Spe%P%radius*Spe%ref%c_omegapi
   inv_r_exo_km = one/(Spe%P%r_exo*Spe%ref%c_omegapi)
       density_H(:,:,:) = zero
       density_O(:,:,:) = zero
       density_CO2(:,:,:) = zero
   nrm=1._dp/729._dp
   !--Main loop
   do kk = 1,ncm(3)-1
    do jj = 1,ncm(2)-1
     do ii = 1,ncm(1)-1
      do kkk=-4,4
      do jjj=-4,4
      do iii=-4,4
      ss(3) = (real((kk-1),dp)+real(kkk,dp)*0.1111111)*gstep(3) + s_min_loc(3) 
      ss(2) = (real((jj-1),dp)+real(jjj,dp)*0.1111111)*gstep(2) + s_min_loc(2)
      ss(1) = (real((ii-1),dp)+real(iii,dp)*0.1111111)*gstep(1) + s_min_loc(1)      
      radius = sqrt(dot_product(ss-Spe%P%centr,ss-Spe%P%centr))
      if (radius >= Spe%P%r_lim) then
       !--Altitude in Kilometers
       altitude_km = (radius-Spe%P%r_lim)*(Spe%ref%c_omegapi)
       
       !--Pour l'hydrogene
       density_H(ii,jj,kk) = density_H(ii,jj,kk)+nrm*(d0_H*exp(e0_H*(&
            &   one/(altitude_km+r_planet_km)-inv_r_exo_km))+&
            &   d1_H*exp(e1_H*(&
            &   one/(altitude_km+r_planet_km)-inv_r_exo_km)))
 
       !--Pour l'Oxygene
       density_O(ii,jj,kk) = density_O(ii,jj,kk)+nrm*(d0_O*exp(-altitude_km*e0_O)+&
            &                d1_O*exp(-altitude_km*e1_O)+&
            &                d2_O*exp(-altitude_km*e2_O)+&
            &                d3_O*exp(-altitude_km*e3_O)+&
            &                d4_O*exp(-altitude_km*e4_O))
       !--For CO2
 !      density_CO2(ii,jj,kk) = density_CO2(ii,jj,kk)+nrm*(d0_CO2*exp(-(altitude_km-140._dp)*e0_CO2))
       density_CO2(ii,jj,kk) = density_CO2(ii,jj,kk)+nrm*(d0_CO2*exp(&
            &-altitude_km*e0_CO2)+d1_CO2*exp(-altitude_km*e1_CO2))
      endif
         enddo
         enddo
         enddo
     enddo
    enddo
   enddo
 
   viscosity(:,:,:)=Spe%viscosity*density_O
   
    if (trim(atmospherename) /= "") then
                write(msg,'(2a)')ch10,&
                "____________Loading atmosphere from external file_________"
                call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
                call load_atmosphere_mars_LMD(Spe,ncm,gstep,s_min_loc,resistivity,density_O,density_CO2)
                !call load_atmosphere_mars_Michigan(Spe,ncm,gstep,s_min_loc,resistivity,density_O,density_CO2)
     endif
   
       if (trim(exospherename) /= "") then
             write(msg,'(2a)')ch10,&
             "____________Loading exosphere from external file_________"
             call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
             call load_exosphere_mars(Spe,ncm,gstep,s_min_loc,resistivity)
     endif
  __WRT_DEBUG_OUT("exosphere_mars")
  end subroutine exosphere_mars


 !!=============================================================
 !!routine: env_mars/photoproduction_mars
 !!
 !! FUNCTION
 !!  computes photoproduction arrays for Mars
 !!   
 subroutine photoproduction_mars(Spe,ncm,gstep,s_min_loc,atmosphere)
  use atm_photoproduction
  use atm_sections_efficaces
  integer,intent(in) :: ncm(3)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type), intent(in) :: Spe
  type(atmosphere_type),intent(inout) ::atmosphere
  real(dp) :: F107,F107a,dist_conv

  atmosphere%F107     = 105.0_dp  ! daily F10.7 (e.g. 74)
  atmosphere%F107_Avg = 105.0_dp  ! 81 day average F10.7 (F10.7A) (e.g. 86)
  dist_conv=2.25_dp  ! Flux solaire a l'orbite de l'objet--(1.e4 conversion to photons/(m2*s))

  call flux_solaire_generic(atmosphere%nb_lo,atmosphere,atmosphere%F107,atmosphere%F107_Avg,dist_conv)
  call section_efficace_H_abs(atmosphere%nb_lo,atmosphere%species(Hn)%ion_abs)
  call section_efficace_O_abs(atmosphere%nb_lo,atmosphere%species(O)%ion_abs)
  call section_efficace_CO2_abs(atmosphere%nb_lo, atmosphere%species(CO2)%ion_abs)
  call section_efficace_CO2_ion(atmosphere%nb_lo,atmosphere%photo_reactions(1)%ion_react)
  call section_efficace_O_ion(atmosphere%nb_lo, atmosphere%photo_reactions(2)%ion_react)
  call section_efficace_H_ion(atmosphere%nb_lo,atmosphere%photo_reactions(3)%ion_react)
  call section_efficace_O_CO2_ion(atmosphere%nb_lo, atmosphere%photo_reactions(4)%ion_react)
  call photoproduction_generic(Spe,ncm,gstep,s_min_loc,atmosphere)! does all the work
 end subroutine photoproduction_mars

 !!=============================================================
 !!routine: env_mars/chrage_exchange_mars
 !!
 !! FUNCTION
 !!  computes charge exchange for Mars
 !!   
subroutine charge_exchange_mars(nn,kpickup,&
                 qsm,irand,ijk,v_p,ww,Spe,particule,atmosphere)
 use defs_parametre
 use defs_particletype
  integer,intent(in) :: nn,ijk(3)
  integer,intent(inout) :: irand,kpickup
  real(dp),intent(in) :: qsm,v_p(3)
  real(dp),intent(in) :: ww(8)
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(in) ::atmosphere
  real(dp) ::Va,conv_fac,vmod
  conv_fac = Spe%ref%c_omegapi*1.e5 ! 1.38e+7
  vmod  = sqrt(sum(v_p*v_p))*dt*conv_fac
   call  charge_exchange_generic(nn,kpickup,qsm,irand,&
      &                     ijk,v_p,ww,Spe,particule,atmosphere)! does all the work
end subroutine

 !!=============================================================
 !!routine: env_mars/create_ionosphere_mars
 !!
 !! FUNCTION
 !!  generate an ionosphere for Mars
 !!   
subroutine create_ionosphere_mars(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
!!  use defs_parametre
  use defs_particletype
  use atm_ionosphere
  use defs_grid,only : ncm

  real(dp),intent(in) :: s_min_loc(3),s_max_loc(3),gstep(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(in) ::atmosphere
  
  integer  ::min_i(3),max_i(3),i,j,k,cnt,npcell
  real(dp) ::s_cen(3),r_lim2,rp2,rb,radius,vol_unit,k4,k5,k6,k7,k8,r_lim
  character(len=100) :: msg



!======V environment dependent code from here V===============
 k4=1.64E-10;k5=9.6E-11;k6=1.1E-9;k7=7.38E-8;k8=2.5E-11*sqrt(300.+2000./16.)
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
end subroutine create_ionosphere_mars

!!=============================================================
 !!routine: env_mars/load_ionosphere_mars
 !!
 !! FUNCTION
 !!  load an ionosphere for Mars
 !!   
subroutine load_ionosphere_mars(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
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

 __WRT_DEBUG_IN("load_ionosphere_mars_LMD")
 
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
__WRT_DEBUG_OUT("load_ionosphere_mars_LMD")  
end subroutine load_ionosphere_mars

!!================================================================
!!routine : env_mars/Iono_mars
!!
subroutine Iono_mars(i,j,k,atmosphere,l,t3,r2,rb,Spe,s_cen)
  use defs_parametre,only :gstep
  type(species_type),intent(in) :: Spe
  type(atmosphere_type),intent(in) ::atmosphere
  real(dp),intent(in) ::r2,rb
  real(dp),intent(in),dimension(3)::s_cen
  real(dp),intent(inout) ::t3
  integer,intent(in) :: i,j,k,l
  real(dp) ::k4,k5,k6,k7,k8
  logical::no_coll

k4=1.64E-10;k5=9.6E-11;k6=1.1E-9;k7=7.38E-8;k8=2.5E-11*sqrt(300.+2000./16.)
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
end subroutine Iono_mars

!!================================================================
!!routine : env_mars/Iono_loaded_mars
!!
subroutine Iono_loaded_mars(i,j,k,atmosphere,l,t3,r2,rb,Spe,s_cen)
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

end subroutine Iono_loaded_mars

!!
!! --------------- end Iono_Mars_loaded --------------------

subroutine feed_ionosphere_mars(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
  use defs_parametre, only : ionospherename 
  use defs_particletype
  use atm_ionosphere
  real(dp),intent(in) :: gstep(3),s_min_loc(3),s_max_loc(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(inout) ::atmosphere
  
  
  if (trim(ionospherename)/="") then
    call Ion_production_generic(Spe,atmosphere,particule,20.,0.1,nptot,irand,gstep,s_min_loc,Iono_loaded_mars)
  else
    call Ion_production_generic(Spe,atmosphere,particule,20.,0.1,nptot,irand,gstep,s_min_loc,Iono_mars)
  endif  

end subroutine feed_ionosphere_mars

subroutine mars_magnetic_field(Bfield,ncm,Spe,gstep,s_min_loc)
 !use defs_parametre
 use defs_arr3Dtype
 use atm_magnetic_fields
 use env_mars_FSU90
  integer, intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type),intent(inout) :: Spe
  real(dp),dimension(1:90,0:90)::g,h
  g(:,:)=0.;h(:,:)=0.
        call FSU_mars(g,h); g=g*1.e-9;h=h*1.e-9;
        call add_multipole_generic(Bfield,ncm,Spe,gstep,s_min_loc,90,g,h)
end subroutine mars_magnetic_field

!!=============================================================
 !!subroutine: env_mars/load_exosphere_mars
 !! NAME
 !!  load_exosphere_mars (RModolo)
 !!
 !! FUNCTION
 !!  load Martian exosphere from external file
 !!
 !! NOTE
 !!  The densities results are in (cm-3) where the
 !!  during the calculation (for the altitude) the units are Km.
 
 subroutine load_exosphere_mars(Spe,ncm,gstep,s_min_loc,resistivity)
 use defs_variable,only : viscosity
 use defs_parametre,only :dt,exospherename
 use defs_mpitype,only     : mpiinfo

  integer, intent(in) :: ncm(3)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  real(dp),intent(inout) :: resistivity(:,:,:)
!!!!!!!! Change by F. Leblanc 01/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(species_type),intent(inout) :: Spe
!  type(species_type),intent(in) :: Spe
!!!!!!!! Change by F. Leblanc 01/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: ii,jj,kk,iii,jjj,kkk
  real(dp) :: radius,altitude_km,r_planet_km,inv_r_exo_km,nrm,ref_alt_exos
  real(dp) :: ss(3),x_cdr,y_cdr,z_cdr,tan_p,tan_t
  character(len=500) :: msg  
  integer :: npt_alt,npt_lat,npt_lon,stId,ncid,n_alt,n_phi,n_theta
  real(dp),dimension(:),allocatable :: altitude,theta,phi
  real(dp),dimension(:,:,:),allocatable :: densO,tempO,densCO2,tempCO2
  real(dp) :: diff_alt,diff_theta,diff_phi,rad_mars
  logical :: find_alt,find_phi,find_theta
  real(dp) :: g,n0_O,n0_CO2,H_scale_height_O,H_scale_height_CO2
  
  !density_O(:,:,:) = 1.e-10
  !density_CO2(:,:,:) = 1.e-10
  
  ! Read external file
 ! if (mpiinfo%me == 4) then
      call wrt_double(6,"Reading file : "//exospherename,wrtscreen,wrtdisk)
 print *,'Read exosphere file :',exospherename 
     stId = nf90_open(exospherename ,  nf90_nowrite, ncid)
     call test_cdf(stId)
     
     call get_simple_variable_cdf(ncid,"npt_alt",npt_alt)
     call get_simple_variable_cdf(ncid,"npt_lat",npt_lat)
     call get_simple_variable_cdf(ncid,"npt_lon",npt_lon)
!!!!!!!! Change by F. Leblanc 01/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(msg,'(2a,3(a,a17,i12))')&
       & ch10," ___________________ Exosphere file Information _________________",&
       & ch10, "   npt_alt  = ",npt_alt,&
       & ch10, "   npt_lat  = ",npt_lat,&
       & ch10, "   npt_lon  = ",npt_lon
     call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
     call get_simple_variable_cdf(ncid,"planetary_radius",rad_mars)
!     call get_simple_variable_cdf(ncid,"Radius",rad_mars)
!     call get_simple_variable_cdf(ncid,"ref_alt",ref_alt_exos)

     Spe%P%r_exo  = Spe%P%radius+220.d0/Spe%ref%c_omegapi
!     Spe%P%r_exo  = Spe%P%radius+ref_alt_exos*1.e-05/Spe%ref%c_omegapi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     allocate(altitude(npt_alt));       altitude(:) = 0.
     allocate(theta(npt_lat));  theta(:) = 0.
     allocate(phi(npt_lon));            phi(:) = 0.
     allocate(densO(npt_alt,npt_lat,npt_lon));  densO(:,:,:) = 0.
     allocate(tempO(npt_alt,npt_lat,npt_lon));  tempO(:,:,:) = 0.
     allocate(densCO2(npt_alt,npt_lat,npt_lon));        densCO2(:,:,:) = 0.
     allocate(tempCO2(npt_alt,npt_lat,npt_lon));        tempCO2(:,:,:) = 0.     
     
     call get_simple_variable_cdf(ncid,"altitude_low",altitude)
     call get_simple_variable_cdf(ncid,"theta_low",theta)
     call get_simple_variable_cdf(ncid,"phi_low",phi)
     call get_simple_variable_cdf(ncid,"density_O",densO)
     call get_simple_variable_cdf(ncid,"TempO",tempO)
     call get_simple_variable_cdf(ncid,"density_CO2",densCO2)
     call get_simple_variable_cdf(ncid,"TempCO2",tempCO2)     
     !--Close the file
     stId = nf90_close(ncid); call test_cdf(stId)

 write(6,'(a,2(1x,e12.5))')' Before Loading Exosphere Min, Max CO2 = ',minval(density_CO2(:,:,:)),maxval(density_CO2(:,:,:))
 write(6,'(a,2(1x,e12.5))')' Before Loading Exosphere Min, Max O   = ',minval(density_O(:,:,:)),maxval(density_O(:,:,:))
     
!!!!!!!! Change by F. Leblanc 01/2015 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     altitude(:) = altitude(:)*1.e-05/(Spe%ref%c_omegapi)
     !altitude(:) = altitude(:)*rad_mars/(Spe%ref%c_omegapi)
     write(msg,'(2a,2(a,a17,e12.5))') &
       & ch10," ___________________ Exosphere Altitude _________________",&
       & ch10, "   Min (km)  = ",minval(altitude(:))*Spe%ref%c_omegapi,&
       & ch10, "   Max       = ",maxval(altitude(:))*Spe%ref%c_omegapi
     call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
     write(msg,'(2a,2(a,a17,e12.5))')&
       & ch10," ___________________ Exosphere nO _________________",&
       & ch10, "   Min (cm-3)  = ",minval(densO(:,:,:)),&
       & ch10, "   Max         = ",maxval(densO(:,:,:))
     call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
     write(msg,'(2a,2(a,a17,e12.5))')&
       & ch10," ___________________ Exosphere nCO2 _________________",&
       & ch10, "   Min (cm-3)  = ",minval(densCO2(:,:,:)),&
       & ch10, "   Max         = ",maxval(densCO2(:,:,:))
     call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
       !--Main loop
       do kk = 1,ncm(3)-1
        do jj = 1,ncm(2)-1
         do ii = 1,ncm(1)-1
         ! do kkk=-4,4
         ! do jjj=-4,4
         ! do iii=-4,4
          ss(3) = (real((kk-1),dp))*gstep(3) + s_min_loc(3) 
          ss(2) = (real((jj-1),dp))*gstep(2) + s_min_loc(2)
          ss(1) = (real((ii-1),dp))*gstep(1) + s_min_loc(1)      
          radius = sqrt(dot_product(ss-Spe%P%centr,ss-Spe%P%centr))
          if (radius <= altitude(npt_alt) .and. (radius >= Spe%P%r_exo)) then

          x_cdr = -(ss(1)-Spe%P%centr(1))
          y_cdr = -(ss(2)-Spe%P%centr(2))
          z_cdr =  (ss(3)-Spe%P%centr(3))
          tan_p = atan(y_cdr/x_cdr)
          tan_t = atan(sqrt((x_cdr**2+y_cdr**2)/z_cdr**2))
! check conditions to obtain the good angles
          if (z_cdr == 0.) tan_t = pi/2._dp
          if (z_cdr < 0.) tan_t = pi-tan_t
          
          if ((x_cdr > 0.).and. (y_cdr < 0.)) tan_p = pi*2._dp+tan_p
          if ((x_cdr > 0.).and. (y_cdr == 0.)) tan_p = 0._dp
          
          if ((x_cdr < 0.).and. (y_cdr < 0.)) tan_p = pi+tan_p
          if ((x_cdr < 0.).and. (y_cdr == 0.)) tan_p = pi
          if ((x_cdr < 0.).and. (y_cdr >0.)) tan_p = pi+tan_p
          
          if ((x_cdr == 0.).and.(y_cdr > 0.)) tan_p = pi/2._dp
          if ((x_cdr == 0.).and.(y_cdr < 0.)) tan_p = pi*3._dp/2._dp
          if ((x_cdr == 0.).and.(y_cdr == 0.)) tan_p = 0._dp

          
          n_alt = 2
          n_phi = 2
          n_theta = 1
          diff_alt = radius-altitude(1)
          diff_phi = tan_p-phi(1)
          !diff_theta = tan_t - theta(1)
          diff_theta = pi
          find_alt = .false.
          find_theta = .false.
          find_phi = .false.
          ! Find the good altitude
          do while (( find_alt.EQV. .FALSE.).and.(n_alt <= npt_alt))
            if (abs(radius-altitude(n_alt)) < abs(diff_alt)) then
              diff_alt = radius-altitude(n_alt)
              n_alt = n_alt+1
            else
              find_alt =.true.
            endif
          enddo
          ! Find the good phi
          do while ((find_phi .EQV. .false.).and.(n_phi <= npt_lon))
            if (abs(tan_p-phi(n_phi)) < abs(diff_phi)) then
              diff_phi = tan_p-phi(n_phi)
              n_phi = n_phi+1
            else
              find_phi =.true.
            endif
          enddo
          ! Find the good theta
          do while ((find_theta .EQV. .false.).and.(n_theta <= npt_lat))
            if (abs(tan_t-theta(n_theta)) < abs(diff_theta)) then
              diff_theta = tan_t-theta(n_theta)
              n_theta = n_theta+1
            else
              find_theta =.true.
            endif
          enddo
         ! if (((find_alt==.true.).and.(find_phi==.true.)).and.(find_theta == .true.)) then
         if (find_phi .EQV. .false.) n_phi = 1
         if (find_theta .EQV. .false.) n_theta = 1
         
         if (find_alt .EQV..true.)  density_O(ii,jj,kk) = densO(n_alt,n_theta,n_phi)
         if (find_alt .EQV..true.)  density_CO2(ii,jj,kk) = densCO2(n_alt,n_theta,n_phi)         
          
!          if ((radius < altitude(1))) then
!            g = M_Mars*G_grav/((altitude(2)*Spe%ref%c_omegapi*1.e3)**2)
!            
!            H_scale_height_O = (kb_JK*tempO(2,n_theta,n_phi))/(g*16.*amu_pmass)
!            H_scale_height_CO2 = (kb_JK*tempCO2(2,n_theta,n_phi))/(g*44.*amu_pmass)
!            
!            n0_CO2 = densCO2(2,n_theta,n_phi)*exp((altitude(2)*Spe%ref%c_omegapi-rad_mars)*1.e3/H_scale_height_CO2)
!            n0_O = densO(2,n_theta,n_phi)*exp((altitude(2)*Spe%ref%c_omegapi-rad_mars)*1.e3/H_scale_height_O)
!            
!            density_O(ii,jj,kk) = n0_O*exp(-(radius*Spe%ref%c_omegapi-rad_mars)*1.e3/H_scale_height_O)
!            density_CO2(ii,jj,kk) = n0_CO2*exp(-(radius*Spe%ref%c_omegapi-rad_mars)*1.e3/H_scale_height_CO2)
!
!                   
!          endif  
          if (radius < Spe%P%radius) density_O(ii,jj,kk) = 1.e-10
          if (radius < Spe%P%radius) density_CO2(ii,jj,kk) = 1.e-10
            
                  
          !if ((radius >= altitude(1)).and.(radius <= altitude(nr))) then
           ! n_alt = minloc(abs(altitude - radius),dim=1)
           ! n_phi = minloc(abs(tan_p-tan(phi)),dim=1)
           ! n_theta = minloc(abs(tan_t-tan(theta)),dim=1)
            !print *,'Altitude :',n_alt,n_theta,n_phi,dens(n_alt,n_theta,n_phi)
           ! density_O(ii,jj,kk) = dens(n_alt,n_theta,n_phi)
          !endif
          endif   !on the altitude
         ! enddo
         ! enddo
         ! enddo
         enddo
       enddo
     enddo
     
    stId = nf90_close(ncid) 
    deallocate(altitude,phi,theta,densO,tempO,densCO2,tempCO2)
 ! endif
  
 write(6,'(a,2(1x,e12.5))')' Loading Exosphere Min, Max CO2 = ',minval(density_CO2(:,:,:)),maxval(density_CO2(:,:,:))
 write(6,'(a,2(1x,e12.5))')' Loading Exosphere Min, Max O   = ',minval(density_O(:,:,:)),maxval(density_O(:,:,:))

  
  end subroutine load_exosphere_mars



end module env_mars

