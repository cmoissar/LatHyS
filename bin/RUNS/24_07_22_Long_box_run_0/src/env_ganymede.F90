!!=============================================================
!!=============================================================
!!module: env_ganymede
!! NAME
!!  env_ganymede (MMancini, RModolo)
!!
!! FUNCTION
!!  Contains definition of type for the environment ganymede
!!
!! NOTE
module env_ganymede

 use defs_basis
 use defs_species
 use defs_particletype
 use defs_mpitype
 use defs_parametre
 use defs_atmospheretype
 use m_writeout
 use atm_sections_efficaces
 use atm_charge_exchange
 use mpi
     
#include "q-p_common.h"

 implicit none
 private
 
 private ::  &
        set_particle_iono 


 !--Pointer towards the density_exo allocated in m_exosphere
 real(dp),pointer ::    density_H(:,:,:), &     ! neutral density of atomic hydogen array
                        density_H2(:,:,:), &! neutral density of molecular hydrogen array
                        density_H2O(:,:,:), &! neutral density of water array
                        density_O2(:,:,:), &    ! neutral density of molecular oxygen
                        density_Hp(:,:,:), &! ion density array of H+
                        density_H2p(:,:,:), & ! ion density array of H2+
                        density_Op(:,:,:),&     ! ion density array of O+
                        density_OHp(:,:,:), &! ion density array of OH+
                        density_H2Op(:,:,:), &! ion density array of H2O+
                        !density_Oi(:,:,:),&    ! ion density array at initialisation
                        prod_H(:,:,:), &! production array of H+
                        prod_H2(:,:,:), &! production array of H2+
                        prod_O(:,:,:), &! production array of O+
                        prod_OH(:,:,:), &! production array of OH+
                        prod_H2O(:,:,:), &! production array of H2O+
                        prod_O2(:,:,:)! production array of O2+

real(dp),allocatable :: density_Oi(:,:,:),prod_temp(:,:,:)
 
 public ::                 &
      init_species_ganymede,  &     !--Intialise the Ganymede
      add_b_dipole_ganymede,  &    ! add dipolar moment for ganymede
      create_ionosphere_ganymede, &
      alloc_ganymede,             &!--allocate pp and density arrays for Ganymede
      dealloc_ganymede,           &!--deallocate pp and density arrays for Ganymede
      exosphere_ganymede,          &
      photoproduction_ganymede,     &
      feed_ionosphere_ganymede,    &
      charge_exchange_ganymede
contains
 !!#####################################################################


 !!=============================================================
 !!routine: m_species/init_species_ganymede
 !!
 !! FUNCTION
 !!  Initialize species for ganymede
 !! IN 
 !! OUT
 !! SIDE EFFECT
 !!  
 !!
 subroutine init_species_ganymede(species,s_centr)

  real(dp),intent(in) :: s_centr(3)
  type(species_type),intent(inout) :: species
  
  integer,parameter :: O=1,H=2
  integer :: is
  real(dp) :: tot_tmp,vitessethermique

  !--Ganymede contains two species of particles

  !--Allocation of species_type
  call alloc_species(2,species)

  !--Name of the planet
  species%planetname = "ganymede"

  !--Intialize Physical parameter used here
  !--Physical density of reference (m-3)
  species%ref%density = 3.47e6

  !--Ions inertial length
  species%ref%c_omegapi = Sp_Lt*sqrt(epsilon0*16*pmasse_e2/species%ref%density)/1.e3

  !--Magnetic Field (Tesla)
  species%ref%mag = sqrt(2._dp)*79.e-9

  !--Alfven Speed
  species%ref%alfvenspeed = species%ref%mag/&
       &                   (sqrt(mu0*16._dp*amu_pmass*species%ref%density)*1.e3)

  !--Inverse of gyrofrequency
  species%ref%inv_gyro = 16*amu_pmass/(e_Cb*species%ref%mag)

  !--Assignation of Planet values
  species%P%centr  = s_centr
#ifndef HAVE_NO_PLANET
  species%P%radius = 2631._dp/species%ref%c_omegapi !/529._dp
  species%P%r_exo  = species%P%radius !TO change!!!!!!!!
  species%P%r_iono = species%P%radius+1000._dp/species%ref%c_omegapi !TO change!!!!!!!!
  species%P%r_lim  = species%P%radius !TO change!!!!!!!
  species%P%speed  = zero 
#endif  

  !--Oxigen and Hydrogen
  species%S(:)%name = (/"O  ","H  "/)
  
  !--Number of particles for cells
  species%S(:)%ng = (/6, 6/) !(/2,2/)

  !--Direct Speeds  (vitesses dirigees?)
  species%S(O)%vxs = 140._dp/species%ref%alfvenspeed !--  O+
  species%S(H)%vxs = species%S(O)%vxs !-- H+
  species%S(:)%vys  = zero
  species%S(:)%vzs  = zero

  !--Percentage of He in the solar wind: n(He++)/(n(He++)+n(H+))
  species%S(H)%percent = 2._dp/15._dp
  species%S(O)%percent = one-species%S(H)%percent

  !--Ratio of charges and masses and Temperatures
  species%S(:)%rcharge = (/one,one/)
  species%S(:)%rmass   = (/one,one/16._dp/)
  species%S(:)%rtemp   = (/one,one/16._dp/)

  !--Betas
  species%S(H)%betas = 0.006_dp
  species%S(O)%betas = 0.04_dp
  species%betae = 0.01_dp

  !--Rapport des vitesses thermque entre parallel et perpendiculair  H+ et He++
  species%S(:)%rspeed = one

  species%tempe_ratio = 0.001_dp !  3 eV

  !--Rapport des vitesse thermiques entre H+ et He++
  species%S(:)%rvth = one
  species%S(:)%rmds = one

  !--Rapport charge sur masse
  species%S(:)%qms  = species%S(:)%rcharge/species%S(:)%rmass

  !--Thermal speed (parallel and perpendicular) 
  vitessethermique = sqrt( three*species%S(O)%betas/(one+two*species%S(O)%rvth**two))

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
 
 end subroutine init_species_ganymede

 !!=============================================================
 !!=============================================================
 !!subroutine: b_dipole/add_b_dipole_ganymede
 !! NAME
 !!  add_b_dipole_mars (RModolo,MMancini,RAllieux)
 !!
 !! FUNCTION
 !!  Contains Dipole moment calculation gfor ganymede
 !!
 !! NOTE
 !! The magnetic field inpu at initialization is derived from Kivelson et Al 2002, Icarus 157, 507522 (2002).
 !! We use the dipolar model superposed with an ambient field and the value calculated for the flyby G2
 !!
 !!
 subroutine add_b_dipole_ganymede(Bfield,ncm,Spe,gstep,s_min_loc)
 use defs_arr3Dtype

  integer, intent(in) :: ncm(3)
  type(arr3Dtype),intent(inout) :: Bfield
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type),intent(inout) :: Spe

  !local
  integer :: ii,jj,kk
  real(dp) :: radius
  real(dp) :: ss(3),moment_dip(3),b_dip(3),moment_dip_u(3)
  real(dp) :: inclinaison, phase,rinclinaison,rphase,bme
  real(dp) :: M_Jia(3)
  
  M_Jia(1) = - 24.7e-9
  M_Jia(2) =   82.5e-9
  M_Jia(3) = -716.8e-9
  
  !bme = 720.e-9 !T equateur
  bme = sqrt(dot_product(M_Jia,M_Jia))
  
  !inclinaison = 180.  !angle axe Z
  !phase = 0.   !tps local en degres
  
  inclinaison = acos(M_Jia(3)/bme) ! angle axe Z
  phase       = acos(M_Jia(1)/sqrt(M_Jia(1)**2+M_Jia(2)**2)) ! tps local
  
  !rinclinaison = inclinaison*deg_to_rad
  !rphase = phase*deg_to_rad
  rinclinaison = inclinaison
  rphase = phase
  

  __WRT_DEBUG_IN("add_b_dipole_ganymede")
  !--Initialisation
  !--Dipolar Moment in planetary units (Spe%ref%mag) mu0/4pi*M
  !moment_dip = mu0/four_pi*1.83e17*(/11.,66.,-728./)/(Spe%ref%mag*(Spe%ref%c_omegapi*1.e3)**three)
   moment_dip_u(1) = sin(rinclinaison)*cos(rphase)
   moment_dip_u(2) = sin(rinclinaison)*sin(rphase)
   moment_dip_u(3) = cos(rinclinaison)
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
          !b_dip = bme*(three*ss*  &
                  !dot_product(ss,moment_dip_u)/        &
                           !(radius*radius)-moment_dip_u)/ &
                           !(radius**three*Spe%ref%mag)
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

  __WRT_DEBUG_OUT("add_b_dipole_ganymede")
 end subroutine add_b_dipole_ganymede
!!=============================================================


 !!=============================================================
 !!subroutine m_ionosphere/create_ionosphere_ganymede
 !!
 !! NAME
 !!  create_ionosphere_ganymede (RModolo,MMancini,RAllioux)
 !!
 !! FUNCTION
 !!  Compute ionosphere profil
 !!
 subroutine create_ionosphere_ganymede(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)

  !use mpi
!  use defs_parametre
  use defs_mpitype,only: mpiinfo
  use defs_atmospheretype
  use atm_ionosphere
  use defs_grid,only : ncm

  real(dp),intent(in) :: gstep(3),s_min_loc(3),s_max_loc(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(in) ::atmosphere

  integer,parameter :: ip_max=100000 ! max number of ionospheric particles
  integer :: ip,ip_fin,ip_tot,ip_tot_all,ioerr,i,j,k,min_i(3),max_i(3),npcell
  real(dp),parameter :: density_iono_max = 520.
  real(dp) :: Hscale,Riono,density_iono_max_H,density_iono_max_O
  real(dp) :: nb_cell_load,r_lim,r_lim2,s_cen(3),vol_unit
  real(dp) :: masschar_fac(2),inv_Hscale,rp2,rb,radius
  character(len=500) :: msg   
  

  ! allocate density ionospheric array
  allocate(density_Oi(ncm(1),ncm(2),ncm(3))); density_Oi(:,:,:) = zero
! allocate temporary production array
  allocate(prod_temp(ncm(1),ncm(2),ncm(3))); prod_temp(:,:,:) = zero

  !--Initialization
  density_iono_max_O = density_iono_max
  density_iono_max_H = density_iono_max
  
  Hscale = 125._dp/(Spe%ref%c_omegapi) !--Scale height of the inosphere=125km
  Riono = ten*Hscale !--Size of the inosphere

  inv_Hscale = one/Hscale

  nb_cell_load = four_thirds*pi*((Spe%P%radius+Riono)**three-(Spe%P%radius)**three)/product(gstep)
  masschar_fac = Spe%S(:)%ng/(real(ip_max,dp)/two)*nb_cell_load

  write(*,'(a,i4.4,a,i8,a)')&
       &' Process ',mpiinfo%me,' has ',nptot,'particles before loading ionosphere'
  !========================================================
  !!rame side

  !--Compute profil 
  !ip_fin = nptot
  !write(*,'(a,i3,a,i8,a)')&
  !     &' Process ',mpiinfo%me,' has ',nptot,'particles before loading O+ rame side'
  !do ip=1,ip_max/2
  ! !--Ram side
  ! !!=====for O+
  ! call set_particle_iono(particule(nptot+1),&
  !      &                 s_min_loc,s_max_loc,&
  !      &                 Riono,inv_Hscale,density_iono_max_O,&
  !      &                 masschar_fac,&
  !      &                 1,0,&
  !      &                 irand,nptot,&
  !      &                 Spe)
!
!
  ! !!=====for H+
  ! call set_particle_iono(particule(nptot+1),&
  !      &                 s_min_loc,s_max_loc,&
  !      &                 Riono,inv_Hscale,density_iono_max_H,&
  !      &                 masschar_fac,&
  !      &                 2,0,&
  !      &                 irand,nptot,&
  !      &                 Spe)  
!
  ! !--Wake side 
!
  ! !!=====for O+
  ! call set_particle_iono(particule(nptot+1),&
  !      &                 s_min_loc,s_max_loc,&
  !      &                 Riono,inv_Hscale,density_iono_max_O,&
  !      &                 masschar_fac,&
  !      &                 1,1,&
  !      &                 irand,nptot,&
  !      &                 Spe)
!
  ! !!=====for H+
  ! call set_particle_iono(particule(nptot+1),&
  !      &                 s_min_loc,s_max_loc,&
  !      &                 Riono,inv_Hscale,density_iono_max_H,&
  !      &                 masschar_fac,&
  !      &                 2,1,&
  !      &                 irand,nptot,&
  !      &                 Spe)  
!
  !enddo
  !ip_fin = nptot-ip_fin
  !ip_tot = ip_fin
  !write(*,'(a,i3,a,i8,a)')&
  !     &' Process ',mpiinfo%me,' has ',ip_fin,'O+ and H+ loaded in the ionosphere(rame)'
  !
  !call MPI_REDUCE(ip_tot,ip_tot_all,1,MPI_INTEGER,MPI_SUM,0,&
  !     &          mpiinfo%comm,ioerr)
  !write(msg,'(3a,3x,a,i12,a)')ch10,&
  !     &" ___________________ Initialization Ionosphere _____________",ch10,&
  !     &"Loaded in Ionosphere ",ip_tot_all," particles O+ and H+"
  !call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
    r_lim=2000._dp !altitude maximum a remplir
     npcell=150 !nombre max de particules par cellule et par espece
    
    r_lim=(r_lim/Spe%ref%c_omegapi+Spe%P%radius)
    r_lim2=(r_lim+sqrt(gstep(1)**2+gstep(2)**2+gstep(3)**2))**2
    rp2=(Spe%P%r_lim-sqrt(gstep(1)**2+gstep(2)**2+gstep(3)**2))**2
    s_cen = Spe%P%centr-s_min_loc
    vol_unit=1E6/Spe%ref%density 
    min_i=max(int((s_cen-r_lim)/gstep)-1,1) !ou commence l'ionosphere, le -1 corrige INT pour des valeurs negatives
    max_i=min(int((s_cen+r_lim)/gstep),ncm-2) !ou finit l'ionosphere
    if(all(min_i.lt.max_i)) then !y a t-il de l'ionosphere dans la boite?
     
  
  !======V environment dependent code from here V===============
    density_Op(:,:,:) = 0.
    do k=min_i(3),max_i(3)+1
          do j=min_i(2),max_i(2)+1
              rb = (float(j-1)*gstep(2)-s_cen(2))**2+(float(k-1)*gstep(3)-s_cen(3))**2
                  do i=min_i(1),max_i(1)+1
                        radius = (float(i-1)*gstep(1)-s_cen(1))**2+rb
                if (((radius > rp2).and.(radius <= r_lim2))) then !!.or.((rb < rp2).and.(float(i)*gstep(1) > s_cen(1)))
  
                    density_Op(i,j,k) = density_iono_max_O*exp(-(sqrt(radius)-Spe%P%radius)*inv_Hscale)
                    density_Oi(i,j,k) = density_iono_max_O*exp(-(sqrt(radius)-Spe%P%radius)*inv_Hscale)
                endif
        
                enddo
        enddo
 enddo
 
 ! normalisation

   density_Op=density_Op*vol_unit
   density_Oi = density_Oi*vol_unit
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

 end subroutine create_ionosphere_ganymede


 !!=============================================================
 !!=============================================================
 !!subroutine m_ionosphere/set_particle_iono
 !!
 !! NAME
 !!  create_ionosphere_ganymede (RModolo,MMancini,RAllioux)
 !!
 !! FUNCTION
 !!  Compute ionosphere profil single entries
 !!
 subroutine set_particle_iono(particle,s_min_loc,s_max_loc,&
      &                       Riono,inv_Hscale,den_iono_species,&
      &                       masschar_fac,is,wake_or_ram,&
      &                       irand,nptot,Spe)

  use m_rand_gen,only                   : rand_gen1 


  integer,intent(in) :: is,wake_or_ram
  real(dp),intent(in) :: Riono,inv_Hscale,den_iono_species
  integer,intent(inout) :: nptot,irand
  real(dp),intent(in) :: masschar_fac(2)
  real(dp),intent(in) :: s_max_loc(3),s_min_loc(3)
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particle

  real(dp) :: theta,alt,phi,ratio
  real(dp) :: poss(3)

  theta = acos(two*rand_gen1(irand)-one)
  phi   = -half_pi + pi*rand_gen1(irand)
  alt = Riono*rand_gen1(irand)
  if(wake_or_ram == 1) phi = phi + pi
  poss(1) = spe%P%centr(1) - (Spe%P%radius+alt)*sin(theta)*cos(phi) 
  poss(2) = Spe%P%centr(2) + (Spe%P%radius+alt)*sin(theta)*sin(phi)
  poss(3) = Spe%P%centr(3) + (Spe%P%radius+alt)*cos(theta)
  !--test on the positions: same proc
  if(all(poss< s_max_loc).and.(all(poss> s_min_loc)))then
   nptot = nptot+1     
   ratio = den_iono_species*exp(-alt*inv_Hscale)
   particle%pos = poss
   particle%vel = zero
   particle%exc = 0
   particle%orig = 3
   particle%mass = Spe%S(is)%sm*ratio*masschar_fac(is)
   particle%char = Spe%S(is)%sq*ratio*masschar_fac(is)
  endif
 end subroutine set_particle_iono
 
 
!!=============================================================
 !!routine: env_ganymede/alloc_ganymede
 !!
 !! FUNCTION
 !!  allocate photoproduction and exosphere arrays for ganymede
 !! IN
 !! density and photoproduction arrays, allocated in environment
 !!
 subroutine alloc_ganymede(density,prod_pp,atmosphere,ncm)
  integer,dimension(3),intent(in) :: ncm
  real(dp),intent(inout),allocatable,target :: density(:,:,:,:)
  real(dp),intent(inout),allocatable,target :: prod_pp(:,:,:,:)  
  type(atmosphere_type),intent(inout),target :: atmosphere
  !integer,parameter :: H=1 
  !integer,parameter :: H2=2
         !integer,parameter :: H2O=3  

  atmosphere%n_species= 9 !nombre d'especes atmospherique (neutres et ions)
  atmosphere%n_spe_pp=5 !nombre d'especes obtenues par photoproduction
  atmosphere%n_pp=6      !nombre de reactions de photoproduction
  atmosphere%nb_lo=37   !nombre de longueur d'onde du spectre UV
  atmosphere%n_exc=2   !nombre de reactions d'echange de charge

 call allocate_atmosphere(ncm,atmosphere,density,prod_pp)
   density_Hp    => density(:,:,:,1)
                   atmosphere%species(1)%name   ="Hp    "
                   atmosphere%species(1)%mass   = 1._dp/16._dp
                   atmosphere%species(1)%charge = 1._dp
                   atmosphere%species(1)%opaque = .FALSE.
                   atmosphere%species(1)%iono   = .TRUE.
                   atmosphere%species(1)%prod   => prod_pp(:,:,:,1)
 
   density_H2p   => density(:,:,:,2)
                   atmosphere%species(2)%name   = "H2p   "
                   atmosphere%species(2)%mass   = 2._dp/16._dp
                  atmosphere%species(2)%charge = 1._dp
                  atmosphere%species(2)%opaque = .FALSE.
                  atmosphere%species(2)%iono   = .TRUE.
                   atmosphere%species(2)%prod   => prod_pp(:,:,:,2)
 
   density_Op    => density(:,:,:,3)
                   atmosphere%species(3)%name   ="Op    "
                   atmosphere%species(3)%mass   = 1._dp
                   atmosphere%species(3)%charge = 1._dp
                   atmosphere%species(3)%opaque = .FALSE.
                   atmosphere%species(3)%iono   = .TRUE.
                   atmosphere%species(3)%prod   => prod_pp(:,:,:,3)
                   
   density_OHp    => density(:,:,:,4)
                   atmosphere%species(4)%name   ="OHp    "
                   atmosphere%species(4)%mass   = 1._dp + 1._dp/16._dp
                   atmosphere%species(4)%charge = 1._dp
                   atmosphere%species(4)%opaque = .FALSE.
                   atmosphere%species(4)%iono   = .TRUE.
                   atmosphere%species(4)%prod   => prod_pp(:,:,:,4)
 
   density_H2Op    => density(:,:,:,5)
                   atmosphere%species(5)%name   ="H2Op    "
                   atmosphere%species(5)%mass   = 1._dp + 2._dp/16._dp
                   atmosphere%species(5)%charge = 1._dp
                   atmosphere%species(5)%opaque = .FALSE.
                   atmosphere%species(5)%iono   = .TRUE.
                  atmosphere%species(5)%prod   => prod_pp(:,:,:,5)
                  

  !density_Oi    => density(:,:,:,6)
  !                atmosphere%species(6)%name   ="Oi    "
  !                atmosphere%species(6)%mass   = 1._dp
  !                atmosphere%species(6)%charge = 1._dp
  !                atmosphere%species(6)%opaque = .FALSE.
  !                atmosphere%species(6)%iono   = .FALSE.                           
 
  density_H    => density(:,:,:,6)
                  atmosphere%species(6)%name  ="H      "
                  atmosphere%species(6)%mass  = 1._dp/16._dp
                  atmosphere%species(6)%charge= zero
                  atmosphere%species(6)%opaque= .FALSE.
                  atmosphere%species(6)%iono  = .FALSE.
  
 density_H2    => density(:,:,:,7)
                  atmosphere%species(7)%name  ="H2     "
                  atmosphere%species(7)%mass  =  2._dp/16._dp
                  atmosphere%species(7)%charge= zero
                  atmosphere%species(7)%opaque= .FALSE.
                  atmosphere%species(7)%iono  = .FALSE.
                  
 density_H2O    => density(:,:,:,8)
                  atmosphere%species(8)%name  ="H2O    "
                  atmosphere%species(8)%mass  =  2._dp/16._dp+1._dp
                  atmosphere%species(8)%charge= zero
                  atmosphere%species(8)%opaque= .FALSE.
                  atmosphere%species(8)%iono  = .FALSE.



  density_O2    => density(:,:,:,9)
                  atmosphere%species(9)%name   ="O2    "
                  atmosphere%species(9)%mass   = 2._dp
                  atmosphere%species(9)%charge = 1._dp
                  atmosphere%species(9)%opaque = .FALSE.
                  atmosphere%species(9)%iono   = .FALSE.

                    
                  
!--Associate pointers to photoproduction array
  prod_H    => prod_pp(:,:,:,1)
  prod_H2   => prod_pp(:,:,:,2)
  prod_O    => prod_pp(:,:,:,3)
  prod_OH   => prod_pp(:,:,:,4)
  prod_H2O  => prod_pp(:,:,:,5)
  !prod_O2  => prod_pp(:,:,:,6)
  

  !H->H+
  atmosphere%photo_reactions(1)%mother     =>atmosphere%species(6)   ! neutral specie H(index 6 from above)
  atmosphere%photo_reactions(1)%daughter   =>atmosphere%species(1)   ! ion specie Hp (index 1 from above)
  !H2->H2+
  atmosphere%photo_reactions(2)%mother     =>atmosphere%species(7)   ! neutral specie H2 (index 7 from above)
  atmosphere%photo_reactions(2)%daughter   =>atmosphere%species(2)   ! ion specie H2p (index 2 from above)
  !H2O->H+
  atmosphere%photo_reactions(3)%mother     =>atmosphere%species(8)   ! neutral specie H2O (index 8 from above)
  atmosphere%photo_reactions(3)%daughter   =>atmosphere%species(1)   ! ion specie Hp (index 1 from above)
  !H2O -> O+
  atmosphere%photo_reactions(4)%mother     =>atmosphere%species(8)   ! neutral specie H2O (index 8 from above)
  atmosphere%photo_reactions(4)%daughter   =>atmosphere%species(3)   ! ion specie Op (index 3 from above)
  !H2O -> OH+
  atmosphere%photo_reactions(5)%mother     =>atmosphere%species(8)   ! neutral specie H2O (index 8 from above)
  atmosphere%photo_reactions(5)%daughter   =>atmosphere%species(4)   ! ion specie OHp (index 4 from above)  
  !H2O -> H2O+
  atmosphere%photo_reactions(6)%mother     =>atmosphere%species(8)   ! neutral specie H2O (index 8 from above)
  atmosphere%photo_reactions(6)%daughter   =>atmosphere%species(5)   ! ion specie H2Op (index 5 from above)  
  !O2 -> O2+
  !atmosphere%photo_reactions(7)%mother     =>atmosphere%species(10)   ! neutral specie O2 (index 10 from above)
  !atmosphere%photo_reactions(7)%daughter   =>atmosphere%species(6)   ! ion specie O2p (index 6 from above)  
  
  
    !--Associate pointers for charge exchange
  !H + H+ ->H+ + H
  atmosphere%exc_reactions(1)%qsm        = 16.
  atmosphere%exc_reactions(1)%ion        =>atmosphere%species(1)     ! Hp (index 1)
  atmosphere%exc_reactions(1)%neutral    =>atmosphere%species(6)     ! H (index 6)
  atmosphere%exc_reactions(1)%cross_section = 2.E-19
  !H2 + H+ ->H2+ + H
  atmosphere%exc_reactions(2)%qsm        = 8.
  atmosphere%exc_reactions(2)%ion        =>atmosphere%species(2)     ! H2p (index 2)
  atmosphere%exc_reactions(2)%neutral    =>atmosphere%species(7)     ! H2 (index 7)
  atmosphere%exc_reactions(2)%cross_section = 2.E-20
  
 end subroutine alloc_ganymede
 
 
 !!=============================================================
 !!routine: env_ganymede/dealloc_ganymede
 !!
 !! FUNCTION
 !!  deallocate photoproduction and exosphere arrays for ganymede
 !!
 subroutine dealloc_ganymede(dummy)
  integer,intent(in) :: dummy
  density_H => Null()
  prod_H =>Null()
  prod_H2 =>Null()
  density_H2 => Null()
  density_H2O => Null()
  prod_H2O=>Null()
  !density_O2 => Null()
  !prod_O2 => Null  
  deallocate(density_Oi,prod_temp)
  
 end subroutine dealloc_ganymede
 
 
 
 !!=============================================================
 !!=============================================================
 !!subroutine: env_ganymede/exosphere_ganymede
 !! NAME
 !!  exosphere_mars (RModolo,MMancini,RAllioux)
 !!
 !! FUNCTION
 !!  Contains Exosphere calculatio for ganymede
 !!
 !! NOTE
 !!  The densities results are in (cm-3) where the
 !!  during the calculation (for the altitude) the units are Km.
 !!  Written by R. Modolo 11/10/02 inspired from the subroutine of the 2D code

 subroutine exosphere_ganymede(Spe,ncm,&
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
  real(dp),save :: d_H,e_H,d_H2,e_H2,d_H2O, e_H2O,d_H2O_2,e_H2O_2, &
                d_O2,d_O2_2,e_O2,e_O2_2


  __WRT_DEBUG_IN("exosphere_ganymede")



  !--radius of the Planet in Km
  r_planet_km = Spe%P%radius*Spe%ref%c_omegapi
  inv_r_exo_km = one/(Spe%P%r_exo*Spe%ref%c_omegapi)
  ! the following values are average values determined from digitizing Marconi et al (2007) profiles
       d_H2 = 3.e6
       e_H2 = one/600._dp
       d_H = 5.e3 
       e_H = one/1000._dp
       d_H2O=1.e8
       e_H2O=one/50._dp
       d_H2O_2 = 2.e3
       e_H2O_2 = one/1000._dp
       d_O2 = 3.e6
       e_O2 = one/40._dp
       d_O2_2 = 4.e3
       e_O2_2 = one/1500._dp
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
      density_H2(ii,jj,kk) = d_H2*exp(-altitude_km*e_H2)
      density_H(ii,jj,kk) = d_H*exp(-altitude_km*e_H)
      density_H2O(ii,jj,kk) = d_H2O*exp(-altitude_km*e_H2O)+ &
        d_H2O_2*exp(-altitude_km*e_H2O_2)
      density_O2(ii,jj,kk) = d_O2*exp(-altitude_km*e_O2)+ &
        d_O2_2*exp(-altitude_km*e_O2_2) 
      else
      density_H2(ii,jj,kk) = zero
      density_H(ii,jj,kk) = zero
      density_H2O(ii,jj,kk) = zero
      density_O2(ii,jj,kk) = zero
     endif
    enddo
   enddo
  enddo
  __WRT_DEBUG_OUT("exosphere_ganymede")
 end subroutine exosphere_ganymede
 
 !!=============================================================
 !!routine: env_ganymede/photoproduction_ganymede
 !!
  subroutine photoproduction_ganymede(Spe,ncm,gstep,s_min_loc,atmosphere)
  use atm_photoproduction
  integer,intent(in) :: ncm(3)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type), intent(in) :: Spe
  type(atmosphere_type),intent(inout) ::atmosphere
  real(dp) :: F107,F107a,dist_conv
 ! integer,parameter :: H=1 
 ! integer,parameter :: H2=2
        ! integer,parameter :: H2O=3  
  F107  = 74._dp  ! daily F10.7 (e.g. 74)
  F107A = 74._dp  ! 81 day average F10.7 (F10.7A) (e.g. 86)
  dist_conv=27.07_dp  ! Flux solaire a l'orbite de l'objet--(1.e4 conversion to photons/(m2*s))

  call flux_solaire_generic(atmosphere%nb_lo,atmosphere,F107,F107a,dist_conv)
  call section_efficace_H_abs(atmosphere%nb_lo,atmosphere%species(7)%ion_abs)
  call section_efficace_H2_abs(atmosphere%nb_lo,atmosphere%species(8)%ion_abs)
  call section_efficace_H2O_abs(atmosphere%nb_lo,atmosphere%species(9)%ion_abs)
  
  call section_efficace_H_ion(atmosphere%nb_lo,atmosphere%photo_reactions(1)%ion_react)
  call section_efficace_H2_ion(atmosphere%nb_lo,atmosphere%photo_reactions(2)%ion_react)
  call section_efficace_H_H2O_ion(atmosphere%nb_lo,atmosphere%photo_reactions(3)%ion_react)
  call section_efficace_O_H2O_ion(atmosphere%nb_lo,atmosphere%photo_reactions(4)%ion_react)
  call section_efficace_OH_H2O_ion(atmosphere%nb_lo,atmosphere%photo_reactions(5)%ion_react)
  call section_efficace_H2O_ion(atmosphere%nb_lo,atmosphere%photo_reactions(6)%ion_react)
  !call section_efficace_O2_ion(atmosphere%nb_lo,atmosphere%photo_reactions(7)%ion_react)
 
                call photoproduction_generic(Spe,ncm,gstep,s_min_loc,atmosphere)! does all the work
 end subroutine photoproduction_ganymede
 !****************************** END PHOTOPRODUCTION ************************************
   

!!=============================================================
!!routine: env_ganymede/charge_exchange_ganymede
!!
!! FUNCTION
!!  computes charge exchange for Ganymede
!!   
subroutine charge_exchange_ganymede(nn,kpickup,&
                 qsm,irand,ijk,v_p,w,Spe,particule,atmosphere)
 use defs_particletype
   integer,intent(in) :: nn,ijk(3)
  integer,intent(inout) :: irand,kpickup
  real(dp),intent(in) :: qsm,v_p(3),w(8)
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(inout) ::atmosphere
  real(dp) ::Va,conv_fac,vmod
   __WRT_DEBUG_IN("charge_exchange_gany")
   call  charge_exchange_generic(nn,kpickup,qsm,irand,&
      &                     ijk,v_p,w,Spe,particule,atmosphere)! does all the work
       __WRT_DEBUG_OUT("charge_exchange_gany")
end subroutine




!!=============================================================
!!routine: env_ganymede/feed_ionosphere_ganymede
!!
!! FUNCTION
!!  inject photoproduct particule
!!   

subroutine feed_ionosphere_ganymede(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)
  use defs_particletype
  use atm_ionosphere
  real(dp),intent(in) :: gstep(3),s_min_loc(3),s_max_loc(3)
  integer,intent(inout) :: nptot,irand
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(inout) ::atmosphere
   real(dp) :: prod_unit
   
  prod_unit=1E6/Spe%ref%density*dt*Spe%ref%inv_gyro
  
  prod_temp = atmosphere%species(3)%prod
  
  where (density_Oi > atmosphere%species(3)%density)
    atmosphere%species(3)%prod = atmosphere%species(3)%prod!+ &
   !      & (density_Oi-atmosphere%species(3)%density)/prod_unit
  end where   

  call Ion_production_generic(Spe,atmosphere,particule,20.,0.01,nptot,irand,gstep,s_min_loc,dummy)
  
  atmosphere%species(3)%prod = prod_temp

end subroutine feed_ionosphere_ganymede
  
 !!================================================================
 !!routine : env_ganymede/Iono_ganymede
 !!
 subroutine Iono_ganymede(i,j,k,atmosphere,l,t3,r2,rb,Spe,s_cen)
   use defs_parametre,only :gstep,dt
   type(species_type),intent(in) :: Spe
   type(atmosphere_type),intent(in) ::atmosphere
   real(dp),intent(in) ::r2,rb
   real(dp),intent(in),dimension(3)::s_cen
   real(dp),intent(inout) ::t3
   integer,intent(in) :: i,j,k,l
   real(dp) :: prod_unit
 
   prod_unit=1E6/Spe%ref%density*dt*Spe%ref%inv_gyro
   
   if ((atmosphere%species(l)%name.eq.("Op        ")) .and. (density_Oi(i,j,k) > density_Op(i,j,k))) then
     t3 = (density_Oi(i,j,k)-density_Op(i,j,k))/prod_unit
     !print *, 't3 = ',t3*prod_unit
        else 
                t3=0.
        endif   
        
end subroutine Iono_ganymede
 
 !****************************** END pghoto_injection ************************************
end module env_ganymede
