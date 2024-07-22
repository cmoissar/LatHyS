!!=============================================================
!!=============================================================
module initialisation

 use defs_basis
 use defs_variable
 use defs_parametre
 use m_writeout
 use m_timing,only     : time_get
 use m_distribution_function
 use mpi

#include "q-p_common.h"

 implicit none
 private

 private ::                     &
      init3,                    &
      print_procs_distrib,      &!--Prt infos process distribution(debug) 
      print_init_infos           !--Prt infos concerning allocation

 public  ::                     &
      allocation,               &
      deallocation,             &
      h3init

contains
 !!#####################################################################

 !**************************** ALLOCATION DES TABLEAUX ******************************* 
 subroutine allocation()

  use defs_arr3Dtype
  use defs_diag_type,only      : init_diag_type
  use defs_particletype 
  use defs_grid
  use defs_species,only           : alloc_species
  use field_lissage,only             : smth_init_mpi_planes
  use field_cond_limit,only          : bound_init_mpi_planes
 
  integer :: tot_size_arr

  __WRT_DEBUG_IN("allocation")
  __GETTIME(2,1) !--Timer start

  tot_size_arr = 0

  !--Allocation des particules
  allocate(particule(npm)) ;   call set_zero_particle(particule,1,npm)

  !--Allocation field and velocity arrays
  call alloc_arr3D(Efield,ncm)
  call alloc_arr3D(Bfield,ncm)
  call alloc_arr3D(vel,   ncm)  
  call alloc_arr3D(Bfield_0,ncm)

  !--TMP allocation field and velocity arrays
  call alloc_arr3D(Bfield_h,ncm)
  call alloc_arr3D(vela,ncm)
  call alloc_arr3D(velf,ncm)

  !--Allocations des champs
  allocate(pe(ncm(1),ncm(2),ncm(3)))  ; pe  = zero 
  allocate(dn(ncm(1),ncm(2),ncm(3)))  ; dn  = zero 
  allocate(dna(ncm(1),ncm(2),ncm(3))) ; dna = zero 
  allocate(dnf(ncm(1),ncm(2),ncm(3))) ; dnf = zero 
  
  allocate(dn_e_pl(ncm(1),ncm(2),ncm(3))) ; dn_e_pl = zero
  allocate(dn_e_incdt(ncm(1),ncm(2),ncm(3))) ; dn_e_incdt = zero

  allocate(viscosity(ncm(1),ncm(2),ncm(3)))      ; viscosity      = zero
  allocate(pos_choc(ncm(2),ncm(3)))           ; pos_choc    = 0

  !allocate(resistivity(ncm(1),ncm(2),ncm(3))) ; resistivity = zero
  !--resistivity actually is not used

  tot_size_arr = tot_size_arr+21*ncm(1)*ncm(2)*ncm(3)

  !--Allocate species
  call alloc_species(ns,Spe)

  allocate(index_exit(ncross)) ; index_exit = 0

  allocate(iwr(nwrm)); iwr = 0  !--Allocatin array for Lissage

  call init_diag_type(diag_info,nhm,ns) !--Intialize array for diagnostic information
  tot_size_arr = tot_size_arr+23*nhm 

  call print_init_infos(tot_size_arr,npm) 

  !--Initialisation of MPI types defined
  call init_MPI_particle()
  call smth_init_mpi_planes(ncm)
  call bound_init_mpi_planes(ncm)

  __GETTIME(2,2)!--Timer stop
  __WRT_DEBUG_OUT("allocation")
 end subroutine allocation
 !******************* FIN D'ALLOCATION DES TABLEAUX *****************************

 !*******************************DEALLOCATION ******************************************
 subroutine deallocation() 
  
  use defs_arr3Dtype
  use defs_tregister,only         : treg_type,clean_tregister
  use defs_diag_type,only      : clean_diag_type
  use defs_species,only           : dealloc_species
  use environment
  use field_lissage,only             : smth_free_mpi_planes
  use field_cond_limit,only          : bound_free_mpi_planes

  __WRT_DEBUG_IN("deallocation")

  !--Deallocation des particules
  deallocate(particule)

  !--Dallocations des champs
  call dealloc_arr3D(Efield)
  call dealloc_arr3D(Bfield)
  call dealloc_arr3D(Bfield_h)
  call dealloc_arr3D(Bfield_0)
  call dealloc_arr3D(vel)
  call dealloc_arr3D(vela)
  call dealloc_arr3D(velf)

  deallocate(pe)
  deallocate(dn)
  deallocate(dna)
  deallocate(dn_e_pl)
  deallocate(dn_e_incdt)
  deallocate(dnf)
  !deallocate(resistivity)
  deallocate(viscosity)
  deallocate(pos_choc)
  
  !--Deassociate Charge_exchange pointers
  call nullify_environment()

  !--Deallocation des tableaux pour les especes du vent solaire
  call dealloc_species(Spe)

  deallocate(iwr)!--Deallocate smoothing control

  !--Clean variables for diagnostic
  call clean_tregister(r_tm)
  call clean_diag_type(diag_info)

  !--Deallocation des tableaux pour la gestion des sorties de macroparticules
  deallocate(index_exit)

  !--Free of MPI types defined
  call smth_free_mpi_planes()
  call bound_free_mpi_planes()
  call free_MPI_particle()

  __WRT_DEBUG_OUT("deallocation")
 end subroutine deallocation
 !*************************** END DEALLOCATION ***********************************************

 !******************************* H3INIT *****************************************************
 subroutine h3init()
  !Initialistion des parametres plasma pour la simulation
  
 ! use mpi
  use defs_mpitype,only         : mpiinfo,nproc
  use defs_grid
  use defs_species,only            : print_species
  use defs_tregister,only          : treg_type,set_tregister
  use environment
  use m_rand_gen,only             : rand_gen1,bi_max_dib,irand
  use field_pe
  use particle_fluxes
  use field_e,only              : ecalc3
  use particle_init,only        : pldf1 
  use part_moment,only               : momtin,Fmtsp3
  use diagnostique     !all
  use particle
  use field_add_waves
  use defs_parametre, only : pos_plan

#ifdef HAVE_WAVE_TEST
  use wavetest,only           : init_wavetest,Amtsp3_wt
#endif

  integer :: ii,is,ireg,iunit,nptot_all,ioerr
  real(dp) :: densco,rphi,rpsi
  real(dp) :: coupdslo
  real(dp) :: s_lon(3)                  !--Length of the box in any Direction 
  real(dp) :: s_min_par(3),s_max_par(3) !--Particles max and min in any direction
  character(len=500) :: msg

  ! Choose a cell where magnetic field values are recorded every time_step
  character(len=1):: a

  __WRT_DEBUG_IN("h3init")  
  __GETTIME(3,1)!--Timer start

  write(msg,'(2a)')ch10,&
       " _____________ Initialisation of Diagnostics ________________"
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
  
  !--Allocate and set variables for time diagnostics
  call set_tregister(r_tm)

  if(restart==0)then
   !--Box size (units of c/wi)
   s_min = zero
   s_max = gstep*real(nc_tot,dp)
   s_r = s_max-s_min
  endif

  !--Definition de la grille relatif a chaque processus
  s_lon = (s_max-s_min)/real((/1,mpiinfo%dims(1),mpiinfo%dims(2)/),dp)
  if ((abs(s_min_loc(2)-s_min(2))).lt.0.001) s_min_loc(2)=s_min(2)
  if ((abs(s_min_loc(3)-s_min(3))).lt.0.001) s_min_loc(3)=s_min(3)
  if ((abs(s_max_loc(2)-s_max(2))).lt.0.001) s_max_loc(2)=s_max(2)
  if ((abs(s_max_loc(3)-s_max(3))).lt.0.001) s_max_loc(3)=s_max(3)
 
 
  if(restart==0)then
  select case(trim(planetname))
  case("moon","mars","venus","mars3try","titan","mercure")
    pos_plan = int((s_max-s_min)/two+s_min)
  case("ganymede")
    pos_plan(1) = int((s_max(1)-s_min(1))*two/three+s_min(1))
    pos_plan(2:3) = int((s_max(2:3)-s_min(2:3))/two+s_min(2:3))
  case("earth")
    pos_plan(1) = int((s_max(1) - 100))
   ! pos_plan(1) = int((s_max(1)-s_min(1))*14.0/15.0+s_min(1))
   ! write(*,*) "pos_plan(1)=", pos_plan(1)
    pos_plan(2) = int((s_max(2)-s_min(2))/two+s_min(2)) !*two/three+s_min(2))
    pos_plan(3) = int((s_max(3)-s_min(3))/two+s_min(3)) 
   ! pos_plan(2:3) = int((s_max(2:3)-s_min(2:3))/two+s_min(2:3))
   ! pos_plan = int((s_max-s_min)/two+s_min)
  case default
    write(msg,'(a,3x,a,4x,3a)')ch10,&
         "ERROR: Selected Environment:",&
         "Planet '",trim(planetname),"' does Not have planet position"
    call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
    stop
  end select
  
 !--Intialisation of Species 
   call init_species(Spe,pos_plan)
   
   densco = sum((Spe%S(:)%percent)*(Spe%S(:)%rmass))  !--Densite du vent solaire   
   vxmean = sum(Spe%S(:)%vxs*Spe%S(:)%percent*Spe%S(:)%rmass)/densco !--Vitesse moyenne du vent solaire
   vx0 = vxmean

   !--Angle du champ IMF
   rphi = phi*deg_to_rad
   rpsi = psi*deg_to_rad

   !--Valeur de B imposee l'entree de la boite

   b0 = one

   bx0  = b0*cos(rphi)*sin(rpsi)  
   by0  = b0*sin(rphi)*sin(rpsi)
   bz0  = b0*cos(rpsi)

   by1 = by0
   bz1 = bz0

   !--Initialisation du champ magntique
   !--Magnetic field (magnitude of background reference field)
   Bfield%x = bx0 
   Bfield%y = by0 
   Bfield%z = bz0 
   Bfield_0%x = 0._dp!bx0
   Bfield_0%y = 0._dp!by0
   Bfield_0%z = 0._dp!bz0

   !--Parametre de lissage
   iwr( 1) = 1
   iwr( 6) = 1
   !iwr( 7) = 1
   iwr(11) = 2
   iwr(13) = 1
   iwr(15) = 1
  else
   write(msg,'(2a)')ch10,&
       " _______________ Restarted from a file ______________________"
   call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
  endif !endif restart==O

  write(msg,'(3a,3(a,f6.3),2a,f6.3)')ch10,&
       &" ________________ Discretization Step _______________________",&
       &ch10,'   gstep          ',gstep(1),',',gstep(2),',',gstep(3),&
       &ch10,'   time step (dt) ',dt
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  !--Print species informations
  call print_species(Spe)

  !--Parameters in .data file
  if(mpiinfo%me == 0 .and. wrtdisk == 0) then
   iunit = 10
   if(restart == 1)then
    open(UNIT   = iunit,&
         FILE   = fildat//'.restart.data',&
         STATUS = 'unknown')
   else
    open(UNIT   = iunit,&
         FILE   = fildat//'.data',&
         STATUS = 'unknown')
   endif
   write(iunit,'(a15,1pe12.4)') 'dt',dt
   write(iunit,'(a15,i8)')      'nsub',nsub
   write(iunit,'(a15,i8)')      'ntest',ntest
   write(iunit,'(a15,1pe12.4)') 'eps',eps
   write(iunit,'(a15,i8)')      'nhm',nhm

   write(iunit,'(a15,i8)')  'nreg',r_tm%nreg
   write(iunit,'(a15,2a12)') 'filreg','treg'
   write(iunit,'(a22)') " Intitialisation file:"
   write(iunit,'(a10,a15)')  "          ",trim(r_tm%file(0))
   write(iunit,'(a26)') " Diagnostique time files :"
   do ireg =1,r_tm%nreg
    write(iunit,'(a10,a15,f12.3)') "          ",trim(r_tm%file(ireg)),r_tm%time(ireg)
   enddo
   write(iunit,'(a15,i8)') 'nc_tot(1)',nc_tot(1)
   write(iunit,'(a15,i8)') 'nc_tot(2)',nc_tot(2)
   write(iunit,'(a15,i8)') 'nc_tot(3)',nc_tot(3)
   write(iunit,'(a15,1pe12.4)') 'dx',gstep(1)
   write(iunit,'(a15,1pe12.4)') 'dy',gstep(2)
   write(iunit,'(a15,1pe12.4)') 'dz',gstep(3)
   write(iunit,'(a15,f12.4)') 'psi',psi
   write(iunit,'(a15,f12.4)') 'phi',phi
   write(iunit,'(a15,f12.4)') 'betae',Spe%betae
   write(iunit,'(a15,i8)') 'ns',ns

   write(iunit,'(a15,5i8)') 'ng',Spe%S(:)%ng
   write(iunit,'(a15,5(1pe12.4))') 'betas',Spe%S(:)%betas
   write(iunit,'(a15,5(1pe12.4))') 'rvth' ,Spe%S(:)%rvth
   write(iunit,'(a15,5(1pe12.4))') 'rmds' ,Spe%S(:)%rmds
   write(iunit,'(a15,5(1pe12.4))') 'qms'  ,Spe%S(:)%qms
   write(iunit,'(a15,5(1pe12.4))') 'vxs'  ,Spe%S(:)%vxs
   write(iunit,'(a15,5(1pe12.4))') 'vys'  ,Spe%S(:)%vys
   write(iunit,'(a15,5(1pe12.4))') 'vzs'  ,Spe%S(:)%vzs
   write(iunit,'(a15,5i8)') 'idisp',idisp
   write(iunit,'(a15,5i8)') 'ipe',ipe
   write(iunit,'(a15,5i8)') 'iresis',iresis
   write(iunit,'(a15,1pe12.4)') 'resis',resis
   write(iunit,'(a15,1pe12.4)') 'dmin',dmin  
   write(iunit,'(a15)') 'iwr'
   write(iunit,'(15x,5i4)') (iwr(ii),ii=1,nwrm)
   close(iunit)
  endif
  
  write(msg,'(2(2a,f14.7),2a,3(e12.4,a))')&
       &ch10,'   vxmean                             ',vxmean,&
       &ch10,'   resistivity (resis)                ',resis,&
       &ch10,'   bx0,by0,bz0 =      ',bx0,', ',by0,', ',bz0,ch10
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
  
  !--Compute dipole moment contribution to magnetic field
  if (restart==0)then
    call add_b_dipole(Bfield_0,ncm,Spe,gstep,s_min_loc)
   Bfield_0%x = Bfield_0%x/(t_init_dip+0._dp)
   Bfield_0%y = Bfield_0%y/(t_init_dip+0._dp)
   Bfield_0%z = Bfield_0%z/(t_init_dip+0._dp)
    Bfield%x=Bfield%x+Bfield_0%x!*0.01_dp
    Bfield%y=Bfield%y+Bfield_0%y!*0.01_dp
    Bfield%z=Bfield%z+Bfield_0%z!*0.01_dp
   endif
 
  !--Compute densities of neutral species in the exosphere
  
  call add_exosphere(Spe,ncm,gstep,s_min_loc,resistivity)
 ! call MPI_BARRIER(mpiinfo%comm,ioerr) 
!--Compute photoproduction and  frequencies in the exosphere
  call calc_photoproduction(Spe,ncm,gstep,s_min_loc)
  !==============================================================
 
  if (diag_part_imp .ne. 0) call init_flux_part_imp()


  !--Basic Parameters only if NOT RESTART
  if(restart==0)then
   !--create temporal_B/ a bit ahead of call init_temporal_B()
   !--otherwise, some procs with me>0 might get there too fast,
   !--when temporal_B/ does not exist yet
   !if (mpiinfo%me==0) then
   !  write(*,*) "creating temporal_B directory"
   !  call system("mkdir temporal_B")
   !endif
   !--Initialisation des parametres pour la simulation
   call init3(betai,cd0,rho0,rmu0i,eps0,&
        &     te,cs,t,iter,&
        &     vac,Spe%betae,&
        &     Spe%S(:)%betas,&
        &     Spe%S(:)%rmds,Spe%S(:)%qms&
        &     )
  !-- Init waves only at the beginning, let them propagate then
   call init_waves()
   !call init_temporal_B()
  endif

  !--resistivity actually is not used
  !resistivity= resis



#ifdef HAVE_WAVE_TEST
   !--Initialisation when Wave test 
   call init_wavetest()
#endif

  ! Attention, these two next lines will have to go if we use
  ! field_add_waves. TODO: make a compilation option taking care
  ! of this.
  by1 = by0
  bz1 = bz0     

  !--Initialisation du champ electrique de convection pour la
  !   condition au bord gauche de la boite
  e_conv = (/zero,vxmean*bz0,-vxmean*by0/)
  e1_conv = (/zero,vxmean*bz1,-vxmean*by1/)

  !write(*,*) e_conv, e1_conv
#ifdef HAVE_DEBUG   
  call print_procs_distrib()
#endif

  if(restart==0)then
  !--Initialisation of Random seed
  if(nproc==1) then
   irand = -1
  else
   irand = mpiinfo%me+1
   coupdslo = rand_gen1(irand)
   irand = 0
  endif
  else
   irand = 0
  endif
  
  !--Initialisation des particules
  if(restart == 0) then
   call pldf1(irand,particule,Spe,s_min_loc,s_max_loc,&
        &     nc,nptot,npm)

   s_min_par(1) = minval(particule(1:nptot)%pos(1))
   s_min_par(2) = minval(particule(1:nptot)%pos(2))
   s_min_par(3) = minval(particule(1:nptot)%pos(3))
   s_max_par(1) = maxval(particule(1:nptot)%pos(1))
   s_max_par(2) = maxval(particule(1:nptot)%pos(2))
   s_max_par(3) = maxval(particule(1:nptot)%pos(3))

#ifdef HAVE_DEBUG   
    write(msg,'(a,i4.4,a,i8,a,3(2a,i4.4,a,2(f9.4)))')&
        &' Process ',mpiinfo%me,' has ',nptot,'particles at initialisation',ch10,&
        &' Process ',mpiinfo%me,' has min and max on X',&
        & s_max_par(1),s_min_par(1),ch10,&
        &' Process ',mpiinfo%me,' has min and max on Y',&
        & s_max_par(2),s_min_par(2),ch10,&
        &' Process ',mpiinfo%me,' has min and max on Z',&
        & s_max_par(3),s_min_par(3)
   call wrtout(6,msg,'PERS')
#endif

   !--Total nuber of macro-particles in the box
   call MPI_REDUCE(nptot,nptot_all,1,MPI_INTEGER,MPI_SUM,0,mpiinfo%comm,ioerr)  
   write(msg,'(2a,1(a,a12,i15))')&
        & ch10," _________ Particles _________",&
        & ch10, "   Nptot  = ",nptot_all
   call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

   !--Distribution de vitesse bi-Maxwellienne
   do is = 1,ns
    call bi_max_dib(irand,Spe%S(is)%vth1,Spe%S(is)%vth2,Spe%S(is)%n1,Spe%S(is)%n2,&
         &          particule(1:nptot)%vel(1),particule(1:nptot)%vel(2), &
         &          particule(1:nptot)%vel(3),Spe%S(is)%rspeed)
    !--On ajoute la vitesse dirigee
    do ii = Spe%S(is)%n1, Spe%S(is)%n2
     particule(ii)%vel = particule(ii)%vel +&
          &             (/ Spe%S(is)%vxs,Spe%S(is)%vys,Spe%S(is)%vzs/)
    enddo
   enddo
  endif
   
   !--Add Ionosphere particles
   if (restart==0) then
     call add_ionosphere(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot)
   endif


#ifdef HAVE_WAVE_TEST
   call Amtsp3_waves(particule,vel,1,nptot,gstep,s_min_loc,nc1)
   call Amtsp3_wt(particule,vel,1,nptot,gstep,s_min_loc,nc1)
#endif

  call Fmtsp3(particule,dna,dn_e_pl,dn_e_incdt,vela,1,nptot,gstep,s_min_loc,mpiinfo,nc1)
  
  
  call momtin(particule,dn,dna,vel,vela,gstep,s_min_loc,iwr,nptot)

  call pecalc(ipe,te,dna,pe)

  !--Here in ecalc I have to use dn (like in first sequential F90
  ! version) and not dna. The use of dna implies divergence in
  ! magnetic energy and structure creation in the fields
  call ecalc3(Efield,Bfield,vela,dn,pe,resistivity,&
       &      e_conv,e1_conv,gstep, &
       &      dmin,rmu0,resis,&
       &      idisp,nc1, &
       &      mpiinfo)

  do is = 1,ns
   call compute_fluxes(Spe,is,nfl)
  enddo

  call load_temporal_change(Spe)
  
  if (distrib_activated == 1) call initialisation_distribution_function

  if (restart == 0.and.wrtdisk == 0)  call diag_all(r_tm%file(0))

  __GETTIME(3,2)!--Timer stop
  __WRT_DEBUG_OUT("h3init")

 end subroutine h3init
 !******************************** END H3INIT **************************************************

 !*********************************** INIT3 **************************************************
 subroutine init3(betai,cd0,rho0,rmu0i,eps0,&
      &           te,cs,t,iter,&
      &           vac,betae,&
      &           betas,&
      &           rmds,qms&
      &           )

  use defs_parametre,only    : restart

  integer,intent(out) :: iter
  real(dp),intent(in) :: vac,betae
  real(dp),intent(out) :: &                      
       &                  betai,cd0,rho0,&
       &                  rmu0i,eps0,te,&
       &                  cs,t

  real(dp),dimension(:),intent(in) :: &
       &                              betas,& !--Beta pour chaque espece
       &                               rmds,& !--Rapport des masses
       &                                qms   !--Rapport charge sur masse par espece

  __WRT_DEBUG_IN("init3")

  !--Species
  cd0 = sum(qms*rmds) 
  betai = sum(betas*rmds)
  
  !--rho0 (reference density)
  rho0 = one
  !--mu0 and inverse (vA = B0/sqrt(mu0*rho0)= 1 : mu0*rho0 = 1)
  rmu0i= rho0
  rmu0 = one/rmu0i
  !--eps0
  eps0 = rmu0i*(vac**two)
  !--Electron temperature
  te = betae/(two*cd0)
  !--Sound speed
  cs = sqrt(half*(betai+betae))

  !--Simulation time (in units of 1/wg(proton))
  if (restart == 0) then
   !--Time
   t = zero
   !--Iteration
   iter = 0
  endif

  __WRT_DEBUG_OUT("init3")
 end subroutine init3
 !************************************ END INIT3 **********************************************
 


 !!=============================================================
 !!routine: initialisation/print_init_infos
 !!
 !! FUNCTION
 !! Print informations concerning grid and particles
 !! initialisation. And print approximative Memory requirement
 !!         
 !! IN
 !!  tot_size_arr= Number of scalar components of field and particule arrays
 !!  npt= Dimension of the particule array
 subroutine print_init_infos(tot_size_arr,npm) 

  use defs_mpitype
  use defs_particletype,only   : particle_type_size
  use defs_grid

  integer,intent(in) :: tot_size_arr
  integer,intent(in) :: npm

  real(dp) :: mem_bit_proc
  integer(1) :: sizeof(9999)
  character(len=3) :: unit
  character(len=500) :: msg


  mem_bit_proc = real(tot_size_arr,dp)*real(size(transfer(1.0_dp,sizeof)),dp)
  mem_bit_proc = mem_bit_proc+real(npm,dp)*real(particle_type_size(),dp)

  write(msg,'(2a,6(a,a23,3i4),2(a,a23,i12))')&
       & ch10," ___________________ Grid Informations ______________________",&
       & ch10, "   nc_tot(1:3)         = ",nc_tot(1),nc_tot(2),nc_tot(3),&
       & ch10, "   nc1_tot(1:3)        = ",nc1_tot(1),nc1_tot(2),nc1_tot(3),&
       & ch10, "   ncm_tot(1:3)        = ",ncm_tot(1),ncm_tot(2),ncm_tot(3),&
       & ch10, "   nc(1:3)             = ",nc(1),nc(2),nc(3),&
       & ch10, "   nc1(1:3)            = ",nc1(1),nc1(2),nc1(3),&
       & ch10, "   ncm(1:3)            = ",ncm(1),ncm(2),ncm(3),&
       & ch10, "   nxyzm               = ",nxyzm,&
       & ch10, "   ncross              = ",ncross
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  write(msg,'(2a,2(a,a17,i12))')&
       & ch10," ___________________ Particles Informations _________________",&
       & ch10, "   n_part_max  = ",n_part_max,&
       & ch10, "   npm         = ",npm
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  !--Conversion in GB or MB
  mem_bit_proc = b2Gb*mem_bit_proc*real(nproc,dp)
  unit = ' GB'
  if(mem_bit_proc < one) then
   mem_bit_proc = mem_bit_proc*1000.0_dp
   unit = ' MB'   
  endif
  
  write(msg,'(2a,2(a,a17,f8.3,a))')&
       & ch10," ___________________ Memory Required ________________________",&
       & ch10, "   Total Memory  = ",mem_bit_proc,unit,&
       & ch10, "   Mem. per proc = ",mem_bit_proc/real(nproc,dp),unit
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

 end subroutine print_init_infos

 !!=============================================================
 !!routine: initialisation/print_procs_distrib
 !!
 !! FUNCTION
 !! Print informations concerning grid and particles
 !! initialisation. And print approximative Memory requirement
 !!         
 !! IN
 !!  tot_size_arr=Number of scalar components of field and particule arrays
 subroutine print_procs_distrib() 

  use defs_mpitype

  character(len=500) :: msg

  __WRT_DEBUG_IN("print_procs_distrib")
  !****************** ecriture ecran voisange pickup ******************************  
  write(msg,'(4a)')ch10,&
       &' --------------------------',ch10,&
       &'   COORD PROC NON PICKUP'
  call wrt_double(6,msg,wrtscreen,wrtdisk)

  write(msg,'(a,i4.4,a,2i4)')&
       &' Process ',mpiinfo%me,&
       &' has coords ',mpiinfo%coord(1),mpiinfo%coord(2)
  call wrtout(6,msg,'PERS')

  write(msg,'(a)')' --------------------------'
  call wrt_double(6,msg,wrtscreen,wrtdisk)

  write(msg,'(4a)')ch10,&
       &' --------------------------',ch10,&
       &'   COORD PROC PICKUP'
  call wrt_double(6,msg,wrtscreen,wrtdisk)

  write(msg,'(a,i4.4,a,2i4)')&
       &' Process ',mpiinfo_pick%me,&
       &' has coords ',mpiinfo_pick%coord(1),mpiinfo_pick%coord(2)
  call wrtout(6,msg,'PERS')

  write(msg,'(a)')' --------------------------'
  call wrt_double(6,msg,wrtscreen,wrtdisk)

  if (wrtscreen == 0) then 
   write(msg,'(4a)')ch10,&
        &'====================================================' ,ch10,&
        &'  COORD ET VOISINAGE PROC DES PARTICULES STANDARD   '
   call wrt_double(6,msg,wrtscreen,wrtdisk)

   write(msg,'(a,i4.4,a,2i4,a,8i4)')&
        &' Process ',mpiinfo%me,&
        &' has coords ',mpiinfo%coord,&
        &' and neighbours ',mpiinfo%voisin 
   call wrtout(6,msg,'PERS')

   ! call MPI_BARRIER(mpiinfo%comm,ioerr)
   ! call MPI_BARRIER(mpiinfo_pick%comm,ioerr)
   write(msg,'(4a)')ch10,&
        &'====================================================' ,ch10,&
        &'  COORD ET VOISINAGE PROC DES DES PICKUPS   '
   call wrt_double(6,msg,wrtscreen,wrtdisk)


   write(msg,'(a,i4.4,a,2i4,a,8i4)')&
        &' Process ',mpiinfo_pick%me,&
        &' has coords ',mpiinfo_pick%coord,&
        &' and neighbours ',mpiinfo_pick%voisin
   call wrtout(6,msg,'PERS')
  endif !--wrtscreen

  __WRT_DEBUG_OUT("print_procs_distrib")
 end subroutine print_procs_distrib
 !********************* END PRINT_PROCS_DISTRIB ****************************

 subroutine init_temporal_B()

  use defs_parametre, only : pos_plan

  real :: x_gsm, y_gsm, z_gsm
  real :: x_sim, y_sim, z_sim
  character(len=150) :: file_name, fn
  logical :: L_open
  integer :: ix_box, iy_box, iz_box
  integer :: file_count, unit_file

  file_count = 0

  do ix_box = 1, int((s_max_loc(1) - s_min_loc(1))/gstep(1))
    do iy_box = 1, int((s_max_loc(2) - s_min_loc(2))/gstep(2))
      do iz_box = 1, int((s_max_loc(3) - s_min_loc(3))/gstep(3))

        x_sim = s_min_loc(1) + ix_box*gstep(1)
        y_sim = s_min_loc(2) + iy_box*gstep(2)
        z_sim = s_min_loc(3) + iz_box*gstep(3)

        x_gsm =  pos_plan(1)-x_sim
        y_gsm =  pos_plan(2)-y_sim
        z_gsm = -pos_plan(3)+z_sim
        
        if (((sqrt(x_gsm**2+y_gsm**2+z_gsm**2)<120) .and. &
               (mod(int(x_gsm), 10)==0 .and. &
                mod(int(y_gsm), 25)==0 .and. &
                mod(int(z_gsm), 25)==0)) .or. &
             (mod(int(x_gsm), 750)==0 .and. &
              mod(int(y_gsm), 100)==0 .and. &
              mod(int(z_gsm), 100)==0)) then

            write(file_name, '(A6,I5,A2,I5,A2,I5,A2,A4)') &
                              "Btime_",int(x_gsm),"x_", &
                                       int(y_gsm),"y_", &
                                       int(z_gsm),"z_",'.txt'
            !--By creating file_name and file_count at the same time
            !--we ensure they have a one-to-one correspondancy
            file_name = trim(file_name)
            file_count = file_count+1

            !--80+file_count because first units are already used
            open(UNIT=80+file_count,FILE=file_name)
            inquire(FILE=file_name, NUMBER=unit_file, OPENED=L_open)
!            write(*,*)  "file_name", file_name, &
!                      & "position", int(x_gsm), int(y_gsm), int(z_gsm), &
!                      & "isopen", L_open,"has number", unit_file
            write(UNIT=80+file_count,FMT=*) 'position:', &
                                           & int(x_gsm), int(y_gsm), int(z_gsm)
!            inquire(UNIT=80+file_count, NAME=FN, OPENED=L_open)
!            write(*,*)  "unit", 80+file_count, &
!                      & "position", int(x_gsm), int(y_gsm), int(z_gsm), &
!                      & "isopen", L_open,"has name", FN
            close(UNIT=80+file_count)

        endif

      enddo
    enddo
  enddo
 end subroutine init_temporal_B

end module initialisation

