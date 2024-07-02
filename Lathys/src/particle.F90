!!=============================================================
!!=============================================================
module particle

 use defs_basis
 use defs_particletype,only  : iXb,particletype
 use defs_counts_types,only  : count_single_type,&
      &                        count_particle_type
 use defs_species
 use defs_grid
 use defs_parametre
 use particle_com
 use environment
 use atm_ionosphere
 use m_writeout
 use m_timing,only           : time_get
 use diag_flux_part_imp
 use m_distribution_function

#include "q-p_common.h"

 implicit none
 private

 private::               &
      mvsp3r,            &
      vcalc3,            &
      sortie
 
 public::                &
      move,              &
      xcalc3,            &
      init_flux_part_imp

contains
 !!#####################################################################


 !********************************* sortie *******************************************
 subroutine sortie(s_p,s_min_loc,s_max_loc,&
      &            nid,exit_p,index_exit,countpar,&
      &            out_xm,out_xp,out_tot,particule)

  integer,intent(out):: exit_p
  integer,intent(in) :: nid!--Index of the particle
  integer,intent(inout) :: out_xp,out_xm !--Exiting particles in x-dir.
  integer,intent(inout) :: out_tot!--Exiting particles in all dir.
  integer,intent(inout) :: index_exit(:)!--Stocks Indices of exiting particles
  real(dp),intent(in) :: s_p(3),s_min_loc(3),s_max_loc(3)
  type(count_single_type),intent(inout) :: countpar
  type(particletype),intent(inout) :: particule(:)

  integer :: test_min_max(4) 
  integer :: test_res
  integer :: test_abnormal(2)

  exit_p = 0
  !--on liste les particules qui doivent sortir et on les compte
  ! les particules qui sortent dans la direction X ne sont pas envoyes vers des 
  ! processeurs voisins (la topologie de la grille est en Y et Z)

  !--First test for particles which go out of the simulation box in X-direction
  if(s_p(1)>s_max_loc(1))  then
   exit_p = -1
   out_xp = out_xp + 1
   out_tot = out_tot + 1
   index_exit(out_tot) = nid
   return
  endif

  if(s_p(1)<s_min_loc(1)) then
   exit_p = -1
   out_xm = out_xm + 1
   out_tot = out_tot + 1
   index_exit(out_tot) = nid
   return
  endif

  !--Now we test particles for going out particles in the Y-Z plane

  !--Get in a integer (test_res) where the particule have to go
  !  0=>no out; 1=>W; 2=>S; 4=>E; 8=>N; 3=>SW; 6=>SE; 9=>NW; 12=>NE  
  !  the rest is error
  test_min_max = 0
  test_abnormal = 0
  where(s_p(2:3) < s_min_loc(2:3))   test_min_max(1:2) = 1
  where(s_p(2:3) > s_max_loc(2:3))   test_min_max(3:4) = 1

  where(s_p(2:3) < (s_min_loc(2:3)-ncm(2:3)/gstep(2:3)))   test_abnormal(1:2) = test_abnormal(1:2)+1
  where(s_p(2:3) > (s_max_loc(2:3)+ncm(2:3)/gstep(2:3)))   test_abnormal(1:2) = test_abnormal(1:2)+1

  if (sum(test_abnormal).ne.0) then 
  ! THIS PARTICLE SHOULD NOT BE HERE, DROP IT OUT OF THE SIMULATIONS
   exit_p = -1
   out_xm = out_xm + 1
   out_tot = out_tot + 1
   index_exit(out_tot) = nid
#ifdef HAVE_DEBUG   
    print *,"Particle at abnormal position",nid!infompi%me,nn
    print *,"Particule%vel ",particule(nid)%vel
    print *,"Particule%pos ",particule(nid)%pos
    print *,"Particule%mass ",particule(nid)%mass
    print *,"Particule%char ",particule(nid)%char
    print *,"Particule%dir ",particule(nid)%dir
    print *,"Particule%orig ",particule(nid)%orig
    print *,"Particule%exc ",particule(nid)%exc
    print *,"ijk ",s_p
    print *,"smnl ",s_min_loc
    print *,"smxl ",s_max_loc
#endif
   return
  endif


  test_res = sum(test_min_max*(/1,2,4,8/))
  
  !--selecting 
  select case(test_res)
  case(0)
   return

  case(4) !--On copte et liste les particules qui sortent vers l'Est
   countpar%out_dir(E) = countpar%out_dir(E) + 1
   exit_p = E
 
  case(12) !--On compte et liste les particules qui sortent vers le nord-est
   countpar%out_dir(NE) = countpar%out_dir(NE) + 1
   exit_p = NE
 
  case(8) !--On compte et liste les particules qui sortent par le nord
   countpar%out_dir(N) = countpar%out_dir(N) + 1
   exit_p = N
 
  case(6) !--on compte et liste les particules qui sortent vers le sud-est
   countpar%out_dir(SE) = countpar%out_dir(SE) + 1
   exit_p = SE
 
  case(2) !--on compte et liste les particules qui sortent vers le sud
   countpar%out_dir(S) = countpar%out_dir(S) + 1
   exit_p = S
 
  case(3) !--on compte et liste les particules qui sortent vers le sud-ouest
   countpar%out_dir(SW) = countpar%out_dir(SW) + 1
   exit_p = SW
 
  case(1) !--on compte et liste les poarticules qui sortent vers l'ouest
   countpar%out_dir(W)  = countpar%out_dir(W) + 1
   exit_p = W
 
  case(9) !--on compte et liste les particules qui sortent vers le nord-ouest
   countpar%out_dir(NW) = countpar%out_dir(NW) + 1
   exit_p = NW
 
  case default
   stop "error in sortie counting outer particles"
  end select

  out_tot = out_tot+1
  index_exit(out_tot) = nid
  
 end subroutine sortie
 !********************************** sortie *******************************************
 !*******************************Init impacting particles flux**********
 subroutine init_flux_part_imp()
   use defs_variable
   use defs_mpitype
   use defs_parametre,only:planetname

   integer :: i,j

   allocate(esurf_grid(n_inte)) ; allocate(asurf_grid(n_inta))
   allocate(phi_grid(nphi)) ; allocate(theta_grid(ntheta))
   allocate(flux_surf(nphi,ntheta,1))
   allocate(fdde_surf(nphi,ntheta,n_inte,1))
   allocate(fdda_surf(nphi,ntheta,n_inta,1))

   esurf_grid(:)=0.d0 ; asurf_grid(:)=0.d0
   phi_grid(:)=0.d0 ; theta_grid(:)=0.d0
      flux_surf(:,:,:) = 0.d0
      do j=1,n_inte
         fdde_surf(:,:,j,:) = 0.d0
      enddo
      do j=1,n_inta
         fdda_surf(:,:,j,:) = 0.d0
      enddo
   

   !--- Phi grid
   do i=1,nphi
      phi_grid(i) = (i - 1.d0) * 2.d0 * pi / (nphi-1)
   enddo

   !--- Theta grid
   do i=1,ntheta
      theta_grid(i) = (i - 1.d0) * pi / (ntheta-1)
   enddo

   !--- Energy intervals
   esurf_grid(1) = ener_min
   dE=(ener_max-ener_min)/(n_inte-1.d0)
   do i = 2 , n_inte
      esurf_grid(i) = esurf_grid(i-1)*(1 + dE_E)
   enddo
   ener_max = esurf_grid(n_inte)
   !--- Angle intervals
   do i = 1 , n_inta
      asurf_grid(i) = 0.5d0*pi+0.5d0*pi * (i-1.d0) / (n_inta - 1.d0)
   enddo


 end subroutine init_flux_part_imp
!***************************************************************************

! Compute surface of the cell of spherical grid for flux of impacting particles
function surf(ind_phi,ind_theta)
use defs_variable

integer,intent(in) :: ind_phi,ind_theta
real(dp)           :: surf

surf = ((Spe%P%radius+alt_part_imp)*Spe%ref%c_omegapi)*((Spe%P%radius*alt_part_imp)*Spe%ref%c_omegapi)*&
   &(cos(theta_grid(ind_theta))-cos(theta_grid(ind_theta+1)))*&
   &(phi_grid(ind_phi+1)-phi_grid(ind_phi))*1.e+10 ! in cm2
return
end function

 !*<xcalc3>**************************************************************
 !--Position advance (3-D) plus periodic boundary conditions.
 subroutine xcalc3(dt,particule,n1,n2,s_min_loc,s_max_loc, &
      &            Preg,index_exit,Spe)
 
  use defs_mpitype,only     : mpiinfo,mpiinfo_pick
  use defs_counts_types
  use particle_creation,only         : ntot_entree
  !use defs_variable
  integer,intent(in) :: n1,n2
  integer,intent(inout) :: index_exit(:)
  real(dp),intent(in) :: dt
  real(dp),intent(in),dimension(3) :: s_min_loc,s_max_loc
  type(species_type),intent(in) :: Spe
  type(count_particle_type),intent(out) :: Preg
  type(particletype),intent(inout) :: particule(:)

  integer :: nn,exit_p,n_abs
  real(dp) :: dst_planet,dst_planet_av
  real(dp),dimension(3) :: s_p,s_temp,s_temp_av
  character(len=500) :: msg
  character(len=15) :: nom_fichier
#ifdef HAVE_DEBUG
  integer  :: compt_sortie_mask
#endif
  integer :: i_theta,i_phi,i_ang,i_ener,ioerr,i,j
  real(dp) :: ang0,w_p,vol_part,N_ions,Vr0,ener,theta0,phi0,r0
  real(dp) :: V0(3),xp(3),test


  __WRT_DEBUG_IN("xcalc3")

  call init_count_particle_type(Preg)

  do nn = n1,n2
   s_p = particule(nn)%pos+dt*particule(nn)%vel
   exit_p = int(particule(nn)%dir)

   !--On elimine les particules les plus energetiques
   !--Apply PERIODIC boundary conditions

   !--Macroparticules provenant 
   !---   du vent solaire                              typ(n) = 0
   !---   de la photoionisation                        typ(n) = 1
   !---   de l'ionisation par collisions electroniques typ(n) = 2
   !---   de l'echappement ionosphrique                typ(n) = 3
   !---   de la zone initialement liee a la planete    typ(n) = 4
   !--Frontiere planetaire: absorbing particles
#ifndef HAVE_NO_PLANET
   s_temp = s_p-Spe%P%centr
   s_temp_av = particule(nn)%pos-Spe%P%centr   
   dst_planet = sqrt(dot_product(s_temp,s_temp))
   dst_planet_av = sqrt(dot_product(s_temp_av,s_temp_av))   
   if (dst_planet <= Spe%P%r_lim) then

 ! ecriture des particules impactant la planete
!    if (diag_part_imp.ne.0) then
! !    if ((iter+1).ge.4000) then
!      if (dst_planet_av > Spe%P%r_lim) then
!            n_abs = n_abs+1
!                !print*,'-------PARTIMP--------',n_abs
!            write(mpiinfo%me+20,*) iter+1,n_abs,s_p,particule(nn)%vel,particule(nn)%char,particule(nn)%mass,particule(nn)%tbirth
 !               !,particule(nn)%orig,particule(nn)%exc
 !     endif
!!     endif
 !   endif
 !  
    if (iter.eq.(0)) then
     Preg%kabs = Preg%kabs +1
     particule(nn)%vel(1)  = Spe%P%speed
     particule(nn)%vel(2:3) = zero
     particule(nn)%orig = 4_iXb
    else
     if ((particule(nn)%orig.ne.4).and.(absorp_surface.ne.1)) then
      s_p(1)  = s_max_loc(1)+0.1_dp
     else
      particule(nn)%vel(1)  = Spe%P%speed
      particule(nn)%vel(2:3) = zero
     endif
    endif
   endif

  ! If the particle impacts the surface
   if(dst_planet < Spe%P%r_lim + alt_part_imp/Spe%ref%c_omegapi) then
     if (((diag_part_imp.ne.0).and.(abs(particule(nn)%char/particule(nn)%mass-1._dp/16._dp) < 0.01)).and. &
            & (iter*dt > 300.)) then
        ! Compute the flux of impacting particles
        if (dst_planet_av >= Spe%P%r_lim + alt_part_imp/Spe%ref%c_omegapi .and. particule(nn)%orig/= 4) then
            n_part_imp=n_part_imp+1
            !-- Position of the particle in spherical coordinates
            xp(1) = (particule(nn)%pos(1) - Spe%P%centr(1))*Spe%ref%c_omegapi
            xp(2) = (particule(nn)%pos(2) - Spe%P%centr(2))*Spe%ref%c_omegapi
            xp(3) = (particule(nn)%pos(3) - Spe%P%centr(3))*Spe%ref%c_omegapi
            r0 = sqrt(dot_product(xp,xp))
            phi0   = acos(xp(1)/sqrt(xp(1)*xp(1)+xp(2)*xp(2))) 
            if(xp(2)<0.d0) phi0 = 2.d0*pi-phi0
            theta0 = acos(xp(3)/r0)
            V0(:) = particule(nn)%vel * Spe%ref%alfvenspeed                                !!--- Velocity of the particle in km/s
            Vr0   = V0(1)*cos(phi0)*sin(theta0) + V0(2)*sin(theta0)*sin(phi0) +&
                &  V0(3)*cos(theta0)  !!--- Radial velocity km/s
            ener  = 0.5d0  * particule(nn)%mass/particule(nn)%char * &
                & amu_pmass * dot_product(1.d3*V0,1.d3*V0) * JtoeV !!---Energy in eV
            if(Vr0<0.d0 .and. ener<ener_max.and.ener>ener_min) then
            
               w_p=Spe%ref%density * particule(nn)%char
            n_abs=n_abs+1
            !--- Compute the weight of the particle (=to a density in cm-3)
            vol_part=(gstep(1)*Spe%ref%c_omegapi*1d5)**3 !!-- in cm3
            N_ions=w_p*vol_part                          !!-- Physical number of
                                                         !!-- ions in the
                                                         !!-- macro-particule

            !--- Add the contribution of this particle to compute the flux of impacting particles
            i_phi   = int(phi0/(2*pi)*(nphi-1))+1
            if(i_phi==nphi) i_phi=i_phi-1
            i_theta = int(theta0/(pi)*(ntheta-1))+1
            if(i_theta==ntheta) i_theta=i_theta-1           
            if(i_phi<=0 .or. i_theta <=0 .or. i_phi>=nphi .or. i_theta>=ntheta) then
               print*,'ERREUR',i_phi,i_theta
               stop 
            endif
            if(phi0<phi_grid(i_phi) .or. phi0>phi_grid(i_phi+1)) then
             !  i_phi = nphi
               print*,'ERREUR_PHI',phi0,phi_grid(i_phi),phi_grid(i_phi+1)
               !stop
            endif
            if(theta0<theta_grid(i_theta) .or. theta0>theta_grid(i_theta+1)) then
               print*,'ERREUR_THETA',theta0,theta_grid(i_theta),theta_grid(i_theta+1)
              ! stop
            endif
            flux_surf(i_phi,i_theta,1) = flux_surf(i_phi,i_theta,1) + &
                & N_ions/surf(i_phi,i_theta)   !!--- Has to be divided by the surface of 
                                                                                                             !!--- the cell.
           ! i_ener = floor((ener-ener_min)/E)+1
            i_ener = floor(log(ener/ener_min)/log(1+dE_E))+1
            if( ener>esurf_grid(i_ener+1) .or. i_ener<1) then
               
                
                print*,'ERREUR_ENER',ener,esurf_grid(i_ener),esurf_grid(i_ener+1),i_ener,ener_min,dE
                !stop
            endif
            fdde_surf(i_phi,i_theta,i_ener,1) = fdde_surf(i_phi,i_theta,i_ener,1) + N_ions/surf(i_phi,i_theta)
            !--- Angle distribution
            ang0 = acos(Vr0/sqrt(dot_product(V0,V0)))
            if(ang0 > pi  ) ang0 = pi * 0.99999999d0
            if(ang0 < 0.d0) ang0 =      0.00000001d0
            i_ang = floor((ang0-pi*0.5d0)*2.d0/pi * (n_inta-1))+ 1

            if(ang0<asurf_grid(i_ang) .or. ang0>asurf_grid(i_ang+1) .or. i_ang<1 .or. ang0<pi*0.5d0) then
               print*,'ERREUR_ANGLE',ang0,asurf_grid(i_ang),asurf_grid(i_ang+1)
               print*,particule(nn)%vel*Spe%ref%alfvenspeed
               print*,particule(nn)%pos-s_p
               print*,dst_planet,dst_planet_av
               print*,Spe%P%r_lim + alt_part_imp/Spe%ref%c_omegapi
               print*,Vr0
               stop
            endif
            fdda_surf(i_phi,i_theta,i_ang,1) = fdda_surf(i_phi,i_theta,i_ang,1) + N_ions/surf(i_phi,i_theta)

               if(minval(fdde_surf)<0.d0 .or. minval(fdda_surf)<0.d0) then
                   print*,'ERREUR NEGATIF',mpiinfo%me,N_ions,surf(i_phi,i_theta)
                   stop
               endif
           endif
         endif
      endif
   endif

#endif

   if((particule(nn)%exc == zero).and.(particule(nn)%orig==0_iXb)) then
    !--Normal particles
    call sortie(s_p,s_min_loc,s_max_loc,&
         &      nn,exit_p,index_exit,&
         &      Preg%np,Preg%out_xm,&
         &      Preg%out_xp,Preg%out_tot,particule)
   else
    !--Pickups 
    call sortie(s_p,s_min_loc,s_max_loc,&
         &      nn,exit_p,index_exit,&
         &      Preg%pp,Preg%out_xm,&
         &      Preg%out_xp,Preg%out_tot,particule)
   endif


   particule(nn)%pos = s_p
   particule(nn)%dir = int(exit_p,iXb) 
  end do

  !--Compute the total number of exiting particles and pickups
  Preg%np%out = sum(Preg%np%out_dir)
  Preg%pp%out = sum(Preg%pp%out_dir)

  !--The total number of exiting particles is given by the sum of
  !  exiting particles and pickups

#ifdef HAVE_DEBUG
  compt_sortie_mask = count((&
       ((particule(n1:n2)%pos(1)>s_max_loc(1)).or.&
       &(particule(n1:n2)%pos(1)<s_min_loc(1))) .or.&
       ((particule(n1:n2)%pos(2)>s_max_loc(2)).or.&
       &(particule(n1:n2)%pos(2)<s_min_loc(2))))  .or.&
       ((particule(n1:n2)%pos(3)>s_max_loc(3)).or.&
       &(particule(n1:n2)%pos(3)<s_min_loc(3))))
  
  if(Preg%out_tot /= compt_sortie_mask) then
   print *,"Probleme particule"
   print *,"counted particles",Preg%out_tot
   print *,"masked particles",compt_sortie_mask
   stop
  endif
 
  write(msg,'(2a,i4.4,a,8i5)')&
       &" EXC PARTICULES =>",&
       &" Proc ",mpiinfo%me,&
       &" is exchanging ",Preg%np%out_dir
  call wrtout(6,msg,"PERS")

  write(msg,'(2a,i4.4,a,8i5)')&
       &" EXC PICKUPS =>",&
       &" Process ",mpiinfo_pick%me,&
       &" is exchanging ",Preg%pp%out_dir
  call wrtout(6,msg,"PERS")
#endif

  write(msg,'(2a,11(2a,i15))')&
       &ch10," _____ Sortie dans xcalc3 (proc 0) ____",&
       &ch10,"   ntot_in  = ", ntot_entree,&
       &ch10,"   ntot_out = ", Preg%out_xm+Preg%out_xp,&
       &ch10,"   out-in   = ", Preg%out_xm+Preg%out_xp-ntot_entree,&
       &ch10,"   ksr      = ", Preg%out_xm,&
       &ch10,"   ksl      = ", Preg%out_xp,&
       &ch10,"   ksd      = ", Preg%np%out_dir(N)+Preg%np%out_dir(NW),&
       &ch10,"   ksd      = ", Preg%np%out_dir(S)+Preg%np%out_dir(SE),&
       &ch10,"   ksf      = ", Preg%np%out_dir(W)+Preg%np%out_dir(SW),&
       &ch10,"   ksb      = ", Preg%np%out_dir(E)+Preg%np%out_dir(NE),&
       &ch10,"   absorbed = ", Preg%kabs,&
       &ch10,"   tot out  = ", Preg%out_tot
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  !!------ Dumping the sputtering flux every XXX time steps
  if (((diag_part_imp.ne.0).and.(mod(iter,300)==0)).and.(iter*dt > 300)) then
   write(nom_fichier,'(a8,a2,i5.5)')trim(fildat),'_t',int(iter*dt)
     
     call wrt_flux_part_imp(nom_fichier)
  endif

  __WRT_DEBUG_OUT("xcalc3")
 end subroutine xcalc3
 !********************************** XCALC3 ******************************************

 !******************************* MOVE ******************************************
 subroutine move(Preg,index_exit)

  !--Particle push : x(new) = x + dt * vh
  ! Calculation of velocities from Lorentz force : dv/dt = (q/m)(E + v x B)
  ! "Raw" moment collection (no normalisation or smoothing)
  ! of v(1) at x(1/2) [ub] and x(3/2) [uf], as well as arrays for
  ! moment advance at x(3/2) [u]

  use field_pe
  use field_lissage,only  : smth_func
  use field_e,only        : ecalc3
  use defs_mpitype,only   : mpiinfo
  use defs_parametre,only : gstep,resis,dmin
  use defs_variable,only  : s_min_loc,&
       &                    Efield,Bfield,&
       &                    te,pe,&
       &                    dn,dnf,dna,velf,&
       &                    resistivity,&
       &                    e_conv,e1_conv,&
       &                    vel,vela,rmu0,&
       &                    particule,nptot,iwr

  integer,intent(inout) :: index_exit(:)
  type(count_particle_type),intent(out) :: Preg

  !*************************************************
  !--IMPORTANT: E = E(vel%coor(1),vel%coor(2),vel%coor(3)) : ADVANCED CURRENTS *
  !*************************************************
  __WRT_DEBUG_IN("move")
  __GETTIME(20,1)!--Timer start

  call pecalc(ipe,te,dnf,pe)
 
  __GETTIME(22,1) !--Timer start ecalc3
  call ecalc3(Efield,Bfield,vel,dnf,pe,resistivity,e_conv,&
       &e1_conv,gstep,dmin,rmu0,resis,idisp,nc1,mpiinfo)
  __GETTIME(22,2)!--Timer stop ecalc3

  call mvsp3r(dt,particule,Efield,Bfield,dna,vela,dnf,velf,dn,vel, &
       1,nptot,gstep,s_min_loc,resis,resistivity,Preg,index_exit)

  !--Multiply by the species charge, and q**two/m for charge advance
  !--4-current at xp (1/2)   (dna,vela)
  !--4-current at xnew (3/2) (dnf,velf)
  !--Arrays for current advance
  !--RAW MOMENTS HAVE BEEN COLLECTED WITH NO PROCESSING
  !--"A(1)" moments are the average of "B(1/2)" and "F(3/2)" moments
  dna    = half*(dna    + dnf)
  vela%x = half*(vela%x + velf%x)
  vela%y = half*(vela%y + velf%y)
  vela%z = half*(vela%z + velf%z)

  !--Smoothing (optional)
  if(iwr(6) == 1)then
   __GETTIME(24,1)  !--Timer start smth_func
   call smth_func(dnf,   nc,mpiinfo)
   call smth_func(dna,   nc,mpiinfo)
   call smth_func(vela%x,nc,mpiinfo)
   call smth_func(vela%y,nc,mpiinfo)
   call smth_func(vela%z,nc,mpiinfo)
   __GETTIME(24,2)  !--Timer stop smth_func
  end if

  __GETTIME(20,2)!--Timer stop
  __WRT_DEBUG_OUT("move")
 end subroutine move
 !************************************** END MOVE ***********************************

!************************************ MVSP3R ***************************************
 subroutine mvsp3r(dt,particule,Efield,Bfield, &
      &            dna,vela,dnf,velf,dn,vel,n1,n2,gstep, &
      &            s_min_loc,resis,resistivity,&
      &            Preg,index_exit)

  !-- (trilineaire moment collection, identical E and B grids)
  !      include "h3_param.inc"

  use defs_arr3Dtype
  use defs_mpitype,only       : mpiinfo
  use defs_variable,only      : s_max_loc,s_min,s_max,s_r,Spe,dn_e_pl,dn_e_incdt
  use part_moment,only             : Bmtsp3,Fmtsp3
  use particle_creation,only  : new_particles
  use particle_sort

  integer,intent(in) :: n1
  integer,intent(inout) :: n2
  integer,intent(inout) :: index_exit(:)
  real(dp),intent(in) :: dt,resis
  real(dp),dimension(3),intent(in) :: gstep,s_min_loc
  type(count_particle_type),intent(out) :: Preg
  real(dp),dimension(:,:,:),intent(inout) :: resistivity
  !--resistivity is intent(inout) because is not used 
  !  if used put it intent(in)
  real(dp),dimension(:,:,:),intent(inout) :: dna,dnf,dn
  type(arr3Dtype),intent(in) :: Efield,Bfield
  type(arr3Dtype),intent(inout) ::vela,velf,vel
  type(particletype),intent(inout):: particule(:)

  __WRT_DEBUG_IN("mvsp3r")
  __GETTIME(23,1) !--Timer start 

  call vcalc3(n1,n2,dt,gstep,s_min_loc,resis,Spe,resistivity,Efield,Bfield,particule)

  __GETTIME(28,1) !--Timer start Fmtsp3
  call Fmtsp3(particule,dna,dn_e_pl,dn_e_incdt,vela,n1,n2,gstep,s_min_loc,mpiinfo,nc1)
  __GETTIME(28,2) !--Timer stop Fmtsp3
  
  __GETTIME(29,1) !--Timer start xcalc3
  call xcalc3(dt,particule,n1,n2,s_min_loc,s_max_loc,Preg,index_exit,Spe)
  __GETTIME(29,2) !--Timer stop xcalc3

  !--Exchange particles
  __GETTIME(30,1) !--Timer start pack_com
  call pack_com(particule,Preg,index_exit,s_min,s_r,n2,s_min_loc,s_max_loc)
  __GETTIME(30,2) !--Timer stop pack_com

  !--Exit and creation of new particles
  __GETTIME(31,1) !--Timer start new_particle
  call new_particles(2,particule)
  __GETTIME(31,2) !--Timer stop new_particles

  __GETTIME(33,1) !--Timer start new_particle
  call p_sort(particule,n2)
  __GETTIME(33,1) !--Timer start new_particle

  __GETTIME(32,1) !--Timer start bmtsp3
  call Bmtsp3(particule,dnf,velf,dn,vel,n1,n2,gstep,s_min_loc,mpiinfo,nc1)
  __GETTIME(32,2) !--Timer stop bmtsp3
 
  __GETTIME(23,2) !--Timer stop
  __WRT_DEBUG_OUT("mvsp3r")
 end subroutine mvsp3r
 !******************************* MVSP3R *********************************************


 !******************************* VCALC3 *********************************************
 subroutine vcalc3(n1,n2,dt,gstep,s_min_loc,resis,Spe, &
      &            resistivity,Efield,Bfield,particule)

  use defs_arr3Dtype
  use defs_variable,only      : s_max_loc,viscosity,t
  use m_rand_gen,only           : rand_gen1,irand

  integer,intent(in) :: n1
  integer,intent(inout) :: n2
  real(dp),intent(in) :: dt,resis
  real(dp),intent(in),dimension(3) :: gstep,s_min_loc
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(arr3Dtype),intent(in) :: Efield,Bfield
  real(dp),intent(inout),dimension(:,:,:) :: resistivity 
  !--resistivity is intent(inout) because is not used 
  !  if used put it intent(in)

  integer :: kpickup,nn,ii,nptot,knew
  integer :: ijk(3),ijk_old(3)
  real(dp) :: sm,sq,qsm,dta,dth
  real(dp) :: w(8)
  real(dp) :: dist2,viscloc
  logical :: part_is_out
  character(len=200) :: msg
  real(dp),dimension(3) ::  s_p,s_e,s_f,s_a,s_b
  real(dp),dimension(3) ::  s_tmp,v_p,v_h
  real(dp),dimension(3) ::  e_local,b_local
  real(dp),dimension(3) ::  gstep_inv,f_lorentz
  real(dp),dimension(8,3) :: b_eight

  __WRT_DEBUG_IN("vcalc3")
  __GETTIME(27,1)!--Timer start

  gstep_inv = one/gstep
  ijk_old=(/-1,-1,-1/)
  kpickup  = 0

  nptot=n2
  do  nn = n1,n2
   !--Assignments
   s_p = particule(nn)%pos
   v_p = particule(nn)%vel
   sm  = particule(nn)%mass
   sq  = particule(nn)%char

   !--Distance from the planet
   s_tmp = s_p-Spe%P%centr
   dist2 = (dot_product(s_tmp,s_tmp))
   part_is_out=(dist2.gt.Spe%P%r_lim**2)

#ifndef HAVE_NO_PLANET
   !--No advance for particles absorbed
   if(.not.part_is_out) then
    particule(nn)%vel = (/Spe%P%speed,zero,zero/)
    cycle
   endif
#endif

   qsm = sq/sm
   dta = qsm*dt
   dth = half*dta

   !------------------
   !--Velocity advance
   !------------------
   !--E and B fields interpolated from same grid
   !--Relative position in cell    i = int(s_e(1));s_f=(xf,yf,zf) center of the particule
   s_e = (1.5_dp + ((s_p-s_min_loc)*gstep_inv))
   ijk = int(s_e)
   s_f = s_e-real(ijk,dp)

   !--Sequence of indices at cell corners
   !-- (i,j,k)
   !--Trilinear weights
   s_a = 1._dp - s_f
   
   w(1) = s_a(1)*s_a(2)*s_a(3) !xa*ya*za
   w(2) = s_f(1)*s_a(2)*s_a(3) !xf*ya*za
   w(3) = s_a(1)*s_f(2)*s_a(3) !xa*yf*za
   w(4) = s_f(1)*s_f(2)*s_a(3) !xf*yf*za
   w(5) = s_a(1)*s_a(2)*s_f(3)
   w(6) = s_f(1)*s_a(2)*s_f(3)
   w(7) = s_a(1)*s_f(2)*s_f(3)
   w(8) = s_f(1)*s_f(2)*s_f(3)
   
   !--Linear interpolation
!-- do not try to optimize E as B, because of the different grid, it doen't work
   e_local(1) = w(1)*Efield%x(ijk(1)  ,ijk(2)  ,ijk(3)  )&
        &     + w(2)*Efield%x(ijk(1)+1,ijk(2)  ,ijk(3)  )&
        &     + w(3)*Efield%x(ijk(1)  ,ijk(2)+1,ijk(3)  )&
        &     + w(4)*Efield%x(ijk(1)+1,ijk(2)+1,ijk(3)  )&
        &     + w(5)*Efield%x(ijk(1)  ,ijk(2)  ,ijk(3)+1)&
        &     + w(6)*Efield%x(ijk(1)+1,ijk(2)  ,ijk(3)+1)&
        &     + w(7)*Efield%x(ijk(1)  ,ijk(2)+1,ijk(3)+1)&
        &     + w(8)*Efield%x(ijk(1)+1,ijk(2)+1,ijk(3)+1)

   e_local(2) = w(1)*Efield%y(ijk(1)  ,ijk(2)  ,ijk(3)  )&
        &     + w(2)*Efield%y(ijk(1)+1,ijk(2)  ,ijk(3)  )&
        &     + w(3)*Efield%y(ijk(1)  ,ijk(2)+1,ijk(3)  )&
        &     + w(4)*Efield%y(ijk(1)+1,ijk(2)+1,ijk(3)  )&
        &     + w(5)*Efield%y(ijk(1)  ,ijk(2)  ,ijk(3)+1)&
        &     + w(6)*Efield%y(ijk(1)+1,ijk(2)  ,ijk(3)+1)&
        &     + w(7)*Efield%y(ijk(1)  ,ijk(2)+1,ijk(3)+1)&
        &     + w(8)*Efield%y(ijk(1)+1,ijk(2)+1,ijk(3)+1)

   e_local(3) = w(1)*Efield%z(ijk(1)  ,ijk(2)  ,ijk(3)  )&
        &     + w(2)*Efield%z(ijk(1)+1,ijk(2)  ,ijk(3)  )&
        &     + w(3)*Efield%z(ijk(1)  ,ijk(2)+1,ijk(3)  )&
        &     + w(4)*Efield%z(ijk(1)+1,ijk(2)+1,ijk(3)  )&
        &     + w(5)*Efield%z(ijk(1)  ,ijk(2)  ,ijk(3)+1)&
        &     + w(6)*Efield%z(ijk(1)+1,ijk(2)  ,ijk(3)+1)&
        &     + w(7)*Efield%z(ijk(1)  ,ijk(2)+1,ijk(3)+1)&
        &     + w(8)*Efield%z(ijk(1)+1,ijk(2)+1,ijk(3)+1)

   s_b = (s_e - half)
   ijk = int(s_b)
   s_f = s_b-real(ijk,dp)
   if (sum(abs(ijk-ijk_old)).ne.0) then
   b_eight(:,1) = (/&
        &         Bfield%x(ijk(1)  ,ijk(2)  ,ijk(3)  ),&
        &         Bfield%x(ijk(1)+1,ijk(2)  ,ijk(3)  ),&
        &         Bfield%x(ijk(1)  ,ijk(2)+1,ijk(3)  ),&
        &         Bfield%x(ijk(1)+1,ijk(2)+1,ijk(3)  ),& 
        &         Bfield%x(ijk(1)  ,ijk(2)  ,ijk(3)+1),&
        &         Bfield%x(ijk(1)+1,ijk(2)  ,ijk(3)+1),&
        &         Bfield%x(ijk(1)  ,ijk(2)+1,ijk(3)+1),&
        &         Bfield%x(ijk(1)+1,ijk(2)+1,ijk(3)+1) /)
   b_eight(:,2) = (/&
        &         Bfield%y(ijk(1)  ,ijk(2)  ,ijk(3)  ),&
        &         Bfield%y(ijk(1)+1,ijk(2)  ,ijk(3)  ),&
        &         Bfield%y(ijk(1)  ,ijk(2)+1,ijk(3)  ),&
        &         Bfield%y(ijk(1)+1,ijk(2)+1,ijk(3)  ),& 
        &         Bfield%y(ijk(1)  ,ijk(2)  ,ijk(3)+1),&
        &         Bfield%y(ijk(1)+1,ijk(2)  ,ijk(3)+1),&
        &         Bfield%y(ijk(1)  ,ijk(2)+1,ijk(3)+1),&
        &         Bfield%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) /)
   b_eight(:,3) = (/&
        &         Bfield%z(ijk(1)  ,ijk(2)  ,ijk(3)  ),&
        &         Bfield%z(ijk(1)+1,ijk(2)  ,ijk(3)  ),&
        &         Bfield%z(ijk(1)  ,ijk(2)+1,ijk(3)  ),&
        &         Bfield%z(ijk(1)+1,ijk(2)+1,ijk(3)  ),& 
        &         Bfield%z(ijk(1)  ,ijk(2)  ,ijk(3)+1),&
        &         Bfield%z(ijk(1)+1,ijk(2)  ,ijk(3)+1),&
        &         Bfield%z(ijk(1)  ,ijk(2)+1,ijk(3)+1),&
        &         Bfield%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) /)
   endif
   ijk_old=ijk
   b_local =  matmul(w,b_eight)
   !--Bilinear weights

   viscloc =viscosity(ijk(1),ijk(2),ijk(3))*dta
   if (viscloc.gt.(1e-6)) then
   s_a = one-s_f
    w(1) = s_a(1)*s_a(2)*s_a(3)
    w(2) = s_f(1)*s_a(2)*s_a(3)
    w(3) = s_a(1)*s_f(2)*s_a(3)
    w(4) = s_f(1)*s_f(2)*s_a(3)
    w(5) = s_a(1)*s_a(2)*s_f(3)
    w(6) = s_f(1)*s_a(2)*s_f(3)
    w(7) = s_a(1)*s_f(2)*s_f(3)
    w(8) = s_f(1)*s_f(2)*s_f(3)

   viscloc = (w(1)*viscosity(ijk(1)  ,ijk(2)  ,ijk(3)  ) &
        &     +w(2)*viscosity(ijk(1)+1,ijk(2)  ,ijk(3)  ) &
        &     +w(3)*viscosity(ijk(1)  ,ijk(2)+1,ijk(3)  ) &
        &     +w(4)*viscosity(ijk(1)+1,ijk(2)+1,ijk(3)  ) &
        &     +w(5)*viscosity(ijk(1)  ,ijk(2)  ,ijk(3)+1) &
        &     +w(6)*viscosity(ijk(1)+1,ijk(2)  ,ijk(3)+1) &
        &     +w(7)*viscosity(ijk(1)  ,ijk(2)+1,ijk(3)+1) &
        &     +w(8)*viscosity(ijk(1)+1,ijk(2)+1,ijk(3)+1) ) 
   viscloc=exp(-viscloc*dta)
   else
   viscloc=1._dp
   endif
   
   !--vh = v + (dt/2)(E + v_p x B)-------
   !--v_p x B
   f_lorentz = (/v_p(2)*b_local(3)-v_p(3)*b_local(2),&
        &        v_p(3)*b_local(1)-v_p(1)*b_local(3),&
        &        v_p(1)*b_local(2)-v_p(2)*b_local(1) /)
   v_h = v_p + dth*(e_local + f_lorentz)
   !--v(new) = v + dt*(E + vh x B)-------
   !--v_h x B
   f_lorentz = (/v_h(2)*b_local(3)-v_h(3)*b_local(2),&
        &        v_h(3)*b_local(1)-v_h(1)*b_local(3),&
        &        v_h(1)*b_local(2)-v_h(2)*b_local(1) /)

   v_p = v_p + dta*(e_local + f_lorentz)
   v_p = v_p*viscloc

        !Check for Nan
        if ((sum(v_p)-sum(v_p)).lt.0.001) then 
        v_p=v_p
        else
        v_p=particule(nn)%vel
        print *,'Vcalc found NaN'
        endif

   !l'ionosphere est bloquee a l'initialisation
   if ((t.gt.t_iono_release).or.(particule(nn)%orig.ne.3))  particule(nn)%vel = v_p 
 
   !******************************* ECHANGE DE CHARGE *********  
   if (part_is_out) then
    call calc_charge_exchange(nn,kpickup,qsm,irand,ijk,v_p,w,Spe,particule,atmosphere)
    call split_particle_iono(Spe,particule,nn,nptot,gstep,dist2,irand, &
        &       s_min_loc,s_max_loc)
   endif
   
   !-- Compute distribution function
   if (distrib_activated .eq. 1) call compute_distribution_function(s_p,v_p,sm,sq,Spe)
  enddo

  knew=nptot-n2
  n2=nptot
  call feed_ionosphere(Spe,particule,gstep,s_min_loc,s_max_loc,irand,nptot,atmosphere)  
  write(msg,'(2a,3(2a,i5),2a,i10)')&
       & ch10," ******************************",&
       & ch10," nombre de charges echangees ",kpickup,&
       & ch10," nombre de pick-up a creer   ",nptot-n2,&
       & ch10," nombre de pick-up a spliter ",knew,&
       & ch10," nombre de particules        ",nptot
  call wrtout(6,msg,'PERS')
  n2 = nptot  

  __GETTIME(27,2)!--Timer stop
  __WRT_DEBUG_OUT("vcalc3")
 end subroutine vcalc3
 !***************************** END VCALC3 **************************************

end module particle
