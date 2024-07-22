!!=============================================================
!!=============================================================
!!module: cration
!! NAME
!!  cration (RModolo,MMancini)
!!
!! FUNCTION
!!  This module contains creation of particles
!!  and Ionization
!!  
module particle_creation

 use defs_basis
 use defs_mpitype
 use defs_particletype
 use defs_variable,only   : dn,nptot,     &
      &                     s_min_loc,s_max_loc,Spe,t
 use defs_grid
 use defs_parametre,only       : dt,npm,nfl,gstep
 use m_rand_gen,only        : rand_gen1,irand
 use m_writeout
 use m_timing
#ifdef HAVE_WAVE_TEST
 use wavetest,only      : by0_b,bz0_b
#endif
 
#include "q-p_common.h"

 implicit none
 private

 integer,public,protected :: ntot_entree !--Total Entering particles 
 
 public ::                &
      new_particles      !--Create new particles

 
contains
 !!#####################################################################
 
 !!=============================================================
 !!routine: creation/new_particles
 !!
 !! FUNCTION
 !!  Create new particles
 !! IN 
 !! OUT
 !!   
 !! SIDE EFFECT
 !!  particles
 subroutine new_particles(kval,particule)
 use defs_parametre,only : ROI

  integer,intent(in) :: kval
  type(particletype),intent(inout) :: particule(:)

  integer,save :: irand_i = 0
  integer :: nlast,ie,ip_inj,ind_vel,n,is
  real(dp) :: dtp,test_0  !--Pas de temps temporaire
  real(dp) :: entree,dec_entree,vmag,theta,test_spec,coupdslo
  real(dp) :: sin_cos(2)
  real(dp) :: mul,imul
#ifdef HAVE_DEBUG
  character(len=500) :: msg
#endif

  __WRT_DEBUG_IN("new_particles")
  !--assumption : the new particles are created on the left

  mul=1._dp
  if ((nc(2)*nc(3)).lt.(nc_tot(2)*nc_tot(3)/real(nproc,dp))) mul=(1._dp+ROI%excess_part)
  imul=1._dp/mul

  nlast = nptot

  if(kval==1) dtp = half*dt
  if(kval==2) dtp = dt

  entree = sum(real(Spe%S(:)%ng,dp)*Spe%S(:)%flux_tot)
!  entree = entree*real(nc_tot(2)*nc_tot(3),dp)*dtp/(gstep(1)*real(nproc,dp)) 
  entree = entree*real(nc(2)*nc(3),dp)*dtp/(gstep(1))*mul 
  ntot_entree = int(entree)
  dec_entree = entree - real(ntot_entree,dp)

  if(rand_gen1(irand) <= dec_entree) ntot_entree = ntot_entree+1 

  ip_inj = 0
  Spe%S(:)%ip_inj = 0

  !--On remplace les particules qui sortent par des particules qui rentrent
  !  si on fait entrer plus de particule qu'on en fait sortir
  entreloop: do ie = 1,ntot_entree
   ip_inj = ip_inj+1
   ind_vel = int(real(nfl-1,dp)*rand_gen1(irand_i))+1
   n = nptot+ie
   vmag = sqrt(-log(one-0.99999_dp*rand_gen1(irand_i)))

   theta = two_pi*rand_gen1(irand_i)
   sin_cos = (/cos(theta),sin(theta)/)

   test_spec = rand_gen1(irand)
   particule(n)%pos(2) = s_min_loc(2) + (s_max_loc(2)-s_min_loc(2))*rand_gen1(irand_i)
   particule(n)%pos(3) = s_min_loc(3) + (s_max_loc(3)-s_min_loc(3))*rand_gen1(irand_i)

   !--Select the species which corresponds to the random generation (test_spec).
   !  The selected species is has S(is-1)%prob< test_spec <=S(is)%prob
   test_0 = zero
   do is = 1,Spe%ns 
    if(test_spec<= Spe%S(is)%prob.and.test_spec>test_0) exit 
    test_0 = Spe%S(is)%prob
   enddo
   !is = mod(is,2)+1

   particule(n)%vel(1)   = Spe%vplus(ind_vel,is)
   particule(n)%vel(2:3) = (Spe%S(is)%vth2)*(vmag*sin_cos)
   particule(n)%mass     = Spe%S(is)%sm*imul
   particule(n)%char     = Spe%S(is)%sq*imul
   particule(n)%tbirth   = t
   Spe%S(is)%ip_inj = Spe%S(is)%ip_inj + 1
  
   particule(n)%vel(2:3) = particule(n)%vel(2:3)+(/Spe%S(is)%vys,Spe%S(is)%vzs/)
 
#ifdef HAVE_WAVE_TEST
   particule(n)%vel(2:3) = particule(n)%vel(2:3)-(/by0_b,bz0_b/)
#endif
   
   !--For debug to control input speeds
   !print *,"SPEEDBO_in",particule(n)%vel(:) 

   !--Comment the next line to obtain the same result than F90 sequential
   coupdslo = rand_gen1(irand)
   particule(n)%exc  = zero
   particule(n)%orig = 0_iXb
   particule(n)%dir  = 0_iXb
   
   if(particule(n)%vel(1)>zero)then
    particule(n)%pos(1) = s_min_loc(1) + rand_gen1(irand_i)*particule(n)%vel(1)*dtp
   else
    particule(n)%pos(1) = s_max_loc(1) + rand_gen1(irand_i)*particule(n)%vel(1)*dtp
   endif
  enddo entreloop

  nlast = n
  if(nlast > npm) stop "ntpar > npm"  
  nptot = nlast       

#ifdef HAVE_DEBUG
  write(msg,'(a,2i6)')&
       &" Injected=  ",Spe%S(:)%ip_inj
  call wrt_double(qp_out,msg,0,0)
  write(msg,'(a,i18,a,i4)')&
       &" Total particles ",nptot," for process ",mpiinfo%me
  call wrtout(6,msg,"PERS")
  write(msg,'(a,i6)')" going in particles   ",ntot_entree
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
#endif

  __WRT_DEBUG_OUT("new_particles")
 end subroutine new_particles


end module particle_creation
