!!=============================================================
!!=============================================================
!!module: atm_ionosphere
!! NAME
!!  atm_ionosphere (SHess)
!!
!!
!! FUNCTION
!!  Contains generic subroutines for the ionosphere computation
!!
!! NOTE
module atm_ionosphere

 use defs_basis
 use defs_species
 use defs_atmospheretype
 use defs_particletype
 use defs_parametre
 use m_writeout
 use m_rand_gen,only                   : rand_gen1 

#include "q-p_common.h"

 implicit none
 private
 
 public ::                 &
      create_ionosphere_generic,  & 
      split_particle_iono,&
      iono_densities_generic,&
      Ion_production_generic,&
      dummy
private ::&
      set_particle_ionosphere,&
      add_particle_ionosphere,&
      compute_ei

contains

! Dummy procedure, used if no fee_ionopshere procedure deffined
subroutine dummy(i,j,k,atmosphere,l,t3,r2,rb,Spe,s_cen)
  type(species_type),intent(in) :: Spe
  type(atmosphere_type),intent(in) ::atmosphere
  real(dp),intent(in) ::r2,rb
  real(dp),intent(in),dimension(3)::s_cen
  real(dp),intent(inout) ::t3
  integer,intent(in) :: i,j,k,l
end subroutine dummy


! routine adding ionosphere (basically only loops on spatial variables and the particle species)
subroutine create_ionosphere_generic(Spe,s_cen,s_min_loc,particule,gstep,s_min_loc_i,s_max_loc_i,irand,nptot,npcell,atmosphere)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  real(dp),dimension(3),intent(in) ::s_cen
  integer,intent(inout) :: nptot,irand
  integer,intent(in) :: npcell,s_min_loc_i(3),s_max_loc_i(3)
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(in) ::atmosphere
  real(dp) ::r_lim2,rp2,rb,radius
  integer::i,j,k,lp_reac
  character(len=100) :: msg

__WRT_DEBUG_IN("create_ionosphere_generic")
  r_lim2=(Spe%P%r_iono+sqrt(gstep(1)**2+gstep(2)**2+gstep(3)**2))**2
  rp2=(Spe%P%r_lim-sqrt(gstep(1)**2+gstep(2)**2+gstep(3)**2))**2
  do k=s_min_loc_i(3),s_max_loc_i(3)
          do j=s_min_loc_i(2),s_max_loc_i(2)
              rb = (float(j-1)*gstep(2)-s_cen(2))**2+(float(k-1)*gstep(3)-s_cen(3))**2
                  do i=s_min_loc_i(1),s_max_loc_i(1)
                radius = (float(i-1)*gstep(1)-s_cen(1))**2+rb
                if (((radius < rp2).or.(radius > r_lim2))) CYCLE 
                        do lp_reac=1,atmosphere%n_species
                if(atmosphere%species(lp_reac)%iono) &
                call set_particle_ionosphere(i,j,k,gstep,s_min_loc,&
                &  atmosphere%species(lp_reac)%density,atmosphere%species(lp_reac)%mass,irand,&
                &  nptot,npcell,particule,Spe)
                enddo
                  enddo
          enddo
  enddo
__WRT_DEBUG_OUT("create_ionosphere_generic")
end subroutine create_ionosphere_generic

! routine effectively adding ionosphere particles
subroutine set_particle_ionosphere(i,j,k,gstep,s_min_loc,density,mass,irand,nptot,npcell,particle,Spe)

  type(species_type),intent(in) :: Spe
  integer,intent(in)::i,j,k,npcell
  real(dp),intent(in) :: gstep(3),mass,s_min_loc(3)
  integer,intent(inout) :: nptot,irand
  real(dp),dimension(:,:,:),intent(in) ::density
  type(particletype),intent(inout) :: particle(:)
 
  real(dp) ::s(3),w(8),ratio,d_cell,pos(3),rp2,rl2
  integer::n,n_cell,sz
  character(len=100) :: msg

  rp2=Spe%P%r_lim**2;  rl2=Spe%P%r_iono**2; sz=size(particle)-100

 d_cell=0.125_dp*(density(i,j,k)    + density(i+1,j+1,k+1)+&
                 density(i+1,j,k)  + density(i+1,j+1,k)+&
                 density(i+1,j,k+1)+ density(i,j+1,k)+&
                 density(i,j+1,k+1)+ density(i,j,k+1))

 n_cell=min(npcell,max(int(d_cell),1))
 d_cell=1._dp/float(n_cell)

 do n=1,n_cell
        if(nptot.ge.sz) CYCLE
        s(1)=rand_gen1(irand) ; s(2)=rand_gen1(irand) ; s(3)=rand_gen1(irand)
        w(1)=(1._dp-s(1))*(1._dp-s(2))*(1._dp-s(3)); w(2)=s(1)*s(2)*s(3)
        w(3)=s(1)*(1._dp-s(2))*(1._dp-s(3)); w(4)=s(1)*s(2)*(1._dp-s(3))
        w(5)=s(1)*(1._dp-s(2))*s(3); w(6)=(1._dp-s(1))*s(2)*(1._dp-s(3))
        w(7)=(1._dp-s(1))*s(2)*s(3); w(8)=(1._dp-s(1))*(1._dp-s(2))*s(3)
        ratio   =d_cell*(w(1)*density(i,j,k)    + w(2)*density(i+1,j+1,k+1)+&
                 w(3)*density(i+1,j,k)  + w(4)*density(i+1,j+1,k)+&
                 w(5)*density(i+1,j,k+1)+ w(6)*density(i,j+1,k)+&
                 w(7)*density(i,j+1,k+1)+ w(8)*density(i,j,k+1))
        if((ratio.ge.0.001)) then 
        pos(1)=(float(i-1)+s(1))*gstep(1); pos(2)=(float(j-1)+s(2))*gstep(2); pos(3)=(float(k-1)+s(3))*gstep(3)
        pos=pos+s_min_loc
        if(((sum((pos-Spe%P%centr)**2)).le.rp2).or.((sum((pos-Spe%P%centr)**2)).ge.rl2)) CYCLE
        
        nptot=nptot+1
        particle(nptot)%pos(1:3) = pos(1:3)
        particle(nptot)%vel(1:3) = zero
        particle(nptot)%mass   = mass*ratio
        particle(nptot)%char   = ratio
        particle(nptot)%exc    = 0
        particle(nptot)%orig   = 3
        endif
 enddo
end subroutine set_particle_ionosphere


! routine splitting heavy macro-particles above a given distance from the planet
subroutine split_particle_iono(Spe,particle,n,nptot,gstep,dist2,irand, &
        &       s_min_loc,s_max_loc)
  use m_timing,only           : time_get

  type(particletype),intent(inout) :: particle(:)
  type(species_type),intent(in) :: Spe
  real(dp),dimension(3),intent(in) :: s_min_loc,s_max_loc
  integer,intent(inout) :: nptot,irand
  integer,intent(in) ::n
  real(dp),intent(in) :: gstep(3),dist2
  integer::sz
  real(dp) ::s(3)
  real(dp),parameter :: separation=0.01,rlim=700.,qlim=3.,mlim=6.

  __GETTIME(77,1)!--Timer start
 
  sz=size(particle)

  if((particle(n)%char.gt.qlim).AND.(nptot <int(0.8*npm))) then 
  if(particle(n)%mass.gt.mlim) then
  if((nptot+1).lt.sz) then
  if(dist2.gt.((rlim/Spe%ref%c_omegapi+Spe%P%r_lim)**2)) then
!__WRT_DEBUG_IN("split_particle_iono")
! check thta the splitting will not introduce a new particle outside
! of the subdomain
  if ((particle(n)%pos(2) >s_min_loc(2)+2*separation*gstep(2)).AND.&
      & (particle(n)%pos(2) < s_max_loc(2)-2*separation*gstep(2)).AND. &
      & (particle(n)%pos(3) >s_min_loc(3)+2*separation*gstep(3)).AND.&
      & (particle(n)%pos(3) < s_max_loc(3)-2*separation*gstep(3))) then
      
          nptot=nptot+1
          s(1)=rand_gen1(irand) ; s(2)=rand_gen1(irand) ; s(3)=rand_gen1(irand)
          s=s/sqrt(sum(s**2))*separation*0.5*gstep
          particle(n)%mass       = particle(n)%mass*0.5
          particle(n)%char       = particle(n)%char*0.5
          particle(nptot)        = particle(n)
          particle(nptot)%pos(1:3) = particle(n)%pos(1:3)-s
          particle(n)%pos(1:3)     = particle(n)%pos(1:3)+s
   endif
!__WRT_DEBUG_OUT("split_particle_iono")
   endif
   endif
   endif
   endif
  __GETTIME(77,2)!--Timer stops
 
end subroutine split_particle_iono

subroutine add_particle_ionosphere(s,gstep,s_min_loc,mass,irand,nptot,particle,Spe,ratio,origine)
  type(species_type),intent(in) :: Spe
  integer,intent(in)::origine
  real(dp),intent(in) :: gstep(3),mass,s_min_loc(3),ratio,s(3)
  integer,intent(inout) :: nptot,irand
  type(particletype),intent(inout) :: particle(:)
  real(dp)::pos(3)
        pos=(s)*gstep+s_min_loc
        nptot=nptot+1
        particle(nptot)%pos(1:3) = pos(1:3)
        particle(nptot)%vel(1:3) = zero
        particle(nptot)%mass   = mass*ratio
        particle(nptot)%char   = ratio
        particle(nptot)%exc    = 0
        particle(nptot)%orig   = origine
end subroutine add_particle_ionosphere

subroutine iono_densities_generic(particle,nptot,atmosphere,Spe,gstep,s_min_loc)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  integer,intent(inout) :: nptot
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particle(:)
  type(atmosphere_type),intent(inout) ::atmosphere
  integer ::i,j,k,ijk(3)
  real(dp) :: qsm,qsm_ind(atmosphere%n_species,2),s(3),pos(3),w(8),nrm

__WRT_DEBUG_IN("iono_densities_generic")  
j=0
do i=1,atmosphere%n_species
if (atmosphere%species(i)%iono) then 
        qsm=0.01*real(int(atmosphere%species(i)%charge/atmosphere%species(i)%mass*100._dp+0.5_dp),dp)
        j=j+1;  qsm_ind(j,1)=qsm;       qsm_ind(j,2)=real(i,dp)
        atmosphere%species(i)%density(:,:,:)=0._dp
endif
enddo
do i=1,nptot
!if ((particle(i)%pos(1)>s_min_loc(1)).AND.(particle(i)%pos(2)>s_min_loc(2)).AND. &
!        & (particle(i)%pos(3)>s_min_loc(3))) then

qsm=0.01*real(int(particle(i)%char/particle(i)%mass*100._dp+0.5_dp),dp)
        do k=1,j
                if (qsm.eq.qsm_ind(k,1)) then
                pos=(particle(i)%pos-s_min_loc)/gstep+1._dp
                ijk=int(pos);   s=pos-real(ijk,dp)
                if ((ijk(1) < 1).OR.(ijk(2) < 1).OR. (ijk(3)<1))  print *, 'pos',particle(i)%pos,s_min_loc
                nrm=particle(i)%char/atmosphere%species(int(qsm_ind(k,2)))%charge*Spe%ref%density/1E6
                w(1)=(1._dp-s(1))*(1._dp-s(2))*(1._dp-s(3)); w(2)=s(1)*s(2)*s(3)
                w(3)=s(1)*(1._dp-s(2))*(1._dp-s(3)); w(4)=s(1)*s(2)*(1._dp-s(3))
                w(5)=s(1)*(1._dp-s(2))*s(3); w(6)=(1._dp-s(1))*s(2)*(1._dp-s(3))
                w(7)=(1._dp-s(1))*s(2)*s(3); w(8)=(1._dp-s(1))*(1._dp-s(2))*s(3)
                w=w*nrm
                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1),ijk(2),ijk(3))       =&
&                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1),ijk(2),ijk(3))+w(1)
                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1)+1,ijk(2)+1,ijk(3)+1) =&
&                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1)+1,ijk(2)+1,ijk(3)+1)+w(2)
                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1)+1,ijk(2),ijk(3))     =&
&                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1)+1,ijk(2),ijk(3))+w(3)
                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1)+1,ijk(2)+1,ijk(3))   =&
&                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1)+1,ijk(2)+1,ijk(3))+w(4)
                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1)+1,ijk(2),ijk(3)+1)   =&
&                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1)+1,ijk(2),ijk(3)+1)+w(5)
                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1),ijk(2)+1,ijk(3))     =&
&                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1),ijk(2)+1,ijk(3))+w(6)
                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1),ijk(2)+1,ijk(3)+1)   =&
&                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1),ijk(2)+1,ijk(3)+1)+w(7)
                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1),ijk(2),ijk(3)+1)     =&
&                atmosphere%species(int(qsm_ind(k,2)))%density(ijk(1),ijk(2),ijk(3)+1)+w(8)
                endif
        enddo
!   endif

enddo
__WRT_DEBUG_OUT("iono_densities_generic")  
end subroutine iono_densities_generic

subroutine Ion_production_generic(Spe,atmosphere,particule,w_max,w_min,nptot,irand,gstep,s_min_loc,iono_func)
 use defs_grid
 use defs_parametre,only      : dt

  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(inout) ::atmosphere
  real(dp),intent(in) ::w_max,w_min
  integer,intent(inout) :: nptot,irand
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  procedure(dummy),optional ::iono_func

  real(dp)::weight,s(3),t1,t2,t3,r2,radius,rb,r_iono2,tirage,s_cen(3),w_min_adapt,prod_unit
  integer ::i,j,k,l,sz

 call iono_densities_generic(particule,nptot,atmosphere,Spe,gstep,s_min_loc)
 radius=Spe%P%r_lim**2;  r_iono2=Spe%P%r_iono**2;  sz=size(particule)-100
 s_cen = Spe%P%centr-s_min_loc 
 prod_unit=1E6/Spe%ref%density*dt*Spe%ref%inv_gyro
         do l=1,atmosphere%n_species
                        if(atmosphere%species(l)%iono) then
 do k=1,ncm(3)-1
        do j=1,ncm(2)-1
                do i=1,ncm(1)-1
                                if(nptot.ge.sz) CYCLE
                                        s(1)=i+rand_gen1(irand)-1.5_dp 
                                        s(2)=j+rand_gen1(irand)-1.5_dp 
                                        s(3)=k+rand_gen1(irand)-1.5_dp
                                        if (any(s.lt.(0._dp)).or.any((s-ncm).ge.0)) CYCLE
                                        rb = (s(2)*gstep(2)-s_cen(2))**2+(s(3)*gstep(3)-s_cen(3))**2
                                        r2 = (s(1)*gstep(1)-s_cen(1))**2+rb
                                        if (r2.lt.radius) CYCLE
                                        t1=0.; t2=0.; t3=0.
                                        if (r2.lt.r_iono2) then
                                                if (present(iono_func)) call iono_func(i,j,k,atmosphere,l,t3,r2,rb,Spe,s_cen)
                                        else
                                                if (associated(atmosphere%species(l)%prod)) t1=atmosphere%species(l)%prod(i,j,k)
                                                call compute_ei(i,j,k,atmosphere,l,t2,r2,rb,Spe,s_cen)
                                                endif
                                        t1=t1*prod_unit; t2=t2*prod_unit;t3=t3*prod_unit
                                        w_min_adapt=5.*(t1+t2+t3)
                                        w_min_adapt=min(max(w_min,w_min_adapt),w_max)
                                        tirage=rand_gen1(irand)*w_min_adapt
                                        weight=max(min(t1+t2+t3,w_max),w_min_adapt)
                                        if (tirage.lt.t1) then 
                                                call add_particle_ionosphere(s,gstep,s_min_loc,&
                                                atmosphere%species(l)%mass,irand,nptot,particule,Spe,weight,1)
                                        else 
                                           if (tirage.lt.(t1+t2)) then 
                                                call add_particle_ionosphere(s,gstep,s_min_loc,&
                                           atmosphere%species(l)%mass,irand,nptot,particule,Spe,weight,2)
                                           else
                                                if (tirage.lt.(t1+t2+t3)) then 
                                                call add_particle_ionosphere(s,gstep,s_min_loc,&
                                                atmosphere%species(l)%mass,irand,nptot,particule,Spe,weight,3)
                                                 endif!test ionosphere
                                           endif!test electron impact
                                        endif!test photoproduction
                enddo!i
        enddo!j
   enddo!k
                                endif!test ionospheric species
                         enddo!loop on species
end subroutine Ion_production_generic


subroutine compute_ei(i,j,k,atmosphere,l,t2,r2,rb,Spe,s_cen)
use defs_variable, only : dn_e_incdt,te
use defs_parametre, only : ipe

  integer,intent(in)::i,j,k,l
  real(dp),intent(inout)::t2
  real(dp),intent(in)::s_cen(3),r2,rb
  type(species_type),intent(in) :: Spe
  type(atmosphere_type),intent(inout) ::atmosphere

   real(dp) :: t,arg,gamma1,rate,conv_t
   integer  :: ii,jj,n

t2=0._dp
!conv_t=2.88E28 !1./(2mu_0k_b)
!conv_t=conv_t*Spe%betae*(Spe%ref%mag**2)/Spe%ref%density
!if (mod(ipe,2) == 0) then   
!        gammam1=2._dp/3._dp  
!else    
!        gammam1=0._dp
!endif
!t=log(max((te*Spe%ref%conv_t*dn_e_incdt(i,j,k)**gammam1),1._dp))
conv_t = 2.*sum(Spe%S(:)%qms*Spe%S(:)%rmds)*(Spe%ref%mag**2)/(Spe%ref%density*2._dp*mu0*e_Cb)

gamma1 = 5._dp/3._dp
t = Spe%betae/e_Cb * dn_e_incdt(i,j,k)**(gamma1-1) * Spe%ref%mag**2 / (2*mu0*Spe%ref%density)
t = log(t*e_Cb/kb_JK)

do ii=1,atmosphere%n_ei
 if ((atmosphere%species(l)%name).eq.atmosphere%ei_reactions(ii)%ion%name) then
    n=atmosphere%ei_reactions(ii)%n-1
    arg=0._dp
        do jj=0,n 
                arg=arg+atmosphere%ei_reactions(ii)%coeff(jj+1)*(t**jj)
        enddo
    if (arg.gt.(-10._dp)) then
    rate=exp(arg)
    t2 = t2+rate*atmosphere%ei_reactions(ii)%neutral%density(i,j,k)*dn_e_incdt(i,j,k)*Spe%ref%density!divided by ref%density later
    endif
 endif
enddo
t2=t2*1.e-6

end subroutine compute_ei

end module atm_ionosphere
