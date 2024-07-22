!!=============================================================
!!=============================================================
!!module: atm_charge_exchange
!! NAME
!!  atm_charge_exchange (SHess)
!!
!! Contains the generic routine in charge of the particle charge exchange
!!
!!
!! NOTE
module atm_charge_exchange

 use defs_basis
 use defs_species
 use defs_atmospheretype
 use defs_particletype
 use defs_atmospheretype
 use defs_parametre
 use m_writeout
 use m_timing
 use m_rand_gen,only            : rand_gen1

#include "q-p_common.h"

contains
!!=============================================================
 !!routine: env_mars/charge_exchange_mars
 !!
 !! FUNCTION
 !!  compute charge-exchange
 !! IN 
 !! OUT
 !!   
 !! SIDE EFFECT
 !!  
 subroutine charge_exchange_generic(nn,pickup,qsm,irand,ijk,v_p,w,Spe,particule,atmosphere)
   
  integer,intent(in) :: nn,ijk(3)
  integer,intent(inout) :: irand,pickup
  real(dp),intent(in) :: qsm,v_p(3),w(8)
  type(species_type),intent(in) :: Spe
  type(particletype),intent(inout) :: particule(:)
  type(atmosphere_type),intent(in) ::atmosphere

  integer :: lp_reac
  real(dp) :: density,vmod,proba_coll,test_val,conv_mass,opa

  __GETTIME(74,1)!--Timer start

  vmod=sqrt(sum(v_p*v_p))
 do lp_reac=1,atmosphere%n_exc
        if(qsm /= atmosphere%exc_reactions(lp_reac)%qsm) CYCLE !particle specie identified by it q/m ratio
   opa = atmosphere%exc_reactions(lp_reac)%cross_section*vmod
   ! interpolation takes time, we estimate first if it is necessary.
   proba_coll = opa*atmosphere%exc_reactions(lp_reac)%neutral%density(ijk(1),ijk(2),ijk(3))
  if(proba_coll>1e-6) then 
   proba_coll =   proba_coll*1E4
   test_val = rand_gen1(irand)
  if(test_val <= proba_coll) then

  !--Neutral density
  density = w(1)*atmosphere%exc_reactions(lp_reac)%neutral%density(ijk(1)  ,ijk(2)  ,ijk(3)  )&
         & +w(2)*atmosphere%exc_reactions(lp_reac)%neutral%density(ijk(1)+1,ijk(2)  ,ijk(3)  )&
         & +w(3)*atmosphere%exc_reactions(lp_reac)%neutral%density(ijk(1)  ,ijk(2)+1,ijk(3)  )&
         & +w(4)*atmosphere%exc_reactions(lp_reac)%neutral%density(ijk(1)+1,ijk(2)+1,ijk(3)  )&
         & +w(5)*atmosphere%exc_reactions(lp_reac)%neutral%density(ijk(1)  ,ijk(2)  ,ijk(3)+1)&
         & +w(6)*atmosphere%exc_reactions(lp_reac)%neutral%density(ijk(1)+1,ijk(2)  ,ijk(3)+1)&
         & +w(7)*atmosphere%exc_reactions(lp_reac)%neutral%density(ijk(1)  ,ijk(2)+1,ijk(3)+1)&
         & +w(8)*atmosphere%exc_reactions(lp_reac)%neutral%density(ijk(1)+1,ijk(2)+1,ijk(3)+1) 

   proba_coll   = 1._dp-exp(-opa*density)
   if(test_val <= proba_coll) then
  !--Charge exchange with a planetary neutral
    conv_mass=atmosphere%exc_reactions(lp_reac)%neutral%mass/atmosphere%exc_reactions(lp_reac)%ion%mass
    particule(nn)%vel(1:3)   = (/Spe%P%speed,0._dp,0._dp/)
    particule(nn)%mass = conv_mass*particule(nn)%mass
    particule(nn)%char = particule(nn)%char
    particule(nn)%exc  = particule(nn)%exc+1._dp
    pickup = pickup+1   ! count the number of pickups
    EXIT  
    endif
  endif
  endif
 enddo
  __GETTIME(74,2)!--Timer stop
 end subroutine charge_exchange_generic

end module atm_charge_exchange

