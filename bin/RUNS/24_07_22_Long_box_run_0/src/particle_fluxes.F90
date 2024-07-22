!!=============================================================
!!=============================================================
!!module: m_fluxes
!! NAME
!!  m_fluxes (RModolo,MMancini)
!!
!! FUNCTION
!!  Contains routine to compute fluxes of injected particles
!!
!! NOTE
module particle_fluxes

 use defs_basis
 use defs_species
 use m_writeout
 use m_timing,only       : time_get
#include "q-p_common.h"

 implicit none
 private
 
 !--Parameter for the speed distribution function
 integer,parameter :: kv = 15
 integer,parameter :: nv = 2**kv

 public  ::              &
      compute_fluxes      !--Compute fluxes of injected particles
contains
 !!#########################################################

 !!=============================================================
 !!routine: m_fluxes/compute_fluxes
 !!
 !! FUNCTION
 !!  Compute fluxes of injected particles
 !! IN 
 !!  is=species considered
 !!  nfl=size of array of velocity distribution
 !!
 !! SIDE EFFECT
 !!  specie(species_type)= %S(is)%vplus velocity distribution of
 !!  injected particles
 !! NOTE
 !!
 subroutine compute_fluxes(specie,is,nfl)
  !--detemination of the fluxes of particles and velocities associated

  integer,intent(in) :: nfl,is
  type(species_type),intent(inout) :: specie

  integer :: iter,ii
  integer :: vup,vdown,vtry
  real(dpo) :: fplus,fmoins,distrib
  real(dpo) :: vmin,vmax,cnorm,dv,vx_flux,sumf
  real(dpo) :: f0,dflux
  real(dp),dimension(nv) :: flux
  character(len=500) :: msg
! #ifdef HAVE_DEBUG
!   real(dp) :: check
! #endif

  __WRT_DEBUG_IN("compute_fluxes")

  vmin = specie%S(is)%vxs - 5.0d0*specie%S(is)%vth1
  vmax = specie%S(is)%vxs + 5.0d0*specie%S(is)%vth1
 
  cnorm = 1.0d0/(specie%S(is)%vth1*real(sqrt(pi),dpo))
  
  write(msg,'(a,a60,(a,a11,i12),4(a,a11,f12.6))')ch10,&
       &" _______ Computation and inversion of the fluxes of "&
       & //trim(specie%S(is)%name)//' ________', &
       & ch10,"   nv    = ",nv   ,&
       & ch10,"   vmean = ",specie%S(is)%vxs,&
       & ch10,"   vther = ",specie%S(is)%vth1,&
       & ch10,"   vmin  = ",vmin ,&
       & ch10,"   vmax  = ",vmax 
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  !--------------------------------------------
  !--Flux (or accumulation function)
  !--------------------------------------------
  !--Step of sampling for speed
  dv = (vmax-vmin)/real(nv,dpo)
  
  !--Integration of the distribution function between (0,vmax)
  sumf    = 0.0d0
  fplus   = 0.0d0
  fmoins  = 0.0d0
  vx_flux = -dv+vmin

  do ii = 1,nv
   vx_flux  = vx_flux + dv 
   distrib  = cnorm*abs(vx_flux)*exp(-((vx_flux-specie%S(is)%vxs)/specie%S(is)%vth1)**2.0d0)
   sumf     = sumf + distrib
   flux(ii) = sumf*dv

   if(vx_flux>0.0d0)then
    fplus  = fplus + distrib
   else
    fmoins = fmoins - distrib
   endif
  enddo

  fplus  = fplus*dv
  fmoins = fmoins*dv
  specie%S(is)%flux_tot = flux(nv)
  
  !--------------------------------------------
  !--Inversion of flux
  !  For each value of sampled flux: find 
  !  the associated speed by dichotomy
  !--------------------------------------------  

  !--Step of sampling of the flux
  dflux = specie%S(is)%flux_tot/real(nfl,dpo)

  loop_on_flux: do ii=1,nfl
   f0 = real(ii,dpo)*dflux !--Target

   !--Dichotomy to find index vtry such that 
   !  f0==fluxp(vtry)
   vtry  = nv/2
   vdown = 1
   vup   = nv
   iter  = 0

   dichotomy: do while ((vtry/=vup).and.(vtry/=vdown))
    iter = iter+1
    if (flux(vtry) >= f0) then
     vup = vtry
    else
     vdown = vtry
    endif
    vtry = (vup+vdown)/2
   enddo dichotomy
   specie%vplus(ii,is) = real(vtry,dpo)*dv+vmin

#ifdef HAVE_DEBUG
   ! check = (f0-fluxp(n2))/(f0+tol10)
   ! write(msg,'(3(a,e15.4),a,i2)')&
   !      &'fluxp =',f0,&
   !      &'   vx =',specie%vplus(ii,is),&
   !      &'   check=',check,'   iter=',iter
   ! call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
#endif
  enddo loop_on_flux

  write(msg,'(9(a20,g15.7,a))')&
       &"   total flux     = ",specie%S(is)%flux_tot  ,ch10,&
       &"   positive flux  = ",fplus     ,ch10,&
       &"   negative flux  = ",fmoins    ,ch10,&
       &"   vplus(nfl/2)   = ",specie%vplus(nfl/2,is),ch10
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

!!!DEBUG
!specie%vplus(:,is)=10.
!!!ENDDEBUG
! do ii =1, nfl
!  print *,'vplus__',specie%vplus(ii,is)
! enddo
! stop

  __WRT_DEBUG_OUT("compute_fluxes")
 end subroutine compute_fluxes
 !*************************** END COMPUTE_FLUXES ***********************************  


end module particle_fluxes
