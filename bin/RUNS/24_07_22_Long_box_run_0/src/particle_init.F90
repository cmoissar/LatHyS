!!=============================================================
!!=============================================================
module particle_init

 use defs_basis
 use defs_particletype
 use m_writeout

#include "q-p_common.h"

 implicit none
 private

 private::    &
      pldf1s  !--initialize all particles and any type (espece)

 public::     &
      pldf1   !--initialize particles of a type (espece)


contains
 !!#####################################################################



 !******************************** PLDF1S ***********************************************
 subroutine pldf1s(irand,particule,s_min_loc,s_max_loc,nc,ng,n1,n2, &
      &            np,npmax,spe,is)
  ! Cette routine initialise les particules
  ! d'une espece ainsi que les indices des especes
  use defs_species
  use m_rand_gen,only :     unif_dist2
  use defs_parametre,only :  ROI
  use defs_variable,only : s_max

  integer,intent(in) :: ng,n1,npmax,is
  integer,intent(inout)   :: irand
  integer,intent(in) :: nc(3)
  type(particletype),intent(inout)  :: particule(:)
  real(dp),dimension(3),intent(in)  :: s_min_loc,s_max_loc
  integer,intent(out) :: n2,np
  type(species_type),intent(inout) :: Spe
 
  integer  :: nprev,next,i,j,k,ng1,ng2
  real(dp) :: dx,dy,dz,z2,z1,y2,y1,x2,x1,xp,ng1sng2
  real(dp) :: s_r(3)
  real(dp),allocatable :: provax(:),provaz(:),provay(:)
  real(dp),allocatable :: provax2(:),provaz2(:),provay2(:)

  s_r = s_max_loc - s_min_loc

  dx = s_r(1)/real(nc(1),dp)
  dy = s_r(2)/real(nc(2),dp)
  dz = s_r(3)/real(nc(3),dp)

  ng1 = nint(real(ng,dp)*one)
  ng2 = nint(real(ng,dp)*(one+ROI%excess_part))
  ng1sng2=real(ng1,dp)/real(ng2,dp)
  nprev = n1
  allocate(provax(ng1),provay(ng1),provaz(ng1))
  allocate(provax2(ng2),provay2(ng2),provaz2(ng2))
  ng1 = nint(real(ng,dp)*one)

  do k = 1,nc(3)
   z2 = s_min_loc(3) + real(k,dp)*dz
   z1 = z2-dz
   do j = 1,nc(2)
    y2 = s_min_loc(2) + real(j,dp)*dy
    y1 = y2-dy
        if (((y1.ge.(ROI%ymin*s_max(2))).and.&
        & (y1.le.(ROI%ymax*s_max(2)))).and.&
        & ((z1.ge.(ROI%zmin*s_max(3))).and.(z1.le.(ROI%zmax*s_max(3))))) then

    do i = 1,nc(1)
     x2 = s_min_loc(1) + real(i,dp)*dx
     x1 = x2 - dx
     xp = half*(x2+x1)-s_min_loc(1)   
     next = nprev-1+ng2
     if (next > npmax) stop ' PLDF1S : Too many particles'
     call unif_dist2(irand,x1,x2,provax2(:))
     call unif_dist2(irand,y1,y2,provay2(:))
     call unif_dist2(irand,z1,z2,provaz2(:))
     particule(nprev:next)%pos(1) = provax2
     particule(nprev:next)%pos(2) = provay2
     particule(nprev:next)%pos(3) = provaz2
     particule(nprev:next)%mass = Spe%S(is)%sm*ng1sng2
     particule(nprev:next)%char = Spe%S(is)%sq*ng1sng2
     nprev = next+1
    enddo

        else

    do i = 1,nc(1)
     x2 = s_min_loc(1) + real(i,dp)*dx
     x1 = x2 - dx
     xp = half*(x2+x1)-s_min_loc(1)   
     next = nprev-1+ng1
     if (next > npmax) stop ' PLDF1S : Too many particles'
     call unif_dist2(irand,x1,x2,provax(:))
     call unif_dist2(irand,y1,y2,provay(:))
     call unif_dist2(irand,z1,z2,provaz(:))
     particule(nprev:next)%pos(1) = provax
     particule(nprev:next)%pos(2) = provay
     particule(nprev:next)%pos(3) = provaz
     particule(nprev:next)%mass = Spe%S(is)%sm
     particule(nprev:next)%char = Spe%S(is)%sq
     nprev = next+1
    enddo

        endif

   enddo
  enddo

  deallocate(provax,provay,provaz)
  deallocate(provax2,provay2,provaz2)
  n2 = next
  np = n2-n1+1

 end subroutine pldf1s
 !**************************** END PLDF1S *****************************************


 !*********************************** PLDF1 ***************************************************
 subroutine pldf1(irand,particule,&
      &           Spe,&
      &           s_min_loc,s_max_loc,&
      &           nc,&
      &           nptot,npm)
  !=====================
  !initialise la position des particules et l'indice des espces
  
  use defs_species

  integer,intent(in) :: npm
  integer,intent(inout) :: irand
  integer,intent (out) :: nptot
  integer,intent(in) :: nc(3)
  real(dp),dimension(3),intent(in)   :: s_min_loc,s_max_loc
  type(species_type),intent(inout) :: Spe
  type(particletype),intent(inout) :: particule(:)

  integer :: npmax,n0,is,ns
  real(dp)  :: psq_min,psq_max,psm_min,psm_max
  character(len=500) :: msg

  __WRT_DEBUG_IN("pldf1")


  !--Dimension de sortie (npm)
  npmax = npm
  n0    = 0
  nptot = 0
  ns = Spe%ns

  !--Loading Spe
  do is = 1,ns
   Spe%S(is)%n1 = n0+1
   call pldf1s(irand,particule,s_min_loc,s_max_loc,&
        &      nc,&
        &      Spe%S(is)%ng,&
        &      Spe%S(is)%n1,&
        &      Spe%S(is)%n2,&
        &      Spe%S(is)%np,&
        &      npmax,spe,is)

   nptot = nptot + Spe%S(is)%np
   particule(Spe%S(is)%n1:Spe%S(is)%n2)%exc  = zero
   particule(Spe%S(is)%n1:Spe%S(is)%n2)%orig = 0_iXb
   particule(Spe%S(is)%n1:Spe%S(is)%n2)%dir  = 0_iXb
   particule(Spe%S(is)%n1:Spe%S(is)%n2)%tbirth = 0._dp

   if (nptot > npm) stop 'PLDF1 : Too many particles'

   write(msg,'(3x,a3,a,i10)') trim(Spe%S(is)%name),&
        &' macro-particles in the box: ',Spe%S(is)%n2-Spe%S(is)%n1+1
   call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
   
   n0 = Spe%S(is)%n2

  enddo

  psm_min = minval(particule(Spe%S(1)%n1:Spe%S(ns)%n2)%mass)
  psm_max = maxval(particule(Spe%S(1)%n1:Spe%S(ns)%n2)%mass)
  psq_min = minval(particule(Spe%S(1)%n1:Spe%S(ns)%n2)%char)
  psq_max = maxval(particule(Spe%S(1)%n1:Spe%S(ns)%n2)%char)

  write(msg,'(2(2(a12,e12.4),a))')&
       &'psq_min  = ',psq_min,'psq_max = ',psq_max,ch10,&
       &'psm_min  = ',psm_min,'psm_max = ',psm_max,ch10
  call wrt_double(6,msg,wrtscreen,wrtdisk)

  __WRT_DEBUG_OUT("pldf1")
 end subroutine pldf1
 !******************************** END PLDF1 **************************************************

end module particle_init
