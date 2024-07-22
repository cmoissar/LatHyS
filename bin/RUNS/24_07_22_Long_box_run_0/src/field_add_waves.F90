!!=============================================================
!!=============================================================
!!module: hybrid_model/m_wave
!! NAME
!! field_add_waves
!!
!! FUNCTION
!!  Module containing waves continually added to the solar wind
!!  at the entrance of the box
module field_add_waves

 use defs_basis
 use defs_particletype
 use defs_parametre
 use m_writeout,only          : wrt_debug
#include "q-p_common.h"

 implicit none
 private

 integer,save,public  :: nwave = 50     !--Wave harmonic
 real(dp),save,public :: awave = 0.0_dp !--Wave amplitude
 real(dp),save,public :: kwave          !--Wave number
 real(dp),save,private :: omega         !--Pulsation
 
 public ::                 &
      &     init_waves, &!--Intialisation waves
      &     add_waves,  &!--continual addition of waves  
      &     Amtsp3_waves !--Initialisation velocities

contains
 
 subroutine init_waves()

  use defs_variable,only    : Bfield,vel,s_min,s_max,&
       &                      b0,vx0,e_conv,e1_conv,vxmean,&
                              bx0,by0,bz0,by1,bz1
  use defs_parametre, only  : gstep
  use defs_grid,only        : nc1
  use m_writeout

  integer :: ii
  real :: by_sw, bz_sw

  kwave = two_pi*real(nwave,dp)/(s_max(1)-s_min(1))*gstep(1)

  Bfield%x(:,:,:) = bx0
  
  do ii=1,nc1(1)+1
   Bfield%y(ii,:nc1(2),:nc1(3)) =  by0! + b0*awave*cos(real(ii,dp)*kwave)
   Bfield%z(ii,:nc1(2),:nc1(3)) =  bz0 - b0*awave*sin(real(ii,dp)*kwave) 
  enddo

  !--Open Boundary condition in X-dirextion
  by_sw = by0
  by0 = by0! + awave*cos(0.0_dp)
  bz_sw = bz0
  bz0 = bz0 - awave*sin(0.0_dp)

  by1 = by0! + awave*cos(gstep(1)*kwave)
  bz1 = bz0 - awave*sin(gstep(1)*kwave)

  e_conv(2)  =  vxmean*bz0
  e_conv(3)  = -vxmean*by0
  e1_conv(2) =  vxmean*bz1
  e1_conv(3) = -vxmean*by1

  Bfield%x(nc1(1),:,:) = Bfield%x(nc1(1)-1,:,:)
  Bfield%y(nc1(1),:,:) = Bfield%y(nc1(1)-1,:,:)
  Bfield%z(nc1(1),:,:) = Bfield%z(nc1(1)-1,:,:)
  
  !--Initiatialisation of the velocity field
  vel%x(:nc1(1),:nc1(2),:nc1(3)) = vx0
  vel%y(:nc1(1),:nc1(2),:nc1(3)) = -(Bfield%y(:nc1(1),:nc1(2),:nc1(3)) - by_sw)
  vel%z(:nc1(1),:nc1(2),:nc1(3)) = -(Bfield%z(:nc1(1),:nc1(2),:nc1(3)) - bz_sw)

 end subroutine init_waves


 !!=============================================================
 !!routine: field_add_waves/Amtsp3_waves

 !! FUNCTION initiate particle velocities in the box,
 !!          using the velocities defined just above         
 !!
 !! IN
 !! OUT
 !! SIDE EFFECT
 subroutine Amtsp3_waves(particule,vel,n1,n2,gstep,s_min,nc1)
  !--Calcul de la densite de charge (rho) au temps 0    -> dn
  !--et des composantes du courant (Ji) au temps 0      -> vel%x,vel%y,vel%z

  use defs_arr3Dtype
  use defs_mpitype,only     : mpitype

  integer,intent(in) :: n1,n2
  integer,intent(in) :: nc1(3)
  real(dp),intent(in) :: gstep(3),s_min(3)
  type(particletype),intent(inout) :: particule(:)
  type(arr3Dtype),intent(in) ::  vel

  integer :: n
  integer :: ijk(3)
  real(dp) :: w(nb_voisins)
  real(dp),dimension(3) :: s_f,s_a,s_m,gstep_inv

  __WRT_DEBUG_IN("Amtsp3_waves")

  gstep_inv = one/gstep

  do n = n1,n2
   !--On collecte la composante X,Y,Z de la vitesse

   !--Relative position in cell s_f=(xf,yf,zf) center of the particule
   s_m = one + (particule(n)%pos-s_min)*gstep_inv

   !--Sequence of indices of B at cell corners
   !   (i,j,k)
   ijk = int(s_m,dp)

   s_f = s_m-real(ijk,dp)

#ifdef HAVE_DEBUG   
   if(any(ijk>nc1)) then
    print *,"ERROR amtsp3_waves"
    print *,"n ",n
    print *,"part ",particule(n)%pos
    print *,"s_min ",s_min
    print *,"ijk",ijk
    stop
   endif
#endif

   !--Trilinear weight
   s_a = one-s_f

   w(1) = s_a(1)*s_a(2)*s_a(3)
   w(2) = s_f(1)*s_a(2)*s_a(3)
   w(3) = s_a(1)*s_f(2)*s_a(3)
   w(4) = s_f(1)*s_f(2)*s_a(3)
   w(5) = s_a(1)*s_a(2)*s_f(3)
   w(6) = s_f(1)*s_a(2)*s_f(3)
   w(7) = s_a(1)*s_f(2)*s_f(3)
   w(8) = s_f(1)*s_f(2)*s_f(3)

   particule(n)%vel(2) = particule(n)%vel(2) +&
        &  vel%y(ijk(1)  ,ijk(2)  ,ijk(3)  ) * w(1) + & 
        &  vel%y(ijk(1)+1,ijk(2)  ,ijk(3)  ) * w(2) + &
        &  vel%y(ijk(1)  ,ijk(2)+1,ijk(3)  ) * w(3) + &
        &  vel%y(ijk(1)+1,ijk(2)+1,ijk(3)  ) * w(4) + &
        &  vel%y(ijk(1)  ,ijk(2)  ,ijk(3)+1) * w(5) + &
        &  vel%y(ijk(1)+1,ijk(2)  ,ijk(3)+1) * w(6) + &
        &  vel%y(ijk(1)  ,ijk(2)+1,ijk(3)+1) * w(7) + &
        &  vel%y(ijk(1)+1,ijk(2)+1,ijk(3)+1) * w(8)

   particule(n)%vel(3) = particule(n)%vel(3) +&
        &  vel%z(ijk(1)  ,ijk(2)  ,ijk(3)  ) * w(1) + & 
        &  vel%z(ijk(1)+1,ijk(2)  ,ijk(3)  ) * w(2) + &
        &  vel%z(ijk(1)  ,ijk(2)+1,ijk(3)  ) * w(3) + &
        &  vel%z(ijk(1)+1,ijk(2)+1,ijk(3)  ) * w(4) + &
        &  vel%z(ijk(1)  ,ijk(2)  ,ijk(3)+1) * w(5) + &
        &  vel%z(ijk(1)+1,ijk(2)  ,ijk(3)+1) * w(6) + &
        &  vel%z(ijk(1)  ,ijk(2)+1,ijk(3)+1) * w(7) + &
        &  vel%z(ijk(1)+1,ijk(2)+1,ijk(3)+1) * w(8)
  enddo

  __WRT_DEBUG_OUT("Amtsp3_waves")
 end subroutine Amtsp3_waves
 !**************************** AMTSP3 ***********************************************


 subroutine add_waves(i, By_alfven, Bz_alfven, By1_alfven, Bz1_alfven, Vy_alfven, Vz_alfven)

 use defs_variable , only : b0, vxmean
 use defs_parametre, only : gstep, dt

  __WRT_DEBUG_IN("add_waves")
  
  integer,  intent(in)  :: i
  real(dp), intent(out) :: By_alfven, Bz_alfven, By1_alfven, Bz1_alfven
  real(dp), intent(out) :: Vy_alfven, Vz_alfven
 
  omega = -kwave*(vxmean+one)

  !--Open Boundary condition in X-dirextion

  !!  Alfven Waves
  By_alfven = zero !b0*awave*cos(omega*i*dt)
  Bz_alfven = -b0*awave*sin(omega*i*dt)

  ! need 2 grid points to initiate maxwell calculations
  By1_alfven = zero !b0*awave*cos(omega*i*dt+gstep(1)*kwave)
  Bz1_alfven = -b0*awave*sin(omega*i*dt+gstep(1)*kwave)

  Vy_alfven = -By_alfven
  Vz_alfven = -Bz_alfven
 
  __WRT_DEBUG_OUT("add_waves")
 
 end subroutine add_waves

end module field_add_waves
