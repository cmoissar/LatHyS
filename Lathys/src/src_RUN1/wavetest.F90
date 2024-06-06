!!=============================================================
!!=============================================================
!!module: hybrid_model/m_wavetest
!! NAME
!!  m_wavetest (MMancini)
!!
!! FUNCTION
!!  Module containing variables and function of the WAVE TEST
module wavetest

 use defs_basis
 use defs_particletype
 use defs_parametre
 use m_writeout,only          : wrt_debug
#include "q-p_common.h"

 implicit none
 private


#ifdef HAVE_WAVE_TEST
 integer,save,public  :: nwave = 2      !--Wave harmonic
 real(dp),save,public :: awave = 0.5_dp !--Wave amplitude
 real(dp),save,public :: kwave          !--Wave number
 real(dp),save,private :: omega0        !--Pulsation
 real(dp),save,public  :: by0_b,bz0_b,by1_b,bz1_b 

 public ::                 &
      &     init_wavetest, &!--Intialisation variable wave test
      &     b_boundary_wt, &!--Boundary Cond. for Magnetic field  
      &     Amtsp3_wt       !--Velocity Intialisation

contains
 !!############################################################


 !!=============================================================
 !!routine: m_wavetest/init_wavetest
 !!
 !! FUNCTION
 !!         intialise varibles for wave test
 !! OUT
 subroutine init_wavetest()

  use defs_variable,only    : Bfield,vel,gstep,s_min,s_max,&
       &                      b0,vx0,e_conv,e1_conv,vxmean
  use defs_grid,only        : nc1
  use m_writeout

  integer :: ii
  character(len=200) :: msg

  !--Initialisation when Wave test 
  kwave = two_pi*real(nwave,dp)/(s_max(1)-s_min(1))*gstep(1)

  write(msg,'(2a,i3,3(2a,f14.7),a)')&
       &ch10,'   nwave ',nwave,&
       &ch10,'   awave ',awave,&
       &ch10,'   kwave ',kwave,ch10
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  Bfield%x(:,:,:) = b0
  
  do ii=1,nc1(1)+1
   Bfield%y(ii,:nc1(2),:nc1(3)) =  b0*awave*cos(real(ii,dp)*kwave)
   Bfield%z(ii,:nc1(2),:nc1(3)) = -b0*awave*sin(real(ii,dp)*kwave) 
  enddo

  call b_boundary_wt(Bfield,nc1(1),e_conv(2),e_conv(3),e1_conv(2),e1_conv(3),vxmean)
  
  vel%x(:nc1(1),:nc1(2),:nc1(3)) = vx0
  vel%y(:nc1(1),:nc1(2),:nc1(3)) = -Bfield%y(:nc1(1),:nc1(2),:nc1(3))
  vel%z(:nc1(1),:nc1(2),:nc1(3)) = -Bfield%z(:nc1(1),:nc1(2),:nc1(3))

 end subroutine init_wavetest


 !!=============================================================
 !!routine: m_wavetest/b_boundary_wt
 !!
 !! FUNCTION
 !!         Boundary Conditions for Magnetic field 
 !! IN
 !!  vxmean= mean velocity
 !!  nc1(3)=cell size+1
 !! OUT
 !!  ey_dm_conv,ez_dm_conv,ey_dm_conv1,ez_dm_conv1=
 !!  Electic field at the boundary
 !!
 !! SIDE EFFECT
 !!  Afield=magnetic field where boundary conditions are setted
 subroutine b_boundary_wt(Afield,ncx1,&
      &                   ey_dm_conv,ez_dm_conv,&
      &                   ey_dm_conv1,ez_dm_conv1,vxmean)

  use defs_arr3Dtype
  use defs_variable,only        : t,gstep


  type(arr3Dtype),intent(inout) :: Afield
  integer,intent(in) :: ncx1
  real(dp),intent(in) :: vxmean
  real(dp),intent(inout) :: ey_dm_conv,ez_dm_conv,ey_dm_conv1,ez_dm_conv1
  !  real(dp) :: omega0

  __WRT_DEBUG_IN("b_boundary_wt")
  !__GETTIME(63,1)


  !--Open Boundary condition in X-dirextion
  omega0 = -kwave*(vxmean+one)*t
  by0_b =  awave*cos(omega0)
  bz0_b = -awave*sin(omega0)

  omega0 = omega0+gstep(1)*kwave
  by1_b =  awave*cos(omega0)
  bz1_b = -awave*sin(omega0)

  ey_dm_conv  =  vxmean*bz0_b
  ez_dm_conv  = -vxmean*by0_b
  ey_dm_conv1 =  vxmean*bz1_b
  ez_dm_conv1 = -vxmean*by1_b

  !ax( 1,:,:) = bx0
  !ay( 1,:,:) = by0
  !az( 1,:,:) = bz0
  Afield%x(ncx1,:,:) = Afield%x(ncx1-1,:,:)
  Afield%y(ncx1,:,:) = Afield%y(ncx1-1,:,:)
  Afield%z(ncx1,:,:) = Afield%z(ncx1-1,:,:)

  !__GETTIME(63,2)
  __WRT_DEBUG_OUT("b_boundary_wt")
 end subroutine b_boundary_wt

 !!=============================================================
 !!routine: m_wavetest/Amtsp3_wt
 !!
 !! FUNCTION
 !!         
 !! IN
 !! OUT
 !! SIDE EFFECT
 subroutine Amtsp3_wt(particule,vel,n1,n2,gstep,s_min,nc1)
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

  __WRT_DEBUG_IN("Amtsp3_wt")

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
    print *,"ERROR amtsp3_wt"
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

  __WRT_DEBUG_OUT("Amtsp3_wt")
 end subroutine Amtsp3_wt
 !**************************** AMTSP3 ***********************************************

#endif
end module wavetest
