!!=============================================================
!!=============================================================
!!file rand_gen.f90
!!
!! Contains modules for generation of random number (RM,MM)
!!
module m_rand_gen

 use defs_basis

#include "q-p_common.h"

 implicit none
 private

 !--Random Variables :
 integer,private :: ix11,ix21,ix31,j1,iff1
 !--Tableaux pour la gestion des nombres aléatoires
 !  utilises lors de la fonction "ran1" et les restart 
 real(dp),private :: R1(97) 

 integer,public,save :: irand = -5  !--Valeur de gsnsration de nombre aleatoire

 
 public::            &
      rand_vars_put, &  !--Assign random variables
      rand_vars_get, &  !--Extract random variables
      rand_gen1,     &  !--generate a random number (scalar)
      unif_dist,     &  !--uniform distribution in a range (array)
      unif_dist2,    &  !--uniform distribution in a range (array)
      bi_max_dib        !--Bi-maxwellian distribution (array)
contains
 !!############################################################

 !=============================================================
 !!subroutine: rand_gen/rand_vars_get
 subroutine rand_vars_get(rstate,istate)
  real(dp),intent(out) :: rstate(97)
  integer,intent(out) :: istate(4)
  istate(1) = ix11
  istate(2) = ix21
  istate(3) = ix31
  istate(4) = j1
  rstate    = R1 
 end subroutine rand_vars_get

 !=============================================================
 !!subroutine: rand_gen/rand_vars_put
 subroutine rand_vars_put(rstate,istate,iff)
  real(dp),intent(in) :: rstate(97)
  integer,intent(in) :: istate(4)
  integer,intent(in) :: iff
  ix11 = istate(1)
  ix21 = istate(2)
  ix31 = istate(3)
  j1   = istate(4)
  R1   = rstate
  iff1 = iff
 end subroutine rand_vars_put

 !=============================================================
 !!function: rand_gen/rand_gen1
 !!génération de nombre aléatoire à partir du Numerical recipes
 real(dp) function rand_gen1(idum)

  integer,intent(inout) :: idum
  integer,parameter :: m1 = 259200,ia1 = 7141,ic1 = 54773
  integer,parameter :: m2 = 134456,ia2 = 8121,ic2 = 28411
  integer,parameter :: m3 = 243000,ia3 = 4561,ic3 = 51349
  real(dp),parameter :: rm1 = one/real(m1,dp),rm2 = one/real(m2,dp)

  if ((idum < 0).or.(iff1 == 0)) then
   iff1 = 1
   ix11 = mod(ic1-idum,m1)
   ix11 = mod(ia1*ix11+ic1,m1)
   ix21 = mod(ix11,m2)
   ix11 = mod(ia1*ix11+ic1,m1)
   ix31 = mod(ix11,m3)
   do j1 = 1,97
    ix11 = mod(ia1*ix11+ic1,m1)
    ix21 = mod(ia2*ix21+ic2,m2)
    R1(j1) = (real(ix11,dp)+real(ix21,dp)*rm2)*rm1
   enddo
   idum = 1
  end if

  ix11 = mod(ia1*ix11+ic1,m1)
  ix21 = mod(ia2*ix21+ic2,m2)
  ix31 = mod(ia3*ix31+ic3,m3)
  j1 = 1+(97*ix31)/m3
  if(j1 > 97.or.j1 < 1) stop 'ran1'
  rand_gen1 = R1(j1)
  R1(j1) = (real(ix11,dp)+real(ix21,dp)*rm2)*rm1

 end function rand_gen1


 !******************************** UNIF_DIST *******************************************
 subroutine unif_dist(irand,x1,x2,n1,n2,x)
  ! Génère une distribution aléatoire uniforme
  ! de nombre réels dans l'intervalle
  ! [x1,x2], retourné en éléments n1:n2 du tableau x
  use m_writeout

  integer,intent(in) :: n1,n2
  integer,intent(inout) :: irand
  real(dp) ,intent(in)   :: x1,x2
  real(dp),intent(inout) :: x(:)
  integer :: n
  real(dp) :: xr

  xr = x2-x1
  if(xr<zero) call wrt_double(qp_out,'error unif_dist : x2 < x1',wrtscreen,wrtdisk)

  do n=n1,n2
   x(n) = x1 + xr*rand_gen1(irand)
  enddo

 end subroutine unif_dist
 !***************************** END UNIF_DIST *******************************************

 !******************************** UNIF_DIST *******************************************
 subroutine unif_dist2(irand,x1,x2,x)
  ! Génère une distribution aléatoire uniforme
  ! de nombre réels dans l'intervalle
  ! [x1,x2], retourné en éléments n1:n2 du tableau x
  use m_writeout

  integer,intent(inout) :: irand
  real(dp) ,intent(in)   :: x1,x2
  real(dp),intent(inout) :: x(:)
  integer :: n
  real(dp) :: xr

  xr = x2-x1
  if(xr<zero) call wrt_double(qp_out,'error unif_dist : x2 < x1',wrtscreen,wrtdisk)

  do n=lbound(x,dim=1),ubound(x,dim=1)
   x(n) = x1 + xr*rand_gen1(irand)
  enddo

 end subroutine unif_dist2
 !***************************** END UNIF_DIST *******************************************




 !******************************** BI_MAX_DIB ********************************************
 subroutine bi_max_dib(irand,vth1,vth2,n1,n2,vx,vy,vz,vtHesvtH)
  ! Generates a Bi-maxwellian distribution
  ! In : irand - seed random integer, 
  !      vth - v(thermal), 
  !      n1:n2 - index range
  ! Out: arrays vx(:),vy(:),vz(:) of generated values

#ifdef HAVE_DEBUG
  use m_writeout,only   : wrt_debug
#endif

  !integer,intent(in) :: is
  integer,intent(inout) :: irand,n1,n2
  real(dp),intent(in) :: vth1,vth2,vtHesvtH
  real(dp),dimension(:),intent(inout) :: vx,vy,vz
  integer :: nn,irand1
  real(dp) :: vmag,theta

  __WRT_DEBUG_IN("bi_max_dib")

  if(irand<=0)  irand = 98354391

  ! if (is==1) then
  !  vtHesvtH_loc = one 
  ! else
  !  vtHesvtH_loc = vtHesvtH
  ! endif

  irand1 = 0
  do nn = n1,n2
   vmag  = sqrt(-log(one-0.99999_dp*rand_gen1(irand1)))
   theta = two_pi*rand_gen1(irand1)
   vy(nn) = vth2*vmag*cos(theta)
   vz(nn) = vth2*vmag*sin(theta)
   vmag = vtHesvtH*sqrt(-log(one-0.99999_dp*rand_gen1(irand1)))
   theta = two_pi*rand_gen1(irand1)
   vx(nn) = vth1*vmag*cos(theta)
  enddo

  __WRT_DEBUG_OUT("bi_max_dib")
 end subroutine bi_max_dib
 !******************************** END BMWDN **************************************

end module m_rand_gen
