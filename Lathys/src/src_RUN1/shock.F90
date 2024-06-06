!!=============================================================
!!=============================================================
module shock

 use defs_basis
 use m_writeout,only  : wrt_debug
 use defs_parametre
 use m_rand_gen,only    : rand_gen1

#include "q-p_common.h"

 implicit none
 private

 public::                &
      shock_position,    &
      refresh

contains
 !!#####################################################################


 !
 !****************************** SHOCK_POSITION **************************************
 subroutine shock_position()
  !========================
  ! Determination de la position du choc
  use defs_grid,only        : nc1_tot
  use defs_variable,only    : pos_choc ,vela,dna,Spe

  integer :: j,k,i
  real(dp) :: saut
  integer :: found

  __WRT_DEBUG_IN("shock_position")

  found = 0
  i = 1
  ! initialisation
  pos_choc(:,:) = nc1_tot(1) !nx

  saut = Spe%S(1)%vxs-one

  do k = 1,nc1_tot(3)
   do j = 1,nc1_tot(2)
    do while ((found == 0).and.(i < nc1_tot(1)))
     if ((vela%x(i,j,k)/dna(i,j,k)) <= (saut)) then
      found = 1
      pos_choc(j,k) = i  ! on releve l'indice X pour laquelle le saut de vitesse
     else      ! est remplie
      i = i+1
     endif
    enddo
    found = 0
    i = 1
   enddo
  enddo

  __WRT_DEBUG_OUT("shock_position")
 end subroutine shock_position
 !*********************************** END SHOCK_POSITION ***********************

 !********************************  REFRESH *********************************
 subroutine refresh
  !========================

  ! written by R. Modolo 07/01/03
  ! On passe sur toutes les particules et on regarde si la particule se trouve 
  ! à l'intérieur d'une pseudo courbe de choc. On re-numérote   les particules
  use defs_variable,only : pos_choc ,vela,dna,nptot,&
       &                 particule,gstep,Spe!vxs,&
!       &                 vth1,vth2,vys,vzs,vtHesvtH

  integer  :: i,j,k,n,compt,is,iran
  real(dp) :: vmag,theta

  __WRT_DEBUG_IN("refresh")

  write(*,*) '****************'
  write(*,*) '     refresh '
  write(*,*) 
  write(qp_out,*) '****************'
  write(qp_out,*) '     refresh '
  write(qp_out,*) 

  write(*,*) 'max(pos_choc) ',maxval(pos_choc),minval(pos_choc), &
       minloc(pos_choc)
  write(qp_out,*) 'max(pos_choc) ',maxval(pos_choc),minval(pos_choc), &
       minloc(pos_choc)

  iran = 0
  compt = 0
  do n = 1,nptot
   i = int(particule(n)%pos(1)/gstep(1))+1
   j = int(particule(n)%pos(2)/gstep(2))+1
   k = int(particule(n)%pos(3)/gstep(3))+1
   if((pos_choc(j,k) <= (i+10)) .or. &
        & particule(n)%exc*real(particule(n)%orig,dp)/=zero) continue

   ! if ((pos_choc(j,k) > (i+10)).and.(particule(9,n) == 0).and. &
   !      (particule(10,n) == 0)) then

   ! On se trouve en amont du choc et on retire les vitesses
   ! si la particule est une particule issu du vent solaire qui n'a pas subit 
   ! d'echange de charge
   ! on regarde si il s'agit de proton ou d'hélium
   if (particule(n)%mass/particule(n)%char == one) then 
    is = 1
   else
    is = 2
   endif
   compt = compt+1
   vmag  = sqrt(-log(one-0.99999_dp*rand_gen1(iran)))
   theta = two_pi*rand_gen1(iran)
   particule(n)%vel(2) = Spe%S(is)%vth2*vmag*cos(theta)+Spe%S(is)%vys
   particule(n)%vel(3) = Spe%S(is)%vth2*vmag*sin(theta)+Spe%S(is)%vzs

   vmag = Spe%S(is)%rspeed*sqrt(-log(one-0.99999_dp*rand_gen1(iran)))
   ! if (is == 1) then
   !  vmag  = sqrt(-log(one-0.99999_dp*rand_gen1(iran)))
   ! else
   !  vmag = vtHesvtH*sqrt(-log(one-0.99999_dp*rand_gen1(iran)))
   ! endif
   theta = two_pi*rand_gen1(iran)
   particule(n)%vel(1) = Spe%S(is)%vth1*vmag*cos(theta)+ Spe%S(is)%vxs
   !endif
  enddo

  write(*,*) 'Le nombre de particule rafraichi est  : ',compt
  write(qp_out,*) 'Le nombre de particule rafraichi est  : ',compt

  write(*,*) '****************'
  write(qp_out,*) '****************'

  __WRT_DEBUG_OUT("refresh")
 end subroutine refresh
 !*********************************** END REFRESH **********************************

end module shock

