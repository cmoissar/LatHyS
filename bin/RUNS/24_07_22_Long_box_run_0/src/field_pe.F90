!!=============================================================
!!=============================================================
!!module: m_pecalc
!! NAME
!!  m_pecalc (RModolo,MMancini)
!!
!! FUNCTION
!!  Contains routine to compute electronic pressure
!!
!! NOTE
module field_pe

 use defs_basis
 use m_writeout
 use m_timing,only       : time_get
#include "q-p_common.h"

 implicit none
 private
 

 public  ::              &
      pecalc              !--Compute electronic pressure

contains
 !!#########################################################


 !!=============================================================
 !!routine: m_pecalc/pecalc
 !!
 !! FUNCTION
 !!  Compute electronic pressure
 !! IN 
 !!  ipe=
 !!  te=Normalised temperature
 !!  dn(:,:,:)=density field
 !!
 !! OUT
 !!
 !! SIDE EFFECT
 !!  pe(:,:,:)=electronic pressure field  
 !! NOTE
 !!
 subroutine pecalc(ipe,te,dn,pe)
 use defs_variable,only : dn_e_pl,dn_e_incdt,Spe,s_min_loc
 use defs_grid,only     : ncm
 use defs_parametre,only :gstep
  !--Electron pressure : pe = nkTe (Te = constant)
  !--Here TE is a normalised temperature (kTe) in simulation units.

  integer,intent(in) ::ipe
  real(dp),intent(in) :: te
  real(dp),dimension(:,:,:),intent(in) :: dn
  real(dp),dimension(:,:,:),intent(inout) :: pe
  real(dp)::rb,radius,s_cen(3)
  integer :: i,j,k  
  real(dp) :: gamma,gamma_iono

  __WRT_DEBUG_IN("pecalc")
  __GETTIME(21,1)!--Timer start
!--0 = adiabatic, 1 = isothermal, 2= adiabatic + hydrostatic ionosphere, 3= isothermal + hydrostatic ionosphere
if (mod(ipe,2) == 0) then
  gamma=5._dp/3._dp 
else 
  gamma=1._dp
endif
  do k=1,ncm(3)
          do j=1,ncm(2)
              rb = (float(j-1)*gstep(2)-s_cen(2))**2+(float(k-1)*gstep(3)-s_cen(3))**2
                  do i=1,ncm(1)
                        radius = (float(i-1)*gstep(1)-s_cen(1))**2+rb-(Spe%P%radius)**2
                        if (((radius.le.0).or.(dn_e_pl(i,j,k).le.(1.1_dp))).or.(ipe.lt.2)) then 
                                gamma_iono=gamma
                        else  
                                gamma_iono = 1._dp/log10(dn_e_pl(i,j,k))
                                gamma_iono = gamma_iono/sqrt(sqrt(1._dp+gamma_iono**4/(gamma)**4))
                        endif
        if ((dn_e_pl(i,j,k)+dn_e_incdt(i,j,k)).ne.0) then
        pe(i,j,k) = te*dn_e_incdt(i,j,k)**gamma + te*Spe%tempe_ratio*dn_e_pl(i,j,k)**gamma_iono
        else
        pe(i,j,k) = 0._dp
        endif
                enddo
        enddo
 enddo
  __GETTIME(21,2)!--Timer stop
  __WRT_DEBUG_OUT("pecalc")
 end subroutine pecalc
 !******************************** END PECALC *************************************

 end module field_pe
