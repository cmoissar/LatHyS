module time_variation

use defs_basis
use defs_parametre
use defs_species
use defs_variable
use m_writeout
#ifdef DIntelFortran
use ifport
#endif

#include "q-p_common.h"

implicit none
#ifdef PGIFortran
#include "lib3f.h"
#endif

private 

public :: &
        init_change_CME,      &!--initialization of temporal change in the incident plasma
        temporal_change_CME    !--apply change in the time loop


contains

subroutine define_VARIABLES_CME(species)


        real(dp)    :: t_start_MC, t_trans_start, Dt_MC, t_end, t_trans_end
        integer     :: n_start_MC, n_trans_start, Dn_MC, n_end, n_trans_end  
        real(dp)    :: V_SW, V_MC
        real(dp)    :: B_SW, Bx_SW, By_SW, Bz_SW, B_MC, By_MC, Bz_MC, B_CME_i_2
        real(dp)    :: a, b, rate
        real(dp)    :: Ng_SW, Ng_MC, Ng_H0, Ng_He0
        real(dp)    :: vth_SW, vth_MC
        real(dp)    :: Pth_SW, Pth_CME_i
        real(dp)    :: rphi,rpsi
        integer     :: nwave
        real(dp)    :: kwave, omega
        type(species_type), intent(inout) :: species
        integer :: i
        integer, parameter :: H=1
        integer, parameter :: He=2

 __WRT_DEBUG_IN("define_VARIABLES_CME")

        t_start_MC = 300.0_dp
        n_start_MC = int(t_start_MC / dt)
        t_trans_start = 4.0_dp
        n_trans_start = real(int(t_trans_start / dt)) ! Avoids integer division in the tanh

        !!! MAGNETIC FIELD !!!
        B_SW   = 1.0_dp   !ref%mag is defined in env_<planet>.F90. 
        Bx_SW  = bx0      !B_SW = 1 essentially means B_SW = ref%mag
        By_SW  = by0
        Bz_SW  = bz0
 
        B_MC   = bx0
        By_MC  = -by0 
        Bz_MC  = bz0
 
        rpsi = psi*deg_to_rad
        rphi = phi*deg_to_rad     

        ! 1.*integer -> real(dp), otherwise tanh(.) is unhappy
        do i=1,nhm
          By_CME(i) = By_SW+(By_MC-By_SW)*(1./2)*(1+tanh(1.*((i-(n_start_MC-1.25*n_trans_start/2.))*5./n_trans_start)))
          Bz_CME(i) = Bz_SW+(Bz_MC-Bz_SW)*(1./2)*(1+tanh(1.*((i-(n_start_MC-1.25*n_trans_start/2.))*5./n_trans_start)))
        enddo

        !!! VELOCITY !!!
        V_SW = species%S(H)%vxs  
        V_MC = V_SW

        ! 1.0_dp*integer -> real(dp), otherwise tanh(.) is unhappy
        ! the factor 1.25 is fine_tuning so that the full velocity is reached by n_start_MC
        do i=1,nhm
          V_CME(i) = V_SW+(V_MC-V_SW)*(1./2)*(1+tanh(1.*((i-(n_start_MC-1.25*n_trans_start/2.))*5./n_trans_start)))
        enddo

        !!! DENSITY !!!
        Ng_SW = 1
        Ng_MC = Ng_SW
        Ng_H0 = species%S(H)%ng   !This is ugly. TODO: Use Ng_SW = species%(:)%ng , or something like that instead
        Ng_He0 = species%S(He)%ng


        ! the factor 1.25 is fine_tuning so that the full density is reached by n_start_MC
        do i=1,nhm
          Ng_CME(i) = Ng_SW+(Ng_MC-Ng_SW)*(1./2)*(1+tanh(1.*((i-(n_start_MC-1.25*n_trans_start/2.))*5./n_trans_start))) 
          Ng_H(i) = Ng_CME(i)*Ng_H0
          Ng_He(i) = Ng_CME(i)*Ng_He0
        enddo
        
        !!! TEMPERATURE !!!
        vth_SW = species%S(H)%vth1
        vth_MC = vth_SW

        ! the factor 1.25 is fine_tuning so that the full temperature is reached by n_start_MC
        do i=1,nhm
          vth1_CME(i) = vth_SW +(vth_MC-vth_SW)*(1./2)*(1+tanh(1.*((i-(n_start_MC-1.25*n_trans_start/2.))*5./n_trans_start)))
          vth2_CME(i) = vth1_CME(i)
        enddo

 __WRT_DEBUG_OUT("define_VARIABLES_CME")


end subroutine define_VARIABLES_CME

subroutine Init_change_CME(species)

type(species_type),intent(inout) :: species

 __WRT_DEBUG_IN("Init_change_CME")

        allocate(Ng_CME(nhm))
        allocate(Ng_H(nhm))
        allocate(Ng_He(nhm))
        allocate(vth1_CME(nhm))
        allocate(vth2_CME(nhm))
        allocate(By_CME(nhm))
        allocate(Bz_CME(nhm))
        allocate(By1_CME(nhm))
        allocate(Bz1_CME(nhm))
        allocate(V_CME(nhm))    
        allocate(Vy_alfven(nhm))
        allocate(Vz_alfven(nhm))   
        call define_VARIABLES_CME(species)
 

 __WRT_DEBUG_OUT("Init_change_CME")

end subroutine Init_change_CME


subroutine temporal_change_CME(by0,bz0,by1,bz1,vxmean,species,iter)

use particle_fluxes
use field_add_waves, only :  add_waves

real(dp) :: By_alfven, Bz_alfven, By1_alfven, Bz1_alfven
real(dp) :: Vy_alfven, Vz_alfven
real(dp)           , intent(inout) :: by0,bz0,by1,bz1,vxmean
type(species_type) , intent(inout) :: species
integer            , intent(in)    :: iter
integer                            :: is
integer, parameter :: H=1
integer, parameter :: He=2


 __WRT_DEBUG_IN("temporal_change_CME")

 call add_waves(iter,By_alfven,Bz_alfven,By1_alfven,Bz1_alfven,Vy_alfven,Vz_alfven)

 by0 = By_CME(iter) + By_Alfven
 bz0 = Bz_CME(iter) + Bz_Alfven

 by1 = By_CME(iter) + By1_Alfven
 bz1 = Bz_CME(iter) + Bz1_Alfven

 species%S(:)%vxs = V_CME(iter)   
 species%S(:)%vys = Vy_alfven
 species%S(:)%vzs = Vz_alfven

 species%S(H)%ng  = Ng_H(iter)
 species%S(He)%ng = Ng_He(iter)
 
 species%S(:)%vth1 = vth1_CME(iter)
 species%S(:)%vth2 = vth2_CME(iter)

 do is = 1,ns
   call compute_fluxes(species,is,nfl)
 enddo

 !--Used to update econv (E at the left boundary of the box)
 vxmean = V_CME(iter)

 __WRT_DEBUG_OUT("temporal_change_CME")

end subroutine temporal_change_CME
end module time_variation
