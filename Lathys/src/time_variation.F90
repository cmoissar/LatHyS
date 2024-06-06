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
        integer     :: i, n_cut
        real(dp)    :: V_SW, V_MC
        real(dp)    :: B_SW, Bx_SW, By_SW, Bz_SW, B_MC, By_MC, Bz_MC, B_CME_i_2
        real(dp)    :: a, b, rate
        real(dp)    :: Ng_SW, Ng_MC
        real(dp)    :: vth_SW, vth_MC
        real(dp)    :: Pth_SW, Pth_CME_i
        real(dp)    :: rphi,rpsi
        integer     :: nwave
        real(dp)    :: kwave, omega
        type(species_type), intent(inout) :: species
        integer,parameter :: H=1

 __WRT_DEBUG_IN("define_VARIABLES_CME")

        t_start_MC = 10.0_dp
        n_start_MC = int(t_start_MC / dt)
        Dt_MC = 70.0_dp
        Dn_MC = int(Dt_MC / dt)
        t_trans_start = 5.0_dp
        n_trans_start = real(int(t_trans_start / dt)) ! Avoids integer division in the tanh
        t_end = 20.0_dp
        n_end = int(t_end / dt)
        t_trans_end = t_end + Dt_MC
        n_trans_end = real(int(t_trans_end / dt))

        !!! MAGNETIC FIELD !!!
        B_SW   = 1.
        Bx_SW  = bx0
        By_SW  = by0
        Bz_SW  = bz0
 
        B_MC   = (50.0_dp / 10.0_dp) * B_SW !B_SW !(12.0_dp / 3.75_dp) * B_SW
        By_MC  = B_MC 
        Bz_MC  = sqrt(B_MC**2 - By_MC**2)
 
        rpsi = psi*deg_to_rad
        rphi = phi*deg_to_rad     

        ! 1.*integer -> real(dp), otherwise tanh(.) is unhappy
        do i=1,n_start_MC
          By_CME(i) = By_SW+(By_MC*bessel_j1(-2.4)-By_SW) &
                                    *(1./2)*(1+tanh(1.*((i-(n_start_MC-1.25*n_trans_start/2.))*5./n_trans_start)))
          Bz_CME(i) = Bz_SW+(Bz_MC*bessel_j0(-2.4)-Bz_SW) &
                                    *(1./2)*(1+tanh(1.*((i-(n_start_MC-1.25*n_trans_start/2.))*5./n_trans_start)))
        enddo
        
        a = 2.4/Dn_MC
        b = -a*(Dn_MC+n_start_MC)

        do i=n_start_MC,nhm
          By_CME(i) = B_MC*bessel_j1(a*i+b) 
          Bz_CME(i) = B_MC*bessel_j0(a*i+b)
        enddo

        !!! VELOCITY !!!
        V_SW = species%S(H)%vxs  
        V_MC = 750.0_dp/species%ref%alfvenspeed

        ! 1.0_dp*integer -> real(dp), otherwise tanh(.) is unhappy
        ! the factor 1.25 is fine_tuning so that the full velocity is reached by
        ! n_start_MC
        do i=1,nhm
          V_CME(i) = V_SW+(V_MC-V_SW)*(1./2)*(1+tanh(1.*((i-(n_start_MC-1.25*n_trans_start/2.))*5./n_trans_start))) &
                         +(V_SW-V_MC)*(1./2)*(1+tanh(1.*((i-(n_end  +n_trans_end  /2.))*5./n_trans_end)))
        enddo

        !do i=1,nhm
        !   B = sqrt(Bx_SW**2 + By_CME(i)**2 + Bz_CME(i)**2)
        !   rate = (B-B_SW) / (B_MC-B_SW)
        !   V_CME(i) = V_SW + rate*(V_MC-V_SW)
        !enddo

        !!! DENSITY !!!
        Ng_SW = species%S(H)%ng
        Ng_MC = (1.7_dp / 1.7_dp) * Ng_SW

        ! the factor 1.25 is fine_tuning so that the full density is reached by
        do i=1,nhm
          Ng_CME(i) = Ng_SW+(Ng_MC-Ng_SW)*(1./2)*(1+tanh(1.*((i-(n_start_MC-1.25*n_trans_start/2.))*5./n_trans_start))) &
                           +(Ng_SW-Ng_MC)*(1./2)*(1+tanh(1.*((i-(n_end  +n_trans_end  /2.))*5./n_trans_end)))
        enddo
        
        !!! TEMPERATURE !!!
        vth_SW = species%S(H)%vth1
        ! The MC should be much colder than the solar wind
        vth_MC = vth_SW * sqrt(1100.0_dp/11000.0_dp)

        ! the factor 1.25 is fine_tuning so that the full temperature is reached
        do i=1,nhm
          vth1_CME(i) = vth_SW +(vth_MC-vth_SW)*(1./2)*(1+tanh(1.*((i-(n_start_MC-1.25*n_trans_start/2.))*5./n_trans_start))) & 
                               +(vth_SW-vth_MC)*(1./2)*(1+tanh(1.*((i-(n_end+n_trans_end  /2.))*5./n_trans_end)))
          vth2_CME(i) = vth1_CME(i)
        enddo

        !!! FLAT !!!
        n_cut = int(120.0_dp / dt)
        do i=1,nhm
          if (i>n_cut) then
            V_CME(i) = V_CME(n_cut)
            By_CME(i) = By_CME(n_cut)
            Bz_CME(i) = Bz_CME(n_cut)
            Ng_CME(i) = Ng_CME(n_cut)
            vth1_CME(i) = vth1_CME(n_cut)
            vth2_CME(i) = vth2_CME(n_cut)
         endif
       enddo
        
 __WRT_DEBUG_OUT("define_VARIABLES_CME")


end subroutine define_VARIABLES_CME

subroutine Init_change_CME(species)

type(species_type),intent(inout) :: species

 __WRT_DEBUG_IN("Init_change_CME")

        allocate(Ng_CME(nhm))
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
integer,parameter :: H=1


 __WRT_DEBUG_IN("temporal_change_CME")

 call add_waves(iter,By_alfven,Bz_alfven,By1_alfven,Bz1_alfven,Vy_alfven,Vz_alfven)

 by0 = By_CME(iter) + By_Alfven
 bz0 = Bz_CME(iter) + Bz_Alfven

 by1 = By_CME(iter) + By1_Alfven
 bz1 = Bz_CME(iter) + Bz1_Alfven

 species%S(:)%vxs = V_CME(iter)   
 species%S(:)%vys = Vy_alfven
 species%S(:)%vzs = Vz_alfven

 species%S(:)%ng  = Ng_CME(iter)
 
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
