!!=============================================================
!!=============================================================
!!module: atm_production
!! NAME
!!  atm_production (SHess)
!!
!!
!! FUNCTION
!!  Contains generic routines for the phtoproduction computations
!!
!! NOTE
module atm_photoproduction

 use defs_basis
 use defs_species
 use defs_atmospheretype
! use defs_parametre
 use m_writeout

#include "q-p_common.h"

 private::                   &
      &    prod_ionisation , & 
      & finalize_absorption, &
      & absorption_EUV,      &
      & shadow     

 public::                     &
      photoproduction_generic,&
      flux_solaire_generic
contains
!==========================================================================

 subroutine photoproduction_generic(Spe,ncm,gstep,s_min_loc,atmosphere)
! use defs_parametre
   
  integer,intent(in) :: ncm(3)
  real(dp),intent(in) :: gstep(3),s_min_loc(3)
  type(species_type), intent(in) :: Spe
  type(atmosphere_type),intent(inout) ::atmosphere

  !--Definition des parametres de la simulations
  integer :: nb_lo ! nombre de longueurs d'ondes
  real(dp) :: s_cen(3),freq_ionis
  character(len=300) :: msg
  !--Profondeur optique, depend de la position (x,y,z) et de la longueur d'ondde
  real(dp),dimension(:,:,:,:),allocatable :: optical_depth 

  __WRT_DEBUG_IN("photoproduction")
  !--Relative (to proc) distance from the planet
  
  write(msg,'(2a)')ch10,&
       " ___________________ Photo-Production  _______________"
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk) 

  nb_lo=atmosphere%nb_lo
  s_cen = Spe%P%centr-s_min_loc

  !--------------------------------------------------------------------------------------
  ! Calcul de l'absoprtion du flux solaire
  !--Optical Depth
  allocate(optical_depth(nb_lo,ncm(1),ncm(2),ncm(3))); optical_depth = zero
  do i=1,atmosphere%n_species
        if(atmosphere%species(i)%opaque) then 
               call absorption_EUV(nb_lo,ncm,s_cen,gstep,Spe,atmosphere,i,optical_depth)
        endif
  enddo
  call shadow(nb_lo,ncm,s_cen,gstep,Spe,optical_depth)
  call finalize_absorption(nb_lo,ncm,optical_depth)  
  !--------------------------------------
  !--Calcul des taux de photoproduction
  ! les taux de photoproduction sont en cm^-3.s^-1
  do i=1,atmosphere%n_pp
        atmosphere%photo_reactions(i)%daughter%prod=zero
  enddo
  do i=1,atmosphere%n_pp-atmosphere%n_pp_freq_fixed
  call prod_ionisation(nb_lo,ncm,atmosphere%photo_reactions(i)%mother%density,optical_depth,&
             &                   atmosphere%photo_reactions(i)%ion_react,atmosphere%EUVFLX,&
             &                   atmosphere%photo_reactions(i)%daughter%prod,freq_ionis)
  atmosphere%photo_reactions(i)%frequency=freq_ionis
  !print *,'Freq ionisation',atmosphere%photo_reactions(i)%daughter%name,freq_ionis
  enddo
  !-- If the photoionisation frequency is not calculated but given
  if (atmosphere%n_pp_freq_fixed /= 0) then
    do i=atmosphere%n_pp-atmosphere%n_pp_freq_fixed+1,atmosphere%n_pp !A VERIFIER
       call prod_ionisation_freq_fixed(nb_lo,ncm,atmosphere%photo_reactions(i)%mother%density,optical_depth,&
             &                   atmosphere%photo_reactions(i)%daughter%prod,freq_ionis)
    enddo
  endif
  
  !--Desallocation des tableaux de photo absorption, de photo sionisation et de flux solaire
  deallocate(optical_depth)
  
  __WRT_DEBUG_OUT("photoproduction")
 end subroutine photoproduction_generic

 !---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Calcul de cone d'ombre
!---------------------------------------------------------------------------
subroutine shadow(nb_lo,ncm,s_cen,gstep,Spe,optical_depth)
! use defs_parametre
  
implicit none
integer,intent(in) :: nb_lo
integer,dimension(3),intent(in)::ncm
real(dp),dimension(3),intent(in) ::gstep,s_cen !necessary to compute the shadow
type(species_type),intent(in) :: Spe
real(dp),dimension(nb_lo,ncm(1),ncm(2),ncm(3)),intent(inout) :: optical_depth 

real :: cst,radius,rb,rp ! used for internal computations
integer :: i,j,k ! used for internal computations
  __WRT_DEBUG_IN("shadow")

 rp=Spe%P%r_lim**2 !better not to compute the square roots      
 cst=0.5*gstep(1)*Spe%ref%c_omegapi*1.e+5 !km to cm
  do k=1,ncm(3)
    do j = 1,ncm(2)
      rb = (float(j-1)*gstep(2)-s_cen(2))**2+(float(k-1)*gstep(3)-s_cen(3))**2
      do i = 2,ncm(1)  
        ! we estimate the optical depth only if we are outside of the obstacle region, except in the shadow of the planet
          radius = (float(i-1)*gstep(1)-s_cen(1))**2+rb
          if (.NOT.((radius > rp).and.((rb > rp).or.(float(i-1)*gstep(1) <= s_cen(1))))) then
                optical_depth(1:nb_lo,i,j,k) = 1. ! means a 0 optical_depth, see finalize_absorption
        endif
      enddo
    enddo
  enddo
  __WRT_DEBUG_OUT("shadow")

return
end subroutine shadow


!---------------------------------------------------------------------------
! Calcul de l'absorption du flux solaire
!---------------------------------------------------------------------------
subroutine absorption_EUV(nb_lo,ncm,s_cen,gstep,Spe,atmosphere,n_spe,optical_depth)
! use defs_parametre
  
implicit none
integer,intent(in) :: nb_lo,n_spe
integer,dimension(3),intent(in)::ncm
real(dp),dimension(3),intent(in) ::gstep,s_cen !necessary to compute the shadow
type(species_type),intent(in) :: Spe
type(atmosphere_type),intent(in) ::atmosphere

real(dp),dimension(nb_lo,ncm(1),ncm(2),ncm(3)),intent(inout) :: optical_depth 

real(dp) :: cst,radius,rb,rp,tmp(1:nb_lo) ! used for internal computations
integer :: i,j,k ! used for internal computations
  __WRT_DEBUG_IN("absorption_EUV")

 rp=Spe%P%r_lim**2 !better not to compute the square roots      
 cst=0.5*gstep(1)*Spe%ref%c_omegapi*1.e+5 !km to cm
  do k=1,ncm(3)
    do j = 1,ncm(2)
      rb = (float(j-1)*gstep(2)-s_cen(2))**2+(float(k-1)*gstep(3)-s_cen(3))**2
      tmp=0._dp
      do i = 2,ncm(1)  
        ! we estimate the optical depth only if we are outside of the obstacle region, except in the shadow of the planet
          radius = (float(i-1)*gstep(1)-s_cen(1))**2+rb
          if ((radius > rp).and.((rb > rp).or.(float(i-1)*gstep(1) <= s_cen(1)))) then
                tmp(1:nb_lo)=tmp(1:nb_lo)+cst*atmosphere%species(n_spe)%ion_abs(:)*sum(atmosphere%species(n_spe)%density(i-1:i,j,k))
                optical_depth(1:nb_lo,i,j,k) = optical_depth(1:nb_lo,i,j,k)-tmp 
! minus is not an error, must be negative,see finalize_absorption
        else
                optical_depth(1:nb_lo,i,j,k) = 1. ! means a 0 optical_depth, see finalize_absorption
        endif
      enddo
    enddo
  enddo
  __WRT_DEBUG_OUT("absorption_EUV")
return
end subroutine absorption_EUV

subroutine finalize_absorption(nb_lo,ncm,optical_depth)
! use defs_parametre
 use defs_species
 integer,intent(in) :: nb_lo
 integer,dimension(3),intent(in)::ncm
 real(dp),dimension(nb_lo,ncm(1),ncm(2),ncm(3)),intent(inout) :: optical_depth 
  __WRT_DEBUG_IN("finalize_absorption")
 where (optical_depth.gt.(-0.05)) !compute the exponantial only if 0.95<exp(-opt_depth)<0 otherwise set to 1-opt_depth or 0.
         optical_depth=1.-abs(optical_depth) ! much faster than exponential
 elsewhere 
        optical_depth=exp(optical_depth)
 end where 
  __WRT_DEBUG_OUT("finalize_absorption")
end subroutine finalize_absorption
!------------------- fin ABSORPTION_EUV ---------------------------

 !---------------------------------------------------------------------------
 ! Calcul du taux de production pour une espece pour une espece neutre donnee
 !---------------------------------------------------------------------------
 !!=============================================================
 !!routine: pp_generic_func/prod_ionisation
 !!
 !! FUNCTION
 !!  Compute Production Rate for a given pair species-neutral species.
 !!  Used for any Planet environnement.
 !!  Production rate are computed in (cm^-3/s).
 !!  
subroutine prod_ionisation(nb_lo,ncm,density,optical_depth,ion,EUVFLX,prod,freq_ionis)
! use defs_parametre
 !========================
implicit none

integer,intent(in) :: nb_lo
integer,dimension(3),intent(in) :: ncm
real(dp),dimension(nb_lo,ncm(1),ncm(2),ncm(3)),intent(in) :: optical_depth
real(dp),dimension(nb_lo),intent(in) :: EUVFLX,ion
real(dp),dimension(ncm(1),ncm(2),ncm(3)),intent(in):: density
real(dp),dimension(ncm(1),ncm(2),ncm(3)),intent(inout) :: prod
real(dp),intent(inout) :: freq_ionis
real(dp),dimension(nb_lo) :: temp

integer :: j,k,i
  __WRT_DEBUG_IN("prod_ionisation")

        
        temp(:) = zero
        temp(1:nb_lo) = EUVFLX(1:nb_lo)*ion(1:nb_lo)
        freq_ionis = sum(temp)
  do j = 1,ncm(2)
    do k = 1,ncm(3)
!compute the exact ionization frequency only if different by more than 1% from the optically thin one at x=nx
        if (sum(optical_depth(1:nb_lo,ncm(1),j,k)).lt.0.99) then 
            prod(1:ncm(1),j,k) = prod(1:ncm(1),j,k)+MATMUL(temp,optical_depth(:,1:ncm(1),j,k))*density(1:ncm(1),j,k)
         else
            prod(1:ncm(1),j,k) = prod(1:ncm(1),j,k)+freq_ionis*density(1:ncm(1),j,k)
        endif
     enddo
  enddo
  __WRT_DEBUG_OUT("prod_ionisation")

return
end subroutine prod_ionisation
!------------------------ fin PROD_IONISATION --------------------------------

 !---------------------------------------------------------------------------
 ! Calcul du taux de production pour une espece pour une espece neutre donnee
 !---------------------------------------------------------------------------
 !!=============================================================
 !!routine: pp_generic_func/prod_ionisation
 !!
 !! FUNCTION
 !!  Compute Production Rate for a given pair species-neutral species.
 !!  Used for any Planet environnement.
 !!  Production rate are computed in (cm^-3/s).
 !!  
subroutine prod_ionisation_freq_fixed(nb_lo,ncm,density,optical_depth,prod,freq_ionis)
 use defs_parametre
 !========================
implicit none

integer,intent(in) :: nb_lo
integer,dimension(3),intent(in) :: ncm
real(dp),dimension(nb_lo,ncm(1),ncm(2),ncm(3)),intent(in) :: optical_depth
real(dp),dimension(ncm(1),ncm(2),ncm(3)),intent(in):: density
real(dp),dimension(ncm(1),ncm(2),ncm(3)),intent(inout) :: prod
real(dp),intent(in) :: freq_ionis


integer :: i,j,k

do i=1,ncm(1)
  do j = 1,ncm(2)
    do k = 1,ncm(3)
!compute the exact ionization frequency only if different by more than 1% from the optically thin one at x=nx
            prod(i,j,k) = prod(i,j,k)+freq_ionis*density(i,j,k)*sum(optical_depth(1:nb_lo,i,j,k))
     enddo
  enddo
enddo  

return
end subroutine prod_ionisation_freq_fixed
!------------------------ fin PROD_IONISATION_FREQ_FIXED --------------------------------
 !---------------------------------------------------------------------------
 ! Calcul du flux solaire (photons/cm2/s) modele EUVAC
 !---------------------------------------------------------------------------
 subroutine flux_solaire_generic(nb_lo,atmosphere,F107,F107A,dist_conv)
! use defs_parametre
 
  ! This EUV flux model uses the F74113, solar reference spectrum and ratios determined
  ! from Hintegger's SERF1 model. It uses the daily F10.7 flux (F107) and the 81 day mean
  ! (F107A) as a proxy for scaling? The fluxes are returned in EUVFLX and correspond to the 
  ! 37 wavelength bins of Torr et al. [1979] Geophys. Res. Lett. p771
  ! see Richards et al. [1994] J. Geophys. Res. p8981 for details

  ! F107   = input daily 10.7 cm flux index
  ! F107A  = input 81 day average of daily F10.7 centered on current day
  ! EUVFLX = output array for EUV flux in units of photons/cm2/s
  ! Changed units in photons/(cm^2*s)

  integer,intent(in) :: nb_lo
  type(atmosphere_type),intent(inout) ::atmosphere
  real(dp),intent(in) :: F107,F107a,dist_conv

  character(len=300) :: msg
  real(dp),dimension(nb_lo) :: AFAC,F74113,EUVFLX
  integer :: ii
  real(dp) :: flxfac

  
  __WRT_DEBUG_IN("flux_solaire")

  write(msg,'(3x,a)')"-Solar Flux"
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  AFAC = &
       [ 1.0017e-002, 7.1250e-003, 1.3375e-002, 1.9450e-002, 2.7750e-003 &
       , 1.3768e-001, 2.6467e-002, 2.5000e-002, 3.3333e-003, 2.2450e-002 &
       , 6.5917e-003, 3.6542e-002, 7.4083e-003, 7.4917e-003, 2.0225e-002 &
       , 8.7583e-003, 3.2667e-003, 5.1583e-003, 3.6583e-003, 1.6175e-002 &
       , 3.3250e-003, 1.1800e-002, 4.2667e-003, 3.0417e-003, 4.7500e-003 &
       , 3.8500e-003, 1.2808e-002, 3.2750e-003, 4.7667e-003, 4.8167e-003 &
       , 5.6750e-003, 4.9833e-003, 3.9417e-003, 4.4167e-003, 5.1833e-003 &
       , 5.2833e-003, 4.3750e-003]

  F74113 = &
       [  1.200,  0.450,  4.800,  3.100,  0.460 &
       ,  0.210,  1.679,  0.800,  6.900,  0.965 & 
       ,  0.650,  0.314,  0.383,  0.290,  0.285 & 
       ,  0.452,  0.720,  1.270,  0.357,  0.530 &
       ,  1.590,  0.342,  0.230,  0.360,  0.141 &
       ,  0.170,  0.260,  0.702,  0.758,  1.625 &
       ,  3.537,  3.000,  4.400,  1.475,  3.500 &
       ,  2.100,  2.467]

  do ii = 1,nb_lo
   flxfac = (one+AFAC(ii)*(half*(F107+F107A) - 80.0))
   if (flxfac < 0.8) flxfac = 0.8
   EUVFLX(ii) = F74113(ii)*flxfac*1.e+9
  end do

  ! Flux solaire a l'orbite de l'objet
  !--(1.e4 conversion to photons/(m2*s))
  atmosphere%EUVFLX = EUVFLX/dist_conv
  
  __WRT_DEBUG_OUT("flux_solaire")
 end subroutine flux_solaire_generic
 !------------------- fin de FLUX_SOLAIRE ---------------------------

end module atm_photoproduction
