!!=============================================================
!!=============================================================
!!module: atm_sections_efficaces
!! NAME
!!  atm_sections_efficaces
!!
!!Define cross section for photonionisation
!!
!! FUNCTION
!! Define cross section for photonionisation
!!
!! NOTE
module atm_sections_efficaces

 use defs_basis
 use m_writeout
 
#include "q-p_common.h"
 
 
 implicit none
 private
 
 
 public ::section_efficace_H_abs,  &
      section_efficace_O_abs,  &
      section_efficace_CO2_abs, &
      section_efficace_H_ion,  &
      section_efficace_O_ion,  &
      section_efficace_CO2_ion, &
      section_efficace_O_CO2_ion, &
      section_efficace_H2_abs, &
      section_efficace_H2_ion, &
      section_efficace_H2O_abs, &
      section_efficace_H2O_ion, &
      section_efficace_H_H2O_ion, &
      section_efficace_O_H2O_ion, &
      section_efficace_OH_H2O_ion, &
      section_efficace_O2_abs, &
      section_efficace_O2_ion, &
      section_efficace_N2_ion, &
      section_efficace_CH4_ion, &
      section_efficace_N2_abs, &
      section_efficace_CH4_abs
      
 contains
 
!---------------------------------------------------------------------------
! Definition of H cross section  for absorbtion
!---------------------------------------------------------------------------
 subroutine section_efficace_H_abs(nb_lo,H_abs)
!--Definition des valeurs des sections efficaces pour l'Hydrogen
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: H_abs
   real(dp),dimension(nb_lo) :: H_ab
   __WRT_DEBUG_IN("section_efficace_H_abs")
!--Section efficace pour l'absorption in (cm^2)
  H_ab = &
       [  0.0024 , 0.0169 , 0.0483 , 0.1007 , 0.1405 ,&
       &  0.1913 , 0.1676 , 0.2324 , 0.2334 , 0.3077 ,&
       &  0.4152 , 0.3984 , 0.6163 , 0.8387 , 0.9739 ,&
       &  1.199  , 1.419  , 1.662  , 1.62   , 1.888  ,&
       &  2.079  , 2.076  , 2.641  , 2.897  , 3.173  ,&
       &  3.73   , 3.807  , 4.093  , 3.868  , 4.784  ,&
       &  5.67   , 3.469  , 0.0    , 0.0    , 0.0    ,&
       &  0.0    , 0.0    ]
       H_abs = H_ab*1.e-18


  __WRT_DEBUG_OUT("section_efficace_H_abs")
  end subroutine section_efficace_H_abs
!------------------- end SECTION_EFFICACE_H_abs ------------------------

!---------------------------------------------------------------------------
! Definition of H2 cross section  for absorbtion
!---------------------------------------------------------------------------
 subroutine section_efficace_H2_abs(nb_lo,H2_abs)
!--Definition des valeurs des sections efficaces pour l'Hydrogen moleculaire
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: H2_abs
   real(dp),dimension(nb_lo) :: H2_ab
   __WRT_DEBUG_IN("section_efficace_H2_abs")
!--Section efficace pour l'absorption in (cm^2)
H2_ab = &
       [  0.011,  0.080,  0.208,  0.433,  0.604 &
       &,  0.839,  0.730,  1.018,  1.022,  1.417 &
       &,  1.942,  1.901,  3.025,  3.870,  4.502 &
       &,  5.356,  6.168,  7.021,  6.864,  7.811 &
       &,  8.464,  8.445,  9.900, 10.731, 11.372 &
       &, 10.755,  8.640,  7.339,  8.748,  8.253 &
       &,  0.476,  0.185,  0.000,  0.046,  0.000 &
       &,  0.000,  0.000]
       H2_abs = H2_ab*1.e-18

  __WRT_DEBUG_OUT("section_efficace_H2_abs")
  end subroutine section_efficace_H2_abs
!------------------- end SECTION_EFFICACE_H2_abs ------------------------


!---------------------------------------------------------------------------
! Definition of H2O cross section  for absorbtion
!---------------------------------------------------------------------------
 subroutine section_efficace_H2O_abs(nb_lo,H2O_abs)
!--Definition des valeurs des sections efficaces pour l'eau moleculaire
!  section efficace in (cm^2
   integer,intent(in) :: nb_lo
   real(dp),dimension(nb_lo),intent(inout) :: H2O_abs
   real(dp),dimension(nb_lo) :: H2O_ab
   __WRT_DEBUG_IN("section_efficace_H2O_abs")
!--Section efficace pour l'absorption in (cm^2)
H2O_ab = &
       [  0.699,  1.971,  4.069,  6.121,  7.52 &
       &,  8.934,  8.113, 9.907,  9.93,  11.35 &
       &,  13.004,  12.734,  16.032,  18.083,  18.897 &
       &,  20.047,  21.159,  21.908,  21.857,  22.446 &
       &,  22.487,  22.502,  22.852, 22.498, 22.118 &
       &, 19.384,  20.992,  16.975, 18.151 , 16.623 &
       &,  19.837,  20.512,  15.072,  15.176,  18.069 &
       &,  15.271,  8.001]
       H2O_abs = H2O_ab*1.e-18

  __WRT_DEBUG_OUT("section_efficace_H2O_abs")
  end subroutine section_efficace_H2O_abs
!------------------- end SECTION_EFFICACE_H2O_abs ------------------------



!---------------------------------------------------------------------------
! Definition of O2 cross section  for absorbtion
!---------------------------------------------------------------------------
 subroutine section_efficace_O2_abs(nb_lo,O2_abs)
!--Definition des valeurs des sections efficaces pour l'oxygen moleculaire
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: O2_abs
   real(dp),dimension(nb_lo) :: O2_ab
   __WRT_DEBUG_IN("section_efficace_O2_abs")
!--Section efficace pour l'absorption in (cm^2)
O2_ab = &
       [    1.316,   3.806,   7.509,  10.900,  13.370 &
       &,  15.790,  14.387,  16.800,  16.800,  17.438 &
       &,  18.320,  18.118,  20.310,  21.910,  23.101 &
       &,  24.606,  26.040,  22.720,  26.610,  28.070 &
       &,  32.060,  26.017,  21.919,  27.440,  28.535 &
       &,  20.800,  18.910,  26.668,  22.145 , 16.631 &
       &,   8.562,  12.817,  18.730,  21.108,   1.630 &
       &,   1.050,   1.346]
       O2_abs = O2_ab*1.e-18

  __WRT_DEBUG_OUT("section_efficace_O2_abs")
end subroutine section_efficace_O2_abs
!------------------- end SECTION_EFFICACE_O2_abs ------------------------



!---------------------------------------------------------------------------
! Definition of N2 cross section  for absorbtion
!---------------------------------------------------------------------------
 subroutine section_efficace_N2_abs(nb_lo,N2_abs)
!--Definition des valeurs des sections efficaces pour l'azote moleculaire
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
    real(dp),dimension(nb_lo),intent(inout) :: N2_abs
   real(dp),dimension(nb_lo) :: N2_ab
   __WRT_DEBUG_IN("section_efficace_N2_abs")
!--Section efficace pour l'absorption in (cm^2)
N2_ab = &
       [    0.720,   2.261,   4.958,   8.392,  10.210 &
       &,  10.900,  10.493,  11.670,  11.700,  13.857 &
       &,  16.910,  16.395,  21.675,  23.160,  23.471 &
       &,  24.501,  24.130,  22.400,  22.787,  22.790 &
       &,  23.370,  23.339,  31.755,  26.540,  24.622 &
       &,  120.490, 14.180,  16.487,  33.578,  16.992 &
       &,   20.249,  9.680,   2.240,  50.988,   0.000 &
       &,    0.000,  0.000]
       N2_abs = N2_ab*1.e-18

  __WRT_DEBUG_OUT("section_efficace_N2_abs")
  end subroutine section_efficace_N2_abs
!------------------- end SECTION_EFFICACE_N2_abs ------------------------



!---------------------------------------------------------------------------
! Definition of CH4 cross section  for absorbtion
!---------------------------------------------------------------------------
 subroutine section_efficace_CH4_abs(nb_lo,CH4_abs)
!--Definition des valeurs des sections efficaces pour la méthane moleculaire
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
    real(dp),dimension(nb_lo),intent(inout) :: CH4_abs
   real(dp),dimension(nb_lo) :: CH4_ab
   __WRT_DEBUG_IN("section_efficace_CH4_abs")
!--Section efficace pour l'absorption in (cm^2)
CH4_ab = &
       [    0.204,   0.593,   1.496,   2.794,   3.857 &
       &,   5.053,   4.360,   6.033,   6.059,   7.829 &
       &,  10.165,   9.776,  14.701,  18.770,  21.449 &
       &,  24.644,  27.924,  31.052,  30.697,  33.178 &
       &,  35.276,  34.990,  39.280,  41.069,  42.927 &
       &,  45.458,  45.716,  46.472,  45.921,  48.327 &
       &,  48.968,  48.001,  41.154,  38.192,  32.700 &
       &,  30.121,  29.108]
       CH4_abs = CH4_ab*1.e-18

  __WRT_DEBUG_OUT("section_efficace_CH4_abs")
  end subroutine section_efficace_CH4_abs
!------------------- end SECTION_EFFICACE_CH4_abs ------------------------


!---------------------------------------------------------------------------
! Definition of O2 cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_O2_ion(nb_lo,O2_ion)
!--Definition des valeurs des sections efficaces pour l'eau
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: O2_ion
   real(dp),dimension(nb_lo) :: O2_io
   __WRT_DEBUG_IN("section_efficace_O2_ion")

O2_io = [   1.316 ,  2.346 , 4.139 , 6.619 , 8.460 , 9.890 , 9.056 , &
         & 10.860 , 10.880 ,12.229 ,13.760 ,13.418 ,15.490 ,16.970 , &
         & 17.754 , 19.469 ,21.600 ,18.840 ,22.789 ,24.540 ,30.070 , &
         & 23.974 , 21.116 ,23.750 ,23.805 ,11.720 , 8.470 ,10.191 , &
         & 10.597 ,  6.413 , 5.494 , 9.374 ,15.540 ,13.940 , 1.050 , &
            0.0   ,  0.259   ]
      O2_ion = O2_io*1.e-18



  __WRT_DEBUG_OUT("section_efficace_O2_ion")
  end subroutine section_efficace_O2_ion
!------------------- end SECTION_EFFICACE_O2_ion ------------------------  


!---------------------------------------------------------------------------
! Definition of H cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_H_ion(nb_lo,H_ion)
!--Definition des valeurs des sections efficaces pour l'Hydrogen
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: H_ion
   real(dp),dimension(nb_lo) :: H_io
   __WRT_DEBUG_IN("section_efficace_H_ion")
!--Section efficace pour l'ionisation in (cm^2)
!H_io = &
!       [  0.0011, 0.004,  0.0075, 0.0305, 0.0527, &
!       &  0.0773, 0.0661, 0.1005, 0.1011, 0.12  , &
!       &  0.1577, 0.1594, 0.1255, 0.0925, 0.0944, &
!       &  0.1020, 0.1184, 0.1208, 0.1237, 0.1429, &
!       &  0.1573, 0.1524, 0.0287, 0.0, 0.0, 0.0 , &
!       &  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
!       &  0.0, 0.0, 0.0]
H_io = &
       [  0.171 , 0.404 , 0.781 , 1.105 , 1.116 , &
       &  1.262 , 1.206 , 1.347 , 1.347 , 1.347 , &
       &  1.327 , 1.352 , 1.354 , 1.458 , 1.529 , &
       &  1.66  , 1.811 , 1.844 , 1.795 , 1.672 , &
       &  1.282 , 1.371 , 0.174 , 0.008 , 0.002, 0.0 , &
       &  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
       &  0.0, 0.0, 0.0]

       H_ion = H_io*1.e-18


  __WRT_DEBUG_OUT("section_efficace_H_ion")
  end subroutine section_efficace_H_ion
!------------------- end SECTION_EFFICACE_H_ion ------------------------
  
!---------------------------------------------------------------------------
! Definition of H2 cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_H2_ion(nb_lo,H2_ion)
!--Definition des valeurs des sections efficaces pour l'Hydrogen
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: H2_ion
   real(dp),dimension(nb_lo) :: H2_io
   __WRT_DEBUG_IN("section_efficace_H2_ion")
!--Section efficace pour l'ionisation in (cm^2)
!H_io = &
!       [  0.0011, 0.004,  0.0075, 0.0305, 0.0527, &
!       &  0.0773, 0.0661, 0.1005, 0.1011, 0.12  , &
!       &  0.1577, 0.1594, 0.1255, 0.0925, 0.0944, &
!       &  0.1020, 0.1184, 0.1208, 0.1237, 0.1429, &
!       &  0.1573, 0.1524, 0.0287, 0.0, 0.0, 0.0 , &
!       &  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
!       &  0.0, 0.0, 0.0]
H2_io = [0.0108 , 0.0798 , 0.2084, 0.4333 , 0.6036 , 0.8227, &
       & 0.7199 , 1.0004 , 1.018 , 1.0052 , 1.416  , 1.9417, & 
       & 1.9014 , 3.0155 , 3.8705, 4.1414 , 5.356  , 6.1684, &
       & 7.028  , 6.8647 , 7.8109, 8.4563 , 8.4404 , 9.7307, &
       &10.731  , 9.761  , 8.6240, 7.071  , 5.072  , 6.629 , &
       & 0.0889 , 0.0    , 0.0   , 0.0    , 0.0    ,0.0    , &
       & 0.0    ]
      H2_ion = H2_io*1.e-18



  __WRT_DEBUG_OUT("section_efficace_H2_ion")
  end subroutine section_efficace_H2_ion
!------------------- end SECTION_EFFICACE_H2_ion ------------------------  

!---------------------------------------------------------------------------
! Definition of H2O cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_H2O_ion(nb_lo,H2O_ion)
!--Definition des valeurs des sections efficaces pour l'eau
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: H2O_ion
   real(dp),dimension(nb_lo) :: H2O_io
   __WRT_DEBUG_IN("section_efficace_H2O_ion")

H2O_io = [ 0.385 , 1.153 , 2.366 , 3.595 , 4.563 , 5.552 , 4.974 , &
         & 6.182 , 6.198 , 7.237 , 8.441 , 8.218 ,10.561 ,11.908 , &
         &12.356 ,12.99  ,13.559 ,13.968 ,13.972 ,14.392 ,14.464 , &
         &14.558 ,17.443 ,18.283 ,17.557 ,13.08  ,13.512 ,10.636 , &
         &11.625 , 9.654 , 9.567 , 8.736 , 6.188 , 4.234 , 0.0   , &
           0.0   , 0.0   ]
      H2O_ion = H2O_io*1.e-18



  __WRT_DEBUG_OUT("section_efficace_H2O_ion")
  end subroutine section_efficace_H2O_ion
!------------------- end SECTION_EFFICACE_H2O_ion ------------------------  

!---------------------------------------------------------------------------
! Definition of H_H2O cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_H_H2O_ion(nb_lo,H_H2O_ion)
!--Definition des valeurs des sections efficaces pour l'eau
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: H_H2O_ion
   real(dp),dimension(nb_lo) :: H_H2O_io
   __WRT_DEBUG_IN("section_efficace_H_H2O_ion")

H_H2O_io = [ 0.171, 0.404, 0.781, 1.105, 1.166, 1.262, 1.206, 1.347, &
           & 1.347, 1.347, &
           & 1.327, 1.352, 1.354, 1.458, 1.529, 1.66 , 1.811, 1.844, &
           & 1.795, 1.672, 1.282, 1.371, 0.174, 0.008, 0.002, 0.0  , &
           & 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           & 0.0  , 0.0  , 0.0]

H_H2O_ion = H_H2O_io*1.e-18



  __WRT_DEBUG_OUT("section_efficace_H_H2O_ion")
  end subroutine section_efficace_H_H2O_ion
!------------------- end SECTION_EFFICACE_H_H2O_ion ------------------------  

!---------------------------------------------------------------------------
! Definition of O_H2O cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_O_H2O_ion(nb_lo,O_H2O_ion)
!--Definition des valeurs des sections efficaces pour l'eau
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: O_H2O_ion
   real(dp),dimension(nb_lo) :: O_H2O_io
   __WRT_DEBUG_IN("section_efficace_O_H2O_ion")

O_H2O_io = [ 0.05 , 0.107, 0.189, 0.223, 0.23 , 0.23 , 0.23 , 0.231, &
           & 0.23 , 0.207, 0.171, 0.18 , 0.075, 0.028, 0.031, 0.024, &
           & 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           & 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           & 0.0  , 0.0  , 0.0  , 0.0  , 0.0]


O_H2O_ion = O_H2O_io*1.e-18



  __WRT_DEBUG_OUT("section_efficace_O_H2O_ion")
  end subroutine section_efficace_O_H2O_ion
!------------------- end SECTION_EFFICACE_O_H2O_ion ------------------------  

!---------------------------------------------------------------------------
! Definition of OH_H2O cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_OH_H2O_ion(nb_lo,OH_H2O_ion)
!--Definition des valeurs des sections efficaces pour l'eau
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: OH_H2O_ion
   real(dp),dimension(nb_lo) :: OH_H2O_io
   __WRT_DEBUG_IN("section_efficace_OH_H2O_ion")

OH_H2O_io = [0.093, 0.306, 0.733, 1.197, 1.56 , 1.889, 1.704, 2.148, &
           & 2.154, 2.559, 3.065, 2.984, 4.042, 4.688, 4.981, 5.374, &
           & 5.789, 6.096, 6.09 , 6.383, 6.279, 6.368, 3.118, 1.364, &
           & 0.386, 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           & 0.0  , 0.0  , 0.0  , 0.0  , 0.0]


OH_H2O_ion = OH_H2O_io*1.e-18



  __WRT_DEBUG_OUT("section_efficace_OH_H2O_ion")
  end subroutine section_efficace_OH_H2O_ion
!------------------- end SECTION_EFFICACE_OH_H2O_ion ------------------------  



!---------------------------------------------------------------------------
! Definition of O cross section  for absorbtion
!---------------------------------------------------------------------------
 subroutine section_efficace_O_abs(nb_lo,O_abs)
!--Definition des valeurs des sections efficaces pour l'Hydrogen
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: O_abs
   real(dp),dimension(nb_lo) :: O_ab
   __WRT_DEBUG_IN("section_efficace_O_abs")
!--Section efficace pour l'absorption in (cm^2)
  O_ab = &
       [  0.730 ,  1.839  ,  3.732 ,  5.202 ,  6.05  ,&
       &  7.08  ,  6.461  ,  7.68  ,  7.70  ,  8.693 ,&
       &  9.84  ,  9.687  , 11.469 , 11.93  , 12.127 ,&
       & 12.059 , 12.59   , 13.09  , 13.024 , 13.4   ,&
       & 13.40  , 13.365  , 17.245 , 11.46  , 10.736 ,&
       &  4.00  ,  3.89   ,  3.749 ,  5.091 ,  3.498 ,&
       &  4.554 ,  1.1315 ,  0.0   ,  0.0   ,  0.000 ,&
       &  0.000 ,  0.000  ]
       O_abs = O_ab*1.e-18

  __WRT_DEBUG_OUT("section_efficace_O_abs")
  end subroutine section_efficace_O_abs
!------------------- end SECTION_EFFICACE_O_abs ------------------------




!---------------------------------------------------------------------------
! Definition of O cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_O_ion(nb_lo,O_ion)
!--Definition des valeurs des sections efficaces pour l'Hydrogen
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: O_ion
   real(dp),dimension(nb_lo) :: O_io
   __WRT_DEBUG_IN("sections_efficaces_O_ion")
  !--Section efficace pour l'ionisation in (cm^2)
O_io = &
       [  0.729, 1.839 , 3.732, 5.201, 6.05 , &
       &  7.080, 6.461 , 7.680, 7.70 , 8.693, &
       &  9.840, 9.707 ,11.496,11.93 ,12.127, &
       & 12.059,12.59  ,13.09 ,13.024,13.4  , &
       & 13.4  ,13.365 ,17.245,11.46 ,10.736, &
       &  4.00 , 3.89  , 3.749, 5.091, 3.498, &
       &  4.554, 1.1315, 0.000, 0.000, 0.000, &
       &  0.000,  0.000]
  O_ion = O_io*1.e-18


  __WRT_DEBUG_OUT("section_efficace_O_ion")
  end subroutine section_efficace_O_ion
!------------------- end SECTION_EFFICACE_O_ion ------------------------
  
  






!---------------------------------------------------------------------------
! Definition of CO2 cross section  for absorbtion
!---------------------------------------------------------------------------
 subroutine section_efficace_CO2_abs(nb_lo,CO2_abs)
!--Definition des valeurs des sections efficaces pour le CO2
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: CO2_abs
   real(dp),dimension(nb_lo) :: CO2_ab
   __WRT_DEBUG_IN("section_efficace_CO2_abs")
!--Section efficace pour l'absorption in (cm^2)
  CO2_ab = &
       [   1.55  ,  4.616 , 9.089  , 14.361 , 16.505 , &
       & 19.016  , 17.518 , 21.492 , 21.594 , 23.574 , &
       & 25.269  , 24.871 , 28.271 , 29.526 , 30.254 , &
       & 31.491  , 33.202 , 34.2   , 34.913 , 35.303 , &
       & 34.3    , 34.447 , 33.699 , 23.518 , 32.832 , &
       & 93.839  , 61.939 , 26.493 , 39.831 , 13.98  , &
       & 44.673  , 52.081 , 42.869 , 50.311 , 15.1   , &
       & 14.2    , 18.241 ]
       CO2_abs = CO2_ab*1.e-18


  __WRT_DEBUG_OUT("section_efficace_CO2_abs")
  end subroutine section_efficace_CO2_abs
!------------------- end SECTION_EFFICACE_CO2_abs ------------------------




!---------------------------------------------------------------------------
! Definition of CO2 cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_CO2_ion(nb_lo,CO2_ion)
!--Definition des valeurs des sections efficaces pour le CO2
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: CO2_ion
   real(dp),dimension(nb_lo) :: CO2_io
   __WRT_DEBUG_IN("section_efficace_CO2_ion")
!--Section efficace pour l'ionisation in (cm^2)
 CO2_io = &
      [  0.447 ,  2.083 ,  4.96  ,  8.515 , 11.113 &
       , 13.004 , 11.906 , 14.39  , 14.414 , 15.954 &
       , 18.271 , 17.982 , 21.082 , 25.378 , 27.163 &
       , 30.138 , 31.451 , 32.382 , 33.482 , 34.418 &
       , 33.795 , 34.003 , 32.387 , 20.856 , 27.49  &
       , 86.317 , 51.765 , 21.676 , 34.094 , 10.93  &
       ,  7.135 ,  0.000 ,  0.000 ,  0.000 ,  0.000 &
       ,  0.000 ,  0.000]
       CO2_ion = CO2_io*1.e-18


  __WRT_DEBUG_OUT("section_efficace_CO2_ion")
  end subroutine section_efficace_CO2_ion
!------------------- end SECTION_EFFICACE_CO2_ion ------------------------
  
  
  
!---------------------------------------------------------------------------
! Definition of O_CO2 cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_O_CO2_ion(nb_lo,O_CO2_ion)
!--Definition des valeurs des sections efficaces pour le O_CO2
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: O_CO2_ion
   real(dp),dimension(nb_lo) :: O_CO2_io
   __WRT_DEBUG_IN("section_efficace_O_CO2_ion")
!--Section efficace pour l'ionisation in (cm^2)
 O_CO2_io = &
       [  0.626 ,  1.32  ,  1.929 ,  2.622 ,  2.26  &
       ,  2.572 ,  2.382 ,  3.271 ,  3.28  ,  3.426 &
       ,  3.128 ,  3.224 ,  2.597 ,  2.13  ,  1.911 &
       ,  1.636 ,  1.351 ,  1.17  ,  1.171 ,  0.85  &
       ,  0.468 ,  0.527 ,  0.000 ,  0.000 ,  0.000 &
       ,  0.000 ,  0.000 ,  0.000 ,  0.000 ,  0.000 &
       ,  0.000 ,  0.000 ,  0.000 ,  0.000 ,  0.000 &
       ,  0.000 ,  0.000]
       O_CO2_ion = O_CO2_io*1.e-18


  __WRT_DEBUG_OUT("section_efficace_O_CO2_ion")
  end subroutine section_efficace_O_CO2_ion
!------------------- end SECTION_EFFICACE_O_CO2_ion ------------------------




!---------------------------------------------------------------------------
! Definition of N2 cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_N2_ion(nb_lo,N2_ion)
!--Definition des valeurs des sections efficaces pour le N2
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: N2_ion
   real(dp),dimension(nb_lo) :: N2_io
   __WRT_DEBUG_IN("section_efficace_N2_ion")
!--Section efficace pour l'ionisation in (cm^2)
 N2_io = &
       [  0.443 ,  1.479 ,  3.153 ,  5.226 ,  6.781 &
       ,  8.100 ,  7.347 ,  9.180 ,  9.210 ,  11.600&
       ,  15.350,  14.669,  20.692,  22.100,  22.772&
       ,  24.468,  24.130,  22.240,  22.787,  22.790&
       ,  23.370,  23.339,  29.235,  25.480,  15.060&
       ,  65.800,  8.500,   8.860,   14.274,  0.000 &
       ,  0.000 ,  0.000 ,  0.000 ,  0.000 ,  0.000 &
       ,  0.000 ,  0.000]
       N2_ion = N2_io*1.e-18


  __WRT_DEBUG_OUT("section_efficace_N2_ion")
  end subroutine section_efficace_N2_ion
!------------------- end SECTION_EFFICACE_N2_ion ------------------------




!---------------------------------------------------------------------------
! Definition of CH4 cross section  for ionisation
!---------------------------------------------------------------------------
 subroutine section_efficace_CH4_ion(nb_lo,CH4_ion)
!--Definition des valeurs des sections efficaces pour le CH4
!  section efficace in (cm^2)
  integer,intent(in) :: nb_lo
  real(dp),dimension(nb_lo),intent(inout) :: CH4_ion
   real(dp),dimension(nb_lo) :: CH4_io
   __WRT_DEBUG_IN("section_efficace_CH4_ion")
!--Section efficace pour l'ionisation in (cm^2)
 CH4_io = &
       [  0.103 ,  0.299 ,  0.796 ,  1.732 ,  2.482 &
       ,  3.505 ,  2.912 ,  4.382 ,  4.405 ,  6.065 &
       ,  8.277 ,  7.921 ,  13.088,  17.218,  20.049&
       ,  23.406,  27.251,  30.302,  30.014,  32.453&
       ,  34.596,  34.305,  38.527,  40.314,  42.164&
       ,  44.035,  44.079,  44.097,  44.031,  40.258&
       ,  25.527,  13.863,  0.136 ,  0.475 ,  0.000 &
       ,  0.000 ,  0.000]
       CH4_ion = CH4_io*1.e-18


  __WRT_DEBUG_OUT("section_efficace_CH4_ion")
  end subroutine section_efficace_CH4_ion
!------------------- end SECTION_EFFICACE_CH4_ion -----------------------

end module atm_sections_efficaces
