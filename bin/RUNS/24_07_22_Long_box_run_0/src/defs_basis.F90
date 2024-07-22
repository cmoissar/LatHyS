!!=============================================================
!!=============================================================
!!module: defs_basis
!! NAME
!!  defs_basis
!!
!! FUNCTION
!!  This module contains definitions for a number of named constants and
!!  physical constants.


module defs_basis

 implicit none
 public

 !Keyword 'integer' stands for default integer type
 !and may be used whenever integer are presumed to be small

 !nb of bytes related to an integer subtype n such as -10^(argument) < n < 10^(argument) (this is standard F90)
 integer, parameter :: i0b=selected_int_kind(1) !--integer 0 bytes
 integer, parameter :: i1b=selected_int_kind(2) !--integer 1 bytes
 integer, parameter :: i2b=selected_int_kind(4) !--integer 2 bytes
 integer, parameter :: i4b=selected_int_kind(9) !--integer 4 bytes default
 integer, parameter :: i8b=selected_int_kind(18)!--integer 8 bytes

 !nb of bytes related to default simple-precision real/complex subtypes
 !(= 4 for many machine architectures, = 8 for e.g. Cray)
 integer, parameter :: sp=kind(1.0)          ! Single precision should not be used
 integer, parameter :: spc=kind((1.0,1.0))

 !nb of bytes related to default double-precision real/complex subtypes
 !(= 8 for many machine architectures)
#ifdef HAVE_DOUBLE_PRECISION
 integer, parameter :: dp=kind(1.0d0)
#else
 integer, parameter :: dp=kind(1.0)
#endif
 integer, parameter :: dpc=kind((1.0_dp,1.0_dp))  ! Complex should not be used presently
 ! except for use of libraries

 !--Double precision obliged (for some subroutines simple precision is
 !  not sufficient so this type is indroduced)
 integer, parameter :: dpo=kind(1.0d0)  

 !To modify sp/spc and / or dp/dpc, insert instructions such as 'dp='
 ! but do not modify the other declarations in this module

 !Default logical type
 integer, parameter :: lgt=kind(.true.)

 !The default lengths
 integer, parameter :: fnlen=264     ! maximum length of file name variables
 integer, parameter :: strlen=2000000 ! maximum length of input string

 !Some constants:
 !integer, parameter :: integer_not_used=0
 !logical, parameter :: logical_not_used=.true.

 !UNIX unit numbers : standard input, standard output, ab_out, and a number
 !for temporary access to a file.
 integer, parameter :: std_in=5,ab_in=5  ! generally, the number 5 is directly used
 integer, parameter :: std_out=6         ! generally, the number 6 is directly used
 integer, parameter :: qp_out=101

 !Real constants
 real(dp), parameter :: zero=0._dp
 real(dp), parameter :: one=1._dp
 real(dp), parameter :: two=2._dp
 real(dp), parameter :: three=3._dp
 real(dp), parameter :: four=4._dp
 real(dp), parameter :: five=5._dp
 real(dp), parameter :: six=6._dp
 real(dp), parameter :: seven=7._dp
 real(dp), parameter :: eight=8._dp
 real(dp), parameter :: nine=9._dp
 real(dp), parameter :: ten=10._dp
 real(dp), parameter :: hundred=100._dp

 !Fractionary real constants
 real(dp), parameter :: half=0.50_dp
 real(dp), parameter :: onehalf=1.50_dp
 real(dp), parameter :: third=one/three
 real(dp), parameter :: quarter=0.25_dp
 real(dp), parameter :: fifth=0.20_dp
 real(dp), parameter :: sixth=one/six
 real(dp), parameter :: seventh=one/seven
 real(dp), parameter :: eighth=0.125_dp
 real(dp), parameter :: ninth=one/nine
 real(dp), parameter :: two_thirds=two*third
 real(dp), parameter :: four_thirds=four*third
 real(dp), parameter :: five_thirds=five*third
 real(dp), parameter :: three_quarters=0.75_dp
 real(dp), parameter :: three_fifth=three/five

 !Real constants related to the golden number
 real(dp), parameter :: gold=1.618033988749894848204586834365638117720309179_dp
 real(dp), parameter :: goldenratio=two-gold

 !Real constants derived from pi
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
 real(dp), parameter :: two_pi=two*pi
 real(dp), parameter :: four_pi=four*pi
 real(dp), parameter :: piinv=one/pi
 !The following are not used
 real(dp), parameter :: rad_to_deg=180._dp/pi
 real(dp), parameter :: deg_to_rad=one/rad_to_deg
 real(dp), parameter :: half_pi=pi*half
 !real(dp), parameter :: third_pi=pi*third
 !real(dp), parameter :: quarter_pi=pi*quarter
 !real(dp), parameter :: two_thirds_pi=two_thirds*pi

 !Real precision
 real(dp), parameter :: greatest_real = huge(one)
 real(dp), parameter :: smallest_real = -greatest_real
 real(dp), parameter :: tol3= 0.001_dp
 real(dp), parameter :: tol4= 0.0001_dp
 real(dp), parameter :: tol5= 0.00001_dp
 real(dp), parameter :: tol6= 0.000001_dp
 real(dp), parameter :: tol7= 0.0000001_dp
 real(dp), parameter :: tol8= 0.00000001_dp
 real(dp), parameter :: tol9= 0.000000001_dp
 real(dp), parameter :: tol10=0.0000000001_dp
 real(dp), parameter :: tol11=0.00000000001_dp
 real(dp), parameter :: tol12=0.000000000001_dp
 real(dp), parameter :: tol13=0.0000000000001_dp
 real(dp), parameter :: tol14=0.00000000000001_dp
 real(dp), parameter :: tol15=0.000000000000001_dp
 real(dp), parameter :: tol16=0.0000000000000001_dp

 !real(dpconstants derived from sqrt(n.)
 real(dp), parameter :: sqrt2=1.4142135623730950488016887242096939_dp
 real(dp), parameter :: half_sqrt2=0.70710678118654752440084436210484697_dp
 real(dp), parameter :: sqrt3=1.7320508075688772935274463415058739_dp
 real(dp), parameter :: half_sqrt3=0.86602540378443864676372317075293693_dp
 real(dp), parameter :: sqrthalf=0.70710678118654752440084436210484697_dp

 !Conversion factors of common use, not directly related to physical quantities.
 real(dp), parameter :: b2Mb=one/1024.0_dp**2  ! conversion factor bytes --> Mbytes
 real(dp), parameter :: b2Gb=b2Mb/1000.0_dp    ! conversion factor bytes --> Gbytes

 !Real physical constants
 !Revised fundamental constants from http://physics.nist.gov/cuu/Constants/index.html
 !(from 2006 least squares adjustment)
 !real(dp), parameter :: Bohr_Ang=0.52917720859_dp    ! 1 Bohr, in Angstrom
 !real(dp), parameter :: Ha_cmm1=219474.6313705_dp  ! 1 Hartree, in cm^-1
 !real(dp), parameter :: Ha_eV=27.21138386_dp ! 1 Hartree, in eV
 !real(dp), parameter :: Ha_K=315774.65_dp ! 1Hartree, in Kelvin
 !real(dp), parameter :: Ha_THz=6579.683920722_dp ! 1 Hartree, in THz
 !real(dp), parameter :: Ha_J=4.35974394d-18    !1 Hartree, in J
 real(dpo), parameter :: e_Cb=1.602176487d-19 ! minus the electron charge, in Coulomb
 !real(dp), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
 real(dpo), parameter :: amu_pmass=1.660538782d-27!/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
 real(dpo), parameter :: kb_JK = 1.38064813d-23   ! Boltzmann constant in J/K
 real(dpo), parameter :: M_Mars = 0.64185d24       ! Mass of Mars in kg
 real(dpo), parameter :: G_grav = 6.6738480d-11   ! Gravitational constant in m3.kg-1.s-2
 !This value is 1Ha/bohr^3 in 1d9 J/m^3
 !real(dp), parameter :: HaBohr3_GPa=29421.033_dp ! 1 Ha/Bohr^3, in GPa
 !real(dp), parameter :: HaBohr3_GPa=Ha_eV/Bohr_Ang**3*e_Cb*1.0d+21 ! 1 Ha/Bohr^3, in GPa
 !real(dp), parameter :: Avogadro=6.02214179d23 ! per mole
 !This value is 1 Ohm.cm in atomic units
 !real(dp), parameter :: Ohmcm=two*pi*Ha_THz*ninth*ten
 real(dpo), parameter :: epsilon0=8.854187817d-12 ! permittivity of free space in F/m
 !real(dp), parameter :: eps0=one/(four_pi*0.0000001_dp*299792458.0_dp**2)
 !real(dp), parameter :: AmuBohr2_Cm2=e_Cb*1.0d20/(Bohr_Ang*Bohr_Ang)
 !real(dp), parameter :: InvFineStruct=137.035999679_dp  ! Inverse of fine structure constant
 real(dp), parameter :: Sp_Lt=2.99792458d8!/2.1876912633d6 ! speed of light in atomic units
 real(dp), parameter :: mu0=four_pi*0.0000001_dp!one/(eps0*Sp_Lt**two) ! permittivity of free space in F/m
 real(dp), parameter :: pmasse_e2=amu_pmass/e_Cb**2! m_proton/e_electron**2
 !real(dp), parameter :: Time_Sec=2.418884326505D-17 !  Atomic unit of time, in seconds
 real(dp),parameter :: JtoeV = 6.24150636309d18 ! 1 Joule in eV
 !Complex constants
 !complex(dpc), parameter :: czero=(0._dp,0._dp)
 !complex(dpc), parameter :: cone =(1._dp,0._dp)
 !complex(dpc) ,parameter :: j_dpc=(0._dp,1.0_dp)

 !Character constants
 character(len=1), parameter :: ch10 = char(10)

 !--Characteristic of 2D MPI-decomposition of space
 integer,parameter :: ndims = 2    !--Dimension of the Cartesian parallelisation
 integer,parameter :: nb_voisins=8 !--Number of Influent neighbours
 integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8 !--neighbours directions

 !Define fake communicator for sequential abinit
 !integer, parameter :: abinit_comm_serial = -12345
 !----------------------------------------------------------------------
end module defs_basis
!!***
