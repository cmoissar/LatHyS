!!=============================================================
!!=============================================================
module defs_parametre

 use defs_basis

 implicit none
 ! Ce module permet l'initialisation de tous les paramètres "variables" pour des simulations

!***************** SPACE ***********************
 !--Number of cells in the three dimensions X,Y,Z (Initial try)
 integer,save :: ncx0 = 200
 integer,save :: ncy0 = 50
 integer,save :: ncz0 = 50
 !--Spatial step (/dx,dy,dz/)
 real(dp),save,dimension(3) :: gstep = one !--grid step in any direction
!***********************************************

!****************** TIME ***********************
 real(dp),save :: dt = 0.05_dp   !--Time step
 integer,save :: nhm = 500 !1000   !--Nombre de pas de temps maximum

! Subcylcings
 integer,save  :: nsub = 4        !--nombre de sous intervalle dans l'integration de B
 integer,save  :: ntest = 5       !--nombre de pas de temps entre le teste des deux version de B (B & Bh)
 real(dp),save :: eps = tol4  !--precision
!***********************************************


!**************** PARTICLES ********************
 integer,save :: n_part_max = 150 !--Max number of particles in per cell at t=0
 integer,save :: npm            !--Max number of points in the box

 integer,parameter :: nfl = 1000 !--dimension of speed distribution function
 integer,save :: nwrm = 15  !--dimension max du tableau de parametre de lissage
!***********************************************

!**************** ENVIRONMENT ******************
 !--Default planet and Number of species
 integer,save :: ns = 2         ! Is overwrite later, hence its value here has no meaning
 character(len=8),save :: planetname = "earth"
 character(len=120) :: exospherename ="" 
 character(len=120) :: atmospherename ="" 
 character(len=120) :: ionospherename=""
 real(dp) :: pos_plan(3)       !--obstacle position
!***********************************************

!**************** SETTINGS *********************
! field and ionosphere settings
 integer,save  ::   idisp  = 1         !--Dispersion or not
 integer,save  ::   iresis = 1         !--Resistivity or not
 real(dp),save ::   resis  = 5.e-02_dp !--Valeur de la resistivité
 real(dp),save ::   dmin   = 1.e-03_dp !--Densité minimale pour le calcul du champ électrique
 integer,save  ::   ipe    = 2         !--0 = adiabatic, 1 = isothermal
               !2= adiabatic + hydrostatic ionosphere, 3= isothermal + hydrostatic ionosphere
 real(dp),save ::   dn_lim_inf_conduc = 300. ! density (in sim. units) above which the
               ! conductivity is considered infinity (vi-ve=0)

 integer,save  ::   absorp_surface    = 0   ! if =1 then the particle are absorbed at the planet surface, otherwise they are lost
 real(dp),save ::   t_iono_release    = 0. ! 50 !time during which the ionosphere is frozen

! IMF direction
! the next three lines are just here for information. 
! The calculation is done in "initialisation.F90"
!   bx0  = b0*cos(rphi)*sin(rpsi)  
!   by0  = b0*sin(rphi)*sin(rpsi)
!   bz0  = b0*cos(rpsi)
! IMF direction
 real(dp),save::   psi=90_dp
 real(dp),save::   phi=-95_dp

 !Planet settings
  integer,save  :: no_env_mag=0  !if 1 then no planetary magnetic field
  real(dp),save :: planet_ssl=-1._dp !sub solar_longitude
  real(dp),save :: planet_sslat = 0._dp ! sub solar latitude
  integer,save  :: t_init_dip = 100  ! number of time steps to initialize the dipole
!***********************************************

 !--Activates restart
 ! restart = 0, non actif; 1, restart simple; 2, restart avec ajout de cellules en Y; 3, restart avec ajout de cellule en Z
 integer,save :: restart = 0

 !--Date for diagnostic (dd_mm_yy format)
 character(len=8),save :: fildat = "" 
 !--Output file name
 character(len=40),save :: qp_out_name = "sortie"
 
 integer,save :: diag_part_imp = 0 
 integer,save :: n_part_imp=0,n_part_del=0
 integer,save :: distrib_activated = 0

 !--Defines a cell where we get the magnetic field at every time_step
 real(dp) :: x_out_sim, y_out_sim, z_out_sim

!-- defines the region of interest in fraction of the total size
 type :: ROI_type
  real(dp) :: ymin=0.4
  real(dp) :: ymax=0.6
  real(dp) :: zmin=0.4
  real(dp) :: zmax=0.6
  real(dp) :: excess=2.              !-- defines by how much more the cells in the ROI weights for the grid dim
  real(dp) :: excess_part=1.  !-- defines how many times more particles are injected in the ROI
 end type ROI_type

 type(ROI_type),SAVE :: ROI

end module defs_parametre


