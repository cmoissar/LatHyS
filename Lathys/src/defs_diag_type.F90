!!=============================================================
!!=============================================================
!!module: defs_diag_type
!! NAME
!!  defs_diag_type (MMancini)
!!
!! FUNCTION
!! Defines type for diagnostic,
!! Initialisation and Clean
module defs_diag_type

 use defs_basis
 
 implicit none
 private

 !!=============================================================
 !!--Type for diagnostic informations. Contains all-procs mean
 !!  quantities
 type,public :: diag_type 
  integer,allocatable,dimension(:) ::  &
       &         n_part ,&
       &         n_out ,&
       &         n_out_g , &
       &         n_in 

  real(dp),allocatable,dimension(:) ::  &
       &         ener_tot , &
       &         ener_therm ,&
       &         ener_mag ,&
       &         ener_mag_moy ,&
       &         ener_kin ,&
       &         ener_kins , &
       &         ener_kin_moy ,&
       &         dn_moy 

  integer,allocatable,dimension(:,:) ::  &
       &         n_in_spec 

  real(dp),allocatable,dimension(:,:) :: &
       &         champ_bmoy , &
       &         champ_bqua , &
       &         vitesse_moy , &
       &         vitesse_qua , &
       &         sv 
 end type diag_type
 
 public::              &
      init_diag_type,  & !--Intialisation ov variables
      clean_diag_type

contains
 !!#####################################################################

 !!--Initialization of the diag_type variable
 !!--nhm=number of diagnostic to do
 subroutine init_diag_type(diag_info,nhm,ns)
  
  type(diag_type),intent(inout) :: diag_info
  integer,intent(in) :: nhm,ns


  allocate(diag_info%sv(3,nhm))
  allocate(diag_info%champ_bmoy(3,nhm))
  allocate(diag_info%champ_bqua(3,nhm))
  allocate(diag_info%vitesse_moy(3,nhm))
  allocate(diag_info%vitesse_qua(3,nhm))
  allocate(diag_info%ener_tot(nhm))
  allocate(diag_info%ener_therm(nhm))
  allocate(diag_info%ener_mag(nhm))
  allocate(diag_info%ener_mag_moy(nhm))
  allocate(diag_info%ener_kin(nhm))
  allocate(diag_info%ener_kins(nhm))
  allocate(diag_info%ener_kin_moy(nhm))
  allocate(diag_info%dn_moy(nhm))
  allocate(diag_info%n_part(nhm))
  allocate(diag_info%n_out(nhm))
  allocate(diag_info%n_out_g(nhm))
  allocate(diag_info%n_in(nhm))
  allocate(diag_info%n_in_spec(nhm,ns))

  diag_info%sv = zero          
  diag_info%champ_bmoy = zero  
  diag_info%champ_bqua = zero  
  diag_info%vitesse_moy = zero 
  diag_info%vitesse_qua = zero 
  diag_info%ener_tot = zero    
  diag_info%ener_therm = zero  
  diag_info%ener_mag = zero    
  diag_info%ener_mag_moy = zero
  diag_info%ener_kin = zero    
  diag_info%ener_kins = zero   
  diag_info%ener_kin_moy = zero
  diag_info%dn_moy = zero      
  diag_info%n_part = 0      
  diag_info%n_out = 0       
  diag_info%n_out_g = 0     
  diag_info%n_in = 0        
  diag_info%n_in_spec = 0     

 end subroutine init_diag_type

 subroutine clean_diag_type(diag_info)

  type(diag_type),intent(inout) :: diag_info
  
  deallocate(diag_info%sv);           
  deallocate(diag_info%champ_bmoy);   
  deallocate(diag_info%champ_bqua);   
  deallocate(diag_info%vitesse_moy);  
  deallocate(diag_info%vitesse_qua);  
  deallocate(diag_info%ener_tot);     
  deallocate(diag_info%ener_therm);   
  deallocate(diag_info%ener_mag);     
  deallocate(diag_info%ener_mag_moy); 
  deallocate(diag_info%ener_kin);     
  deallocate(diag_info%ener_kins);    
  deallocate(diag_info%ener_kin_moy); 
  deallocate(diag_info%dn_moy);       
  deallocate(diag_info%n_part);       
  deallocate(diag_info%n_out);        
  deallocate(diag_info%n_out_g);      
  deallocate(diag_info%n_in);               
  deallocate(diag_info%n_in_spec);      

 end subroutine clean_diag_type

end module defs_diag_type
