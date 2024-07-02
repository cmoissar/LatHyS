module m_distribution_function
use defs_basis
use defs_parametre
use defs_species
use defs_variable
use m_writeout
 use defs_mpitype
 use defs_grid
 use m_writeout
 use m_timing,only         : time_get 
#ifdef HAVE_NETCDF
 use defs_basic_cdf
 use diag_wrt_common_cdf
 use netcdf
#endif 

#include "q-p_common.h"

implicit none

private :: &
        determine_nb_species_distrib,&! determine the number of species for distribution function (can be a subset of all ions)
        determine_species_distrib       ! determine the 
        

public :: &

        initialisation_distribution_function    !Initialization of distribution function
        
 type,public :: Species_characteristic
    character(len=10) :: Ion_label
    real(dp)          :: qsm_value
    integer       :: origin_value
 !   integer       :: CE_value
 end type Species_characteristic 
 
 real(dp) :: Ener_min_distrib = 10.,dE_E_distrib = 0.3 ! First level of energy (energy max) and energy resolution


contains

!!---------------------------------
!! initialisation of distribution function arrays
!! RModolo, May 2014
subroutine initialisation_distribution_function
   real(dp) :: mem_bit_proc
   integer(1) :: sizeof(9999),i
   integer :: nb_species,tot_size_arr
  character(len=3) :: unit
  character(len=500) :: msg   
   character(len=10),dimension(:),allocatable :: Ion_label_tab
   type(Species_characteristic),dimension(:), allocatable :: species_info
   
 
   __WRT_DEBUG_IN("initialisation_distribution_function")
   __GETTIME(9,1) !--Timer start

   !--Determine the number of ion species
   call determine_nb_species_distrib(planetname,nb_species)
   
   if (mpiinfo%me == 0) print *,'Nb ion species for ',trim(planetname),' is ',nb_species
   


  tot_size_arr = ncm(1)*ncm(2)*ncm(3)*n_Ene_distrib*n_theta_distrib*n_phi_distrib/8._dp
  mem_bit_proc = real(tot_size_arr,dp)*real(size(transfer(1.0_dp,sizeof)),dp)
  !--Conversion in GB or MB
  mem_bit_proc = b2Gb*mem_bit_proc*real(nproc,dp)
  unit = ' GB'
  if(mem_bit_proc < one) then
   mem_bit_proc = mem_bit_proc*1000.0_dp
   unit = ' MB'   
  endif
  
  write(msg,'(2a,2(a,a17,f8.3,a))')&
       & ch10," ____________ Memory Required for distribution function _________________",&
       & ch10, "   Total Memory  = ",mem_bit_proc,unit,&
       & ch10, "   Mem. per proc = ",mem_bit_proc/real(nproc,dp),unit
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)  
   
   allocate(species_info(nb_species))
      
  ! determine the charateristic of the different ion species
  call determine_species_distrib(planetname,nb_species,species_info)

! correct a bug with some netCDF libraries 
  !allocate(Ion_label_tab(nb_species))
  !do jj=1,nb_species
!	Ion_label_tab(jj)=species_info(jj)%Ion_label
  !enddo
  
  !-- record the restart iteration
  iter_start_distrib = iter
  
 !-- allocation of distribution function arrays
 !-- The first 3 columns indicate the position x,y,z
 !-- the 4th column is the Energy
 !-- the 5th column indicate the theta angle (V_theta/Vtot) which is the co-latitude of the velocity vector
 !--    theta angle is defined as follow : theta = acos(Vz/sqrt(Vx^2+Vy^2+Vz^2))
 !--	theta = 90 means velocity vector in the equatorial plane, theta = 0° means North pole
 !--	theta = 180° means South pole
 !-- the 6th column indicate the phi angle (V_phi/Vtot) which is the longitude of the velocity vector
 !--	phi angle is defined as follow : phi = atan(Vy/Vx)
 !--	phi=0° means velcoty vector aligned with -X_simu axis
 !-- 	phi=90° means velocity vector aligned with -Y_simu axis
 
 allocate(distrib_func(int(ncm(1)/2)+1,int(ncm(2)/2)+1,int(ncm(3)/2)+1,n_Ene_distrib,n_theta_distrib,n_phi_distrib))
 distrib_func(:,:,:,:,:,:) = zero
 
 !-- allocation of Energy , Theta and Phi array
 allocate(Energy_distrib(n_Ene_distrib));       Energy_distrib(:) = zero
 allocate(Theta_distrib(n_theta_distrib));      Theta_distrib(:) = zero
 allocate(Phi_distrib(n_phi_distrib));          Phi_distrib(:) = zero
 
 !-- Initialisation of Energy, Theta and Phi
 do i=1,n_phi_distrib
        Phi_distrib(i) = (i - 1._dp) * 2._dp * pi / (n_phi_distrib-1)
 enddo
  
 do i=1,n_theta_distrib
        Theta_distrib(i) = (i - 1._dp) * pi / (n_theta_distrib-1)
 enddo
  
     !--- Energy intervals
     Energy_distrib(1) = Ener_min_distrib
     do i = 2 , n_Ene_distrib
        Energy_distrib(i) = Energy_distrib(i-1)*(1 + dE_E_distrib)
     enddo
 
  
   __GETTIME(9,2) !--Timer stop
   __WRT_DEBUG_OUT("initialisation_distribution_function")
end subroutine initialisation_distribution_function

!!================================================================
!! subroutine : compute_distribution_function
!! FUNCTION
!! 	from the particle information (position, velocity, mass, charge)
!!	determine the contribution of this particle to the total distribution function
!! INPUTS : position of the particle (3), velocity of the particle (3), mass (1), charge (1)
!! OUTPUTS : none (add contribution to the distribution function array)

subroutine compute_distribution_function(xp,vp,sm,sq,Spe)
real(dp),intent(in) :: xp(3),vp(3),sq,sm
type(species_type),intent(in) :: Spe

real(dp) :: vtot,E_part,phi0,theta0,dvphi,dvtheta
integer  :: ijk(3),i_ener,iV_phi,iV_theta
real(dp) :: w(8),s_b(3),s_a(3),s_f(3)
real(dp) :: weight,gstep_inv(3)
 character(len=500) :: msg

   __WRT_DEBUG_IN("compute_distribution_function")
   __GETTIME(9,1) !--Timer start

dvphi = two_pi/n_phi_distrib
dvtheta = pi/n_theta_distrib
vtot = sqrt(dot_product(vp,vp))*Spe%ref%alfvenspeed*1.e3 !(in m/s)
E_part = 0.5*amu_pmass*(sm/sq)*vtot**2*JtoeV    ! energy of the particule in eV
!if (abs(sm/sq-2._dp) < 0.1_dp) E_part = E_part*2. ! to take into account He++ SW 
!if (mpiinfo%me == 0) print *,sm,sq,E_part
weight = sq
if (abs(sm/sq-2._dp) < tol3) weight = weight/2.

 !   write(msg,'(2a,f9.7)')ch10,&
 !      " dvphi = ",dvphi
 !  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 !   write(msg,'(2a,f9.7)')ch10,&
 !      " dvtheta = ",dvtheta
 !  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)   
 !      write(msg,'(2a,f9.2)')ch10,&
 !         " Ener = ",E_part
 !  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

!-- Determination of the Energy index

i_ener = int(log(E_part/Ener_min_distrib)/log(1.+dE_E_distrib))+2 ! it corresponds to the upper energy value
if (i_ener <1) i_ener=1
   if (i_ener <= n_ene_distrib) then
! if (i_ener <1 .or. i_ener > n_ene_distrib) print *, 'Pb energy :',E_part,i_ener,Ener_min_distrib,dE_E_distrib

!-- Determinition of the index position and relative weight with neighbours
   gstep_inv = one/(2.*gstep)
   s_b = (1._dp + ((xp-s_min_loc)*gstep_inv))
   ijk = int(s_b)
   if (any(ijk(:) < 1).or.any(ijk(:) > int(ncm(:)/2.))) print *,'Pb position :',xp,s_min_loc,gstep_inv
   
   s_f = s_b-real(ijk,dp)

   !--Sequence of indices at cell corners
   !-- (i,j,k)
   !--Trilinear weights
   s_a = 1._dp - s_f
   
   w(1) = s_a(1)*s_a(2)*s_a(3) !xa*ya*za
   w(2) = s_f(1)*s_a(2)*s_a(3) !xf*ya*za
   w(3) = s_a(1)*s_f(2)*s_a(3) !xa*yf*za
   w(4) = s_f(1)*s_f(2)*s_a(3) !xf*yf*za
   w(5) = s_a(1)*s_a(2)*s_f(3)
   w(6) = s_f(1)*s_a(2)*s_f(3)
   w(7) = s_a(1)*s_f(2)*s_f(3)
   w(8) = s_f(1)*s_f(2)*s_f(3)
   
   !-- Be CAREFULL WORK ONLY FOR MSO coordinate
   phi0 = atan(vp(2)/vp(1))
   if (-vp(1) < 0._dp) then
     phi0 = phi0+pi*0.99999_dp
   else
     if ((-vp(2)<0._dp).and.(vp(1).ne.0._dp)) phi0 = phi0+two_pi*0.99999_dp
   endif
   iV_phi = int(phi0/dvphi)+1
   if ((iV_phi <1).or.(iV_phi > n_phi_distrib)) print *,'Pb phi',vp,phi0,dvphi,n_phi_distrib
   
   theta0 = acos(vp(3)/sqrt(dot_product(vp,vp)))
   if (theta0 >= pi) theta0 = pi*0.99999
   if (theta0 <0._dp) theta0 = tol5
   iV_theta = int(theta0/dvtheta)+1
   
   if ((iV_theta <1).or.(iV_theta > n_theta_distrib)) print *,'Pb theta',vp,theta0,dvtheta,n_theta_distrib
   
   
   distrib_func(ijk(1)  ,ijk(2)  ,ijk(3)  , i_ener,iV_theta,iV_phi) = &
     distrib_func(ijk(1)  ,ijk(2)  ,ijk(3)  , i_ener,iV_theta,iV_phi) + w(1)*weight
   distrib_func(ijk(1)+1,ijk(2)  ,ijk(3)  , i_ener,iV_theta,iV_phi) = &
     distrib_func(ijk(1)+1,ijk(2)  ,ijk(3)  , i_ener,iV_theta,iV_phi) + w(2)*weight     
   distrib_func(ijk(1)  ,ijk(2)+1,ijk(3)  , i_ener,iV_theta,iV_phi) = &
     distrib_func(ijk(1)  ,ijk(2)+1,ijk(3)  , i_ener,iV_theta,iV_phi) + w(3)*weight   
   distrib_func(ijk(1)+1,ijk(2)+1,ijk(3)  , i_ener,iV_theta,iV_phi) = &
     distrib_func(ijk(1)+1,ijk(2)+1,ijk(3)  , i_ener,iV_theta,iV_phi) + w(4)*weight 
   distrib_func(ijk(1)  ,ijk(2)  ,ijk(3)+1, i_ener,iV_theta,iV_phi) = &
     distrib_func(ijk(1)  ,ijk(2)  ,ijk(3)+1, i_ener,iV_theta,iV_phi) + w(5)*weight
   distrib_func(ijk(1)+1,ijk(2)  ,ijk(3)+1, i_ener,iV_theta,iV_phi) = &
     distrib_func(ijk(1)+1,ijk(2)  ,ijk(3)+1, i_ener,iV_theta,iV_phi) + w(6)*weight     
   distrib_func(ijk(1)  ,ijk(2)+1,ijk(3)+1, i_ener,iV_theta,iV_phi) = &
     distrib_func(ijk(1)  ,ijk(2)+1,ijk(3)+1, i_ener,iV_theta,iV_phi) + w(7)*weight   
   distrib_func(ijk(1)+1,ijk(2)+1,ijk(3)+1, i_ener,iV_theta,iV_phi) = &
     distrib_func(ijk(1)+1,ijk(2)+1,ijk(3)+1, i_ener,iV_theta,iV_phi) + w(8)*weight 
     
     endif
     
   __GETTIME(9,2) !--Timer stop
   __WRT_DEBUG_OUT("compute_distribution_function")     

end subroutine compute_distribution_function


!!================================================================
!! subroutine : wrt_distribution_function
!!  FUNCTION
!! write the distribution function and information concerning it in
!! a netcdf file
subroutine wrt_distribution_function(filwrt)
 character(len=*),intent(in) :: filwrt
 
   integer :: ncid, stId,ii,jj
   integer :: dimid(11), varid(100)
   integer :: dim_distrib(9)
   real(dp) :: integration_time,fac,vol,surf
   character(len=40) :: name_file 
   integer,allocatable :: dimspec(:)   
   
!-- Conversion in physical unit (ions.m-2.s-1.st-1.eV-1)
integration_time = (iter-iter_start_distrib)*dt*Spe%ref%inv_gyro
vol = gstep(1)*gstep(2)*gstep(3)*(Spe%ref%c_omegapi)**3
surf = 6*(2*gstep(1))*(2*gstep(2))*(Spe%ref%c_omegapi)**2
fac = vol/surf/integration_time/Energy_distrib(n_ene_distrib)*Spe%ref%density
   
 
   __WRT_DEBUG_IN("wrt_distribution_function")
   __GETTIME(9,1) !--Timer start
   !--Create the file name
   call create_file_name(name_file,filwrt(1:8),mpiinfo%me) 

   stId = nf90_create(name_file , nf90_clobber, ncid)
   call test_cdf(stId)
 
    !--Define Speces dimensions
   call species_def_dim_cdf(ncid,dimspec)
 
   !--Define Ion species dimension
  stId = nf90_def_dim(ncid, "dim_scalar"      , 1         , dim_distrib(1))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "bidimensional"   , 2         , dim_distrib(2))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "space_dimension" , 3         , dim_distrib(3))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_x"          , int(ncm(1)/2)+1    , dim_distrib(4))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_y"          , int(ncm(2)/2)+1    , dim_distrib(5))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_z"          , int(ncm(3)/2)+1    , dim_distrib(6))
  call test_cdf(stId)
   stId = nf90_def_dim(ncid, "size_Energy"          , n_ene_distrib    , dim_distrib(7))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_Theta"          , n_theta_distrib    , dim_distrib(8))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_Phi"          , n_phi_distrib    , dim_distrib(9))
  call test_cdf(stId) 
   
   
   
   !--Define Speces dimensions
!   call species_def_dim_cdf(ncid,dimspec)

   ii = 1
      !--Define Speces variables
   call species_def_var_cdf(ncid,varid,dimspec,ii)
!   call common_def_var_cdf(ncid,varid,dimid,ii) 
   stId = nf90_def_var(ncid, "ncm",    nf90_int,dim_distrib(3),varid(ii))
   call test_cdf(stId); ii = ii+1 
   stId = nf90_def_var(ncid, "gstep",    QP_NF90_DP,dim_distrib(3),varid(ii))
   call test_cdf(stId); ii = ii+1 
   stId = nf90_def_var(ncid, "s_min_loc",  QP_NF90_DP,dim_distrib(3),  varid(ii))
   call test_cdf(stId); ii = ii+1 
   stId = nf90_def_var(ncid, "s_max_loc",  QP_NF90_DP,dim_distrib(3),  varid(ii))
   call test_cdf(stId); ii = ii+1    
   !--number of proc
   stId = nf90_def_var(ncid, "nproc", nf90_int,dim_distrib(1),varid(ii))
   call test_cdf(stId); ii = ii+1  
   !-- number of energy level
   stId = nf90_def_var(ncid, "nEnergy", nf90_int,dim_distrib(1),varid(ii))
   call test_cdf(stId); ii = ii+1  
   !-- number of Theta values
   stId = nf90_def_var(ncid, "nTheta", nf90_int,dim_distrib(1),varid(ii))
   call test_cdf(stId); ii = ii+1   
   !-- number of phi values
   stId = nf90_def_var(ncid, "nPhi", nf90_int,dim_distrib(1),varid(ii))
   call test_cdf(stId); ii = ii+1   
   !-- Energy range
   stId = nf90_def_var(ncid, "Energy", QP_NF90_DP,dim_distrib(7),varid(ii))
   call test_cdf(stId); ii = ii+1  
   !-- Theta Range
   stId = nf90_def_var(ncid, "Theta", QP_NF90_DP,dim_distrib(8),varid(ii))
   call test_cdf(stId); ii = ii+1   
   !-- Phi Range
   stId = nf90_def_var(ncid, "Phi", QP_NF90_DP,dim_distrib(9),varid(ii))
   call test_cdf(stId); ii = ii+1    
   !-- Distribution function array
   stId = nf90_def_var(ncid, "Distribution_function", QP_NF90_DP,dim_distrib(4:9),varid(ii))
   call test_cdf(stId); ii = ii+1    
   
   
 !--Switch to write mode
   stId = nf90_enddef(ncid); call test_cdf(stId)
 
   ii = 1   
   !--Write common variables into the file
   !--Write Species infos into the file
   call species_put_var_cdf(Spe,ncid,varid,dimspec,ii)
  
   !--Write the other variables into the file
   stId = nf90_put_var(ncid, varid(ii), ncm)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), gstep)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), s_min_loc)
   call test_cdf(stId); ii = ii+1              
   stId = nf90_put_var(ncid, varid(ii), s_max_loc)
   call test_cdf(stId); ii = ii+1                           
   stId = nf90_put_var(ncid, varid(ii), mpiinfo%nproc)
   call test_cdf(stId); ii = ii+1   
   stId = nf90_put_var(ncid, varid(ii), n_ene_distrib)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), n_theta_distrib)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), n_phi_distrib)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), Energy_distrib)
   call test_cdf(stId); ii = ii+1   
   stId = nf90_put_var(ncid, varid(ii), Theta_distrib)
   call test_cdf(stId); ii = ii+1     
   stId = nf90_put_var(ncid, varid(ii), Phi_distrib)
   call test_cdf(stId); ii = ii+1   
   stId = nf90_put_var(ncid, varid(ii), distrib_func*fac)
   call test_cdf(stId); ii = ii+1    
   
   !--Close the file
   stId = nf90_close(ncid); call test_cdf(stId)   
   
   __GETTIME(9,2) !--Timer stop
   __WRT_DEBUG_OUT("wrt_distribution_function")
end subroutine wrt_distribution_function

 !!=============================================================
 !!subroutine: diag_particles/create_file_name
 !! FUNCTION 
 !!  Create the name of the file containing particles information
 !! INPUT
 !!  filwrt=suffix containgt data
 !!  me=processus identity
 !! OUTPUT
 !!  name_file=name of the file where particles will be recorded
 subroutine create_file_name(name_file,filwrt,me)

  integer,intent(in) :: me
  character(len=40),intent(out) :: name_file 
  character(len=*),intent(in) :: filwrt 
  
  write(name_file,'(a8,i4.4,a1,a)')"Distrib_",me,'_',trim(filwrt)
#ifdef HAVE_NETCDF
  name_file = trim(name_file)//".nc"
#endif

 end subroutine create_file_name


 !!==============================================================
 !! subroutine: m_split_cdf/determine_nb_species
 !!  FUNCTION
 !! returns the number of different ion species present in the simulation
 !! and create a characteris profile for each ion species
 !! INPUT 
 !!  Spe  =  information concerning the simulation
 !! OUTPUT 
 !!  nb_species : number of Ion species present
 !! species_info = charactersitic information for each ion species
    
   subroutine determine_nb_species_distrib(planetname,nb_species)
   
   character(len=*),intent(in) :: planetname
   integer,intent(inout)           :: nb_species  
   character(len=500) :: msg,planet
 __WRT_DEBUG_IN("determine_nb_species_distrib")
 
   
   select case(trim(planetname))
   case("mars")
      nb_species = 6
   case("venus")
      nb_species = 6
   case("titan")
      nb_species = 3
   case("mercure")
      nb_species = 2
   case("ganymede")
      nb_species = 5
   case("earth")
      nb_species=1
   case default
      write(*,*) &
          "ERROR: Selected Species determination:",&
          "Planet '",trim(planetname),"' does Not Exist"
      stop
    end select    
      
   __WRT_DEBUG_OUT("determine_nb_species_distrib")
 end subroutine determine_nb_species_distrib
 
 !!########################################################
  !! subroutine: diag_moment_species/determine_species
  !!  FUNCTION
  !! returns the number of different ion species present in the simulation
  !! and create a characteris profile for each ion species
  !! INPUT 
  !!  Spe  =  information concerning the simulation
  !! OUTPUT 
  !!  nb_species : number of Ion species present
  !! species_info = charactersitic information for each ion species
   
  subroutine determine_species_distrib(planetname,nb_species,species_info)
  
  character(len=*),intent(in)                              :: planetname
  integer,intent(in)                                       :: nb_species
  type(Species_characteristic), dimension(nb_species),intent(inout) :: species_info
  character(len=500) :: msg
 __WRT_DEBUG_IN("determine_species_distrib")
  
  select case(trim(planetname))
   case("mars")
     !  H+ sw
     species_info(1)%Ion_label = "Hsw"
     species_info(1)%qsm_value = 1._dp
     species_info(1)%origin_value = 0
 !    species_info(1)%CE_value = 0
     
     !He++ sw
     species_info(2)%Ion_label = "Hesw"
     species_info(2)%qsm_value = 2._dp/4._dp
     species_info(2)%origin_value = 0
 !    species_info(2)%CE_value = 0
     
     !O+ planetary
     species_info(3)%Ion_label = "Opl"
     species_info(3)%qsm_value = 1._dp/16._dp
     species_info(3)%origin_value = 1
 !    species_info(3)%CE_value = 0
     
     !O2+ planetary
     species_info(4)%Ion_label = "O2pl"
     species_info(4)%qsm_value = 1._dp/32._dp
     species_info(4)%origin_value = 1
 !    species_info(4)%CE_value = 0
 
     !CO2+ planetary
     species_info(5)%Ion_label = "CO2pl"
     species_info(5)%qsm_value = 1._dp/44._dp
     species_info(5)%origin_value = 1
 !    species_info(5)%CE_value = 0
 
     !H+ planetary
     species_info(6)%Ion_label = "Hpl"
     species_info(6)%qsm_value = 1._dp
     species_info(6)%origin_value = 1
     
    case("mercure")

     !  H+ sw
     species_info(1)%Ion_label = "Hsw"
     species_info(1)%qsm_value = 1._dp/1._dp
     species_info(1)%origin_value = 0
 !    species_info(1)%CE_value = 0
     
     !He++ sw
     species_info(2)%Ion_label = "Hesw"
     species_info(2)%qsm_value = 2._dp/4._dp
     species_info(2)%origin_value = 0
 !    species_info(2)%CE_value = 0
     
 
   case("ganymede")
     !  O+ jv
     species_info(1)%Ion_label = "Ojv"
     species_info(1)%qsm_value = 1._dp/1._dp
     species_info(1)%origin_value = 0
 !    species_info(1)%CE_value = 0
     
     !  H+ jv
     species_info(2)%Ion_label = "Hjv"
     species_info(2)%qsm_value = 16._dp/1._dp
     species_info(2)%origin_value = 0
 !    species_info(2)%CE_value = 0
     
     !  H+ pl
     species_info(3)%Ion_label = "Hpl"
     species_info(3)%qsm_value = 16._dp/1._dp
     species_info(3)%origin_value = 1
 !    species_info(3)%CE_value = 0
     
     !  O+ pl
     species_info(4)%Ion_label = "Opl"
     species_info(4)%qsm_value = 1._dp/1._dp
     species_info(4)%origin_value = 1
 !    species_info(4)%CE_value = 0     
 
      !  H2+ pl
      species_info(5)%Ion_label = "H2pl"
      species_info(5)%qsm_value = 16._dp/2._dp
      species_info(5)%origin_value = 1
 !    species_info(5)%CE_value = 0   

 
    case("earth")
      !  H+ sw
      species_info(1)%Ion_label = "Hsw"
      species_info(1)%qsm_value = 1._dp
      species_info(1)%origin_value = 0
 !    species_info(1)%CE_value = 0
    case default
    write(msg,'(a,3x,a,4x,3a)')ch10,&
         "ERROR: Selected Species determination:",&
         "Planet '",trim(Spe%planetname),"' does Not Exist"
    stop
   end select    
 __WRT_DEBUG_OUT("determine_species_distrib")
  end subroutine determine_species_distrib
end module m_distribution_function
