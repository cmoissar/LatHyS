!!=============================================================
!!=============================================================
module diag_moment_species

 use defs_basis
 use defs_mpitype
 use defs_species
 use defs_diag_type,only   : diag_type
 use defs_parametre,only   :planetname,gstep,npm
 use defs_variable
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
 private

 private ::              &
      create_file_name,  &    !--create the file name for particles record
      species_moments
 
 public  ::              &
      wrt_moment_species, &       !--ions moments for each species
      determine_nb_species, &      
      determine_species

 type,public :: Species_characteristic
    character(len=10) :: Ion_label
    real(dp)          :: qsm_value
    integer       :: origin_value
 !   integer       :: CE_value
 end type Species_characteristic      

contains
 !!############################################################

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
  
  write(name_file,'(a5,i4.4,a1,a)')"Mom3_",me,'_',trim(filwrt)
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
    
   subroutine determine_nb_species(planetname,nb_species)
   
   character(len=*),intent(in) :: planetname
   integer,intent(inout)           :: nb_species  
   character(len=500) :: msg,planet
 __WRT_DEBUG_IN("determine_nb_species")
 
   
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
      
   __WRT_DEBUG_OUT("determine_nb_species")
 end subroutine determine_nb_species
 
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
   
  subroutine determine_species(planetname,nb_species,species_info)
  
  character(len=*),intent(in)                              :: planetname
  integer,intent(in)                                       :: nb_species
  type(Species_characteristic), dimension(nb_species),intent(inout) :: species_info
  character(len=500) :: msg
 __WRT_DEBUG_IN("determine_species")
  
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

  case("venus")
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
     
     !  H+ pl
!     species_info(3)%Ion_label = "Hpl"
!     species_info(3)%qsm_value = 1._dp/1._dp
!     species_info(3)%origin_value = 1
 !    species_info(1)%CE_value = 0    
 
      !  Na+ pl
!     species_info(4)%Ion_label = "Napl"
!     species_info(4)%qsm_value = 1._dp/23._dp
!     species_info(4)%origin_value = 1
 !    species_info(1)%CE_value = 0    
 
     !  He+ pl
     !species_info(4)%Ion_label = "Hepl"
     !species_info(4)%qsm_value = 1._dp/4._dp
     !species_info(4)%origin_value = 1
 
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

    !  H2O+ pl
  !    species_info(6)%Ion_label = "H2Opl"
  !    species_info(6)%qsm_value = 16._dp/18._dp
  !    species_info(6)%origin_value = 1
 !!    species_info(6)%CE_value = 0    
 !
  !  !  OH+ pl
  !     species_info(7)%Ion_label = "OHpl"
  !     species_info(7)%qsm_value = 16._dp/17._dp
  !     species_info(7)%origin_value = 1
 !!    species_info(7)%CE_value = 0   
!
  !  !  O2+ pl
  !     species_info(8)%Ion_label = "O2pl"
  !     species_info(8)%qsm_value = 2._dp
  !     species_info(8)%origin_value = 1
 !!    species_info(8)%CE_value = 0    
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
 __WRT_DEBUG_OUT("determine_species")
  end subroutine determine_species
  
!!########################################################
!! subroutine: diag_moment_species/species_moments
!!  FUNCTION
!! returns the number of particle for a given species
!! INPUT 
!!  particle  =  particle information
!!  qsm = charge over mass criterai
!!  org = origin criteria
!! OUTPUT 
!!  density_species = density array for each ion species
 
 subroutine species_moments(particule,nptot,species_info,ncm,nc1,gstep,s_min,nb_species,density,Vx,Vy,Vz,Temp)
 
 integer,intent(in) :: nptot,nb_species
 integer,dimension(3),intent(in) :: ncm,nc1
 type(particletype),dimension(nptot),intent(in) :: particule
 real(dp),intent(in)           :: gstep(3),s_min(3)
 real(dp),dimension(ncm(1),ncm(2),ncm(3),nb_species),intent(inout) :: density,Vx,Vy,Vz,Temp
 type(species_characteristic),intent(in) :: species_info(nb_species)
 ! local variable
 integer :: nn,ijk(3),ii
 real(dp) :: sqp,smp,qsmp,vxp,vyp,vzp
 real(dp) :: w1,w2,w3,w4,w5,w6,w7,w8
 real(dp),dimension(:,:,:),allocatable :: dn_temp
 real(dp),dimension(3) :: s_f,s_a,v_p,s_m,gstep_inv
 __WRT_DEBUG_IN("species_moments")
 
 !allocate(dn_temp(ncm(1),ncm(2),ncm(3)));	
 
 gstep_inv = one/gstep
 
 do nn = 1,nptot
!   print *,'Particle n:',nn
   !dn_temp(:,:,:) = 0.
   sqp = particule(nn)%char !--On collecte la charge de la particule
   smp = particule(nn)%mass !--On collecte la masse de la particule
   qsmp = sqp/smp
   vxp = particule(nn)%vel(1) !-- On collecte la composante Vx de la vitesse
   vyp = particule(nn)%vel(2) !-- On collecte la composante Vy de la vitesse
   vzp = particule(nn)%vel(3) !-- On collecte la composante Vz de la vitesse
      
      !--Relative position in cell s_f=(xf,yf,zf) center of the particule
      s_m = one + (particule(nn)%pos-s_min)*gstep_inv
      
      ijk = int(s_m)
      
#ifdef HAVE_DEBUG   
   if(any(ijk>nc1)) then
    print *,"ERROR Moment",particule(nn)%pos
    stop
   endif
#endif   
      s_f = s_m-real(ijk,dp)
      
      !--Sequence of indices of B at cell corners
      !--Trilinear weight
      s_a = one-s_f
   
      w1 = s_a(1)*s_a(2)*s_a(3)*sqp
      w2 = s_f(1)*s_a(2)*s_a(3)*sqp
      w3 = s_a(1)*s_f(2)*s_a(3)*sqp
      w4 = s_f(1)*s_f(2)*s_a(3)*sqp
      w5 = s_a(1)*s_a(2)*s_f(3)*sqp
      w6 = s_f(1)*s_a(2)*s_f(3)*sqp
      w7 = s_a(1)*s_f(2)*s_f(3)*sqp
      w8 = s_f(1)*s_f(2)*s_f(3)*sqp
      
      !--Number density (in DN)
      !dn_temp(ijk(1)  ,ijk(2)  ,ijk(3)  ) =  w1
      !dn_temp(ijk(1)+1,ijk(2)  ,ijk(3)  ) =  w2
      !dn_temp(ijk(1)  ,ijk(2)+1,ijk(3)  ) =  w3
      !dn_temp(ijk(1)+1,ijk(2)+1,ijk(3)  ) =  w4
      !dn_temp(ijk(1)  ,ijk(2)  ,ijk(3)+1) =  w5
      !dn_temp(ijk(1)+1,ijk(2)  ,ijk(3)+1) =  w6
      !dn_temp(ijk(1)  ,ijk(2)+1,ijk(3)+1) =  w7
      !dn_temp(ijk(1)+1,ijk(2)+1,ijk(3)+1) =  w8
      
      do ii = 1,nb_species
        if (species_info(ii)%qsm_value == qsmp) then
          if (species_info(ii)%origin_value == 0) then
            if ((particule(nn)%orig == 0).and.(particule(nn)%exc == 0)) then
               density(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) = density(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii)   + w1
               density(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) = density(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii)   + w2    
               density(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) = density(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii)   + w3
               density(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) = density(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii)   + w4    
               density(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) = density(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii)   + w5
               density(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) = density(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii)   + w6    
               density(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) = density(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii)   + w7
               density(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) = density(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii)   + w8

               Vx(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) = Vx(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii)   + w1*vxp
               Vx(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) = Vx(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii)   + w2*vxp    
               Vx(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) = Vx(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii)   + w3*vxp
               Vx(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) = Vx(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii)   + w4*vxp    
               Vx(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) = Vx(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii)   + w5*vxp
               Vx(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) = Vx(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii)   + w6*vxp    
               Vx(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) = Vx(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii)   + w7*vxp
               Vx(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) = Vx(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii)   + w8*vxp
               
               Vy(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) = Vy(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii)   + w1*vyp
               Vy(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) = Vy(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii)   + w2*vyp    
               Vy(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) = Vy(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii)   + w3*vyp
               Vy(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) = Vy(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii)   + w4*vyp    
               Vy(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) = Vy(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii)   + w5*vyp
               Vy(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) = Vy(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii)   + w6*vyp    
               Vy(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) = Vy(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii)   + w7*vyp
               Vy(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) = Vy(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii)   + w8*vyp               
               
               Vz(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) = Vz(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii)   + w1*vzp
               Vz(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) = Vz(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii)   + w2*vzp    
               Vz(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) = Vz(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii)   + w3*vzp
               Vz(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) = Vz(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii)   + w4*vzp    
               Vz(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) = Vz(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii)   + w5*vzp
               Vz(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) = Vz(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii)   + w6*vzp    
               Vz(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) = Vz(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii)   + w7*vzp
               Vz(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) = Vz(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii)   + w8*vzp               

             endif  
          else
            if ((particule(nn)%orig /= 0).or.(particule(nn)%exc /= 0)) then 
               density(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) = density(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii)   + w1
               density(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) = density(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii)   + w2    
               density(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) = density(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii)   + w3
               density(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) = density(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii)   + w4    
               density(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) = density(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii)   + w5
               density(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) = density(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii)   + w6    
               density(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) = density(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii)   + w7
               density(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) = density(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii)   + w8

               Vx(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) = Vx(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii)   + w1*vxp
               Vx(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) = Vx(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii)   + w2*vxp    
               Vx(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) = Vx(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii)   + w3*vxp
               Vx(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) = Vx(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii)   + w4*vxp    
               Vx(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) = Vx(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii)   + w5*vxp
               Vx(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) = Vx(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii)   + w6*vxp    
               Vx(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) = Vx(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii)   + w7*vxp
               Vx(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) = Vx(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii)   + w8*vxp
               
               Vy(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) = Vy(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii)   + w1*vyp
               Vy(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) = Vy(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii)   + w2*vyp    
               Vy(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) = Vy(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii)   + w3*vyp
               Vy(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) = Vy(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii)   + w4*vyp    
               Vy(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) = Vy(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii)   + w5*vyp
               Vy(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) = Vy(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii)   + w6*vyp    
               Vy(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) = Vy(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii)   + w7*vyp
               Vy(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) = Vy(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii)   + w8*vyp               
               
               Vz(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) = Vz(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii)   + w1*vzp
               Vz(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) = Vz(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii)   + w2*vzp    
               Vz(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) = Vz(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii)   + w3*vzp
               Vz(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) = Vz(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii)   + w4*vzp    
               Vz(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) = Vz(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii)   + w5*vzp
               Vz(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) = Vz(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii)   + w6*vzp    
               Vz(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) = Vz(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii)   + w7*vzp
               Vz(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) = Vz(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii)   + w8*vzp               
            endif 
          endif
        endif  
      enddo
 
 enddo
 
 Vx = Vx/density
 Vy = Vy/density
 Vz = Vz/density
 
#ifdef HAVE_DEBUG   
 print *,'Calculation of Temperature'
#endif
 
 !-- Now that the bulk speed for each ion species is known, we determine the temperature for each species
  do nn = 1,nptot
 !   print *,'Particle n:',nn
    !dn_temp(:,:,:) = 0.
    sqp = particule(nn)%char !--On collecte la charge de la particule
    smp = particule(nn)%mass !--On collecte la charge de la particule
    qsmp = sqp/smp
    vxp = particule(nn)%vel(1) !-- On collecte la composante Vx de la vitesse
    vyp = particule(nn)%vel(2) !-- On collecte la composante Vy de la vitesse
    vzp = particule(nn)%vel(3) !-- On collecte la composante Vz de la vitesse
       
       !--Relative position in cell s_f=(xf,yf,zf) center of the particule
       s_m = one + (particule(nn)%pos-s_min)*gstep_inv
       
       ijk = int(s_m)
       
#ifdef HAVE_DEBUG   
    if(any(ijk>nc1)) then
     print *,"ERROR Moment temperature",particule(nn)%pos
     stop
    endif
#endif   
       s_f = s_m-real(ijk,dp)
       
       !--Sequence of indices of B at cell corners
       !--Trilinear weight
       s_a = one-s_f
    
       w1 = s_a(1)*s_a(2)*s_a(3)*sqp
       w2 = s_f(1)*s_a(2)*s_a(3)*sqp
       w3 = s_a(1)*s_f(2)*s_a(3)*sqp
       w4 = s_f(1)*s_f(2)*s_a(3)*sqp
       w5 = s_a(1)*s_a(2)*s_f(3)*sqp
       w6 = s_f(1)*s_a(2)*s_f(3)*sqp
       w7 = s_a(1)*s_f(2)*s_f(3)*sqp
       w8 = s_f(1)*s_f(2)*s_f(3)*sqp
       
       !--Number density (in DN)
       !dn_temp(ijk(1)  ,ijk(2)  ,ijk(3)  ) =  w1
       !dn_temp(ijk(1)+1,ijk(2)  ,ijk(3)  ) =  w2
       !dn_temp(ijk(1)  ,ijk(2)+1,ijk(3)  ) =  w3
       !dn_temp(ijk(1)+1,ijk(2)+1,ijk(3)  ) =  w4
       !dn_temp(ijk(1)  ,ijk(2)  ,ijk(3)+1) =  w5
       !dn_temp(ijk(1)+1,ijk(2)  ,ijk(3)+1) =  w6
       !dn_temp(ijk(1)  ,ijk(2)+1,ijk(3)+1) =  w7
       !dn_temp(ijk(1)+1,ijk(2)+1,ijk(3)+1) =  w8
       
       do ii = 1,nb_species
         if (species_info(ii)%qsm_value == qsmp) then
           if (species_info(ii)%origin_value == 0) then
             if ((particule(nn)%orig == 0).and.(particule(nn)%exc == 0)) then
               Temp(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) = Temp(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii)   + w1*( &
                & (Vx(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) - vxp)**2  + &
                & (Vy(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) - vyp)**2  + &
                & (Vz(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) - vzp)**2)  
               Temp(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) = Temp(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii)   + w2*( &
                & (Vx(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) - vxp)**2  + &
                & (Vy(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) - vyp)**2  + &
                & (Vz(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) - vzp)**2)                  
               Temp(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) = Temp(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii)   + w3*( &
                & (Vx(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) - vxp)**2  + &
                & (Vy(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) - vyp)**2  + &
                & (Vz(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) - vzp)**2)             
               Temp(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) = Temp(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii)   + w4*(&
                & (Vx(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) - vxp)**2  + &
                & (Vy(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) - vyp)**2  + &
                & (Vz(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) - vzp)**2)                  
                Temp(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) = Temp(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii)   + w5*(&
                 & (Vx(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) - vxp)**2  + &
                 & (Vy(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) - vyp)**2  + &
                 & (Vz(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) - vzp)**2)  
                Temp(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) = Temp(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii)   + w6*(&
                 & (Vx(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) - vxp)**2  + &
                 & (Vy(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) - vyp)**2  + &
                 & (Vz(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) - vzp)**2)                  
                Temp(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) = Temp(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii)   + w7*(&
                 & (Vx(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) - vxp)**2  + &
                 & (Vy(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) - vyp)**2  + &
                 & (Vz(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) - vzp)**2)             
                Temp(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) = Temp(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii)   + w8*( &
                 & (Vx(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) - vxp)**2  + &
                 & (Vy(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) - vyp)**2  + &
                & (Vz(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) - vzp)**2)
               
 !               Temp(:,:,:,ii) = Temp(:,:,:,ii) + (Vx(:,:,:,ii) - dn_temp(:,:,:)*vxp)**2 + &
 !               	& (Vy(:,:,:,ii) - dn_temp(:,:,:)*vyp)**2 + &
 !               	& (Vz(:,:,:,ii) - dn_temp(:,:,:)*vzp)**2
                
              endif  
           else
             if ((particule(nn)%orig /= 0).or.(particule(nn)%exc /= 0)) then 
!                Temp(:,:,:,ii) = Temp(:,:,:,ii) + (Vx(:,:,:,ii) - dn_temp(:,:,:)*vxp)**2 + &
!                	& (Vy(:,:,:,ii) - dn_temp(:,:,:)*vyp)**2 + &
!                	& (Vz(:,:,:,ii) - dn_temp(:,:,:)*vzp)**2
               Temp(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) = Temp(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii)   + w1*( &
                & (Vx(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) - vxp)**2  + &
                & (Vy(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) - vyp)**2  + &
                & (Vz(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii) - vzp)**2)  
               Temp(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) = Temp(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii)   + w2*( &
                & (Vx(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) - vxp)**2  + &
                & (Vy(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) - vyp)**2  + &
                & (Vz(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii) - vzp)**2)                  
               Temp(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) = Temp(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii)   + w3*( &
                & (Vx(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) - vxp)**2  + &
                & (Vy(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) - vyp)**2  + &
                & (Vz(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii) - vzp)**2)             
               Temp(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) = Temp(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii)   + w4*(&
                & (Vx(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) - vxp)**2  + &
                & (Vy(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) - vyp)**2  + &
                & (Vz(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii) - vzp)**2)                 
                Temp(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) = Temp(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii)   + w5*(&
                 & (Vx(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) - vxp)**2  + &
                 & (Vy(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) - vyp)**2  + &
                 & (Vz(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii) - vzp)**2)  
                Temp(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) = Temp(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii)   + w6*(&
                 & (Vx(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) - vxp)**2  + &
                 & (Vy(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) - vyp)**2  + &
                 & (Vz(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii) - vzp)**2)                  
                Temp(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) = Temp(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii)   + w7*(&
                 & (Vx(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) - vxp)**2  + &
                 & (Vy(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) - vyp)**2  + &
                 & (Vz(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii) - vzp)**2)             
                Temp(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) = Temp(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii)   + w8*( &
                 & (Vx(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) - vxp)**2  + &
                 & (Vy(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) - vyp)**2  + &
                & (Vz(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii) - vzp)**2)
             endif
          endif 
        endif  
       enddo
  
 enddo
 
 


 where (density/=0) 
   Temp = Temp/density
 end where  

#ifdef HAVE_DEBUG   
do ii = 1,nb_species
  print *,' Temp Min & Max  espece ',ii,minval(Temp(:,:,:,ii)),maxval(Temp(:,:,:,ii))
enddo
#endif
 
 !deallocate(dn_temp)
 __WRT_DEBUG_OUT("species_moments")
 
 end subroutine species_moments  

 

!#ifdef HAVE_NETCDF 
 !!=============================================================
 !!subroutine: diag_moment_species/wrt_moment_species
 !! FUNCTION 
 !!  Write information about Particles moments on CDF files for any proc
 !!
 !! OUTPUT
 !!  Only write
 subroutine wrt_moment_species(filwrt)
 character(len=*),intent(in) :: filwrt
 
   integer :: ncid, stId,ii,nb_species,jj
   integer :: dimid(12), varid(100)
   integer,allocatable :: dimspec(:)
   character(len=10),dimension(:),allocatable :: Ion_label_tab
   integer :: dimion(5)
   character(len=40) :: name_file 
   type(Species_characteristic),dimension(:), allocatable :: species_info
   real(dp),dimension(:,:,:,:),allocatable :: density_species,Vx_species,Vy_species,Vz_species,Temp_species
 
   __WRT_DEBUG_IN("wrt_moment_species")
   __GETTIME(9,1) !--Timer start
   
   !--Determine the number of ion species
   call determine_nb_species(planetname,nb_species)
   
   if (mpiinfo%me == 0) print *,'Nb ion species for ',trim(planetname),' is ',nb_species
   
   allocate(species_info(nb_species))
      
  ! determine the charateristic of the different ion species
  call determine_species(planetname,nb_species,species_info)

! correct a bug with some netCDF libraries 
  allocate(Ion_label_tab(nb_species))
  do jj=1,nb_species
        Ion_label_tab(jj)=species_info(jj)%Ion_label
  enddo
      
  allocate(density_species(ncm(1),ncm(2),ncm(3),nb_species));   density_species(:,:,:,:) = zero   
  allocate(Vx_species(ncm(1),ncm(2),ncm(3),nb_species));        Vx_species(:,:,:,:) = zero
  allocate(Vy_species(ncm(1),ncm(2),ncm(3),nb_species));        Vy_species(:,:,:,:) = zero
  allocate(Vz_species(ncm(1),ncm(2),ncm(3),nb_species));        Vz_species(:,:,:,:) = zero
  allocate(Temp_species(ncm(1),ncm(2),ncm(3),nb_species));      Temp_species(:,:,:,:) = zero
  
   call species_moments(particule,nptot,species_info,ncm,nc1,gstep,s_min_loc,nb_species,density_species, &
        & Vx_species,Vy_species,Vz_species,Temp_species)
   
 
   !--Create the file name
   call create_file_name(name_file,filwrt,mpiinfo%me)
 
   stId = nf90_create(name_file , nf90_clobber, ncid)
   call test_cdf(stId)
 
   !--Define the dimensions that will define the size of the file.
   call  create_wrt_dimensions_cdf(ncid,dimid,ncm)
 
   !--Define Speces dimensions
   call species_def_dim_cdf(ncid,dimspec)
   
   !--Define Ion species dimension
   stId = nf90_def_dim(ncid, "name_ion_length" , 10 , dimion(1))
   call test_cdf(stId)
   stId = nf90_def_dim(ncid, "cell_x"          , ncm(1)    , dimion(2))
   call test_cdf(stId)
   stId = nf90_def_dim(ncid, "cell_y"          , ncm(2)    , dimion(3))
   call test_cdf(stId)
   stId = nf90_def_dim(ncid, "cell_z"          , ncm(3)    , dimion(4))
   call test_cdf(stId)
   stId = nf90_def_dim(ncid, "nb_species" , nb_species , dimion(5))
   call test_cdf(stId)
   
   
   !--Set the global attributes
   call set_global_attribute_cdf(ncid,"Moments")
 
   ii = 1
   !--Define commons variables
   call common_def_var_cdf(ncid,varid,dimid,ii)
 
   !--Define MPI_INFO
   call mpiinfo_def_var_cdf(ncid,varid,dimid,ii)
 
   !--Define Speces variables
   call species_def_var_cdf(ncid,varid,dimspec,ii)
 
   !--Define other variables
   !--GRID_INFO
   stId = nf90_def_var(ncid, "nxyzm",  nf90_int,dimid(1), varid(ii))
   call test_cdf(stId); ii = ii+1 
   stId = nf90_def_var(ncid, "ncm_tot",nf90_int,dimid(3), varid(ii))
   call test_cdf(stId); ii = ii+1 
   stId = nf90_def_var(ncid, "nc_tot", nf90_int,dimid(3),varid(ii))
   call test_cdf(stId); ii = ii+1 
   stId = nf90_def_var(ncid, "ncm",    nf90_int,dimid(3),varid(ii))
   call test_cdf(stId); ii = ii+1 
   stId = nf90_def_var(ncid, "s_min",  QP_NF90_DP,dimid(3),  varid(ii))
   call test_cdf(stId); ii = ii+1 
   stId = nf90_def_var(ncid, "s_max",  QP_NF90_DP,dimid(3),  varid(ii))
   call test_cdf(stId); ii = ii+1 
 
   ! !--SPECIES_INFO
   ! stId = nf90_def_var(ncid, "nfl", nf90_int, dimid(1), varid(ii))
   ! call test_cdf(stId); ii = ii+1 
 
   !--PARTICLES_INFO
   stId = nf90_def_var(ncid, "npm",               nf90_int,dimid(1),varid(ii))
   call test_cdf(stId); ii = ii+1 
   
   !--Nb_ion species
   stId = nf90_def_var(ncid, "nb_ion_species",nf90_int,dimid(1),varid(ii))
   call test_cdf(stId); ii = ii+1 
   !--ion species label
   stId = nf90_def_var(ncid, "ion_label",nf90_char,(/dimion(1),dimion(5)/), varid(ii))
   call test_cdf(stId); ii = ii+1
   !-- density for each species
   stId = nf90_def_var(ncid, "density_species",QP_NF90_DP,dimion(2:5), varid(ii))
   call test_cdf(stId); ii = ii+1   
   !-- Vx for each species
   stId = nf90_def_var(ncid, "Vx_species",QP_NF90_DP,dimion(2:5), varid(ii))
   call test_cdf(stId); ii = ii+1   
   !-- Vy for each species
   stId = nf90_def_var(ncid, "Vy_species",QP_NF90_DP,dimion(2:5), varid(ii))
   call test_cdf(stId); ii = ii+1
   !-- Vz for each species
   stId = nf90_def_var(ncid, "Vz_species",QP_NF90_DP,dimion(2:5), varid(ii))
   call test_cdf(stId); ii = ii+1      
   !-- Temperature for each species
   stId = nf90_def_var(ncid, "Temp_species",QP_NF90_DP,dimion(2:5), varid(ii))
   call test_cdf(stId); ii = ii+1      

  
   !--PARTICLES_DEF
   !call def_var_particle_cdf(ncid,varid,dimid,ii)
                                               
   !stId = nf90_put_att(ncid, varid(2), "units", "atomic units")
   !call test_cdf(stId)
 
   !--Switch to write mode
   stId = nf90_enddef(ncid); call test_cdf(stId)
 
   ii = 1   
   !--Write common variables into the file
   call common_put_var_cdf(ncid,varid,ii)
 
   !--Write mpi variables into the file
   call mpiinfo_put_var_cdf(ncid,varid,ii)
 
   !--Write Species infos into the file
   call species_put_var_cdf(Spe,ncid,varid,dimspec,ii)
 
   !--Write the other variables into the file
   stId = nf90_put_var(ncid, varid(ii), nxyzm)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), ncm_tot)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), nc_tot)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), ncm)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), s_min)
   call test_cdf(stId); ii = ii+1              
   stId = nf90_put_var(ncid, varid(ii), s_max)
   call test_cdf(stId); ii = ii+1              
   ! stId = nf90_put_var(ncid, varid(ii), nfl)
   ! call test_cdf(stId); ii = ii+1              
   stId = nf90_put_var(ncid, varid(ii), npm)
   call test_cdf(stId); ii = ii+1        
   stId = nf90_put_var(ncid, varid(ii), nb_species)
   call test_cdf(stId); ii = ii+1      
   stId = nf90_put_var(ncid, varid(ii), Ion_label_tab)
   call test_cdf(stId); ii = ii+1
   stId = nf90_put_var(ncid, varid(ii), density_species)
   call test_cdf(stId); ii = ii+1   
   stId = nf90_put_var(ncid, varid(ii), Vx_species)
   call test_cdf(stId); ii = ii+1   
   stId = nf90_put_var(ncid, varid(ii), Vy_species)
   call test_cdf(stId); ii = ii+1      
   stId = nf90_put_var(ncid, varid(ii), Vz_species)
   call test_cdf(stId); ii = ii+1      
   stId = nf90_put_var(ncid, varid(ii), Temp_species)
   call test_cdf(stId); ii = ii+1   
   
   !--Put particles
  ! call put_var_particle_cdf(particule,nptot,ncid,varid,ii)
 
   !--Close the file
   stId = nf90_close(ncid); call test_cdf(stId)
   
   deallocate(density_species,Vx_species,Vy_species,Vz_species,Temp_species,Ion_label_tab)
 
   __GETTIME(9,2) !--Timer stop
   __WRT_DEBUG_OUT("wrt_moment_species")
  end subroutine wrt_moment_species
 
! #else
! 
!  subroutine wrt_moment_species(filwrt)
! 
!   character(len=*),intent(in) :: filwrt
!   character(len = 40) :: name_file
!   character(len=500) :: msg
! 
!   __WRT_DEBUG_IN("wrt_particles")
!   __GETTIME(9,1) !--Timer start
! 
!   !--Particles
!   !--Create the file name
!   call create_file_name(name_file,filwrt,mpiinfo%me)
! 
!   write(msg,'(a,i4,a,a)')&
!        &' ..Writing process ',mpiinfo%me,&
!        &' on file ',name_file
!   call wrt_double(6,msg,wrtscreen,0)
! 
!   open(UNIT   =  1, &
!        FILE   =  name_file, &
!        ACTION = 'write', &
!        FORM   = 'unformatted', &
!        STATUS = 'unknown', &
!        ACCESS = 'sequential')
! 
!   !--On ecrit en tete du fichier les caracteristiques de la simulation
!   write(unit = 1) mpiinfo%me   
!   write(unit = 1) nxyzm
!   write(unit = 1) nproc,nb_voisins,ncm,ndims
!   write(unit = 1) mpiinfo%dims
!   write(unit = 1) nptot
!   write(unit = 1) ncm_tot
!   write(unit = 1) nc_tot,gstep,s_min,s_max,iter,dt,t,nptot,ns,nfl,npm
!   write(unit = 1) mpiinfo%voisin
! 
!   !--Write Species infos
!   call species_put_var_bin(Spe,unit = 1)
! 
!   !--Write Particles
!   call put_var_particle_bin(particule,nptot,unit = 1)
! 
!   close(unit = 1)
! 
!   write(msg,'(a,i4,a,a)')&
!        &' ..Process ',mpiinfo%me,&
!        &' written on file ',name_file
!   call wrt_double(6,msg,wrtscreen,0)
! 
!   __GETTIME(9,2) !--Timer stop
!   __WRT_DEBUG_OUT("wrt_particles")
!  end subroutine wrt_particles
!  !******************************** END WRT_MOMENT_SPECIES_PROCESSUS **************
! 
! #endif
 
 end module diag_moment_species

