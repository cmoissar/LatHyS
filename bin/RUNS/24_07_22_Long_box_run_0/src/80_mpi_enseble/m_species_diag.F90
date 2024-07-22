!=============================================================
!=============================================================
module m_species_diag

 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype
 use defs_particletype
 use defs_species
 use m_writeout
 use defs_parametre
! use m_verify_reassemble_mpi
#ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
#endif
#include "q-p_common.h"

 implicit none
 private 


 public ::                &
      determine_nb_species, &
      determine_species, &
      species_part_count, &
      species_moments

 type,public :: Species_characteristic
    character(len=10) :: Ion_label
    real(dp)          :: qsm_value
    integer       :: origin_value
    integer       :: CE_value
 end type Species_characteristic
 
contains
 !!############################################################

 !********************************************************************
 ! Auteur				:	 MMancini RModolo
 ! Date					:	 27/05/11
 ! Institution				:	CETP/CNRS/IPSL
 ! Derniere modification		:	27/05/10	
 ! Resume	
 ! 
 !********************************************************************
 !!########################################################
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


  
  select case(trim(planetname))
  case("mars")
    nb_species = 6 
  case("titan","mercure")
     nb_species = 3
  case("ganymede")
     nb_species = 5
  case default
     write(*,*) &
         "ERROR: Selected Species determination:",&
         "Planet '",trim(planetname),"' does Not Exist"
     stop
   end select    
     
  
 end subroutine determine_nb_species
 
 !!########################################################
 !! subroutine: m_split_cdf/determine_species
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

 
 select case(trim(planetname))
  case("mars")
    !  H+ sw
    species_info(1)%Ion_label = "Hsw"
    species_info(1)%qsm_value = 1.
    species_info(1)%origin_value = 0
    species_info(1)%CE_value = 0
    
    !He++ sw
    species_info(2)%Ion_label = "Hesw"
    species_info(2)%qsm_value = 2./4.
    species_info(2)%origin_value = 0
    species_info(2)%CE_value = 0
    
    !O+ planetary
    species_info(3)%Ion_label = "Opl"
    species_info(3)%qsm_value = 1./16.
    species_info(3)%origin_value = 1
    species_info(3)%CE_value = 0
    
    !O2+ planetary
    species_info(4)%Ion_label = "O2pl"
    species_info(4)%qsm_value = 1./32.
    species_info(4)%origin_value = 1
    species_info(4)%CE_value = 0

    !CO2+ planetary
    species_info(5)%Ion_label = "CO2pl"
    species_info(5)%qsm_value = 1./44.
    species_info(5)%origin_value = 1
    species_info(5)%CE_value = 0

    !H+ planetary
    species_info(6)%Ion_label = "Hpl"
    species_info(6)%qsm_value = 1./1.
    species_info(6)%origin_value = 1
    species_info(6)%CE_value = 0    
   case("mercure")
    !  H+ sw
    species_info(1)%Ion_label = "Hsw"
    species_info(1)%qsm_value = 1./1.
    species_info(1)%origin_value = 0
    species_info(1)%CE_value = 0
    
    !He++ sw
    species_info(2)%Ion_label = "Hesw"
    species_info(2)%qsm_value = 1./2.
    species_info(2)%origin_value = 0
    species_info(2)%CE_value = 0
    
    !  H+ pl
    species_info(3)%Ion_label = "Hpl"
    species_info(3)%qsm_value = 1./1.
    species_info(3)%origin_value = 1
    species_info(3)%CE_value = 0    
  case("ganymede")
    !  O+ jv
    species_info(1)%Ion_label = "Ojv"
    species_info(1)%qsm_value = 1./1.
    species_info(1)%origin_value = 0
    species_info(1)%CE_value = 0
    
    !  H+ jv
    species_info(2)%Ion_label = "Hjv"
    species_info(2)%qsm_value = 16./1.
    species_info(2)%origin_value = 0
    species_info(2)%CE_value = 0
    
    !  H+ pl
    species_info(3)%Ion_label = "Hjv"
    species_info(3)%qsm_value = 16./1.
    species_info(3)%origin_value = 1
    species_info(3)%CE_value = 0    
    
    !  O+ pl
    species_info(4)%Ion_label = "Opl"
    species_info(4)%qsm_value = 1./1.
    species_info(4)%origin_value = 1
    species_info(4)%CE_value = 0  
    
    !  H2+ pl
    species_info(5)%Ion_label = "H2pl"
    species_info(5)%qsm_value = 16./2.
    species_info(5)%origin_value = 1
    species_info(5)%CE_value = 0     
    
    !  OH+ pl
    !species_info(6)%Ion_label = "OHjv"
    !species_info(6)%qsm_value = 16./17.
    !species_info(6)%origin_value = 1
    !species_info(6)%CE_value = 0  
    
   case default
   write(msg,'(a,3x,a,4x,3a)')ch10,&
        "ERROR: Selected Species determination:",&
        "Planet '",trim(Spe%planetname),"' does Not Exist"
   stop
  end select    

 end subroutine determine_species

!!########################################################
!! subroutine: m_species_diag/species_part_count
!!  FUNCTION
!! returns the number of particle for a given species
!! INPUT 
!!  particle  =  particle information
!!  qsm = charge over mass criterai
!!  org = origin criteria
!! OUTPUT 
!!  nb_species : number of Ion species present
!! species_info = charactersitic information for each ion species
  
 subroutine species_part_count(particle,nptot,qsm,org,CE,count_nb)
 
 integer,intent(in)            :: nptot
 type(particletype),dimension(nptot),intent(in) :: particle
 real(dp),intent(in)           :: qsm
 integer,intent(in)            :: org,CE
 integer,intent(inout)         :: count_nb

 if (org == 0) then
   count_nb = count(((particle(:)%char/particle(:)%mass == qsm) .and. (particle(:)%orig == org)).and. (particle(:)%exc == 0))
 else
   count_nb = count(((particle(:)%char/particle(:)%mass == qsm) .and. (particle(:)%orig == org)).or. (particle(:)%exc /= 0))
 endif  
 
 
 end subroutine species_part_count


!!########################################################
!! subroutine: m_species_diag/species_index
!!  FUNCTION
!! returns the number of particle for a given species
!! INPUT 
!!  particle  =  particle information
!!  qsm = charge over mass criterai
!!  org = origin criteria
!! OUTPUT 
!!  index_part = index of selected particle

  
! subroutine species_index(particle,nptot,qsm,org,CE,nb_part,index_selec)
 
! integer,inten(in)             :: nptot,nb_part
! type(particletype),dimension(nptot),intent(in) :: particule
! integer,dimension(nb_part)    :: index_selec
! real(dp),intent(in)           :: qsm
! integer,intent(in)            :: org,CE
 
 
! mask_selec = 
! PACK(
! 
! end subroutine species_index

!!########################################################
!! subroutine: m_species_diag/species_moments
!!  FUNCTION
!! returns the number of particle for a given species
!! INPUT 
!!  particle  =  particle information
!!  qsm = charge over mass criterai
!!  org = origin criteria
!! OUTPUT 
!!  density_species = density array for each ion species
 
 subroutine species_moments(particule,nptot,species_info,ncm,gstep,s_min,nb_species,density)
 
 integer,intent(in) :: nptot,nb_species
 integer,dimension(3),intent(in) :: ncm
 type(particletype),dimension(nptot),intent(in) :: particule
 real(dp),intent(in)           :: gstep(3),s_min(3)
 real(dp),dimension(ncm(1),ncm(2),ncm(3),nb_species),intent(inout) :: density
 type(species_characteristic),intent(in) :: species_info(nb_species)
 ! local variable
 integer :: nn,ijk(3),ii
 real(dp) :: sqp,smp,qsmp
 real(dp) :: w1,w2,w3,w4,w5,w6,w7,w8
 real(dp),dimension(:,:,:),allocatable :: dn_temp
 real(dp),dimension(3) :: s_f,s_a,v_p,s_m,gstep_inv
 
 allocate(dn_temp(ncm(1),ncm(2),ncm(3)));	
 
 gstep_inv = one/gstep
 
 do nn = 1,nptot
!   print *,'Particle n:',nn
   dn_temp(:,:,:) = 0.
   sqp = particule(nn)%char !--On collecte la charge de la particule
   smp = particule(nn)%mass !--On collecte la masse de la particule
   qsmp = sqp/smp
      
      !--Relative position in cell s_f=(xf,yf,zf) center of the particule
      s_m = one + (particule(nn)%pos-s_min)*gstep_inv
      
      ijk = int(s_m)
   
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
      dn_temp(ijk(1)  ,ijk(2)  ,ijk(3)  ) =  w1
      dn_temp(ijk(1)+1,ijk(2)  ,ijk(3)  ) =  w2
      dn_temp(ijk(1)  ,ijk(2)+1,ijk(3)  ) =  w3
      dn_temp(ijk(1)+1,ijk(2)+1,ijk(3)  ) =  w4
      dn_temp(ijk(1)  ,ijk(2)  ,ijk(3)+1) =  w5
      dn_temp(ijk(1)+1,ijk(2)  ,ijk(3)+1) =  w6
      dn_temp(ijk(1)  ,ijk(2)+1,ijk(3)+1) =  w7
      dn_temp(ijk(1)+1,ijk(2)+1,ijk(3)+1) =  w8
      
      do ii = 1,nb_species
        if (species_info(ii)%qsm_value == qsmp) then
          if (species_info(ii)%origin_value == 0) then
            if ((particule(nn)%orig == 0).and.(particule(nn)%exc == 0)) density(:,:,:,ii) = density(:,:,:,ii) + dn_temp(:,:,:)
          else
            if ((particule(nn)%orig /= 0).or.(particule(nn)%exc /= 0)) density(:,:,:,ii) = density(:,:,:,ii) + dn_temp(:,:,:)
          endif
        endif  
      enddo
 
 enddo
 
 deallocate(dn_temp)
 
 end subroutine species_moments
  



end module m_species_diag
