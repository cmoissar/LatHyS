!=============================================================
!=============================================================
module m_flux_computation
 
 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype
 use defs_particletype
 use defs_species
 use defs_parametre, only : fildat,gstep,ns,planetname
 use m_writeout
 use diag_moment_species
#ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
#endif
#include "q-p_common.h"

 implicit none
 private :: &
   Read_p3_compute_flux
 public :: extract_fluxes  
   
contains
 !!############################################################

 !********************************************************************
 ! Auteur				:	 RModolo
 ! Date					:	 18/01/14
 ! Institution				:	LATMOS/UVSQ
 ! Derniere modification		:	18/01/14	
 ! Resume	
 ! 
 !********************************************************************   

subroutine extract_fluxes(run_name,Esep)
! This routine call a reading routine for particle file
! creates arrays and stores ion fluxes for each planetary species
  character(len=*),intent(in) :: run_name
  real(dp),intent(inout) :: Esep
!--Reading variables
  integer :: rang,nb_procs,ndims,nb_voisin
  integer,allocatable :: dims(:),dimspec(:)
  integer :: nptot_proc
  integer,allocatable :: voisin_proc(:)
  integer,allocatable :: coord_proc(:)


  !--Others
  integer :: iproc,nb_species
  integer  :: stId,ncid,ii,is,nxyzm,iter
  integer  :: varid(47),dimid(11)
  integer  :: itmp(3),nc_tot(3),ncm_tot(3),ncm(3)
  real(dp) :: rtmp(3),dt,t,s_min(3),s_max(3),s_min_loc(3),s_max_loc(3),gstep(3),phys_speed,n0,V0,x0,r_planet,s_centr(3)
  logical :: file_e
  character(len=64) :: filename
  character(len=32) :: global_attribute,diag_split_name
  character(len=100) :: test_name
  character(len=500) :: msg  
  real(dp),dimension(:,:,:,:),allocatable :: ion_flux, ion_flux_Energy_separation
  
  
  type(Species_characteristic),dimension(:), allocatable :: species_info
  
  
  !--Create file name for Field (READ, proc 0) 
  write(filename,'(a3,i4.4,a1,2a)')"p3_",0,'_',trim(run_name),".nc"
!--Inquire if the file exists
  inquire( file=trim(filename), exist=file_e )
  !--No file: return
  if(.not.(file_e)) then
   call wrtout(6,"File: "//trim(filename)//" does not exits","PERS")
   return
  endif

  !--Open NetCDF file
  stId = nf90_open(trim(filename), nf90_nowrite, ncid)
  call test_cdf(stId)

  !--Control the Date
  stId = nf90_get_att(ncid, nf90_global, "Date", test_name)
  call test_cdf(stId)
  if(fildat /=trim(test_name))  stop "Problem of Dates"

  !--Get the number of processus which generated the file
  call get_simple_variable_cdf(ncid,"nproc",nb_procs)
 ! write(*,*) 'nb_procs : ',nb_procs

  !--Get the rang of the open file 
  call get_simple_variable_cdf(ncid,"mpiinfo_me",rang)

  !--Get topology dimension
  call get_simple_dimens_cdf(ncid,"mpi_ndims",ndims)

  !--Get number of neighbours
  call get_simple_dimens_cdf(ncid,"mpi_nb_voisin",nb_voisin)
  
  !--Get physical speed value
  call get_simple_variable_cdf(ncid,"phys_speed",phys_speed)  

  !--Get number of cells in X,Y,Z
  call get_simple_variable_cdf(ncid,"nc_tot",nc_tot(:))
  call set_grid(nc_tot,5)

  !--Get number of point for fields in Y,Z
  call get_simple_variable_cdf(ncid,"ncm_tot",ncm_tot(:))
  call set_grid(ncm_tot,3) 

  !--Get number point per processus
  call get_simple_variable_cdf(ncid,"ncm",ncm(:))
  call set_grid(ncm,2)
  
  !--Get total number of points
  call get_simple_variable_cdf(ncid,"nxyzm",nxyzm)
  itmp = ncm
  itmp(1) = nxyzm
  call set_grid(itmp,4)

  !--Get the grid step
  call get_simple_variable_cdf(ncid,"gstep",gstep)

  !--Get the dimensions of the box
  !--MIN
  call get_simple_variable_cdf(ncid,"s_min",s_min)
  !--MAX
  call get_simple_variable_cdf(ncid,"s_max",s_max)

  !--Get Time Informations
  !--Iteration
  call get_simple_variable_cdf(ncid,"iter",iter)
  !--Time interval
  call get_simple_variable_cdf(ncid,"dt_t",rtmp(:2))
  dt = rtmp(1);   t = rtmp(2)

  !--Get Number of Particle of this proc
  call get_simple_variable_cdf(ncid,"nptot",nptot_proc)

  !--Get Number of Species
  call get_simple_variable_cdf(ncid,"ns",ns)

  !--Allocate species and get informations
  call alloc_species(ns,Spe)
  call species_get_var_cdf(Spe,ncid)
  
      ! get normalization value
      !--Get normalisation values
      call get_simple_variable_cdf(ncid,"phys_density",n0)
      call get_simple_variable_cdf(ncid,"phys_speed",V0)  
      call get_simple_variable_cdf(ncid,"phys_length",x0)
      
  !-- Get Planet position and radius
  call get_simple_variable_cdf(ncid,"r_planet",r_planet)
  call get_simple_variable_cdf(ncid,"s_centr",s_centr)
  print *,'Normalisation density value ',n0
  print *,'Normalisation speed value ',V0
  print *,'Normalisation length value ',x0
  print *,' Planet radius ',r_planet
  print *,'Planet position ',s_centr
  
  !--Get planetname
  call get_simple_variable_cdf(ncid,"planetname",planetname)


  stId = nf90_close(ncid);  call test_cdf(stId)
  
  
     !-- Convert Energy separation from eV to speed simulation values (for O ions)
   Esep = sqrt(2.*Esep*e_Cb/(16.*amu_pmass))/(phys_speed*1000.)
   print *,'speed separation ',Esep
  
 !! FOR MARS ONLY
 !! To adapt after for other planets
 nb_species =6
 allocate(species_info(nb_species))
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
     
 !! allocation of Ion fluxes array
 allocate(ion_flux(ncm_tot(1),ncm_tot(2),ncm_tot(3),nb_species-ns))
 ion_flux(:,:,:,:) = zero
 allocate(ion_flux_Energy_separation(ncm_tot(1),ncm_tot(2),ncm_tot(3),2))
 ion_flux_Energy_separation(:,:,:,:) = zero

 !-- Reading p3 files and extracting ion fluxes for each species
 do ii=0,nb_procs-1
 !--Create file name 
   write(filename,'(a3,i4.4,a1,2a)')"p3_",ii,'_',trim(run_name),".nc"
 !--Inquire if the file exists
   inquire( file=trim(filename), exist=file_e )
   !--No file: return
   if(.not.(file_e)) then
    call wrtout(6,"File: "//trim(filename)//" does not exits","PERS")
    return
  endif
 
 
 call Read_p3_compute_flux(filename,nb_species,ns,species_info,Esep,ion_flux,ion_flux_Energy_separation,ii)
 
 enddo
 
 !--Normalizing ion fluxes
 ion_flux = ion_flux*n0*(V0*1.e3)
 ion_flux_Energy_separation = ion_flux_Energy_separation*n0*(V0*1.e3)
 
 print *,'Maximum ion flux Energy separation',maxval(ion_flux_Energy_separation)
 
 !-- Compute the total escaping flux from simulation border box
 call Ion_escape_outer_box(nb_species,species_info,ns,ncm_tot,gstep,x0,ion_flux)
 
  !-- Compute the total escaping flux for O+ and for different energy
 !call Ion_escape_wake_slices(nb_species,species_info,ns,ncm_tot,gstep,x0,ion_flux_Energy_separation,s_centr(1),r_planet)
 

 deallocate(species_info)
 deallocate(ion_flux,ion_flux_Energy_separation)

end subroutine extract_fluxes

  !!=============================================================
  !!subroutine: m_flux_computation/Read_p3_compute_flux
  !! FUNCTION 
  !!  Read information about particle on CDF files for a given proc
  !! and collect flux on the ion flux array
  !!
  !! OUTPUT
  !!  Only read
subroutine Read_p3_compute_flux(filename,nb_species,ns,species_info,Esep,ion_flux,ion_flux_Energy_separation,rank)
character(len=*),intent(in) :: filename
integer,intent(in)::nb_species,ns,rank
type(Species_characteristic),dimension(:),intent(in) :: species_info
real(dp),dimension(:,:,:,:),intent(inout) :: ion_flux,ion_flux_Energy_separation
real(dp),intent(in) :: Esep
character(len=27) :: write_name
integer ::iunit




! local variable
integer :: nptot_proc,nn,ii
real(dp),dimension(:),allocatable :: posx,posy,posz,vx,vy,vz,pmass,pcharge,pexch,porig
integer :: stId,ncid,ijk(3),ncm_tot(3)
real(dp) :: gstep(3),sqp,smp,qsmp,vxp,vyp,vzp,s_m(3),gstep_inv(3),s_f(3),s_a(3)
real(dp) :: w1,w2,w3,w4,w5,w6,w7,w8,vtot


!!-diagnostic of recording O+ ions in the plume
  !--Create file name for prticle 
  write(write_name,'(a19,i4.4,a5)')"Particle_selection_",rank,"_.dat"
iunit=10
open(unit=iunit,file = trim(write_name),form = "formatted",status = "unknown")


  !--Open NetCDF file
  stId = nf90_open(trim(filename), nf90_nowrite, ncid)
  call test_cdf(stId)
  print *,'===  Treating file === ',filename
  
  !--Get Number of Particle of this proc
  call get_simple_variable_cdf(ncid,"nptot",nptot_proc)
! allocate particle arrays for reading
   allocate(posx(nptot_proc),posy(nptot_proc),posz(nptot_proc))
   posx(:) = zero;	posy(:)=zero;	posz(:)=zero
   allocate(vx(nptot_proc),vy(nptot_proc),vz(nptot_proc))
   vx(:) = zero; vy(:) = zero; vz(:) = zero
   allocate(pmass(nptot_proc),pcharge(nptot_proc),pexch(nptot_proc),porig(nptot_proc))
   pmass(:) = zero;	pcharge(:) = zero;	pexch(:) = zero;	porig(:) = zero
   
     !--Get the grid step
  call get_simple_variable_cdf(ncid,"gstep",gstep)
  !--Get number of point for fields in Y,Z
  call get_simple_variable_cdf(ncid,"ncm_tot",ncm_tot(:))  
   
   !-- get particle position
   call get_simple_variable_cdf(ncid,"particule_x",posx)
   call get_simple_variable_cdf(ncid,"particule_y",posy)
   call get_simple_variable_cdf(ncid,"particule_z",posz)
   !-- get particle velocity
   call get_simple_variable_cdf(ncid,"particule_vx",vx)
   call get_simple_variable_cdf(ncid,"particule_vy",vy)
   call get_simple_variable_cdf(ncid,"particule_vz",vz)
   !-- get particle mass,charge
   call get_simple_variable_cdf(ncid,"particule_mass",pmass)
   call get_simple_variable_cdf(ncid,"particule_char",pcharge)
   !-- get particle origin and charge exchange number
   call get_simple_variable_cdf(ncid,"particule_exc",pexch)
   call get_simple_variable_cdf(ncid,"particule_orig",porig)
   

 
 stId = nf90_close(ncid);  call test_cdf(stId)

gstep_inv = one/gstep

do nn = 1,nptot_proc
!   print *,'Particle n:',nn
   !dn_temp(:,:,:) = 0.
   sqp = pcharge(nn) !--On collecte la charge de la particule
   smp = pmass(nn) !--On collecte la masse de la particule
   qsmp = sqp/smp
   vxp = vx(nn) !-- On collecte la composante Vx de la vitesse
   vyp = vy(nn) !-- On collecte la composante Vy de la vitesse
   vzp = vz(nn) !-- On collecte la composante Vz de la vitesse
   vtot = sqrt(vxp**2+vyp**2+vzp**2)
      
   
   !-- position in cell s_f=(xf,yf,zf) center of the particule
      s_m(1) = one + posx(nn)*gstep_inv(1)
      s_m(2) = one + posy(nn)*gstep_inv(2)
      s_m(3) = one + posz(nn)*gstep_inv(3)
      
      ijk = int(s_m)
      
#ifdef HAVE_DEBUG   
   if(any(ijk>ncm_tot-1)) then
    print *,"ERROR Moment",posx(nn),posy(nn),posz(nn)
    stop
   endif
#endif   
      s_f = s_m-real(ijk,dp)
      
      !--Sequence of indices of B at cell corners
      !--Trilinear weight
      s_a = one-s_f
   
      w1 = s_a(1)*s_a(2)*s_a(3)*sqp*vtot
      w2 = s_f(1)*s_a(2)*s_a(3)*sqp*vtot
      w3 = s_a(1)*s_f(2)*s_a(3)*sqp*vtot
      w4 = s_f(1)*s_f(2)*s_a(3)*sqp*vtot
      w5 = s_a(1)*s_a(2)*s_f(3)*sqp*vtot
      w6 = s_f(1)*s_a(2)*s_f(3)*sqp*vtot
      w7 = s_a(1)*s_f(2)*s_f(3)*sqp*vtot
      w8 = s_f(1)*s_f(2)*s_f(3)*sqp*vtot   


      do ii = ns+1,nb_species
        if (species_info(ii)%qsm_value == qsmp) then  
            if ((porig(nn)/= 0).or.(pexch(nn) /= 0)) then 
               ion_flux(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii-ns) = ion_flux(ijk(1)  ,ijk(2)  ,ijk(3)  ,ii-ns)   + w1
               ion_flux(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii-ns) = ion_flux(ijk(1)+1,ijk(2)  ,ijk(3)  ,ii-ns)   + w2    			
               ion_flux(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii-ns) = ion_flux(ijk(1)  ,ijk(2)+1,ijk(3)  ,ii-ns)   + w3
               ion_flux(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii-ns) = ion_flux(ijk(1)+1,ijk(2)+1,ijk(3)  ,ii-ns)   + w4    			
               ion_flux(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii-ns) = ion_flux(ijk(1)  ,ijk(2)  ,ijk(3)+1,ii-ns)   + w5
               ion_flux(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii-ns) = ion_flux(ijk(1)+1,ijk(2)  ,ijk(3)+1,ii-ns)   + w6    			
               ion_flux(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii-ns) = ion_flux(ijk(1)  ,ijk(2)+1,ijk(3)+1,ii-ns)   + w7
               ion_flux(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii-ns) = ion_flux(ijk(1)+1,ijk(2)+1,ijk(3)+1,ii-ns)   + w8
               ! -- Energy separation for O ions
               if ((qsmp==1._dp/16._dp).and.(vtot <= Esep)) then
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,1) = ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,1)   + w1
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,1) = ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,1)   + w2    			
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,1) = ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,1)   + w3
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,1) = ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,1)   + w4    			
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,1) = ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,1)   + w5
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,1) = ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,1)   + w6    			
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,1) = ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,1)   + w7
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,1) = ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,1)   + w8
               endif
               if ((qsmp==1._dp/16._dp).and.(vtot > Esep)) then
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,2) = ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,2)   + w1
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,2) = ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,2)   + w2    			
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,2) = ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,2)   + w3
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,2) = ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,2)   + w4    			
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,2) = ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,2)   + w5
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,2) = ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,2)   + w6    			
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,2) = ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,2)   + w7
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,2) = ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,2)   + w8
               endif  
              ! if (((qsmp==1._dp/1._dp).and.(((posz(nn) > 165).and.(posz(nn)<225.)).and.(posx(nn)>200))).and.((porig(nn)/=0).or.(pexch(nn)/=0))) then
               if ((qsmp==1._dp/16._dp).and.(((posz(nn) > 90).and.(posz(nn)<125.)).and.(posx(nn)>105))) then
                 write(iunit,*) posx(nn),posy(nn),posz(nn),vx(nn),vy(nn),vz(nn),pcharge(nn),porig(nn),pexch(nn)
               endif  
               
            endif
        endif
      enddo 
   
enddo  

!print *,'Max and Min of ion flux :',maxval(ion_flux),minval(ion_flux)
close(iunit)
   deallocate(posx,posy,posz)
   deallocate(vx,vy,vz)
   deallocate(pmass,pcharge,pexch,porig)
   
end subroutine Read_p3_compute_flux

!!===============================================================
!!subroutine: m_flux_computation/Ion_flux_outer_box
!! FUNCTION 
!!  compute ion escape over the different simulation faces

subroutine Ion_escape_outer_box(nb_species,species_info,ns,ncm_tot,gstep,x0,ion_flux)
integer,intent(in) :: nb_species,ns,ncm_tot(3)
type(Species_characteristic),dimension(:),intent(in) :: species_info
real(dp),intent(in) :: x0,gstep(3),ion_flux(:,:,:,:)
!! Local Variables
real(dp),dimension(:),allocatable :: total_escape
integer :: ii,jj,kk,interval,species
real(dp) :: sum_flux

!-- alocate total escape fluxes
allocate(total_escape(nb_species-ns));	total_escape(:)=zero

interval=5

!-- Computation in the wake X=xmax-5*delta -> X=xmax
write(*,*) ' '
write(*,*) '=============================='
write(*,*) ' Escaping flux for X=Xmax'
write(*,*) ' '
do species=1,nb_species-ns
sum_flux = 0._dp
do ii=ncm_tot(1)-1-interval,ncm_tot(1)-1
  do jj=1,ncm_tot(2)-1
    do kk=1,ncm_tot(3)-1
       sum_flux = sum_flux + ion_flux(ii,jj,kk,species)      
    enddo
  enddo
enddo  
sum_flux = sum_flux/(interval+1)*(gstep(2)*x0*1.e3)*(gstep(3)*x0*1.e3) 
print *,'Ion escape for ',trim(species_info(species+ns)%Ion_label),' ions'
print *,sum_flux,' ions/s'
total_escape(species) = total_escape(species) + sum_flux
enddo

write(*,*) ' '
write(*,*) '=============================='
write(*,*) ' Escaping flux for Z=Zmax'
write(*,*) ' '
do species=1,nb_species-ns
sum_flux = 0._dp
do ii=1,ncm_tot(1)-1
  do jj=1,ncm_tot(2)-1
    do kk=ncm_tot(3)-1-interval,ncm_tot(3)-1
       sum_flux = sum_flux + ion_flux(ii,jj,kk,species)      
    enddo
  enddo
enddo  
sum_flux = sum_flux/(interval+1)*(gstep(2)*x0*1.e3)*(gstep(3)*x0*1.e3) 
print *,'Ion escape for ',trim(species_info(species+ns)%Ion_label),' ions'
print *,sum_flux,' ions/s'
total_escape(species) = total_escape(species) + sum_flux
enddo

write(*,*) ' '
write(*,*) '=============================='
write(*,*) ' Escaping flux for Y=Ymax'
write(*,*) ' '
do species=1,nb_species-ns
sum_flux = 0._dp
do ii=1,ncm_tot(1)-1
  do jj=ncm_tot(2)-1-interval,ncm_tot(2)-1
    do kk=1,ncm_tot(3)-1
       sum_flux = sum_flux + ion_flux(ii,jj,kk,species)      
    enddo
  enddo
enddo  
sum_flux = sum_flux/(interval+1)*(gstep(2)*x0*1.e3)*(gstep(3)*x0*1.e3) 
print *,'Ion escape for ',trim(species_info(species+ns)%Ion_label),' ions'
print *,sum_flux,' ions/s'
total_escape(species) = total_escape(species) + sum_flux
enddo


write(*,*) ' '
write(*,*) '=============================='
write(*,*) ' Escaping flux for Z=Zmin'
write(*,*) ' '
do species=1,nb_species-ns
sum_flux = 0._dp
do ii=1,ncm_tot(1)-1
  do jj=1,ncm_tot(2)-1
    do kk=1,interval+1
       sum_flux = sum_flux + ion_flux(ii,jj,kk,species)      
    enddo
  enddo
enddo  
sum_flux = sum_flux/(interval+1)*(gstep(2)*x0*1.e3)*(gstep(3)*x0*1.e3) 
print *,'Ion escape for ',trim(species_info(species+ns)%Ion_label),' ions'
print *,sum_flux,' ions/s'
total_escape(species) = total_escape(species) + sum_flux
enddo

write(*,*) ' '
write(*,*) '=============================='
write(*,*) ' Escaping flux for Y=Ymin'
write(*,*) ' '
do species=1,nb_species-ns
sum_flux = 0._dp
do ii=1,ncm_tot(1)-1
  do jj=1,interval+1
    do kk=1,ncm_tot(3)-1
       sum_flux = sum_flux + ion_flux(ii,jj,kk,species)      
    enddo
  enddo
enddo  
sum_flux = sum_flux/(interval+1)*(gstep(2)*x0*1.e3)*(gstep(3)*x0*1.e3) 
print *,'Ion escape for ',trim(species_info(species+ns)%Ion_label),' ions'
print *,sum_flux,' ions/s'
total_escape(species) = total_escape(species) + sum_flux
enddo


write(*,*) ' '
write(*,*) '=============================='
write(*,*) '=============================='
write(*,*) ' Total escape'
write(*,*) ' '
do species=1,nb_species-ns
print *,'Ion escape for ',trim(species_info(species+ns)%Ion_label),' ions'
print *,total_escape(species),' ions/s'
enddo



end subroutine Ion_escape_outer_box


!!===============================================================
!!subroutine: m_flux_computation/Ion_flux_wake_slices
!! FUNCTION 
!!  compute ion escape over the different slices in the wake
subroutine Ion_escape_wake_slices(nb_species,species_info,ns,ncm_tot,gstep,x0,ion_flux_Energy_separation,xcentr,r_planet)
integer,intent(in) :: nb_species,ns,ncm_tot(3)
type(Species_characteristic),dimension(:),intent(in) :: species_info
real(dp),intent(in) :: x0,gstep(3),ion_flux_Energy_separation(:,:,:,:),xcentr,r_planet
!! Local Variables
integer :: ii,jj,kk,interval,Energylevel,icut
real(dp) :: sum_flux,xcut,xmax


interval=6

xcut = xcentr+1.*r_planet
xmax = (ncm_tot(1)-1)*gstep(1)

write(*,*) '###############################'
write(*,*) ' '
do while (xcut<xmax)
write(*,*) ' '
write(*,*) '=============================='
write(*,*) ' Escaping flux for X='
write(*,*) ' '
icut = int(xcut/gstep(1)-interval/2.)
do Energylevel=1,2
sum_flux = 0._dp
do ii=icut-int(interval/2.),icut+int(interval/2.)
  do jj=1,ncm_tot(2)-1
    do kk=1,ncm_tot(3)-1
       sum_flux = sum_flux + ion_flux_Energy_separation(ii,jj,kk,Energylevel)      
    enddo
  enddo
enddo  
sum_flux = sum_flux/(interval+1)*(gstep(2)*x0*1.e3)*(gstep(3)*x0*1.e3) 
print *,'Ion escape for Energie ',Energylevel,' at X=',(xcut-xcentr)/r_planet,' Rm'
print *,sum_flux,' ions/s'
enddo
xcut = xcut + 0.5*r_planet

enddo !while
end subroutine Ion_escape_wake_slices
 

end module m_flux_computation
