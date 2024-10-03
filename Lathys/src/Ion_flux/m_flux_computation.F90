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
   Read_p3_compute_flux, &
   Ion_escape_outer_box, &
   Ion_escape_wake_slices, &
   Ion_escape_spherical_shell, &
   Ion_flux_shell_map
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

subroutine extract_fluxes(run_name,Esep,Shell_alt,Map)
! This routine call a reading routine for particle file
! creates arrays and stores ion fluxes for each planetary species
  character(len=*),intent(in) :: run_name
  real(dp),intent(inout) :: Esep
  real(dp),intent(inout) :: Shell_alt
  integer,intent(in) :: Map
!--Reading variables
  integer :: rang,nb_procs,ndims,nb_voisin
  integer,allocatable :: dims(:),dimspec(:)
  integer :: nptot_proc
  integer,allocatable :: voisin_proc(:)
  integer,allocatable :: coord_proc(:)


  !--Others
  integer :: iproc,nb_species
  integer  :: stId,ncid,ii,is,nxyzm,iteri,sgn
  integer  :: varid(47),dimid(12)
  integer  :: itmp(3),nc_tot(3),ncm_tot(3),ncm(3)
  real(dp) :: rtmp(3),dt,t,s_min(3),s_max(3),s_min_loc(3),s_max_loc(3),gstep(3),phys_speed,n0,V0,x0,r_planet,s_centr(3)
  logical :: file_e
  character(len=64) :: filename
  character(len=32) :: global_attribute,diag_split_name
  character(len=100) :: test_name
  character(len=500) :: msg
  !character(len=8),intent(inout) :: planet
  real(dp),dimension(:,:,:,:),allocatable :: ion_flux 
  real(dp), dimension(:, :, :, :, :), allocatable:: ion_flux_Energy_separation
  real(dp),dimension(:,:,:),allocatable :: Density_Esup,Ux_Esup,Uy_Esup,Uz_Esup, &
& Temp_Esup
  real(dp),dimension(36,72) :: Flux_shell_in=0._dp,Flux_shell_out=0._dp


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
 allocate(ion_flux_Energy_separation(ncm_tot(1),ncm_tot(2),ncm_tot(3),2, 2))
 ion_flux_Energy_separation(:,:,:,:, :) = zero
 allocate(Density_Esup(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 Density_Esup(:,:,:) = zero
 allocate(Ux_Esup(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 Ux_Esup(:,:,:) = zero
 allocate(Uy_Esup(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 Uy_Esup(:,:,:) = zero
 allocate(Uz_Esup(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 Uz_Esup(:,:,:) = zero
 allocate(Temp_Esup(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 Temp_Esup(:,:,:) = zero


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

 call Read_p3_compute_flux(filename,nb_species,ns,species_info,Esep,&
 ion_flux,ion_flux_Energy_separation,ii,Density_Esup,Ux_Esup,Uy_Esup,Uz_Esup)

 enddo

!**** saving the ion flux map
Ux_Esup = Ux_Esup/Density_Esup
Uy_Esup = Uy_Esup/Density_Esup
Uz_Esup = Uz_Esup/Density_Esup

sgn = -1

Density_Esup = Density_Esup*n0/1.e6
Ux_Esup = sgn*Ux_Esup*v0
Uy_Esup = sgn*Uy_Esup*v0
Uz_Esup = Uz_Esup*v0
!!******* creation of ion flux file
stId = nf90_create('Op_ion_flux_Esup.nc',nf90_clobber,ncid)
call test_cdf(stId)

call create_wrt_dimensions_cdf(ncid,dimid,ncm_tot)
stId = nf90_def_dim(ncid,"dimplanet",20,dimid(12))
call test_cdf(stId)

call set_global_attribute_cdf(ncid,"Ion flux")
ii=1
call common_def_var_cdf(ncid,varid,dimid,ii)

 
stId = nf90_def_var(ncid,"s_centr",QP_NF90_DP,dimid(3),varid(ii))
call test_cdf(stId); ii = ii + 1

stId = nf90_def_var(ncid,"r_planet",QP_NF90_DP,dimid(1),varid(ii))
call test_cdf(stId); ii = ii + 1

stId = nf90_def_var(ncid,"phys_length",QP_NF90_DP,dimid(1),varid(ii))
call test_cdf(stId); ii = ii + 1

stId = nf90_def_var(ncid,"planetname",nf90_char,dimid(12),varid(ii))
call test_cdf(stId); ii = ii + 1

  stId  = nf90_def_var(ncid,"Density",QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1

  stId  = nf90_def_var(ncid,"Ux",QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId  = nf90_def_var(ncid,"Uy",QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId  = nf90_def_var(ncid,"Uz",QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId  = nf90_def_var(ncid,"Temperature",QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1
!
  stId = nf90_enddef(ncid); call test_cdf(stId)
  ii=1
  call common_put_var_cdf(ncid,varid,ii)

  stId = nf90_put_var(ncid,varid(ii),s_centr)
  call test_cdf(stId); ii = ii + 1
 
  stId = nf90_put_var(ncid,varid(ii),r_planet)
  call test_cdf(stId); ii = ii +1

  stId = nf90_put_var(ncid,varid(ii),x0)
  call test_cdf(stId); ii = ii +1

  stId = nf90_put_var(ncid,varid(ii),planetname)
  call test_cdf(stId); ii = ii +1

  stId =  nf90_put_var(ncid,varid(ii),Density_Esup(:,:,:))
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid,varid(ii),Ux_Esup(:,:,:))
  call test_cdf(stId); ii = ii+1
  stId =  nf90_put_var(ncid,varid(ii),Uy_Esup(:,:,:))
  call test_cdf(stId); ii = ii+1
  stId =  nf90_put_var(ncid,varid(ii),Uz_Esup(:,:,:))
  call test_cdf(stId); ii = ii+1
  stId =  nf90_put_var(ncid,varid(ii),Temp_Esup(:,:,:))
  call test_cdf(stId); ii = ii+1
!
  stId = nf90_close(ncid); call test_cdf(stId)


 if (Map == 1) then
   call Ion_flux_shell_map(run_name,nb_species,ns,species_info,Flux_shell_in,Flux_shell_out, &
   Shell_alt,Esep,n0,V0,x0,s_centr,gstep,r_planet,nb_procs)
 endif

 !--Normalizing ion fluxes
 ion_flux = ion_flux*n0*(V0*1.e3)
 ion_flux_Energy_separation = ion_flux_Energy_separation*n0*(V0*1.e3)

 print *,'Maximum ion flux Energy separation',maxval(ion_flux_Energy_separation)

 !-- Compute the total escaping flux from simulation border box
 call Ion_escape_outer_box(nb_species,species_info,ns,ncm_tot,gstep,x0,ion_flux)


 !-- Compute the total escaping flux from a shel at Shell_alt
 call Ion_escape_spherical_shell(nb_species,species_info,ns,ncm_tot,&
        & gstep,x0,ion_flux_Energy_separation,s_centr,r_planet,Shell_alt)
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
subroutine Read_p3_compute_flux(filename,nb_species,ns,species_info,Esep,ion_flux,ion_flux_Energy_separation,rank,Density,Ux,Uy,Uz)
character(len=*),intent(in) :: filename
integer,intent(in)::nb_species,ns,rank
type(Species_characteristic),dimension(:),intent(in) :: species_info
real(dp),dimension(:,:,:,:),intent(inout):: ion_flux
real(dp),dimension(:,:,:,:, :),intent(inout) :: ion_flux_Energy_separation
real(dp),dimension(:,:,:),intent(inout) :: Density,Ux,Uy,Uz
real(dp),intent(in) :: Esep
character(len=27) :: write_name
integer ::iunit




! local variable
integer :: nptot_proc,nn,ii
real(dp),dimension(:),allocatable :: posx,posy,posz,vx,vy,vz,pmass,pcharge,pexch,porig
integer :: stId,ncid,ijk(3),ncm_tot(3)
real(dp) :: gstep(3),sqp,smp,qsmp,vxp,vyp,vzp,s_m(3),gstep_inv(3),s_f(3),s_a(3)
real(dp) ::w1,w2,w3,w4,w5,w6,w7,w8,vtot,w01,w02,w03,w04,w05,w06,w07,w08,nrm_n,nrm_U

integer :: dimid(11),varid(200),sgn

!!-diagnostic of recording O+ ions in the plume
  !--Create file name for prticle
!!  write(write_name,'(a19,i4.4,a5)')"Particle_selection_",rank,"_.dat"
!!iunit=10
!!open(unit=iunit,file = trim(write_name),form = "formatted",status = "unknown")


  !--Open NetCDF file
  stId = nf90_open(trim(filename), nf90_nowrite, ncid)
  call test_cdf(stId)
  print *,'===  Treating file === ',filename

  !--Get Number of Particle of this proc
  call get_simple_variable_cdf(ncid,"nptot",nptot_proc)
! allocate particle arrays for reading
   allocate(posx(nptot_proc),posy(nptot_proc),posz(nptot_proc))
   posx(:) = zero; posy(:)=zero; posz(:)=zero
   allocate(vx(nptot_proc),vy(nptot_proc),vz(nptot_proc))
   vx(:) = zero; vy(:) = zero; vz(:) = zero
   allocate(pmass(nptot_proc),pcharge(nptot_proc),pexch(nptot_proc),porig(nptot_proc))
   pmass(:) = zero; pcharge(:) = zero; pexch(:) = zero; porig(:) = zero

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
  
   call get_simple_variable_cdf(ncid,"phys_density",nrm_n)

   call get_simple_variable_cdf(ncid,"phys_speed",nrm_U)


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


      w01 = s_a(1)*s_a(2)*s_a(3)*sqp
      w02 = s_f(1)*s_a(2)*s_a(3)*sqp
      w03 = s_a(1)*s_f(2)*s_a(3)*sqp
      w04 = s_f(1)*s_f(2)*s_a(3)*sqp
      w05 = s_a(1)*s_a(2)*s_f(3)*sqp
      w06 = s_f(1)*s_a(2)*s_f(3)*sqp
      w07 = s_a(1)*s_f(2)*s_f(3)*sqp
      w08 = s_f(1)*s_f(2)*s_f(3)*sqp

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
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,1, 1) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,1, 1)   + w1
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,1, 1) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,1, 1)   + w2
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,1, 1) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,1, 1)   + w3
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,1, 1) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,1, 1)   + w4
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,1, 1) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,1, 1)   + w5
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,1, 1) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,1, 1)   + w6
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,1, 1) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,1, 1)   + w7
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,1, 1) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,1, 1)   + w8
               else if ((qsmp==1._dp/1._dp).and.(vtot <= Esep)) then
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,1, 2) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,1, 2)   + w1
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,1, 2) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,1, 2)   + w2
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,1, 2) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,1, 2)   + w3
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,1, 2) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3) ,1, 2)   + w4
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,1, 2) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,1, 2)   + w5
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,1, 2) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,1, 2)   + w6
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,1, 2) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,1, 2)   + w7
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,1, 2) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,1, 2)   + w8
               endif
               if ((qsmp==1._dp/16._dp).and.(vtot > Esep)) then
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,2, 1) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,2, 1)   + w1
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,2, 1) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,2, 1)   + w2
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,2, 1) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,2, 1)   + w3
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,2, 1) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,2, 1)   + w4
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,2, 1) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,2, 1)   + w5
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,2, 1) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,2, 1)   + w6
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,2, 1) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,2, 1)   + w7
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,2, 1) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,2, 1)   + w8
               else if ((qsmp==1._dp/1._dp).and.(vtot > Esep)) then
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,2, 2) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)  ,2, 2)   + w1
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,2, 2) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)  ,2, 2)   + w2
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,2, 2) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)  ,2, 2)   + w3
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,2, 2) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)  ,2, 2)   + w4
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,2, 2) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)  ,ijk(3)+1,2, 2)   + w5
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,2, 2) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)  ,ijk(3)+1,2, 2)   + w6
               ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,2, 2) = &
                ion_flux_Energy_separation(ijk(1)  ,ijk(2)+1,ijk(3)+1,2, 2)   + w7
               ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,2, 2) = &
                ion_flux_Energy_separation(ijk(1)+1,ijk(2)+1,ijk(3)+1,2, 2)   + w8
                endif
              ! if (((qsmp==1._dp/1._dp).and.(((posz(nn) > 165).and.(posz(nn)<225.))&
              !      .and.(posx(nn)>200))).and.((porig(nn)/=0).or.(pexch(nn)/=0))) then
!!               if ((qsmp==1._dp/16._dp).and.(((posz(nn) > 90).and.(posz(nn)<125.)).and.(posx(nn)>105))) then
!!                 write(iunit,*) posx(nn),posy(nn),posz(nn),vx(nn),vy(nn),vz(nn),pcharge(nn),porig(nn),pexch(nn)
 !!              endif

            endif
        endif
      enddo

     if ((qsmp==1._dp/16._dp).and.(vtot >= Esep)) then
!-- saving density and velcoity
               Density(ijk(1)  ,ijk(2)  ,ijk(3)  ) = Density(ijk(1)  ,ijk(2)  ,ijk(3)  )   + w01
               Density(ijk(1)+1,ijk(2)  ,ijk(3)  ) = Density(ijk(1)+1,ijk(2)  ,ijk(3)  )   + w02
               Density(ijk(1)  ,ijk(2)+1,ijk(3)  ) = Density(ijk(1)  ,ijk(2)+1,ijk(3)  )   + w03
               Density(ijk(1)+1,ijk(2)+1,ijk(3)  ) = Density(ijk(1)+1,ijk(2)+1,ijk(3)  )   + w04
               Density(ijk(1)  ,ijk(2)  ,ijk(3)+1) = Density(ijk(1)  ,ijk(2)  ,ijk(3)+1)   + w05
               Density(ijk(1)+1,ijk(2)  ,ijk(3)+1) = Density(ijk(1)+1,ijk(2)  ,ijk(3)+1)   + w06
               Density(ijk(1)  ,ijk(2)+1,ijk(3)+1) = Density(ijk(1)  ,ijk(2)+1,ijk(3)+1)   + w07
               Density(ijk(1)+1,ijk(2)+1,ijk(3)+1) = Density(ijk(1)+1,ijk(2)+1,ijk(3)+1)   + w08

               Ux(ijk(1)  ,ijk(2)  ,ijk(3)  ) = Ux(ijk(1)  ,ijk(2)  ,ijk(3)  )   + w01*vxp
               Ux(ijk(1)+1,ijk(2)  ,ijk(3)  ) = Ux(ijk(1)+1,ijk(2)  ,ijk(3)  )   + w02*vxp
               Ux(ijk(1)  ,ijk(2)+1,ijk(3)  ) = Ux(ijk(1)  ,ijk(2)+1,ijk(3)  )   + w03*vxp
               Ux(ijk(1)+1,ijk(2)+1,ijk(3)  ) = Ux(ijk(1)+1,ijk(2)+1,ijk(3)  )   + w04*vxp
               Ux(ijk(1)  ,ijk(2)  ,ijk(3)+1) = Ux(ijk(1)  ,ijk(2)  ,ijk(3)+1)   + w05*vxp
               Ux(ijk(1)+1,ijk(2)  ,ijk(3)+1) = Ux(ijk(1)+1,ijk(2)  ,ijk(3)+1)   + w06*vxp
               Ux(ijk(1)  ,ijk(2)+1,ijk(3)+1) = Ux(ijk(1)  ,ijk(2)+1,ijk(3)+1)   + w07*vxp
               Ux(ijk(1)+1,ijk(2)+1,ijk(3)+1) = Ux(ijk(1)+1,ijk(2)+1,ijk(3)+1)   + w08*vxp

               Uy(ijk(1)  ,ijk(2)  ,ijk(3)  ) = Uy(ijk(1)  ,ijk(2)  ,ijk(3)  )   + w01*vyp
               Uy(ijk(1)+1,ijk(2)  ,ijk(3)  ) = Uy(ijk(1)+1,ijk(2)  ,ijk(3)  )   + w02*vyp
               Uy(ijk(1)  ,ijk(2)+1,ijk(3)  ) = Uy(ijk(1)  ,ijk(2)+1,ijk(3)  )   + w03*vyp
               Uy(ijk(1)+1,ijk(2)+1,ijk(3)  ) = Uy(ijk(1)+1,ijk(2)+1,ijk(3)  )   + w04*vyp
               Uy(ijk(1)  ,ijk(2)  ,ijk(3)+1) = Uy(ijk(1)  ,ijk(2)  ,ijk(3)+1)   + w05*vyp
               Uy(ijk(1)+1,ijk(2)  ,ijk(3)+1) = Uy(ijk(1)+1,ijk(2)  ,ijk(3)+1)   + w06*vyp
               Uy(ijk(1)  ,ijk(2)+1,ijk(3)+1) = Uy(ijk(1)  ,ijk(2)+1,ijk(3)+1)   + w07*vyp
               Uy(ijk(1)+1,ijk(2)+1,ijk(3)+1) = Uy(ijk(1)+1,ijk(2)+1,ijk(3)+1)   + w08*vyp

               Uz(ijk(1)  ,ijk(2)  ,ijk(3)  ) = Uz(ijk(1)  ,ijk(2)  ,ijk(3)  )   + w01*vzp
               Uz(ijk(1)+1,ijk(2)  ,ijk(3)  ) = Uz(ijk(1)+1,ijk(2)  ,ijk(3)  )   + w02*vzp
               Uz(ijk(1)  ,ijk(2)+1,ijk(3)  ) = Uz(ijk(1)  ,ijk(2)+1,ijk(3)  )   + w03*vzp
               Uz(ijk(1)+1,ijk(2)+1,ijk(3)  ) = Uz(ijk(1)+1,ijk(2)+1,ijk(3)  )   + w04*vzp
               Uz(ijk(1)  ,ijk(2)  ,ijk(3)+1) = Uz(ijk(1)  ,ijk(2)  ,ijk(3)+1)   + w05*vzp
               Uz(ijk(1)+1,ijk(2)  ,ijk(3)+1) = Uz(ijk(1)+1,ijk(2)  ,ijk(3)+1)   + w06*vzp
               Uz(ijk(1)  ,ijk(2)+1,ijk(3)+1) = Uz(ijk(1)  ,ijk(2)+1,ijk(3)+1)   + w07*vzp
               Uz(ijk(1)+1,ijk(2)+1,ijk(3)+1) = Uz(ijk(1)+1,ijk(2)+1,ijk(3)+1)   + w08*vzp
endif

enddo

!Ux = Ux/Density
!Uy = Uy/Density
!Uz = Uz/Density
!
!sgn = -1
!
!Density = Density*nrm_n
!Ux = sgn*Ux*nrm_U
!Uy = sgn*Uy*nrm_U
!Uz = Uz*nrm_n
!!******* creation of ion flux file
!stId = nf90_create('Op_ion_flux_Esup.nc',nf90_clobber,ncid)
!call test_cdf(stId)
!
!call create_wrt_dimensions_cdf(ncid,dimid,ncm)
!call set_global_attribute_cdf(ncid,"Ion flux")
!ii=1
!call common_def_var_cdf(ncid,varid,dimid,ii)
! 
!  stId  = nf90_def_var(ncid,"Density",QP_NF90_DP,dimid(6:8),varid(ii))
!  call test_cdf(stId); ii = ii+1
!
!  stId  = nf90_def_var(ncid,"Ux",QP_NF90_DP,dimid(6:8),varid(ii))
!  call test_cdf(stId); ii = ii+1
!  stId  = nf90_def_var(ncid,"Uy",QP_NF90_DP,dimid(6:8),varid(ii))
!  call test_cdf(stId); ii = ii+1
!  stId  = nf90_def_var(ncid,"Uz",QP_NF90_DP,dimid(6:8),varid(ii))
!  call test_cdf(stId); ii = ii+1
!
!  stId = nf90_enddef(ncid); call test_cdf(stId)
!  ii=1
!  call common_put_var_cdf(ncid,varid,ii)
!
!  stId =  nf90_put_var(ncid,varid(ii),Density(:,:,:))
!  call test_cdf(stId); ii = ii+1
!  stId = nf90_put_var(ncid,varid(ii),Ux(:,:,:))
!  call test_cdf(stId); ii = ii+1
!  stId =  nf90_put_var(ncid,varid(ii),Uy(:,:,:))
!  call test_cdf(stId); ii = ii+1
!  stId =  nf90_put_var(ncid,varid(ii),Uz(:,:,:))
!  call test_cdf(stId); ii = ii+1
!
!  stId = nf90_close(ncid); call test_cdf(stId)

!print *,'Max and Min of ion flux :',maxval(ion_flux),minval(ion_flux)
!!close(iunit)
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
allocate(total_escape(nb_species-ns)); total_escape(:)=zero

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
real(dp),intent(in) :: x0,gstep(3),ion_flux_Energy_separation(:,:,:,:, :),xcentr,r_planet
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
       sum_flux = sum_flux + ion_flux_Energy_separation(ii,jj,kk,Energylevel, 1)
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


!
! ==== Computation of Ion escape on spherical shell with
! ==== Energy separation
subroutine Ion_escape_spherical_shell(nb_species,species_info,ns,ncm_tot,&
        & gstep,x0,ion_flux_Energy_separation,s_centr,r_planet,Shell_alt)
integer,intent(in) :: nb_species,ns,ncm_tot(3)
type(Species_characteristic),dimension(:),intent(in) :: species_info
real(dp),intent(in) :: x0,gstep(3),ion_flux_Energy_separation(:,:,:,:,:),s_centr(3),r_planet, Shell_alt
!! Local Variables
integer :: ii,jj,kk,interval,Energylevel,icut
real(dp) :: sum_flux,x,y,z,r,flux_pE,flux_mE,flux_day,flux_night
real(dp) :: dtheta,dphi,surf,thetaval
integer :: nphi,ntheta,itheta,iphi,nr,ir,sp
real(dp),dimension(:),allocatable :: phishell,thetashell
real(dp) :: min_x,min_y,min_z,max_x,max_y,max_z
! -  Determination of the angular resolution depending on the altitude and the
! grid size
dtheta = gstep(1)*x0/(r_planet*x0+Shell_alt)
print *,' Angular resolution of the Shell [rad]  : ',dtheta
print *,'Shell Altitude [km] :',Shell_alt
print *,'r_planet, x0 gstep(1) : ',r_planet,x0,gstep(1)
dphi = dtheta
ntheta = int(pi/dtheta)+1
nphi = int(2*pi/dphi)+1
!averaging fluxes of a thickness shell
nr = 5

 allocate(phishell(nphi));   phishell(:) = 0._dp
 allocate(thetashell(ntheta)); thetashell(:) = 0._dp
do ii=1,nphi
   phishell(ii) = (ii-1)*dphi
enddo
do ii=1, ntheta
   thetashell(ii) = (ii-1)*dtheta
enddo

r = r_planet + Shell_alt/x0

do sp = 1, 2
 sum_flux = 0._dp
 flux_pE = 0._dp
 flux_mE = 0._dp
 flux_day = 0._dp
 flux_night = 0._dp
 min_x = ncm_tot(1)
 min_y = ncm_tot(2)
 min_z = ncm_tot(3)
 max_x = 1
 max_y = 1
 max_z = 1

 do ir = 1,nr
   do itheta = 1,ntheta
     do iphi = 1,nphi
      x = (r+(ir-1)*gstep(1))*cos(phishell(iphi))*sin(thetashell(itheta))+s_centr(1)
      y = (r+(ir-1)*gstep(1))*sin(phishell(iphi))*sin(thetashell(itheta))+s_centr(2)
      z = (r+(ir-1)*gstep(1))*cos(thetashell(itheta))+s_centr(3)
      ii = int(x/gstep(1))+1
      jj = int(y/gstep(2))+1
      kk = int(z/gstep(3))+1
      if (ii < min_x) min_x = ii
      if (jj < min_y) min_y = jj
      if (kk < min_z) min_z = kk
      if (ii>max_x) max_x = ii
      if (jj>max_y) max_y = jj
      if (kk>max_z) then
                max_z = kk
                thetaval = thetashell(itheta)
      endif
     surf =sin(thetashell(itheta))*dtheta*dphi*((r+(ir-1)*gstep(1))*x0*1.e3)**2
     sum_flux = sum_flux+ion_flux_Energy_separation(ii,jj,kk,2, sp)*surf
     if (thetashell(itheta) < 45.*pi/180.) then
         flux_pE = flux_pE + ion_flux_Energy_separation(ii,jj,kk,2, sp)*surf
     elseif (thetashell(itheta) < 135.*pi/180.) then
         if (x>s_centr(1)) then
         flux_night = flux_night+ion_flux_Energy_separation(ii,jj,kk,2, sp)*surf
         else
         flux_day = flux_day+ion_flux_Energy_separation(ii,jj,kk,2, sp)*surf
         endif
       else
         flux_mE = flux_mE + ion_flux_Energy_separation(ii,jj,kk,2, sp)*surf
     endif
    enddo
   enddo
 enddo
if (sp == 1) then
sum_flux = sum_flux/nr
flux_pE = flux_pE/nr
flux_mE = flux_mE/nr
flux_day = flux_day/nr
flux_night = flux_night/nr
!sum_flux = sum_flux/(ntheta*nphi)*4*pi*((r_planet*x0+Shell_alt)*1.e3)**2
print *,' Nb of grid point used', nphi*ntheta
print *,'O+ escape above Esep at a shperical shell: ',sum_flux
!print *,'min (x,y,z) :',min_x,min_y,min_z
print *,'max z, theta :',max_z,thetaval
print *,'O+ escape along +Z_MSO (>45 lat) :',flux_pE
print *,'O+ escape along -Z_MSO (>45 lat) :',flux_mE
print *,'O+ escape along dayside :',flux_day
print *,'O+ escape nightside :',flux_night
else
sum_flux = sum_flux/nr
flux_pE = flux_pE/nr
flux_mE = flux_mE/nr
flux_day = flux_day/nr
flux_night = flux_night/nr
!sum_flux = sum_flux/(ntheta*nphi)*4*pi*((r_planet*x0+Shell_alt)*1.e3)**2
print *,' Nb of grid point used', nphi*ntheta
print *,'H+ escape above Esep at a shperical shell: ',sum_flux
!print *,'min (x,y,z) :',min_x,min_y,min_z
print *,'max z, theta :',max_z,thetaval
print *,'H+ escape along +Z_MSO (>45 lat) :',flux_pE
print *,'H+ escape along -Z_MSO (>45 lat) :',flux_mE
print *,'H+ escape along dayside :',flux_day
print *,'H+ escape nightside :',flux_night
endif
enddo

deallocate(phishell,thetashell)

end subroutine Ion_escape_spherical_shell

!
! ==== Computation of Ion flux on spherical shell
! ====  at a given altitude : separation between
! ==== inward and outward flux
subroutine Ion_flux_shell_map(run_name,nb_species,ns,species_info,Flux_shell_in, &
           Flux_shell_out,Shell_alt,Esep,n0,V0,x0,s_centr,gstep,r_planet,nbprocs)
  character(len=*),intent(in) :: run_name
  real(dp),intent(in) :: Esep,n0,V0,s_centr(3),gstep(3),r_planet,x0
  real(dp),intent(in) :: Shell_alt
  real(dp),dimension(:,:),intent(inout) :: Flux_shell_in, Flux_shell_out
  integer,intent(in)::nb_species,ns,nbprocs
  type(Species_characteristic),dimension(:),intent(in) :: species_info
  !! local Variable
  integer :: ii
  character(len=64) :: filename
  logical :: file_e
  integer :: nptot_proc,nn
  real(dp),dimension(:),allocatable :: posx,posy,posz,vx,vy,vz,pmass,pcharge,pexch,porig
  integer :: stId,ncid,ijk(3),ncm_tot(3)
  real(dp) :: sqp,smp,qsmp,vxp,vyp,vzp,s_m(3),gstep_inv(3),s_f(3),s_a(3)
  real(dp) :: w1,w2,w3,w4,w5,w6,w7,w8,vtot
  real(dp) :: rval,thetaval,phival,xp,yp,zp,rmin,rmax,Ur,volp,dV
  integer :: itheta,iphi
  integer :: dimid(2),varid(2)

  !-- Reading p3 files and extracting ion fluxes for each species
   do ii=0,nbprocs-1
   !--Create file name
     write(filename,'(a3,i4.4,a1,2a)')"p3_",ii,'_',trim(run_name),".nc"
   !--Inquire if the file exists
     inquire( file=trim(filename), exist=file_e )
     !--No file: return
     if(.not.(file_e)) then
      call wrtout(6,"File: "//trim(filename)//" does not exits","PERS")
      return
    endif


!!-diagnostic of recording O+ ions in the plume
  !--Create file name for prticle
!!  write(write_name,'(a19,i4.4,a5)')"Particle_selection_",rank,"_.dat"
!!iunit=10
!!open(unit=iunit,file = trim(write_name),form = "formatted",status = "unknown")


  !--Open NetCDF file
  stId = nf90_open(trim(filename), nf90_nowrite, ncid)
  call test_cdf(stId)
  print *,'===  Treating file === ',filename

  !--Get Number of Particle of this proc
  call get_simple_variable_cdf(ncid,"nptot",nptot_proc)
! allocate particle arrays for reading
   allocate(posx(nptot_proc),posy(nptot_proc),posz(nptot_proc))
   posx(:) = zero; posy(:)=zero; posz(:)=zero
   allocate(vx(nptot_proc),vy(nptot_proc),vz(nptot_proc))
   vx(:) = zero; vy(:) = zero; vz(:) = zero
   allocate(pmass(nptot_proc),pcharge(nptot_proc),pexch(nptot_proc),porig(nptot_proc))
   pmass(:) = zero; pcharge(:) = zero; pexch(:) = zero; porig(:) = zero

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
rmin = 1.+Shell_alt/6052.
rmax = 1.+(Shell_alt+400.)/6052.

do nn = 1,nptot_proc
!   print *,'Particle n:',nn
   !dn_temp(:,:,:) = 0.
   sqp = pcharge(nn) !--On collecte la charge de la particule
   smp = pmass(nn) !--On collecte la masse de la particule
   qsmp = sqp/smp
   ! looking only for O+
   if (qsmp==1._dp/16._dp) then
   ! looking if the particle is in the shell thickness
     xp = -(posx(nn)/r_planet-s_centr(1)/r_planet)
     yp = -(posy(nn)/r_planet-s_centr(2)/r_planet)
     zp =  (posz(nn)/r_planet-s_centr(3)/r_planet)
     rval = sqrt(xp**2+yp**2+zp**2)
     if ((rval .ge. rmin) .and. (rval .le. rmax)) then
     thetaval = acos(zp/rval)
     phival = atan(yp/xp)
     if ((yp >=0._dp) .and. (xp <0.)) phival = pi+phival
     if ((yp >0._dp) .and. (xp==0._dp)) phival = pi/2._dp
     if ((yp <0._dp) .and. (xp ==0._dp)) phival = 3.*pi/2._dp
     if ((yp <0._dp) .and. (xp <0._dp)) phival = phival + pi
     if ((yp <0._dp) .and. (xp >0._dp)) phival = 2.*pi+phival
     itheta = int(thetaval/(5./180.*pi))+1.
     if (itheta > 36 ) itheta=36
     iphi = int(phival/(5./180.*pi))+1
     if (iphi > 72) iphi=72
     vxp = vx(nn) !-- On collecte la composante Vx de la vitesse
     vyp = vy(nn) !-- On collecte la composante Vy de la vitesse
     vzp = vz(nn) !-- On collecte la composante Vz de la vitesse
     Ur = -vxp*cos(phival)*sin(thetaval) - vyp*sin(phival)*sin(thetaval) + &
      vzp*cos(thetaval)
     volp = (gstep(1)*x0*1.e3)**3
     dV = (rmin*6052.e5)**2*sin((itheta-1)*5./180.*pi)*(5.*pi/180.)**2*400.e5
     if (Ur >0._dp) then
        Flux_shell_out(itheta,iphi) = Flux_shell_out(itheta,iphi) + &
        sqp*n0*volp*Ur*V0*1.e5/dV
     else
        Flux_shell_in(itheta,iphi) = Flux_shell_in(itheta,iphi) + &
        sqp*n0*volp*Ur*V0*1.e5/dV
     endif
     !print *,'itheta,iphi,Ur,dV',itheta,iphi,Ur,dV
     endif
   endif

enddo

   deallocate(posx,posy,posz)
   deallocate(vx,vy,vz)
   deallocate(pmass,pcharge,pexch,porig)
 enddo
 !--Create file name for Flux (WRITE)
  write(filename,'(a5,2a)')"MapS_",trim(run_name),".nc"
  !--Create the file
  stId = nf90_create(filename , nf90_clobber, ncid)
  call test_cdf(stId)

  stId = nf90_def_dim(ncid, "nTheta"          , 36    , dimid(1))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "nPhi"            , 72    , dimid(2))
  call test_cdf(stId)

  ii=1
  stId = nf90_def_var(ncid, "Flux_Op_In",       QP_NF90_DP,dimid(1:2),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "Flux_Op_Out",       QP_NF90_DP,dimid(1:2),varid(ii))
  call test_cdf(stId); ii = ii+1


  !--Switch to write mode
  stId = nf90_enddef(ncid); call test_cdf(stId)

  ii = 1
  stId = nf90_put_var(ncid, varid(ii), Flux_Shell_in)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), Flux_Shell_out)
  call test_cdf(stId); ii = ii+1

    !--Close the file
  stId = nf90_close(ncid); call test_cdf(stId)

end subroutine Ion_flux_shell_map

end module m_flux_computation
