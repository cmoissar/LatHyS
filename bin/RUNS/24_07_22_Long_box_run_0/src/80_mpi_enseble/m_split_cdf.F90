!=============================================================
!=============================================================
module m_split_cdf

 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype
 use defs_particletype
 use defs_species
 use defs_parametre, only : fildat,gstep,ns
 use m_writeout
 use m_verify_reassemble_mpi
#ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
#endif
#include "q-p_common.h"

 implicit none
 private 

#ifdef HAVE_NETCDF
 private ::               &
      wrt_field_diag,     &
      create_file_diag,   &
      wrt_field_diag_1_array

 public ::                &
      split_cdf
contains
 !!############################################################

 !********************************************************************
 ! Auteur				:	 MMancini RModolo
 ! Date					:	 23/06/03
 ! Institution				:	CETP/CNRS/IPSL
 ! Derniere modification		:	19/07/10	
 ! Resume	
 ! 
 !********************************************************************

 subroutine split_cdf(run_name)
! this routine split the c3 file created during the simulation to 4 files
! 1 file mag : including only magnetic arrays
! 1 file elec : including only electric arrays
! 1 file vel  : including only velocity (bulk) arrays
! 1 file dens : including only density (electron number) array
! after the creation for each process it recombines in one unique file
  character(len=*),intent(in) :: run_name
  !--Reading variables
  integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
  integer :: rang,nb_procs,ndims,nb_voisin
  integer,allocatable :: dims(:),dimspec(:)
  integer :: nptot_proc
  integer,allocatable :: voisin_proc(:)
  integer,allocatable :: coord_proc(:)


  !--Variables to read fields diagnostic
  real(dp),dimension(:,:,:),allocatable :: bx_proc,by_proc,bz_proc, &
       ex_proc,ey_proc,ez_proc,dn_proc,ux_proc,uy_proc,uz_proc,uxa_proc,uya_proc, &
       uza_proc,dna_proc,pe_proc

  !--Others
  integer :: ind,deby,debz,finy,finz,iproc
  integer :: deb,fin,cumul
  integer  :: stId,ncid,ii,is,nxyzm
  integer  :: varid(47),dimid(11)
  integer  :: itmp(3),nc_tot(3),ncm_tot(3),ncm(3)
  real(dp) :: rtmp(3),dt,t,s_min(3),s_max(3),s_min_loc(3),s_max_loc(3)
  logical :: file_e
  character(len=50) :: write_name
  character(len=64) :: filename,var_name1,var_name2,var_name3
  character(len=32) :: global_attribute,diag_split_name
  character(len=100) :: test_name
  character(len=40)  :: diag_name_file
  character(len=500) :: msg

  !--Create file name for Field (READ, proc 0) 
  write(filename,'(a3,i4.4,a1,2a)')"c3_",0,'_',trim(run_name),".nc"
 

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

  allocate(dims(ndims)) ; dims = 0

  !--Define coordinate vector for proc
  allocate(coord_proc(ndims)) ; coord_proc =0

  !--Get the grid step
  call get_simple_variable_cdf(ncid,"gstep",gstep)

  !--Get the dimensions of the box
  !--MIN
  call get_simple_variable_cdf(ncid,"s_min",s_min)
  !--MAX
  call get_simple_variable_cdf(ncid,"s_max",s_max)
  !--Get the dimensions of the process box
  !--MIN_LOC
  call get_simple_variable_cdf(ncid,"s_min_loc",s_min_loc)
  !--MAX_LOC
  call get_simple_variable_cdf(ncid,"s_max_loc",s_max_loc)

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

  allocate(voisin_proc(nb_voisin)) ; voisin_proc = 0

  stId = nf90_close(ncid);  call test_cdf(stId)

  ! Maintenat que l'on connait le nombre de processus qui a tourne
  ! on lit les autres fichiers
  
  do iproc=1,nb_procs 

   !--Create file name
   write(filename,'(a3,i4.4,a1,a)')"c3_",iproc-1,'_',trim(run_name)
   stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
   write(*,*) "Reading file :",filename
   call test_cdf(stId)
  !--allocation des tableaux de lectures du fichier diag de champ
  !--Get number point per processus
  call get_simple_variable_cdf(ncid,"ncm",ncm(:))
  call set_grid(ncm,2)

  allocate(bx_proc(ncm(1),ncm(2),ncm(3)))
  allocate(by_proc(ncm(1),ncm(2),ncm(3)))
  allocate(bz_proc(ncm(1),ncm(2),ncm(3)))
  allocate(ex_proc(ncm(1),ncm(2),ncm(3)))
  allocate(ey_proc(ncm(1),ncm(2),ncm(3)))
  allocate(ez_proc(ncm(1),ncm(2),ncm(3)))
  allocate(dn_proc(ncm(1),ncm(2),ncm(3)))
  allocate(ux_proc(ncm(1),ncm(2),ncm(3)))
  allocate(uy_proc(ncm(1),ncm(2),ncm(3)))
  allocate(uz_proc(ncm(1),ncm(2),ncm(3)))
  allocate(dna_proc(ncm(1),ncm(2),ncm(3)))
  allocate(uxa_proc(ncm(1),ncm(2),ncm(3)))
  allocate(uya_proc(ncm(1),ncm(2),ncm(3)))
  allocate(uza_proc(ncm(1),ncm(2),ncm(3)))
  allocate(pe_proc(ncm(1),ncm(2),ncm(3)))
  bx_proc = zero ; by_proc = zero ; bz_proc = zero
  ex_proc = zero ; ey_proc = zero ; ez_proc = zero
  ux_proc = zero ; uy_proc = zero ; uz_proc = zero
  uxa_proc = zero ; uya_proc = zero ; uza_proc = zero
  dn_proc = zero ; dna_proc = zero; pe_proc = zero

   !--Get Total Number of Particles
   call get_simple_variable_cdf(ncid,"nptot",nptot_proc)
   call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(:) )
   call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(:ndims))
   call get_simple_variable_cdf(ncid,"Bfield_x" ,bx_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"Bfield_y" ,by_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"Bfield_z" ,bz_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"Efield_x" ,ex_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"Efield_y" ,ey_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"Efield_z" ,ez_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"dn"       ,dn_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"dna"     ,dna_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"vel_x"    ,ux_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"vel_y"    ,uy_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"vel_z"    ,uz_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"vela_x"  ,uxa_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"vela_y"  ,uya_proc(:,:,:) )
   call get_simple_variable_cdf(ncid,"vela_z"  ,uza_proc(:,:,:) )
  ! added to compute electronic pressure
   call get_simple_variable_cdf(ncid,"Pe"  ,pe_proc(:,:,:) )

   stId = nf90_close(ncid);  call test_cdf(stId)
  ! writing the magnetic field components in a file
  diag_split_name = 'Mag3'
  global_attribute = "Bfield"
  var_name1 = "Bx"
  var_name2 = "By"
  var_name3 = "Bz"
  
  call create_file_diag(diag_name_file,run_name,diag_split_name,iproc-1)

   write(*,*) diag_name_file
  call wrt_field_diag(diag_name_file,bx_proc,by_proc,bz_proc, &
  trim(var_name1),trim(var_name2),trim(var_name3),& 
  nb_procs,nc_tot,ncm_tot,ncm,nxyzm,dt,t,iter, &
  nptot_proc,ns,coord_proc,voisin_proc,rang,& 
  ndims,nb_voisin,s_min,s_max,s_min_loc,s_max_loc,Spe)
        
  
  ! writing the bulk speed components in a file
  diag_split_name = 'Vel3'
  global_attribute = "Vbulk"
  var_name1 = "Ux"
  var_name2 = "Uy"
  var_name3 = "Uz"
  
  call create_file_diag(diag_name_file,run_name,diag_split_name,iproc-1)

  write(*,*) diag_name_file
  call wrt_field_diag(diag_name_file,uxa_proc/dna_proc,uya_proc/dna_proc,uza_proc/dna_proc, &
  trim(var_name1),trim(var_name2),trim(var_name3),& 
  nb_procs,nc_tot,ncm_tot,ncm,nxyzm,dt,t,iter, &
  nptot_proc,ns,coord_proc,voisin_proc,rang,& 
  ndims,nb_voisin,s_min,s_max,s_min_loc,s_max_loc,Spe)
  
  
! writing the electric field components in a file
  diag_split_name = 'Ele3'
  global_attribute = "Efield"
  var_name1 = "Ex"
  var_name2 = "Ey"
  var_name3 = "Ez"
  
  call create_file_diag(diag_name_file,run_name,diag_split_name,iproc-1)

  write(*,*) diag_name_file
  call wrt_field_diag(diag_name_file,ex_proc,ey_proc,ez_proc, &
  trim(var_name1),trim(var_name2),trim(var_name3),& 
  nb_procs,nc_tot,ncm_tot,ncm,nxyzm,dt,t,iter, &
  nptot_proc,ns,coord_proc,voisin_proc,rang,& 
  ndims,nb_voisin,s_min,s_max,s_min_loc,s_max_loc,Spe)  
  

  ! writing the total charge density components in a file
  diag_split_name = 'Den3'
  global_attribute = "Density"
  var_name1 = "Density"
  
  call create_file_diag(diag_name_file,run_name,diag_split_name,iproc-1)

  write(*,*) diag_name_file
  call wrt_field_diag_1_array(diag_name_file,dna_proc, &
  trim(var_name1),& 
  nb_procs,nc_tot,ncm_tot,ncm,nxyzm,dt,t,iter, &
  nptot_proc,ns,coord_proc,voisin_proc,rang,& 
  ndims,nb_voisin,s_min,s_max,s_min_loc,s_max_loc,Spe)

  ! writing the electron pressure in a file
  diag_split_name = 'Pre3'
  global_attribute = "Pressure"
  var_name1 = "Pe_tot"
  
  call create_file_diag(diag_name_file,run_name,diag_split_name,iproc-1)

  write(*,*) diag_name_file
  call wrt_field_diag_1_array(diag_name_file,pe_proc, &
  trim(var_name1),& 
  nb_procs,nc_tot,ncm_tot,ncm,nxyzm,dt,t,iter, &
  nptot_proc,ns,coord_proc,voisin_proc,rang,& 
  ndims,nb_voisin,s_min,s_max,s_min_loc,s_max_loc,Spe)
  
  deallocate(bx_proc,by_proc,bz_proc)
  deallocate(ex_proc,ey_proc,ez_proc)
  deallocate(dn_proc,dna_proc,pe_proc)
  deallocate(ux_proc,uy_proc,uz_proc)
  deallocate(uxa_proc,uya_proc,uza_proc)
  enddo

  
  !nptot = sum(nptot_proc)
  
  end subroutine split_cdf
  
 
  !!=============================================================
  !!subroutine: m_split_reassemble_cdf/wrt_field_diag
  !! FUNCTION 
  !!  Write information about Fields on CDF files for any proc
  !!
  !! OUTPUT
  !!  Only write
  subroutine wrt_field_diag(filwrt,array1,array2,array3,var1,var2,var3,nproc,nc_tot, &
                            ncm_tot,ncm,nxyzm,dt,t,iter,nptot,ns, &
                            coord_proc,voisin_proc,rang,ndims,nb_voisin,s_min,s_max, &
                            s_min_loc,s_max_loc,Spe)
 
   character(len=*),intent(in) :: filwrt
   real(dp),dimension(:,:,:),intent(in) :: array1,array2,array3
   character(len=*),intent(in) :: var1,var2,var3
   integer,intent(in) :: nproc,nxyzm,iter,nptot,ns
   integer,intent(in) :: nc_tot(3),ncm_tot(3),ncm(3)
   real(dp),intent(in) :: s_min(3),s_max(3),s_min_loc(3),s_max_loc(3),dt,t
   type(species_type) :: Spe
   
   
   integer,intent(in) :: rang,ndims,nb_voisin
   integer,allocatable :: dims(:),dimspec(:)
   !integer :: nptot_proc
   integer,intent(in) :: voisin_proc(:)
   integer,intent(in) :: coord_proc(:)
   
   integer :: ncid, stId,ii
   integer :: dimid(12), varid(75)
   character(len=40) :: name_file 
   real(dp) :: dt_t(2)
   
 
  
   !--Create the file
   stId = nf90_create(trim(filwrt) , nf90_clobber, ncid)
   call test_cdf(stId)

     
   !--Define the dimensions that will define the size of the file.
   call create_wrt_dimensions_cdf(ncid,dimid,ncm)
   
   
   !--Set the global attributes
   call set_global_attribute_cdf(ncid,"Fields")
   
   
    ii=1  
  !--Set number of cells in X,Y,Z
  stId = nf90_def_var(ncid,"nc_tot",      nf90_int,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii + 1  
  !--Get number of point for fields in Y,Z
  stId = nf90_def_var(ncid,"ncm_tot",   nf90_int,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii + 1  
  !--Set number point per processus
  stId = nf90_def_var(ncid,"ncm",   nf90_int,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii + 1  
  !--Set total number of points
  stId = nf90_def_var(ncid,"nxyzm",   nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1
  !-- Set time step
  stId = nf90_def_var(ncid,"dt_t",     QP_NF90_DP,dimid(2),varid(ii))
  call test_cdf(stId);  ii = ii + 1
  !--Set spatial step
  stId = nf90_def_var(ncid,"gstep",   QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii + 1
  !-- Set iteration
  stId = nf90_def_var(ncid,"iter",   nf90_int,dimid(1),varid(ii))
  call test_cdf(stId) ; ii = ii + 1
  !-- Set total number of particle for the opened file
  stId = nf90_def_var(ncid,"nptot",  nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1
  !-- Set number of incoming plasma species
  !stId = nf90_def_var(ncid,"ns",  nf90_int,dimid(1),varid(ii))
  !call test_cdf(stId); ii = ii + 1
  
  !--Define MPI_INFO
  !call mpiinfo_def_var_cdf(ncid,varid,dimid,ii)
  stId = nf90_def_var(ncid,"nproc",     nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1
 
  !stId = nf90_def_var(ncid, "mpiinfo_dims"  , nf90_int,dimid(4),varid(ii))
  !call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "mpiinfo_coord" , nf90_int,dimid(4),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "mpiinfo_voisin", nf90_int,dimid(5), varid(ii))
  call test_cdf(stId); ii = ii+1 
  !--Set the rank of the opened file 
  stId = nf90_def_var(ncid,"mpiinfo_me",  nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1  
  !--Set topology dimension
  stId = nf90_def_var(ncid,"mpi_ndims",  nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1 
  !--Set number of neighbours
  stId = nf90_def_var(ncid,"mpi_nb_voisin",  nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1

    
   !--Define other variables
   !--GRID_INFO
   stId = nf90_def_var(ncid, "s_min",    QP_NF90_DP,dimid(3),varid(ii))
   call test_cdf(stId); ii = ii+1                   
   stId = nf90_def_var(ncid, "s_max",    QP_NF90_DP,dimid(3),varid(ii))
   call test_cdf(stId); ii = ii+1 
   stId = nf90_def_var(ncid, "s_min_loc",  QP_NF90_DP,dimid(3),varid(ii))
   call test_cdf(stId); ii = ii+1                   
   stId = nf90_def_var(ncid, "s_max_loc",  QP_NF90_DP,dimid(3),varid(ii))
   call test_cdf(stId); ii = ii+1 
 
   !--FIELD_DEF
   stId = nf90_def_var(ncid, var1, QP_NF90_DP,dimid(6:8),varid(ii))
   call test_cdf(stId); ii = ii+1                                                
   stId = nf90_def_var(ncid, var2, QP_NF90_DP,dimid(6:8),varid(ii))
   call test_cdf(stId); ii = ii+1                                                
   stId = nf90_def_var(ncid, var3, QP_NF90_DP,dimid(6:8),varid(ii))
   call test_cdf(stId); ii = ii + 1      
   
     !--Define Speces dimensions
   call species_def_dim_cdf(ncid,dimspec)
   
   !--Define Speces variables
   call species_def_var_cdf(ncid,varid,dimspec,ii)
 
   !stId = nf90_put_att(ncid, varid(2), "units", "atomic units")
   !call test_cdf(stId)
 
   write(*,*) 'Switch to writing mode'
   !--Switch to write mode
   stId = nf90_enddef(ncid); call test_cdf(stId)
 
   ii = 1
  !--Write the  variables into the file
   stId = nf90_put_var(ncid, varid(ii), nc_tot)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), ncm_tot)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), ncm)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), nxyzm)
   call test_cdf(stId); ii = ii+1   
   dt_t(1) = dt
   dt_t(2) = t
   stId = nf90_put_var(ncid, varid(ii), dt_t)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), gstep)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), iter)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), nptot)
   call test_cdf(stId); ii = ii+1    
   !stId = nf90_put_var(ncid, varid(ii), ns)
   !call test_cdf(stId); ii = ii+1   
   stId = nf90_put_var(ncid, varid(ii), nproc)
   call test_cdf(stId); ii = ii+1    
   !stId = nf90_put_var(ncid, varid(ii), mpiinfo_dims)
   !call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), coord_proc)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid,varid(ii),voisin_proc)
   call test_cdf(stId); ii = ii + 1
   stId = nf90_put_var(ncid,varid(ii),rang)
   call test_cdf(stId); ii = ii + 1
   stId = nf90_put_var(ncid, varid(ii), ndims)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), nb_voisin)
   call test_cdf(stId); ii = ii+1       
   stId = nf90_put_var(ncid, varid(ii), s_min)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), s_max)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), s_min_loc)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), s_max_loc)
   call test_cdf(stId); ii = ii+1    
   
   stId = nf90_put_var(ncid, varid(ii), array1)
   call test_cdf(stId); ii = ii+1              
   stId = nf90_put_var(ncid, varid(ii), array2)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), array3)
   call test_cdf(stId); ii = ii+1
   
  !--Write Species infos into the file
  call species_put_var_cdf(Spe,ncid,varid,dimspec,ii)
 
   !--Close the file
   stId = nf90_close(ncid); call test_cdf(stId)
 

   end subroutine wrt_field_diag

  !!=============================================================
  !!subroutine: m_split_reassemble_cdf/wrt_field_diag_1_array
  !! FUNCTION 
  !!  Write information about Fields on CDF files for any proc
  !!
  !! OUTPUT
  !!  Only write
  subroutine wrt_field_diag_1_array(filwrt,array1,var1,nproc,nc_tot, &
                            ncm_tot,ncm,nxyzm,dt,t,iter,nptot,ns, &
                            coord_proc,voisin_proc,rang,ndims,nb_voisin,s_min,s_max, &
                            s_min_loc,s_max_loc,Spe)
 
   character(len=*),intent(in) :: filwrt
   real(dp),dimension(:,:,:),intent(in) :: array1
   character(len=*),intent(in) :: var1
   integer,intent(in) :: nproc,nxyzm,iter,nptot,ns
   integer,intent(in) :: nc_tot(3),ncm_tot(3),ncm(3)
   real(dp),intent(in) :: s_min(3),s_max(3),s_min_loc(3),s_max_loc(3),dt,t
   type(species_type) :: Spe
   
   
   integer,intent(in) :: rang,ndims,nb_voisin
   integer,allocatable :: dims(:),dimspec(:)
   !integer :: nptot_proc
   integer,intent(in) :: voisin_proc(:)
   integer,intent(in) :: coord_proc(:)
   
   integer :: ncid, stId,ii
   integer :: dimid(12), varid(75)
   character(len=40) :: name_file 
   real(dp) :: dt_t(2)
   
 
  
   !--Create the file
   stId = nf90_create(trim(filwrt) , nf90_clobber, ncid)
   call test_cdf(stId)

     
   !--Define the dimensions that will define the size of the file.
   call create_wrt_dimensions_cdf(ncid,dimid,ncm)
   
   
   !--Set the global attributes
   call set_global_attribute_cdf(ncid,"Fields")
   
   
    ii=1  
  !--Set number of cells in X,Y,Z
  stId = nf90_def_var(ncid,"nc_tot",      nf90_int,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii + 1  
  !--Get number of point for fields in Y,Z
  stId = nf90_def_var(ncid,"ncm_tot",   nf90_int,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii + 1  
  !--Set number point per processus
  stId = nf90_def_var(ncid,"ncm",   nf90_int,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii + 1  
  !--Set total number of points
  stId = nf90_def_var(ncid,"nxyzm",   nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1
  !-- Set time step
  stId = nf90_def_var(ncid,"dt_t",     QP_NF90_DP,dimid(2),varid(ii))
  call test_cdf(stId);  ii = ii + 1
  !--Set spatial step
  stId = nf90_def_var(ncid,"gstep",   QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii + 1
  !-- Set iteration
  stId = nf90_def_var(ncid,"iter",   nf90_int,dimid(1),varid(ii))
  call test_cdf(stId) ; ii = ii + 1
  !-- Set total number of particle for the opened file
  stId = nf90_def_var(ncid,"nptot",  nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1
  !-- Set number of incoming plasma species
  !stId = nf90_def_var(ncid,"ns",  nf90_int,dimid(1),varid(ii))
  !call test_cdf(stId); ii = ii + 1
  
  !--Define MPI_INFO
  !call mpiinfo_def_var_cdf(ncid,varid,dimid,ii)
  stId = nf90_def_var(ncid,"nproc",     nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1
 
  !stId = nf90_def_var(ncid, "mpiinfo_dims"  , nf90_int,dimid(4),varid(ii))
  !call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "mpiinfo_coord" , nf90_int,dimid(4),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "mpiinfo_voisin", nf90_int,dimid(5), varid(ii))
  call test_cdf(stId); ii = ii+1 
  !--Set the rank of the opened file 
  stId = nf90_def_var(ncid,"mpiinfo_me",  nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1  
  !--Set topology dimension
  stId = nf90_def_var(ncid,"mpi_ndims",  nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1 
  !--Set number of neighbours
  stId = nf90_def_var(ncid,"mpi_nb_voisin",  nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii + 1

    
   !--Define other variables
   !--GRID_INFO
   stId = nf90_def_var(ncid, "s_min",    QP_NF90_DP,dimid(3),varid(ii))
   call test_cdf(stId); ii = ii+1                   
   stId = nf90_def_var(ncid, "s_max",    QP_NF90_DP,dimid(3),varid(ii))
   call test_cdf(stId); ii = ii+1 
   stId = nf90_def_var(ncid, "s_min_loc",  QP_NF90_DP,dimid(3),varid(ii))
   call test_cdf(stId); ii = ii+1                   
   stId = nf90_def_var(ncid, "s_max_loc",  QP_NF90_DP,dimid(3),varid(ii))
   call test_cdf(stId); ii = ii+1 
 
   !--FIELD_DEF
   stId = nf90_def_var(ncid, var1, QP_NF90_DP,dimid(6:8),varid(ii))
   call test_cdf(stId); ii = ii+1                                                
     
   
     !--Define Speces dimensions
   call species_def_dim_cdf(ncid,dimspec)
   
   !--Define Speces variables
   call species_def_var_cdf(ncid,varid,dimspec,ii)
 
   !stId = nf90_put_att(ncid, varid(2), "units", "atomic units")
   !call test_cdf(stId)
 
   write(*,*) 'Switch to writing mode'
   !--Switch to write mode
   stId = nf90_enddef(ncid); call test_cdf(stId)
 
   ii = 1
  !--Write the  variables into the file
   stId = nf90_put_var(ncid, varid(ii), nc_tot)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), ncm_tot)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), ncm)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), nxyzm)
   call test_cdf(stId); ii = ii+1   
   dt_t(1) = dt
   dt_t(2) = t
   stId = nf90_put_var(ncid, varid(ii), dt_t)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), gstep)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), iter)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), nptot)
   call test_cdf(stId); ii = ii+1    
   !stId = nf90_put_var(ncid, varid(ii), ns)
   !call test_cdf(stId); ii = ii+1   
   stId = nf90_put_var(ncid, varid(ii), nproc)
   call test_cdf(stId); ii = ii+1    
   !stId = nf90_put_var(ncid, varid(ii), mpiinfo_dims)
   !call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), coord_proc)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid,varid(ii),voisin_proc)
   call test_cdf(stId); ii = ii + 1
   stId = nf90_put_var(ncid,varid(ii),rang)
   call test_cdf(stId); ii = ii + 1
   stId = nf90_put_var(ncid, varid(ii), ndims)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), nb_voisin)
   call test_cdf(stId); ii = ii+1       
   stId = nf90_put_var(ncid, varid(ii), s_min)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), s_max)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), s_min_loc)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), s_max_loc)
   call test_cdf(stId); ii = ii+1    
   
   stId = nf90_put_var(ncid, varid(ii), array1)
   call test_cdf(stId); ii = ii+1              

   
  !--Write Species infos into the file
  call species_put_var_cdf(Spe,ncid,varid,dimspec,ii)
 
   !--Close the file
   stId = nf90_close(ncid); call test_cdf(stId)
 

   end subroutine wrt_field_diag_1_array


  !!############################################################
  
   !!=============================================================
   !!subroutine: diag/create_file_diag
   !! FUNCTION 
   !!  Create the name of the file containing fields information
   !! INPUT
   !!  filwrt=suffix containgt data
   !!  me=processus identity
   !! OUTPUT
   !!  name_file=name of the file where fields will be recorded
   subroutine create_file_diag(name_file,filwrt,prefix,me)
  
    integer,intent(in) :: me
    character(len=40),intent(out) :: name_file 
    character(len=*),intent(in) :: filwrt,prefix 
  
    write(name_file,'(a4,a1,i4.4,a1,a)')trim(prefix),"_",me,'_',trim(filwrt)
  
#ifdef HAVE_NETCDF
    name_file = trim(name_file)//".nc"
#endif
  
  
 end subroutine create_file_diag


#endif

end module m_split_cdf
