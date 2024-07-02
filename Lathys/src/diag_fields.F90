!!=============================================================
!!=============================================================
module diag_fields

 use defs_basis
 use defs_mpitype
 use defs_species
 use defs_diag_type,only   : diag_type
 use defs_variable
 use defs_grid
 use environment
 use m_timing,only         : time_get
 use m_writeout
#ifdef HAVE_NETCDF 
 use defs_basic_cdf
 use diag_wrt_common_cdf
 use netcdf
#endif

#include "q-p_common.h"

 implicit none
 private

 private ::              &
      create_file_name    !--create the file name for field record
 
 public  ::              &
      wrt_fields          !--field record

contains
 !!############################################################

 !!=============================================================
 !!subroutine: diag_fields/create_file_name
 !! FUNCTION 
 !!  Create the name of the file containing fields information
 !! INPUT
 !!  filwrt=suffix containgt data
 !!  me=processus identity
 !! OUTPUT
 !!  name_file=name of the file where fields will be recorded
 subroutine create_file_name(name_file,filwrt,me)

  integer,intent(in) :: me
  character(len=40),intent(out) :: name_file 
  character(len=*),intent(in) :: filwrt 

  write(name_file,'(a3,i4.4,a1,a)')"c3_",me,'_',trim(filwrt)
#ifdef HAVE_NETCDF
  name_file = trim(name_file)//".nc"
#endif

 end subroutine create_file_name

#ifdef HAVE_NETCDF 
 !!=============================================================
 !!subroutine: diag_fields/wrt_fields
 !! FUNCTION 
 !!  Write information about Fields on CDF files for any proc
 !!
 !! OUTPUT
 !!  Only write
 subroutine wrt_fields(filwrt)

  character(len=*),intent(in) :: filwrt

  integer :: ncid, stId,ii
  integer :: dimid(12), varid(75)
  integer,allocatable :: dimspec(:)
  character(len=40) :: name_file 

  __WRT_DEBUG_IN("wrt_fields")
  __GETTIME(8,1) !--Timer start

  !--Create the file name
  call create_file_name(name_file,filwrt,mpiinfo%me)

  !--Create the file
  stId = nf90_create(name_file , nf90_clobber, ncid)
  call test_cdf(stId)

  !--Define the dimensions that will define the size of the file.
  call create_wrt_dimensions_cdf(ncid,dimid,ncm)

  !--Define Speces dimensions
  call species_def_dim_cdf(ncid,dimspec)

  !--Set the global attributes
  call set_global_attribute_cdf(ncid,"Fields")

  !print *, 'species info ok'
  ii=1  
  !--Define common variables
  call common_def_var_cdf(ncid,varid,dimid,ii)

  !--Define MPI_INFO
  call mpiinfo_def_var_cdf(ncid,varid,dimid,ii)

  !--Define Speces variables
  call species_def_var_cdf(ncid,varid,dimspec,ii)

  !--Define other variables
  !--GRID_INFO
  stId = nf90_def_var(ncid, "nxyzm",      nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "ncm_tot",    nf90_int,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "nc_tot ",    nf90_int,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "ncm"    ,    nf90_int,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "s_min",    QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1                   
  stId = nf90_def_var(ncid, "s_max",    QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "s_min_loc",  QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1                   
  stId = nf90_def_var(ncid, "s_max_loc",  QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1 

  !--FIELD_DEF
  stId = nf90_def_var(ncid, "Bfield_x", QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "Bfield_y", QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "Bfield_z", QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "Efield_x", QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "Efield_y", QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "Efield_z", QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "vel_x",    QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                             
  stId = nf90_def_var(ncid, "vel_y",    QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                             
  stId = nf90_def_var(ncid, "vel_z",    QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                                 
  stId = nf90_def_var(ncid, "vela_x",   QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                              
  stId = nf90_def_var(ncid, "vela_y",   QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "vela_z",   QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "dn"   ,    QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "dna"  ,    QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1   
  stId = nf90_def_var(ncid, "Pe"  , QP_NF90_DP,dimid(6:8),varid(ii))
!!  call test_cdf(stId); ii = ii+1   
!!  stId = nf90_def_var(ncid, "Prod_O"  , QP_NF90_DP,dimid(6:8),varid(ii))
!!  call test_cdf(stId); ii = ii+1   
!!  stId = nf90_def_var(ncid, "Prod_CO2", QP_NF90_DP,dimid(6:8),varid(ii))
  call test_cdf(stId);  

  !stId = nf90_put_att(ncid, varid(2), "units", "atomic units")
  !call test_cdf(stId)

  !--Switch to write mode
  stId = nf90_enddef(ncid); call test_cdf(stId)
  !print * ,'switching to writing mode'

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
  stId = nf90_put_var(ncid, varid(ii), ncm )
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), s_min)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), s_max)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), s_min_loc)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), s_max_loc)
  call test_cdf(stId); ii = ii+1              
  ! stId = nf90_put_var(ncid, varid(ii), nfl)
  ! call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), Bfield%x)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), Bfield%y)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), Bfield%z)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), Efield%x)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), Efield%y)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), Efield%z)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), vel%x)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), vel%y)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), vel%z)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), vela%x)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), vela%y)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), vela%z)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), dn)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), dna)   
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), pe)   
!!  call test_cdf(stId); ii = ii+1 
!!  stId = nf90_put_var(ncid, varid(ii), density_Hp)   
!!  call test_cdf(stId); ii = ii+1 
!!  stId = nf90_put_var(ncid, varid(ii), density_Hp)   
  call test_cdf(stId)

  !--Close the file
  stId = nf90_close(ncid); call test_cdf(stId)

  __GETTIME(8,2) !--Timer stop
  __WRT_DEBUG_OUT("wrt_fields")
 end subroutine wrt_fields

#else

 !******************************* WRT_FIELDS ***************************
 subroutine wrt_fields(filwrt)
 use defs_parametre
  character(len=*),intent(in) :: filwrt

  character(len=40)  :: name_file
  character(len=500) :: msg

  __WRT_DEBUG_IN("wrt_fields")
  __GETTIME(8,1) !--Timer start
  !--Fields

  !--Create the file name
  call create_file_name(name_file,filwrt,mpiinfo%me)

  write(msg,'(a,i4.4,a,a)')&
       &' ..Writing process ',mpiinfo%me,&
       &' on file ',name_file
  call wrt_double(6,msg,wrtscreen,0)

  open(UNIT   = 1, &
       FILE   = name_file, &
       ACTION = 'write', &
       FORM   = 'unformatted', &
       STATUS = 'unknown', &
       ACCESS = 'sequential')

  !--On écrit en tete du fichier les caractéristiques de la simulation
  write(unit = 1) mpiinfo%me   
  write(unit = 1) nxyzm
  write(unit = 1) nproc,nb_voisins,ncm,ndims
  write(unit = 1) mpiinfo%dims 
  write(unit = 1) mpiinfo%coord
  write(unit = 1) ncm_tot
  write(unit = 1) nc_tot,gstep,s_min,s_max,s_min_loc,s_max_loc
  write(unit = 1) iter,dt,t,nptot,ns,nfl
  write(unit = 1) mpiinfo%voisin

  !--Write Species informations
  call species_put_var_bin(Spe,unit = 1)

  !--Ecriture des champs
  write(unit = 1) &
       Bfield%x(:,:,:), &
       Bfield%y(:,:,:), &
       Bfield%z(:,:,:)
  write(unit = 1) &
       Efield%x(:,:,:), &
       Efield%y(:,:,:), &
       Efield%z(:,:,:)
  write(unit = 1) dn(:,:,:), &
       vel%x(:,:,:), &
       vel%y(:,:,:), &
       vel%z(:,:,:)
  write(unit = 1) dna(:,:,:), &
       vela%x(:,:,:), &
       vela%y(:,:,:), &
       vela%z(:,:,:)
  close(unit = 1)

  write(msg,'(a,i4.4,a,a)')&
       &' ..Process ',mpiinfo%me,&
       &' written on file ',name_file
  call wrt_double(6,msg,wrtscreen,0)

  __GETTIME(8,2) !--Timer stop
  __WRT_DEBUG_OUT("wrt_fields")
 end subroutine wrt_fields
 !******************************* END WRT_FIELDS ***********************

#endif

end module diag_fields
