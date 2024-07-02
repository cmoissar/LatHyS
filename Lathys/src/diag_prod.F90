!!=============================================================
!!=============================================================
module diag_prod

 use defs_basis
 use defs_mpitype
 use defs_species
 use defs_diag_type,only   : diag_type
 use defs_variable
 use defs_grid
 use environment
 use atm_ionosphere
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
      create_file_name_io    !--create the file name for field record
 
 public  ::              &
      wrt_prod          !--field record

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
 subroutine create_file_name_io(name_file,filwrt,me)

  integer,intent(in) :: me
  character(len=40),intent(out) :: name_file 
  character(len=*),intent(in) :: filwrt 

  write(name_file,'(a5,i4.4,a1,a)')"Pro3_",me,'_',trim(filwrt)
#ifdef HAVE_NETCDF
  name_file = trim(name_file)//".nc"
#endif

 end subroutine create_file_name_io

#ifdef HAVE_NETCDF 
 !!=============================================================
 !!subroutine: diag_fields/wrt_fields
 !! FUNCTION 
 !!  Write information about Fields on CDF files for any proc
 !!
 !! OUTPUT
 !!  Only write
 subroutine wrt_prod(filwrt)
 use environment,only : atmosphere

  character(len=*),intent(in) :: filwrt

  integer :: ncid, stId,ii,i, compt
  integer :: dimid(12), varid(75)
  integer,allocatable :: dimspec(:)
  character(len=40) :: name_file
  character(len=11) :: nam

  __WRT_DEBUG_IN("wrt_iono")
  __GETTIME(8,1) !--Timer start

  !--Create the file name
  call create_file_name_io(name_file,filwrt,mpiinfo%me)

  !--Create the file
  stId = nf90_create(name_file , nf90_clobber, ncid)
  call test_cdf(stId)

  !--Define the dimensions that will define the size of the file.
  call create_wrt_dimensions_cdf(ncid,dimid,ncm)

  !--Define Speces dimensions
  call species_def_dim_cdf(ncid,dimspec)

  !--Set the global attributes
  call set_global_attribute_cdf(ncid,"Production")

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
  stId = nf90_def_var(ncid, "n_species",  QP_NF90_DP,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "n_spe_pp",  QP_NF90_DP,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "n_reac_pp",  QP_NF90_DP,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "n_reac_exc",  QP_NF90_DP,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1 

  !--prod_DEF
  compt = 0
  do i=1,atmosphere%n_species
  ! on regarde s'il s'agit d'une espece ionospherique
    if (atmosphere%species(i)%iono .EQV. .TRUE.) then
      compt = compt+1
      write(nam,'(a,i2)')"Spe_",compt
      stId = nf90_def_var(ncid,nam,nf90_char,dimid(12),varid(ii))
      call test_cdf(stId); ii = ii+1
      write(nam,'(2a)')"Prod_",TRIM(atmosphere%species(i)%name)
      stId = nf90_def_var(ncid, nam , QP_NF90_DP,dimid(6:8),varid(ii))
      call test_cdf(stId); ii = ii+1
    endif
  enddo                                            
  ii=ii-1


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
  stId = nf90_put_var(ncid, varid(ii), atmosphere%n_species)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), atmosphere%n_spe_pp)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), atmosphere%n_pp)
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), atmosphere%n_exc)
  call test_cdf(stId); ii = ii+1 

 compt = 0
 do i=1,atmosphere%n_species
  ! on regarde s'il s'agit d'une espece ionospherique
    if (atmosphere%species(i)%iono .EQV. .TRUE.) then
      compt = compt+1
      write(nam,'(2a)')"Prod_",TRIM(atmosphere%species(i)%name)
      stId = nf90_put_var(ncid, varid(ii), nam)
      write(*,*) nam
      call test_cdf(stId); ii = ii+1
      stId = nf90_put_var(ncid, varid(ii), atmosphere%species(i)%prod)
      call test_cdf(stId); ii = ii+1
    endif  
  enddo                                            
  ii=ii-1

  !--Close the file
  stId = nf90_close(ncid); call test_cdf(stId)

  __GETTIME(8,2) !--Timer stop
  __WRT_DEBUG_OUT("wrt_prod")
 end subroutine wrt_prod

#endif

end module diag_prod
