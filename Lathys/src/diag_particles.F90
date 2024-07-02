!!=============================================================
!!=============================================================
module diag_particles

 use defs_basis
 use defs_mpitype
 use defs_species
 use defs_diag_type,only   : diag_type
 use defs_parametre,only   : npm
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
      create_file_name    !--create the file name for particles record
 
 public  ::              &
      wrt_particles       !--particles record

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
  
  write(name_file,'(a3,i4.4,a1,a)')"p3_",me,'_',trim(filwrt)
#ifdef HAVE_NETCDF
  name_file = trim(name_file)//".nc"
#endif

 end subroutine create_file_name

#ifdef HAVE_NETCDF 
 !!=============================================================
 !!subroutine: diag_particles/wrt_particles
 !! FUNCTION 
 !!  Write information about Particles on CDF files for any proc
 !!
 !! OUTPUT
 !!  Only write
 subroutine wrt_particles(filwrt)

  character(len=*),intent(in) :: filwrt

  integer :: ncid, stId,ii
  integer :: dimid(12), varid(100)
  integer,allocatable :: dimspec(:)
  character(len=40) :: name_file 

  __WRT_DEBUG_IN("wrt_particles")
  __GETTIME(9,1) !--Timer start

  !--Create the file name
  call create_file_name(name_file,filwrt,mpiinfo%me)

  stId = nf90_create(name_file , nf90_clobber, ncid)
  call test_cdf(stId)

  !--Define the dimensions that will define the size of the file.
  call  create_wrt_dimensions_cdf(ncid,dimid,ncm)

  !--Define Speces dimensions
  call species_def_dim_cdf(ncid,dimspec)
  
  !--Set the global attributes
  call set_global_attribute_cdf(ncid,"Particles")

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
 
  !--PARTICLES_DEF
  call def_var_particle_cdf(ncid,varid,dimid,ii)
                                              
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

  !--Put particles
  call put_var_particle_cdf(particule,nptot,ncid,varid,ii)

  !--Close the file
  stId = nf90_close(ncid); call test_cdf(stId)

  __GETTIME(9,2) !--Timer stop
  __WRT_DEBUG_OUT("wrt_particles")
 end subroutine wrt_particles

#else

 subroutine wrt_particles(filwrt)

  character(len=*),intent(in) :: filwrt
  character(len = 40) :: name_file
  character(len=500) :: msg

  __WRT_DEBUG_IN("wrt_particles")
  __GETTIME(9,1) !--Timer start

  !--Particles
  !--Create the file name
  call create_file_name(name_file,filwrt,mpiinfo%me)

  write(msg,'(a,i4.4,a,a)')&
       &' ..Writing process ',mpiinfo%me,&
       &' on file ',name_file
  call wrt_double(6,msg,wrtscreen,0)

  open(UNIT   =  1, &
       FILE   =  name_file, &
       ACTION = 'write', &
       FORM   = 'unformatted', &
       STATUS = 'unknown', &
       ACCESS = 'sequential')

  !--On écrit en tete du fichier les caractéristiques de la simulation
  write(unit = 1) mpiinfo%me   
  write(unit = 1) nxyzm
  write(unit = 1) nproc,nb_voisins,ncm,ndims
  write(unit = 1) mpiinfo%dims
  write(unit = 1) nptot
  write(unit = 1) ncm_tot
  write(unit = 1) nc_tot,gstep,s_min,s_max,iter,dt,t,nptot,ns,nfl,npm
  write(unit = 1) mpiinfo%voisin

  !--Write Species infos
  call species_put_var_bin(Spe,unit = 1)

  !--Write Particles
  call put_var_particle_bin(particule,nptot,unit = 1)

  close(unit = 1)

  write(msg,'(a,i4.4,a,a)')&
       &' ..Process ',mpiinfo%me,&
       &' written on file ',name_file
  call wrt_double(6,msg,wrtscreen,0)

  __GETTIME(9,2) !--Timer stop
  __WRT_DEBUG_OUT("wrt_particles")
 end subroutine wrt_particles
 !******************************** END WRT_PARTICLES_PROCESSUS **************

#endif

end module diag_particles
