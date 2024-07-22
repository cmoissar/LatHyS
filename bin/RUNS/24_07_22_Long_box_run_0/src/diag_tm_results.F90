!!=============================================================
!!=============================================================
module diag_tm_results

 use defs_basis
 use defs_mpitype
 use defs_species
 use defs_diag_type,only   : diag_type
 use defs_parametre,only   : idisp,iresis,nsub,psi,phi,resis,nhm
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
      create_file_name    !--create the file name for time_result record
 
 public  ::              &
      wrt_tm_results      !--time_result record

contains
 !!############################################################

 !!=============================================================
 !!subroutine: diag_tm_results/create_file_name
 !! FUNCTION 
 !!  Create the name of the file containing time results information
 !! INPUT
 !!  filwrt=suffix containgt data
 !! OUTPUT
 !!  name_file=name of the file where time results will be recorded
 subroutine create_file_name(name_file,filwrt)

  character(len=40),intent(out) :: name_file 
  character(len=*),intent(in) :: filwrt 
    
  write(name_file,'(a3,a)')"b3_",trim(filwrt)
#ifdef HAVE_NETCDF
  name_file = trim(name_file)//".nc"
#endif

 end subroutine create_file_name

#ifdef HAVE_NETCDF 
 !!=============================================================
 !!subroutine: diag_tm_results/wrt_tm_results
 !! FUNCTION 
 !!  Write Time Results on CDF files for processus 0
 !!
 !! OUTPUT
 !!  Only write
 subroutine wrt_tm_results(filwrt)

  character(len=*),intent(in) :: filwrt

  integer :: ncid, stId,ii
  integer :: dimid(12), varid(80)
  integer,allocatable :: dimspec(:)
  character(len=40) :: name_file 

  __WRT_DEBUG_IN("wrt_tm_results")
  __GETTIME(10,1)  !--Timer start

  !--Select Proc 0
  if(mpiinfo%me/=0) return

  !--Create the file name 
  call create_file_name(name_file,filwrt)

  stId = nf90_create(name_file , nf90_clobber, ncid)
  call test_cdf(stId)

  !--Define the dimensions that will define the size of the file.
  call create_wrt_dimensions_cdf(ncid,dimid,ncm)

  !--Define Speces dimensions
  call species_def_dim_cdf(ncid,dimspec)

  !--Global attributes.
  call set_global_attribute_cdf(ncid,"Time results")

  ii = 1
  !--Define commons variables
  call common_def_var_cdf(ncid,varid,dimid,ii)

  !--Define other variables
  stId = nf90_def_var(ncid, "nc_tot",            nf90_int,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "s_min",           QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1             
  stId = nf90_def_var(ncid, "s_max",           QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1            
  stId = nf90_def_var(ncid, "idisp",            nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1            
  stId = nf90_def_var(ncid, "iresis",           nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1            
  stId = nf90_def_var(ncid, "nsub",             nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1            
  stId = nf90_def_var(ncid, "psi",            QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1            
  stId = nf90_def_var(ncid, "phi",            QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1            
  stId = nf90_def_var(ncid, "b0",             QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1            
  stId = nf90_def_var(ncid, "resis",          QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1           
  stId = nf90_def_var(ncid, "nhm",              nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1           
  stId = nf90_def_var(ncid, "dinfo_sv",    QP_NF90_DP,(/dimid(3), dimid(11)/),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_B_mean",QP_NF90_DP,(/dimid(3), dimid(11)/),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_B_quad",QP_NF90_DP,(/dimid(3), dimid(11)/),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_v_mean",QP_NF90_DP,(/dimid(3), dimid(11)/),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_v_quad",QP_NF90_DP,(/dimid(3), dimid(11)/),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_etot",    QP_NF90_DP,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1             
  stId = nf90_def_var(ncid, "dinfo_etherm",  QP_NF90_DP,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1             
  stId = nf90_def_var(ncid, "dinfo_mag"  ,   QP_NF90_DP,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_mag_mean",QP_NF90_DP,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_kin  ",   QP_NF90_DP,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1             
  stId = nf90_def_var(ncid, "dinfo_kins ",   QP_NF90_DP,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_kin_mean",QP_NF90_DP,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_dn_mean", QP_NF90_DP,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_n_part",  nf90_int,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_n_out  ", nf90_int,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_n_out_g", nf90_int,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_n_in   ", nf90_int,dimid(11),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dinfo_n_in_spec ", nf90_int,(/dimid(11),dimspec(2)/),varid(ii))
  call test_cdf(stId); ii = ii+1

  !--Define Speces variables
  call species_def_var_cdf(ncid,varid,dimspec,ii)

  !stId = nf90_put_att(ncid, varid(2), "units", "atomic units")
  !call test_cdf(stId)

  !--Switch to write mode.
  stId = nf90_enddef(ncid); call test_cdf(stId)

  ii = 1
  !--Write common variables into the file
  call common_put_var_cdf(ncid,varid,ii)

  !--Write the other variables into the file
  stId = nf90_put_var(ncid, varid(ii), nc_tot)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), s_min)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), s_max)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), idisp)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), iresis)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), nsub)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), psi)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), phi)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), b0)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), resis)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), nhm)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), diag_info%sv)
  call test_cdf(stId); ii = ii+1    
  stId = nf90_put_var(ncid, varid(ii), diag_info%champ_bmoy)
  call test_cdf(stId); ii = ii+1    
  stId = nf90_put_var(ncid, varid(ii), diag_info%champ_bqua)
  call test_cdf(stId); ii = ii+1    
  stId = nf90_put_var(ncid, varid(ii), diag_info%vitesse_moy)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%vitesse_qua)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%ener_tot)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%ener_therm)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%ener_mag)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%ener_mag_moy)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%ener_kin)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%ener_kins)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%ener_kin_moy)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%dn_moy)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%n_part)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%n_out)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%n_out_g)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%n_in)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), diag_info%n_in_spec)
  call test_cdf(stId); ii = ii+1

  !--Write Species infos into the file
  call species_put_var_cdf(Spe,ncid,varid,dimspec,ii)

  !--Close the file
  stId = nf90_close(ncid); call test_cdf(stId)

  __GETTIME(10,1)  !--Timer stop
  __WRT_DEBUG_OUT("wrt_tm_results")
 end subroutine wrt_tm_results

#else

 !********************************* WRT_TM_RESULTS **************************************
 subroutine wrt_tm_results(filwrt)

  character(len=*),intent(in) :: filwrt

  integer :: iunit
  character(len = 40) :: name_file
  character(len = 500):: msg

  __WRT_DEBUG_IN("wrt_tm_results")
  __GETTIME(10,1)  !--Timer start
  !--Bilans en fonction du temps

  iunit = 7

  !--Create the file name
  call create_file_name(name_file,filwrt)

  !--Select processus 0
  if(mpiinfo%me == 0) then
   write(msg,'(a,a)')' ..Creating file ',trim(name_file)
   call wrt_double(6,msg,wrtscreen,wrtdisk)

   open(UNIT   = iunit,&
        FILE   = name_file,&
        FORM   = 'unformatted',&
        ACTION = 'write', &
        STATUS = 'unknown')

   write(iunit) nc_tot,gstep,s_min,s_max,iter,dt,t,nptot
   
   !--Write Species infos
   call species_put_var_bin(Spe,iunit)

   write(iunit) nhm
   write(iunit) diag_info%sv(1,:),diag_info%sv(2,:),diag_info%sv(3,:)
   write(iunit) diag_info%ener_kins,diag_info%ener_kin,diag_info%ener_mag,diag_info%ener_tot,diag_info%ener_therm
   write(iunit) diag_info%ener_mag_moy,diag_info%ener_kin_moy
   write(iunit) diag_info%champ_bmoy(1,:),diag_info%champ_bmoy(2,:),diag_info%champ_bmoy(3,:)
   write(iunit) diag_info%champ_bqua(1,:),diag_info%champ_bqua(2,:),diag_info%champ_bqua(3,:)
   write(iunit) diag_info%vitesse_moy(1,:),diag_info%vitesse_moy(2,:),diag_info%vitesse_moy(3,:)
   write(iunit) diag_info%vitesse_qua(1,:),diag_info%vitesse_qua(2,:),diag_info%vitesse_qua(3,:)
   write(iunit) psi,phi,b0
   write(iunit) idisp,iresis,resis
   write(iunit) nsub
   write(iunit) diag_info%n_part,diag_info%n_in,diag_info%n_out,diag_info%n_out_g
   write(iunit) diag_info%n_in_spec

   close(iunit)
  endif

  __GETTIME(10,2)!--Timer stop
  __WRT_DEBUG_OUT("wrt_tm_results")
 end subroutine wrt_tm_results
 !*********************************** END WRT_TM_RESULTS ********************************

#endif

end module diag_tm_results
