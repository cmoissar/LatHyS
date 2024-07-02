!!=============================================================
!!=============================================================
module diag_wrt_common_cdf

 use defs_basis
 use defs_parametre,only :nhm,dt,gstep,fildat
 use defs_variable
 use defs_mpitype
#ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
#endif
#include "q-p_common.h"

 implicit none
 private

#ifdef HAVE_NETCDF
 public  ::                        &
      create_wrt_dimensions_cdf,   &!--Create Dimensions Array
      set_global_attribute_cdf,    &
      common_def_var_cdf,          &
      common_put_var_cdf,          &
      mpiinfo_def_var_cdf,         &
      mpiinfo_put_var_cdf

contains
 !!############################################################


 !!=============================================================
 !!subroutine: m_wrt_common_cdf/create_wrt_dimensions_cdf
 !! FUNCTION 
 !!  Create dimid(11) array contaninif dimension Identity with respect
 !!  to ncid (cdf file, yet opened)
 !!
 !! INPUT
 !!  ncid=Identity of the cdf file
 !! OUTPUT
 !!  dimid(11)=Dimensions array
 subroutine create_wrt_dimensions_cdf(ncid,dimid,ncm)

  integer,intent(in) :: ncid
  integer,intent(in) :: ncm(3)
  integer,intent(out) :: dimid(12)

  integer :: stId

  !--Define the dimensions that will define the size of the file.
  stId = nf90_def_dim(ncid, "dim_scalar"      , 1         , dimid(1))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "bidimensional"   , 2         , dimid(2))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "space_dimension" , 3         , dimid(3))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "mpi_ndims"       , ndims     , dimid(4))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "mpi_nb_voisin"   , nb_voisins, dimid(5))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_x"          , ncm(1)    , dimid(6))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_y"          , ncm(2)    , dimid(7))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "size_z"          , ncm(3)    , dimid(8))
  call test_cdf(stId) 
  stId = nf90_def_dim(ncid, "particles coords", 11        , dimid(9))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "nptot"           , abs(nptot), dimid(10))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "number of diag. ", nhm       , dimid(11))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "string size ", 14       , dimid(12))
  call test_cdf(stId)

 end subroutine create_wrt_dimensions_cdf


 !!=============================================================
 !!subroutine: m_wrt_common_cdf/set_global_attribute_cdf
 !! FUNCTION 
 !!  Set the global attribute in a NetCDF file
 !!
 !! INPUT
 !!  ncid=Identity of the cdf file
 !!  title_wrt=Title(Metadata) of the file
 subroutine set_global_attribute_cdf(ncid,title_wrt)

  integer,intent(in) :: ncid
  character(len=*),intent(in) :: title_wrt

  integer :: stId

  !--Set the global attributes
  stId = nf90_put_att(ncid, nf90_global, "file_format", "NetCDF")
  call test_cdf(stId)
  stId = nf90_put_att(ncid, nf90_global, "file_format_version", 1.4)
  call test_cdf(stId)
  stId = nf90_put_att(ncid, nf90_global, "Conventions", "http://www.unidata.ucar.edu/software/netcdf/")
  call test_cdf(stId)
  stId = nf90_put_att(ncid, nf90_global, "Program", "Quiet Plasma")
  call test_cdf(stId)
  stId = nf90_put_att(ncid, nf90_global, "Title", title_wrt)
  call test_cdf(stId)
  stId = nf90_put_att(ncid, nf90_global, "Date", fildat)
  call test_cdf(stId)

 end subroutine set_global_attribute_cdf


 !!=============================================================
 !!subroutine: m_wrt_common_cdf/common_def_var_cdf
 !! FUNCTION 
 !!
 !! INPUT
 subroutine common_def_var_cdf(ncid,varid,dimid,ii)

  integer,intent(in)   :: ncid
  integer,intent(inout):: ii
  integer,intent(inout),dimension(:) :: varid,dimid

  integer :: stId

  stId = nf90_def_var(ncid, "gstep",QP_NF90_DP,dimid(3), varid(ii))
  call test_cdf(stId); ii = ii+1                     
  stId = nf90_def_var(ncid, "iter",   nf90_int,dimid(1), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "dt_t", QP_NF90_DP,dimid(2), varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "nptot",  nf90_int,dimid(1), varid(ii))
  call test_cdf(stId); ii = ii+1 

 end subroutine common_def_var_cdf


 !!=============================================================
 !!subroutine: m_wrt_common_cdf/common_put_var_cdf
 !! FUNCTION 
 !!
 !! INPUT
 subroutine common_put_var_cdf(ncid,varid,ii)

  integer,intent(in)   :: ncid
  integer,intent(inout):: ii
  integer,intent(in),dimension(:) :: varid

  integer :: stId

  stId = nf90_put_var(ncid, varid(ii), gstep)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), iter)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), (/dt,t/))
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), nptot)
  call test_cdf(stId); ii = ii+1

 end subroutine common_put_var_cdf

 !!=============================================================
 !!subroutine: m_wrt_common_cdf/mpiinfo_def_var_cdf
 !! FUNCTION 
 !!
 !! INPUT
 subroutine mpiinfo_def_var_cdf(ncid,varid,dimid,ii)

  integer,intent(in)   :: ncid
  integer,intent(inout):: ii
  integer,intent(inout),dimension(:) :: varid,dimid

  integer :: stId

  stId = nf90_def_var(ncid, "nproc"         , nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "mpiinfo_me"    , nf90_int,dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "mpiinfo_dims"  , nf90_int,dimid(4),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "mpiinfo_coord" , nf90_int,dimid(4),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "mpiinfo_voisin", nf90_int,dimid(5), varid(ii))
  call test_cdf(stId); ii = ii+1 

 end subroutine mpiinfo_def_var_cdf

 !!=============================================================
 !!subroutine: m_wrt_common_cdf/mpiinfo_put_var_cdf
 !! FUNCTION 
 !!
 !! INPUT
 subroutine mpiinfo_put_var_cdf(ncid,varid,ii)

  integer,intent(in)   :: ncid
  integer,intent(inout):: ii
  integer,intent(in),dimension(:) :: varid

  integer :: stId

  stId = nf90_put_var(ncid, varid(ii), nproc)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), mpiinfo%me)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), mpiinfo%dims)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), mpiinfo%coord)
  call test_cdf(stId); ii = ii+1 
  stId = nf90_put_var(ncid, varid(ii), mpiinfo%voisin)
  call test_cdf(stId); ii = ii+1 

 end subroutine mpiinfo_put_var_cdf

#endif
end module diag_wrt_common_cdf
