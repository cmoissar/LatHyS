!!=============================================================
!!=============================================================
module defs_basic_cdf

 use defs_basis
#ifdef HAVE_NETCDF 
 use netcdf
#endif
#include "q-p_common.h"

 implicit none
 private

#ifdef HAVE_NETCDF
 public  ::                        &
      test_cdf,                    &
      get_simple_dimens_cdf,       &
      get_simple_variable_cdf

 interface get_simple_variable_cdf
  module procedure get_simple_variable_cdf_int
  module procedure get_simple_variable_cdf_int2_kind1
  module procedure get_simple_variable_cdf_int2_kind2
  module procedure get_simple_variable_cdf_int2
  module procedure get_simple_variable_cdf_dp
  module procedure get_simple_variable_cdf_dparr
  module procedure get_simple_variable_cdf_dparr2
  module procedure get_simple_variable_cdf_dparr3
  module procedure get_simple_variable_cdf_dparr4
  module procedure get_simple_variable_cdf_dparr6
  module procedure get_simple_variable_cdf_string
  module procedure get_simple_variable_cdf_string2
 end interface

contains
 !!############################################################

!  !!=============================================================
!  !!function: m_basic_cdf/get_version_cdf
!  !! FUNCTION 
!  !!  Get the NetCDF version if used 
!  !!
!  !! OUTPUT
!  !!  cdf_version (string)
!  subroutine get_version_cdf()
! #ifdef HAVE_NETCDF
!   print *,"Compiled with NetCDF : ",trim(nf90_inq_libvers())
! #else 
!   print *,"NetCDF not used: Unformatted Output "
! #endif
!  end subroutine get_version_cdf

 !!=============================================================
 !!subroutine: m_basic_cdf/test_cdf
 !! FUNCTION 
 !!  Test the error flag of netcdf function
 !!
 !! OUTPUT
 !!  Only write
 subroutine test_cdf(stId)

  integer,intent(in) :: stId

  if (stId /= nf90_noerr) then
   write(6,*)  nf90_strerror(stId)
   stop
  endif

 end subroutine test_cdf

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_dimens_cdf
 !! FUNCTION 
 !!  Get dimension from an open Netcdf file
 !!
 !! OUTPUT
 subroutine get_simple_dimens_cdf(ncid,namedim,dimens)

  integer,intent(in)  :: ncid
  integer,intent(out) :: dimens
  character(len=*),intent(in)  :: namedim
  
  integer :: stId,dimid

  stId = nf90_inq_dimid(ncid, namedim, dimid)
  call test_cdf(stId)
  stId = nf90_inquire_dimension(ncid,dimid,len=dimens) 
  call test_cdf(stId)

 end subroutine get_simple_dimens_cdf

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_int2_kind2
 !! FUNCTION 
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_int2_kind2(ncid,namevar,values)

  integer,intent(in)  :: ncid
  integer(i2b),intent(out) :: values(:)
  character(len=*),intent(in)  :: namevar
  
  integer :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_int2_kind2

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_int2_kind1
 !! FUNCTION 
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_int2_kind1(ncid,namevar,values)

  integer,intent(in)  :: ncid
  integer(i1b),intent(out) :: values(:)
  character(len=*),intent(in)  :: namevar
  
  integer :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_int2_kind1


 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_int
 !! FUNCTION 
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_int(ncid,namevar,values)

  integer,intent(in)  :: ncid
  integer,intent(out) :: values
  character(len=*),intent(in)  :: namevar
  
  integer :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_int

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_int2
 !! FUNCTION 
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_int2(ncid,namevar,values)

  integer,intent(in)  :: ncid
  integer,intent(out) :: values(:)
  character(len=*),intent(in)  :: namevar
  
  integer :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_int2

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_dp
 !! FUNCTION 
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_dp(ncid,namevar,values)

  integer,intent(in)   :: ncid
  real(dp),intent(out) :: values
  character(len=*),intent(in)  :: namevar
  
  integer :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_dp

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_int
 !! FUNCTION 
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_dparr(ncid,namevar,values)

  integer,intent(in)   :: ncid
  real(dp),intent(out) :: values(:)
  character(len=*),intent(in)  :: namevar
  
  integer :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_dparr

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_dparr2
 !! FUNCTION 
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_dparr2(ncid,namevar,values)

  integer,intent(in)   :: ncid
  real(dp),intent(out) :: values(:,:)
  character(len=*),intent(in)  :: namevar
  
  integer :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_dparr2

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_dparr3
 !! FUNCTION 
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_dparr3(ncid,namevar,values)

  integer,intent(in)   :: ncid
  real(dp),intent(out) :: values(:,:,:)
  character(len=*),intent(in)  :: namevar
  
  integer :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_dparr3
 
  !!=============================================================
  !!subroutine: m_basic_cdf/get_simple_variable_cdf_dparr4
  !! FUNCTION 
  !!  Get variable from an open Netcdf file
  !!
  !! OUTPUT
  !!  Only write
  subroutine get_simple_variable_cdf_dparr4(ncid,namevar,values)
 
   integer,intent(in)   :: ncid
   real(dp),intent(out) :: values(:,:,:,:)
   character(len=*),intent(in)  :: namevar
   
   integer :: stId,varid
 
   stId = nf90_inq_varid(ncid, namevar, varid)
   call test_cdf(stId)
   stId = nf90_get_var(ncid,varid,values)
   call test_cdf(stId)
 
 end subroutine get_simple_variable_cdf_dparr4
 
  !!=============================================================
  !!subroutine: m_basic_cdf/get_simple_variable_cdf_dparr6
  !! FUNCTION 
  !!  Get variable from an open Netcdf file
  !!
  !! OUTPUT
  !!  Only write
  subroutine get_simple_variable_cdf_dparr6(ncid,namevar,values)
 
   integer,intent(in)   :: ncid
   real(dp),intent(out) :: values(:,:,:,:,:,:)
   character(len=*),intent(in)  :: namevar
   
   integer :: stId,varid
 
   stId = nf90_inq_varid(ncid, namevar, varid)
   call test_cdf(stId)
   stId = nf90_get_var(ncid,varid,values)
   call test_cdf(stId)
 
 end subroutine get_simple_variable_cdf_dparr6 

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_string
 !! FUNCTION 
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_string(ncid,namevar,values)

  integer,intent(in)   :: ncid
  character(len=*),intent(out) :: values
  character(len=*),intent(in)  :: namevar
  
  integer :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_string

 !!=============================================================
 !!subroutine: m_basic_cdf/get_simple_variable_cdf_string2
 !! FUNCTION 
 !!  Get variable from an open Netcdf file
 !!
 !! OUTPUT
 !!  Only write
 subroutine get_simple_variable_cdf_string2(ncid,namevar,values)

  integer,intent(in)   :: ncid
  character(len=*),intent(out) :: values(:)
  character(len=*),intent(in)  :: namevar
  
  integer :: stId,varid

  stId = nf90_inq_varid(ncid, namevar, varid)
  call test_cdf(stId)
  stId = nf90_get_var(ncid,varid,values)
  call test_cdf(stId)

 end subroutine get_simple_variable_cdf_string2


#endif
end module defs_basic_cdf


 
