!!=============================================================
!!=============================================================
!!module: defs_particletype
!! NAME
!!  defs_particletype (MMancini)
!!
!! FUNCTION
!!  Definition of type for particles
!! NOTE
module defs_particletype

 use defs_basis
#ifndef NOTHAVE_MPI 
use mpi
#endif

 implicit none
 private

 integer,public,parameter :: iXb = &
#ifdef HAVE_DOUBLE_PRECISION
      &      i4b
#else
      &      i2b
#endif
 !--This parameter iXb is needed because the heterogeneous type
 !--particles have to occupe a number of bytes multiple of 8.
 !--(gcc, ifort do that). If it is not, MPI create a type which size
 !--is smaller than particletype.

 integer,public,save :: MPI_particle !--MPI planes


 !!=============================================================
 !!--Type for three-dimensional vector
 type,public :: particletype

  real(dp) :: pos(3) !--Position of macro-particle
  real(dp) :: vel(3) !--Velocity 
  real(dp) :: mass   !--Mass
  real(dp) :: char   !--Charge
  real(dp) :: exc    !--Charge exchange
  integer(iXb):: orig!--Origin of the macro-particle:
                     !  orig = 0 => Solar wind
                     !  orig = 1 => Photionisation
                     !  orig = 2 => Ionisation by impact
                     !  orig = 3 => Ionospheric loading
  integer(iXb):: dir !--dir = -1 if going out of system
                     !  if it is going in the neighbours:
                     !  N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
                     !  dir = 0 otherwise
  real(dp):: tbirth ! birth time of the particle 
  
  ! integer(i1b)::bo,ba!--this two vatiables are not used for the moment
  !                    !but they are used (in memory) in any case so I
  !could put some information here. Note that I could use a only
  !variable of size integer(i2b)

 end type particletype

 public ::                  &
      particle_type_size,   &!--Gives the size of the particletype
      init_MPI_particle,    &!--Create the MPI_type for particles
      free_MPI_particle,    &!--Free the MPI_type for particles
      set_zero_particle      !--Set to zero all fields of particletype
 
#ifdef HAVE_NETCDF
 public ::                  &
      def_var_particle_cdf, &!--Define cdf variable for particletype
      put_var_particle_cdf, &!--Put particletype in a cdf file
      get_var_particle_cdf   !--Get particletype from a cdf file
#else
 public ::                  &
      put_var_particle_bin, &!--Put particletype in a binary file
      get_var_particle_bin   !--Get particletype from a binary file
#endif

 interface set_zero_particle
  module procedure set_zero_particle_1d
  module procedure set_zero_particle_2d
 end interface

contains
 !!#####################################################################


 !!=============================================================
 !!routine: defs_particletype/particle_type_size
 !!
 !! FUNCTION
 !!  Gives the memory size (in bit) of the particletype 
 function particle_type_size()
  integer :: particle_type_size
  type(particletype) :: part 
  integer(1) :: sizeof(9999)

  particle_type_size = size(transfer(part,sizeof))

 end function  particle_type_size

!!=============================================================
 !!routine: defs_particletype/init_MPI_particle
 !!
 !! FUNCTION
 !!  Create the MPI types for particles
 !!         
 !! IN
 !! OUT
 subroutine init_MPI_particle()
#include "q-p_common.h"
#ifndef NOTHAVE_MPI 
!use mpi



  integer :: types(8)
  integer :: length_block(8)
  integer(kind=MPI_ADDRESS_KIND) :: address(8),displacement(8)
  integer :: ioerr,ii,MPI_iXb,MPI_particle_size
  type(particletype) :: parts

#ifdef HAVE_DOUBLE_PRECISION
  MPI_iXb = MPI_INTEGER4
#else
  MPI_iXb = MPI_INTEGER2
#endif

  !--Construction of type of data
  types = (/QP_MPI_DP,QP_MPI_DP,QP_MPI_DP,QP_MPI_DP,QP_MPI_DP,&
       &    MPI_iXb,MPI_iXb,QP_MPI_DP/)

  length_block = (/3,3,1,1,1,1,1,1/)

  call MPI_GET_ADDRESS(parts%pos ,address(1),ioerr)
  call MPI_GET_ADDRESS(parts%vel ,address(2),ioerr)
  call MPI_GET_ADDRESS(parts%mass,address(3),ioerr)
  call MPI_GET_ADDRESS(parts%char,address(4),ioerr)
  call MPI_GET_ADDRESS(parts%exc ,address(5),ioerr)
  call MPI_GET_ADDRESS(parts%orig,address(6),ioerr)
  call MPI_GET_ADDRESS(parts%dir ,address(7),ioerr)
  call MPI_GET_ADDRESS(parts%tbirth ,address(8),ioerr)
  ! call MPI_GET_ADDRESS(parts%bo  ,address(8),ioerr)
  ! call MPI_GET_ADDRESS(parts%ba  ,address(9),ioerr)

  !--Compute the displacement from the starting address
  do ii = 1,8
   displacement(ii) = address(ii) - address(1)
  enddo
    
  !--Creation of the MPI-type
  call MPI_TYPE_CREATE_STRUCT(8,length_block,displacement,types,MPI_particle,ioerr)

  !--Compute the size of the MPI_particle type
  call MPI_TYPE_SIZE(MPI_particle,MPI_particle_size,ioerr)

  !--Control on the size of types related to particles
  if(MPI_particle_size /= particle_type_size())then
   print *,'ERROR: sizes of particletype MPI_particle differs ',MPI_particle_size,particle_type_size()
   stop
  endif

  !--Commit of the type
  call MPI_TYPE_COMMIT(MPI_particle,ioerr)
#endif
 end subroutine init_MPI_particle

 !!=============================================================
 !!routine: defs_particletype/free_MPI_particle
 !!
 !! FUNCTION
 !!  Free the MPI types for particles 
 !!         
 !! IN
 !! OUT
 subroutine free_MPI_particle()
#ifndef NOTHAVE_MPI
!use mpi
  integer :: ioerr
  call MPI_TYPE_FREE(MPI_particle,ioerr)
#endif 
 end subroutine free_MPI_particle

 !!=============================================================
 !!routine: defs_particletype/set_zero_particle_1d
 !!
 !! FUNCTION
 !!  Set to zero to the field of a particles 1d-array
 !!         
 !! IN
 !! OUT
 subroutine set_zero_particle_1d(particle,n1,n2)
  integer,intent(in) :: n1 !--Start index
  integer,intent(in) :: n2 !--Stop index
  type(particletype),intent(inout) :: particle(:)

  integer :: ii

  do ii = n1,n2
   particle(ii)%pos  = zero
   particle(ii)%vel  = zero
   particle(ii)%mass = zero
   particle(ii)%char = zero
   particle(ii)%exc  = zero
   particle(ii)%orig = 0_iXb
   particle(ii)%dir  = 0_iXb
   particle(ii)%tbirth = 0_iXb
  enddo
 end subroutine set_zero_particle_1d

 !!=============================================================
 !!routine: defs_particletype/set_zero_particle_2d
 !!
 !! FUNCTION
 !!  Set to zero to the field of a particles 1d-array
 !!         
 !! IN
 !! OUT
 subroutine set_zero_particle_2d(particle)
  type(particletype),intent(inout) :: particle(:,:)

  integer :: ii,jj
  do jj = 1,size(particle,dim=2)
   do ii = 1,size(particle,dim=1)
    particle(ii,jj)%pos  = zero
    particle(ii,jj)%vel  = zero
    particle(ii,jj)%mass = zero
    particle(ii,jj)%char = zero
    particle(ii,jj)%exc  = zero
    particle(ii,jj)%orig = 0_iXb
    particle(ii,jj)%dir  = 0_iXb
    particle(ii,jj)%tbirth = 0_iXb
   enddo
  enddo
 end subroutine set_zero_particle_2d

#ifdef HAVE_NETCDF
 !!=============================================================
 !!routine: defs_particletype/def_var_particle_cdf
 !!
 !! FUNCTION
 !!  Define variable for write the particletype with netcdf
 !!         
 !! IN
 !! OUT
 subroutine def_var_particle_cdf(ncid,varid,dimid,num)
  use netcdf
  use defs_basic_cdf,only     : test_cdf
  integer,intent(in) :: ncid
  integer,intent(inout) :: num
  integer,intent(in) :: dimid(:)
  integer,intent(inout) :: varid(:)
  integer :: stId,nf90_iXb

#ifdef HAVE_DOUBLE_PRECISION
  nf90_iXb = nf90_short
#else
  nf90_iXb = nf90_int
#endif

  stId = nf90_def_var(ncid, "particule_x", QP_NF90_DP, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_y", QP_NF90_DP, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_z", QP_NF90_DP, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_vx", QP_NF90_DP, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_vy", QP_NF90_DP, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_vz", QP_NF90_DP, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_mass", QP_NF90_DP, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_char", QP_NF90_DP, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_exc ", QP_NF90_DP, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_orig", nf90_iXb, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_dir ", nf90_iXb, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1
  stId = nf90_def_var(ncid, "particule_tbirth ", QP_NF90_DP, dimid(10),varid(num)) 
  call test_cdf(stId); num = num+1  
 end subroutine def_var_particle_cdf

 !!=============================================================
 !!routine: defs_particletype/put_var_particle_cdf
 !!
 !! FUNCTION
 !!  Put variable writing particletype in netcdf
 !!         
 !! IN
 !! OUT
 subroutine put_var_particle_cdf(particule,nptot,ncid,varid,num)
  use netcdf
  use defs_basic_cdf,only     : test_cdf
  type(particletype),intent(in) :: particule(:)
  integer,intent(in) :: ncid,nptot
  integer,intent(inout) :: num
  integer,intent(inout) :: varid(:)
  integer :: stId,ii,nwrite,inum
  integer :: lidx,hidx

  !--QP_MAXSIZEWRITE Size of a real(dp) array which can be written at the same time for
  !  a single proc (it is machine dependent we fix it for ciclad
  !  cluster) defined in q-p_common.h
  nwrite = nptot/QP_MAXSIZEWRITE

  do ii = 0,nwrite
   inum = num
   lidx = 1+ii*QP_MAXSIZEWRITE
   hidx = min(nptot,(ii+1)*QP_MAXSIZEWRITE)
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%pos(1),(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%pos(2),(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%pos(3),(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%vel(1),(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%vel(2),(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%vel(3),(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%mass,(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%char,(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%exc,(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%orig,(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%dir,(/lidx/))
   call test_cdf(stId); inum = inum+1
   stId = nf90_put_var(ncid, varid(inum), particule(lidx:hidx)%tbirth,(/lidx/))
   call test_cdf(stId); inum = inum+1   
  enddo
  num = inum
 end subroutine put_var_particle_cdf

 !!=============================================================
 !!routine: defs_particletype/get_var_particle_cdf
 !!
 !! FUNCTION
 !!  Get variable particletype from a netcdf file
 !!         
 !! IN
 !! OUT
 subroutine get_var_particle_cdf(particule,nptot,ncid)
  use netcdf
  use defs_basic_cdf,only     : test_cdf,get_simple_variable_cdf

  type(particletype),intent(inout) :: particule(:)
  integer,intent(in) :: ncid,nptot

  call get_simple_variable_cdf(ncid,"particule_x"   ,particule(:nptot)%pos(1))
  call get_simple_variable_cdf(ncid,"particule_y"   ,particule(:nptot)%pos(2))
  call get_simple_variable_cdf(ncid,"particule_z"   ,particule(:nptot)%pos(3))
  call get_simple_variable_cdf(ncid,"particule_vx"  ,particule(:nptot)%vel(1))
  call get_simple_variable_cdf(ncid,"particule_vy"  ,particule(:nptot)%vel(2))
  call get_simple_variable_cdf(ncid,"particule_vz"  ,particule(:nptot)%vel(3))
  call get_simple_variable_cdf(ncid,"particule_mass",particule(:nptot)%mass)
  call get_simple_variable_cdf(ncid,"particule_char",particule(:nptot)%char)
  call get_simple_variable_cdf(ncid,"particule_exc" ,particule(:nptot)%exc)
  call get_simple_variable_cdf(ncid,"particule_orig",particule(:nptot)%orig)
  call get_simple_variable_cdf(ncid,"particule_dir" ,particule(:nptot)%dir)
  !call get_simple_variable_cdf(ncid,"particule_tbirth" ,particule(:nptot)%tbirth)

 end subroutine get_var_particle_cdf

#else

 !!=============================================================
 !!routine: defs_particletype/put_var_particle_bin
 !!
 !! FUNCTION
 !!  Put variable of particletype in binary file
 !!         
 !! IN
 !! particule(particletype)=particule variable
 !! nptot= assigned particles
 !! unit= output unit
 !! OUT
 subroutine put_var_particle_bin(particule,nptot,unit)

  type(particletype),intent(in) :: particule(:)
  integer,intent(in) :: unit,nptot
  integer :: ii,nwrite
  integer :: lidx,hidx

  !--QP_MAXSIZEWRITE Size of a real(dp) array which can be written at the same time for
  !  a single proc (it is machine dependent we fix it for ciclad
  !  cluster) defined in q-p_common.h
  nwrite = nptot/QP_MAXSIZEWRITE
  
  do ii = 0,nwrite
   lidx = 1+ii*QP_MAXSIZEWRITE
   hidx = min(nptot,(ii+1)*QP_MAXSIZEWRITE)
   
   write(unit) particule(lidx:hidx)%pos(1)
   write(unit) particule(lidx:hidx)%pos(2)
   write(unit) particule(lidx:hidx)%pos(3)
   write(unit) particule(lidx:hidx)%vel(1)
   write(unit) particule(lidx:hidx)%vel(2)
   write(unit) particule(lidx:hidx)%vel(3)
   write(unit) particule(lidx:hidx)%mass
   write(unit) particule(lidx:hidx)%char
   write(unit) particule(lidx:hidx)%exc
   write(unit) particule(lidx:hidx)%orig 
   write(unit) particule(lidx:hidx)%dir   
  enddo

 end subroutine put_var_particle_bin

 !!=============================================================
 !!routine: defs_particletype/get_var_particle_bin
 !!
 !! FUNCTION
 !!  Get variable particletype from a binary file
 !!         
 !! IN
 !! nptot= assigned particles
 !! unit= output unit
 !! OUT
 !!  particule(particletype)=particule variable
 !!
 subroutine get_var_particle_bin(particule,nptot,unit)

  type(particletype),intent(inout) :: particule(:)
  integer,intent(in) :: unit,nptot

  integer :: ii,nwrite
  integer :: lidx,hidx

  !--QP_MAXSIZEWRITE Size of a real(dp) array which can be written at the same time for
  !  a single proc (it is machine dependent we fix it for ciclad
  !  cluster) defined in q-p_common.h
  nwrite = nptot/QP_MAXSIZEWRITE
  
  do ii = 0,nwrite
   lidx = 1+ii*QP_MAXSIZEWRITE
   hidx = min(nptot,(ii+1)*QP_MAXSIZEWRITE)

   read(unit) particule(lidx:hidx)%pos(1)
   read(unit) particule(lidx:hidx)%pos(2)
   read(unit) particule(lidx:hidx)%pos(3)
   read(unit) particule(lidx:hidx)%vel(1)
   read(unit) particule(lidx:hidx)%vel(2)
   read(unit) particule(lidx:hidx)%vel(3)
   read(unit) particule(lidx:hidx)%mass
   read(unit) particule(lidx:hidx)%char
   read(unit) particule(lidx:hidx)%exc
   read(unit) particule(lidx:hidx)%orig 
   read(unit) particule(lidx:hidx)%dir
  enddo

 end subroutine get_var_particle_bin
 
#endif
end module defs_particletype







