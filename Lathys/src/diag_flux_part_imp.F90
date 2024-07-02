!!=============================================================
!!=============================================================
module diag_flux_part_imp
  
  use mpi
  use defs_mpitype
  use defs_variable
  use defs_basis
  use defs_parametre!,only :nhm,dt,gstep,fildat
  use env_ganymede
  use defs_atmospheretype
  use environment
#ifdef HAVE_NETCDF
  use defs_basic_cdf
  use diag_wrt_common_cdf
  use netcdf
#endif
  
#include "q-p_common.h"
  
  implicit none
  
  private :: &
       create_file_name,&
       compute_global_flux
  public :: wrt_flux_part_imp 
  
contains
  !!###############################################################
  
  !!================================================================
  !!subroutine: create_file_name
  !! FUNCTION
  !!  Create the name of the file containing the flux of impacting
  !!  particles
  !! INPUT
  !!  diag_time = time of the diagnostic (in term of number of 
  !!  time steps)
  !! OUTPUT
  !!  name_file=name of the file where particles will be recorded
  
  subroutine create_file_name(name_file,filwrt)
   character(len=*),intent(in) :: filwrt
   character(len=40),intent(out) :: name_file
   integer                       :: diag_time

   write(name_file,'(a14,a)')"Flux_part_imp_",trim(filwrt)
#ifdef HAVE_NETCDF
  name_file = trim(name_file)//".nc"
#endif

  end subroutine create_file_name

  !!###############################################################
  
  !!================================================================
  !!subroutine: compute_global_flux
  !! FUNCTION
  !! Compute the global flux with the computed flux of each process
  !! INPUT
  !!  
  !! OUTPUT
  !!  Global arrays containing flux of impacting particles
  
  subroutine compute_global_flux(flux_surf_g,fdde_surf_g,fdda_surf_g,n_part_imp_g,n_part_del_g)
   use defs_variable
   use defs_mpitype
   use mpi

   integer,intent(inout)                                         :: n_part_imp_g,n_part_del_g
   real(dp),intent(inout),dimension(nphi,ntheta,1)        :: flux_surf_g
   real(dp),intent(inout),dimension(nphi,ntheta,n_inte,1) :: fdde_surf_g
   real(dp),intent(inout),dimension(nphi,ntheta,n_inta,1) :: fdda_surf_g
   integer                                                       :: code,i,j,k,l

        flux_surf_g(:,:,1)=0.d0
        do j=1,n_inte
            fdde_surf_g(:,:,j,1)=0.d0 
        enddo 
        do j=1,n_inta      
            fdda_surf_g(:,:,j,1)=0.d0
        enddo
      do j=1,nphi
         do k=1,ntheta
          call MPI_REDUCE(flux_surf(j,k,1),flux_surf_g(j,k,1),1,&
               MPI_REAL,MPI_SUM,0,mpiinfo%comm,code)
            do l=1,n_inte
              call MPI_REDUCE(fdde_surf(j,k,l,1),fdde_surf_g(j,k,l,1),1,&
                   MPI_REAL,MPI_SUM,0,mpiinfo%comm,code)
            enddo
            do l=1,n_inta
              call MPI_REDUCE(fdda_surf(j,k,l,1),fdda_surf_g(j,k,l,1),1,&
                   MPI_REAL,MPI_SUM,0,mpiinfo%comm,code)
            enddo
         enddo
      enddo
   call MPI_REDUCE(n_part_imp,n_part_imp_g,1,MPI_INTEGER,MPI_SUM,0,mpiinfo%comm,code)
   call MPI_REDUCE(n_part_del,n_part_del_g,1,MPI_INTEGER,MPI_SUM,0,mpiinfo%comm,code)


  end subroutine compute_global_flux
  
  
  !!###############################################################

  !!================================================================
  !!subroutine: wrt_flux_part_imp
  !! FUNCTION
  !! Compute the global flux with the computed flux of each process
  subroutine wrt_flux_part_imp(filwrt)
    character(len=*),intent(in) :: filwrt
    integer :: ncid, stId,ii
    integer,parameter :: nspecies_flux=1
    integer :: dimid(12), varid(75)
    integer,allocatable :: dimspec(:)
    character(len=40) :: name_file
    real(dp),dimension(nphi,ntheta,nspecies_flux)        :: flux_surf_g
    real(dp),dimension(nphi,ntheta,n_inte,nspecies_flux) :: fdde_surf_g
    real(dp),dimension(nphi,ntheta,n_inta,nspecies_flux) :: fdda_surf_g
    integer               :: n_part_imp_g,n_part_del_g
    integer, dimension(3) :: dimidang
    integer               :: nb_species
    character(len=120) :: filename_Efields,filename_Bfields,filename_IonProd, &
          filename_PhoProd, filename_ImpProd, filename_ChXProd,planet,filename_SputSur
    character(len=120) :: Frame0,Date0,Author0,Reference0,Comment0,Units0,frm
    real(dp) :: rmax0,Rexo,Rm,esurf_min,esurf_max,deltae, &
          esurf_grid_mean, esurf_grid_sigm, asurf_grid_mean, asurf_grid_sigm, &
          dtmin,tmax0,Rot_Exosph,Coefficient,r_gridion_min, Surf0, &
          fdde_max, fdda_max, amu, ergev, sqratio
    integer  :: i,j,ie,ia,nionprod_end, &
 npart_cell, part0, ipart, ipart0, npart, Or0, CX0, eof, iunit, ip0, it0, ie0, ia0
    real(dp),dimension(:),allocatable :: masse
    character(len=120) :: species_name

    ! Compute the global flux of impacting particles from all processes
    call compute_global_flux(flux_surf_g,fdde_surf_g,fdda_surf_g,n_part_imp_g,n_part_del_g)
if(sum(fdde_surf_g)/=sum(fdda_surf_g)) then
print*,'ERREUR DIAG',sum(fdde_surf_g),sum(fdde_surf),sum(fdda_surf_g),sum(fdda_surf)
endif

if(mpiinfo%me==0) then
 !   if(trim(Spe%planetname)=='mars') 
     nb_species=1
    allocate(masse(nb_species))
    !if(trim(Spe%planetname)=='mars') then
      ! Mass of the particles in g
      masse(1) = 16.d0 * Spe%S(1)%rmass * amu_pmass * 1d3 ! O+
    !endif


    Frame0     = 'X Sun/Object, Y dusk,/dawn Z South/North'
    Date0      = '06/04/2014'
    Author0    = 'F. Leblanc/L.Leclercq/R.Modolo'
    Reference0 = ' '
    Units0     = 'CGS'
    Comment0   = ' '
    filename_SputSur = ' '
    filename_Efields = 'Elw_'//filwrt//'.nc'
    filename_Bfields = 'Magw_'//filwrt//'.nc'
    filename_PhoProd = 'Phot_'//filwrt//'.nc'
    filename_ChXProd = 'CX_O_'//filwrt//'.nc'
    filename_impProd = 'Eimp_'//filwrt//'.nc'
    !nionprod_start   = 15231
    nionprod_end     = iter
    npart_cell       = 0
    dtmin            = dt
    tmax0            = iter*dt*Spe%ref%inv_gyro
    r_gridion_min    = (Spe%P%radius+alt_part_imp)*1d+5 !in cm
    rmax0            = r_gridion_min
    Rot_Exosph       = 0.d0
    Rm               = 0.d0

    if(iter/=0) then
    ! Flux in cm-2 : we have to divide by the time in sec
    flux_surf_g(:,:,1)=flux_surf_g(:,:,1)/((nionprod_end-nionprod_start)*dt*Spe%ref%inv_gyro*two_pi)
    fdde_surf_g(:,:,:,1)=fdde_surf_g(:,:,:,1)/((nionprod_end-nionprod_start)*dt*Spe%ref%inv_gyro*two_pi)
    fdda_surf_g(:,:,:,1)=fdda_surf_g(:,:,:,1)/((nionprod_end-nionprod_start)*dt*Spe%ref%inv_gyro*two_pi)
    endif

    ! Creation of the file
    call create_file_name(name_file,filwrt)
    write(6,*)' Writing of file = ',name_file

    !--Create the file
    stId = nf90_create(name_file , nf90_clobber, ncid)
    call test_cdf(stId)

    stId = nf90_def_dim(ncid, "nphi",    nphi , dimid(1))
    call test_cdf(stId)
    stId = nf90_def_dim(ncid, "ntet",    ntheta  , dimid(2))
    call test_cdf(stId) 
    stId = nf90_def_dim(ncid, "nener",   n_inte  , dimid(3))
    call test_cdf(stId)
    stId = nf90_def_dim(ncid, "nangl",   n_inta  , dimid(4))
    call test_cdf(stId)
    stId = nf90_def_dim(ncid, "nspec",   nb_species  , dimid(5))
    call test_cdf(stId)
    stId = nf90_def_dim(ncid, "dim_int"   , 1    , dimid(6))
    call test_cdf(stId)
    stId = nf90_def_dim(ncid, "dim_real"  , 1    , dimid(7)) 
    call test_cdf(stId) 
    stId = nf90_def_dim(ncid, "dim_string"  , 120    , dimid(8))
    call test_cdf(stId)
    ii=1

stId = nf90_def_var(ncid, "species_name", nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "alt_of_flux", nf90_int,dimid(7),varid(ii))
call test_cdf(stId) ; ii = ii + 1 
stId = nf90_def_var(ncid, "nphi", nf90_int,dimid(6),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "ntet", nf90_int,dimid(6),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "nener", nf90_int,dimid(6),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "nangl", nf90_int,dimid(6),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "phi", nf90_real,dimid(1),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "theta", nf90_real,dimid(2),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "energy", nf90_real,dimid(3),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "angle", nf90_real,dimid(4),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "n_part_imp", nf90_int,dimid(6),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "n_part_del", nf90_int,dimid(6),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Flux_totsurf", nf90_real,(/dimid(1:2),dimid(5)/),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Fdd_Energy", nf90_real,(/dimid(1:3),dimid(5)/),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Fdd_Angle", nf90_real,(/dimid(1:2),dimid(4),dimid(5)/),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Object",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Input_Efiles",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Input_Bfiles",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Input_PhoProd_files",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Input_ChXProd_files",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Input_ImpProd_files",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Mass_of_species",nf90_real,dimid(5),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Nb_testpart_cell",nf90_real,dimid(6),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Time_step",nf90_real,dimid(7),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Time_max",nf90_real,dimid(7),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Rm",nf90_real,dimid(7),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Rmin",nf90_real,dimid(7),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Rmax",nf90_real,dimid(7),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Rot_Exo", nf90_real,dimid(7),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Frame",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Date_of_writing",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Authors",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Reference",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Units",nf90_char,dimid(8),varid(ii))
call test_cdf(stId) ; ii = ii + 1
stId = nf90_def_var(ncid, "Comments",nf90_char,dimid(8),varid(ii))
stId = nf90_enddef(ncid)
ii = 1
species_name="1-Opl"
stId = nf90_put_var(ncid, varid(ii), species_name)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), alt_part_imp)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), nphi)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), ntheta)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), n_inte)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), n_inta)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), phi_grid)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), theta_grid)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), esurf_grid)
call test_cdf(stId) ;ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), asurf_grid)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), n_part_imp_g)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), n_part_del_g)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), flux_surf_g)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), fdde_surf_g)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), fdda_surf_g)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), Spe%planetname)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), filename_Efields)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), filename_Bfields)
call test_cdf(stId) ; ii = ii + 1
if (nionprod_start == 1) then
 stId = nf90_put_var(ncid, varid(ii), filename_PhoProd)
else
 stId = nf90_put_var(ncid, varid(ii), ' ')
endif
call test_cdf(stId) ; ii = ii + 1
if (nionprod_start.le.2.and.nionprod_end.ge.2) then
 stId = nf90_put_var(ncid, varid(ii), filename_ChXProd)
else
 stId = nf90_put_var(ncid, varid(ii), ' ')
endif
call test_cdf(stId) ; ii = ii + 1
if (nionprod_end == 3) then
 stId = nf90_put_var(ncid, varid(ii), filename_ImpProd)
else
 stId = nf90_put_var(ncid, varid(ii), ' ')
endif
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), masse)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), npart_cell)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), dtmin)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), tmax0)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), Rm)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), r_gridion_min)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), rmax0)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), Rot_Exosph)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), Frame0)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), Date0)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), Author0)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), Reference0)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), Units0)
call test_cdf(stId) ; ii = ii + 1
stId = nf90_put_var(ncid, varid(ii), Comment0)
stId = nf90_close(ncid)
write(6,*)' End of Writing of Flux Surface file'
endif
! re-initializing fluxes to zero
! local array (for each process)
fdda_surf(:,:,:,:) = zero
fdde_surf(:,:,:,:) = zero
flux_surf(:,:,:) = zero
! global array
fdda_surf_g(:,:,:,:)=zero
fdde_surf_g(:,:,:,:) = zero
flux_surf(:,:,:) = zero
! setting the starting intergration to the new iteration
nionprod_start = iter


  end subroutine wrt_flux_part_imp

end module diag_flux_part_imp
