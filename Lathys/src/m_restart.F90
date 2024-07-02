!!=============================================================
!!=============================================================
module m_restart
 use mpi
 use defs_basis
 use defs_mpitype,only     : mpiinfo
 use defs_parametre
 use defs_variable
 use defs_grid
 use part_moment,only : Emtsp3
 use initialisation,only   : h3init
 use m_writeout
 use m_timing,only         : time_get
 use time_variation
 use environment,only : load_temporal_change
#ifdef HAVE_NETCDF 
 use defs_basic_cdf
 use diag_wrt_common_cdf
 use netcdf
#endif
#ifdef HAVE_WAVE_TEST
 use wavetest,only       : awave,kwave,nwave
#endif

#include "q-p_common.h"

 implicit none
 private

 integer,private,save :: re_st = 0 !--Variable determining the restart file name

 private::           &
      read_restart

 public::            &
      wrt_restart,   &
      re_start

contains
 !!#####################################################################

 !**************************** WRT_RESTART ***************
 subroutine wrt_restart()
  !======================================
  ! ecriture du fichier de restart pour le code parallele
  ! ecrit par R. Modolo le 24/06/03
  !--tous les processus ecrivent dans des ficghiers separes

  use m_rand_gen,only        : rand_vars_get
  use defs_species


  character(len=80)  :: namfilrest
  character(len=200) :: msg
  !--tableaux pour la gestion des nombres alatoires
  integer :: istate_1(4)
  real(dp) :: rstate_1(97)
#ifdef HAVE_NETCDF
  integer :: ncid,stId,ii
  integer :: dimid(14),varid(115),dimidloc(3)
  integer,allocatable :: dimspec(:)
#else
  integer :: iunit
#endif
  __WRT_DEBUG_IN("wrt_restart")
  __GETTIME(11,1)!--Timer start

  !--Get variables for random number
  call rand_vars_get(rstate_1,istate_1)

  namfilrest = trim(fildat)
  if(re_st == 1) then
   namfilrest = "1_"//namfilrest(:8);  re_st = 0
  else
   namfilrest = "0_"//namfilrest(:8);  re_st = 1
  endif

  !--Create the file name
  write(msg,'(a2,i4.4,a1,a)')"r_",mpiinfo%me,'_',trim(namfilrest)
  namfilrest = trim(msg)

  call wrt_double(qp_out,"Restart file name : "//namfilrest,wrtscreen,wrtdisk)

#ifdef HAVE_NETCDF
  !--Create the file
  stId = nf90_create(trim(namfilrest)//".nc" , nf90_clobber, ncid); call test_cdf(stId)

  !--Define the dimensions that will define the size of the file.
  call  create_wrt_dimensions_cdf(ncid,dimid,ncm)

  !--Local Definition
  stId = nf90_def_dim(ncid, "istat_dim", size(istate_1),dimidloc(1))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "rstate_dim",size(rstate_1),dimidloc(2))
  call test_cdf(stId)
  stId = nf90_def_dim(ncid, "nwrm", nwrm, dimidloc(3))
  call test_cdf(stId)

  !--Define Speces dimensions
  call species_def_dim_cdf(ncid,dimspec)

  !--Set the global attributes
  call set_global_attribute_cdf(ncid,"Restart")

  ii = 1

  !--Define common variables
  call common_def_var_cdf(ncid,varid,dimid,ii)

  !--Define MPI_INFO
  call mpiinfo_def_var_cdf(ncid,varid,dimid,ii)

  !--Define Speces variables
  call species_def_var_cdf(ncid,varid,dimspec,ii)

  !--Define other variables
  stId = nf90_def_var(ncid, "s_min", QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "s_max", QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "s_r",    QP_NF90_DP,dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "te",    QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "nwrm",     nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "nhm",  nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "nsub",    nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "ntest",   nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "eps",   QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "idisp",   nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "ipe",     nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "iresis",  nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "eps0",  QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "vac",   QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "cd0",   QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dmin",  QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "resis", QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "betai", QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "cs",    QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "psi",   QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "phi",   QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "rmu0",  QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "rho0",  QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "vxmean",QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "theta", QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
#ifdef HAVE_WAVE_TEST
  stId = nf90_def_var(ncid, "nwave",   nf90_int, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "kwave",  QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "awave",   QP_NF90_DP, dimid(1),varid(ii))
  call test_cdf(stId); ii = ii+1
#endif
  stId = nf90_def_var(ncid, "nc1_tot",       nf90_int, dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "vx0_vy0_vz0", QP_NF90_DP, dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "bx0_by0_bz0", QP_NF90_DP, dimid(3),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "iwr",        nf90_int, dimidloc(3),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "istate_1",   nf90_int, dimidloc(1),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "rstate_1", QP_NF90_DP, dimidloc(2),varid(ii))
  call test_cdf(stId); ii = ii+1

  !--FIELDS_DEF
  stId = nf90_def_var(ncid, "Bfield_x", QP_NF90_DP, dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "Bfield_y", QP_NF90_DP, dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "Bfield_z", QP_NF90_DP, dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1              
    stId = nf90_def_var(ncid, "Bfield0_x", QP_NF90_DP, dimid(6:8),varid(ii))
    call test_cdf(stId); ii = ii+1                                                
    stId = nf90_def_var(ncid, "Bfield0_y", QP_NF90_DP, dimid(6:8),varid(ii))
    call test_cdf(stId); ii = ii+1                                                
    stId = nf90_def_var(ncid, "Bfield0_z", QP_NF90_DP, dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1 
  stId = nf90_def_var(ncid, "Efield_x", QP_NF90_DP, dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "Efield_y", QP_NF90_DP, dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "Efield_z", QP_NF90_DP, dimid(6:8),varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "vel_x", QP_NF90_DP, dimid(6:8),   varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "vel_y", QP_NF90_DP, dimid(6:8),   varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "vel_z", QP_NF90_DP, dimid(6:8),   varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "vela_x", QP_NF90_DP, dimid(6:8),  varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "vela_y", QP_NF90_DP, dimid(6:8),  varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "vela_z", QP_NF90_DP, dimid(6:8),  varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "velf_x", QP_NF90_DP, dimid(6:8),  varid(ii))
  call test_cdf(stId); ii = ii+1                                                
  stId = nf90_def_var(ncid, "velf_y", QP_NF90_DP, dimid(6:8),  varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "velf_z", QP_NF90_DP, dimid(6:8),  varid(ii))
  call test_cdf(stId); ii = ii+1  
  stId = nf90_def_var(ncid, "dn   ", QP_NF90_DP, dimid(6:8),   varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dna  ", QP_NF90_DP, dimid(6:8),   varid(ii))
  call test_cdf(stId); ii = ii+1
  stId = nf90_def_var(ncid, "dnf  ", QP_NF90_DP, dimid(6:8),   varid(ii))
  call test_cdf(stId); ii = ii+1
  !--PARTICLES_DEF
  call def_var_particle_cdf(ncid,varid,dimid,ii)

  !--Switch to write mode
  stId = nf90_enddef(ncid); call test_cdf(stId)

  ii = 1;  
  !--Write Common variables into the file
  call common_put_var_cdf(ncid,varid,ii)
  
  !--Write MPI variables into the file
  call mpiinfo_put_var_cdf(ncid,varid,ii)

  !--Write Species infos into the file
  call species_put_var_cdf(Spe,ncid,varid,dimspec,ii)

  !--Write the other variables into the file
  stId = nf90_put_var(ncid, varid(ii), s_min)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), s_max)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), s_r)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), te)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), nwrm)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), nhm)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), nsub)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), ntest)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), eps)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), idisp)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), ipe)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), iresis)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), eps0)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), vac)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), cd0)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), dmin)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), resis)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), betai)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), cs)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), psi)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), phi)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), rmu0)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), rho0)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), vxmean)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), theta)
  call test_cdf(stId); ii = ii+1
#ifdef HAVE_WAVE_TEST
  stId = nf90_put_var(ncid, varid(ii), nwave)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), kwave)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), awave)
  call test_cdf(stId); ii = ii+1
#endif
  stId = nf90_put_var(ncid, varid(ii), nc1_tot)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), (/vx0,vy0,vz0/))
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), (/bx0,by0,bz0/))
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), iwr)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), istate_1)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), rstate_1)
  call test_cdf(stId); ii = ii+1


  stId = nf90_put_var(ncid, varid(ii), Bfield%x)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), Bfield%y)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), Bfield%z)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), Bfield_0%x)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), Bfield_0%y)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), Bfield_0%z)
  call test_cdf(stId); ii = ii+1  
  stId = nf90_put_var(ncid, varid(ii), Efield%x)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), Efield%y)
  call test_cdf(stId); ii = ii+1               
  stId = nf90_put_var(ncid, varid(ii), Efield%z)
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), vel%x  )
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), vel%y  )
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), vel%z  )
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), vela%x )
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), vela%y )
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), vela%z )
  call test_cdf(stId); ii = ii+1   
  stId = nf90_put_var(ncid, varid(ii), velf%x )
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), velf%y )
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), velf%z )
  call test_cdf(stId); ii = ii+1   
  stId = nf90_put_var(ncid, varid(ii), dn     )
  call test_cdf(stId); ii = ii+1              
  stId = nf90_put_var(ncid, varid(ii), dna    )
  call test_cdf(stId); ii = ii+1
  stId = nf90_put_var(ncid, varid(ii), dnf    )
  call test_cdf(stId); ii = ii+1  

  !--Put particles in netcdf file
  call put_var_particle_cdf(particule,nptot,ncid,varid,ii)
  
  !--Close the file
  stId = nf90_close(ncid); call test_cdf(stId)

#else

  iunit = 1
  open (UNIT = iunit, &
       FILE = namfilrest, &
       ACCESS = "sequential", &
       ACTION = "write", &
       STATUS = "unknown",& 
       FORM = "unformatted")

  write(iunit) gstep,s_min,s_max,iter,dt,t,nptot,te

  !--Write Species infos
  call  species_put_var_bin(Spe,iunit)

  write(iunit) fildat,nwrm
  write(iunit) nhm,nsub,ntest,eps,idisp,ipe,iresis
  write(iunit) dmin,vac,eps0,cd0,resis,betai,cs
  write(iunit) bx0,by0,bz0,vx0,vy0,vz0
  write(iunit) psi,phi,theta
#ifdef HAVE_WAVE_TEST
  write(iunit) kwave,awave,nwave
#endif
  write(iunit) nc1_tot,s_r,rmu0,rho0
  write(iunit) vxmean,iwr
  write(iunit) rstate_1,istate_1

  !--Fields
  write(iunit) Efield%x,Efield%y,Efield%z
  write(iunit) Bfield%x,Bfield%y,Bfield%z
  write(iunit) vel%x,vel%y,vel%z
  write(iunit) vela%x,vela%y,vela%z
  write(iunit) dn,dna

  !--Particles
  call put_var_particle_bin(particule,nptot,iunit)

  close(unit=iunit)
#endif

  if(mpiinfo%me == 0) then 
   write(msg,'(4a,4i11,2a,i2)')&
        &"Random generation parameters ",ch10,&
        &"istate_1(1),istate_1(2),istate_1(3),istate_1(4)",ch10,&
        & istate_1(1),istate_1(2),istate_1(3),istate_1(4),ch10,&
        &" Writing mode (0 or 1) :",re_st
   call wrt_double(6,msg,wrtscreen,wrtdisk)
  endif

  __GETTIME(11,2)!--Timer stop
  __WRT_DEBUG_OUT("wrt_restart")
 end subroutine wrt_restart
 !************************* END WRT_RESTART***************

 !************************** READ_RESTART ******************
 subroutine read_restart(rstate_1,istate_1)
  !--lecture du fichier de restart parallele

  use defs_mpitype,only     : mpiinfo
  use defs_species
  !use mpi

  integer,intent(out) :: istate_1(4)
  real(dp),intent(out) :: rstate_1(97)

  integer :: nhm_loc
  integer :: itmp(3)
  character(len=80)  :: namfilrest
  character(len=200) :: msg
#ifdef HAVE_NETCDF
  integer :: ncid,stId
  real(dp) :: rtmp(3)
#else
  integer :: iunit
#endif
  
  !--Select which file use to restart
  if(re_st==1)then
   namfilrest = "1_"//fildat;  re_st = 0
  else
   namfilrest = "0_"//fildat;  re_st = 1
  endif

  !--Create the file name
  write(msg,'(a2,i4.4,a1,a)')"r_",mpiinfo%me,'_',trim(namfilrest)

#ifdef HAVE_NETCDF
  namfilrest = trim(msg)//".nc"

  call wrt_double(6,"Reading file : "//namfilrest,wrtscreen,wrtdisk)

  !--Create the file
  stId = nf90_open(namfilrest ,  nf90_nowrite, ncid)
  call test_cdf(stId)

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

  !--Get Total Number of Particles for this proc
  call get_simple_variable_cdf(ncid,"nptot",nptot)

 
  !--Get Species informations
  call species_get_var_cdf(Spe,ncid)

  call get_simple_variable_cdf(ncid,"te",te)
  call get_simple_variable_cdf(ncid,"nwrm",nwrm)
  call get_simple_variable_cdf(ncid,"nhm",nhm_loc)
  call get_simple_variable_cdf(ncid,"nsub",nsub)
  call get_simple_variable_cdf(ncid,"ntest",ntest)
  call get_simple_variable_cdf(ncid,"eps",eps)
  call get_simple_variable_cdf(ncid,"idisp",idisp)
  call get_simple_variable_cdf(ncid,"ipe",ipe)
  call get_simple_variable_cdf(ncid,"iresis",iresis)
  call get_simple_variable_cdf(ncid,"dmin",dmin)
  call get_simple_variable_cdf(ncid,"vac",vac)
  call get_simple_variable_cdf(ncid,"eps0",eps0)
  call get_simple_variable_cdf(ncid,"cd0",cd0)
  call get_simple_variable_cdf(ncid,"resis",resis)
  call get_simple_variable_cdf(ncid,"betai",betai)
  call get_simple_variable_cdf(ncid,"cs",cs)

  call get_simple_variable_cdf(ncid,"bx0_by0_bz0",rtmp)
  bx0 = rtmp(1);  by0 = rtmp(2);  bz0 = rtmp(3);
  call get_simple_variable_cdf(ncid,"vx0_vy0_vz0",rtmp)
  vx0 = rtmp(1);  vy0 = rtmp(2);  vz0 = rtmp(3);
  call get_simple_variable_cdf(ncid,"psi",psi)
  call get_simple_variable_cdf(ncid,"phi",phi)
  call get_simple_variable_cdf(ncid,"theta",theta)
#ifdef HAVE_WAVE_TEST
  call get_simple_variable_cdf(ncid,"nwave",nwave)
  call get_simple_variable_cdf(ncid,"kwave",kwave)
  call get_simple_variable_cdf(ncid,"awave",awave)
#endif
  call get_simple_variable_cdf(ncid,"nc1_tot",itmp)
  call set_grid(itmp,1)
  call get_simple_variable_cdf(ncid,"s_r",s_r)
  call get_simple_variable_cdf(ncid,"rmu0",rmu0)
  call get_simple_variable_cdf(ncid,"rho0",rho0)
  call get_simple_variable_cdf(ncid,"vxmean",vxmean)
  call get_simple_variable_cdf(ncid,"iwr",iwr)
  call get_simple_variable_cdf(ncid,"rstate_1",rstate_1)
  call get_simple_variable_cdf(ncid,"istate_1",istate_1)

  call get_simple_variable_cdf(ncid,"Efield_x",Efield%x)
  call get_simple_variable_cdf(ncid,"Efield_y",Efield%y)
  call get_simple_variable_cdf(ncid,"Efield_z",Efield%z)
  call get_simple_variable_cdf(ncid,"Bfield_x",Bfield%x)  
  call get_simple_variable_cdf(ncid,"Bfield_y",Bfield%y)  
  call get_simple_variable_cdf(ncid,"Bfield_z",Bfield%z)  
  call get_simple_variable_cdf(ncid,"Bfield0_x",Bfield_0%x)  
  call get_simple_variable_cdf(ncid,"Bfield0_y",Bfield_0%y)  
  call get_simple_variable_cdf(ncid,"Bfield0_z",Bfield_0%z)  
  call get_simple_variable_cdf(ncid,"vel_x",vel%x)  
  call get_simple_variable_cdf(ncid,"vel_y",vel%y)  
  call get_simple_variable_cdf(ncid,"vel_z",vel%z)  

  call get_simple_variable_cdf(ncid,"vela_x",vela%x)  
  call get_simple_variable_cdf(ncid,"vela_y",vela%y)  
  call get_simple_variable_cdf(ncid,"vela_z",vela%z)  
  call get_simple_variable_cdf(ncid,"velf_x",velf%x)  
  call get_simple_variable_cdf(ncid,"velf_y",velf%y)  
  call get_simple_variable_cdf(ncid,"velf_z",velf%z)   
  call get_simple_variable_cdf(ncid,"dn",dn)  
  call get_simple_variable_cdf(ncid,"dna",dna) 
  call get_simple_variable_cdf(ncid,"dnf",dnf) 

  !--Get particles from netcdf file
  call get_var_particle_cdf(particule,nptot,ncid)

  !--Close the file
  stId = nf90_close(ncid); call test_cdf(stId)

#else  

  namfilrest = trim(msg)
  iunit = 1

  open (UNIT = iunit, &
       FILE = namfilrest, &
       ACCESS = "sequential", &
       ACTION = "read", &
       STATUS = "unknown",& 
       FORM = "unformatted")

  read(iunit) gstep,s_min,s_max,iter,dt,t,nptot,te

  !--Get Species infos
  call species_get_var_bin(Spe,iunit)

  read(iunit) fildat,nwrm
  read(iunit) nhm_loc,nsub,ntest,eps,idisp,ipe,iresis
  read(iunit) dmin,vac,eps0,cd0,resis,betai,cs
  read(iunit) bx0,by0,bz0,vx0,vy0,vz0
  read(iunit) psi,phi,theta
#ifdef HAVE_WAVE_TEST
  read(iunit) kwave,awave,nwave
#endif
  read(iunit) itmp,s_r,rmu0,rho0
  call set_grid(itmp,1)
  read(iunit) vxmean,iwr
  read(iunit) rstate_1,istate_1

  !--Fields
  read(iunit) Efield%x,Efield%y,Efield%z
  read(iunit) Bfield%x,Bfield%y,Bfield%z
  read(iunit) vel%x,vel%y,vel%z
  read(iunit) vela%x,vela%y,vela%z
  read(iunit) dn,dna
  write(msg,'(a,i4.4)')"End of reading process fiels : ",mpiinfo%me
  call wrt_double(6,msg,wrtscreen,wrtdisk)

  !--Particles
  call get_var_particle_bin(particule,nptot,iunit)

  close(unit=iunit)
#endif

  !--Set the max number of iteration as the maximum between the stored
  !value and the actual one in parameter module
  nhm = max(nhm,nhm_loc)

  write(msg,'(a,i4.4)')" End of reading process fiels : ",mpiinfo%me
  call wrt_double(6,msg,wrtscreen,wrtdisk)

 end subroutine read_restart
 !************************* END READ_RESTART ******************



 !*********************************** RE_START *******************************
 subroutine re_start(restart_opt)
  !--subroutine de restart

  use defs_mpitype,only    : mpiinfo
  use field_lissage,only         : smth_func
  use part_moment,only          : momt3d
  use m_rand_gen,only        : rand_vars_put

  integer,intent(in) :: restart_opt  

  integer :: istate_1(4)
  real(dp) :: rstate_1(97)
  character(len=200) :: msg

  __WRT_DEBUG_IN("re_start")
  __GETTIME(4,1) !--Timer start


  re_st = restart_opt

  !--Lecture des donnes archives
  call read_restart(rstate_1,istate_1)

  write(msg,'(a,i8)')"iteration ",iter
  call wrt_double(6,msg,wrtscreen,wrtdisk)

  !--Mise a jour des parametres des nombres aleatoires
  call rand_vars_put(rstate_1,istate_1,1)

  write(msg,'(4a,4i11,2a,i2)')&
       &" Random generation parameter ",ch10,&
       &" istate_1(1),istate_1(2),istate_1(3),istate_1(4)",ch10,&
       &  istate_1(1),istate_1(2),istate_1(3),istate_1(4),ch10,&
       &" Writing mode (0 or 1) :",re_st
  call wrt_double(6,msg,wrtscreen,wrtdisk)

  !--Affichage des parametres du restart
  call h3init()

  !--Collect dn,u,dnf,uf
  !call momt3d(2)
  call Emtsp3(particule,dn,vel,1,nptot, &
          &      gstep,s_min_loc,&
          &      mpiinfo,&
        &      nc1)

  !--Smoothing (optional)
 ! if(iwr(6) == 1)then
 !  call smth_func(dnf,nc,mpiinfo)
 ! end if

  !--Set Bh = B
  Bfield_h = Bfield
  
  ! intialise temporal change of incident plasma
  !call load_temporal_change(Spe)

  __GETTIME(4,2)!--Timer stop
  __WRT_DEBUG_OUT("re_start")
 end subroutine re_start
 !************************************ END RE_START **********************************

end module m_restart
