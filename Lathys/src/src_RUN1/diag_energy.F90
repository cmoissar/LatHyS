!!=============================================================
!!=============================================================
!!module: m_energy
!! NAME
!!  m_energy (RModolo,MMancini)
!!
!! FUNCTION
!!  Contains routine for the diagnostic and controls
!!
module diag_energy

 use defs_basis 
 use defs_mpitype
 use defs_arr3dtype
 use m_timing,only         : time_get
 use m_writeout
 use mpi
#include "q-p_common.h"

 implicit none
 private

 private::                     &
      controls                  !--Control on tranches dividing x-axis

 public::                      &
      energy_proc               !--Compute diagnostic

contains
 !!#####################################################################


 !!=============================================================
 !!=============================================================
 !!module: m_energy/energy_proc
 !! NAME
 !!  energy_proc (RModolo,MMancini)
 !!
 !! FUNCTION
 !!  Compute energies, mean quantities on grids
 !!
 !! IN
 !!  Preg(count_particle_type)=infos manage of particles
 !!
 !! SIDE EFFECT
 !!  diag_info(diag_type)=infos on all procs mean quantities
 !!  
 subroutine energy_proc(Preg,diag_info)

  use defs_diag_type,only     : diag_type
  use defs_grid
  use defs_parametre,only :ns
  use defs_variable,only           : Bfield,Efield,particule,&
       &                        nptot,iter,s_min,s_max,Spe
  use defs_counts_types,only  : count_particle_type
  use particle_creation,only           : ntot_entree

  type(count_particle_type),intent(in) :: Preg
  type(diag_type),intent(inout) :: diag_info

  integer :: nptot_all,ioerr,is
  real(dp) :: ekin,eterm, psmoy,etot
  real(dp) :: ekin_proc,eterm_proc, psmoy_proc,emag_moy,ekin_moy
  real(dp) :: dnp_inv, nptot_inv ,nmacro_inv
  integer,allocatable :: sendinteger(:),sendinteger_all(:)
  real(dp),dimension(3) :: b_moy,b_quad,b_moy_proc,&
       &                   qms_proc,v_moy_proc,v_qua_proc
  real(dp),dimension(3) :: v_moy,v_qua,qmps,emag_proc,emag
  character(len=200) :: msg

  __WRT_DEBUG_IN("energy_proc")
  __GETTIME(43,1)!--Timer start 

  !--Put to zero
  emag_proc = zero
  emag_moy = zero
  ekin_moy = zero
  b_moy_proc = zero
  
  !--Total nuber of macro-particles in the box
  call MPI_ALLREDUCE(nptot,nptot_all,1,MPI_INTEGER,MPI_SUM,mpiinfo%comm,ioerr)  

  nptot_inv = one/real(nptot_all,dp)
  nmacro_inv = real(sum(Spe%S(:)%ng),dp)/(real(nptot_all,dp))
  dnp_inv = one/real(nc_tot(1)*nc_tot(2)*nc_tot(3),dp)
  
  !--Magnetic Energy on any proc
  emag_proc(1) = sum(Bfield%x(1:nc(1),1:nc(2),1:nc(3))**two)
  emag_proc(2) = sum(Bfield%y(1:nc(1),1:nc(2),1:nc(3))**two)
  emag_proc(3) = sum(Bfield%z(1:nc(1),1:nc(2),1:nc(3))**two)

  !--Total Magnetic Energy on proc 0
  call MPI_REDUCE(emag_proc,emag,3,QP_MPI_DP,MPI_SUM,0,mpiinfo%comm,ioerr)
  if(mpiinfo%me == 0) emag = emag*dnp_inv

  !--Mean Magnetic Field on any proc
  b_moy_proc(1) = sum(Bfield%x(1:nc(1),1:nc(2),1:nc(3)))
  b_moy_proc(2) = sum(Bfield%y(1:nc(1),1:nc(2),1:nc(3)))
  b_moy_proc(3) = sum(Bfield%z(1:nc(1),1:nc(2),1:nc(3)))

  !--Total Mean Magnetic Field on all procs 
  call MPI_REDUCE(b_moy_proc,b_moy,3,QP_MPI_DP,MPI_SUM,0,mpiinfo%comm,ioerr)
  if(mpiinfo%me == 0) b_moy = b_moy*dnp_inv

  !--Quadratic deviation of the Magnetic Field 
  b_quad(1) = emag(1)-b_moy(1)**two
  b_quad(2) = emag(2)-b_moy(2)**two
  b_quad(3) = emag(3)-b_moy(3)**two

  !--Kinetic Energy for any proc
  ! E = m*(vx**2+vy**2+vz**2)
  ! This sum is splitted because when nptot is big, for simple
  ! precision the results is not correct with a single sum
  ekin_proc = sum(particule(1:nptot/2)%mass*(      &
       &          particule(1:nptot/2)%vel(1)**two + &
       &          particule(1:nptot/2)%vel(2)**two + &
       &          particule(1:nptot/2)%vel(3)**two))
  ekin_proc = ekin_proc+ sum(particule(nptot/2+1:nptot)%mass*(      &
       &          particule(nptot/2+1:nptot)%vel(1)**two + &
       &          particule(nptot/2+1:nptot)%vel(2)**two + &
       &          particule(nptot/2+1:nptot)%vel(3)**two))


  !--Total Kinetic energy on proc 0
  call MPI_REDUCE(ekin_proc,ekin,1,QP_MPI_DP,MPI_SUM,0,mpiinfo%comm,ioerr)
  if(mpiinfo%me == 0) ekin = ekin*nmacro_inv
 
  !--Mean Mass of the particles
  psmoy_proc = sum(particule(1:nptot)%mass)
  !--Total Mass on all procs
  call MPI_ALLREDUCE(psmoy_proc,psmoy,1,QP_MPI_DP,MPI_SUM,mpiinfo%comm,ioerr)
  psmoy = psmoy*nmacro_inv

  !--Moment of all particles for any proc
  qms_proc(1) = sum(particule(1:nptot)%vel(1)*particule(1:nptot)%mass)
  qms_proc(2) = sum(particule(1:nptot)%vel(2)*particule(1:nptot)%mass)
  qms_proc(3) = sum(particule(1:nptot)%vel(3)*particule(1:nptot)%mass)
  !--Total moment on all proc
  call MPI_ALLREDUCE(qms_proc,qmps,3,QP_MPI_DP,MPI_SUM,mpiinfo%comm,ioerr)
  qmps = qmps*nmacro_inv

  !--Thermal energy (vibrations around the center of mass)
  eterm_proc = sum(particule(1:nptot)%mass*(&
       &          (particule(1:nptot)%vel(1)-qmps(1)/psmoy)**two +&
       &          (particule(1:nptot)%vel(2)-qmps(2)/psmoy)**two +&
       &          (particule(1:nptot)%vel(3)-qmps(3)/psmoy)**two))
  !--Total termal energy on proc 0
  call MPI_REDUCE(eterm_proc,eterm,1,QP_MPI_DP,MPI_SUM,0,mpiinfo%comm,ioerr)
  eterm = eterm*nmacro_inv
 
  !--Velocity of particles for any proc
  v_moy_proc(1) = sum(particule(1:nptot)%vel(1))
  v_moy_proc(2) = sum(particule(1:nptot)%vel(2))
  v_moy_proc(3) = sum(particule(1:nptot)%vel(3))
  !--Total Velocity of particles on all procs
  call MPI_ALLREDUCE(v_moy_proc,v_moy,3,QP_MPI_DP,MPI_SUM,mpiinfo%comm,ioerr)
  v_moy = v_moy*nptot_inv

  !--On calcul l'ecart quadratique à la vitesse myenne suivant chacune des directions
  v_qua_proc(1) = sum((particule(1:nptot)%vel(1)-v_moy(1))**two)
  v_qua_proc(2) = sum((particule(1:nptot)%vel(2)-v_moy(2))**two)
  v_qua_proc(3) = sum((particule(1:nptot)%vel(3)-v_moy(3))**two)

  !--On envoie les contributions des différents processeurs vers le processeur 0
  v_qua = zero
  call MPI_REDUCE(v_qua_proc,v_qua,3,QP_MPI_DP,MPI_SUM,0,mpiinfo%comm,ioerr)
  if(mpiinfo%me == 0) v_qua = v_qua*nptot_inv

  !--Calcul de l'energie totale
  etot = ekin + sum(emag)

  write(msg,'(2a,4(a,a13,f12.5))')&
       & ch10," _________ Energies __________",&
       & ch10, "   Etotal   = ",etot, &
       & ch10, "   Emagne   = ",sum(emag), &
       & ch10, "   Ekinet   = ",ekin, &
       & ch10, "   Etermi   = ",eterm
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  !--Sendinteger to send many indices at the same time
  allocate(sendinteger(3+Spe%ns))
  allocate(sendinteger_all(3+Spe%ns))

  sendinteger(:3) =(/  ntot_entree, &
       &               Preg%out_xp,&
       &               Preg%out_xm/)

  do is=1,Spe%ns
   sendinteger(3+is) = Spe%S(is)%ip_inj
  enddo

  !--Total nuber of macro-particles in the box
  call MPI_REDUCE(sendinteger,sendinteger_all,3+Spe%ns,MPI_INTEGER, &
       &          MPI_SUM,0,mpiinfo%comm,ioerr)

  !--Ecriture dans les tableaux temporels 
  if (mpiinfo%me == 0) then
   diag_info%ener_mag(iter)      = sum(emag)
   diag_info%ener_kin(iter)      = ekin
   diag_info%ener_tot(iter)      = etot
   diag_info%ener_therm(iter)    = eterm
   diag_info%ener_mag_moy(iter)  = emag_moy
   diag_info%ener_kin_moy(iter)  = ekin_moy
   diag_info%champ_bqua(:,iter)  = b_quad
   diag_info%champ_bmoy(:,iter)  = b_moy 
   diag_info%vitesse_moy(:,iter) = v_moy
   diag_info%vitesse_qua(:,iter) = v_qua
   diag_info%sv(:,iter)          = qmps
   diag_info%ener_kins(iter)     = ekin
   diag_info%n_part(iter)        = nptot_all
   diag_info%n_in(iter)          = sendinteger_all(1)
   diag_info%n_out(iter)         = sendinteger_all(2)+sendinteger_all(3)
   diag_info%n_out_g(iter)       = sendinteger_all(3)
   do is = 1, ns
    diag_info%n_in_spec(iter,is) = sendinteger_all(3+ns)
   enddo
  endif

  !--Print particle balancing
  write(msg,'(2a,5(2a,i13))')&
       &ch10," ______ Particle Balance _____",&
       &ch10,"   Nptot    = ",nptot_all,&
       &ch10,"   ntot_in  = ", sendinteger_all(1),&
       &ch10,"   ntot_out = ", sendinteger_all(2)+sendinteger_all(3),&
       &ch10,"   out-in   = ", sendinteger_all(2)+sendinteger_all(3)-sendinteger_all(1),&
       &ch10,"   ksr      = ", sendinteger_all(3)
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  deallocate(sendinteger,sendinteger_all)

  !--Control of many quantities dividing x-axis in tranche
  !  useful for debugging. Comment it, does not delete it!
  ! call controls(particule,Efield,Bfield,Spe,nc1,nptot,s_min(1),s_max(1),mpiinfo)

  __GETTIME(43,2)!--Timer stop
  __WRT_DEBUG_OUT("energy_proc")
  !*********************************** END ENERGY_PROC *****************************
 end subroutine energy_proc


 !!=============================================================
 !!=============================================================
 !!module: m_energy/controls
 !! NAME
 !!  controls (MMancini)
 !!
 !! FUNCTION
 !!  Compute energies, mean quantities on grids
 !!
 !! IN
 !!  particule(particletype)=particule vector
 !!  Efield,Bfield(arr3Dtype)=Electric and Magnetic fields
 !!  nc1(3)=Bfield points
 !!  nptot=total number of particles (for proc)
 !!  x_min=lower bound of x-axis simulation box
 !!  x_max=upper bound of x-axis simulation box
 !!  infompi(mpitype)=mpi informations
 !!  
 !! NOTE
 !!  Normally this subroutine is not called. It is usefull only
 !!  for complicated debugs
 subroutine controls(particule,Efield,Bfield,Spe,nc1,nptot,x_min,x_max,infompi)

  use defs_particletype,only     : particletype
  use defs_species

  integer,intent(in) :: nptot
  real(dp),intent(in) :: x_min,x_max
  integer,intent(in) :: nc1(3)
  type(mpitype),intent(in) :: infompi
  type(species_type),intent(in) :: Spe
  type(arr3Dtype),intent(in) :: Bfield,Efield
  type(particletype),intent(in) :: particule(:)

  integer :: ioerr,rindex
  integer :: ntranche,itranche,npt_tranche,npt_tranche_tot
  real(dp) :: rtranche,min_tranche,max_tranche,norm
  real(dp) :: vmeantranche,vmeantranche_tot
  real(dp) :: Bmeantranche(3),Emeantranche(3)
  real(dp) :: Emeantranche_tot(3),Bmeantranche_tot(3)
  character(len=200) :: msg
  logical,allocatable :: mask(:)

  !--Allocate the array containing the particles in the tranche
  allocate(mask(nptot))

  !--Number of tranches
  ntranche = 10
  rtranche = (x_max-x_min)/real(ntranche,dp)
  rindex = nc1(1)/ntranche

  write(msg,'(4a,a3,a10,7a18)') ch10," ______ Control tranches _____",&
       &ch10,"           ","n","npt","<vx>",&
       &      "<Ex>","<Ey>","<Ez>",&
       &      "<Bx>","<By>","<Bz>"
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  !--Loop on traches
  do itranche = 1, ntranche
   !--Bounds of the tranche
   min_tranche = x_min+real(itranche-1,dp)*rtranche
   max_tranche = min_tranche+rtranche

   !--Mask selecting the particles in the tranche (1 proc)
   mask = (particule(:nptot)%pos(1)> min_tranche .and. particule(:nptot)%pos(1)<=max_tranche)

   !--Number of particles in the tranche (1 proc)
   npt_tranche = count(mask)
   
   !--Sum of speeds (1 proc)
   vmeantranche = sum(particule(:nptot)%vel(1),mask=mask)

   !--Sum to proc 0 the number of points in the tranche and the sum of
   !  speeds of these
   call MPI_REDUCE(npt_tranche,npt_tranche_tot,1,MPI_INTEGER,MPI_SUM,0,infompi%comm,ioerr)
   call MPI_REDUCE(vmeantranche,vmeantranche_tot,1,QP_MPI_DP,MPI_SUM,0,infompi%comm,ioerr)

   !--The mean of speed on all procs in the tranche
   if(infompi%me == 0) vmeantranche_tot = vmeantranche_tot/real(npt_tranche_tot,dp)

   !--Sum Electric Field
   emeantranche(1) = sum(Efield%x(rindex*(itranche-1):rindex*itranche-1,:,:))
   emeantranche(2) = sum(Efield%y(rindex*(itranche-1):rindex*itranche-1,:,:))
   emeantranche(3) = sum(Efield%z(rindex*(itranche-1):rindex*itranche-1,:,:))
   call MPI_REDUCE(emeantranche,emeantranche_tot,3,QP_MPI_DP,MPI_SUM,0,infompi%comm,ioerr)   
   if(infompi%me == 0)  emeantranche_tot = emeantranche_tot/&
        &               real(rindex*(nc1(2)+1)*(nc1(3)+1)*infompi%nproc,dp)

   !--Sum Electric Field
   bmeantranche(1) = sum(Bfield%x(rindex*(itranche-1):rindex*itranche-1,:nc1(2),:nc1(3)))
   bmeantranche(2) = sum(Bfield%y(rindex*(itranche-1):rindex*itranche-1,:nc1(2),:nc1(3)))
   bmeantranche(3) = sum(Bfield%z(rindex*(itranche-1):rindex*itranche-1,:nc1(2),:nc1(3)))
   call MPI_REDUCE(bmeantranche,bmeantranche_tot,3,QP_MPI_DP,MPI_SUM,0,infompi%comm,ioerr)   
     if(infompi%me == 0)  bmeantranche_tot = bmeantranche_tot/&
        &               real(rindex*nc1(2)*nc1(3)*infompi%nproc,dp)

   !--Write control on tranches
   write(msg,'(a,i3,i10,7f18.6)')&
        & "   tranche ",itranche,npt_tranche_tot,vmeantranche_tot,&
        &               emeantranche_tot,bmeantranche_tot
   call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  enddo

  deallocate(mask)

  norm = real(sum(Spe%S(:)%ng),dp)/real(nptot,dp)

  write(msg,'(2a,2(a,a12,f12.5))')&
       & ch10," _____ Mass and charge _____",&
       & ch10, "   Mass   = ",sum(particule(:nptot)%mass)*norm,& 
       & ch10, "   Charge = ",sum(particule(:nptot)%char)*norm
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

 end subroutine controls

end module diag_energy

