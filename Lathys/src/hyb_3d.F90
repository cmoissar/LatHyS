program Mars3d_f90
 ! =================
 
 !*********************************************************************************************
 !
 ! Code numerique hybride s'appyant sur l'algorithme CAM-CL d'Alan Matthews ( Matthews A.,
 ! Comp.Jour. Physics, 1994)
 !
 ! Ce code permet de decrire l'interaction entre le vent solaire et une atmosphere planetaire
 ! pour ce code la planete est Mars.
 !
 ! Code ecrit en fortran 90 par Ronan Modolo  le 20 février 2003
 !
 !********************************************************************************************

 use defs_basis
 use defs_parametre
 use defs_mpitype
 use defs_variable
 use defs_grid
 use environment
 use initialisation  !all
 use diag_energy,only   : energy_proc
 use diagnostique    !all
 use part_moment,only             : Dmtsp3,Emtsp3
 use time_schedule   !all
 use m_restart       !all
 use m_cmdline
 use m_logo
 use m_writeout
 use m_timing
 use diag_impex_xml

 
 implicit none
#include "q-p_common.h"

 integer :: ireg     !--iteration pour le nombre de diagnostique
 integer :: ioerr,ii,restart_opt,iunits
 character(len=500) :: msg,name_file_imp

 __WRT_DEBUG_IN("quiet_plasma")

 !--Initialisation of MPI, and variables relative to parallelisation
 call init_all_mpiinfo(mpiinfo,mpiinfo_pick)

 !--Timers initialisation
 call time_init
 call time_get(1,1)

 !--Get command line
 call cmd_line(restart_opt)

 !--Define Output file
 if(wrtdisk==0.and.mpiinfo%me==0) &
      &  open(unit=qp_out,file = trim(qp_out_name),form = "formatted",status = "unknown")

 !--Computation of the cell size (if to change)
 call ncell_compute(ncx0,ncy0,ncz0)

 !--Print log of quiet_plasma
 call logo(qp_out,date_ok=1)

 !--Select the Planet Environment
 call select_environment(planetname,ns)

 !--Print MPI information
 call print_mpiinfo(qp_out)

 !--Initialization of files name, for output is the DATE OF DAY
 if(trim(fildat)=="") fildat = print_date()

 !--Allocation of arrays and global variables
 call allocation()

 !--Initialisation of 3D data
 Select case(restart)
 case(0)
  call h3init()
 case(1)
  call re_start(restart_opt)
 case default
  stop "Restart Option: Not yet implemented"
 end Select

 if(mpiinfo%me==0) call write_impex_xml()

 !--Check time
 !--The system will be advanced from T to TMAX = TREG(IREG)
 if(r_tm%nreg<=0) stop " No time advance nreg <= 0"
 !--Check size array
 if(product(ncm_tot) > nxyzm) stop "ERROR : GRID TOO LARGE"
 if(Spe%S(ns)%n2 > npm) stop "ERROR : TOO MANY PARTICLES"
 
 
  !--------------fichiers particules impactantes-------------
 !diag_part_imp = 1
 !if (diag_part_imp.ne.0) then
 !  iunits = mpiinfo%me+20
 ! write(name_file_imp,'(a19,i4.4,a1,a)')"part_imp_formatted_",mpiinfo%me,"_",trim(fildat)
 !
 !   open(UNIT   = iunits,&
 !        FILE   = name_file_imp,&
 !        FORM   = 'formatted',&
 !        ACTION = 'write', &
 !        STATUS = 'unknown')
 ! print*,'iunits open',iunits
 !endif
 !----------------------------------------------------------


 if(restart==0) call first()

 allocate(diag(nhm+1),diag_pickup(nhm+1))
 diag = 0;   diag(1) = mpiinfo%me
 diag_pickup = 0; diag_pickup(1) = mpiinfo%me

 do ireg = 1,r_tm%nreg
  !--CAM3 advances the system to the next data dump
  write(msg,'(64a,i4,2(a,f8.2),62a)')ch10,&
       &" ",('=',ii=1,60),ch10,&
       & "   Data dump:", ireg,&
       & "    from t= ",t,&
       & "    to tmax= ", r_tm%time(ireg),ch10,&
       &" ",('=',ii=1,60)
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  call cam3(r_tm%time(ireg))

  call energy_proc(Preg,diag_info)

  if(wrtdisk == 0) then
   !--Particle data is selectively written out
   !iwr(12) = ipfile(ireg)
   call Dmtsp3(particule,dn,vel,1,nptot, &
        &      gstep,s_min_loc,&
        &      mpiinfo,&
        &      nc1)

   call diag_all(r_tm%file(ireg))

   !--Restart record
   call wrt_restart()

   call Emtsp3(particule,dn,vel,1,nptot, &
        &      gstep,s_min_loc,&
        &      mpiinfo,&
        &      nc1)
  endif
 enddo

 deallocate(diag,diag_pickup)

 call last()


!if (diag_part_imp.ne.0) then
! close(iunits)
!endif


 call deallocation()
 call time_get(1,2)
 call time_clean()
 call MPI_FINALIZE(ioerr)
 __WRT_DEBUG_OUT("quiet_plasma")
end program Mars3d_f90
