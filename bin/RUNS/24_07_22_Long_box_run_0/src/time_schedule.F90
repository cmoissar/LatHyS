!!=============================================================
!!=============================================================
module time_schedule

 use defs_basis
 use defs_mpitype
 use defs_grid
 use part_moment,only                : momt3d,momad
 use particle            !all
 use particle_creation,only     : new_particles
 use particle_com
 use defs_parametre,only             : dt,nsub,ntest,ipe,idisp,t_init_dip
 use defs_variable,only : by0,bz0,by1,bz1,vxmean
 use m_timing,only              : time_get
 use m_writeout
 use time_variation
 use environment,only : set_time_variation
#include "q-p_common.h"

 implicit none 
 private

 public::        &
      first,     &
      cam3,      &
      last

contains
 !!#####################################################################


 !****************************** FIRST *********************************************
 subroutine first()
  !--Half step of x; moment collection; set Bh = B

  use field_lissage,only       : smth_func
  use defs_variable,only      : iter,nptot,&
       &                   iwr,particule,&
       &                   Bfield,Bfield_h,&
       &                   dnf,dna,pe,&
       &                   s_min_loc,s_max_loc,&
       &                   Preg,index_exit,s_r,s_min,&
       &                   Spe
 use particle_sort

  real(dp) :: dth
  character(len=100) :: msg

  __WRT_DEBUG_IN("first")
  __GETTIME(5,1)!--Timer start

  !--Check variables
  if(nsub<=1)then
   write(msg,'(a)')' WARNING : nsub <= 1 : set to 2'
   call wrt_double(qp_out,msg,0,wrtdisk)
   nsub = 2
  else if(modulo(nsub,2)/=0)then
   write(msg,'(a)')' WARNING : nsub odd : set even'
   call wrt_double(qp_out,msg,0,wrtdisk)
   nsub = nsub + 1
  end if
  
  !--Collect dna,ua
  call momt3d(1)

  !--Half-step of position
  dth = half*dt

  call xcalc3(dth,particule,1,nptot,s_min_loc,s_max_loc, &
       &Preg,index_exit,Spe)

  !--echange des informations
  call pack_com(particule,Preg,index_exit,s_min,s_r,nptot,s_min_loc,s_max_loc)

  !--exit and creation of particles
  call new_particles(1,particule)
  
  !--Collect dn,u,dnf,uf
  call momt3d(2)

  call p_sort(particule,nptot)
  
  !--Smoothing (optional)
  if(iwr(6)==1) call smth_func(dnf,nc,mpiinfo)

  !--Set Bh = B
  Bfield_h= Bfield
 
  __GETTIME(5,2)!--Timer stop
  __WRT_DEBUG_OUT("first")
 end subroutine first
 !********************************* END FIRST ****************************************

 !********************************** CAM3 *******************************************
 subroutine cam3(tmax)
!  use diag_energy,only   : energy_proc
  use defs_parametre,only : planetname,gstep,pos_plan
  use particle,only      : move 
  use field_b,only       : testBfield
  use field,only         : calc_field
  use defs_variable,only : iter,nptot,&
       &                   iwr,particule,&
       &                   Bfield,Bfield_h,vela,&
       &                   dnf,dna,pe,&
       &                   diag,diag_pickup,&
       &                   s_min_loc,s_max_loc,&
       &                   t,diag_info,&
       &                   Preg,index_exit,Bfield_0,&
       &                   Spe, vel
  use field_lissage

  real(dp),intent(in) :: tmax
  integer :: nn
  character(len=200) :: msg
  real(dp)           :: alpha !multiplying factor to initialize step by step the dipolar field (Bfield_0)

  __WRT_DEBUG_IN("cam3")
  __GETTIME(6,1)!--Timer start

  !--The system will be advanced from T to TMAX, with a first and last step.
  !--Every NTEST steps, average the two B solutions if they have separated.
  
  alpha = (t_init_dip+0._dp)**(1._dp/(t_init_dip+1._dp))

  do while(t < (tmax-0.01_dp*dt))    
   if (iter.le.t_init_dip) then 
!        Bfield%x=Bfield%x+Bfield_0%x*0.0466513_dp
!        Bfield%y=Bfield%y+Bfield_0%y*0.0466513_dp
!        Bfield%z=Bfield%z+Bfield_0%z*0.0466513_dp
! 	Bfield_0%x=Bfield_0%x*1.0466513_dp
!        Bfield_0%y=Bfield_0%y*1.0466513_dp
!        Bfield_0%z=Bfield_0%z*1.0466513_dp
        Bfield%x=Bfield%x+Bfield_0%x*(alpha-one)
        Bfield%y=Bfield%y+Bfield_0%y*(alpha-one)
        Bfield%z=Bfield%z+Bfield_0%z*(alpha-one)
        Bfield_0%x=Bfield_0%x*alpha
        Bfield_0%y=Bfield_0%y*alpha
        Bfield_0%z=Bfield_0%z*alpha
   endif
   iter = iter+1
   t = real(iter,dp)*dt
   write(msg,'(64a,i6,a,f8.3,63a)')ch10,&
        &" ", ('=',nn=1,60),ch10,&
        &"   ITERATION: ",iter,"     time: ",t, ch10,& 
        &" ", ('=',nn=1,60),ch10
   call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

   !--Step of length dt, advancing x, v and B
   !--(1)................ B(0) to B(1/2) using dna(0)
   call calc_field(mpiinfo)

   !--(2)................ Current half step
   call momad()

   !--(3)................ Particle step and moments
   call move(Preg,index_exit)

   !--(4)................ B(1/2) to B(1) using dna(1)
   call calc_field(mpiinfo)

!   call energy_proc(Preg,diag_info)

   if(mod(iter, 40*ntest) == 0) then
    call testBfield(Bfield,Bfield_h,t,iter,mpiinfo)
    !--Smoothing Bfield in the last part of the box (along the x-axis)
    call smth_func(Bfield,nc,mpiinfo)
    call smth_func(Bfield_h,nc,mpiinfo)
    !--write magnetic field at specific location
    !call write_temporal_B()
   endif

   diag(iter+1) = nptot
   diag_pickup(iter+1) = count(particule(1:nptot)%exc/=zero)
   
   ! re-initiliase plasma incident info if temporal variation are loaded
   call set_time_variation(by0,bz0,by1,bz1,vxmean,Spe,iter)

  enddo
  
  __GETTIME(6,2)!--Timer stop
  __WRT_DEBUG_OUT("cam3")



 end subroutine cam3

 subroutine write_temporal_B()

  use defs_variable, only : Bfield, vel, dn
  use defs_variable, only : s_min_loc, s_max_loc, t
  use defs_parametre, only : pos_plan, gstep
  

  real :: x_gsm, y_gsm, z_gsm
  real :: x_sim, y_sim, z_sim
  character(len=150):: file_name, FN
  integer :: ix_box, iy_box, iz_box
  logical :: L_open
  integer :: file_count

  file_count = 0  

 ! write(*,*) "pos_plan:", pos_plan
 ! write(*,*) "s_min_loc(2):", s_min_loc(2)

  do ix_box = 1, int((s_max_loc(1) - s_min_loc(1))/gstep(1))
    do iy_box = 1, int((s_max_loc(2) - s_min_loc(2))/gstep(2))
      do iz_box = 1, int((s_max_loc(3) - s_min_loc(3))/gstep(3))

        x_sim = s_min_loc(1) + ix_box*gstep(1)
        y_sim = s_min_loc(2) + iy_box*gstep(2)
        z_sim = s_min_loc(3) + iz_box*gstep(3)

!        x_sim = s_min_loc(1) + (ix_box-1)*gstep(1)
!        y_sim = s_min_loc(2) + (iy_box-1)*gstep(2)
!        z_sim = s_min_loc(3) + (iz_box-1)*gstep(3)

        x_gsm =  pos_plan(1)-x_sim
        y_gsm =  pos_plan(2)-y_sim
        z_gsm = -pos_plan(3)+z_sim

        if (((sqrt(x_gsm**2+y_gsm**2+z_gsm**2)<120) .and. &
               (mod(int(x_gsm), 10)==0 .and. &
                mod(int(y_gsm), 25)==0 .and. &
                mod(int(z_gsm), 25)==0)) .or. &
             (mod(int(x_gsm), 750)==0 .and. &
              mod(int(y_gsm), 100)==0 .and. &
              mod(int(z_gsm), 100)==0)) then
           
            write(file_name, '(A6,I5,A2,I5,A2,I5,A2,A4)') &
                              "Btime_",int(x_gsm),"x_", &
                                       int(y_gsm),"y_", &
                                       int(z_gsm),"z_",'.txt'
            !--By creating file_name and file_count at the same time
            !--we ensure they have a one-to-one correspondancy
            file_name = trim(file_name)
            file_count = file_count+1
            open(UNIT=80+file_count,FILE=file_name, &
                & status="old", position="append", action="write")

            inquire(UNIT=80+file_count, NAME=FN, OPENED=L_open)

!            write(*,*)  "unit", 80+file_count, &
!                      & "position", int(x_gsm), int(y_gsm), int(z_gsm), &
!                      & "isopen", L_open,"has name", FN

            write(UNIT=80+file_count,FMT=*) 'time', t,  &
                              & 'B_field%xyz', &
                              & Bfield%x(ix_box,iy_box,iz_box), &
                              & Bfield%y(ix_box,iy_box,iz_box), &
                              & Bfield%z(ix_box,iy_box,iz_box), &
                              & 'velocity%xyz', &
                              & vel%x(ix_box,iy_box,iz_box), &
                              & vel%y(ix_box,iy_box,iz_box), &
                              & vel%z(ix_box,iy_box,iz_box), &
                              & 'density', &
                              & dn(ix_box,iy_box,iz_box)
        endif

      enddo
    enddo
  enddo
 end subroutine write_temporal_B

 !***************************** END CAM3 ********************************************

 !******************************* LAST **********************************************
 subroutine last()
  !--Backward half-step of x; average B and Bh; moment collection; E-field

  use defs_arr3dtype
  use field_pe
  use field_e,only       : ecalc3
  use defs_parametre,only :gstep,dmin,resis
  use defs_variable,only : iter,nptot,index_exit,&
       &                   iwr,particule,&
       &                   Bfield,Bfield_h,Efield,vela,&
       &                   resistivity,dnf,dna,pe,&
       &                   s_min_loc,s_max_loc,&
       &                   e_conv,e1_conv,&
       &                   rmu0,te,t,&
       &                   Preg,index_exit,s_r,s_min,&
       &                   Spe

  real(dp) :: dthm

  __WRT_DEBUG_IN("last")
  __GETTIME(13,1)!--Timer start

  dthm = half*dt

  call xcalc3(dthm,particule,1,nptot,s_min_loc,s_max_loc, &
       &      Preg,index_exit,Spe)

  !--echange des informations 
  call pack_com(particule,Preg,index_exit,s_min,s_r,nptot,s_min_loc,s_max_loc)

  !--exit and creation of new particles
  call new_particles(1,particule)

  call momt3d(3)

  Bfield%x = half*(Bfield_h%x+Bfield%x)
  Bfield%y = half*(Bfield_h%y+Bfield%y)
  Bfield%z = half*(Bfield_h%z+Bfield%z)

  call pecalc(ipe,te,dna,pe)

  call ecalc3(Efield,Bfield,vela,dna,pe,resistivity,e_conv,&
       &e1_conv,gstep,dmin,rmu0,resis,idisp,nc1,mpiinfo)

  __GETTIME(13,2)!--Timer stop
  __WRT_DEBUG_OUT("last")
 end subroutine last
 !********************************* END LAST
 !******************************************

end module time_schedule
