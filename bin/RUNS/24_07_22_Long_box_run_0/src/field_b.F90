!!=============================================================
!!=============================================================
!!module: m_bcalc
!! NAME
!!  m_bcalc (MMancini,RModolo)
!!
!! FUNCTION
!!  Contains the routine for the Efield computation
!! NOTE
module field_b

 use defs_basis
 use defs_arr3Dtype
 use defs_mpitype,only     : mpitype
 use m_writeout
 use m_timing,only         : time_get
 use mpi

#include "q-p_common.h"

 implicit none
 private

 public::           &
      bcalc3,       &!--Compute Magnetic Fields
      testBfield     !--Test differences between 2 solution of Bfield

contains
 !!#####################################################################

 !!=============================================================
 !!routine: m_bcalc/bcalc3
 !!
 !! FUNCTION
 !!  Compute Magnetic field
 !! IN 
 !!  dt=time step
 !!  by0,bz0=Magnetic Field at x=0
 !!  vxmean=velocity at x=0
 !!  nc1(3)=grid size +1
 !!  gstep(3)=Grid step
 !!  Efield(arr3Dtype)=Electric Field
 !!
 !! OUT
 !!  ey_conv,ez_conv,e1y_conv,e1z_conv=Electric Field at x=0,x=1
 !! SIDE EFFECT
 !!  Bfield(arr3Dtype)=Magntic Field
 !! NOTES
 !!  Calculation of B from Maxwell's equation dB/dt = - curl(E)
 !!  B(new) = B - dtb * curl(E)
 subroutine bcalc3(dt,Bfield,Efield,nc1,gstep,&
      &            by0,bz0,by1,bz1,ey_conv,ez_conv,&
      &            e1y_conv,e1z_conv,vxmean)

  use field_cond_limit,only       : cond_limit_func
  use field_lissage,only          : smth_func
  use defs_mpitype
#ifdef HAVE_WAVE_TEST
  use m_wavetest,only       : b_boundary_wt
#endif

  integer,intent(in) :: nc1(3)
  real(dp),intent(in) :: dt
  real(dp),intent(in) :: by0,bz0,by1,bz1,vxmean
  real(dp),intent(out) :: ey_conv,ez_conv,e1y_conv,e1z_conv
  real(dp),dimension(3),intent(in) :: gstep  !--Grid step
  type(arr3Dtype),intent(in) :: Efield
  type(arr3Dtype),intent(inout) :: Bfield

  real(dp),dimension(3) :: dt_gstep_invq

  __WRT_DEBUG_IN("bcalc3")
  __GETTIME(47,1)!--Timer start

  ! dt/(4*dx), dt/(4*dy) and dt/(4*dz) for differential operators
  dt_gstep_invq = (quarter)/gstep

  !--Bfield has size (:nc1(1)+1,:nc1(2)+1,:nc1(3)+1) but its significative values are
  !--given by (:nc1(1),:nc1(2),:nc1(3)). For electric Field all (:nc1(1)+1,:nc1(2)+1,:nc1(3)+1)
  !--are significative!!!!!!!!!!!!!!!
  Bfield%x(1:nc1(1),1:nc1(2),1:nc1(3)) = Bfield%x(1:nc1(1),1:nc1(2),1:nc1(3)) - dt*(&       
       + dt_gstep_invq(2)*(&
       &   - Efield%z(1:nc1(1),1:nc1(2),1:nc1(3)) &
       &   - Efield%z(2:      ,1:nc1(2),1:nc1(3)) &
       &   + Efield%z(1:nc1(1),2:      ,1:nc1(3)) &
       &   + Efield%z(2:      ,2:      ,1:nc1(3)) &
       &   - Efield%z(1:nc1(1),1:nc1(2),2:      ) &
       &   - Efield%z(2:      ,1:nc1(2),2:      ) &
       &   + Efield%z(1:nc1(1),2:      ,2:      ) &
       &   + Efield%z(2:      ,2:      ,2:      ))&
       -dt_gstep_invq(3)*(&                  
       &   - Efield%y(1:nc1(1),1:nc1(2),1:nc1(3)) &
       &   - Efield%y(2:      ,1:nc1(2),1:nc1(3)) &
       &   - Efield%y(1:nc1(1),2:      ,1:nc1(3)) &
       &   - Efield%y(2:      ,2:      ,1:nc1(3)) &
       &   + Efield%y(1:nc1(1),1:nc1(2),2:      ) &
       &   + Efield%y(2:      ,1:nc1(2),2:      ) &
       &   + Efield%y(1:nc1(1),2:      ,2:      ) &
       &   + Efield%y(2:      ,2:      ,2:      ))&
       )

  Bfield%y(1:nc1(1),1:nc1(2),1:nc1(3)) = Bfield%y(1:nc1(1),1:nc1(2),1:nc1(3)) - dt*(&
       + dt_gstep_invq(3)*(&
       &   - Efield%x(1:nc1(1),1:nc1(2),1:nc1(3)) &
       &   - Efield%x(2:      ,1:nc1(2),1:nc1(3)) &
       &   - Efield%x(1:nc1(1),2:      ,1:nc1(3)) &
       &   - Efield%x(2:      ,2:      ,1:nc1(3)) &
       &   + Efield%x(1:nc1(1),1:nc1(2),2:      ) &
       &   + Efield%x(2:      ,1:nc1(2),2:      ) &
       &   + Efield%x(1:nc1(1),2:      ,2:      ) &
       &   + Efield%x(2:      ,2:      ,2:      ))&
       - dt_gstep_invq(1)*(&                     
       &   - Efield%z(1:nc1(1),1:nc1(2),1:nc1(3)) &
       &   + Efield%z(2:      ,1:nc1(2),1:nc1(3)) &
       &   - Efield%z(1:nc1(1),2:      ,1:nc1(3)) &
       &   + Efield%z(2:      ,2:      ,1:nc1(3)) &
       &   - Efield%z(1:nc1(1),1:nc1(2),2:      ) &
       &   + Efield%z(2:      ,1:nc1(2),2:      ) &
       &   - Efield%z(1:nc1(1),2:      ,2:      ) &
       &   + Efield%z(2:      ,2:      ,2:      ))&
       )

  Bfield%z(1:nc1(1),1:nc1(2),1:nc1(3)) = Bfield%z(1:nc1(1),1:nc1(2),1:nc1(3)) - dt*(& 
       + dt_gstep_invq(1)*(&
       &   - Efield%y(1:nc1(1),1:nc1(2),1:nc1(3)) &    
       &   + Efield%y(2:      ,1:nc1(2),1:nc1(3)) &
       &   - Efield%y(1:nc1(1),2:      ,1:nc1(3)) &
       &   + Efield%y(2:      ,2:      ,1:nc1(3)) &
       &   - Efield%y(1:nc1(1),1:nc1(2),2:      ) &
       &   + Efield%y(2:      ,1:nc1(2),2:      ) &
       &   - Efield%y(1:nc1(1),2:      ,2:      ) &
       &   + Efield%y(2:      ,2:      ,2:      ))&
       - dt_gstep_invq(2)*(&                       
       &   - Efield%x(1:nc1(1),1:nc1(2),1:nc1(3)) &
       &   - Efield%x(2:      ,1:nc1(2),1:nc1(3)) &
       &   + Efield%x(1:nc1(1),2:      ,1:nc1(3)) &
       &   + Efield%x(2:      ,2:      ,1:nc1(3)) &
       &   - Efield%x(1:nc1(1),1:nc1(2),2:      ) &
       &   - Efield%x(2:      ,1:nc1(2),2:      ) &
       &   + Efield%x(1:nc1(1),2:      ,2:      ) &
       &   + Efield%x(2:      ,2:      ,2:      ))&
       )

#ifdef HAVE_WAVE_TEST
  call b_boundary_wt(Bfield,nc1(1),ey_conv,ez_conv,e1y_conv,e1z_conv,vxmean)
#else
  call cond_limit_func(Bfield,nc1(1),by0,bz0,by1,bz1,ey_conv,ez_conv,e1y_conv,e1z_conv,vxmean)
#endif

  !--Smoothing Bfield in the last part of the box (along the x-axis)
  !call smth_func(Bfield,nc1-1,infompi)

  __GETTIME(47,2)!--Timer stop
  __WRT_DEBUG_OUT("bcalc3")
 end subroutine bcalc3



 !!=============================================================
 !!routine: m_bcalc/testBfield
 !!
 !! FUNCTION
 !!  Tests the system to see whether the two solutions have separated.
 !!  If so, set ICLAVE = 0
 !!
 !! IN 
 !!  t=time step
 !!  iter=iteration
 !!  infompi(mpitype)=mpi informations
 !!
 !! OUT
 !! SIDE EFFECT
 !!  Bfield(arr3Dtype)=Magntic Field 1
 !!  Bfield_h(arr3Dtype)=Magntic Field 2
 !! NOTES
 subroutine testBfield(Bfield,Bfield_h,t,iter,infompi)

  use defs_parametre,only      : eps
  !use mpi

  integer,intent(in) :: iter
  real(dp),intent(in) :: t
  type(mpitype),intent(in) :: infompi
  type(arr3Dtype),intent(inout) :: Bfield_h,Bfield

  integer :: iclave,ii,jj,kk,ioerr
  real(dp) :: err,err_tout_proc
  real(dp) :: b_a_tot,e_b_tot
  integer :: Bbnd(3)
  real(dp) :: rtmp(3)
  character(len=200) :: msg

  __WRT_DEBUG_IN("test")
  __GETTIME(50,1)!--Timer start

  Bbnd = ubound(Bfield%x)

  err = zero
  do kk=1,Bbnd(3)-1
   do jj=1,Bbnd(2)-1
    do ii=1,Bbnd(1)-1
     rtmp(1) = (Bfield%x(ii,jj,kk) + Bfield_h%x(ii,jj,kk))
     rtmp(2) = (Bfield%y(ii,jj,kk) + Bfield_h%y(ii,jj,kk))
     rtmp(3) = (Bfield%z(ii,jj,kk) + Bfield_h%z(ii,jj,kk))
     b_a_tot = dot_product(rtmp,rtmp)!*quarter
     if(b_a_tot< tol10) cycle
     rtmp(1) = (Bfield%x(ii,jj,kk) - Bfield_h%x(ii,jj,kk))
     rtmp(2) = (Bfield%y(ii,jj,kk) - Bfield_h%y(ii,jj,kk))
     rtmp(3) = (Bfield%z(ii,jj,kk) - Bfield_h%z(ii,jj,kk))
     e_b_tot = dot_product(rtmp,rtmp)
     err = max(err,e_b_tot/b_a_tot)
    enddo
   enddo
  enddo
  err = sqrt(err)

  !--On cherhce le maximum de err sur tous les processus
  call MPI_ALLREDUCE(err,err_tout_proc,1,QP_MPI_DP,MPI_MAX,infompi%comm,ioerr)

  !--Si le max des "err" est suprieurs  la precision alors tous les processus
  !  on fait la moyenne entre la copie t et le champ magntique
  iclave = 0
  if(err_tout_proc > eps) then
   iclave = 1

   !--Average of Bh and B
   Bfield%x = half*(Bfield_h%x+Bfield%x)
   Bfield%y = half*(Bfield_h%y+Bfield%y)
   Bfield%z = half*(Bfield_h%z+Bfield%z)

   !--Set Bh = B
   Bfield_h = Bfield

   write(msg,'(a5,i8,a3,f10.2,a7,e16.7,a12)')&
        &     'CAM3',iter,' t=',t,' err = ',err_tout_proc
   call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
  end if

  __GETTIME(50,2)!--Timer stop
  __WRT_DEBUG_OUT("test")
 end subroutine testBfield
end module field_b
