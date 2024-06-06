!=============================================================
!=============================================================
module field

 use defs_basis
 use defs_arr3Dtype
 use defs_mpitype,only     : mpitype
 use m_writeout
 use m_timing,only         : time_get
#include "q-p_common.h"

 implicit none
 private

 public::           &
      calc_field     !--Compute Electric and Magnetic Fields

contains
 !!#####################################################################


 !********************************** CALC_FIELD ******************************************
 subroutine calc_field(infompi)

  use defs_parametre,only : dt,nsub,idisp,ipe,resis,dmin,gstep
  use defs_variable,only  : Efield,Bfield,Bfield_h,vela,&
       &                    dna,pe,resistivity,&
       &                    e_conv,e1_conv,&
       &                    rmu0,te,by0,bz0,by1,bz1,vxmean
  use defs_grid,only      : nc1
  use field_pe
  use field_e,only        : ecalc3
  use field_b,only        : bcalc3

  type(mpitype),intent(in) :: infompi
  real(dp) :: dtb,hs
  integer  :: ms,ig

  __WRT_DEBUG_IN("calc_field")  
  __GETTIME(41,1)!--Timer start

  !----------------------------------------------------------------
  dtb = dt/real(nsub,dp)
  !----------------------------------------------------------------
  !--Advance B to t=t+dt/2
  !----------------------------------------------------------------
  !--En utilsant que la methode trilineaire
  !--Pour une utilisation NGP change ecalc3 en e3r
  ms=0
  ig = 1

  !--ajout du calcul%z de la pression electronique a t0 et t1
  call pecalc(ipe,te,dna,pe)

  do while(ms<=nsub)!  do ms=0,nsub
   hs = dtb
   if(ms==0 .or. ms==nsub) hs = hs*half
   if(ig == 1)then
    ig = 0

    __GETTIME(46,1)!--Timer start ecalc3 in calc_field 
    call ecalc3(Efield,Bfield,vela,&
         &      dna,pe,resistivity,&
         &      e_conv,e1_conv,gstep,&
         &      dmin,rmu0,resis,&
         &      idisp,nc1,&
         &      infompi)
    __GETTIME(46,2)!--Timer stop ecalc3 in calc_field 

    !if(iwr(7) == 1) call smth_func(Efield,nx1,ny1,nz1)

    call bcalc3(hs,Bfield_h,Efield,nc1,gstep,&
         &      by0,bz0,by1,bz1,e_conv(2),e_conv(3),&
         &      e1_conv(2),e1_conv(3),vxmean)

   else 
    ig = 1

    __GETTIME(46,1) !--Timer start ecalc3 in calc_field 
    call ecalc3(Efield,Bfield_h,vela,&
         &      dna,pe,resistivity,&
         &      e_conv,e1_conv,gstep, &
         &      dmin,rmu0,resis,&
         &      idisp,nc1,&
         &      infompi)
    __GETTIME(46,2) !--Timer stop ecalc3 in calc_field 

    !if(iwr(7) == 1) call smth_func(Efield,nx1,ny1,nz1)

    call bcalc3(hs,Bfield,Efield,nc1,gstep,&
         &      by0,bz0,by1,bz1,e_conv(2),e_conv(3),&
         &      e1_conv(2),e1_conv(3),vxmean)

   end if

   ms = ms+1
  enddo

  __GETTIME(41,2)!--Timer stop

  __WRT_DEBUG_OUT("calc_field")
 end subroutine calc_field
 !************************************ END CALC_FIELD *********************************


end module field
