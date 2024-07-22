!!=============================================================
!!=============================================================
!!module: m_grid
!! NAME
!!  cration (RModolo,MMancini)
!!
!! FUNCTION
!!  Contains variables concernig the spatial grid of simulation
!!  Starting from the total grid the grid on any proc are calculated
!!  taking into account the number of procs.
!!  
module defs_grid
 use defs_basis
#include "q-p_common.h"

 implicit none
 private

 !--------------- GRID --------------------

 !--Number of cells in the three dimensions X,Y,Z
 integer,public,save,protected :: nc_tot(3)

 !--N of Points (internal) in any dimension (global)
 integer,public,save,protected :: nc1_tot(3)   

 !--N of Points (external) in any dimension (Electric Field) (global)
 integer,public,save,protected :: ncm_tot(3) != nc_tot+2   

 !--Max number of grid points for the Fields
 integer,public,save,protected :: nxyzm

 !--N Cells for any sub-region (proc)
 integer,public,save,protected :: nc(3)

 !--N (interal) Points  to any sub-region (proc)
 integer,public,save,protected :: nc1(3)

 !--N (external) Points (electric field)  to any sub-region (proc)
 integer,public,save,protected :: ncm(3)

 integer,public,save,protected :: ncross  !--Nombre maximum de macroparticule qui touchent un bord

 public ::              &
      set_grid,         &
      ncell_compute      !--Compute the grid and cell size (Total and Procs)
contains
 !!#####################################################################

 subroutine set_grid(itmp,opt)
  integer,intent(in) :: itmp(3)
  integer,intent(in) :: opt
  select case(opt)
  case(1)
   nc1_tot = itmp(:)
  case(2)
   ncm(:) = itmp(:)
  case(3)
   ncm_tot = itmp(:)
  case(4)
   nxyzm = itmp(1) 
  case(5)
   nc_tot = itmp(:)
   ! case(7)
   !  ncross = itmp(1);
  end select
 end subroutine set_grid

 !!=============================================================
 !!routine: m_grid/ncell_compute
 !!
 !! FUNCTION
 !! Compute the grid and cell sizes (Total and per Proc)
 !!         
 !! IN
 !! OUT
 !!   nc1_tot(3)=The grid where fields (Excepted Efield) are defined (total)
 !!   ncm_tot(3)=The grid where the fields are allocated (total)
 !!   nc1(3)=The grid where fields (Excepted Efield) are defined (proc)
 !!   ncm(3)=The grid where the fields are allocated are defined (proc)
 !!   nc(3)=Number of cells defining the space of simulation (proc)
 !!   nxyzm=total (scalar) size of the field allocation grid
 !!   ncross=.......
 !!   npm=Max number of particles per proc
 !! SIDE EFFECT
 !!   nc_tot(3)=Number of cells defining the space of simulation
 !!    (total) could be changed (only nc_tot(2:3) if not compatible with
 !!    MPI Cartesian grid
 subroutine ncell_compute(ncx0,ncy0,ncz0)
#ifndef NOTHAVE_MPI
  use m_writeout
  use defs_parametre,only    : npm,n_part_max,ROI,gstep
  use defs_mpitype
  use defs_variable
#endif
  integer,intent(in) :: ncx0,ncy0,ncz0

#ifndef NOTHAVE_MPI
  character(len=500) :: msg
  real(dp),dimension(ncy0+100)::ct_y
  real(dp),dimension(ncz0+100)::ct_z
  integer,dimension(mpiinfo%dims(1))::cl_y
  integer,dimension(mpiinfo%dims(2))::cl_z
  integer ::ii,jj,w0,w1,ct0y,ct0z
  !/////////////////////////////////////////////////////////////////////
  !//        PRIMORDIAL                                               //
  !//  il faut que le nombre de processus en X soit un multiple du    //
  !//  du nombre de cellule en X                                      //
  !//  il faut que le nombre de processus en Y soit un multiple du    //
  !//  du nombre de cellule en Y                                      //
  !/////////////////////////////////////////////////////////////////////

  __WRT_DEBUG_IN("ncell_compute")

  !--Setting cell number equal to defined in parameter
  nc_tot = (/ ncx0,ncy0,ncz0 /)

  !--Control if the number of cells is compatible with the cartesian grid
  if (mod(nc_tot(2),mpiinfo%dims(1)) /= 0) then
   nc_tot(2) = int(nc_tot(2)/mpiinfo%dims(1))*mpiinfo%dims(1)+mpiinfo%dims(1)
   write(msg,'(a,i6)')' WARNING : Y NUMBER CELL CHANGE : ',nc_tot(2)
   call wrt_double(6,msg,wrtscreen,wrtdisk)
  endif

  if (mod(nc_tot(3),mpiinfo%dims(2)) /= 0) then
   nc_tot(3) = int(nc_tot(3)/mpiinfo%dims(2))*mpiinfo%dims(2)+mpiinfo%dims(2)
   write(msg,'(a,i6)')' WARNING : Z NUMBER CELL CHANGE : ',nc_tot(3) 
   call wrt_double(6,msg,wrtscreen,wrtdisk)
  endif
  ct_y(:)=0._dp
  ct_z(:)=0._dp
  ct_y(1:nc_tot(2))=1.
  ct_y(int(ROI%ymin*nc_tot(2)):int(ROI%ymax*nc_tot(2)))=ct_y(int(ROI%ymin*nc_tot(2)):int(ROI%ymax*nc_tot(2)))+ROI%excess
  ct_y(1:nc_tot(2))=real(mpiinfo%dims(1),dp)*ct_y(1:nc_tot(2))/(sum(ct_y(1:nc_tot(2)))+1._dp)
  do jj=2,nc_tot(2)
     ct_y(jj)=ct_y(jj)+ct_y(jj-1)
  enddo

  do ii=1,mpiinfo%dims(1)
    do jj=2,nc_tot(2)
            if (int(ct_y(jj)).lt.ii) cl_y(ii)=jj
    enddo
  enddo
  do ii=mpiinfo%dims(1),2,-1
     cl_y(ii)=cl_y(ii)-cl_y(ii-1)
  enddo

  ct0y=count(cl_y.eq.0)
  do ii=1,ct0y
   w0=minloc(cl_y,dim=1)
   if (cl_y(w0).eq.0) then
           w1=maxloc(cl_y,dim=1)
           if (cl_y(w1).le.1) then
             print *,mpiinfo%me,' has a problem, too much process for grid size? cly ',cl_y
             stop
           endif
           cl_y(w1)=cl_y(w1)-1
          cl_y(w0)=cl_y(w0)+1
   endif
  enddo

  ct_z(1:nc_tot(3))=1.
  ct_z(int(ROI%zmin*nc_tot(3)):int(ROI%zmax*nc_tot(3)))=ct_z(int(ROI%zmin*nc_tot(3)):int(ROI%zmax*nc_tot(3)))+ROI%excess
  ct_z(1:nc_tot(3))=real(mpiinfo%dims(2),dp)*ct_z(1:nc_tot(3))/(sum(ct_z(1:nc_tot(3)))+1._dp)
  do jj=2,nc_tot(3)
     ct_z(jj)=ct_z(jj)+ct_z(jj-1)
  enddo
  do ii=1,mpiinfo%dims(2)
    do jj=2,nc_tot(3)
      if (int(ct_z(jj)).lt.ii) cl_z(ii)=jj
    enddo
  enddo
  do ii=mpiinfo%dims(2),2,-1
    cl_z(ii)=cl_z(ii)-cl_z(ii-1)
  enddo

  ct0z=count(cl_z.eq.0)
  do ii=1,ct0z
   w0=minloc(cl_z,dim=1)
   if (cl_z(w0).eq.0) then
     w1=maxloc(cl_z,dim=1)
     if (cl_z(w1).le.1) then
      print *,mpiinfo%me,&
      ' has a problem, too much process for grid size? clz ',cl_z
      stop
     endif
     cl_z(w1)=cl_z(w1)-1
     cl_z(w0)=cl_z(w0)+1
   endif
  enddo
  !--On definit le nombre particule par cellule maximum
  npm = int(product(real(nc_tot,dp))*real(n_part_max,dp)/real(nproc,dp))
  !--Déclaration du nombre maximal de points de grille
  ncm_tot = nc_tot + 2
  nxyzm = product(ncm_tot)
  ncross = 2*32*int((ncm_tot(2)*ncm_tot(3) + ncm_tot(2)*ncm_tot(1) + ncm_tot(3)*ncm_tot(1))/nproc)
  !--Dimensions of "full" grid (BX,BY,BZ,DEN,UX,UY,UZ,PE)
  nc1_tot = nc_tot + 1
  !--Nombre de cellule en X et en Y pour chaque processeur
  nc(1) = nc_tot(1)
  nc(2) = cl_y(mpiinfo%coord(1)+1)!--int(nc_tot(2)/mpiinfo%dims(1))
  nc(3) = cl_z(mpiinfo%coord(2)+1)!--int(nc_tot(3)/mpiinfo%dims(2))
  ncm = nc + 2
  nc1 = nc + 1
#ifdef HAVE_DEBUG   
print *,'me: ',mpiinfo%me,' nc: ',nc
#endif

  s_min_loc(:) =0._dp
  s_max_loc(:) =0._dp
  if (mpiinfo%coord(1).gt.0) then
    do ii=1,mpiinfo%coord(1)
      s_min_loc(2)=s_min_loc(2)+real(cl_y(ii),dp)*gstep(2)
    enddo
  endif
  if (mpiinfo%coord(2).gt.0) then
    do ii=1,mpiinfo%coord(2)
      s_min_loc(3)=s_min_loc(3)+real(cl_z(ii),dp)*gstep(3)
    enddo
  endif
  s_max_loc=s_min_loc+real(nc,dp)*gstep
#ifdef HAVE_DEBUG   
print *,'me: ',mpiinfo%me,mpiinfo%coord(:),s_min_loc,s_max_loc
#endif
  __WRT_DEBUG_OUT("ncell_compute")
#endif !NOTHAVE_MPI
 end subroutine ncell_compute

end module defs_grid


