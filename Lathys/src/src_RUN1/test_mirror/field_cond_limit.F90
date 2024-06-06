!!=============================================================
!!=============================================================
module field_cond_limit

 use defs_basis
 use mpi
 use m_writeout,only     : wrt_debug
 use m_timing,only       : time_get
#include "q-p_common.h"

 implicit none
 private

 integer,private,save :: planXY_type_D,planXZ_type_D,planXZ_type_E,planXY_type_E


 public ::                    &
      bound_init_mpi_planes,  &!--Inti. mpi_type_vector for bound. condition
      cond_limit_func,        &!--Generic interface for Boundary conditions
      bound_free_mpi_planes    !--Free mpi_type_vector

 private ::                   &
      mtpd3,                  &!--Bound. cond. for density vectors (cells)
      b_boundary,             &!--X_limits  for Magnetic fields
      apdh3d_arr3d,           &!--Bound. Cond. for Electric Field
      mtpd_four                !--Bound. Cond. for density and velocity vectors

 interface cond_limit_func
  module procedure mtpd3,mtpd_four,apdh3d_arr3d,b_boundary !,apdh3d
 end interface

contains
 !!#####################################################################

 !!=============================================================
 !!routine: cond_limit/bound_init_mpi_planes
 !!
 !! FUNCTION
 !!  Initialize planes for MPI sendrecv for cond_limit
 !!         
 !! IN
 !!  The bounds of the arrays
 !! OUT
 subroutine bound_init_mpi_planes(ncm)

  integer,intent(in) :: ncm(3)
  integer :: ioerr

  call MPI_TYPE_VECTOR(ncm(3)-1-1+1,ncm(1)-1-1+1,ncm(1)*ncm(2),QP_MPI_DP,planXZ_type_D,ioerr)
  call MPI_TYPE_COMMIT(planXZ_type_D,ioerr)
  call MPI_TYPE_VECTOR(ncm(2)-1-1+1,ncm(1)-1-1+1,ncm(1), QP_MPI_DP,planXY_type_D,ioerr)
  call MPI_TYPE_COMMIT(planXY_type_D,ioerr)
  call MPI_TYPE_VECTOR(ncm(3)-1-1-1+1,ncm(1)-1+1,ncm(1)*ncm(2),QP_MPI_DP,planXZ_type_E,ioerr)
  call MPI_TYPE_COMMIT(planXZ_type_E,ioerr)
  call MPI_TYPE_VECTOR(ncm(2)-1-1-1+1,ncm(1)-1+1,ncm(1),QP_MPI_DP,planXY_type_E,ioerr)
  call MPI_TYPE_COMMIT(planXY_type_E,ioerr)

 end subroutine bound_init_mpi_planes

 !!=============================================================
 !!routine: lissage/bound_free_mpi_planes
 !!
 !! FUNCTION
 !!  Free the types for MPI
 !!         
 !! IN
 !! OUT
 subroutine bound_free_mpi_planes()

  integer :: ioerr
  call MPI_TYPE_FREE(planXZ_type_D,ioerr)
  call MPI_TYPE_FREE(planXY_type_D,ioerr)
  call MPI_TYPE_FREE(planXZ_type_E,ioerr)
  call MPI_TYPE_FREE(planXY_type_E,ioerr)
 end subroutine bound_free_mpi_planes


 !!=============================================================
 !!routine: cond_limit/mtpd3
 !!
 !! FUNCTION
 !!   Periodic condition for densities
 !!
 !! INPUTS
 !!   infompi=type(mpitype) information about MPI
 !!   nc1(3)= dimension of the physical array (meanful)
 !!
 !! SIDE EFFECT
 !!  dn(nx+1,ny+1,nz+1)=array where put periodic conditions
 subroutine mtpd3(dn,infompi,nc1)

  use defs_mpitype

  integer,intent(in) :: nc1(3)
  type(mpitype),intent(in) :: infompi
  real(dp),dimension(:,:,:),intent(inout)::dn

  integer :: ioerr
  integer :: status(MPI_STATUS_SIZE)
!  real(dp),allocatable :: r_line(:)
  real(dp),allocatable :: r_plan(:,:)

  __WRT_DEBUG_IN("mtpd3")
  __GETTIME(63,1)!--start timer

  !---------------------------------------------------------------  
  allocate(r_plan(nc1(1),nc1(3))) 

  !--Tout le monde envoie vers l'est sa derniere colonne et recoit de l'ouest 
  !  les valeurs que l'on doit ajuter  la premiere colonne

  !-- E: dn(:nc1(1),nc1(2),:nc1(3)) ==> W: dn(:nc1(1),1,:nc1(3))+dn(:nc1(1),nc1(2),:nc1(3)) (no contiguous)
  call MPI_SENDRECV(dn(1,nc1(2),1),1,planXZ_type_D,infompi%voisin(E),etiquette, &
       &            r_plan,nc1(1)*nc1(3),QP_MPI_DP,infompi%voisin(W),etiquette, &
       &            infompi%comm,status,ioerr)

  !--on ajoute les moments du processeur voisn E  la premeire colonnne de moment
  dn(:nc1(1),1,:nc1(3)) = dn(:nc1(1),1,:nc1(3)) + r_plan

  !---------------------------------------------------------------
  !--tout le monde envoie vers l'ouest la colonne que l'on vient de modifier
  !-- W: dn(:nc1(1),1,:nc1(3)) ==> E: dn(:nc1(1),nc1(2),:nc1(3)) (no contiguous)
  call MPI_SENDRECV(dn(1,1,1),     1,planXZ_type_D,infompi%voisin(W),etiquette, &
       &            dn(1,nc1(2),1),1,planXZ_type_D,infompi%voisin(E),etiquette, &
       &            infompi%comm,status,ioerr)

  if(nc1(2)/=nc1(3)) then
   deallocate(r_plan)
   allocate(r_plan(nc1(1),nc1(2))) 
  endif

  !---------------------------------------------------------------
  !--Tout le monde envoie vers le nord sa derniere lign et recoit du sud!
  !  les valeurs que l'on doit ajouter  la premiere ligne
  
  !-- N: dn(:nc1(1),:nc1(2),nc1(3)) ==> S: dn(:nc1(1),:nc1(2),nc1(3)) + dn(:nc1(1),:nc1(2),1) (no contiguous)
  call MPI_SENDRECV(dn(1,1,nc1(3)),1,planXY_type_D,infompi%voisin(N),etiquette, &
       &            r_plan,nc1(1)*nc1(2),QP_MPI_DP,infompi%voisin(S),etiquette, &
       &            infompi%comm,status,ioerr)

  !--On ajoute les les moments du processeur voisin S  la premire ligne de moment
  if (infompi%coords(2).ne.(infompi%dims(2)-1)) then 
  dn(:nc1(1),:nc1(2),1) = dn(:nc1(1),:nc1(2),1) + r_plan
  else
  dn(:nc1(1),:nc1(2),1) = dn(:nc1(1),:nc1(2),2)
  endif
 
  !---------------------------------------------------------------
  !--tout le monde envoie vers le sud la ligne que l'on vient de modifier!
  !-- S: dn(:nc1(1),:nc1(2),1) ==> N: dn(:nc1(1),:nc1(2),nc1(3))  (no contiguous)
  call MPI_SENDRECV(dn(1,1,1),     1,planXY_type_D,infompi%voisin(S),etiquette, &
       &            dn(1,1,nc1(3)),1,planXY_type_D,infompi%voisin(N),etiquette, &
       &            infompi%comm,status,ioerr)
 
  if (infompi%coords(2).eq.0) then dn(:nc1(1),:nc1(2),nc1(3))=dn(:nc1(1),:nc1(2),nc1(3)-1)

  deallocate(r_plan)

  !---------------------------------------------------------------
  !--Adding contribution to the first plain
  dn(1,:,:)    = two*dn(1,:,:)   
  
  !--Conditions aux limites en X 
  dn(nc1(1),:,:) = dn(nc1(1)-1,:,:)

  __GETTIME(63,2)!--stop timer
  __WRT_DEBUG_OUT("mtpd3")
 end subroutine mtpd3
 !******************************* END MTPD3 ***************************************


 !!=============================================================
 !!routine: cond_limit/apdh3d_arr3d
 !!
 !! FUNCTION
 !!   Periodic Condition for Electric Field array
 !!
 !! INPUTS
 !!   infompi=type(mpitype) information about MPI
 !!   ncm(3)= dimension of the Electric Field
 !!   e_conv(3),e1_conv(3)=imposed Electric Field ad x,y,z=1 and =2, respect
 !!
 !! SIDE EFFECT
 !!  Efield(ncm(1),ncm(2),ncm(3))=arr3Dtype where put periodic conditions
 subroutine apdh3d_arr3d(Efield,infompi,ncm,e_conv,e1_conv)

  use defs_mpitype
  use defs_arr3Dtype

  integer,intent(in)  :: ncm(3)
  real(dp),intent(in) :: e_conv(3),e1_conv(3)
  type(mpitype),intent(in) :: infompi
  type(arr3Dtype),intent(inout) :: Efield

  integer :: ioerr
  integer :: status(MPI_STATUS_SIZE)

  __WRT_DEBUG_IN("apdh3d_arr3d")
  __GETTIME(63,1)

  !--Envoie vers le NORD et le SUD ----------------------------------

  !-- S: Efield%x(1:ncm(1),2:ncm(2)-1,2) ==>  N: Efield%x(1:ncm(1),2:ncm(2)-1,ncm(3))
  call MPI_SENDRECV(Efield%x(1,2,2)     ,1,planXY_type_E,infompi%voisin(S),etiquette, &
       &            Efield%x(1,2,ncm(3)),1,planXY_type_E,infompi%voisin(N),etiquette, &
       &            infompi%comm,status,ioerr)
  if (infompi%coords(2).eq.0) then Efield%x(1:ncm(1),2:ncm(2)-1,ncm(3))=Efield%x(1:ncm(1),2:ncm(2)-1,ncm(3)-1)

  !-- N: Efield%x(1:ncm(1),2:ncm(2)-1,ncm(3)-1) ==> Efield%x(1:ncm(1),2:ncm(2)-1,1)
  call MPI_SENDRECV(Efield%x(1,2,ncm(3)-1),1,planXY_type_E,infompi%voisin(N),etiquette, &
       &            Efield%x(1,2,1)       ,1,planXY_type_E,infompi%voisin(S),etiquette, &
       &            infompi%comm,status,ioerr)
  if (infompi%coords(2).eq.(infompi%dims(2)-1)) then Efield%x(1:ncm(1),2:ncm(2)-1,1:)=Efield%x(1:ncm(1),2:ncm(2)-1,2)
  !--Envoie vers l'EST et l'OUEST ----------------------------------

  !-- W: Efield%x(1:ncm(1),2,2:ncm(3)-1) ==> E: Efield%x(1:ncm(1),ncm(2),2:ncm(3)-1) (no contiguous)
  call MPI_SENDRECV(Efield%x(1,2,2)     ,1,planXZ_type_E,infompi%voisin(W),etiquette, &
       &            Efield%x(1,ncm(2),2),1,planXZ_type_E,infompi%voisin(E),etiquette, &
       &            infompi%comm,status,ioerr)

  !-- E: Efield%x(1:ncm(1),ncm(2)-1,2:ncm(3)-1) ==> Z: Efield%x(1:ncm(1),1,2:ncm(3)-1) (no contiguous)
  call MPI_SENDRECV(Efield%x(1,ncm(2)-1,2),1,planXZ_type_E,infompi%voisin(E),etiquette, &
       &            Efield%x(1,1,2),       1,planXZ_type_E,infompi%voisin(W),etiquette, &
       &            infompi%comm,status,ioerr)
  
  
  !--The Corners -------
  !--Send the line Efield%x(:,2,2) from SW ==> Efield%x(:,ncm(2),ncm(3)) on NE 
  call MPI_SENDRECV(Efield%x(1,2,2)          ,ncm(1),QP_MPI_DP,infompi%voisin(SW),etiquette, &
       &            Efield%x(1,ncm(2),ncm(3)),ncm(1),QP_MPI_DP,infompi%voisin(NE),etiquette, &
       &            infompi%comm,status,ioerr)

  !--Send the line Efield%x(:,ncm(2)-1,ncm(3)-1) from NE ==> Efield%x(:,1,1) on SW
  call MPI_SENDRECV(Efield%x(1,ncm(2)-1,ncm(3)-1),ncm(1),QP_MPI_DP,infompi%voisin(NE),etiquette, &
       &            Efield%x(1,1,1)              ,ncm(1),QP_MPI_DP,infompi%voisin(SW),etiquette, &
       &            infompi%comm,status,ioerr) 

  !--Send the line Efield%x(:,2,ncm(3)-1) from SE ==> Efield%x(:,ncm(2),1) on NW
  call MPI_SENDRECV(Efield%x(1,2,ncm(3)-1),ncm(1),QP_MPI_DP,infompi%voisin(NW),etiquette, &
       &            Efield%x(1,ncm(2),1)  ,ncm(1),QP_MPI_DP,infompi%voisin(SE),etiquette, &
       &            infompi%comm,status,ioerr)

  !--Send the line Efield%x(:,ncm(2)-1,2) from NW ==> Efield%x(:,1,ncm(3)) on SE
  call MPI_SENDRECV(Efield%x(1,ncm(2)-1,2),ncm(1),QP_MPI_DP,infompi%voisin(SE),etiquette, &
       &            Efield%x(1,1,ncm(3))  ,ncm(1),QP_MPI_DP,infompi%voisin(NW),etiquette, &
       &            infompi%comm,status,ioerr)

  !--Le plan de sortie est egale au plan prcdent
  Efield%x(ncm(1),:,:) = Efield%x(ncm(1)-1,:,:)
  Efield%x(1,:,:) = e_conv(1)
  Efield%x(2,:,:) = e1_conv(1)

  !--Envoie vers le NORD et le SUD ----------------------------------

  !-- S: Efield%y(1:ncm(1),2:ncm(2)-1,2) ==>  N: Efield%y(1:ncm(1),2:ncm(2)-1,ncm(3))
  call MPI_SENDRECV(Efield%y(1,2,2)     ,1,planXY_type_E,infompi%voisin(S),etiquette, &
       &            Efield%y(1,2,ncm(3)),1,planXY_type_E,infompi%voisin(N),etiquette, &
       &            infompi%comm,status,ioerr)

  if (infompi%coords(2).eq.0) then Efield%y(1:ncm(1),2:ncm(2)-1,ncm(3))=Efield%y(1:ncm(1),2:ncm(2)-1,ncm(3)-1)

  !-- N: Efield%y(1:ncm(1),2:ncm(2)-1,ncm(3)-1) ==> Efield%y(1:ncm(1),2:ncm(2)-1,1)
  call MPI_SENDRECV(Efield%y(1,2,ncm(3)-1),1,planXY_type_E,infompi%voisin(N),etiquette, &
       &            Efield%y(1,2,1)       ,1,planXY_type_E,infompi%voisin(S),etiquette, &
       &            infompi%comm,status,ioerr)
  if (infompi%coords(2).eq.(infompi%dims(2)-1)) then Efield%y(1:ncm(1),2:ncm(2)-1,1:)=Efield%y(1:ncm(1),2:ncm(2)-1,2)
  !--Envoie vers l'EST et l'OUEST ----------------------------------

  !-- W: Efield%y(1:ncm(1),2,2:ncm(3)-1) ==> E: Efield%y(1:ncm(1),ncm(2),2:ncm(3)-1) (no contiguous)
  call MPI_SENDRECV(Efield%y(1,2,2)     ,1,planXZ_type_E,infompi%voisin(W),etiquette, &
       &            Efield%y(1,ncm(2),2),1,planXZ_type_E,infompi%voisin(E),etiquette, &
       &            infompi%comm,status,ioerr)

  !-- E: Efield%y(1:ncm(1),ncm(2)-1,2:ncm(3)-1) ==> Z: Efield%y(1:ncm(1),1,2:ncm(3)-1) (no contiguous)
  call MPI_SENDRECV(Efield%y(1,ncm(2)-1,2),1,planXZ_type_E,infompi%voisin(E),etiquette, &
       &            Efield%y(1,1,2),       1,planXZ_type_E,infompi%voisin(W),etiquette, &
       &            infompi%comm,status,ioerr)
  
  
  !--The Corners -------
  !--Send the line Efield%y(:,2,2) from SW ==> Efield%y(:,ncm(2),ncm(3)) on NE 
  call MPI_SENDRECV(Efield%y(1,2,2)          ,ncm(1),QP_MPI_DP,infompi%voisin(SW),etiquette, &
       &            Efield%y(1,ncm(2),ncm(3)),ncm(1),QP_MPI_DP,infompi%voisin(NE),etiquette, &
       &            infompi%comm,status,ioerr)

  !--Send the line Efield%y(:,ncm(2)-1,ncm(3)-1) from NE ==> Efield%y(:,1,1) on SW
  call MPI_SENDRECV(Efield%y(1,ncm(2)-1,ncm(3)-1),ncm(1),QP_MPI_DP,infompi%voisin(NE),etiquette, &
       &            Efield%y(1,1,1)              ,ncm(1),QP_MPI_DP,infompi%voisin(SW),etiquette, &
       &            infompi%comm,status,ioerr) 

  !--Send the line Efield%y(:,2,ncm(3)-1) from SE ==> Efield%y(:,ncm(2),1) on NW
  call MPI_SENDRECV(Efield%y(1,2,ncm(3)-1),ncm(1),QP_MPI_DP,infompi%voisin(NW),etiquette, &
       &            Efield%y(1,ncm(2),1)  ,ncm(1),QP_MPI_DP,infompi%voisin(SE),etiquette, &
       &            infompi%comm,status,ioerr)                               
                                                                             
  !--Send the line Efield%y(:,ncm(2)-1,2) from NW ==> Efield%y(:,1,ncm(3)) on SE
  call MPI_SENDRECV(Efield%y(1,ncm(2)-1,2),ncm(1),QP_MPI_DP,infompi%voisin(SE),etiquette, &
       &            Efield%y(1,1,ncm(3))  ,ncm(1),QP_MPI_DP,infompi%voisin(NW),etiquette, &
       &            infompi%comm,status,ioerr)

  !--Le plan de sortie est egale au plan prcdent
  Efield%y(ncm(1),:,:) = Efield%y(ncm(1)-1,:,:)
  Efield%y(1,:,:) = e_conv(2)
  Efield%y(2,:,:) = e1_conv(2)

  !--Envoie vers le NORD et le SUD ----------------------------------
  !-- S: Efield%z(1:ncm(1),2:ncm(2)-1,2) ==>  N: Efield%z(1:ncm(1),2:ncm(2)-1,ncm(3))
  call MPI_SENDRECV(Efield%z(1,2,2)     ,1,planXY_type_E,infompi%voisin(S),etiquette, &
       &            Efield%z(1,2,ncm(3)),1,planXY_type_E,infompi%voisin(N),etiquette, &
       &            infompi%comm,status,ioerr)
  if (infompi%coords(2).eq.0) then Efield%z(1:ncm(1),2:ncm(2)-1,ncm(3))=-Efield%z(1:ncm(1),2:ncm(2)-1,ncm(3)-1)

  !-- N: Efield%z(1:ncm(1),2:ncm(2)-1,ncm(3)-1) ==> Efield%z(1:ncm(1),2:ncm(2)-1,1)
  call MPI_SENDRECV(Efield%z(1,2,ncm(3)-1),1,planXY_type_E,infompi%voisin(N),etiquette, &
       &            Efield%z(1,2,1)       ,1,planXY_type_E,infompi%voisin(S),etiquette, &
       &            infompi%comm,status,ioerr)
  if (infompi%coords(2).eq.(infompi%dims(2)-1)) then Efield%z(1:ncm(1),2:ncm(2)-1,1:)=-Efield%z(1:ncm(1),2:ncm(2)-1,2)
  !--Envoie vers l'EST et l'OUEST ----------------------------------
  !-- W: Efield%z(1:ncm(1),2,2:ncm(3)-1) ==> E: Efield%z(1:ncm(1),ncm(2),2:ncm(3)-1) (no contiguous)
  call MPI_SENDRECV(Efield%z(1,2,2)     ,1,planXZ_type_E,infompi%voisin(W),etiquette, &
       &            Efield%z(1,ncm(2),2),1,planXZ_type_E,infompi%voisin(E),etiquette, &
       &            infompi%comm,status,ioerr)

  !-- E: Efield%z(1:ncm(1),ncm(2)-1,2:ncm(3)-1) ==> Z: Efield%z(1:ncm(1),1,2:ncm(3)-1) (no contiguous)
  call MPI_SENDRECV(Efield%z(1,ncm(2)-1,2),1,planXZ_type_E,infompi%voisin(E),etiquette, &
       &            Efield%z(1,1,2),       1,planXZ_type_E,infompi%voisin(W),etiquette, &
       &            infompi%comm,status,ioerr)
  
  
  !--The Corners -------
  !--Send the line Efield%z(:,2,2) from SW ==> Efield%z(:,ncm(2),ncm(3)) on NE 
 if (infompi%coords(2).eq.0) then
  call MPI_SENDRECV(Efield%z(1,2,2)          ,ncm(1),QP_MPI_DP,infompi%voisin(SW),etiquette, &
       &            Efield%z(1,ncm(2),ncm(3)),ncm(1),QP_MPI_DP,infompi%voisin(NE),etiquette, &
       &            infompi%comm,status,ioerr)
 else
  call MPI_SENDRECV(Efield%z(1,2,ncm(3)-1)   ,ncm(1),QP_MPI_DP,infompi%voisin(W),etiquette, &
       &            Efield%z(1,ncm(2),ncm(3)),ncm(1),QP_MPI_DP,infompi%voisin(E),etiquette, &
       &            infompi%comm,status,ioerr)
 endif

  !--Send the line Efield%z(:,ncm(2)-1,ncm(3)-1) from NE ==> Efield%z(:,1,1) on SW
  call MPI_SENDRECV(Efield%z(1,ncm(2)-1,ncm(3)-1),ncm(1),QP_MPI_DP,infompi%voisin(NE),etiquette, &
       &            Efield%z(1,1,1)              ,ncm(1),QP_MPI_DP,infompi%voisin(SW),etiquette, &
       &            infompi%comm,status,ioerr) 

  !--Send the line Efield%z(:,2,ncm(3)-1) from SE ==> Efield%z(:,ncm(2),1) on NW
  call MPI_SENDRECV(Efield%z(1,2,ncm(3)-1),ncm(1),QP_MPI_DP,infompi%voisin(NW),etiquette, &
       &            Efield%z(1,ncm(2),1)  ,ncm(1),QP_MPI_DP,infompi%voisin(SE),etiquette, &
       &            infompi%comm,status,ioerr)                               
                                                                             
  !--Send the line Efield%z(:,ncm(2)-1,2) from NW ==> Efield%z(:,1,ncm(3)) on SE
 if (infompi%coords(2).eq.0) then
  call MPI_SENDRECV(Efield%z(1,ncm(2)-1,2),ncm(1),QP_MPI_DP,infompi%voisin(SE),etiquette, &
       &            Efield%z(1,1,ncm(3))  ,ncm(1),QP_MPI_DP,infompi%voisin(NW),etiquette, &
       &            infompi%comm,status,ioerr)
 else
  call MPI_SENDRECV(Efield%z(1,ncm(2)-1,ncm(3)-1),ncm(1),QP_MPI_DP,infompi%voisin(SE),etiquette, &
       &            Efield%z(1,1,ncm(3))  ,ncm(1),QP_MPI_DP,infompi%voisin(NW),etiquette, &
       &            infompi%comm,status,ioerr)
 endif
  !--Le plan de sortie est egale au plan prcdent
  Efield%z(ncm(1),:,:) = Efield%z(ncm(1)-1,:,:)
  Efield%z(1,:,:) = e_conv(3)
  Efield%z(2,:,:) = e1_conv(3)

  __GETTIME(63,2)
  __WRT_DEBUG_OUT("apdh3d_arr3d")
 end subroutine apdh3d_arr3d
 !****************************** END APDH3D_ARR3D ***************************************

 !****************************** B_BOUNDARY ***************************************
 subroutine b_boundary(Afield,ncx1,&
      &                by0,bz0,ey_dm_conv,ez_dm_conv,&
      &                ey_dm_conv1,ez_dm_conv1,vxmean)

  use defs_arr3Dtype

  type(arr3Dtype),intent(inout) :: Afield
  integer,intent(in) :: ncx1
  real(dp),intent(in) :: by0,bz0,vxmean
  real(dp),intent(inout) :: ey_dm_conv,ez_dm_conv,ey_dm_conv1,ez_dm_conv1

  __WRT_DEBUG_IN("b_boundary")
  __GETTIME(63,1)

  !--Open Boundary condition in X-dirextion
  ey_dm_conv  =  vxmean*bz0
  ez_dm_conv  = -vxmean*by0
  ey_dm_conv1 =  vxmean*bz0
  ez_dm_conv1 = -vxmean*by0

  !ax( 1,:,:) = bx0
  !ay( 1,:,:) = by0
  !az( 1,:,:) = bz0
  Afield%x(ncx1,:,:) = Afield%x(ncx1-1,:,:)
  Afield%y(ncx1,:,:) = Afield%y(ncx1-1,:,:)
  Afield%z(ncx1,:,:) = Afield%z(ncx1-1,:,:)

  __GETTIME(63,2)
  __WRT_DEBUG_OUT("b_boundary")
 end subroutine b_boundary

 !*************************** END B_BOUNDARY ************************************


 !!=============================================================
 !!routine: cond_limit/mtpd_four
 !!
 !! FUNCTION
 !!   Periodic condition for densities AND velocity arrays
 !!
 !! INPUTS
 !!   infompi=type(mpitype) information about MPI
 !!   nc1(3)= dimension of the physical array (meanful)
 !!
 !! SIDE EFFECT
 !!  dn(nx+1,ny+1,nz+1)=array density where put periodic conditions
 !!  vel=arr(3d)= velocity  where put periodic conditions
 subroutine mtpd_four(dn,vel,infompi,nc1)

  use defs_mpitype
  use defs_arr3Dtype

  integer,intent(in) :: nc1(3)
  type(mpitype),intent(in) :: infompi
  real(dp),intent(inout)::dn(:,:,:)
  type(arr3Dtype),intent(inout) :: vel

  integer :: ioerr
  integer :: status(MPI_STATUS_SIZE)
!  real(dp),allocatable :: r_line(:)
  real(dp),allocatable :: r_plan(:,:)

  __WRT_DEBUG_IN("mtpd_four")
  __GETTIME(63,1)!--start timer

  !---------------------------------------------------------------  
  allocate(r_plan(nc1(1),nc1(3))) 

  !--Tout le monde envoie vers l'est sa derniere colonne et recoit de l'ouest 
  !  les valeurs que l'on doit ajuter  la premiere colonne

  !--DN-----------------------------------------------------------
  !-- E: dn(:nc1(1),nc1(2),:nc1(3)) ==> W: dn(:nc1(1),1,:nc1(3))+dn(:nc1(1),nc1(2),:nc1(3)) (no contiguous)
  call MPI_SENDRECV(dn(1,nc1(2),1),1,planXZ_type_D,infompi%voisin(E),etiquette, &
       &            r_plan,nc1(1)*nc1(3),QP_MPI_DP,infompi%voisin(W),etiquette, &
       &            infompi%comm,status,ioerr)

  !--on ajoute les moments du processeur voisn E  la premeire colonnne de moment
  dn(:nc1(1),1,:nc1(3)) = dn(:nc1(1),1,:nc1(3)) + r_plan

  !--tout le monde envoie vers l'ouest la colonne que l'on vient de modifier
  !-- W: dn(:nc1(1),1,:nc1(3)) ==> E: dn(:nc1(1),nc1(2),:nc1(3)) (no contiguous)
  call MPI_SENDRECV(dn(1,1,1),     1,planXZ_type_D,infompi%voisin(W),etiquette, &
       &            dn(1,nc1(2),1),1,planXZ_type_D,infompi%voisin(E),etiquette, &
       &            infompi%comm,status,ioerr)

  !--VEL%X-----------------------------------------------------------
  call MPI_SENDRECV(vel%x(1,nc1(2),1),1,planXZ_type_D,infompi%voisin(E),etiquette, &
       &            r_plan,nc1(1)*nc1(3),QP_MPI_DP,   infompi%voisin(W),etiquette, &
       &            infompi%comm,status,ioerr)

  !--on ajoute les moments du processeur voisn E  la premeire colonnne de moment
  vel%x(:nc1(1),1,:nc1(3)) = vel%x(:nc1(1),1,:nc1(3)) + r_plan

  !--tout le monde envoie vers l'ouest la colonne que l'on vient de modifier
  !-- W: dn(:nc1(1),1,:nc1(3)) ==> E: dn(:nc1(1),nc1(2),:nc1(3)) (no contiguous)
  call MPI_SENDRECV(vel%x(1,1,1),     1,planXZ_type_D,infompi%voisin(W),etiquette, &
       &            vel%x(1,nc1(2),1),1,planXZ_type_D,infompi%voisin(E),etiquette, &
       &            infompi%comm,status,ioerr)

  !--VEL%Y-----------------------------------------------------------
  call MPI_SENDRECV(vel%y(1,nc1(2),1),1,planXZ_type_D,infompi%voisin(E),etiquette, &
       &            r_plan,nc1(1)*nc1(3),QP_MPI_DP,   infompi%voisin(W),etiquette, &
       &            infompi%comm,status,ioerr)

  !--on ajoute les moments du processeur voisn E  la premeire colonnne de moment
  vel%y(:nc1(1),1,:nc1(3)) = vel%y(:nc1(1),1,:nc1(3)) + r_plan

  !--tout le monde envoie vers l'ouest la colonne que l'on vient de modifier
  !-- W: dn(:nc1(1),1,:nc1(3)) ==> E: dn(:nc1(1),nc1(2),:nc1(3)) (no contiguous)
  call MPI_SENDRECV(vel%y(1,1,1),1,planXZ_type_D,infompi%voisin(W),etiquette, &
       &            vel%y(1,nc1(2),1),1,planXZ_type_D,infompi%voisin(E),etiquette, &
       &            infompi%comm,status,ioerr)


  !--VEL%Z-----------------------------------------------------------
  call MPI_SENDRECV(vel%z(1,nc1(2),1),1,planXZ_type_D,infompi%voisin(E),etiquette, &
       &            r_plan,nc1(1)*nc1(3),QP_MPI_DP,   infompi%voisin(W),etiquette, &
       &            infompi%comm,status,ioerr)

  !--on ajoute les moments du processeur voisn E  la premeire colonnne de moment
  vel%z(:nc1(1),1,:nc1(3)) = vel%z(:nc1(1),1,:nc1(3)) + r_plan
  
  !--tout le monde envoie vers l'ouest la colonne que l'on vient de modifier
  !-- W: dn(:nc1(1),1,:nc1(3)) ==> E: dn(:nc1(1),nc1(2),:nc1(3)) (no contiguous)
  call MPI_SENDRECV(vel%z(1,1,1)     ,1,planXZ_type_D,infompi%voisin(W),etiquette, &
       &            vel%z(1,nc1(2),1),1,planXZ_type_D,infompi%voisin(E),etiquette, &
       &            infompi%comm,status,ioerr)


  if(nc1(2)/=nc1(3)) then
   deallocate(r_plan)
   allocate(r_plan(nc1(1),nc1(2))) 
  endif

  !--DN-----------------------------------------------------------
  !--Tout le monde envoie vers le nord sa derniere lign et recoit du sud!
  !  les valeurs que l'on doit ajouter  la premiere ligne
  
  !-- N: dn(:nc1(1),:nc1(2),nc1(3)) ==> S: dn(:nc1(1),:nc1(2),nc1(3)) + dn(:nc1(1),:nc1(2),1) (no contiguous)
  call MPI_SENDRECV(dn(1,1,nc1(3)),1,planXY_type_D,infompi%voisin(N),etiquette, &
       &            r_plan,nc1(1)*nc1(2),QP_MPI_DP,infompi%voisin(S),etiquette, &
       &            infompi%comm,status,ioerr)

  !--On ajoute les les moments du processeur voisin S  la premire ligne de moment
  dn(:nc1(1),:nc1(2),1) = dn(:nc1(1),:nc1(2),1) + r_plan

  !--tout le monde envoie vers le sud la ligne que l'on vient de modifier!
  !-- S: dn(:nc1(1),:nc1(2),1) ==> N: dn(:nc1(1),:nc1(2),nc1(3))  (no contiguous)
  call MPI_SENDRECV(dn(1,1,1),1,planXY_type_D,infompi%voisin(S),etiquette, &
       &            dn(1,1,nc1(3)),1,planXY_type_D,infompi%voisin(N),etiquette, &
       &            infompi%comm,status,ioerr)

  !--VEL%X-----------------------------------------------------------
  !--Tout le monde envoie vers le nord sa derniere lign et recoit du sud!
  !  les valeurs que l'on doit ajouter  la premiere ligne
  
  !-- N: vel%x(:nc1(1),:nc1(2),nc1(3)) ==> S: vel%x(:nc1(1),:nc1(2),nc1(3)) + vel%x(:nc1(1),:nc1(2),1) (no contiguous)
  call MPI_SENDRECV(vel%x(1,1,nc1(3)),1,planXY_type_D,infompi%voisin(N),etiquette, &
       &            r_plan,nc1(1)*nc1(2),QP_MPI_DP,infompi%voisin(S),etiquette, &
       &            infompi%comm,status,ioerr)

  !--On ajoute les les moments du processeur voisin S  la premire ligne de moment
  vel%x(:nc1(1),:nc1(2),1) = vel%x(:nc1(1),:nc1(2),1) + r_plan


  !--tout le monde envoie vers le sud la ligne que l'on vient de modifier!
  !-- S: vel%x(:nc1(1),:nc1(2),1) ==> N: vel%x(:nc1(1),:nc1(2),nc1(3))  (no contiguous)
  call MPI_SENDRECV(vel%x(1,1,1),1,planXY_type_D,infompi%voisin(S),etiquette, &
       &            vel%x(1,1,nc1(3)),1,planXY_type_D,infompi%voisin(N),etiquette, &
       &            infompi%comm,status,ioerr)

  !--VEL%Y-----------------------------------------------------------
  !--Tout le monde envoie vers le nord sa derniere lign et recoit du sud!
  !  les valeurs que l'on doit ajouter  la premiere ligne
  
  !-- N: vel%y(:nc1(1),:nc1(2),nc1(3)) ==> S: vel%y(:nc1(1),:nc1(2),nc1(3)) + vel%y(:nc1(1),:nc1(2),1) (no contiguous)
  call MPI_SENDRECV(vel%y(1,1,nc1(3)),1,planXY_type_D,infompi%voisin(N),etiquette, &
       &            r_plan,nc1(1)*nc1(2),QP_MPI_DP,infompi%voisin(S),etiquette, &
       &            infompi%comm,status,ioerr)

  !--On ajoute les les moments du processeur voisin S  la premire ligne de moment
  vel%y(:nc1(1),:nc1(2),1) = vel%y(:nc1(1),:nc1(2),1) + r_plan

  !--tout le monde envoie vers le sud la ligne que l'on vient de modifier!
  !-- S: vel%y(:nc1(1),:nc1(2),1) ==> N: vel%y(:nc1(1),:nc1(2),nc1(3))  (no contiguous)
  call MPI_SENDRECV(vel%y(1,1,1),1,planXY_type_D,infompi%voisin(S),etiquette, &
       &            vel%y(1,1,nc1(3)),1,planXY_type_D,infompi%voisin(N),etiquette, &
       &            infompi%comm,status,ioerr)

  !--VEL%Z-----------------------------------------------------------
  !--Tout le monde envoie vers le nord sa derniere lign et recoit du sud!
  !  les valeurs que l'on doit ajouter  la premiere ligne
  
  !-- N: vel%z(:nc1(1),:nc1(2),nc1(3)) ==> S: vel%z(:nc1(1),:nc1(2),nc1(3)) + vel%z(:nc1(1),:nc1(2),1) (no contiguous)
  call MPI_SENDRECV(vel%z(1,1,nc1(3)),1,planXY_type_D,infompi%voisin(N),etiquette, &
       &            r_plan,nc1(1)*nc1(2),QP_MPI_DP,infompi%voisin(S),etiquette, &
       &            infompi%comm,status,ioerr)

  !--On ajoute les les moments du processeur voisin S  la premire ligne de moment
  vel%z(:nc1(1),:nc1(2),1) = vel%z(:nc1(1),:nc1(2),1) + r_plan

  !--tout le monde envoie vers le sud la ligne que l'on vient de modifier!
  !-- S: vel%z(:nc1(1),:nc1(2),1) ==> N: vel%z(:nc1(1),:nc1(2),nc1(3))  (no contiguous)
  call MPI_SENDRECV(vel%z(1,1,1),1,planXY_type_D,infompi%voisin(S),etiquette, &
       &            vel%z(1,1,nc1(3)),1,planXY_type_D,infompi%voisin(N),etiquette, &
       &            infompi%comm,status,ioerr)


  deallocate(r_plan)

  dn(1,:,:)    = two*dn(1,:,:)   
  vel%x(1,:,:) = two*vel%x(1,:,:)
  vel%y(1,:,:) = two*vel%y(1,:,:)
  vel%z(1,:,:) = two*vel%z(1,:,:)

  !--Conditions aux limites en X 
  dn(nc1(1),:,:)    = dn(nc1(1)-1,:,:)
  vel%x(nc1(1),:,:) = vel%x(nc1(1)-1,:,:)
  vel%y(nc1(1),:,:) = vel%y(nc1(1)-1,:,:)
  vel%z(nc1(1),:,:) = vel%z(nc1(1)-1,:,:)

  __GETTIME(63,2)!--stop timer
  __WRT_DEBUG_OUT("mtpd_four")
 end subroutine mtpd_four
 !******************************* END MTPD_FOUR ***************************************


end module field_cond_limit
