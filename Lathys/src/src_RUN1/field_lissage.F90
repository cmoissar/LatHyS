!!=============================================================
!!=============================================================
module field_lissage

 use mpi
 use defs_basis
 use m_writeout,only     : wrt_debug
 use m_timing,only       : time_get
#include "q-p_common.h"


 implicit none
 private

 integer,save :: size_smooth_B=5!--Size of the smoothing done on the tail of B
 integer,save :: planXZ_type,planXZ_type_B,planXY_type_B !--MPI planes
 real(dp),parameter :: b1=0.25_dp, b2=0.5_dp, b3=0.25_dp!--Parameters for smoothing

 public::                     &
      smth_init_mpi_planes,   &!--Intialize mpi_type_vector
      smth_func,              &!--Compute smoothing using MPI
      smth_free_mpi_planes     !--Free mpi_type_vector

 interface smth_func
  module procedure smth3r_mpi !,smth_four
  module procedure smth_B !,smth_four
 end interface

contains
 !!#####################################################################

 !!=============================================================
 !!routine: lissage/smth_init_mpi_planes
 !!
 !! FUNCTION
 !!  Initialize planes for MPI sendrecv for smoothing of fields
 !!         
 !! IN
 !!  The bounds of the arrays
 !! OUT
 subroutine smth_init_mpi_planes(ncm)

  integer,intent(in) :: ncm(3)
  integer :: ioerr

  size_smooth_B=ncm(1)-1
  call MPI_TYPE_VECTOR(ncm(3)-1+1,ncm(1)-1+1,ncm(1)*ncm(2), QP_MPI_DP,planXZ_type,ioerr)
  call MPI_TYPE_COMMIT(planXZ_type,ioerr)
  call MPI_TYPE_VECTOR(ncm(3)-1+1,size_smooth_B+1,ncm(1)*ncm(2), QP_MPI_DP,planXZ_type_B,ioerr)
  call MPI_TYPE_COMMIT(planXZ_type_B,ioerr)
  call MPI_TYPE_VECTOR(ncm(2)-1+1,size_smooth_B+1,ncm(1), QP_MPI_DP,planXY_type_B,ioerr)
  call MPI_TYPE_COMMIT(planXY_type_B,ioerr)

 end subroutine smth_init_mpi_planes

 !!=============================================================
 !!routine: lissage/smth_free_mpi_planes
 !!
 !! FUNCTION
 !!  Free the types for MPI
 !!         
 !! IN
 !! OUT
 subroutine smth_free_mpi_planes()

  integer :: ioerr
  call MPI_TYPE_FREE(planXZ_type,ioerr)
  call MPI_TYPE_FREE(planXZ_type_B,ioerr)
  call MPI_TYPE_FREE(planXY_type_B,ioerr)
 end subroutine smth_free_mpi_planes

 !*********************************** SMTH3R_MPI *************************************
 subroutine smth3r_mpi(as,nc,infompi)

  use defs_mpitype

  integer,intent(in) :: nc(3)!--Cells numbers
  type(mpitype),intent(in) :: infompi
  real(dp),dimension(:,:,:),intent(inout)  :: as

  integer :: i,j,k,ioerr
  integer :: status(MPI_STATUS_SIZE)
  integer :: ncm(3),nc1(3)
  real(dp) :: as1,as2

  __GETTIME(64,1)
  __WRT_DEBUG_IN("smth3r_mpi")

  !--Physical size arrays
  nc1 = nc+1

  !--Memory size arrays 
  ncm = nc+2

  !--On suppose que les conditions aux limites sont deja appliquees
  !-x-smoothing
  do k = 1,nc1(3)
   do j = 1,nc1(2)
    !--Applications des conditions priodiques en X 
    as(nc1(1),j,k) = as(nc(1),j,k)!!!!!!!!!!!!!!!!!!!!!!!!!!a quoi serve?
    as1 = as(1,j,k)
    do i = 2,nc(1)
     as2 = as(i,j,k)
     as(i,j,k) = b1*as1 + b2*as2 + b3*as(i+1,j,k)
     as1 = as2
    enddo
    !--conditions aux limites
    as(nc1(1),j,k) = as(nc(1),j,k)
   enddo
  enddo

  !--Pour lisser les moments en Y et Z nous avons besoins de connaitre des lignes et 
  !des colonnes des processus voisins.

  !--Envoi des colonnes vers les voisins E
 
  !--Reception par le voisin W de la colonne 
  !-- E: as(1:ncm(1),nc1(2)-1,1:ncm(3)) ==> W: plan_r (no contiguous)
  !--We use the void plan as(:,ncm(2),:)  to receive
  call MPI_SENDRECV(as(1,nc(2),1),1,planXZ_type,infompi%voisin(E),etiquette3, &
       &            as(1,ncm(2),1),1,planXZ_type,infompi%voisin(W),etiquette3, &
       &            infompi%comm,status,ioerr)
 
  !--Envoi des lignes vers le N Reception de la lign par le voisin S
  !-- N: as(1:ncm(1),1:ncm(2),nc(3)) ==> S: plan_r
  !--We use the void plan as(:,ncm(3),:)  to receive
  call MPI_SENDRECV(as(1,1,nc(3)),ncm(1)*ncm(2),QP_MPI_DP,infompi%voisin(N),etiquette3, &
       &            as(1,1,ncm(3)),ncm(1)*ncm(2),QP_MPI_DP,infompi%voisin(S),etiquette3,& 
       &            infompi%comm,status,ioerr)

 
  !--y-smoothing
  as = b1*cshift(as,shift=-1,dim=2)+b2*as+b3*cshift(as,shift=1,dim=2)
  !--z-smoothing
  as = b1*cshift(as,shift=-1,dim=3)+b2*as+b3*cshift(as,shift=1,dim=3)

  !--On applique les conditions periodiques en Y
  !-- W: as(1:ncm(1),1,1:ncm(3)) ==> E: as(1:ncm(1),nc1(2),1:ncm(1))
  call MPI_SENDRECV(as(1,1,1)     ,1,planXZ_type,infompi%voisin(W),etiquette3, &
       &            as(1,nc1(2),1),1,planXZ_type,infompi%voisin(E),etiquette3,  &
       &            infompi%comm,status,ioerr)

  !--On applique les conditions periodiques en Z
  !--S: as(:ncm(1),:ncm(2),1) ==> as(:ncm(1),:ncm(2),nc1(3))
  call MPI_SENDRECV(as(1,1,1),     ncm(1)*ncm(2),QP_MPI_DP,infompi%voisin(S),etiquette3, &
       &            as(1,1,nc1(3)),ncm(1)*ncm(2),QP_MPI_DP,infompi%voisin(N),etiquette3,&
       &            infompi%comm,status,ioerr)


  !--Applications des conditions priodiques en X
  as(nc1(1),:,:) = as(nc(1),:,:) 
 
  __GETTIME(64,2)
  __WRT_DEBUG_OUT("smth3r_mpi")
 end subroutine smth3r_mpi
 !******************************** END SMTH3R_MPI ****************************************



 !*********************************** SMTH_B *************************************
 subroutine smth_B(Bfield,nc,infompi)

  use defs_mpitype
  use defs_arr3Dtype
  use defs_variable,only : Bfield_0,iter
  integer,intent(in) :: nc(3)!--Cells numbers
  type(mpitype),intent(in) :: infompi
  type(arr3Dtype),intent(inout) :: Bfield

  integer :: i,j,k,ioerr,smth0
  integer :: status(MPI_STATUS_SIZE)
  integer :: ncm(3),nc1(3)
  real(dp) :: Bfield1(3),Bfield2(3)

  __GETTIME(64,1)
  __WRT_DEBUG_IN("smth_B")
  Bfield%x=Bfield%x-Bfield_0%x
  Bfield%y=Bfield%y-Bfield_0%y
  Bfield%z=Bfield%z-Bfield_0%z

  size_smooth_B=nc(1)

  !--Physical size arrays
  nc1 = nc+1

  !--Memory size arrays 
  ncm = nc+2

  smth0 = ncm(1)-size_smooth_B

  !--On suppose que les conditions aux limites sont deja appliquees
  !-x-smoothing
  do k = 1,nc1(3)
   do j = 1,nc1(2)
    !--Applications des conditions priodiques en X 
    Bfield%x(nc1(1),j,k) = Bfield%x(nc1(1)-1,j,k)
    Bfield%y(nc1(1),j,k) = Bfield%y(nc1(1)-1,j,k)
    Bfield%z(nc1(1),j,k) = Bfield%z(nc1(1)-1,j,k)
    Bfield1(1) = Bfield%x(smth0,j,k)
    Bfield1(2) = Bfield%y(smth0,j,k)
    Bfield1(3) = Bfield%z(smth0,j,k)
    do i = smth0+1,nc(1)
     Bfield2(1) = Bfield%x(i,j,k)
     Bfield2(2) = Bfield%y(i,j,k)
     Bfield2(3) = Bfield%z(i,j,k)
     Bfield%x(i,j,k) = b1*Bfield1(1) + b2*Bfield2(1) + b3*Bfield%x(i+1,j,k)
     Bfield%y(i,j,k) = b1*Bfield1(2) + b2*Bfield2(2) + b3*Bfield%y(i+1,j,k)
     Bfield%z(i,j,k) = b1*Bfield1(3) + b2*Bfield2(3) + b3*Bfield%z(i+1,j,k)        
     Bfield1 = Bfield2
    enddo
    !--conditions aux limites
    Bfield%x(nc1(1),j,k) = Bfield%x(nc1(1)-1,j,k)
    Bfield%y(nc1(1),j,k) = Bfield%y(nc1(1)-1,j,k)
    Bfield%z(nc1(1),j,k) = Bfield%z(nc1(1)-1,j,k)
   enddo
  enddo

  !--Pour lisser les moments en Y et Z nous avons besoins de connaitre des lignes et 
  !des colonnes des processus voisins.

  !--Envoi des colonnes vers les voisins E
 
  !--Reception par le voisin W de la colonne 
  !-- E: Bfield(1:ncm(1),nc1(2)-1,1:ncm(3)) ==> W: plan_r (no contiguous)
  !--We use the void plan Bfield(:,ncm(2),:)  to receive
  call MPI_SENDRECV(Bfield%x(smth0,nc(2) ,1),1,planXZ_type_B,infompi%voisin(E),etiquette3, &
       &            Bfield%x(smth0,ncm(2),1),1,planXZ_type_B,infompi%voisin(W),etiquette3, &
       &            infompi%comm,status,ioerr)
 
  !--Envoi des lignes vers le N Reception de la lign par le voisin S
  !-- N: Bfield(smth0:ncm(1),:,nc(3)) ==> S: plan_r
  !--We use the void plan Bfield(smth0:ncm(1),:,ncm(3))  to receive
  call MPI_SENDRECV(Bfield%x(smth0,1,nc(3) ),1,planXY_type_B,infompi%voisin(N),etiquette3, &
       &            Bfield%x(smth0,1,ncm(3)),1,planXY_type_B,infompi%voisin(S),etiquette3,& 
       &            infompi%comm,status,ioerr)

  !--Reception par le voisin W de la colonne 
  !-- E: Bfield(1:ncm(1),nc1(2)-1,1:ncm(3)) ==> W: plan_r (no contiguous)
  !--We use the void plan Bfield(:,ncm(2),:)  to receive
  call MPI_SENDRECV(Bfield%y(smth0,nc(2) ,1),1,planXZ_type_B,infompi%voisin(E),etiquette3, &
       &            Bfield%y(smth0,ncm(2),1),1,planXZ_type_B,infompi%voisin(W),etiquette3, &
       &            infompi%comm,status,ioerr)
 
  !--Envoi des lignes vers le N Reception de la lign par le voisin S
  !-- N: Bfield(1:ncm(1),1:ncm(2),nc(3)) ==> S: plan_r
  !--We use the void plan Bfield(:,ncm(3),:)  to receive
  call MPI_SENDRECV(Bfield%y(smth0,1,nc(3) ),1,planXY_type_B,infompi%voisin(N),etiquette3, &
       &            Bfield%y(smth0,1,ncm(3)),1,planXY_type_B,infompi%voisin(S),etiquette3,& 
       &            infompi%comm,status,ioerr)

  !--Reception par le voisin W de la colonne 
  !-- E: Bfield(1:ncm(1),nc1(2)-1,1:ncm(3)) ==> W: plan_r (no contiguous)
  !--We use the void plan Bfield(:,ncm(2),:)  to receive
  call MPI_SENDRECV(Bfield%z(smth0,nc(2) ,1),1,planXZ_type_B,infompi%voisin(E),etiquette3, &
       &            Bfield%z(smth0,ncm(2),1),1,planXZ_type_B,infompi%voisin(W),etiquette3, &
       &            infompi%comm,status,ioerr)
 
  !--Envoi des lignes vers le N Reception de la lign par le voisin S
  !-- N: Bfield(1:ncm(1),1:ncm(2),nc(3)) ==> S: plan_r
  !--We use the void plan Bfield(:,ncm(3),:)  to receive
  call MPI_SENDRECV(Bfield%z(smth0,1,nc(3)),1,planXY_type_B,infompi%voisin(N),etiquette3, &
       &            Bfield%z(smth0,1,ncm(3)),1,planXY_type_B,infompi%voisin(S),etiquette3,& 
       &            infompi%comm,status,ioerr)

  !--y-smoothing
  Bfield%x(smth0:,:,:) = b1*cshift(Bfield%x(smth0:,:,:),shift=-1,dim=2) +&
       &                 b2*Bfield%x(smth0:,:,:) +&
       &                 b3*cshift(Bfield%x(smth0:,:,:),shift=1,dim=2)
  Bfield%y(smth0:,:,:) = b1*cshift(Bfield%y(smth0:,:,:),shift=-1,dim=2) +&
       &                 b2*Bfield%y(smth0:,:,:) +&
       &                 b3*cshift(Bfield%y(smth0:,:,:),shift=1,dim=2)
  Bfield%z(smth0:,:,:) = b1*cshift(Bfield%z(smth0:,:,:),shift=-1,dim=2) +&
       &                 b2*Bfield%z(smth0:,:,:) +&
       &                 b3*cshift(Bfield%z(smth0:,:,:),shift=1,dim=2)

  !--z-smoothing
  Bfield%x(smth0:,:,:) = b1*cshift(Bfield%x(smth0:,:,:),shift=-1,dim=3) +&
       &                 b2*Bfield%x(smth0:,:,:) +&
       &                 b3*cshift(Bfield%x(smth0:,:,:),shift=1,dim=3)
  Bfield%y(smth0:,:,:) = b1*cshift(Bfield%y(smth0:,:,:),shift=-1,dim=3) +&
       &                 b2*Bfield%y(smth0:,:,:) +&
       &                 b3*cshift(Bfield%y(smth0:,:,:),shift=1,dim=3)
  Bfield%z(smth0:,:,:) = b1*cshift(Bfield%z(smth0:,:,:),shift=-1,dim=3) +&
       &                 b2*Bfield%z(smth0:,:,:) +&
       &                 b3*cshift(Bfield%z(smth0:,:,:),shift=1,dim=3)
 
 !--On applique les conditions periodiques en Y
  !-- W: Bfield(1:ncm(1),1,1:ncm(3)) ==> E: Bfield(1:ncm(1),nc1(2),1:ncm(1))
  call MPI_SENDRECV(Bfield%x(smth0,1,1)     ,1,planXZ_type_B,infompi%voisin(W),etiquette3, &
       &            Bfield%x(smth0,nc1(2),1),1,planXZ_type_B,infompi%voisin(E),etiquette3,  &
       &            infompi%comm,status,ioerr)

  !--On applique les conditions periodiques en Z
  !--S: Bfield(:ncm(1),:ncm(2),1) ==> Bfield(:ncm(1),:ncm(2),nc1(3))
  call MPI_SENDRECV(Bfield%x(smth0,1,1),     1,planXY_type_B,infompi%voisin(S),etiquette3, &
       &            Bfield%x(smth0,1,nc1(3)),1,planXY_type_B,infompi%voisin(N),etiquette3,&
       &            infompi%comm,status,ioerr)

  !--On applique les conditions periodiques en Y
  !-- W: Bfield(1:ncm(1),1,1:ncm(3)) ==> E: Bfield(1:ncm(1),nc1(2),1:ncm(1))
  call MPI_SENDRECV(Bfield%y(smth0,1,1)     ,1,planXZ_type_B,infompi%voisin(W),etiquette3, &
       &            Bfield%y(smth0,nc1(2),1),1,planXZ_type_B,infompi%voisin(E),etiquette3,  &
       &            infompi%comm,status,ioerr)

  !--On applique les conditions periodiques en Z
  !--S: Bfield(:ncm(1),:ncm(2),1) ==> Bfield(:ncm(1),:ncm(2),nc1(3))
  call MPI_SENDRECV(Bfield%y(smth0,1,1),     1,planXY_type_B,infompi%voisin(S),etiquette3, &
       &            Bfield%y(smth0,1,nc1(3)),1,planXY_type_B,infompi%voisin(N),etiquette3,&
       &            infompi%comm,status,ioerr)

  !--On applique les conditions periodiques en Y
  !-- W: Bfield(1:ncm(1),1,1:ncm(3)) ==> E: Bfield(1:ncm(1),nc1(2),1:ncm(1))
  call MPI_SENDRECV(Bfield%z(smth0,1,1)     ,1,planXZ_type_B,infompi%voisin(W),etiquette3, &
       &            Bfield%z(smth0,nc1(2),1),1,planXZ_type_B,infompi%voisin(E),etiquette3,  &
       &            infompi%comm,status,ioerr)

  !--On applique les conditions periodiques en Z
  !--S: Bfield(:ncm(1),:ncm(2),1) ==> Bfield(:ncm(1),:ncm(2),nc1(3))
  call MPI_SENDRECV(Bfield%z(smth0,1,1),     1,planXY_type_B,infompi%voisin(S),etiquette3, &
       &            Bfield%z(smth0,1,nc1(3)),1,planXY_type_B,infompi%voisin(N),etiquette3,&
       &            infompi%comm,status,ioerr)

  !--Applications des conditions priodiques en X
  Bfield%x(nc1(1),:,:) = Bfield%x(nc1(1)-1,:,:) 
  Bfield%y(nc1(1),:,:) = Bfield%y(nc1(1)-1,:,:) 
  Bfield%z(nc1(1),:,:) = Bfield%z(nc1(1)-1,:,:) 


  Bfield%x=Bfield%x+Bfield_0%x
  Bfield%y=Bfield%y+Bfield_0%y
  Bfield%z=Bfield%z+Bfield_0%z


  __GETTIME(64,2)
  __WRT_DEBUG_OUT("smth_B")
 end subroutine smth_B
 !******************************** END SMTH_B ****************************************
 
 

end module field_lissage
