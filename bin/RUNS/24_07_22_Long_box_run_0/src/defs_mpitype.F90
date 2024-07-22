!!=============================================================
!!=============================================================
!!module: defs_mpitype
!! NAME
!!  defs_mpitype (MMancini)
!!
!! FUNCTION
!!  This module contains definitions of types 
!!  for MPI variables and relied functions
module defs_mpitype

#ifndef NOTHAVE_MPI
 use mpi
#endif
 implicit none
 private

 integer,public,save :: rang  = 0 !--Absolute rank of the proc
 integer,public,save :: nproc = 1 !--Overall number of procs

 !--Mapping neighbours of process
 integer,private,parameter :: ndims = 2    !--Dimension of the Cartesian parallelisation
 integer,private,parameter :: nb_voisins = 8 !--Number of Influent neighbours
 integer,private,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8 !--neighbours directions
 !  -----------------
 ! |NW=8 | N=1 | NE=2|
 !  -----------------
 ! | W=7 |     |  E=3|
 !  -----------------
 ! |SW=6 | S=5 | SE=4|
 !  -----------------

 integer,public,save       :: etiquette=100,etiquette1=200,etiquette3=303

 !!=============================================================
 !!--Type for MPI informations
 type,public :: mpitype  
  integer :: me = 0    !--numero du process (me) 
  integer :: nproc = 1 !--nombre de processus (old nb_procs)
  integer :: comm    !--communicator
  integer :: dims(ndims) = 0   !--dimensions of the cartesian grid of parallelisation
  integer :: coord(ndims) = 0  !--Coordinates in the grid
  integer :: voisin(nb_voisins) = 0 !--Neighbours of this proc in the sub-domain
 end type mpitype

 type(mpitype),public,save :: mpiinfo       !--mpi type for normal particles
 type(mpitype),public,save :: mpiinfo_pick  !--mpi type for normal particles

#ifndef NOTHAVE_MPI
 private ::               &
      init_MPI,           &!--Intialisation of MPI (MPI_COMM_WORLD)
      init_mpiinfo,       &!--Intialisation of type mpitype
      voisinage            !--Compute all neighbours for any processus

 public  ::               &
      init_all_mpiinfo,   &!--Intialise two mpitype, normal and pick_up
      print_mpiinfo        !-Print information about parallelisation
contains
 !!#####################################################################

 !!=============================================================
 !!routine: defs_mpitype/init_MPI
 !!
 !! FUNCTION
 !!  Initialize MPI, and MPI_COMM_WORLD
 !! IN 
 !! OUT
 !!   
 !! SIDE EFFECT
 !!  MPI_COMM_WORLD,rang,nproc are initialized
 subroutine init_MPI()
  integer :: ioerr

  !--Initialisation of MPI
  call MPI_INIT(ioerr)

  !--Who are me?
  call MPI_COMM_RANK(MPI_COMM_WORLD,rang,ioerr)

  !--How many procs?
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ioerr)

 end subroutine init_MPI

 !!=============================================================
 !!routine: defs_mpitype/init_mpiinfo
 !!
 !! FUNCTION 
 !!  Initalize cartesian grid MPI
 !!
 !! IN 
 !!  period=if True use periodic grid else no periodic grid
 !! OUT
 !!
 !! SIDE EFFECT
 !!  <type(mpitype)>infompi= Contains cartesian grid MPI information
 !!                           It is initialized only on the MPI part
 subroutine init_mpiinfo(infompi,period)

  logical,intent(in) :: period
  type(mpitype),intent(inout) :: infompi

  integer :: ioerr
  logical :: reorganisation
  logical,dimension(ndims) :: periods  

  !-Initialisation of work variables
  periods(:) = period
  reorganisation = .true.

  !--Write the number of proc in info
  infompi%nproc = nproc

  !--Connaitre le nombre de processus suivant X et Y 
  !  en fonction du nombre total de procesus nproc
  call MPI_DIMS_CREATE(nproc,ndims,infompi%dims,ioerr)

  !--Creation de la grille du processus 2D
  call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,infompi%dims,&
       &               periods,reorganisation,infompi%comm,ioerr)

  !--Connaitre mon rang dans la topologie
  call MPI_COMM_RANK(infompi%comm,infompi%me,ioerr)

  !--Connaitre les coordonées du processus dans la topologie
  call MPI_CART_COORDS(infompi%comm,infompi%me,ndims,infompi%coord,ioerr)

 end subroutine init_mpiinfo

 !!=============================================================
 !!routine: defs_mpitype/voisinage
 !!
 !! FUNCTION
 !!  given a mpitype, computes its coordinates 
 !!  in the cartesian grid associated with its communcator
 !! IN 
 !! OUT
 !!
 !! SIDE EFFECT
 !!  <type(mpiglobaltype)>%coord=coordinates in the grid
 subroutine voisinage(infompi)

  type(mpitype),intent(inout) :: infompi

  integer :: ioerr
  integer,dimension(ndims):: loc_coords

  !--Initialisation du tableau de local coords
  infompi%voisin(:) = MPI_PROC_NULL
  loc_coords(:) = 0

  ! ----- ----- -----
  !|NW=8 | N=1 | NE=2|
  ! ----- ----- -----
  !| W=7 |     | E=3 |
  ! ----- ----- -----
  !|SW=6 | S=5 | SE=4|
  ! ----- ----- -----

  !--Recherhce des voisin Nords, Sud Est et Ouest
  call MPI_CART_SHIFT(infompi%comm,1,1,infompi%voisin(S),infompi%voisin(N),ioerr)
  call MPI_CART_SHIFT(infompi%comm,0,1,infompi%voisin(W),infompi%voisin(E),ioerr)

  !--Recherche des voisins NW 
  if(  (infompi%voisin(N) /= MPI_PROC_NULL) .and. &
       (infompi%voisin(W) /= MPI_PROC_NULL)) then
   loc_coords(:) = infompi%coord(:)+(/-1,1/)
   call MPI_CART_RANK(infompi%comm,loc_coords,infompi%voisin(NW),ioerr)
  endif

  !--Recherches des voisins NE
  if(  (infompi%voisin(N) /= MPI_PROC_NULL) .and. &
       (infompi%voisin(E) /= MPI_PROC_NULL)) then
   loc_coords(:) = infompi%coord(:)+(/1,1/)
   call MPI_CART_RANK(infompi%comm,loc_coords,infompi%voisin(NE),ioerr)
  endif

  !--recherche des voisins SW
  if(  (infompi%voisin(S) /= MPI_PROC_NULL) .and. &
       (infompi%voisin(W) /= MPI_PROC_NULL)) then
   loc_coords(:) = infompi%coord(:)+(/-1,-1/)
   call MPI_CART_RANK(infompi%comm,loc_coords,infompi%voisin(SW),ioerr)
  endif

  !--recherche des voisins SE
  if(  (infompi%voisin(S) /= MPI_PROC_NULL) .and. &
       (infompi%voisin(E) /= MPI_PROC_NULL)) then
   loc_coords(:) = infompi%coord(:)+(/1,-1/)
   call MPI_CART_RANK(infompi%comm,loc_coords,infompi%voisin(SE),ioerr)
  endif

 end subroutine voisinage

 !!=============================================================
 !!routine: defs_mpitype/init_all_mpiinfo
 !!
 !! FUNCTION
 !!  Initialize MPI, Initialize mpi information for normal
 !!  end pickups particles. Compute the neighbours of any 
 !!  cell.
 !!         
 !! IN 
 !! OUT
 !!
 !! SIDE EFFECT
 !!  infompi<type(mpitype)> information mpi about
 !!   normal particles
 !!  infompi_pick<type(mpitype)>information mpi about
 !!   pickups particles
 subroutine init_all_mpiinfo(infompi,infompi_pick)

  type(mpitype),intent(inout) :: infompi
  type(mpitype),intent(inout) :: infompi_pick

  call init_MPI()

  !--NB pickups and normal particles dont have 
  !  the same periodicity so two communicator are needed

  !---------------------------------------------------------
  !--Normal particles with periodic boundary condition (Y,Z)

  !--Create MPI grid and get info in mpiinfo
  call init_mpiinfo(infompi,.true.)

  !--Find the height neighbours of normals particles
  call voisinage(infompi)

  !---------------------------------------------------------
  !--Pickup particles with free boundary condition

  !--Create MPI grid and get info in mpiinfo_pick
  call init_mpiinfo(infompi_pick,.false.)

  !--Find the height neighbours of pickups
  call voisinage(infompi_pick)

 end subroutine init_all_mpiinfo

 !!=============================================================
 !!routine: defs_mpitype/print_mpiinfo
 !!
 !! FUNCTION
 !! Print MPI informations on standard outputs and unit file
 !!         
 !! IN 
 !! OUT
 subroutine print_mpiinfo(unit)

  integer,intent(in) :: unit
  integer :: name_length,ioerr
  integer :: vers,subvers
  character(len=30)  :: name_process
  character(len=200) :: msg

  call MPI_GET_VERSION(vers,subvers,ioerr)
  call MPI_GET_PROCESSOR_NAME(name_process,name_length,ioerr)
  if(rang==0)then
   write(msg,'(2a,2a,i1,a,i2.2,5a,i4,a,2(a,i2),a)')char(10),&
        &" ___________________ MPI informations _______________________",&
        & char(10),"   MPI version: ",vers,".",subvers,&
        & char(10),"   Machine: ",name_process,&
        & char(10),"   MPI processus  = ",nproc,&
        & char(10),"   MPI (Y,Z)-dims = (",&
        & mpiinfo%dims(1),',',mpiinfo%dims(2),')'
   write(6,'(a)')msg
   write(unit,'(a)')msg
  end if
 end subroutine print_mpiinfo

#endif
end module defs_mpitype
