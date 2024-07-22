!!=============================================================
!!=============================================================
module particle_com

 use defs_basis
 use defs_mpitype
 use defs_particletype
 use defs_counts_types,only  : count_particle_type
 use defs_variable,only           : s_max
 use m_writeout
 use m_timing,only           : time_get
 use mpi
#include "q-p_common.h"

 implicit none
 private

 private::               &
      pack_part,         &
      communication,     &
      pre_communication, &
      rangement

 public:: pack_com

contains
 !!#####################################################################

 !******************************** PACK_PART ******************************************
 subroutine pack_part(Preg,index_exit,&
      &               p_out_N,p_out_NE,p_out_E,p_out_SE,& 
      &               p_out_S,p_out_SW,p_out_W,p_out_NW,&
      &               pick_out_N,pick_out_NE,pick_out_E,pick_out_SE,&
      &               pick_out_S,pick_out_SW,pick_out_W,pick_out_NW,&
      &               particule)

  integer,intent(in) :: index_exit(:)
  integer,intent(out) :: p_out_N(:)
  integer,intent(out) :: p_out_NE(:)
  integer,intent(out) :: p_out_E(:)
  integer,intent(out) :: p_out_SE(:)
  integer,intent(out) :: p_out_S(:)
  integer,intent(out) :: p_out_SW(:)
  integer,intent(out) :: p_out_W(:)
  integer,intent(out) :: p_out_NW(:)
  integer,intent(out) :: pick_out_N(:)
  integer,intent(out) :: pick_out_NE(:)
  integer,intent(out) :: pick_out_E(:)
  integer,intent(out) :: pick_out_SE(:)
  integer,intent(out) :: pick_out_S(:)
  integer,intent(out) :: pick_out_SW(:)
  integer,intent(out) :: pick_out_W(:)
  integer,intent(out) :: pick_out_NW(:)
  type(count_particle_type),intent(inout) :: Preg
  type(particletype),intent(in) ::particule(:)

  logical :: part_test  
  integer :: jj,ii
  integer :: index_pick_out(8),index_p_out(8)
   
  __WRT_DEBUG_IN("pack_part")
  __GETTIME(35,1)!--Timer start

  !--Intialize------
  index_p_out = 1
  index_pick_out = 1

  !--Loop on the exiting particles
  do ii = 1,Preg%out_tot
   !--Here the index of particles which going out
   jj = index_exit(ii)

   !--If particle goes out on x-direction the cycle   
   if(particule(jj)%dir== -1_iXb) cycle
   
   part_test = (particule(jj)%exc==zero).and.(particule(jj)%orig==0_iXb)
   if(part_test) then
    select case(int(particule(jj)%dir))
    case(N)
     p_out_N(index_p_out(N)) = jj-1
     index_p_out(N) = index_p_out(N)+1
    case(NE) 
     p_out_NE(index_p_out(NE)) = jj-1
     index_p_out(NE) = index_p_out(NE)+1
    case(E) 
     p_out_E(index_p_out(E)) = jj-1
     index_p_out(E) = index_p_out(E)+1
    case(SE)
     p_out_SE(index_p_out(SE)) = jj-1
     index_p_out(SE) = index_p_out(SE)+1
    case(S) 
     p_out_S(index_p_out(S)) = jj-1
     index_p_out(S) = index_p_out(S)+1
    case(SW)
     p_out_SW(index_p_out(SW)) = jj-1
     index_p_out(SW) = index_p_out(SW)+1
    case(W) 
     p_out_W(index_p_out(W)) = jj-1
     index_p_out(W) = index_p_out(W)+1
    case(NW)
     p_out_NW(index_p_out(NW)) = jj-1
     index_p_out(NW) = index_p_out(NW)+1
    end select
   else
    select case(int(particule(jj)%dir))
    case(N)
     pick_out_N(index_pick_out(N)) = jj-1
     index_pick_out(N) = index_pick_out(N)+1
    case(NE) 
     pick_out_NE(index_pick_out(NE)) = jj-1
     index_pick_out(NE) = index_pick_out(NE)+1
    case(E) 
     pick_out_E(index_pick_out(E)) = jj-1
     index_pick_out(E) = index_pick_out(E)+1
    case(SE)
     pick_out_SE(index_pick_out(SE)) = jj-1
     index_pick_out(SE) = index_pick_out(SE)+1
    case(S) 
     pick_out_S(index_pick_out(S)) = jj-1
     index_pick_out(S) = index_pick_out(S)+1
    case(SW)
     pick_out_SW(index_pick_out(SW)) = jj-1
     index_pick_out(SW) = index_pick_out(SW)+1
    case(W) 
     pick_out_W(index_pick_out(W)) = jj-1
     index_pick_out(W) = index_pick_out(W)+1
    case(NW)
     pick_out_NW(index_pick_out(NW)) = jj-1
     index_pick_out(NW) = index_pick_out(NW)+1
    end select
   end if
  enddo

  __GETTIME(35,2)!--Timer stop
  __WRT_DEBUG_OUT("pack_part")
 end subroutine pack_part
 !******************************** END PACK_PART *************************************

 !***************************** PRE_COMMUNICATION **********************************
 subroutine pre_communication(np_out,np_in,infompi, etiquette1)

  integer,intent(in) :: etiquette1
  integer,dimension(nb_voisins),intent(in) :: np_out
  integer,dimension(nb_voisins),intent(out):: np_in
  type(mpitype),intent(in) :: infompi
  integer :: ioerr
  integer :: statut(MPI_STATUS_SIZE)

  __WRT_DEBUG_IN("pre_communication")
  __GETTIME(36,1)!--Timer start

  np_in = 0

  !--on envoie le nombre de particule d'change vers le nord et on le recoit par le sud
  call MPI_SENDRECV(np_out(N),1,MPI_INTEGER,infompi%voisin(N),etiquette1, &
       np_in(S),1,MPI_INTEGER,infompi%voisin(S),etiquette1,infompi%comm,statut,ioerr)

  !--on envoie le nombre de particule d'change vers le nord-est et on le recoit par le sud-ouest
  call MPI_SENDRECV(np_out(NE),1,MPI_INTEGER,infompi%voisin(NE),etiquette1, &
       np_in(SW),1,MPI_INTEGER,infompi%voisin(SW),etiquette1,infompi%comm,statut,ioerr)

  !--on envoie le nombre de particule d'change vers l'est et on le recoit par l'ouest
  call MPI_SENDRECV(np_out(E),1,MPI_INTEGER,infompi%voisin(E),etiquette1, &
       np_in(W),1,MPI_INTEGER,infompi%voisin(W),etiquette1,infompi%comm,statut,ioerr)

  !--on envoie le nombre de particule d'change vers le sud-est et on le recoit par le nord-ouest
  call MPI_SENDRECV(np_out(SE),1,MPI_INTEGER,infompi%voisin(SE),etiquette1, &
       np_in(NW),1,MPI_INTEGER,infompi%voisin(NW),etiquette1,infompi%comm,statut,ioerr)

  !--on envoie le nombre de particule d'change vers le sud et on le recoit par le nord
  call MPI_SENDRECV(np_out(S),1,MPI_INTEGER,infompi%voisin(S),etiquette1, &
       np_in(N),1,MPI_INTEGER,infompi%voisin(N),etiquette1,infompi%comm,statut,ioerr)

  !--on envoie le nombre de particule d'change vers le sud-ouest et on le recoit par le nord-est
  call MPI_SENDRECV(np_out(SW),1,MPI_INTEGER,infompi%voisin(SW),etiquette1, &
       np_in(NE),1,MPI_INTEGER,infompi%voisin(NE),etiquette1,infompi%comm,statut,ioerr)

  !--on envoie le nombre de particule d'change vers l'ouest et on le recoit par l'est
  call MPI_SENDRECV(np_out(W),1,MPI_INTEGER,infompi%voisin(W),etiquette1, &
       np_in(E),1,MPI_INTEGER,infompi%voisin(E),etiquette1,infompi%comm,statut,ioerr)

  !--on envoie le nombre de particule d'change vers le nord-ouest et on le recoit par le sud-est
  call MPI_SENDRECV(np_out(NW),1,MPI_INTEGER,infompi%voisin(NW),etiquette1, &
       np_in(SE),1,MPI_INTEGER,infompi%voisin(SE),etiquette1,infompi%comm,statut,ioerr)

  __GETTIME(36,2)!--Timer start
  __WRT_DEBUG_OUT("pre_communication")
 end subroutine pre_communication
 !********************************** END PRE_COMMUNICATION *********************************


 !***************************** COMMUNICATION **********************************
 subroutine communication(&
      &                   particule,&
      &                   p_out_N,p_out_NE,p_out_E,p_out_SE,&
      &                   p_out_S,p_out_SW,p_out_W,p_out_NW, &
      &                   p_in,np_out,np_in,np_trace,tot_in, &
      &                   nb_voisinss,&
      &                   comm2d,voisin)

  integer,intent(in) :: nb_voisinss,comm2d,tot_in
  type(particletype),intent(in) :: particule(:)
  type(particletype),intent(out) :: p_in(tot_in)
  integer,intent(in) :: p_out_N(:)
  integer,intent(in) :: p_out_NE(:)
  integer,intent(in) :: p_out_E(:) 
  integer,intent(in) :: p_out_SE(:)
  integer,intent(in) :: p_out_S(:) 
  integer,intent(in) :: p_out_SW(:)
  integer,intent(in) :: p_out_W(:) 
  integer,intent(in) :: p_out_NW(:)
  integer,dimension(nb_voisinss),intent(in) :: voisin,np_out,np_in,np_trace

  integer :: ioerr,MPI_particle_block
  integer :: statut(MPI_STATUS_SIZE)
  integer,allocatable :: block_len(:)

  __WRT_DEBUG_IN("communication")
  __GETTIME(37,1)!--Timer start

  !--Array to containt blocks length (=1)
  allocate(block_len(maxval(np_out)))
  block_len = 1

  !--On envoie les particules  vers le nord et on recoit par le sud
  if(np_out(N)/=0)then
   call MPI_TYPE_INDEXED(np_out(N), block_len(:np_out(N)), p_out_N, MPI_particle,MPI_particle_block,ioerr)
  else
   call MPI_TYPE_INDEXED(np_out(N), (/0/), (/0/), MPI_particle,MPI_particle_block,ioerr)
  endif
  call MPI_TYPE_COMMIT(MPI_particle_block,ioerr)
  call MPI_SENDRECV(particule,1,MPI_particle_block,voisin(N),etiquette1, &
       &            p_in(np_trace(S))%pos,np_in(S),MPI_particle,voisin(S),etiquette1,&
       &            comm2d,statut,ioerr)
  call MPI_TYPE_FREE(MPI_particle_block,ioerr)

  !--on envoie les particules par le nord-est et on recoit par le sud-ouest
  if(np_out(NE)/=0)then
   call MPI_TYPE_INDEXED(np_out(NE), block_len(:np_out(NE)), p_out_NE, MPI_particle,MPI_particle_block,ioerr)
  else
   call MPI_TYPE_INDEXED(np_out(NE), (/0/), (/0/), MPI_particle,MPI_particle_block,ioerr)
  endif
  call MPI_TYPE_COMMIT(MPI_particle_block,ioerr)   
  call MPI_SENDRECV(particule,1,MPI_particle_block,voisin(NE),etiquette1, &
       &            p_in(np_trace(SW))%pos,np_in(SW),MPI_particle,voisin(SW),etiquette1,&
       &            comm2d,statut,ioerr)
  call MPI_TYPE_FREE(MPI_particle_block,ioerr)

  !--on envoie les particules par l'est et on recoit les particules par l'ouest 
  if(np_out(E)/=0)then
   call MPI_TYPE_INDEXED(np_out(E), block_len(:np_out(E)), p_out_E, MPI_particle,MPI_particle_block,ioerr)
  else
   call MPI_TYPE_INDEXED(np_out(E), (/0/), (/0/), MPI_particle,MPI_particle_block,ioerr)
  endif
  call MPI_TYPE_COMMIT(MPI_particle_block,ioerr)
  call MPI_SENDRECV(particule,1,MPI_particle_block,voisin(E),etiquette1, &
       &            p_in(np_trace(W))%pos,np_in(W),MPI_particle,voisin(W),&
       &            etiquette1,comm2d,statut,ioerr)
  call MPI_TYPE_FREE(MPI_particle_block,ioerr)

  !--on envoie les particules par le sud-est et on recoit les particules par le nord-ouest
  if(np_out(SE)/=0)then
    call MPI_TYPE_INDEXED(np_out(SE), block_len(:np_out(SE)), p_out_SE, MPI_particle,MPI_particle_block,ioerr)
  else
   call MPI_TYPE_INDEXED(np_out(SE), (/0/), (/0/), MPI_particle,MPI_particle_block,ioerr)
  endif
  call MPI_TYPE_COMMIT(MPI_particle_block,ioerr)
  call MPI_SENDRECV(particule,1,MPI_particle_block,voisin(SE),etiquette1, &
       &            p_in(np_trace(NW))%pos,np_in(NW),MPI_particle,voisin(NW),etiquette1,&
       &            comm2d,statut,ioerr)
  call MPI_TYPE_FREE(MPI_particle_block,ioerr)

  !--on envoie les particules par le le sud et on recoit les particules par le nord
  if(np_out(S)/=0)then
   call MPI_TYPE_INDEXED(np_out(S), block_len(:np_out(S)), p_out_S, MPI_particle,MPI_particle_block,ioerr)
  else
   call MPI_TYPE_INDEXED(np_out(S), (/0/), (/0/), MPI_particle,MPI_particle_block,ioerr)
  endif
  call MPI_TYPE_COMMIT(MPI_particle_block,ioerr)
  call MPI_SENDRECV(particule,1,MPI_particle_block,voisin(S),etiquette1, &
       &            p_in(np_trace(N))%pos,np_in(N),MPI_particle,voisin(N),etiquette1,&
       &           comm2d,statut,ioerr)
  call MPI_TYPE_FREE(MPI_particle_block,ioerr)

  !--on envoie les particules par le sud-ouest et on recoit les particules par le nord-est
  if(np_out(SW)/=0)then
   call MPI_TYPE_INDEXED(np_out(SW), block_len(:np_out(SW)), p_out_SW, MPI_particle,MPI_particle_block,ioerr)
  else
   call MPI_TYPE_INDEXED(np_out(SW), (/0/), (/0/), MPI_particle,MPI_particle_block,ioerr)
  endif
  call MPI_TYPE_COMMIT(MPI_particle_block,ioerr)
  call MPI_SENDRECV(particule,1,MPI_particle_block,voisin(SW),etiquette1, &
       &            p_in(np_trace(NE))%pos,np_in(NE),MPI_particle,voisin(NE),etiquette1,&
       &            comm2d,statut,ioerr)
  call MPI_TYPE_FREE(MPI_particle_block,ioerr)

  !--on envoie les particules par l'ouest et on recoit par l'est
  if(np_out(W)/=0)then
   call MPI_TYPE_INDEXED(np_out(W), block_len(:np_out(W)), p_out_W, MPI_particle,MPI_particle_block,ioerr)
  else
   call MPI_TYPE_INDEXED(np_out(W), (/0/), (/0/), MPI_particle,MPI_particle_block,ioerr)
  endif
  call MPI_TYPE_COMMIT(MPI_particle_block,ioerr)
  call MPI_SENDRECV(particule,1,MPI_particle_block,voisin(W),etiquette1, &
       &            p_in(np_trace(E))%pos,np_in(E),MPI_particle,voisin(E),etiquette1,&
       &            comm2d,statut,ioerr)
  call MPI_TYPE_FREE(MPI_particle_block,ioerr)

  !--on envoie les particules par le nord-ouest et on recoit par le sud-est
  if(np_out(NW)/=0)then
   call MPI_TYPE_INDEXED(np_out(NW), block_len(:np_out(NW)), p_out_NW, MPI_particle,MPI_particle_block,ioerr)
  else
   call MPI_TYPE_INDEXED(np_out(NW), (/0/), (/0/), MPI_particle,MPI_particle_block,ioerr)
  endif
  call MPI_TYPE_COMMIT(MPI_particle_block,ioerr)
  call MPI_SENDRECV(particule,1,MPI_particle_block,voisin(NW),etiquette1, &
       &            p_in(np_trace(SE))%pos,np_in(SE),MPI_particle,voisin(SE),etiquette1,&
       &            comm2d,statut,ioerr)
  call MPI_TYPE_FREE(MPI_particle_block,ioerr)

  deallocate(block_len)
  
  __GETTIME(37,2)!--Timer stop
  __WRT_DEBUG_OUT("communication")
 end subroutine communication
 !********************************** END COMMUNICATION *********************************


 !********************************** RANGEMENT *******************************************
 subroutine rangement(Preg,index_exit,p_in,&
    &                 particule,np_out,np_out_pick,np, &
    &                 s_min,s_r,s_min_loc,s_max_loc &

#ifdef HAVE_DEBUG
    &                 ,xmax &
#endif
    &                 )

  
#ifdef HAVE_DEBUG
  use defs_mpitype
#endif

#ifdef HAVE_DEBUG
  real(dp),intent(in) :: xmax
#endif
  integer,intent(inout) :: np
  integer,intent(inout) :: index_exit(:)
  integer,dimension(nb_voisins),intent(in) :: np_out,np_out_pick
  real(dp),intent(in) :: s_min(3),s_r(3),s_min_loc(3),s_max_loc(3)
  type(count_particle_type),intent(inout) :: Preg
  type(particletype),intent(in) :: p_in(:)
  type(particletype),intent(inout) :: particule(:)

  integer :: j,nlast,incr = 0
  integer :: compt,indice,np_ini
  character(len=200) :: msg

  __WRT_DEBUG_IN("rangement")
  __GETTIME(38,1)!--Timer start

  np_ini = np
  compt = 1

  !--premiere etape
  ! on rempli tous les "trous" des particules qui sont parties vers un autre processeur
  ! par les particules qui viennent d'un autre processeur
  ! et s'il est ncessaire (dans le cas ou nombre entree > nombre sortie), on rajoute 
  ! les particules a la suite

#ifdef HAVE_DEBUG 
  if(wrtscreen == 0) then 
   if(mpiinfo%me == 0) then
    write(*,*) 
    write(*,*) "-----------------------------"
    write(*,*) "         RANGEMENT"
    write(*,*) " processus   sortie X&Y pickup sortie et non rempl.   sortie Z   sortie   npart"
   endif

   print '(i10,i13,i22,i11,i10,i10)', &
        mpiinfo%me,count((particule(:)%dir/=0_iXb).and.(particule(:)%pos(1) <= xmax)), &
        abs(Preg%pp%in-sum(np_out_pick)), &
        count(particule(:)%pos(1)>xmax),count(particule(:)%dir/=0_iXb),np
  endif !--wrtscreen
#endif
  !*********************************** on remplace les particules **********************

  compt = 1

  !--On remplace les particules qui sortent vers les autres processeurs 
  !par ceux qui sont recus par les autres processeurs
  nlast = np

  write(msg,'(a,4i8)')"rangement ",&
       &  Preg%in_tot,Preg%out_xm,&
       &  Preg%out_xp,Preg%out_tot
  call wrt_double(qp_out,msg,wrtscreen,1)

  do j = 1,min(Preg%out_tot,Preg%in_tot)
   indice = index_exit(j)
   compt = compt+1
   particule(indice)%pos(1) = p_in(j)%pos(1)
   particule(indice)%pos(2) = modulo(p_in(j)%pos(2)-s_min(2),s_r(2))
   particule(indice)%pos(3) = modulo(p_in(j)%pos(3)-s_min(3),s_r(3))
   particule(indice)%vel    = p_in(j)%vel
   particule(indice)%mass   = p_in(j)%mass
   particule(indice)%char   = p_in(j)%char
   particule(indice)%exc    = p_in(j)%exc
   particule(indice)%orig   = p_in(j)%orig
   particule(indice)%dir    = 0_iXb
if (any((particule(indice)%pos-s_min_loc).lt.0)) then
        if ((s_max_loc(1).eq.s_r(1)).and.(abs(particule(indice)%pos(1)).lt.0.01)) particule(indice)%pos(1)=s_min_loc(1)+1E-10
        if ((s_max_loc(2).eq.s_r(2)).and.(abs(particule(indice)%pos(2)).lt.0.01)) particule(indice)%pos(2)=s_max_loc(2)-1E-10
        if ((s_max_loc(3).eq.s_r(3)).and.(abs(particule(indice)%pos(3)).lt.0.01)) particule(indice)%pos(3)=s_max_loc(3)-1E-10
endif
if (any((particule(indice)%pos-s_max_loc).gt.0)) then
        if ((s_min_loc(1).eq.0.).and.(abs(particule(indice)%pos(1)-s_r(1)).lt.0.01)) particule(indice)%pos(1)=s_max_loc(1)-1E-10
        if ((s_min_loc(2).eq.0.).and.(abs(particule(indice)%pos(2)-s_r(2)).lt.0.01)) particule(indice)%pos(2)=s_min_loc(2)+1E-10
        if ((s_min_loc(3).eq.0.).and.(abs(particule(indice)%pos(3)-s_r(3)).lt.0.01)) particule(indice)%pos(3)=s_min_loc(3)+1E-10
endif
if (any((particule(indice)%pos-s_min_loc).lt.0)) then
        if ((s_min_loc(1)-particule(indice)%pos(1)).lt.0.001) particule(indice)%pos(1)=s_min_loc(1)+1E-10
        if ((s_min_loc(2)-particule(indice)%pos(2)).lt.0.001) particule(indice)%pos(2)=s_min_loc(2)+1E-10
        if ((s_min_loc(3)-particule(indice)%pos(3)).lt.0.001) particule(indice)%pos(3)=s_min_loc(3)+1E-10
endif
if (any((particule(indice)%pos-s_max_loc).gt.0)) then
        if ((particule(indice)%pos(1)-s_max_loc(1)).lt.0.001) particule(indice)%pos(1)=s_max_loc(1)-1E-10
        if ((particule(indice)%pos(2)-s_max_loc(2)).lt.0.001) particule(indice)%pos(2)=s_max_loc(2)-1E-10
        if ((particule(indice)%pos(3)-s_max_loc(3)).lt.0.001) particule(indice)%pos(3)=s_max_loc(3)-1E-10
endif

if ((any((particule(indice)%pos-s_min_loc).lt.0)).or.(any((particule(indice)%pos-s_max_loc).gt.0))) then

print *,"Rangement problem: particle out of bounds case 1"
print *,'Particle #: ',compt
print *,'Position ',particule(indice)%pos
print *,'Velocity ',particule(indice)%vel
print *,'Origin ',particule(indice)%orig
print *,'Ch. exc ',particule(indice)%exc
print *,'s_min_loc ',s_min_loc
print *,'s_max_loc ',s_max_loc
particule(indice)%pos(1:3) = s_min_loc(1:3)+1E-10
particule(indice)%vel(1:3) = zero
particule(indice)%mass   = 0.0001
particule(indice)%char   = 0.0001
particule(indice)%orig   = -1
endif
  enddo
  
!if (min(Preg%out_tot,Preg%in_tot) == 0) then
! print *, 'proc, Min de sortie - entree = 0 ',mpiinfo%me,Preg%out_tot,Preg%in_tot
!endif

  !--si on fait rentrer plus de particules que l'on en fait sortir
  if(Preg%in_tot > Preg%out_tot) then
   do j = Preg%out_tot+1,Preg%in_tot
    np = np_ini + j-Preg%out_tot
    particule(np)%pos(1) = p_in(j)%pos(1)
    particule(np)%pos(2) = modulo(p_in(j)%pos(2)-s_min(2),s_r(2))
    particule(np)%pos(3) = modulo(p_in(j)%pos(3)-s_min(3),s_r(3))
    particule(np)%vel    = p_in(j)%vel
    particule(np)%mass   = p_in(j)%mass
    particule(np)%char   = p_in(j)%char
    particule(np)%exc    = p_in(j)%exc
    particule(np)%orig   = p_in(j)%orig
    particule(np)%dir    = 0_iXb
if (any((particule(np)%pos-s_min_loc).lt.0)) then
        if ((s_max_loc(1).eq.s_r(1)).and.(abs(particule(np)%pos(1)).lt.0.01)) particule(np)%pos(1)=s_min_loc(1)+1E-10
        if ((s_max_loc(2).eq.s_r(2)).and.(abs(particule(np)%pos(2)).lt.0.01)) particule(np)%pos(2)=s_max_loc(2)-1E-10
        if ((s_max_loc(3).eq.s_r(3)).and.(abs(particule(np)%pos(3)).lt.0.01)) particule(np)%pos(3)=s_max_loc(3)-1E-10
endif
if (any((particule(np)%pos-s_max_loc).gt.0)) then
        if ((s_min_loc(1).eq.0.).and.(abs(particule(np)%pos(1)-s_r(1)).lt.0.01)) particule(np)%pos(1)=s_max_loc(1)-1E-10
        if ((s_min_loc(2).eq.0.).and.(abs(particule(np)%pos(2)-s_r(2)).lt.0.01)) particule(np)%pos(2)=s_min_loc(2)+1E-10
        if ((s_min_loc(3).eq.0.).and.(abs(particule(np)%pos(3)-s_r(3)).lt.0.01)) particule(np)%pos(3)=s_min_loc(3)+1E-10
endif
if (any((particule(np)%pos-s_min_loc).lt.0)) then
        if ((s_min_loc(1)-particule(np)%pos(1)).lt.0.001) particule(np)%pos(1)=s_min_loc(1)+1E-10
        if ((s_min_loc(2)-particule(np)%pos(2)).lt.0.001) particule(np)%pos(2)=s_min_loc(2)+1E-10
        if ((s_min_loc(3)-particule(np)%pos(3)).lt.0.001) particule(np)%pos(3)=s_min_loc(3)+1E-10
endif
if (any((particule(np)%pos-s_max_loc).gt.0)) then
        if ((particule(np)%pos(1)-s_max_loc(1)).lt.0.001) particule(np)%pos(1)=s_max_loc(1)-1E-10
        if ((particule(np)%pos(2)-s_max_loc(2)).lt.0.001) particule(np)%pos(2)=s_max_loc(2)-1E-10
        if ((particule(np)%pos(3)-s_max_loc(3)).lt.0.001) particule(np)%pos(3)=s_max_loc(3)-1E-10
endif
if ((any((particule(np)%pos-s_min_loc).lt.0)).or.(any((particule(np)%pos-s_max_loc).gt.0))) then

print *,"Rangement problem: particle out of bounds case 2"
print *,'Particle #: ',compt
print *,'Position ',particule(np)%pos
print *,'Velocity ',particule(np)%vel
print *,'Origin ',particule(np)%orig
print *,'Ch. exc ',particule(np)%exc
print *,'s_min_loc ',s_min_loc
print *,'s_max_loc ',s_max_loc
particule(np)%pos(1:3) = s_min_loc(1:3)+1E-10
particule(np)%vel(1:3) = zero
particule(np)%mass   = 0.0001
particule(np)%char   = 0.0001
particule(np)%orig   = -1

endif
   enddo
   nlast = np
  endif

  incr = 0
  !--If out > in: fill all the holes 
  if(Preg%out_tot > Preg%in_tot) then
   nlast = np
   do j = Preg%in_tot+1,Preg%out_tot
    do while (particule(nlast)%dir /= 0_iXb)
     nlast = nlast-1
     incr = incr+1
    enddo
    indice  = index_exit(j)
    if(indice == 0)  then
     print *,"PROBLEME INDICE ",indice,j,index_exit(j-1); stop
    endif
    compt = compt + 1
    particule(indice) = particule(nlast)
    nlast = nlast -1
   enddo
  endif
  nlast = nlast + incr

  np = nlast
  !--If particles are less than initially then set to zero
  !  older particles
  if(np<np_ini) call set_zero_particle(particule,np+1,np_ini)

  !--Set to zero the index_exit
  index_exit(:Preg%out_tot) = 0

  !--Control
  if(np /= (np_ini-Preg%np%out+Preg%np%in-Preg%pp%out+&
       &    Preg%pp%in-Preg%out_xm-Preg%out_xp)) stop "probleme particule"

#ifdef HAVE_DEBUG
  if(wrtscreen == 0) then 
   print '(" Moi processus ",i4," il reste ",i3," particules a sortir sur ",i10)',&
        mpiinfo%me,count(particule(1:np)%dir/=0_iXb),np

   print '("Moi processus ",i2," max & min ",6(f8.3))',mpiinfo%me,&
        &    minval(particule(1:np)%pos(1)),maxval(particule(1:np)%pos(1)), &
        &    minval(particule(1:np)%pos(2)),maxval(particule(1:np)%pos(2)), &
        &    minval(particule(1:np)%pos(3)),maxval(particule(1:np)%pos(3))
   if(count(particule(1:np)%pos(1)>xmax) > 0) then
    print '("Moi processus ",i5," j ai ",i5," particules tq x>xmax son indice est",i7," et np ",i7)', &
         mpiinfo%me,count((particule(1:np)%pos(1)>xmax).and.(particule(1:np)%exc==zero)),maxloc(particule(1:np)%pos(1)),np
    print '("Moi processus ",i5," j ai ",i5," pickups tq z>zmax son indice est",i7," et np ",i7)', &
         mpiinfo%me,count((particule(1:np)%pos(1)>xmax).and.(particule(1:np)%exc/=zero)),maxloc(particule(1:np)%pos(1)),np
   endif
   if(mpiinfo%me == 0) then
    write(*,*) "          VERIFICATION"
   endif
   print '("Moi processus ",i4,"j ai maintenant ",i10," part alors que j en attends ",i10)', &
        mpiinfo%me,np,np_ini-sum(np_out)+Preg%np%in&
        &  -sum(np_out_pick)+Preg%pp%in-Preg%out_xp-Preg%out_xm
  endif
#endif

  __GETTIME(38,2)!--Timer stop
  __WRT_DEBUG_OUT("rangement")
 end subroutine rangement
 !*********************** END RANGEMENT ******************************************

 !************************** PACK_COM ********************************************
 !--ce module package les particules a echange et les echanges
 subroutine pack_com(particule,Preg,index_exit,s_min,s_r,nptot,s_min_loc,s_max_loc)

  integer,intent(inout) :: nptot
  integer,intent(inout) :: index_exit(:)
  real(dp),intent(in) :: s_min(3),s_r(3),s_min_loc(3),s_max_loc(3)
  type(count_particle_type),intent(inout) :: Preg
  type(particletype),intent(inout) :: particule(:)


  integer,dimension(:),allocatable :: &
       &       p_out_N,p_out_NE,p_out_E,p_out_SE, &
       &       p_out_S,p_out_SW,p_out_W,p_out_NW, & 
       &       pick_out_N,pick_out_NE,pick_out_E,pick_out_SE, &
       &       pick_out_S,pick_out_SW,pick_out_W,pick_out_NW
  type(particletype),allocatable ::  p_in(:)

  integer :: jj
  integer,dimension(nb_voisins) :: pin_trace,pickin_trace
  character(len=500) :: msg

  __WRT_DEBUG_IN("pack_com")
  __GETTIME(34,1)!--Timer start

  !--on envoie le nombre de particule a echange aux processeurs voisins et on le receptionne
  call pre_communication(Preg%np%out_dir,Preg%np%in_dir,mpiinfo,etiquette1)

  !--on envoie le nombre de particule a echange aux processeurs voisins et on le receptionne
  call pre_communication(Preg%pp%out_dir,Preg%pp%in_dir,mpiinfo_pick,etiquette1)

#ifdef HAVE_DEBUG 
  write(msg,'(2(2a,i4.4,a,8i5,a))')&
       &" RCP  PARTICLES <=",&
       &" Process ",mpiinfo%me,&
       &" is exchanging ",Preg%np%in_dir,ch10,&
       &" RCP  PICKUPS <=",&
       &" Process ",mpiinfo%me,&
       &" is exchanging ",Preg%pp%in_dir,ch10
  call wrtout(6,msg,"PERS")
#endif

  Preg%np%in = sum(Preg%np%in_dir)
  Preg%pp%in = sum(Preg%pp%in_dir)
  Preg%in_tot = Preg%np%in + Preg%pp%in

  !--Allocation de tableaux temporaire pour l'envoi des particules
  allocate(p_out_N(Preg%np%out_dir(N))); allocate(p_out_NE(Preg%np%out_dir(NE)))
  allocate(p_out_E(Preg%np%out_dir(E))); allocate(p_out_SE(Preg%np%out_dir(SE)))
  allocate(p_out_S(Preg%np%out_dir(S))); allocate(p_out_SW(Preg%np%out_dir(SW))) 
  allocate(p_out_W(Preg%np%out_dir(W))); allocate(p_out_NW(Preg%np%out_dir(NW)))

  !--Allocation de tableaux temporaire pour l'envoi des particules
  allocate(pick_out_N(Preg%pp%out_dir(N))); allocate(pick_out_NE(Preg%pp%out_dir(NE)))
  allocate(pick_out_E(Preg%pp%out_dir(E))); allocate(pick_out_SE(Preg%pp%out_dir(SE)))
  allocate(pick_out_S(Preg%pp%out_dir(S))); allocate(pick_out_SW(Preg%pp%out_dir(SW))) 

  allocate(pick_out_W(Preg%pp%out_dir(W))); allocate(pick_out_NW(Preg%pp%out_dir(NW)))

  !--On compacte toutes les particules  et pickups qui doivent sortir dans des tabeaux temporaire
  call pack_part(Preg,index_exit,&
       &         p_out_N,p_out_NE,p_out_E,p_out_SE,& 
       &         p_out_S,p_out_SW,p_out_W,p_out_NW,&
       &         pick_out_N,pick_out_NE,pick_out_E,pick_out_SE,&
       &         pick_out_S,pick_out_SW,pick_out_W,pick_out_NW,&
       &         particule)

  !--allocation des tableaux temporaires pour a reception des
  !  particules et pickups
  allocate(p_in(Preg%in_tot)) ;  call set_zero_particle(p_in,1,Preg%in_tot)

  !--Construction of the reference array to indicate where the
  !  particules and pickups have to be receveid from other procs
  !  in the array p_in
  pin_trace(1) = 1
  pin_trace(2:) = (/(1+sum(Preg%np%in_dir(:jj)),jj=1,nb_voisins-1)/)
  pin_trace = min(Preg%in_tot,pin_trace)
  pickin_trace(1) = 1+Preg%np%in
  pickin_trace(2:) = (/(pickin_trace(1)+sum(Preg%pp%in_dir(:jj)),jj=1,nb_voisins-1)/)
  pickin_trace = min(Preg%in_tot,pickin_trace)

  !--On envoie les particules vers les processeur voisins et on les receptionne
  call communication(particule,&
       &             p_out_N,p_out_NE,p_out_E,p_out_SE,&
       &             p_out_S,p_out_SW,p_out_W,p_out_NW, &
       &             p_in,Preg%np%out_dir,Preg%np%in_dir,pin_trace,Preg%in_tot, &
       &             nb_voisins,mpiinfo%comm,mpiinfo%voisin)

  !--On envoie les particules vers les processeur voisins et on les receptionne
  call communication(particule,&
       &             pick_out_N,pick_out_NE,pick_out_E,pick_out_SE,&
       &             pick_out_S,pick_out_SW,pick_out_W,pick_out_NW, &
       &             p_in,Preg%pp%out_dir,Preg%pp%in_dir,pickin_trace,Preg%in_tot, &
       &             nb_voisins,mpiinfo_pick%comm,mpiinfo_pick%voisin)

  call rangement(Preg,index_exit,p_in, &
       &         particule,Preg%np%out_dir,Preg%pp%out_dir,nptot, &
       &         s_min,s_r,s_min_loc,s_max_loc &
#ifdef HAVE_DEBUG
       &         ,s_max(1) &
#endif
       )

  deallocate(p_out_N) 
  deallocate(p_out_NE)
  deallocate(p_out_E,p_out_SE)
  deallocate(p_out_S,p_out_SW) 
  deallocate(p_out_W,p_out_NW)
  deallocate(pick_out_N) 
  deallocate(pick_out_NE)
  deallocate(pick_out_E,pick_out_SE)
  deallocate(pick_out_S,pick_out_SW) 
  deallocate(pick_out_W,pick_out_NW)
  deallocate(p_in) 

  __GETTIME(34,2)!--Timer stop
  __WRT_DEBUG_OUT("pack_com")
 end subroutine pack_com
 !******************************* END PACK_COM *************************************


end module particle_com
