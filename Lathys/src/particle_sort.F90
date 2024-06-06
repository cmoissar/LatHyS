module particle_sort

use defs_particletype
use defs_grid,only : ncm
use defs_variable,only : s_min_loc
use defs_parametre,only : gstep
use m_writeout

#include "q-p_common.h"

implicit none
private

real,public,allocatable,save :: sort_tmp(:,:,:,:)
integer,public,allocatable,save :: sort_indice(:)

public :: p_sort

contains

subroutine p_sort(particle,nptot)
type(particletype),dimension(:),intent(inout) :: particle
integer,intent(in)::nptot
integer::POS,I,J,k,ip,lp,tmp,tmp2
integer,dimension(1:nptot)::ind
integer,dimension(1:(ncm(2)*ncm(3)+1))::first_p,last_p
integer :: y,z
  character(len=100) :: msg

type(particletype)::BUF1,BUF2

! write(msg,'(a,i10)')&
!       & " start ",int(0)
!    call wrtout(6,msg,'PERS')

if (.NOT.(Allocated(sort_tmp))) then
        allocate(sort_tmp(1:8,1:2,1:2,1:ncm(1)))
        allocate(sort_indice(1:ncm(2)*ncm(3)))
endif

sort_indice(:)=0
do I=1,nptot
y=INT((particle(i)%pos(2)-s_min_loc(2))/gstep(2))+1
z=INT((particle(i)%pos(3)-s_min_loc(3))/gstep(3))
J=z*ncm(2)+y
 if ((j.le.0).or.(j.gt.(ncm(2)*ncm(3)))) then
 write(msg,'(a,i10)')&
       & " J ",int(J)
    call wrtout(6,msg,'PERS')
write(msg,'(a,i10)')&
       & " y ",INT((particle(i)%pos(2)-s_min_loc(2))/gstep(2))+1
    call wrtout(6,msg,'PERS')
write(msg,'(a,i10)')&
       & " ny ",int(ncm(2))
    call wrtout(6,msg,'PERS')
write(msg,'(a,i10)')&
       & " z ",INT((particle(i)%pos(3)-s_min_loc(3))/gstep(3))
    call wrtout(6,msg,'PERS')
write(msg,'(a,i10)')&
       & " nz ",int(ncm(3))
    call wrtout(6,msg,'PERS')
write(msg,'(a,i10)')&
       & " origin ",particle(i)%orig
    call wrtout(6,msg,'PERS')
                stop
 endif
 sort_indice(J)=sort_indice(J)+1 ! sort_indice(J,1)= nombre de particules dans la case J
 ind(I)=J  ! ind(I)= case dans laquelle se trouve la particule I
ENDDO
 first_p(1)=1 !first_p(I) = premiere position dans la liste des particules de la case I
DO I=2,ncm(2)*ncm(3)+1
 first_p(I)=first_p(I-1)+sort_indice(I-1)
ENDDO
 last_p(:)=first_p(:)  !last_p(I) = position dans la liste de la derniere particule classee parmis
            ! les particules de la case I
 POS=1     ! on classe la POSieme particule
 IP=1      ! IP = case en cours
 LP=first_p(2) ! LP = derniere position+1 des particules de la case en cours
DO
  if (POS.ge.LP) then   !changement de case on recalcul IP et LP
    DO  ! on parcours les cases jusqu'a ce qu'on en trouve une pas classee
        IP=IP+1
        LP=first_p(IP+1)
        POS=last_p(IP)          ! POS prend la valeur de la prochaine particule pas encore classee
        if (POS.gt.nptot) EXIT
        if (POS.lt.LP) EXIT
    ENDDO
  ENDIF
  if (POS.gt.nptot) EXIT
  I=ind(POS)        ! i = case dans laquelle se trouve la POSieme particule
  IF (I.EQ.IP) then    ! si elle est deja dans la bonne case on y touche pas
     POS=POS+1
  else               !sinon on la classe
     tmp=ind(pos)
     BUF1=particle(pos)
     K=POS+1               !K parcours la liste pour trouver la bonne place
     DO
       DO
           J=ind(last_p(I))    ! case dans laquelle se trouve la premiere particule pas classee 
                                        ! du domaine reserve a la classe I
           IF(J.NE.I) EXIT      ! si cette particule n'est pas dans la bonne case on continue
           last_p(I)=last_p(I)+1    ! sinon on regarde la particules pas classee suivante
       ENDDO
       BUF2=particle(last_p(i))! on intervertit les particules
       tmp2=ind(last_p(i))
       ind(last_p(i))=tmp
       tmp=tmp2
       particle(last_p(I))=BUF1
       last_p(I)=last_p(I)+1
       if (last_p(I).EQ.K) EXIT
       BUF1=BUF2
       I=J
    ENDDO
  endif
ENDDO
end subroutine p_sort

end module particle_sort
