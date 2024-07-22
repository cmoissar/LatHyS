!=============================================================
!=============================================================
module m_verify_reassemble_mpi

 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype

 implicit none
 private

 public ::                &
      verify_reassemble_mpi
contains
 !!############################################################

 !********************************************************************
 ! Auteur				:	 MMancini RModolo
 ! Date					:	 23/06/03
 ! Instituion				:	CETP/CNRS/IPSL
 ! Dernière modification		:	19/07/10	
 ! Résumé	
 ! 
 !********************************************************************

 subroutine verify_reassemble_mpi(&
      &                    bx_proc,by_proc,bz_proc,&
      &                    ex_proc,ey_proc,ez_proc,&
      &                    ux_proc,uy_proc,uz_proc,&
      &                    uxa_proc,uya_proc,uza_proc,&
      &                    dn_proc,dna_proc,voisin_proc,&
      &                    nb_procs)


  !--Variables to read fields diagnostic
  integer,intent(in) :: nb_procs
  integer,dimension(:,:),intent(in) :: voisin_proc
  real(dp),dimension(:,:,:,:),intent(in) :: bx_proc,by_proc,bz_proc,&
       &                                    ex_proc,ey_proc,ez_proc,&
       &                                    ux_proc,uy_proc,uz_proc,&
       &                                    uxa_proc,uya_proc,uza_proc,&
       &                                    dn_proc,dna_proc

  integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8

  !--Others
  integer :: iproc
  integer :: indx,indy,indz
  real(dp) :: diff

  !--Verification de la zone de recouvrement Nord/Sud
  do iproc = 1,nb_procs
   do indx = 1,ncm(1)-1
    do indy = 1,ncm(2)-1
     
     diff = abs(bx_proc(iproc,indx,indy,1)-bx_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Bx  N/S',iproc,indx,indy,diff
      !if(diff > tol6) stop 4
     endif

     diff = abs(by_proc(iproc,indx,indy,1)-by_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob By  N/S',iproc,indx,indy,diff
      !if(diff > tol7) stop 4
     endif

     diff = abs(bz_proc(iproc,indx,indy,1)-bz_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Bz  N/S',iproc,indx,indy,diff
      !if(diff > tol7) stop 4
     endif


     diff = abs(dn_proc(iproc,indx,indy,1)-dn_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Dn  N/S',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ux_proc(iproc,indx,indy,1)-ux_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Ux  N/S',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(uy_proc(iproc,indx,indy,1)-uy_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Uy  N/S',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(uz_proc(iproc,indx,indy,1)-uz_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Uz  N/S',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(dna_proc(iproc,indx,indy,1)-dna_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Dna  N/S',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(uxa_proc(iproc,indx,indy,1)-uxa_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Uxa  N/S',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(uya_proc(iproc,indx,indy,1)-uya_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Uya  N/S',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(uza_proc(iproc,indx,indy,1)-uza_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Uza  N/S',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

    enddo
   enddo
  enddo

  ! Verification de la zone de recouvrement Est/Ouest
  do iproc = 1,nb_procs
   do indx = 1,ncm(1)-1
    do indz = 1,ncm(3)-1
     
     diff = abs(bx_proc(iproc,indx,1,indz)-bx_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Bx  E/W',iproc,indx,indz,diff
      !if(diff > tol6) stop 4
     endif

     diff = abs(by_proc(iproc,indx,1,indz)-by_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob By  E/W',iproc,indx,indz,diff
      !if(diff > tol7) stop 4
     endif

     diff = abs(bz_proc(iproc,indx,1,indz)-bz_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Bz  E/W',iproc,indx,indz,diff
      !if(diff > tol7) stop 4
     endif

     diff = abs(dn_proc(iproc,indx,1,indz)-dn_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Dn  E/W',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ux_proc(iproc,indx,1,indz)-ux_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Ux  E/W',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(uy_proc(iproc,indx,1,indz)-uy_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Uy  E/W',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(uz_proc(iproc,indx,1,indz)-uz_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Uz  E/W',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(dna_proc(iproc,indx,1,indz)-dna_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Dna  E/W',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(uxa_proc(iproc,indx,1,indz)-uxa_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Uxa  E/W',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(uya_proc(iproc,indx,1,indz)-uya_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Uya  E/W',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(uza_proc(iproc,indx,1,indz)-uza_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Uza  E/W',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

    enddo
   enddo
  enddo

  ! Pour le champ Ã©lectrique

  ! Verification Interface Nord/Sud
  do iproc = 1,nb_procs
   do indx = 1,ncm(1)
    do indy = 1,ncm(2)


     diff = abs(ex_proc(iproc,indx,indy,2)-ex_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)))
     if(diff/=zero) then
      print *,'Prob Ex  N/S 1',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ex_proc(iproc,indx,indy,1)-ex_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Ex  N/S 2',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ey_proc(iproc,indx,indy,2)-ey_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)))
     if(diff/=zero) then
      print *,'Prob Ey  N/S 1',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ey_proc(iproc,indx,indy,1)-ey_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Ey  N/S 2',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ez_proc(iproc,indx,indy,2)-ez_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)))
     if(diff/=zero) then
      print *,'Prob Ez  N/S 1',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ez_proc(iproc,indx,indy,1)-ez_proc(voisin_proc(iproc,S)+1,indx,indy,ncm(3)-1))
     if(diff/=zero) then
      print *,'Prob Ez  N/S 2',iproc,indx,indy,diff
      if(diff > tol7) stop 4
     endif

    enddo
   enddo
  enddo

  ! Verification Interface Es/Ouest

  do iproc = 1,nb_procs
   do indx = 1,ncm(1)
    do indz = 1,ncm(3)
     
     diff = abs(ex_proc(iproc,indx,2,indz)-ex_proc(voisin_proc(iproc,W)+1,indx,ncm(2),indz))
     if(diff/=zero) then
      print *,'Prob Ex  E/W 1',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ex_proc(iproc,indx,1,indz)-ex_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Ex  E/W 2',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ey_proc(iproc,indx,2,indz)-ey_proc(voisin_proc(iproc,W)+1,indx,ncm(2),indz))
     if(diff/=zero) then
      print *,'Prob Ey  E/W 1',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ey_proc(iproc,indx,1,indz)-ey_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Ey  E/W 2',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ez_proc(iproc,indx,2,indz)-ez_proc(voisin_proc(iproc,W)+1,indx,ncm(2),indz))
     if(diff/=zero) then
      print *,'Prob Ez  E/W 1',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

     diff = abs(ez_proc(iproc,indx,1,indz)-ez_proc(voisin_proc(iproc,W)+1,indx,ncm(2)-1,indz))
     if(diff/=zero) then
      print *,'Prob Ez  E/W 2',iproc,indx,indz,diff
      if(diff > tol7) stop 4
     endif

    enddo
   enddo
  enddo


 end subroutine verify_reassemble_mpi

end module m_verify_reassemble_mpi
