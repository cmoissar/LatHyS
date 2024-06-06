!=============================================================
!=============================================================
module m_reassemble_unform

 use defs_basis
 use defs_variable
 use defs_grid
 use defs_species
 use defs_arr3Dtype
 use defs_particletype
 use m_writeout
 use m_verify_reassemble_mpi

 implicit none
 private

#ifndef HAVE_NETCDF

 public ::                &
      reassemble_unform
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

 subroutine reassemble_unform(run_name)

  character(len=*),intent(in) :: run_name
  !--Reading variables
  integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
  integer :: rang,nb_procs,ndims,nb_voisin,is
  integer,allocatable :: dims(:),nptot_proc(:)
  integer,allocatable :: voisin_proc(:,:)
  integer,allocatable :: coord_proc(:,:)

  !--Variables to read particles diagnostic
  type(particletype),allocatable  :: particule_proc(:,:)

  !--Variables to read fields diagnostic
  real(dp),dimension(:,:,:,:),allocatable :: bx_proc,by_proc,bz_proc,&
       &                                     ex_proc,ey_proc,ez_proc,&
       &                                     dn_proc,ux_proc,uy_proc,&
       &                                     uz_proc,uxa_proc,uya_proc,&
       &                                     uza_proc,dna_proc

  !--Others
  integer :: iunit,nfl_w,ii
  integer :: ind,deby,debz,finy,finz,iproc
  integer :: deb,fin,cumul
  integer :: itmp(3)
  logical :: file_e
  character(len=50) :: write_name
  character(len=64) :: filename
  character(len=500) :: msg

  !--Create file name for Field (READ, proc 0) 
  write(filename,'(a3,i4.4,a1,a)')"c3_",0,'_',trim(run_name)

  !--Inquire if the file exists
  inquire( file=trim(filename), exist=file_e )
  !--No file: return
  if(.not.(file_e)) then
   call wrtout(6,"File: "//trim(filename)//" does not exits","PERS")
   return
  endif

  iunit = 1
  open(unit=iunit,file = filename,action = 'read', &
       status = 'old',form = 'unformatted')

  read(unit = iunit) rang
  read(unit = iunit) itmp(1)
  call set_grid(itmp,4)!--nxyzm
  read(unit = iunit) nb_procs,nb_voisin,itmp(1:3),ndims
  call set_grid(itmp,2)!--ncxm,ncym_proc,nczm_proc

  allocate(dims(ndims)) ; dims = 0
  allocate(coord_proc(nb_procs,ndims)) ; coord_proc =0
  allocate(nptot_proc(nb_procs)) ; nptot_proc = 0

  read(unit = iunit) dims
  read(unit = iunit) coord_proc(1,:)

  !--Ecriture des coordonnees des processus
  read(unit = iunit) itmp(:)
  call set_grid(itmp,3)!--ncxm,ncym,nczm 

  read(unit = iunit) itmp(:),gstep,s_min,s_max,s_min_loc,s_max_loc
  call set_grid(itmp,5)!--ncx,ncy,ncz
  read(unit = iunit) iter,dt,t,nptot_proc(1),ns,nfl_w

  allocate(voisin_proc(nb_procs,nb_voisin)) ; voisin_proc = 0
  read(unit = iunit) voisin_proc(1,1:nb_voisin)

  !--Allocate species and get informations
  call alloc_species(ns,Spe)  
  call species_get_var_bin(Spe,unit = iunit)

  !--allocation des tableaux de lectures du fichier diag de champ
  allocate(bx_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(by_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(bz_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(ex_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(ey_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(ez_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(dn_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(ux_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(uy_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(uz_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(dna_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(uxa_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(uya_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  allocate(uza_proc(nb_procs,ncm(1),ncm(2),ncm(3)))
  bx_proc = zero ; by_proc = zero ; bz_proc = zero
  ex_proc = zero ; ey_proc = zero ; ez_proc = zero
  ux_proc = zero ; uy_proc = zero ; uz_proc = zero
  uxa_proc = zero ; uya_proc = zero ; uza_proc = zero
  dn_proc = zero ; dna_proc = zero

  ! Lecture des dfonnées champ du processus 0
  read(unit = iunit) bx_proc(1,:,:,:), &
       &             by_proc(1,:,:,:), &
       &             bz_proc(1,:,:,:)

  read(unit = iunit) ex_proc(1,:,:,:), &
       &             ey_proc(1,:,:,:), &
       &             ez_proc(1,:,:,:)
  read(unit = iunit) dn_proc(1,:,:,:), &
       &             ux_proc(1,:,:,:), &
       &             uy_proc(1,:,:,:), &
       &             uz_proc(1,:,:,:)
  read(unit = iunit) dna_proc(1,:,:,:), &
       &             uxa_proc(1,:,:,:), &
       &             uya_proc(1,:,:,:), &
       &             uza_proc(1,:,:,:)

  close(unit = iunit)

  ! Maintenat que l'on connait le nombre de processus qui a tourne
  ! on lit les autres fichiers
  if (nb_procs>1) then
   do iproc=2,nb_procs 
    !--On convertit le rang_du processus en charactere    
    write(filename,'(a3,i4.4,a1,a)')"c3_",iproc-1,'_',trim(run_name)
    ! write(*,*) '************************************************'
    ! write(*,*) 'Lecture du fichier de simulations du processus ',iproc
    ! write(*,*)
    iunit = 1
    open(unit=iunit,file = filename, action = 'read', &
         status = 'old',form = 'unformatted')


    read(unit = iunit) rang

    read(unit = iunit) itmp(1)
    call set_grid(itmp,4)!--nxyzm
    read(unit = iunit) nb_procs,nb_voisin,itmp(:),ndims
    call set_grid(itmp,2)!--ncxm,ncym_proc,nczm_proc

    read(unit = iunit) dims
    read(unit = iunit) coord_proc(iproc,:)
    read(unit = iunit) itmp(:)
    call set_grid(itmp,3)!--ncxm,ncym,nczm 

    read(unit = iunit) itmp(:),gstep,s_min,s_max,s_min_loc,s_max_loc
    call set_grid(itmp,5)!--ncx,ncy,ncz
    read(unit = iunit) iter,dt,t,nptot_proc(1),ns,nfl_w

    read(unit = iunit) voisin_proc(iproc,1:nb_voisin)

    !--Get species informations
    call species_get_var_bin(Spe,unit = iunit)

    !--Lecture des dfonnées champ du processus 0
    read(unit = iunit) bx_proc(iproc,:,:,:), &
         by_proc(iproc,:,:,:), &
         bz_proc(iproc,:,:,:)
    read(unit = iunit) ex_proc(iproc,:,:,:), &
         ey_proc(iproc,:,:,:), &
         ez_proc(iproc,:,:,:)
    read(unit = iunit) dn_proc(iproc,:,:,:), &
         ux_proc(iproc,:,:,:), &
         uy_proc(iproc,:,:,:), &
         uz_proc(iproc,:,:,:)
    read(unit = iunit) dna_proc(iproc,:,:,:), &
         uxa_proc(iproc,:,:,:), &
         uya_proc(iproc,:,:,:), &
         uza_proc(iproc,:,:,:)

    close(unit = iunit)
   enddo
  endif

  nptot = sum(nptot_proc)


  !-- Verification des INTERFACES
  call verify_reassemble_mpi(bx_proc,by_proc,bz_proc,&
       &                     ex_proc,ey_proc,ez_proc,&
       &                     ux_proc,uy_proc,uz_proc,&
       &                     uxa_proc,uya_proc,uza_proc,&
       &                     dn_proc,dna_proc,voisin_proc,&
       &                     nb_procs)


  !**************************** On recolle les morceaux ************************
  ! a partir des tableaux de chaque processeur, et connaisant les coordonnées des 
  ! différents processuers dans la topologie de processus, on reconstruit les 
  ! tableaux globaux

  !--Allocation of global arrays
  call alloc_arr3D(Bfield,ncm_tot)
  call alloc_arr3D(Efield,ncm_tot)
  call alloc_arr3D(vel,ncm_tot)
  call alloc_arr3D(vela,ncm_tot)
  allocate(dn(ncm_tot(1),ncm_tot(2),ncm_tot(3)));   dn = zero 
  allocate(dna(ncm_tot(1),ncm_tot(2),ncm_tot(3))); dna = zero

  !--On reconstruit les tableaux ayant la meme grille , celle de B
  do ind = 1,nb_procs
   deby = coord_proc(ind,ndims-1)*(ncm(2)-2) + 1
   debz = coord_proc(ind,ndims)*(ncm(3)-2) + 1
   finy = deby + ncm(2) - 2
   finz = debz + ncm(3) - 2
   ! print '("P ",i3, " deby,finy ",2(i5)," debz,finz ",2(i5))',ind,deby,finy,debz,finz
   Bfield%x(1:ncm_tot(1),deby:finy,debz:finz) = bx_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
   Bfield%y(1:ncm_tot(1),deby:finy,debz:finz) = by_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
   Bfield%z(1:ncm_tot(1),deby:finy,debz:finz) = bz_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
   vel%x(1:ncm_tot(1),deby:finy,debz:finz) = ux_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
   vel%y(1:ncm_tot(1),deby:finy,debz:finz) = uy_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
   vel%z(1:ncm_tot(1),deby:finy,debz:finz) = uz_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
   vela%x(1:ncm_tot(1),deby:finy,debz:finz) = uxa_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
   vela%y(1:ncm_tot(1),deby:finy,debz:finz) = uya_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
   vela%z(1:ncm_tot(1),deby:finy,debz:finz) = uza_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
   dn(1:ncm_tot(1),deby:finy,debz:finz) = dn_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
   dna(1:ncm_tot(1),deby:finy,debz:finz) = dna_proc(ind,1:ncm(1),1:ncm(2)-1,1:ncm(3)-1)
  enddo

  !--on reconstruit les tableaux du champ électrique
  do ind = 1,nb_procs
   deby = coord_proc(ind,ndims-1)*(ncm(2)-2) + 1
   debz = coord_proc(ind,ndims)*(ncm(3)-2) + 1
   finy = deby + ncm(2) - 1
   finz = debz + ncm(3) - 1
   ! print '("P ",i3, " deby,finy ",2(i5)," debz,finz ",2(i5))', &
   !      ind,deby,finy,debz,finz
   Efield%x(1:ncm_tot(1),deby:finy,debz:finz) = ex_proc(ind,:,:,:)
   Efield%y(1:ncm_tot(1),deby:finy,debz:finz) = ey_proc(ind,:,:,:)
   Efield%z(1:ncm_tot(1),deby:finy,debz:finz) = ez_proc(ind,:,:,:)
  enddo

  !--Creation du fichier de diagnostique de champ a acces sequentiel
  write_name = "cw_"//trim(run_name)
  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")


  !--Open the file for write Particles
  open(UNIT= iunit,FILE = write_name,FORM = 'unformatted',ACTION = 'write', &
       STATUS = 'unknown')

  write(iunit) itmp(:),gstep,s_min,s_max &
       &       ,iter,dt,t,nptot,ns
  call set_grid(itmp,5)!--ncx,ncy,ncz

  call species_put_var_bin(Spe,unit=iunit)

  write(iunit) Bfield%x(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  write(iunit) Bfield%y(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  write(iunit) Bfield%z(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  write(iunit) Efield%x(:,:,:)
  write(iunit) Efield%y(:,:,:)
  write(iunit) Efield%z(:,:,:)
  write(iunit)     dn(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  write(iunit)  vel%x(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  write(iunit)  vel%y(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  write(iunit)  vel%z(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  write(iunit)    dna(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  write(iunit) vela%x(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  write(iunit) vela%y(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  write(iunit) vela%z(1:nc_tot(1)+1,1:nc_tot(2)+1,1:nc_tot(3)+1)
  close(iunit)

  !--Desallocation des tableaux
  deallocate(bx_proc,by_proc,bz_proc)
  deallocate(ex_proc,ey_proc,ez_proc)
  deallocate(ux_proc,uy_proc,uz_proc)
  deallocate(uxa_proc,uya_proc,uza_proc)
  deallocate(dn_proc,dna_proc)
  deallocate(dn,dna) ;
  call dealloc_arr3D(Bfield) ; call dealloc_arr3D(Efield)
  call dealloc_arr3D(vel)
  call dealloc_arr3D(vela)

  !================================== lecture du fichier de particules ================
  ! nom du fichier de lec ture des données de la simulation pour les diagnostiques des 
  ! particules

  !--Create file name to read Particles
  write(filename,'(a3,i4.4,a1,a)')"p3_",0,'_',trim(run_name)
  
  !--Inquire if the file exists
  inquire( file=trim(filename), exist=file_e )
  !--No file: return
  if(.not.(file_e)) then
   call wrtout(6,"File: "//trim(filename)//" does not exits","PERS")
   return
  endif

  iunit = 1
  open(unit=iunit,file = filename, action = 'read', &
       status = 'old',form = 'unformatted')

  dims = 0; coord_proc =0; nptot_proc = 0

  read(unit = iunit) rang
  read(unit = iunit) itmp(1)
  call set_grid(itmp,4)!--nxyzm
  read(unit = iunit) nb_procs,nb_voisin,itmp(:),ndims
  call set_grid(itmp,2)!--ncxm,ncym_proc,nczm_proc

  read(unit = iunit) dims

  !--Ecriture du nombre de particules dans  le processus
  read(unit = iunit) nptot_proc(1)

  read(unit = iunit) itmp(:2)
  call set_grid(itmp,3)!--ncym,nczm 

  read(unit = iunit) itmp(:),gstep,s_min,s_max,iter,dt,t,&
       nptot_proc(1),ns,nfl_w,npm
  call set_grid(itmp,5)!--ncx,ncy,ncz
  read(unit = iunit) voisin_proc(1,1:nb_voisin)

  !--Set to zero some infos
  dims = 0; coord_proc =0; nptot_proc = 0

  !--Get species informations
  call species_get_var_bin(Spe,unit = iunit)

  !--Allocation du tableaux de particules pour chaque processuers
  allocate(particule_proc(nb_procs,npm)) 
  call set_zero_particle(particule_proc)

  call get_var_particle_bin(particule_proc(1,:),nptot_proc(1),unit=iunit)

  close (unit = iunit)

  !--Print some information
  write(msg,'(81a,(2a,i8),2(2a,i4),2(2a,3i4),(2a,i4),(2a,i2,a,i2,a))')&
       &ch10," ",(/("-",ii=1,79)/),&
       &ch10," Nmax particles per proc         : ",npm,&
       &ch10," N procs                         : ",nb_procs,&
       &ch10," N neighbours per proc           : ",nb_voisin,&
       &ch10," N cells in the box (X,Y,Z)      : ",ncm_tot(1),ncm_tot(2),ncm_tot(3),&
       &ch10," N cells per proc  (X,Y,Z)       : ",ncm,&
       &ch10," Topology                        : ",ndims,& 
       &ch10," Proc 0 coordinates              : (",coord_proc(1,1),",",coord_proc(1,2),")"
  call wrtout(6,msg,'PERS')                   

  write(msg,'((2a,3f8.3),2a,2(2a,3f8.3),2(2a,f8.3),(2a,i10))')&
       &ch10," grid step (gstep)               : ",gstep,&
       &ch10," Box Dimensions                  : ",&
       &ch10,"                Min              : ",s_min,&
       &ch10,"                Max              : ",s_max,&
       &ch10," Time step                       : ",dt,&
       &ch10," Diagnostic Time                 : ",t,&
       &ch10," N particles                     : ",nptot
  call wrtout(6,msg,'PERS')

  do is = 1,ns
   write(msg,'((81a),a,1x,a,(2a,2i10),2(2a,2f10.6))')&
        &ch10," ",(/("-",ii=1,79)/),&
        &ch10,Spe%S(is)%name,&
        &ch10,"       Start, Stop indices       : ",Spe%S(is)%n1,Spe%S(is)%n2,&
        &ch10,"       Thermal Speed (para,perp) : ",Spe%S(is)%vth1,Spe%S(is)%vth2,&
        &ch10,"       Charge, Masse             : ",Spe%S(is)%sq,Spe%S(is)%sm
   call wrtout(6,msg,'PERS')

   write(msg,'(1(2a,i6),3(2a,f10.6))')&
        &ch10,"       N Protons per cell        : ",Spe%S(is)%ng,&
        &ch10,"       % Alphas in solar wind    : ",Spe%S(is)%percent,& 
        &ch10,"       Charges Ratio             : ",Spe%S(is)%rcharge,& 
        &ch10,"       Masses Ratio              : ",Spe%S(is)%rmass
   call wrtout(6,msg,'PERS')                   

  enddo
  write(msg,'((81a),2a,(2a,3f8.3),(2a,f8.3))')&
       &ch10," ",(/("-",ii=1,79)/),&
       &ch10," Planet",&
       &ch10,"       Center                    : ",Spe%P%centr,&
       &ch10,"       Desplacement Speed        : ",Spe%P%speed
  call wrtout(6,msg,'PERS')

  if (nb_procs > 1) then
   do iproc=2,nb_procs

    !--Create file name
    write(filename,'(a3,i4.4,a1,a)')"p3_",iproc-1,'_',trim(run_name)
    iunit = 1
    open(unit=iunit,file = filename, action = 'read', &
         status = 'old',form = 'unformatted')

    read(unit = iunit) rang
    read(unit = iunit) itmp(1)
    call set_grid(itmp,4)!--nxyzm
    read(unit = iunit) nb_procs,nb_voisin,itmp(:),ndims
    call set_grid(itmp,2)!--ncxm,ncym_proc,nczm_proc

    read(unit = iunit) dims
    read(unit = iunit) nptot_proc(iproc)

    read(unit = iunit) itmp(:2)
    call set_grid(itmp,3)!--ncym,nczm 

    read(unit = iunit) itmp(:),gstep,s_min,s_max,iter,dt,t,&
         nptot_proc(iproc),ns,nfl_w,npm
    call set_grid(itmp,5)!--ncx,ncy,ncz
    read(unit = iunit) voisin_proc(iproc,1:nb_voisin)

    call species_get_var_bin(Spe,unit = iunit)

    call get_var_particle_bin(particule_proc(iproc,:),nptot_proc(iproc),unit=iunit)
 
    close (unit = iunit)
   enddo
  endif

  ! On recolle les morceaux en mettant outes les infos des différents processus
  ! dans un seul tableaux
  ! On boucle sur les processeurs
  allocate(particule(npm*nb_procs))
  call set_zero_particle(particule,1,npm*nb_procs)

  cumul = 0
  do iproc = 1,nb_procs
   deb = cumul + 1
   fin = cumul + nptot_proc(iproc)
   cumul = fin
   ! print *,'deb et fin ',deb,fin
   ! print *, 'Cumul du nombre de particule par processeur',cumul
   particule(deb:fin) = particule_proc(iproc,1:nptot_proc(iproc))
  enddo
  nptot = sum(nptot_proc)

  !Ecriture dans un fichier a aces sequentiel
  write_name = 'pw_'//run_name

  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")

  open(UNIT= iunit,FILE = write_name,FORM = 'unformatted',ACTION = 'write', &
       STATUS = 'unknown')

  write(iunit) nc_tot,gstep,s_min,s_max &
       ,iter,dt,t,nptot,ns,npm

  call species_put_var_bin(Spe,unit = iunit)

  call put_var_particle_bin(particule,nptot,unit=iunit)

  close(iunit)

  call dealloc_species(Spe)
  deallocate(dims,coord_proc)
  deallocate(nptot_proc)
  deallocate(particule_proc,particule)

 end subroutine reassemble_unform
#endif

end module m_reassemble_unform
