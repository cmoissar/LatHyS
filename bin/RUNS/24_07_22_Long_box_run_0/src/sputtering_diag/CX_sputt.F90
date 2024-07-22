!=============================================================
!=============================================================
module CX_sputt

 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype
 use defs_particletype
 use defs_species
 use defs_parametre,only : fildat,dt,ns,npm,gstep
 use m_writeout
 use ieee_arithmetic
 #ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
#endif
#include "q-p_common.h"

 implicit none
 private 

#ifdef HAVE_NETCDF
 private ::		  &
      merge_field_cdf_pro,&
      create_file_diag

 public ::                &
      merge_CXprod
contains
 !!############################################################

 !********************************************************************
 ! Auteur				:	 MMancini RModolo
 ! Date					:	 23/06/03
 ! Institution				:	CETP/CNRS/IPSL
 ! Derniere modification		:	19/07/10	
 ! Resume	
 ! 
 !********************************************************************
 
 
   !!############################################################
   !!=============================================================
   !!subroutine: m_split_reassemble/merge_photoprod
   !! This routine reads one splited file and merge it 
   !! accoriding to processes position
   
    subroutine merge_CXprod(run_name)
   ! this routine merge the file created during the post_treatment to 1 file
   ! 1 file O+ photoproduction : including only O+ photoproduction arrays
   ! after the creation for each process it recombines in one unique file
   character(len=*),intent(in) :: run_name
   !Local Variable
   character(len=64) :: var_name1,var_name2,var_name3
   character(len=4)  :: prefix,prefix_output,prefixb
   
   prefix = "Mom3"
   prefixb = "Atm3"
   prefix_output = "CX_O"
   !var_name1 = "Prod_O+"
 
    call merge_field_cdf_pro(run_name,prefix,prefixb,prefix_output) 
   
   
  end subroutine merge_CXprod
 
    
    
    

!!############################################################
!!****************************************************
 !!-----------------------------------------------
 !! m_split_reassemble_cdf.F90
 !! reads 1 file computed for each process
 !! and merge it to one big file
 
  subroutine merge_field_cdf_pro(run_name,prefix,prefixb,prefix_output)
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: prefix,prefixb,prefix_output
    !--Reading variables
    integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
    integer :: rang,nb_procs,ndims,nb_voisin
    integer,allocatable :: dims(:),nptot_proc(:),dimspec(:)
    integer,allocatable :: voisin_proc(:,:)
    integer,allocatable :: coord_proc(:,:),proc_tab(:,:)
    character(len=10),dimension(:),allocatable:: ion_label    
    
      
    !--Variables to read fields diagnostic
    real(dp),dimension(:,:,:,:),allocatable :: A_proc,Ax_proc,Ay_proc,Az_proc
    real(dp),dimension(:,:,:),allocatable :: Qx,A
    
    !--Others
    integer :: ind,deby,debz,finy,finz,iproc,ish
    integer :: deb,fin,cumul
    integer  :: stId,ncid,ii,is,ny,nz,jj,ind_i,ind_j,ind_k,kk
    integer  :: varid(105),dimid(11)
    integer  :: itmp(3)
    real(dp) :: rtmp(3)
    logical :: file_e
    character(len=50) :: write_name
    character(len=64) :: filename
     character(len=100) :: test_name
    character(len=500) :: msg
    integer :: nb_spe,i_spe,n_spe,n_exc
    character(len=11):: nam,nam_tmp
    character(len=11),dimension(:),allocatable:: name_dens
    real(dp) :: sigma_Hp_O
 
   
    !--Create file name for Field (READ, proc 0) 
    write(filename,'(a4,a1,i3.3,a1,2a)')trim(prefix),"_",0,'_',trim(run_name),".nc"
    write(*,*) 'merging file'  ,filename  
  
    !--Inquire if the file exists
    inquire( file=trim(filename), exist=file_e )
    !--No file: return
    if(.not.(file_e)) then
     call wrtout(6,"File: "//trim(filename)//" does not exits","PERS")
     return
    endif

    !--Open NetCDF file
    stId = nf90_open(trim(filename), nf90_nowrite, ncid)
    call test_cdf(stId)

    !--Control the Date
    stId = nf90_get_att(ncid, nf90_global, "Date", test_name)
    call test_cdf(stId)
    if(fildat /=trim(test_name))  stop "Problem of Dates"

    !--Get the number of processus which generated the file
    call get_simple_variable_cdf(ncid,"nproc",nb_procs)
    write(*,*) 'nb_procs : ',nb_procs
  
    !--Get the rang of the open file 
    call get_simple_variable_cdf(ncid,"mpiinfo_me",rang)
  
    !--Get topology dimension
    call get_simple_dimens_cdf(ncid,"mpi_ndims",ndims)
  
    !--Get number of neighbours
    call get_simple_dimens_cdf(ncid,"mpi_nb_voisin",nb_voisin)
  
    !--Get number of cells in X,Y,Z
    call get_simple_variable_cdf(ncid,"nc_tot",itmp(:))
    call set_grid(itmp,5)
  
    !--Get number of point for fields in Y,Z
    call get_simple_variable_cdf(ncid,"ncm_tot",itmp(:))
    call set_grid(itmp,3) 
  
    !--Get number point per processus
    call get_simple_variable_cdf(ncid,"ncm",itmp(:))
    call set_grid(itmp,2)
    
    !--Get total number of points
    call get_simple_variable_cdf(ncid,"nxyzm",itmp(1))
    call set_grid(itmp,4)
  
    allocate(dims(ndims)) ; dims = 0
  
    !--Define coordinate vector for proc
    allocate(coord_proc(nb_procs,ndims)) ; coord_proc =0
  
    !--Get the grid step
    call get_simple_variable_cdf(ncid,"gstep",gstep)
  
    !--Get the dimensions of the box
    !--MIN
    call get_simple_variable_cdf(ncid,"s_min",s_min)
    !--MAX
    call get_simple_variable_cdf(ncid,"s_max",s_max)
  
    !--Get Time Informations
    !--Iteration
    call get_simple_variable_cdf(ncid,"iter",iter)
    !--Time interval
    call get_simple_variable_cdf(ncid,"dt_t",rtmp(:2))
    dt = rtmp(1);   t = rtmp(2)
  
    !--Get Total Number of Particles for this proc
    allocate(nptot_proc(nb_procs)) ; nptot_proc = 0
    call get_simple_variable_cdf(ncid,"nptot",nptot_proc(1))
  
    !--Get Number of Species
    call get_simple_variable_cdf(ncid,"ns",ns)
  
    !--Allocate species and get informations
    call alloc_species(ns,Spe)
    call species_get_var_cdf(Spe,ncid)
    
    call get_simple_variable_cdf(ncid,"nb_ion_species",nb_spe)
    write(*,*) 'nb_species ',nb_spe
  
  
        
    !-- allocation du tableaux de nom des especes
    allocate(ion_label(nb_spe))	;	ion_label(:) = ''
    call get_simple_variable_cdf(ncid,"ion_label",ion_label)
    write(*,*) 'Ion label',ion_label    
  
     
    allocate(Qx(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
    Qx = zero	  
       
    allocate(voisin_proc(nb_procs,nb_voisin)) ; voisin_proc = 0
  
    stId = nf90_close(ncid);  call test_cdf(stId)
    
 ! Maintenat que l'on connait le nombre de processus qui a tourne
   ! on lit les autres fichiers
   ny=0
   do iproc=1,nb_procs 
     write(filename,'(a4,a1,i3.3,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
     write(*,*) 'read tmp file ',filename
     stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
     	call test_cdf(stId)
     	call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
     	if ((coord_proc(iproc,ndims-1)).gt.ny) ny=coord_proc(iproc,ndims-1)
     stId = nf90_close(ncid);  call test_cdf(stId)
     write(*,*) 'close file ',filename
   enddo
   ny=ny+1
   nz=nb_procs/ny
  allocate(proc_tab(ny,nz))
  do iproc=1,nb_procs 
     write(filename,'(a4,a1,i3.3,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
     write(*,*) 'read tmp file ',filename
     stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
     	call test_cdf(stId)
     	call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
     proc_tab(coord_proc(iproc,ndims-1)+1,coord_proc(iproc,ndims)+1)=iproc
     stId = nf90_close(ncid);  call test_cdf(stId)
     write(*,*) 'close file ',filename
   enddo

  ! Maintenat que l'on connait le nombre de processus qui a tourne
   ! on lit les autres fichiers
        debz = 1
  do jj=1,nz
	deby = 1
  do ii=1,ny
    iproc=proc_tab(ii,jj)
    !--Create file name
    write(filename,'(a4,a1,i3.3,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
    write(*,*) 'read file ',filename
    stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
    	call test_cdf(stId)
    	call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
    	!--Get number point per processus
    	call get_simple_variable_cdf(ncid,"ncm",itmp(:))
    	call set_grid(itmp,2) 
    	
    	nam_tmp = "Prod_CX_O+"
    	

        !--allocation des tableaux de lectures du fichier diag de moment
        allocate(A_proc(ncm(1),ncm(2),ncm(3),nb_spe))
        allocate(Ax_proc(ncm(1),ncm(2),ncm(3),nb_spe))
        allocate(Ay_proc(ncm(1),ncm(2),ncm(3),nb_spe))
        allocate(Az_proc(ncm(1),ncm(2),ncm(3),nb_spe))
        A_proc = zero; Ax_proc = zero;	Ay_proc = zero;	Az_proc = zero 
          
    
    !--Get Total Number of Particles
    call get_simple_variable_cdf(ncid,"nptot",nptot_proc(iproc))
    call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )

      finy = deby + ncm(2) - 2
      finz = debz + ncm(3) - 2
      
 
     call get_simple_variable_cdf(ncid,"density_species" ,A_proc(:,:,:,:) )
     call get_simple_variable_cdf(ncid,"Vx_species" ,Ax_proc(:,:,:,:) )
     call get_simple_variable_cdf(ncid,"Vy_species" ,Ay_proc(:,:,:,:) )
     call get_simple_variable_cdf(ncid,"Vz_species" ,Az_proc(:,:,:,:) )
     ish = 1 ! we are interstering only for H+sw =>O+
     
      Qx(1:ncm_tot(1),deby:finy,debz:finz) = Qx(1:ncm_tot(1),deby:finy,debz:finz)+ &
      &		A_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,ish)*&
      &	sqrt(Ax_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,ish)**2 + &
      &	     Ay_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,ish)**2 + &
      &	     Az_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,ish)**2)
      
      ! we add the production from planetary proton
      ish = 6
      
      do ind_i=1,ncm_tot(1)
        do ind_j = 1,ncm(2)-1
          do ind_k = 1,ncm(3)-1
            if (A_proc(ind_i,ind_j,ind_k,ish) > 1.e-2) then      
              Qx(ind_i,deby+ind_j-1,debz+ind_k-1) = Qx(ind_i,deby+ind_j-1,debz+ind_k-1)+ &
      	&	A_proc(ind_i,ind_j,ind_k,ish)*&
      	&	sqrt(Ax_proc(ind_i,ind_j,ind_k,ish)**2 + &
      	&	     Ay_proc(ind_i,ind_j,ind_k,ish)**2 + &
      	&	     Az_proc(ind_i,ind_j,ind_k,ish)**2)      
      	    endif
      	  enddo
      	enddo
      enddo 	
      	
      
      


   deby=deby+(ncm(2)-2)
   deallocate(A_proc,Ax_proc,Ay_proc,Az_proc)
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename
   enddo  
   debz=debz+(ncm(3)-2)   
   enddo
   nptot = sum(nptot_proc)  
   deallocate(proc_tab)    
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we read neutral density and charge exchange cross-section
   !--Create file name for Field (READ, proc 0) 
    write(filename,'(a4,a1,i3.3,a1,2a)')trim(prefixb),"_",0,'_',trim(run_name),".nc"
    write(*,*) 'merging file'  ,filename  
  
    !--Inquire if the file exists
    inquire( file=trim(filename), exist=file_e )
    !--No file: return
    if(.not.(file_e)) then
     call wrtout(6,"File: "//trim(filename)//" does not exits","PERS")
     return
    endif

    !--Open NetCDF file
    stId = nf90_open(trim(filename), nf90_nowrite, ncid)
    call test_cdf(stId)

    !--Control the Date
    stId = nf90_get_att(ncid, nf90_global, "Date", test_name)
    call test_cdf(stId)
    if(fildat /=trim(test_name))  stop "Problem of Dates"

    !--Get the number of processus which generated the file
    call get_simple_variable_cdf(ncid,"nproc",nb_procs)
    write(*,*) 'nb_procs : ',nb_procs
  
    !--Get the rang of the open file 
    call get_simple_variable_cdf(ncid,"mpiinfo_me",rang)
  
    !--Get topology dimension
    call get_simple_dimens_cdf(ncid,"mpi_ndims",ndims)
  
    !--Get number of neighbours
    call get_simple_dimens_cdf(ncid,"mpi_nb_voisin",nb_voisin)
  
    !--Get number of cells in X,Y,Z
    call get_simple_variable_cdf(ncid,"nc_tot",itmp(:))
    call set_grid(itmp,5)
  
    !--Get number of point for fields in Y,Z
    call get_simple_variable_cdf(ncid,"ncm_tot",itmp(:))
    call set_grid(itmp,3) 
    allocate(A(ncm_tot(1),ncm_tot(2),ncm_tot(3))); A = zero
  
    !--Get number point per processus
    call get_simple_variable_cdf(ncid,"ncm",itmp(:))
    call set_grid(itmp,2)
    
    !--Get total number of points
    call get_simple_variable_cdf(ncid,"nxyzm",itmp(1))
    call set_grid(itmp,4)
  
    !allocate(dims(ndims)) ; 
    dims = 0
  
    !--Define coordinate vector for proc
    !allocate(coord_proc(nb_procs,ndims)) ; 
    coord_proc =0
  
  
    !--Get Total Number of Particles for this proc
    !allocate(nptot_proc(nb_procs)) ; 
        nptot_proc = 0
    call get_simple_variable_cdf(ncid,"nptot",nptot_proc(1))
  
    !--Get Number of Species
    call get_simple_variable_cdf(ncid,"ns",ns)
  
    !--Allocate species and get informations
    call alloc_species(ns,Spe)
    call species_get_var_cdf(Spe,ncid)
  
    call get_simple_variable_cdf(ncid,"n_species",n_spe)
    call get_simple_variable_cdf(ncid,"n_reac_exc",n_exc)
  
    allocate(name_dens(n_spe))
    !--allocation des tableaux de lectures du fichier diag de champ
     
    !allocate(voisin_proc(nb_procs,nb_voisin)) ; 
    voisin_proc = 0
  
    stId = nf90_close(ncid);  call test_cdf(stId)
    
  ! Maintenat que l'on connait le nombre de processus qui a tourne
   ! on lit les autres fichiers
   ny=0
   do iproc=1,nb_procs 
     write(filename,'(a4,a1,i3.3,a1,a)') trim(prefixb),'_',iproc-1,'_',trim(run_name)
     write(*,*) 'read tmp file ',filename
     stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
     	call test_cdf(stId)
     	call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
     	if ((coord_proc(iproc,ndims-1)).gt.ny) ny=coord_proc(iproc,ndims-1)
     stId = nf90_close(ncid);  call test_cdf(stId)
     write(*,*) 'close file ',filename
   enddo
   ny=ny+1
   nz=nb_procs/ny
  allocate(proc_tab(ny,nz))
  do iproc=1,nb_procs 
     write(filename,'(a4,a1,i3.3,a1,a)') trim(prefixb),'_',iproc-1,'_',trim(run_name)
     write(*,*) 'read tmp file ',filename
     stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
     	call test_cdf(stId)
     	call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
     proc_tab(coord_proc(iproc,ndims-1)+1,coord_proc(iproc,ndims)+1)=iproc
     stId = nf90_close(ncid);  call test_cdf(stId)
     write(*,*) 'close file ',filename
   enddo

  ! Maintenat que l'on connait le nombre de processus qui a tourne
   ! on lit les autres fichiers
        debz = 1
  do jj=1,nz
	deby = 1
  do ii=1,ny
    iproc=proc_tab(ii,jj)
    !--Create file name
    write(filename,'(a4,a1,i3.3,a1,a)') trim(prefixb),'_',iproc-1,'_',trim(run_name)
    write(*,*) 'read file ',filename
    stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
    	call test_cdf(stId)
    	call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
    	!--Get number point per processus
    	call get_simple_variable_cdf(ncid,"ncm",itmp(:))
    	call set_grid(itmp,2) 

    allocate(Ax_proc(ncm(1),ncm(2),ncm(3),n_spe))
    Ax_proc = zero 

 
    !--Get Total Number of Particles
    call get_simple_variable_cdf(ncid,"nptot",nptot_proc(iproc))
    call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )

      finy = deby + ncm(2) - 2
      finz = debz + ncm(3) - 2
!	do ish=1,n_spe
      ! we are just interested for O density
      ish = 5
	
       	write(nam,'(a,i2)')"Spe_",ish
!		call get_simple_variable_cdf(ncid,nam,nam_tmp)
        	!name_dens(ish)=nam_tmp
!        	print *,'nam_tmp,nam',nam_tmp,nam
                nam_tmp = "Den_O"
                print *,'nam_tmp,nam',nam_tmp,nam
		call get_simple_variable_cdf(ncid,nam_tmp,Ax_proc(:,:,:,ish))
		write(*,*) 'read file ',nam_tmp
		A(1:ncm_tot(1),deby:finy,debz:finz) =&
		  &	Ax_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,ish)
		
!	        Ax(1:ncm_tot(1),deby:finy,debz:finz,ish) = Ax_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,ish)
!	enddo

   deby=deby+(ncm(2)-2)
   deallocate(Ax_proc)
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename
   enddo  
   debz=debz+(ncm(3)-2)   
   enddo
   nptot = sum(nptot_proc)  
   deallocate(proc_tab)
   
   Qx = Qx*A
   
   
   sigma_Hp_O = 1.e-19 
      
   ! normalization
   Qx = Qx*(Spe%ref%density*1.e-6)*(Spe%ref%alfvenspeed*1.e5)*(sigma_Hp_O*1.e4)
   
   do ii = 1,ncm_tot(1)
     do jj = 1,ncm_tot(2)
       do kk = 1,ncm_tot(3)
         if (ieee_is_nan(Qx(ii,jj,kk))) Qx(ii,jj,kk)=0
       enddo
     enddo
   enddo  
         
   
   !--Creation du fichier de diagnostique de champ a acces sequentiel
     write_name = prefix_output//'_'//trim(run_name)//".nc"
     call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
   
     !--Open the file for write Particles
     stId = nf90_create(trim(write_name), nf90_clobber, ncid)
     call test_cdf(stId)
   
   call set_global_attribute_cdf(ncid,"Field")
   
   !--This is needed by m_wrt_common_cdf to use global number of points
     !and not relative to processus
     call create_wrt_dimensions_cdf(ncid,dimid,ncm_tot)
   
     ii = 1
     call common_def_var_cdf(ncid,varid,dimid,ii)
     
     !--Define Speces dimensions
     call species_def_dim_cdf(ncid,dimspec)
     
     !--Define Speces variables
     call species_def_var_cdf(ncid,varid,dimspec,ii)
   
     stId = nf90_def_var(ncid, "s_min",   QP_NF90_DP,dimid(3), varid(ii))
     call test_cdf(stId); ii = ii+1 
     stId = nf90_def_var(ncid, "s_max",   QP_NF90_DP,dimid(3), varid(ii))
     call test_cdf(stId); ii = ii+1   
   
	!do ish=1,n_spe
	nam_tmp="ProdCX_O+"
	    stId = nf90_def_var(ncid,nam_tmp, QP_NF90_DP,dimid(6:8),varid(ii))
	    call test_cdf(stId); ii = ii+1   
 	!enddo
	ii = ii-1   

     
     !--Switch to write mode
       stId = nf90_enddef(ncid); call test_cdf(stId)
     
       ii = 1
       !--Write common variables into the file
       call common_put_var_cdf(ncid,varid,ii)
       
       !--Write Species infos into the file
       call species_put_var_cdf(Spe,ncid,varid,dimspec,ii)
     
       stId = nf90_put_var(ncid, varid(ii), s_min)
       call test_cdf(stId); ii = ii+1              
       stId = nf90_put_var(ncid, varid(ii), s_max)
       call test_cdf(stId); ii = ii+1              

	!do ish=1,n_spe
	    stId = nf90_put_var(ncid, varid(ii), Qx(:,:,:))
	    call test_cdf(stId); ii = ii+1   
 	!enddo
	ii = ii-1   
       
       !--Close the file
         stId = nf90_close(ncid); call test_cdf(stId)
         
       
   
   !deallocate(Ax_proc)
   deallocate(Qx,A)
   deallocate(voisin_proc,dims,coord_proc,nptot_proc)
   
  
  end subroutine merge_field_cdf_pro
  
    
    

!!############################################################


!!############################################################
!!=============================================================
!!subroutine: diag/create_file_diag
!! FUNCTION 
!!  Create the name of the file containing fields information
!! INPUT
!!  filwrt=suffix containgt data
!!  me=processus identity
!! OUTPUT
!!  name_file=name of the file where fields will be recorded
 subroutine create_file_diag(name_file,filwrt,prefix,me)
     
 integer,intent(in) :: me
 character(len=40),intent(out) :: name_file 
 character(len=*),intent(in) :: filwrt,prefix 
  
  write(name_file,'(a4,a1,i3.3,a1,a)')trim(prefix),"_",me,'_',trim(filwrt)
  
 #ifdef HAVE_NETCDF
    name_file = trim(name_file)//".nc"
 #endif
    
   
 end subroutine create_file_diag
 
 !!!!!################################################################

 
 #endif
 end module CX_sputt
