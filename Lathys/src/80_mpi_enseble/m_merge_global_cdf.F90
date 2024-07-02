!=============================================================
!=============================================================
module m_merge_global_cdf

 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype
 use defs_particletype
 use defs_species
 use defs_parametre,only : fildat,dt,ns,npm,gstep
 use m_writeout
 use m_verify_reassemble_mpi
#ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
#endif
#include "q-p_common.h"

 implicit none
 private 

#ifdef HAVE_NETCDF
 private ::               &
      merge_field_cdf,    &
      merge_field_cdf_atm,&
      merge_field_cdf_pro,&
      create_file_diag,   &
      merge_moment_species_cdf,&
      merge_field_cdf_1_array, &
      merge_field_cdf_Te_array

 public ::                &
      merge_cdf
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
   !!subroutine: m_split_reassemble/Merge_cdf
   !! This routine reads one splited file and merge it 
   !! accoriding to processes position
   
    subroutine merge_cdf(run_name)
   ! this routine merge the file created during the post_treatment to 4 files
   ! 1 file mag : including only magnetic arrays
   ! 1 file elec : including only electric arrays
   ! 1 file vel  : including only velocity (bulk) arrays
   ! 1 file dens : including only density (electron number) array
   ! after the creation for each process it recombines in one unique file
   character(len=*),intent(in) :: run_name
   !Local Variable
   character(len=64) :: var_name1,var_name2,var_name3,var_name4,var_name5
   character(len=4)  :: prefix,prefix_output,prefix1,prefix2,prefix3
   
   prefix = "Mag3"
   var_name1 = "Bx"
   var_name2 = "By"
   var_name3 = "Bz"
   var_name4 = "Btot"  
   prefix_output = "Magw"
   
   call merge_field_cdf(run_name,prefix,var_name1,var_name2,var_name3,var_name4,prefix_output) 
   
  !   prefix = "Vel3"
  !   var_name1 = "Ux"
  !   var_name2 = "Uy"
  !   var_name3 = "Uz"
     
  ! prefix_output = "Velw"
  ! call merge_field_cdf(run_name,prefix,var_name1,var_name2,var_name3,prefix_output) 
   
     prefix = "Ele3"
     var_name1 = "Ex"
     var_name2 = "Ey"
     var_name3 = "Ez"
     var_name4 = "Etot"
   prefix_output = "Elew"
   call merge_field_cdf(run_name,prefix,var_name1,var_name2,var_name3,var_name4,prefix_output)    
   
  !   prefix = "Den3"
  !   var_name1 = "Density"
     
  ! prefix_output = "Denw"
  ! call merge_field_cdf_1_array(run_name,prefix,var_name1,prefix_output)    
   
  !   prefix = "Pre3"
  !   var_name1 = "Pe_tot"
        
  ! prefix_output = "Prew"
  ! call merge_field_cdf_1_array(run_name,prefix,var_name1,prefix_output) 
   
      
   ! create file of electronic temperature
  !   prefix1 = "Den3"
  !   prefix2 = "Pre3"
  !   var_name1 = "Density"
  !   var_name2 = "Pe_tot"
        
  !   prefix_output = "Tew_"
  ! call merge_field_cdf_Te_array(run_name,prefix1,var_name1,prefix2,var_name2,prefix_output)   
    
     prefix = "Atm3"
     prefix_output = "Atmw"
   call merge_field_cdf_atm(run_name,prefix,prefix_output)    
   
   prefix = "Mom3"
   prefix_output = "Momw"
   call merge_moment_species_cdf(run_name,prefix,prefix_output)
 
   prefix = "Pro3"
   prefix_output = "Prow"
  ! call merge_field_cdf_pro(run_name,prefix,prefix_output) 
   
   
   prefix1 = "Den3"
   prefix2 = "Vel3"
   prefix3 = "Pre3"
   var_name1 = "Density"
   var_name2 = "Ux"
   var_name3 = "Uy"
   var_name4 = "Uz"
   var_name5 = "Pe_tot"
   prefix_output = "Thew"
   call merge_moment_thermal_cdf(run_name,prefix1,prefix2,prefix3,var_name1,var_name2,var_name3,var_name4,var_name5,prefix_output)
     
   
   
  end subroutine merge_cdf
 
 !!***********************************************
 !!-----------------------------------------------
 !! m_split_reassemble_cdf.F90
 !! reads 1 file computed for each process
 !! and merge it to one big file
 
  subroutine merge_field_cdf(run_name,prefix,varname1,varname2,varname3,varname4,prefix_output)
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: prefix,prefix_output
    character(len=*),intent(in) :: varname1,varname2,varname3,varname4
    !--Reading variables
    integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
    integer :: rang,nb_procs,ndims,nb_voisin
    integer,allocatable :: dims(:),nptot_proc(:),dimspec(:)
    integer,allocatable :: voisin_proc(:,:)
    integer,allocatable :: coord_proc(:,:),proc_tab(:,:)
    integer :: sgn,dimchar
    character(len=10) :: unit,unit_axis,coord
    real(dp) :: nrm
    
      
    !--Variables to read fields diagnostic
    real(dp),dimension(:,:,:),allocatable :: Ax_proc,Ay_proc,Az_proc
    real(dp),dimension(:,:,:),allocatable :: Ax,Ay,Az,Atot
    real(dp),dimension(:),allocatable :: X_axis,Y_axis,Z_axis
    
    !--Others
    integer :: ind,deby,debz,finy,finz,iproc
    integer :: deb,fin,cumul
    integer  :: stId,ncid,ii,is,ny,nz,jj,i
    integer  :: varid(1000),dimid(12)
    integer  :: itmp(3),var_id
    real(dp) :: rtmp(3)
    logical :: file_e
    character(len=8) :: planet
    character(len=50) :: write_name
    character(len=64) :: filename
     character(len=100) :: test_name
    character(len=500) :: msg
    real(dp) :: n0,B0,V0,x0
   
    !--Create file name for Field (READ, proc 0) 
    write(filename,'(a4,a1,i4.4,a1,2a)')trim(prefix),"_",0,'_',trim(run_name),".nc"
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
       write(*,*) 'open file ',filename
    !--Control the Date
    stId = nf90_get_att(ncid, nf90_global, "Date", test_name)
    call test_cdf(stId)
    
    !--Get Planetname
    stId = nf90_inq_varid(ncid,"planetname",var_id)
    stId = nf90_get_var(ncid,var_id,planet)
   write(*,*) 'Planet :',planet

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
       
    allocate(voisin_proc(nb_procs,nb_voisin)) ; voisin_proc = 0

   ! get normalization value
   !--Get normalisation values
   call get_simple_variable_cdf(ncid,"phys_mag",B0)
   call get_simple_variable_cdf(ncid,"phys_density",n0)
   call get_simple_variable_cdf(ncid,"phys_speed",V0)  
   call get_simple_variable_cdf(ncid,"phys_length",x0)  
  
    stId = nf90_close(ncid);  call test_cdf(stId)
         write(*,*) 'close file ',filename

     ! coordinate transformation for each object
     select case (trim(planet))
     case ("mars","mercury")
       sgn = -1
       coord="MSO"
     case ("ganymede")
       sgn=+1
       coord = "GPhiO"
     case("titan")
       sgn = +1
       coord="TIIS"
     case("earth")
       sgn=-1
       coord="GSM"
     case default 
       planet="mars    "
       sgn = -1
       coord = "MSO"
    end select
  
  
   select case(trim(prefix_output))
   case("Magw")
      nrm = B0*1.e9
     unit ="nT"
   case("Elew")
     nrm = V0*B0
     unit = "mV.m-1"
   case("Velw")
     nrm = V0
     unit="km.s-1"
   case default
     nrm=1
     unit =  ""
   end select  
  
  
  !--Allocation of global arrays
   allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   allocate(Ay(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   allocate(Az(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   allocate(Atot(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   Ax = zero    ; Ay = zero     ; Az = zero    ; Atot = zero
   
   !-- Allocate coordinate axis
   allocate(X_axis(ncm_tot(1)));        allocate(Y_axis(ncm_tot(2)))
   allocate(Z_axis(ncm_tot(3)))
   X_axis = zero;       Y_axis = zero;  Z_axis = zero
   
   X_axis = sgn*((/(i,i=1,ncm_tot(1))/)*gstep(1)-Spe%P%centr(1))*x0
   Y_axis = sgn*((/(i,i=1,ncm_tot(2))/)*gstep(2)-Spe%P%centr(2))*x0
   Z_axis = ((/(i,i=1,ncm_tot(3))/)*gstep(3)-Spe%P%centr(3))*x0
   unit_axis = "km"
   
   if (trim(prefix) == 'Ele3') then
     X_axis = X_axis - sgn*x0*gstep(1)/2._dp
     Y_axis = Y_axis - sgn*x0*gstep(2)/2._dp
     Z_axis = Z_axis - x0*gstep(3)/2._dp
   endif  


   ny=0
   do iproc=1,nb_procs 
     write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
    ! write(*,*) 'read tmp file ',filename
     stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
     call test_cdf(stId)
     call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
     if ((coord_proc(iproc,ndims-1)).gt.ny) ny=coord_proc(iproc,ndims-1)
     stId = nf90_close(ncid);  call test_cdf(stId)
     !write(*,*) 'close file ',filename
   enddo
   ny=ny+1
   nz=nb_procs/ny
  allocate(proc_tab(ny,nz))
  do iproc=1,nb_procs 
     write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
    ! write(*,*) 'read tmp file ',filename
     stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
     call test_cdf(stId)
     call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
     proc_tab(coord_proc(iproc,ndims-1)+1,coord_proc(iproc,ndims)+1)=iproc
     stId = nf90_close(ncid);  call test_cdf(stId)
    ! write(*,*) 'close file ',filename
   enddo

  ! Maintenat que l'on connait le nombre de processus qui a tourne
   ! on lit les autres fichiers
        debz = 1
  do jj=1,nz
        deby = 1
  do ii=1,ny
    iproc=proc_tab(ii,jj)
    !--Create file name
    write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
   ! write(*,*) 'read file ',filename
    stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
    call test_cdf(stId)
    call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
    !--Get number point per processus
    call get_simple_variable_cdf(ncid,"ncm",itmp(:))
    call set_grid(itmp,2)

    !--allocation des tableaux de lectures du fichier diag de champ
    allocate(Ax_proc(ncm(1),ncm(2),ncm(3)))
    allocate(Ay_proc(ncm(1),ncm(2),ncm(3)))
    allocate(Az_proc(ncm(1),ncm(2),ncm(3)))
    Ax_proc = zero ; Ay_proc = zero ; Az_proc = zero

    !--Get Total Number of Particles
    call get_simple_variable_cdf(ncid,"nptot",nptot_proc(iproc))
    call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )
    call get_simple_variable_cdf(ncid,varname1 ,Ax_proc(:,:,:) )
    call get_simple_variable_cdf(ncid,varname2 ,Ay_proc(:,:,:) )
    call get_simple_variable_cdf(ncid,varname3 ,Az_proc(:,:,:) )
    
   if (trim(prefix) /= 'Ele3') then     
   !--On reconstruit les tableaux ayant la meme grille , celle de B
      finy = deby + ncm(2) - 2
      finz = debz + ncm(3) - 2
      if ((finy.gt.ncm_tot(2)).or.(finz.gt.ncm_tot(3))) then
        print *,'finy ',finy,'finz ',finz,'ii ',ii,'jj ',jj
        stop
      endif
      Ax(1:ncm_tot(1),deby:finy,debz:finz) = Ax_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
      Ay(1:ncm_tot(1),deby:finy,debz:finz) = Ay_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
      Az(1:ncm_tot(1),deby:finy,debz:finz) = Az_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
   else
   !--on reconstruit les tableaux du champ lectrique
      finy = deby + ncm(2) - 1
      finz = debz + ncm(3) - 1
      if ((finy.gt.ncm_tot(2)).or.(finz.gt.ncm_tot(3))) then
        print *,'finy ',finy,'finz ',finz,'ii ',ii,'jj ',jj
        stop
      endif
      Ax(1:ncm_tot(1),deby:finy,debz:finz) = Ax_proc(1:ncm_tot(1),1:ncm(2),1:ncm(3))
      Ay(1:ncm_tot(1),deby:finy,debz:finz) = Ay_proc(1:ncm_tot(1),1:ncm(2),1:ncm(3))
      Az(1:ncm_tot(1),deby:finy,debz:finz) = Az_proc(1:ncm_tot(1),1:ncm(2),1:ncm(3))
   endif
   deby=deby+(ncm(2)-2)
   deallocate(Ax_proc,Ay_proc,Az_proc)
    stId = nf90_close(ncid);  call test_cdf(stId)
   ! write(*,*) 'close file ',filename
   enddo  
   debz=debz+(ncm(3)-2)   
   enddo
   nptot = sum(nptot_proc)  
   deallocate(proc_tab)
   
   ! normalization
   Ax = Ax*nrm*sgn
   Ay = Ay*nrm*sgn
   Az = Az*nrm
   Atot = sqrt(Ax**2+Ay**2+Az**2)
   
   !--Creation du fichier de diagnostique de champ a acces sequentiel
     write_name = prefix_output//'_'//trim(run_name)//".nc"
     call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
   
     !--Open the file for write Particles
     stId = nf90_create(trim(write_name), nf90_64BIT_OFFSET, ncid)
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
     stId = nf90_def_dim(ncid, "dimunit",10, dimchar)     
     call test_cdf(stId)
   
     stId = nf90_def_var(ncid, "s_min",   QP_NF90_DP,dimid(3), varid(ii))
     call test_cdf(stId); ii = ii+1 
     stId = nf90_def_var(ncid, "s_max",   QP_NF90_DP,dimid(3), varid(ii))
     call test_cdf(stId); ii = ii+1   
     stId = nf90_def_var(ncid, "unit",   nf90_char,dimchar, varid(ii))
     call test_cdf(stId); ii = ii+1  
     stId = nf90_def_var(ncid, "X_axis",   QP_NF90_DP,dimid(6), varid(ii))
     call test_cdf(stId); ii = ii+1     
     stId = nf90_def_var(ncid, "Y_axis",   QP_NF90_DP,dimid(7), varid(ii))
     call test_cdf(stId); ii = ii+1     
     stId = nf90_def_var(ncid, "Z_axis",   QP_NF90_DP,dimid(8), varid(ii))
     call test_cdf(stId); ii = ii+1     
     stId = nf90_def_var(ncid, "Unit_axis",   nf90_char,dimchar, varid(ii))
     call test_cdf(stId); ii = ii+1     
     stId = nf90_def_var(ncid, "Coordinate_system",  nf90_char,dimchar, varid(ii))
     call test_cdf(stId); ii = ii+1      

     
   
     stId = nf90_def_var(ncid, varname1, QP_NF90_DP,dimid(6:8),varid(ii))
     call test_cdf(stId); ii = ii+1                                                
     stId = nf90_def_var(ncid, varname2, QP_NF90_DP,dimid(6:8),varid(ii))
     call test_cdf(stId); ii = ii+1                                                
     stId = nf90_def_var(ncid, varname3, QP_NF90_DP,dimid(6:8),varid(ii))
     call test_cdf(stId); ii = ii+1
     !stId = nf90_def_var(ncid, varname4, QP_NF90_DP,dimid(6:8),varid(ii))
     !call test_cdf(stId)
     
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
       stId = nf90_put_var(ncid, varid(ii), unit)
       call test_cdf(stId); ii = ii+1    
       stId = nf90_put_var(ncid, varid(ii), X_axis)
       call test_cdf(stId); ii = ii+1        
       stId = nf90_put_var(ncid, varid(ii), Y_axis)
       call test_cdf(stId); ii = ii+1        
       stId = nf90_put_var(ncid, varid(ii), Z_axis)
       call test_cdf(stId); ii = ii+1        
       stId = nf90_put_var(ncid, varid(ii), Unit_axis)
       call test_cdf(stId); ii = ii+1       
       stId = nf90_put_var(ncid, varid(ii), coord)
       call test_cdf(stId); ii = ii+1        
       
       stId = nf90_put_var(ncid, varid(ii), Ax)
       call test_cdf(stId); ii = ii+1              
       stId = nf90_put_var(ncid, varid(ii), Ay)
       call test_cdf(stId); ii = ii+1 
       stId = nf90_put_var(ncid, varid(ii), Az)
       call test_cdf(stId);ii = ii+1
       !stId = nf90_put_var(ncid, varid(ii), Atot)
       !call test_cdf(stId)
       
       !--Close the file
         stId = nf90_close(ncid); call test_cdf(stId)
       
   
   deallocate(Ax,Ay,Az)
   deallocate(voisin_proc,dims,coord_proc,nptot_proc)
   deallocate(X_axis,Y_axis,Z_axis)
   !call dealloc_arr3D(Bfield)
  
  end subroutine merge_field_cdf
  
!!***********************************************
 !!-----------------------------------------------
 !! m_split_reassemble_cdf.F90
 !! reads 1 file computed for each process
 !! and merge it to one big file
 
  subroutine merge_field_cdf_1_array(run_name,prefix,varname1,prefix_output)
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: prefix,prefix_output
    character(len=*),intent(in) :: varname1
    !--Reading variables
    integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
    integer :: rang,nb_procs,ndims,nb_voisin
    integer,allocatable :: dims(:),nptot_proc(:),dimspec(:)
    integer,allocatable :: voisin_proc(:,:)
    integer,allocatable :: coord_proc(:,:),proc_tab(:,:)
    
    
      
    !--Variables to read fields diagnostic
    real(dp),dimension(:,:,:),allocatable :: Ax_proc
    real(dp),dimension(:,:,:),allocatable :: Ax
    
    !--Others
    integer :: ind,deby,debz,finy,finz,iproc
    integer :: deb,fin,cumul
    integer  :: stId,ncid,ii,is,jj,ny,nz
    integer  :: varid(1000),dimid(12)
    integer  :: itmp(3)
    real(dp) :: rtmp(3)
    logical :: file_e
    character(len=50) :: write_name
    character(len=64) :: filename
     character(len=100) :: test_name
    character(len=500) :: msg
    real(dp) :: nrm,n0
    character(len=10) :: unit
    integer :: dimchar
    
   
    !--Create file name for Field (READ, proc 0) 
    write(filename,'(a4,a1,i4.4,a1,2a)')trim(prefix),"_",0,'_',trim(run_name),".nc"
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
  
  
  
  
    !--allocation des tableaux de lectures du fichier diag de champ
   allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   Ax = zero
     
    allocate(voisin_proc(nb_procs,nb_voisin)) ; voisin_proc = 0
  
     !--Get normalisation values
    
     call get_simple_variable_cdf(ncid,"phys_density",n0)
     
     select case(trim(prefix_output))
     case("Denw")
        nrm = n0*1.e-6
       unit ="cm-3"
     case default
       nrm=1
       unit =  ""
   end select  
  
    stId = nf90_close(ncid);  call test_cdf(stId)
    
   ny=0
   do iproc=1,nb_procs 
     write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
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
     write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
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
    write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
    write(*,*) 'read file ',filename
    stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
    call test_cdf(stId)
    call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
    !--Get number point per processus
    call get_simple_variable_cdf(ncid,"ncm",itmp(:))
    call set_grid(itmp,2)

    !--allocation des tableaux de lectures du fichier diag de champ
    allocate(Ax_proc(ncm(1),ncm(2),ncm(3)))
    Ax_proc = zero ; 

    !--Get Total Number of Particles
    call get_simple_variable_cdf(ncid,"nptot",nptot_proc(iproc))
    call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )
    call get_simple_variable_cdf(ncid,varname1 ,Ax_proc(:,:,:) )
    
   if (trim(prefix) /= 'Ele3') then     
   !--On reconstruit les tableaux ayant la meme grille , celle de B
      finy = deby + ncm(2) - 2
      finz = debz + ncm(3) - 2
      if ((finy.gt.ncm_tot(2)).or.(finz.gt.ncm_tot(3))) then
        print *,'finy ',finy,'finz ',finz,'ii ',ii,'jj ',jj
        stop
      endif
      Ax(1:ncm_tot(1),deby:finy,debz:finz) = Ax_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
   else
   !--on reconstruit les tableaux du champ lectrique
      finy = deby + ncm(2) - 1
      finz = debz + ncm(3) - 1
      if ((finy.gt.ncm_tot(2)).or.(finz.gt.ncm_tot(3))) then
        print *,'finy ',finy,'finz ',finz,'ii ',ii,'jj ',jj
        stop
      endif
      Ax(1:ncm_tot(1),deby:finy,debz:finz) = Ax_proc(1:ncm_tot(1),1:ncm(2),1:ncm(3))
   endif
   deby=deby+(ncm(2)-2)
   deallocate(Ax_proc)
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename
   enddo  
   debz=debz+(ncm(3)-2)   
   enddo
   nptot = sum(nptot_proc)  
   deallocate(proc_tab)
   
 
   nptot = sum(nptot_proc)  
   
   
   ! normalization
   Ax = nrm*Ax

   
   !--Creation du fichier de diagnostique de champ a acces sequentiel
     write_name = prefix_output//'_'//trim(run_name)//".nc"
     call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
   
     !--Open the file for write Particles
     stId = nf90_create(trim(write_name), nf90_64BIT_OFFSET, ncid)
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
     stId = nf90_def_dim(ncid, "dimunit",10, dimchar)     
     call test_cdf(stId)
   
     stId = nf90_def_var(ncid, "s_min",   QP_NF90_DP,dimid(3), varid(ii))
     call test_cdf(stId); ii = ii+1 
     stId = nf90_def_var(ncid, "s_max",   QP_NF90_DP,dimid(3), varid(ii))
     call test_cdf(stId); ii = ii+1   
     stId = nf90_def_var(ncid, "unit",   nf90_char,dimchar, varid(ii))
     call test_cdf(stId); ii = ii+1   
   
   
     stId = nf90_def_var(ncid, varname1, QP_NF90_DP,dimid(6:8),varid(ii))
     call test_cdf(stId)
     
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
       stId = nf90_put_var(ncid, varid(ii), unit)
       call test_cdf(stId); ii = ii+1              
     
       stId = nf90_put_var(ncid, varid(ii), Ax)
       call test_cdf(stId)
       
       !--Close the file
         stId = nf90_close(ncid); call test_cdf(stId)
       
   

   deallocate(Ax)
   deallocate(voisin_proc,dims,coord_proc,nptot_proc)
   
  
  end subroutine merge_field_cdf_1_array

!!***********************************************
 !!-----------------------------------------------
 !! m_merge_global_cdf.F90
 !! reads 1 file computed for each process
 !! and merge it to one big file
 
  subroutine merge_field_cdf_Te_array(run_name,prefix1,varname1,prefix2,varname2,prefix_output)
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: prefix1,prefix_output,prefix2
    character(len=*),intent(in) :: varname1,varname2
    character(len=2) :: varname_output
    !--Reading variables
    integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
    integer :: rang,nb_procs,ndims,nb_voisin
    integer,allocatable :: dims(:),nptot_proc(:),dimspec(:)
    integer,allocatable :: voisin_proc(:,:)
    integer,allocatable :: coord_proc(:,:),proc_tab(:,:)
    
    
      
    !--Variables to read fields diagnostic
    real(dp),dimension(:,:,:),allocatable :: A_proc,B_proc
    real(dp),dimension(:,:,:),allocatable :: A,B
    
    !--Others
    integer :: ind,deby,debz,finy,finz,iproc
    integer :: deb,fin,cumul
    integer  :: stId,ncid,ii,is,jj,ny,nz
    integer  :: varid(1000),dimid(12)
    integer  :: itmp(3)
    real(dp) :: rtmp(3)
    logical :: file_e
    character(len=50) :: write_name
    character(len=64) :: filename1,filename2
     character(len=100) :: test_name
    character(len=500) :: msg
   
    !--Create file name for Field (READ, proc 0) 
    write(filename1,'(a4,a1,i4.4,a1,2a)')trim(prefix1),"_",0,'_',trim(run_name),".nc"
    write(filename2,'(a4,a1,i4.4,a1,2a)')trim(prefix2),"_",0,'_',trim(run_name),".nc"
    write(*,*) 'merging file'  ,filename1,filename2  
  
  !---- Reading density file
  
    !--Inquire if the file exists
    inquire( file=trim(filename1), exist=file_e )
    !--No file: return
    if(.not.(file_e)) then
     call wrtout(6,"File: "//trim(filename1)//" does not exits","PERS")
     return
    endif
  
    !--Open NetCDF file
    stId = nf90_open(trim(filename1), nf90_nowrite, ncid)
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
  
  
  
  
  
    !--allocation des tableaux de lectures du fichier diag de champ
   allocate(A(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   A = zero
   allocate(B(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   B = zero
     
     
    allocate(voisin_proc(nb_procs,nb_voisin)) ; voisin_proc = 0
  
    stId = nf90_close(ncid);  call test_cdf(stId)
    
   ny=0
   do iproc=1,nb_procs 
     write(filename1,'(a4,a1,i4.4,a1,a)') trim(prefix1),'_',iproc-1,'_',trim(run_name)
     write(*,*) 'read tmp file ',filename1
     stId = nf90_open(trim(filename1)//".nc", nf90_nowrite, ncid)
     call test_cdf(stId)
     call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
     if ((coord_proc(iproc,ndims-1)).gt.ny) ny=coord_proc(iproc,ndims-1)
     stId = nf90_close(ncid);  call test_cdf(stId)
     write(*,*) 'close file ',filename1
   enddo
   ny=ny+1
   nz=nb_procs/ny
  allocate(proc_tab(ny,nz))
  do iproc=1,nb_procs 
     write(filename1,'(a4,a1,i4.4,a1,a)') trim(prefix1),'_',iproc-1,'_',trim(run_name)
     write(*,*) 'read tmp file ',filename1
     stId = nf90_open(trim(filename1)//".nc", nf90_nowrite, ncid)
     call test_cdf(stId)
     call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
     proc_tab(coord_proc(iproc,ndims-1)+1,coord_proc(iproc,ndims)+1)=iproc
     stId = nf90_close(ncid);  call test_cdf(stId)
     write(*,*) 'close file ',filename1
   enddo

  ! Maintenat que l'on connait le nombre de processus qui a tourne
   ! on lit les autres fichiers
        debz = 1
  do jj=1,nz
        deby = 1
  do ii=1,ny
    iproc=proc_tab(ii,jj)
    !--Create file name
    write(filename1,'(a4,a1,i4.4,a1,a)') trim(prefix1),'_',iproc-1,'_',trim(run_name)
    write(*,*) 'read file ',filename1
    stId = nf90_open(trim(filename1)//".nc", nf90_nowrite, ncid)
    call test_cdf(stId)
    call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
    !--Get number point per processus
    call get_simple_variable_cdf(ncid,"ncm",itmp(:))
    call set_grid(itmp,2)

    !--allocation des tableaux de lectures du fichier diag de champ
    allocate(A_proc(ncm(1),ncm(2),ncm(3)))
    A_proc = zero ; 

    !--Get Total Number of Particles
    call get_simple_variable_cdf(ncid,"nptot",nptot_proc(iproc))
    call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )
    call get_simple_variable_cdf(ncid,varname1 ,A_proc(:,:,:) )
    
     
   !--On reconstruit les tableaux ayant la meme grille , celle de B
      finy = deby + ncm(2) - 2
      finz = debz + ncm(3) - 2
      if ((finy.gt.ncm_tot(2)).or.(finz.gt.ncm_tot(3))) then
        print *,'finy ',finy,'finz ',finz,'ii ',ii,'jj ',jj
        stop
      endif
      A(1:ncm_tot(1),deby:finy,debz:finz) = A_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
  
   deby=deby+(ncm(2)-2)
   deallocate(A_proc)
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename1
   enddo  
   debz=debz+(ncm(3)-2)   
   enddo
   nptot = sum(nptot_proc)  
   !deallocate(proc_tab)
  
   
   
   !--- Reading Pe file
        debz = 1
  do jj=1,nz
        deby = 1
  do ii=1,ny
    iproc=proc_tab(ii,jj)
    !--Create file name
    write(filename2,'(a4,a1,i4.4,a1,a)') trim(prefix2),'_',iproc-1,'_',trim(run_name)
    write(*,*) 'read file ',filename2
    stId = nf90_open(trim(filename2)//".nc", nf90_nowrite, ncid)
    call test_cdf(stId)
    call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
    !--Get number point per processus
    call get_simple_variable_cdf(ncid,"ncm",itmp(:))
    call set_grid(itmp,2)

    !--allocation des tableaux de lectures du fichier diag de champ
    allocate(B_proc(ncm(1),ncm(2),ncm(3)))
    B_proc = zero ; 

    !--Get Total Number of Particles
    call get_simple_variable_cdf(ncid,"nptot",nptot_proc(iproc))
    call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )
    call get_simple_variable_cdf(ncid,varname2 ,B_proc(:,:,:) )
    
     
   !--On reconstruit les tableaux ayant la meme grille , celle de B
      finy = deby + ncm(2) - 2
      finz = debz + ncm(3) - 2
      if ((finy.gt.ncm_tot(2)).or.(finz.gt.ncm_tot(3))) then
        print *,'finy ',finy,'finz ',finz,'ii ',ii,'jj ',jj
        stop
      endif
      B(1:ncm_tot(1),deby:finy,debz:finz) = B_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
  
   deby=deby+(ncm(2)-2)
   deallocate(B_proc)
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename2
   enddo  
   debz=debz+(ncm(3)-2)   
   enddo   
   
   
   deallocate(proc_tab)
   
   ! creating the Te array
   A = B/A
   varname_output = 'Temperature'
   
   !--Creation du fichier de diagnostique de champ a acces sequentiel
     write_name = prefix_output//'_'//trim(run_name)//".nc"
     call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
   
     !--Open the file for write Particles
     stId = nf90_create(trim(write_name), nf90_64BIT_OFFSET, ncid)
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
   
     stId = nf90_def_var(ncid, varname_output, QP_NF90_DP,dimid(6:8),varid(ii))
     call test_cdf(stId)
     
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
       stId = nf90_put_var(ncid, varid(ii), A)
       call test_cdf(stId)
       
       !--Close the file
         stId = nf90_close(ncid); call test_cdf(stId)
       
   

   deallocate(A,B)
   deallocate(voisin_proc,dims,coord_proc,nptot_proc)
   
  
  end subroutine merge_field_cdf_Te_array
  
  
  
  !!***********************************************
   !!-----------------------------------------------
   !! m_split_reassemble_cdf.F90
   !! reads 3 files (density, speed and pressure) computed for each process
   !! and merge it to one big file
   
    subroutine merge_moment_thermal_cdf(run_name,prefix1,prefix2,prefix3,varname1,varname2,varname3,&
    & varname4,varname5,prefix_output)
      character(len=*),intent(in) :: run_name
      character(len=*),intent(in) :: prefix1,prefix2,prefix3,prefix_output
      character(len=*),intent(in) :: varname1,varname2,varname3,varname4,varname5
      !--Reading variables
      integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
      integer :: rang,nb_procs,ndims,nb_voisin
      integer,allocatable :: dims(:),nptot_proc(:),dimspec(:)
      integer,allocatable :: voisin_proc(:,:)
      integer,allocatable :: coord_proc(:,:),proc_tab(:,:)
      
      
        
      !--Variables to read fields diagnostic
      real(dp),dimension(:,:,:),allocatable :: A_proc,B_proc,Ax_proc,Ay_proc,Az_proc
      real(dp),dimension(:,:,:),allocatable :: A,B,Ax,Ay,Az
      real(dp),dimension(:),allocatable :: X_axis,Y_axis,Z_axis
      
      !--Others
      integer :: ind,deby,debz,finy,finz,iproc
      integer :: deb,fin,cumul
      integer  :: stId,ncid,ii,is,jj,ny,nz,i
      integer  :: varid(1000),dimid(12)
      integer  :: itmp(3)
      real(dp) :: rtmp(3)
      logical :: file_e
      character(len=50) :: write_name
      character(len=64) :: filename
       character(len=100) :: test_name
      character(len=500) :: msg
      real(dp) :: nrm_n,nrm_U,nrm_T,n0,V0,x0
      character(len=10) :: unit_n,unit_U,unit_T,coord,unit_axis
      integer :: dimchar,sgn,var_id
      character(len=8) :: planet
      
     
      !--Create file name for Field (READ, proc 0) 
      !-- Reading electron density file
      write(filename,'(a4,a1,i4.4,a1,2a)')trim(prefix1),"_",0,'_',trim(run_name),".nc"
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

      !--Get Planetname
      stId = nf90_inq_varid(ncid,"planetname",var_id)
      stId = nf90_get_var(ncid,var_id,planet)
      write(*,*) 'Planet :',planet        
    
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
    
    
    
    
      !--allocation des tableaux de lectures du fichier diag de champ
     allocate(A(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
     A = zero
     allocate(B(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
     B = zero
     allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
     Ax = zero  
     allocate(Ay(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
     Ay = zero  
     allocate(Az(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
     Az = zero
         ! coordinate transformation for each object
        select case (trim(planet))
        case ("mars","mercury")
          sgn = -1
          coord="MSO"
        case ("ganymede")
          sgn=+1
          coord = "GPHiO"
        case("titan")
          sgn = +1
          coord="TIIS"
        case("earth")
          sgn=-1
          coord="GSM"
        case default 
          planet="mars    "
          sgn = -1
          coord="MSO"
       end select    
     
       
      allocate(voisin_proc(nb_procs,nb_voisin)) ; voisin_proc = 0
    
       !--Get normalisation values
      
       call get_simple_variable_cdf(ncid,"phys_density",n0)
       call get_simple_variable_cdf(ncid,"phys_speed",V0)
       call get_simple_variable_cdf(ncid,"phys_length",x0)
       

         nrm_n = n0*1.e-6
         nrm_U = V0
         !nrm_T = (Spe%ref%mag**2)/(Spe%ref%density*2._dp*4._dp*pi*mu0*e_Cb)
         nrm_T = 2.*sum(Spe%S(:)%qms*Spe%S(:)%rmds)*(Spe%ref%mag**2)/(Spe%ref%density*2._dp*mu0*e_Cb)
         unit_n ="cm-3"
         unit_U = "km.s-1"
         unit_T = "eV"

   !-- Allocate coordinate axis
   allocate(X_axis(ncm_tot(1)));        allocate(Y_axis(ncm_tot(2)))
   allocate(Z_axis(ncm_tot(3)))
   X_axis = zero;       Y_axis = zero;  Z_axis = zero
   
   X_axis = sgn*((/(i,i=1,ncm_tot(1))/)*gstep(1)-Spe%P%centr(1))*x0
   Y_axis = sgn*((/(i,i=1,ncm_tot(2))/)*gstep(2)-Spe%P%centr(2))*x0
   Z_axis = ((/(i,i=1,ncm_tot(3))/)*gstep(3)-Spe%P%centr(3))*x0
   unit_axis = "km"
   
    
      stId = nf90_close(ncid);  call test_cdf(stId)
      
     ny=0
     do iproc=1,nb_procs 
       write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix1),'_',iproc-1,'_',trim(run_name)
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
       write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix1),'_',iproc-1,'_',trim(run_name)
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
      write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix1),'_',iproc-1,'_',trim(run_name)
      write(*,*) 'read file ',filename
      stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
      call test_cdf(stId)
      call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
      !--Get number point per processus
      call get_simple_variable_cdf(ncid,"ncm",itmp(:))
      call set_grid(itmp,2)
  
      !--allocation des tableaux de lectures du fichier diag de champ
      allocate(A_proc(ncm(1),ncm(2),ncm(3)))
      A_proc = zero ; 
  
      call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )
      call get_simple_variable_cdf(ncid,varname1 ,A_proc(:,:,:) )
      
     
     !--On reconstruit les tableaux ayant la meme grille , celle de B
        finy = deby + ncm(2) - 2
        finz = debz + ncm(3) - 2
        if ((finy.gt.ncm_tot(2)).or.(finz.gt.ncm_tot(3))) then
        print *,'finy ',finy,'finz ',finz,'ii ',ii,'jj ',jj
        stop
        endif
        A(1:ncm_tot(1),deby:finy,debz:finz) = A_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
     deby=deby+(ncm(2)-2)
     deallocate(A_proc)
      stId = nf90_close(ncid);  call test_cdf(stId)
      write(*,*) 'close file ',filename
     enddo  
     debz=debz+(ncm(3)-2)   
     enddo  
 !    deallocate(proc_tab)
     
   
     
     ! normalization
!!     A = nrm*A
     
     
     ! Now we read bulk speed file
      !--Create file name for Field (READ, proc 0) 
      !-- Reading electron density file
      write(filename,'(a4,a1,i4.4,a1,2a)')trim(prefix2),"_",0,'_',trim(run_name),".nc"
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
      
          debz = 1
    do jj=1,nz
        deby = 1
    do ii=1,ny
      iproc=proc_tab(ii,jj)
      !--Create file name
      write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix2),'_',iproc-1,'_',trim(run_name)
      write(*,*) 'read file ',filename
      stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
      call test_cdf(stId)
      call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
      !--Get number point per processus
      call get_simple_variable_cdf(ncid,"ncm",itmp(:))
      call set_grid(itmp,2)
  
      !--allocation des tableaux de lectures du fichier diag de champ
      allocate(Ax_proc(ncm(1),ncm(2),ncm(3)))
      Ax_proc = zero ; 
     allocate(Ay_proc(ncm(1),ncm(2),ncm(3)))
      Ay_proc = zero ;
      allocate(Az_proc(ncm(1),ncm(2),ncm(3)))
      Az_proc = zero ;     
  
      !--Get Total Number of Particles
      call get_simple_variable_cdf(ncid,"nptot",nptot_proc(iproc))
      call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )
      call get_simple_variable_cdf(ncid,varname2 ,Ax_proc(:,:,:) )
      call get_simple_variable_cdf(ncid,varname3 ,Ay_proc(:,:,:) )
      call get_simple_variable_cdf(ncid,varname4 ,Az_proc(:,:,:) )
      
     
     !--On reconstruit les tableaux ayant la meme grille , celle de B
        finy = deby + ncm(2) - 2
        finz = debz + ncm(3) - 2
        if ((finy.gt.ncm_tot(2)).or.(finz.gt.ncm_tot(3))) then
        print *,'finy ',finy,'finz ',finz,'ii ',ii,'jj ',jj
        stop
        endif
        Ax(1:ncm_tot(1),deby:finy,debz:finz) = Ax_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
        Ay(1:ncm_tot(1),deby:finy,debz:finz) = Ay_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
        Az(1:ncm_tot(1),deby:finy,debz:finz) = Az_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
     deby=deby+(ncm(2)-2)
     deallocate(Ax_proc,Ay_proc,Az_proc)
      stId = nf90_close(ncid);  call test_cdf(stId)
      write(*,*) 'close file ',filename
     enddo  
     debz=debz+(ncm(3)-2)   
     enddo
      

    ! Now we read electron pressure file
      !--Create file name for Field (READ, proc 0) 
      !-- Reading electron density file
      write(filename,'(a4,a1,i4.4,a1,2a)')trim(prefix3),"_",0,'_',trim(run_name),".nc"
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
      
          debz = 1
    do jj=1,nz
        deby = 1
    do ii=1,ny
      iproc=proc_tab(ii,jj)
      !--Create file name
      write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix3),'_',iproc-1,'_',trim(run_name)
      write(*,*) 'read file ',filename
      stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
      call test_cdf(stId)
      call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
      !--Get number point per processus
      call get_simple_variable_cdf(ncid,"ncm",itmp(:))
      call set_grid(itmp,2)
  
      !--allocation des tableaux de lectures du fichier diag de champ
      allocate(B_proc(ncm(1),ncm(2),ncm(3)))
      B_proc = zero ; 
      
      !--Get Total Number of Particles
      call get_simple_variable_cdf(ncid,"nptot",nptot_proc(iproc))
      call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )
      call get_simple_variable_cdf(ncid,varname5 ,B_proc(:,:,:) )
      
      
     
     !--On reconstruit les tableaux ayant la meme grille , celle de B
        finy = deby + ncm(2) - 2
        finz = debz + ncm(3) - 2
        if ((finy.gt.ncm_tot(2)).or.(finz.gt.ncm_tot(3))) then
        print *,'finy ',finy,'finz ',finz,'ii ',ii,'jj ',jj
        stop
        endif
        B(1:ncm_tot(1),deby:finy,debz:finz) = B_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1)
        deby=deby+(ncm(2)-2)
     deallocate(B_proc)
      stId = nf90_close(ncid);  call test_cdf(stId)
      write(*,*) 'close file ',filename
     enddo  
     debz=debz+(ncm(3)-2)   
     enddo      
      
      
     ! now we derive electron temperature from electron pressure and electron density
     B = B/A
     
     ! normalisation
     A = A*nrm_n
     Ax = Ax*nrm_U*sgn
     Ay = Ay*nrm_U*sgn
     Az = Az*nrm_U
     B = B*nrm_T
     
     
     
  
     
     !--Creation du fichier de diagnostique de champ a acces sequentiel
       write_name = prefix_output//'_'//trim(run_name)//".nc"
       call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
     
       !--Open the file for write Particles
       stId = nf90_create(trim(write_name), nf90_64BIT_OFFSET, ncid)
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
       stId = nf90_def_dim(ncid, "dimunit",10, dimchar)     
       call test_cdf(stId)
     
       stId = nf90_def_var(ncid, "s_min",   QP_NF90_DP,dimid(3), varid(ii))
       call test_cdf(stId); ii = ii+1 
       stId = nf90_def_var(ncid, "s_max",   QP_NF90_DP,dimid(3), varid(ii))
       call test_cdf(stId); ii = ii+1   
       stId = nf90_def_var(ncid, "unit n",   nf90_char,dimchar, varid(ii))
       call test_cdf(stId); ii = ii+1   
       stId = nf90_def_var(ncid, "unit U",   nf90_char,dimchar, varid(ii))
       call test_cdf(stId); ii = ii+1  
       stId = nf90_def_var(ncid, "unit T",   nf90_char,dimchar, varid(ii))
       call test_cdf(stId); ii = ii+1  
       stId = nf90_def_var(ncid, "X_axis",   QP_NF90_DP,dimid(6), varid(ii))
       call test_cdf(stId); ii = ii+1     
       stId = nf90_def_var(ncid, "Y_axis",   QP_NF90_DP,dimid(7), varid(ii))
       call test_cdf(stId); ii = ii+1     
       stId = nf90_def_var(ncid, "Z_axis",   QP_NF90_DP,dimid(8), varid(ii))
       call test_cdf(stId); ii = ii+1     
       stId = nf90_def_var(ncid, "Unit_axis",   nf90_char,dimchar, varid(ii))
       call test_cdf(stId); ii = ii+1     
       stId = nf90_def_var(ncid, "Coordinate_system",  nf90_char,dimchar, varid(ii))
       call test_cdf(stId); ii = ii+1      
      
     
     
       stId = nf90_def_var(ncid, "Density", QP_NF90_DP,dimid(6:8),varid(ii))
       call test_cdf(stId); ii = ii+1 
       stId = nf90_def_var(ncid, "Ux", QP_NF90_DP,dimid(6:8),varid(ii))
       call test_cdf(stId); ii = ii+1 
       stId = nf90_def_var(ncid, "Uy", QP_NF90_DP,dimid(6:8),varid(ii))
       call test_cdf(stId); ii = ii+1 
       stId = nf90_def_var(ncid, "Uz", QP_NF90_DP,dimid(6:8),varid(ii))
       call test_cdf(stId); ii = ii+1 
       stId = nf90_def_var(ncid, "Temperature", QP_NF90_DP,dimid(6:8),varid(ii))
       call test_cdf(stId); ii = ii+1 
       
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
         stId = nf90_put_var(ncid, varid(ii), unit_n)
         call test_cdf(stId); ii = ii+1              
         stId = nf90_put_var(ncid, varid(ii), unit_U)
         call test_cdf(stId); ii = ii+1            
         stId = nf90_put_var(ncid, varid(ii), unit_T)
         call test_cdf(stId); ii = ii+1     
         stId = nf90_put_var(ncid, varid(ii), X_axis)
         call test_cdf(stId); ii = ii+1        
         stId = nf90_put_var(ncid, varid(ii), Y_axis)
         call test_cdf(stId); ii = ii+1        
         stId = nf90_put_var(ncid, varid(ii), Z_axis)
         call test_cdf(stId); ii = ii+1        
         stId = nf90_put_var(ncid, varid(ii), Unit_axis)
         call test_cdf(stId); ii = ii+1       
         stId = nf90_put_var(ncid, varid(ii), coord)
         call test_cdf(stId); ii = ii+1 


       
         stId = nf90_put_var(ncid, varid(ii), A)
         call test_cdf(stId); ii = ii+1
         stId = nf90_put_var(ncid, varid(ii), Ax)
         call test_cdf(stId); ii = ii+1
         stId = nf90_put_var(ncid, varid(ii), Ay)
         call test_cdf(stId); ii = ii+1
         stId = nf90_put_var(ncid, varid(ii), Az)
         call test_cdf(stId); ii = ii+1
         stId = nf90_put_var(ncid, varid(ii), B)
         call test_cdf(stId); ii = ii+1
         
         !--Close the file
           stId = nf90_close(ncid); call test_cdf(stId)
         
     
  
     deallocate(A,Ax,Ay,Az,B)
     deallocate(voisin_proc,dims,coord_proc,nptot_proc)
     deallocate(X_axis,Y_axis,Z_axis)
     
    
  end subroutine merge_moment_thermal_cdf
  
  
  
 !!****************************************************
 !!-----------------------------------------------
 !! m_split_reassemble_cdf.F90
 !! reads 1 file computed for each process
 !! and merge it to one big file
 
  subroutine merge_field_cdf_atm(run_name,prefix,prefix_output)
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: prefix,prefix_output
    !--Reading variables
    integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
    integer :: rang,nb_procs,ndims,nb_voisin
    integer,allocatable :: dims(:),nptot_proc(:),dimspec(:)
    integer,allocatable :: voisin_proc(:,:)
    integer,allocatable :: coord_proc(:,:),proc_tab(:,:)
    
      
    !--Variables to read fields diagnostic
    real(dp),dimension(:,:,:,:),allocatable :: Ax_proc
    real(dp),dimension(:,:,:,:),allocatable :: Ax
    
    !--Others
    integer :: ind,deby,debz,finy,finz,iproc,ish
    integer :: deb,fin,cumul
    integer  :: stId,ncid,ii,is,jj,ny,nz
    integer  :: varid(1000),dimid(12)
    integer  :: itmp(3)
    real(dp) :: rtmp(3)
    logical :: file_e
    character(len=50) :: write_name
    character(len=64) :: filename
     character(len=100) :: test_name
    character(len=500) :: msg
    integer :: n_spe,i_spe
    character(len=11):: nam,nam_tmp
    character(len=11),dimension(:),allocatable:: name_dens
   
    !--Create file name for Field (READ, proc 0) 
    write(filename,'(a4,a1,i4.4,a1,2a)')trim(prefix),"_",0,'_',trim(run_name),".nc"
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
  
    call get_simple_variable_cdf(ncid,"n_species",n_spe)
  
    allocate(name_dens(n_spe))
    allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3),n_spe))
    Ax = zero
    !--allocation des tableaux de lectures du fichier diag de champ
     
    allocate(voisin_proc(nb_procs,nb_voisin)) ; voisin_proc = 0
  
    stId = nf90_close(ncid);  call test_cdf(stId)
    
  ! Maintenat que l'on connait le nombre de processus qui a tourne
   ! on lit les autres fichiers
   ny=0
   do iproc=1,nb_procs 
     write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
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
     write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
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
    write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
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
        do ish=1,n_spe
        write(nam,'(a,i2)')"Spe_",ish
                print *,'nam = ',nam
                call get_simple_variable_cdf(ncid,nam,nam_tmp)
                name_dens(ish)=nam_tmp
                call get_simple_variable_cdf(ncid,nam_tmp,Ax_proc(:,:,:,ish))
                write(*,*) 'read file ',nam_tmp
                Ax(1:ncm_tot(1),deby:finy,debz:finz,ish) = Ax_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,ish)
        enddo

   deby=deby+(ncm(2)-2)
   deallocate(Ax_proc)
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename
   enddo  
   debz=debz+(ncm(3)-2)   
   enddo
   nptot = sum(nptot_proc)  
   deallocate(proc_tab)
   
   !--Creation du fichier de diagnostique de champ a acces sequentiel
     write_name = prefix_output//'_'//trim(run_name)//".nc"
     call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
   
     !--Open the file for write Particles
     stId = nf90_create(trim(write_name), nf90_64BIT_OFFSET, ncid)
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
   
        do ish=1,n_spe
            stId = nf90_def_var(ncid,name_dens(ish), QP_NF90_DP,dimid(6:8),varid(ii))
            call test_cdf(stId); ii = ii+1   
        enddo
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

        do ish=1,n_spe
            stId = nf90_put_var(ncid, varid(ii), Ax(:,:,:,ish))
            call test_cdf(stId); ii = ii+1   
        enddo
        ii = ii-1   
       
       !--Close the file
         stId = nf90_close(ncid); call test_cdf(stId)
         
       
   
   deallocate(Ax)
   deallocate(voisin_proc,dims,coord_proc,nptot_proc)
   deallocate(name_dens) 
  
  end subroutine merge_field_cdf_atm
  
  
!****************************************************
 !!-----------------------------------------------
 !! m_merge_global_cdf/merge_moment_species_cdf
 !! reads 1 file computed for each process
 !! and merge it to one big file
 
  subroutine merge_moment_species_cdf(run_name,prefix,prefix_output)
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: prefix,prefix_output
    !--Reading variables
    integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
    integer :: rang,nb_procs,ndims,nb_voisin
    integer,allocatable :: dims(:),nptot_proc(:),dimspec(:)
    integer,allocatable :: voisin_proc(:,:)
    integer,allocatable :: coord_proc(:,:),proc_tab(:,:)
    
      
    !--Variables to read fields diagnostic
    real(dp),dimension(:,:,:,:),allocatable :: A_proc,Ax_proc,Ay_proc,Az_proc,B_proc
    real(dp),dimension(:,:,:,:),allocatable :: A,Ax,Ay,Az,B
    real(dp),dimension(:),allocatable :: X_axis,Y_axis,Z_axis
    
    !--Others
    integer :: ind,deby,debz,finy,finz,iproc,ish
    integer :: deb,fin,cumul
    integer  :: stId,ncid,ii,is,ny,nz,jj,iii,i
    integer  :: varid(1000),dimid(12)
    integer :: dimion(5),dimchar
    integer  :: itmp(3),var_id
    real(dp) :: rtmp(3)
    logical :: file_e
    character(len=50) :: write_name
    character(len=64) :: filename
     character(len=100) :: test_name
    character(len=500) :: msg
    integer :: nb_spe,i_spe,iunit,sgn
    character(len=11):: nam,nam_tmp
    character(len=10),dimension(:),allocatable:: ion_label
    character(len=50),dimension(:),allocatable:: diag_name
    character(len=10) :: unit_n,unit_U,unit_T,unit_axis,coord
    real(dp) :: nrm_n,nrm_U,nrm_T,v0,n0,x0
    character(len=8) :: planet 
   
    !--Create file name for Field (READ, proc 0) 
    write(filename,'(a4,a1,i4.4,a1,2a)')trim(prefix),"_",0,'_',trim(run_name),".nc"
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
    
    !--Get Planetname
    stId = nf90_inq_varid(ncid,"planetname",var_id)
    stId = nf90_get_var(ncid,var_id,planet)
   write(*,*) 'Planet :',planet    
    
    !--Get the number of processus which generated the file
        call get_simple_variable_cdf(ncid,"nproc",nb_procs)
        write(*,*) 'nb_procs : ',nb_procs
      
        !--Get the rang of the open file 
        call get_simple_variable_cdf(ncid,"mpiinfo_me",rang)
        write(*,*) 'rang ',rang
      
        !--Get topology dimension
        call get_simple_dimens_cdf(ncid,"mpi_ndims",ndims)
        write(*,*) ' ndims ',ndims
      
        !--Get number of neighbours
        call get_simple_dimens_cdf(ncid,"mpi_nb_voisin",nb_voisin)
        write(*,*) 'mpi_nb_voisin ',nb_voisin
      
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
        write(*,*) 'nxyzm ',nxyzm
      
        allocate(dims(ndims)) ; dims = 0
      
        !--Define coordinate vector for proc
        allocate(coord_proc(nb_procs,ndims)) ; coord_proc =0
      
        !--Get the grid step
        call get_simple_variable_cdf(ncid,"gstep",gstep)
        write(*,*) 'gstep ',gstep
      
        !--Get the dimensions of the box
        !--MIN
        call get_simple_variable_cdf(ncid,"s_min",s_min)
        !--MAX
        call get_simple_variable_cdf(ncid,"s_max",s_max)
      
        !--Get Time Informations
        !--Iteration
        call get_simple_variable_cdf(ncid,"iter",iter)
        write(*,*) 'iter ',iter
        !--Time interval
        call get_simple_variable_cdf(ncid,"dt_t",rtmp(:2))
        dt = rtmp(1);   t = rtmp(2)
        
      !--Allocate species and get informations
      call alloc_species(ns,Spe)
      call species_get_var_cdf(Spe,ncid)   
      
          !--Get Number of Species
      call get_simple_variable_cdf(ncid,"ns",ns)
      
        !--Get Total Number of Particles for this proc
        allocate(nptot_proc(nb_procs)) ; nptot_proc = 0
        call get_simple_variable_cdf(ncid,"nptot",nptot_proc(1))
        write(*,*) 'nptot proc (1)', nptot_proc(1)
        
        call get_simple_variable_cdf(ncid,"nb_ion_species",nb_spe)
        write(*,*) 'nb_species ',nb_spe
        
        call get_simple_variable_cdf(ncid,"phys_length",x0)
  
        
        !-- allocation du tableaux de nom des especes
        allocate(ion_label(nb_spe))     ;       ion_label(:) = ''
        call get_simple_variable_cdf(ncid,"ion_label",ion_label)
        write(*,*) 'Ion label',ion_label
     
        allocate(voisin_proc(nb_procs,nb_voisin)) ; voisin_proc = 0
        
        ! coordinate transformation for each object
        select case (trim(planet))
        case ("mars","mercury")
          sgn = -1
          coord="MSO"
        case ("ganymede")
          sgn=+1
          coord="GPHiO"
        case("titan")
          sgn = +1
          coord="TIIS"
        case("earth")  
          sgn = -1
          coord = "GSM"
        case default 
          planet="mars    "
          sgn = -1
       end select        
        
   !-- Allocate coordinate axis
   allocate(X_axis(ncm_tot(1)));        allocate(Y_axis(ncm_tot(2)))
   allocate(Z_axis(ncm_tot(3)))
   X_axis = zero;       Y_axis = zero;  Z_axis = zero
   
   X_axis = sgn*((/(i,i=1,ncm_tot(1))/)*gstep(1)-Spe%P%centr(1))*x0
   Y_axis = sgn*((/(i,i=1,ncm_tot(2))/)*gstep(2)-Spe%P%centr(2))*x0
   Z_axis = ((/(i,i=1,ncm_tot(3))/)*gstep(3)-Spe%P%centr(3))*x0
   unit_axis = "km"
           

           !-- allocation du tableau general
           allocate(A(ncm_tot(1),ncm_tot(2),ncm_tot(3),nb_spe));        A = zero
           allocate(B(ncm_tot(1),ncm_tot(2),ncm_tot(3),nb_spe));        B = zero
           allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3),nb_spe));       Ax = zero
           allocate(Ay(ncm_tot(1),ncm_tot(2),ncm_tot(3),nb_spe));       Ay = zero
           allocate(Az(ncm_tot(1),ncm_tot(2),ncm_tot(3),nb_spe));       Az = zero
   A = zero
   Ax = zero
   Ay = zero
   Az = zero
   B = zero

   allocate(diag_name(nb_spe));diag_name = ""
   

     
      ! get normalization value
      !--Get normalisation values
      call get_simple_variable_cdf(ncid,"phys_density",n0)
      call get_simple_variable_cdf(ncid,"phys_speed",V0)    

     nrm_n = n0*1.e-6!to have it in cm-3
     unit_n = "cm-3"
     nrm_U = V0     ! to have it in km/s
     unit_U = "km.s-1"
     nrm_T = (V0*1.e3)**2*amu_pmass/(3._dp*kb_JK*11600._dp) ! to have it in eV
     !nrm_T = (Spe%ref%mag**2)/(Spe%ref%density*2._dp*mu0*kb_JK)
     unit_T = "eV"
     select case (trim(planet))
     case("ganymede","titan")
       nrm_T = nrm_T*16._dp ! to take into account that normalisation is O+ mass
     end select
      
      
        stId = nf90_close(ncid);  call test_cdf(stId)


   ny=0
   do iproc=1,nb_procs 
     write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
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
     write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
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
    write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
    write(*,*) 'read file ',filename
    stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
    call test_cdf(stId)
    call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))
    !--Get number point per processus
    call get_simple_variable_cdf(ncid,"ncm",itmp(:))
    call set_grid(itmp,2) 
        !--allocation des tableaux de lectures du fichier diag de moment
        allocate(A_proc(ncm(1),ncm(2),ncm(3),nb_spe))
        allocate(B_proc(ncm(1),ncm(2),ncm(3),nb_spe))
        allocate(Ax_proc(ncm(1),ncm(2),ncm(3),nb_spe))
        allocate(Ay_proc(ncm(1),ncm(2),ncm(3),nb_spe))
        allocate(Az_proc(ncm(1),ncm(2),ncm(3),nb_spe))
        A_proc = zero; Ax_proc = zero;Ay_proc = zero;Az_proc = zero 
        B_proc = zero
 
         !--Get Total Number of Particles
         call get_simple_variable_cdf(ncid,"nptot",nptot_proc(iproc))
         call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )
         call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))        
        
         call get_simple_variable_cdf(ncid,"density_species" ,A_proc(:,:,:,:) )
         call get_simple_variable_cdf(ncid,"Vx_species" ,Ax_proc(:,:,:,:) )
         call get_simple_variable_cdf(ncid,"Vy_species" ,Ay_proc(:,:,:,:) )
         call get_simple_variable_cdf(ncid,"Vz_species" ,Az_proc(:,:,:,:) )
         call get_simple_variable_cdf(ncid,"Temp_species" ,B_proc(:,:,:,:) )
         finy = deby + ncm(2) - 2
         finz = debz + ncm(3) - 2
    stId = nf90_close(ncid);  call test_cdf(stId)
         A(1:ncm_tot(1),deby:finy,debz:finz,:) = &
        & A(1:ncm_tot(1),deby:finy,debz:finz,:) + A_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,:)
         Ax(1:ncm_tot(1),deby:finy,debz:finz,:) = &
        & Ax_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,:)     
         Ay(1:ncm_tot(1),deby:finy,debz:finz,:) = &
        & Ay_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,:)
         Az(1:ncm_tot(1),deby:finy,debz:finz,:) = & 
        & Az_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,:)
         B(1:ncm_tot(1),deby:finy,debz:finz,:) = B_proc(1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,:)
         
     deby=deby+(ncm(2)-2)
        deallocate(A_proc,Ax_proc,Ay_proc,Az_proc,B_proc)
    write(*,*) 'close file ',filename
   enddo  
   debz=debz+(ncm(3)-2)   
   enddo
   nptot = sum(nptot_proc)  
   deallocate(proc_tab)

! normalization
A = A*nrm_n
Ax = sgn*Ax*nrm_U
Ay = sgn*Ay*nrm_U
Az = Az*nrm_U
B = B*nrm_T
       
  do i_spe=1,nb_spe
   !--Create file name for Field (READ, proc 0) 
   write(filename,'(a,a1,2a)')trim(ion_label(i_spe)),"_",trim(run_name),".nc"
   write(*,*) " ======= Creation of file : "//trim(filename)
   diag_name(i_spe) = trim(filename)
   
   !--Open the file for write moments
    stId = nf90_create(trim(filename), nf90_64BIT_OFFSET, ncid)
    call test_cdf(stId)
   
    call set_global_attribute_cdf(ncid,"Moments")
    
      !-- Information relative a l'ecriture des fichiers
      !--Define the dimensions that will define the size of the file.
      call  create_wrt_dimensions_cdf(ncid,dimid,ncm_tot)
    
      !--Define Speces dimensions
      call species_def_dim_cdf(ncid,dimspec)
      
       iii = 1
       !--Define commons variables
       call common_def_var_cdf(ncid,varid,dimid,iii)
     
       !--Define MPI_INFO
       call mpiinfo_def_var_cdf(ncid,varid,dimid,iii)
     
       !--Define Speces variables
       call species_def_var_cdf(ncid,varid,dimspec,iii)
       
     stId = nf90_def_dim(ncid, "dimunit",10, dimchar)     
     call test_cdf(stId)       
     
       !--Define other variables
       !--GRID_INFO
       stId = nf90_def_var(ncid, "nxyzm",  nf90_int,dimid(1), varid(iii))
       call test_cdf(stId); iii = iii+1 
       stId = nf90_def_var(ncid, "ncm_tot",nf90_int,dimid(3), varid(iii))
       call test_cdf(stId); iii = iii+1 
       stId = nf90_def_var(ncid, "nc_tot", nf90_int,dimid(3),varid(iii))
       call test_cdf(stId); iii = iii+1 
       stId = nf90_def_var(ncid, "ncm",    nf90_int,dimid(3),varid(iii))
       call test_cdf(stId); iii = iii+1 
       stId = nf90_def_var(ncid, "s_min",  QP_NF90_DP,dimid(3),  varid(iii))
       call test_cdf(stId); iii = iii+1 
       stId = nf90_def_var(ncid, "s_max",  QP_NF90_DP,dimid(3),  varid(iii))
       call test_cdf(stId); iii = iii+1 
       stId = nf90_def_var(ncid, "unit n",  nf90_char,dimchar,  varid(iii))
       call test_cdf(stId); iii = iii+1 
       stId = nf90_def_var(ncid, "unit U",  nf90_char,dimchar,  varid(iii))
       call test_cdf(stId); iii = iii+1
       stId = nf90_def_var(ncid, "unit T",  nf90_char,dimchar,  varid(iii))
       call test_cdf(stId); iii = iii+1   
       stId = nf90_def_var(ncid, "X_axis",   QP_NF90_DP,dimid(6), varid(iii))
       call test_cdf(stId); iii = iii+1     
       stId = nf90_def_var(ncid, "Y_axis",   QP_NF90_DP,dimid(7), varid(iii))
       call test_cdf(stId); iii = iii+1     
       stId = nf90_def_var(ncid, "Z_axis",   QP_NF90_DP,dimid(8), varid(iii))
       call test_cdf(stId); iii = iii+1     
       stId = nf90_def_var(ncid, "Unit_axis",   nf90_char,dimchar, varid(iii))
       call test_cdf(stId); iii = iii+1     
       stId = nf90_def_var(ncid, "Coordinate_system",  nf90_char,dimchar, varid(iii))
       call test_cdf(stId); iii = iii+1       
        

       !--PARTICLES_INFO
       stId = nf90_def_var(ncid, "npm",               nf90_int,dimid(1),varid(iii))
       call test_cdf(stId); iii = iii+1 
       
       !--Nb_ion species
       stId = nf90_def_var(ncid, "nb_ion_species",nf90_int,dimid(1),varid(iii))
       call test_cdf(stId); iii = iii+1 
       
       !-- Density
       stId = nf90_def_var(ncid,"Density",QP_NF90_DP,dimid(6:8),varid(iii))
       call test_cdf(stId); iii = iii +1
       
       !-- Vx component
       stId = nf90_def_var(ncid,"Ux",QP_NF90_DP,dimid(6:8),varid(iii))
       call test_cdf(stId); iii = iii +1
       !-- Vy component
       stId = nf90_def_var(ncid,"Uy",QP_NF90_DP,dimid(6:8),varid(iii))
       call test_cdf(stId); iii = iii +1
       !-- Vz component
       stId = nf90_def_var(ncid,"Uz",QP_NF90_DP,dimid(6:8),varid(iii))
       call test_cdf(stId); iii = iii +1
       !-- Temperature 
       stId = nf90_def_var(ncid,"Temperature",QP_NF90_DP,dimid(6:8),varid(iii))
       call test_cdf(stId); iii = iii +1       
       
       !--Switch to write mode
       stId = nf90_enddef(ncid); call test_cdf(stId)
  
       iii = 1   
       !--Write common variables into the file
       call common_put_var_cdf(ncid,varid,iii)
 
       !--Write mpi variables into the file
       call mpiinfo_put_var_cdf(ncid,varid,iii)
 
       !--Write Species infos into the file
       call species_put_var_cdf(Spe,ncid,varid,dimspec,iii)

       !--Write the other variables into the file
       stId = nf90_put_var(ncid, varid(iii), nxyzm)
       call test_cdf(stId); iii = iii+1 
       stId = nf90_put_var(ncid, varid(iii), ncm_tot)
       call test_cdf(stId); iii = iii+1 
       stId = nf90_put_var(ncid, varid(iii), nc_tot)
       call test_cdf(stId); iii = iii+1 
       stId = nf90_put_var(ncid, varid(iii), ncm)
       call test_cdf(stId); iii = iii+1 
       stId = nf90_put_var(ncid, varid(iii), s_min)
       call test_cdf(stId); iii = iii+1              
       stId = nf90_put_var(ncid, varid(iii), s_max)
       call test_cdf(stId); iii = iii+1              
       stId = nf90_put_var(ncid, varid(iii), unit_n)
       call test_cdf(stId); iii = iii+1     
       stId = nf90_put_var(ncid, varid(iii), unit_U)
       call test_cdf(stId); iii = iii+1            
       stId = nf90_put_var(ncid, varid(iii), unit_T)
       call test_cdf(stId); iii = iii+1        
       stId = nf90_put_var(ncid, varid(iii), X_axis)
       call test_cdf(stId); iii = iii+1        
       stId = nf90_put_var(ncid, varid(iii), Y_axis)
       call test_cdf(stId); iii = iii+1        
       stId = nf90_put_var(ncid, varid(iii), Z_axis)
       call test_cdf(stId); iii = iii+1        
       stId = nf90_put_var(ncid, varid(iii), Unit_axis)
       call test_cdf(stId); iii = iii+1       
       stId = nf90_put_var(ncid, varid(iii), coord)
       call test_cdf(stId); iii = iii+1       
       stId = nf90_put_var(ncid, varid(iii), npm)
       call test_cdf(stId); iii = iii+1        
       stId = nf90_put_var(ncid, varid(iii), nb_spe)
       call test_cdf(stId); iii = iii+1      
       stId = nf90_put_var(ncid, varid(iii), A(:,:,:,i_spe))
       call test_cdf(stId); iii = iii+1 
       stId = nf90_put_var(ncid, varid(iii), Ax(:,:,:,i_spe))
       call test_cdf(stId); iii = iii+1 
       stId = nf90_put_var(ncid, varid(iii), Ay(:,:,:,i_spe))
       call test_cdf(stId); iii = iii+1 
       stId = nf90_put_var(ncid, varid(iii), Az(:,:,:,i_spe))
       call test_cdf(stId); iii = iii+1        
       stId = nf90_put_var(ncid, varid(iii), B(:,:,:,i_spe))
       call test_cdf(stId); iii = iii+1           
       !--Close the file
       stId = nf90_close(ncid); call test_cdf(stId)       
   enddo   


   
   iunit = 10
   write(write_name,'(a,a1,2a)')"Read_moment_species","_",trim(run_name),".dat"
   open(unit = iunit,file = write_name,form = 'formatted',status  = 'unknown')
   write(iunit,*) nb_spe
   write(iunit,*) ion_label
   do ii = 1,nb_spe
     write(iunit,*) diag_name(ii)
   enddo
   close(iunit)
   
   
        
        deallocate(A,B,Ax,Ay,Az)
        deallocate(voisin_proc,dims,coord_proc,nptot_proc)  
        deallocate(diag_name)
        deallocate(X_axis,Y_axis,Z_axis)
        
end subroutine merge_moment_species_cdf    
    
    
    
    

!!############################################################
!!****************************************************
 !!-----------------------------------------------
 !! m_split_reassemble_cdf.F90
 !! reads 1 file computed for each process
 !! and merge it to one big file
 
  subroutine merge_field_cdf_pro(run_name,prefix,prefix_output)
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: prefix,prefix_output
    !--Reading variables
    integer,parameter :: N=1,NE=2,E=3,SE=4,S=5,SW=6,W=7,NW=8
    integer :: rang,nb_procs,ndims,nb_voisin
    integer,allocatable :: dims(:),nptot_proc(:),dimspec(:)
    integer,allocatable :: voisin_proc(:,:)
    integer,allocatable :: coord_proc(:,:)
    
      
    !--Variables to read fields diagnostic
    real(dp),dimension(:,:,:,:,:),allocatable :: Ax_proc
    real(dp),dimension(:,:,:,:),allocatable :: Ax
    
    !--Others
    integer :: ind,deby,debz,finy,finz,iproc,ish
    integer :: deb,fin,cumul
    integer  :: stId,ncid,ii,is
    integer  :: varid(1000),dimid(12)
    integer  :: itmp(3)
    real(dp) :: rtmp(3)
    logical :: file_e
    character(len=50) :: write_name
    character(len=64) :: filename
     character(len=100) :: test_name
    character(len=500) :: msg
    integer :: n_spe_pp,i_spe
    character(len=11):: nam,nam_tmp
    character(len=11),dimension(:),allocatable:: name_dens
   
    !--Create file name for Field (READ, proc 0) 
    write(filename,'(a4,a1,i4.4,a1,2a)')trim(prefix),"_",0,'_',trim(run_name),".nc"
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
  
    call get_simple_variable_cdf(ncid,"n_spe_pp",n_spe_pp)
  
  
    !--allocation des tableaux de lectures du fichier diag de champ
    allocate(Ax_proc(nb_procs,ncm(1),ncm(2),ncm(3),n_spe_pp))
    Ax_proc = zero 
     
    allocate(voisin_proc(nb_procs,nb_voisin)) ; voisin_proc = 0
  
    stId = nf90_close(ncid);  call test_cdf(stId)
    
  ! Maintenat que l'on connait le nombre de processus qui a tourne
   ! on lit les autres fichiers
   do iproc=1,nb_procs 
 
    !--Create file name
    write(filename,'(a4,a1,i4.4,a1,a)') trim(prefix),'_',iproc-1,'_',trim(run_name)
    write(*,*) 'read file ',filename
    !--Open NetCDF file for read fields 
    stId = nf90_open(trim(filename)//".nc", nf90_nowrite, ncid)
    call test_cdf(stId)
 
    !--Get Total Number of Particles
    call get_simple_variable_cdf(ncid,"nptot",nptot_proc(iproc))
    call get_simple_variable_cdf(ncid,"mpiinfo_voisin" ,voisin_proc(iproc,:) )
    call get_simple_variable_cdf(ncid,"mpiinfo_coord",coord_proc(iproc,:ndims))

    allocate(name_dens(n_spe_pp))
        do ish=1,n_spe_pp
        write(nam,'(a,i2)')"Spe_",ish
        call get_simple_variable_cdf(ncid,nam,nam_tmp)
        name_dens(ish)=nam_tmp
        call get_simple_variable_cdf(ncid,nam_tmp,Ax_proc(iproc,:,:,:,ish))
        write(*,*) 'read file ',filename,nam_tmp
        enddo

    
    stId = nf90_close(ncid);  call test_cdf(stId)
    
   end do
 
   nptot = sum(nptot_proc)  
   
   allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3),n_spe_pp))
   Ax = zero
   
     !--Allocation of global arrays
   !call alloc_arr3D(Bfield,ncm_tot)


   !--On reconstruit les tableaux ayant la meme grille , celle de B
   do ish=1,n_spe_pp  
   do ind = 1,nb_procs
      deby = coord_proc(ind,ndims-1)*(ncm(2)-2) + 1
      debz = coord_proc(ind,ndims)*(ncm(3)-2) + 1
      finy = deby + ncm(2) - 2
      finz = debz + ncm(3) - 2
      ! print *, "COORD PROC",ind,coord_proc(ind,:),ndims
      !    print '("P ",i3, " deby,finy ",2(i5)," debz,finz ",2(i5))', ind,deby,finy,debz,finz
      Ax(1:ncm_tot(1),deby:finy,debz:finz,ish) = Ax_proc(ind,1:ncm_tot(1),1:ncm(2)-1,1:ncm(3)-1,ish)
   enddo
   enddo

   
   !--Creation du fichier de diagnostique de champ a acces sequentiel
     write_name = prefix_output//'_'//trim(run_name)//".nc"
     call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
   
     !--Open the file for write Particles
     stId = nf90_create(trim(write_name), nf90_64BIT_OFFSET, ncid)
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
   
        do ish=1,n_spe_pp
            stId = nf90_def_var(ncid,name_dens(ish), QP_NF90_DP,dimid(6:8),varid(ii))
            call test_cdf(stId); ii = ii+1   
        enddo
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

        do ish=1,n_spe_pp
            stId = nf90_put_var(ncid, varid(ii), Ax(:,:,:,ish))
            call test_cdf(stId); ii = ii+1   
        enddo
        ii = ii-1   
       
       !--Close the file
         stId = nf90_close(ncid); call test_cdf(stId)
         
       
   
   deallocate(Ax_proc)
   deallocate(Ax)
   deallocate(voisin_proc,dims,coord_proc,nptot_proc)
   deallocate(name_dens)
  
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
  
  write(name_file,'(a4,a1,i4.4,a1,a)')trim(prefix),"_",me,'_',trim(filwrt)
  
#ifdef HAVE_NETCDF
    name_file = trim(name_file)//".nc"
#endif
    
   
 end subroutine create_file_diag
 
#endif
 end module m_merge_global_cdf
