!!===================================================
!!===================================================
module m_IMPEX_moment_orbit
!! problem with planet case mars

 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype
 use defs_particletype
 use defs_species
 use defs_parametre,only : fildat
 use m_writeout
 use m_VO
 use ImpexTreeXML_generator
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
      value_along_orbit,	  &
      value_species_along_orbit, &
      save_value_orbit, &
      save_value_species_orbit, &
      density_along_orbit, &
      save_density_orbit, &
      header_orbit_VOTABLE
      !extract_XY,	  &
      !extract_XZ
      !create_file_diag,   &
      !read_moment_species_cdf,&
      !read_field_cdf_1_array

 public ::                &
      read_traj_file,	&
      read_file_size,	&
      read_field_cdf,	&
      read_dim_field_cdf, &
      read_density_cdf, &
      read_species_cdf,	&
      extract_1D
contains
 !********************************************************************
 ! Auteur				:	 RModolo, SHess
 ! Date					:	 05/04/12
 ! Institution				:	LATMOS/CNRS/IPSL
 ! Derniere modification		:	05/04/12	
 ! Resume	
 ! 
 !********************************************************************
 
  !!##############################################################################
   !! m_IMPEX_moment_orbit/read_traj_file
   !! This routine reads the trajectory
   !! information. Outputs are S/C position in 
   !! simulation coordinate system
   subroutine read_traj_file(traj_name,X_MSO,Y_MSO,Z_MSO,YYYY,MM,DD,hh,minut,ss)
   character(len=50),intent(in) :: traj_name
   real(dp),dimension(:),intent(inout) :: X_MSO,Y_MSO,Z_MSO
   integer,dimension(:),intent(inout) :: YYYY,MM,DD,hh,minut
   real(dp),dimension(:),intent(inout) :: ss
   ! Local variables
   integer :: iunit,i, n_traj,j
   real(dp) :: read_second,read_x,read_y,read_z
   integer :: read_year,read_month,read_day,read_hour,read_minute
   integer :: Reason
   
   iunit = 1
   n_traj = size(X_MSO,1)
   open(UNIT = iunit, FILE = traj_name, STATUS = 'OLD', ACTION = 'READ', &
        FORM = 'FORMATTED')

   Reason = 0
   do i=1,n_traj
     read(iunit,*) read_year,read_month,read_day,read_hour,read_minute,read_second,read_x,read_y,read_z
     X_MSO(i) = read_x
     Y_MSO(i) = read_y
     Z_MSO(i) = read_z
     YYYY(i)  = read_year
     MM(i)    = read_month
     DD(i)    = read_day
     hh(i)    = read_hour
     minut(i) = read_minute
     ss(i)    = read_second     
   enddo
   close(iunit)


   print *,' Max & Min location of the trajectory'
   print *, 'X MSO',maxval(X_MSO),minval(X_MSO)
   print *, 'Y MSO',maxval(Y_MSO),minval(Y_MSO)
   print *, 'Z MSO',maxval(Z_MSO),minval(Z_MSO)
   
   end subroutine read_traj_file
   
   !!#############################################################
   !!m_IMPEX_moment_orbit/read_file_size
   !! This routine reads the trajectory
   !! file and send the number of lines
   subroutine read_file_size(traj_name,n_traj)
   character(len=50),intent(in) :: traj_name
   integer,intent(out) :: n_traj
   !Local Variables
   integer :: Reason,j,iunit
   real(dp) :: read_second,read_x,read_y,read_z
   integer :: read_year,read_month,read_day,read_hour,read_minute
      
      iunit = 1
      open(UNIT = iunit, FILE = traj_name, STATUS = 'OLD', ACTION = 'READ', &
           FORM = 'FORMATTED')
   
      j = 0
      Reason = 0
      do while (Reason==0)
        read(iunit,*,IOSTAT=Reason) read_year,read_month,read_day,read_hour,read_minute,read_second,read_x,read_y,read_z
        !  print *, read_line,read_time,read_d,read_x,read_y,read_z
        j = j+1
      enddo
      n_traj = j-1   
   close(iunit)
   end subroutine read_file_size
   
   !!#############################################################
   !!m_IMPEX_moment_orbit/extract_1D
   !! This routine reads the global
   !! output simulation file and extract
   !! moments and fields along the trajectory
   subroutine extract_1D(traj_name,run_name,X0_traj,Y0_traj,Z0_traj,YYYY,MM,DD,hh,minut,ss,angle,spacecraft)
   character(len=50),intent(in) :: traj_name
   character(len=20),intent(in) :: run_name
   character(len=*),intent(in) :: spacecraft
   real(dp),dimension(:),intent(inout) :: X0_traj,Y0_traj,Z0_traj
   integer,dimension(:),intent(in) :: YYYY,MM,DD,hh,minut
   real(dp),dimension(:),intent(in) :: ss   
   real(dp),intent(in) :: angle
   ! Local varaibales
    character(len=50) :: write_name
   real(dp),dimension(:,:,:),allocatable :: A,Ax,Ay,Az,B
   real(dp),dimension(:),allocatable :: A_1D,Ax_1D,Ay_1D,Az_1D,B_1D
   character(len=64) :: var_name1,var_name2,var_name3,var_name4,var_name5
   character(len=4)  :: prefix
   integer :: ncm_tot(3)
   real(dp)  :: centr(3),radius,gstep(3),n0,V0,B0,nrm_n,nrm_U,nrm_T,nrm,box_size(6),orient(3),c_wpi
   character(len=8) :: planet
   integer :: nb_point,nb_spe,ii,iunit,sgn
    character(len=5) :: coord   
    character(len=10),dimension(:),allocatable:: ion_label
    character(len=50) :: diag_name_the
    character(len=50),dimension(:),allocatable:: diag_name
   character(len=600)::numdatid
   real(dp),dimension(size(Y0_traj)) :: X_traj,Y_traj,Z_traj
  
  
   if (angle /= -1) then  
   X_traj=X0_traj
   Y_traj=Y0_traj*cos((angle+90)/180.*pi)-Z0_traj*sin((angle+90)/180.*pi)
   Z_traj=Z0_traj*cos((angle+90)/180.*pi)+Y0_traj*sin((angle+90)/180.*pi)
   else
   X_traj = X0_traj
   Y_traj = Y0_traj
   Z_traj = Z0_traj
   endif

   print *,' Max & Min location of the trajectory after rotation'
   print *, 'angle',angle
   print *, 'X ',maxval(X_traj),minval(X_traj)
   print *, 'Y ',maxval(Y_traj),minval(Y_traj)
   print *, 'Z ',maxval(Z_traj),minval(Z_traj)
   
   prefix = "Magw"
   var_name1 = "Bx"
   var_name2 = "By"
   var_name3 = "Bz"
      
   call read_dim_field_cdf(run_name,prefix, var_name1,ncm_tot,radius,centr,gstep,planet,c_wpi)
   
   print *,'Dimension array :',ncm_tot
   print *,'Spatial step :',gstep
   print *,'Obstacle position :',centr
   print *,'Obstacle radius :',radius
   print *,'Planet :',planet
   
   box_size=(/ -centr(1),ncm_tot(1)*gstep(1)-centr(1),-centr(2),ncm_tot(2)*gstep(2)-centr(2),-centr(3),ncm_tot(3)*gstep(3)-centr(3) /)
   box_size=box_size*c_wpi
       
   ! coordinate transformation for each object
   select case (trim(planet))
   case ("mars","mercury")
     sgn = -1
     coord='MSO'
   case ("ganymede")
     sgn=+1
     coord = 'GPhiO'
   case("titan")
     sgn = +1
     coord='TIIS'
     case default 
     planet="mars    "
     sgn = -1
     coord = 'MSO'
    end select
   
   !--Convert MSO trajectory in simulation unit and coordinate
   X_traj = (sgn*X_traj*radius + centr(1))
   Y_traj = (sgn*Y_traj*radius + centr(2))
   Z_traj =  (Z_traj*radius + centr(3))
   
   nb_point = size(X_traj,1)
   
   
   !--Write tree.xml metadata
   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.true.,.false.,.false.,.false.,'',0.,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
   call write_granule_3D_XML(prefix,run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)

   !--Allocating arrays : will read global output files
   allocate(A(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   allocate(Ay(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   allocate(Az(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   allocate(B(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
   A = zero;	Ax = zero;	Ay = zero;	Az = zero
   B = zero
   
   !-- Allocating arrays for values along trajectory
   allocate(Ax_1D(nb_point),Ay_1D(nb_point),Az_1D(nb_point))
   Ax_1D = zero;	Ay_1D = zero;	Az_1D = zero
   allocate(A_1D(nb_point),B_1D(nb_point))
   A_1D = zero;		B_1D = zero
   
   
   !-- Read Magw file
   call read_field_cdf(run_name,prefix,var_name1,var_name2,var_name3, Ax,Ay,Az,n0,V0,B0)
   nrm = 1 ! Magw file is already normalized in nT
   call value_along_orbit(X_traj,Y_traj,Z_traj,gstep,Ax,Ay,Az,Ax_1D,Ay_1D,Az_1D,sgn,nrm)
   A_1D = sqrt(Ax_1D**2+Ay_1D**2+Az_1D**2)
   call save_value_orbit(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss, &
   		A_1D,Ax_1D,Ay_1D,Az_1D,planet,coord,angle)

   !--Write tree.xml metadata
   call write_numdat_XML(fildat,planet,prefix(1:3),.true.,.false.,.false.,.false.,.false.,trim(spacecraft),angle,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
   call write_granule_TS_XML(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss,numdatid,angle)

   		

   !-- Read Thew file
    prefix = "Thew"
   var_name1 = "Density"
   var_name2 = "Ux"
   var_name3 = "Uy"
   var_name4 = "Uz"
   var_name5 = "Temperature"
   !--Write tree.xml metadata
   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.true.,.false.,.false.,.false.,'',0.,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
!	minval(X_traj),maxval(X_traj),minval(Y_traj),maxval(Y_traj),minval(Z_traj),minval(Z_traj))   
   call write_granule_3D_XML(prefix,run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)   
   
    A = zero;	B = zero
    Ax = zero;	Ay = zero;	Az = zero
    A_1D = zero;	B_1D = zero
    Ax_1D = zero;	Ay_1D = zero;	Az_1D = zero
    
    write(diag_name_the,'(a4,a1,2a)')trim(prefix),"_",trim(run_name),".nc"    
     call read_species_cdf(diag_name_the,var_name1,var_name2,var_name3,var_name4,var_name5,A, Ax,Ay,Az,B)
     
    nrm_n =1 !density is already  in cm-3
    nrm_U =1     ! speed is already in km/s
    nrm_T = 1 ! temperature is already in eV
     call value_species_along_orbit(X_traj,Y_traj,Z_traj,gstep,A,Ax,Ay,Az,B,A_1D,Ax_1D,Ay_1D,Az_1D,B_1D, &
     	sgn,nrm_n,nrm_U,nrm_T) 	
    call save_value_species_orbit(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss, &
   		A_1D,Ax_1D,Ay_1D,Az_1D,B_1D,planet,coord,angle)    	
     

   !--Write tree.xml metadata
   print *,'Spacecraft name :',spacecraft
   call write_numdat_XML(fildat,planet,prefix(1:3),.true.,.false.,.false.,.false.,.false.,trim(spacecraft),angle,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
   call write_granule_TS_XML(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss,numdatid,angle)
   
   
    
!   !--Write tree.xml metadata
!   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.true.,.false.,.false.,.false.,'',0.,&
!	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
!!	minval(X_traj),maxval(X_traj),minval(Y_traj),maxval(Y_traj),minval(Z_traj),minval(Z_traj))   
!   call write_granule_3D_XML(prefix,run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
!
!    Ax = zero;	Ay = zero;	Az = zero
!    Ax_1D = zero;	Ay_1D = zero;	Az_1D = zero
!    A_1D = zero
!    
!   call read_field_cdf(run_name,prefix,var_name1,var_name2,var_name3, Ax,Ay,Az,n0,V0,B0)
!   nrm = 1 ! Plasma speed is already normalized in km/s   
!   call value_along_orbit(X_traj,Y_traj,Z_traj,gstep,Ax,Ay,Az,Ax_1D,Ay_1D,Az_1D,sgn,nrm)
!   A_1D = sqrt(Ax_1D**2+Ay_1D**2+Az_1D**2)
!   call save_value_orbit(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss, &
!   		A_1D,Ax_1D,Ay_1D,Az_1D,planet,coord,angle)
   		

!   !--Write tree.xml metadata
!   call write_numdat_XML(fildat,planet,prefix(1:3),.true.,.false.,.false.,.false.,.false.,trim(spacecraft),angle,&
!	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
!   call write_granule_TS_XML(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss,numdatid,angle)
!
!
!    !--3DCubes for electron density
!    prefix = "Thew"
!   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.true.,.false.,.false.,.false.,'',0.,&
!	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
!   call write_granule_3D_XML(prefix,run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)


    !-- Read Elew file
    prefix = "Elew"
    var_name1 = "Ex"
    var_name2 = "Ey"
    var_name3 = "Ez"
    
   !--Write tree.xml metadata
   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.true.,.false.,.false.,.false.,'',0.,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
   call write_granule_3D_XML(prefix,run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)

    Ax = zero;	Ay = zero;	Az = zero
    Ax_1D = zero;	Ay_1D = zero;	Az_1D = zero
    A_1D = zero
    
   call read_field_cdf(run_name,prefix,var_name1,var_name2,var_name3, Ax,Ay,Az,n0,V0,B0)
    nrm=1 ! eletric field is already  in mV/km   
   call value_along_orbit(X_traj,Y_traj,Z_traj,gstep,Ax,Ay,Az,Ax_1D,Ay_1D,Az_1D,sgn,nrm)
   A_1D = sqrt(Ax_1D**2+Ay_1D**2+Az_1D**2)
   call save_value_orbit(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss, &
   		A_1D,Ax_1D,Ay_1D,Az_1D,planet,coord,angle)   		
   		
   !--Write tree.xml metadata
   call write_numdat_XML(fildat,planet,prefix(1:3),.true.,.false.,.false.,.false.,.false.,trim(spacecraft),angle,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
   call write_granule_TS_XML(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss,numdatid,angle)
   !-- Read Denw file
   ! prefix = "Denw"
   ! var_name1 = "Dn_tot"
   ! A = zero
   ! A_1D = zero
   ! call read_density_cdf(run_name,prefix,var_name1, A)
   ! call density_along_orbit(X_traj,Y_traj,Z_traj,gstep,A,A_1D)
   ! call save_density_orbit(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss, &
   !		A_1D)      
   
   !-- Read Ion species information
   iunit = 10
   write(write_name,'(a,a1,2a)')"Read_moment_species","_",trim(run_name),".dat"
   open(unit = iunit,file = write_name,form = 'formatted',status  = 'old',action = 'read')
   read(iunit,*) nb_spe
   allocate(ion_label(nb_spe)); ion_label(:) = ''
   allocate(diag_name(nb_spe));	diag_name(:) =' '
   read(iunit,*) ion_label
   do ii = 1,nb_spe
     read(iunit,*) diag_name(ii)
   !  print *,diag_name(ii),'  ',ion_label(ii)
   enddo
   close(iunit)  
   
   !-- Read Ion species files
   var_name1 = "Density"
   var_name2 = "Ux"
   var_name3 = "Uy"
   var_name4 = "Uz"
   var_name5 = "Temperature"
   do ii = 1,nb_spe
     A = zero;	B = zero
     Ax = zero;	Ay = zero;	Az = zero
     A_1D = zero;	B_1D = zero
     Ax_1D = zero;	Ay_1D = zero;	Az_1D = zero
     call read_species_cdf(diag_name(ii),var_name1,var_name2,var_name3,var_name4,var_name5,A, Ax,Ay,Az,B)
     nrm_n =1 !density is already  in cm-3
     nrm_U =1     ! speed is already in km/s
     nrm_T = 1 ! temperature is already in eV
     select case (trim(planet))
     case("ganymede","titan")
       nrm_T = 1. ! to take into account that normalisation is O+ mass
     end select
     
     call value_species_along_orbit(X_traj,Y_traj,Z_traj,gstep,A,Ax,Ay,Az,B,A_1D,Ax_1D,Ay_1D,Az_1D,B_1D, &
     	sgn,nrm_n,nrm_U,nrm_T)
    prefix = trim(ion_label(ii)) 	
    call save_value_species_orbit(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss, &
   		A_1D,Ax_1D,Ay_1D,Az_1D,B_1D,planet,coord,angle)    	
     
   !--Write tree.xml metadata
   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.true.,.false.,.false.,.false.,'',0.,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
!	minval(X_traj),maxval(X_traj),minval(Y_traj),maxval(Y_traj),minval(Z_traj),minval(Z_traj))   
   call write_granule_3D_XML(prefix,run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
   !--Write tree.xml metadata
   print *,'Spacecraft name :',spacecraft
   call write_numdat_XML(fildat,planet,prefix(1:3),.true.,.false.,.false.,.false.,.false.,trim(spacecraft),angle,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
   call write_granule_TS_XML(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss,numdatid,angle)
   enddo
   
   
   deallocate(A,Ax,Ay,Az,B)
   deallocate(A_1D,Ax_1D,Ay_1D,Az_1D,B_1D)


   end subroutine extract_1D
   
   !!##############################################################
   !! m_IMPEX_diag/read_dim_field_cdf
   !! This routines reads a simulation file
   !! and extract the dimension of a given variabe
   !! associated to a tag
   subroutine read_dim_field_cdf(run_name,prefix,varname1,ncm_tot,radius,centr,gstep,planet,c_wpi)
        
   character(len=*),intent(in) :: run_name
   character(len=*),intent(in) :: prefix
   character(len=*),intent(in) :: varname1
   integer,intent(inout)         :: ncm_tot(3)
   real(dp),intent(out)        :: radius,centr(3),gstep(3),c_wpi
   character(len=8),intent(out) :: planet
   
   !--Others
   integer  :: stId,ncid
   logical :: file_e
   character(len=50) :: write_name
   character(len=64) :: filename
   integer :: varfieldId, numDims,var_id
   integer,dimension(nf90_max_var_dims) :: DimfieldId
       
   !--Create file name for Field 
   write(filename,'(a4,a1,2a)')trim(prefix),"_",trim(run_name),".nc"
   write(*,*) 'Open file'  ,filename  
     
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
                
         
   !--Get number of point for fields in X,Y,Z
   StId = nf90_inq_varId(ncid,varname1,varfieldId)
   StId = nf90_Inquire_Variable(ncid,varfieldId,dimids = DimfieldId)
   StId = nf90_Inquire_Dimension(ncid,DimfieldId(1), len=ncm_tot(1))
   StId = nf90_Inquire_Dimension(ncid,DimfieldId(2), len=ncm_tot(2))
   StId = nf90_Inquire_Dimension(ncid,DimfieldId(3), len=ncm_tot(3))
   
    !--Get the grid step
    call get_simple_variable_cdf(ncid,"gstep",gstep)

   !--Get the obstacle radius
   call get_simple_variable_cdf(ncid,"r_planet",radius)
   
   !--Get the obstacle position
   call get_simple_variable_cdf(ncid,"s_centr",centr)  
   
    !--Get the refernce length (ion inertial length)
    call get_simple_variable_cdf(ncid,"phys_length",c_wpi)  

   !--Get Planetname
   stId = nf90_inq_varid(ncid,"planetname",var_id)
   stId = nf90_get_var(ncid,var_id,planet)
   write(*,*) 'Planet :',planet
   

    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename   
   
   
   end subroutine read_dim_field_cdf
   
   
   !!############################################################
   !! m_IMPEX_diag/read_field_cdf
   !! This routines reads the global output
   !! simulation files. Concern files are :
   !!		- Magw
   !!		- Elew
   !!		- Velw    
    subroutine read_field_cdf(run_name,prefix,varname1,varname2,varname3,Ax,Ay,Az,n0,v0,B0)
        
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: prefix
    character(len=*),intent(in) :: varname1,varname2,varname3
    real(dp),dimension(:,:,:),intent(inout) :: Ax,Ay,Az
    real(dp),intent(inout) :: n0,B0,V0
        
    !--Local variables
    character(len=64) :: varname    
    integer  :: stId,ncid
    logical :: file_e
    character(len=50) :: write_name
    character(len=64) :: filename

    
    !--Create file name for Field 
    write(filename,'(a4,a1,2a)')trim(prefix),"_",trim(run_name),".nc"
    write(*,*) 'Open file'  ,filename  
  
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
             
    
    call get_simple_variable_cdf(ncid,varname1 ,Ax(:,:,:) )
    call get_simple_variable_cdf(ncid,varname2 ,Ay(:,:,:) )
    call get_simple_variable_cdf(ncid,varname3 ,Az(:,:,:) )  

   !--Get normalisation values
   call get_simple_variable_cdf(ncid,"phys_mag",B0)
   call get_simple_variable_cdf(ncid,"phys_density",n0)
   call get_simple_variable_cdf(ncid,"phys_speed",V0)    
     
    
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename
    
    end subroutine read_field_cdf
    
!!############################################################
   !! m_IMPEX_diag/read_species_cdf
   !! This routines reads the global output
   !! of ion species simulation files. Concern the quantities are :
   !!		- Density
   !!		- Vx
   !!		- Vy
   !!		- Vz
   !!		- Temperature
    subroutine read_species_cdf(run_name,varname1,varname2,varname3,varname4,varname5,A,Ax,Ay,Az,B)
        
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: varname1,varname2,varname3,varname4,varname5
    real(dp),dimension(:,:,:),intent(inout) :: Ax,Ay,Az,A,B
        
    !--Local variables
    character(len=64) :: varname    
    integer  :: stId,ncid
    logical :: file_e
    character(len=50) :: write_name
    character(len=64) :: filename

    
    !--Create file name for Field 
    !write(filename,'(a4,a1,2a)')trim(prefix),"_",trim(run_name),".nc"
    !write(*,*) 'Open file'  ,filename  
  
    !--Inquire if the file exists
    inquire( file=trim(run_name), exist=file_e )
    !--No file: return
    if(.not.(file_e)) then
     call wrtout(6,"File: "//trim(run_name)//" does not exits","PERS")
     return
    endif
  
    !--Open NetCDF file
    stId = nf90_open(trim(run_name), nf90_nowrite, ncid)
    call test_cdf(stId)
             
    call get_simple_variable_cdf(ncid,varname1 ,A(:,:,:) )
    call get_simple_variable_cdf(ncid,varname2 ,Ax(:,:,:) )
    call get_simple_variable_cdf(ncid,varname3 ,Ay(:,:,:) )
    call get_simple_variable_cdf(ncid,varname4 ,Az(:,:,:) )  
    call get_simple_variable_cdf(ncid,varname5 ,B(:,:,:) )
     
    
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',run_name
    
    end subroutine read_species_cdf  
    
!!############################################################
   !! m_IMPEX_diag/read_density_cdf
   !! This routines reads the global output
   !! simulation files. Concern files are :
   !!		- Denw
    subroutine read_density_cdf(run_name,prefix,varname,A)
        
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: prefix
    character(len=*),intent(in) :: varname
    real(dp),dimension(:,:,:),intent(inout) :: A
        
    !--Local variables
    integer  :: stId,ncid
    logical :: file_e
    character(len=50) :: write_name
    character(len=64) :: filename

    
    !--Create file name for Field 
    write(filename,'(a4,a1,2a)')trim(prefix),"_",trim(run_name),".nc"
    write(*,*) 'Open file'  ,filename  
  
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
             
    
    call get_simple_variable_cdf(ncid,varname ,A(:,:,:) )
     
    
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename
    
    end subroutine read_density_cdf    
    
    !!#############################################################
    !! m_IMPEX_diag/value_along_orbit
    !! This routine extracts value from the simulation
    !! along a trajectory
    subroutine value_along_orbit(X_traj,Y_traj,Z_traj,gstep,Ax,Ay,Az,Ax_1D,Ay_1D,Az_1D,sgn,nrm)
    real(dp),dimension(:),intent(in) :: X_traj,Y_traj,Z_traj
    real(dp),dimension(:),intent(inout) :: Ax_1D,Ay_1D,Az_1D
    real(dp),dimension(:,:,:),intent(in) :: Ax,Ay,Az
    real(dp) :: gstep(3),nrm
    integer,intent(in) :: sgn
    !Local variables
    integer :: n,nb_traj
    integer :: i,j,k
    real(dp) :: w1,w2,w3,w4,w5,w6,w7,w8
    real(dp) :: xa,xf,ya,yf,za,zf,x,y,z
    integer::ncm(3)
    
    !-- Get array dimension
    ncm(1) = size(Ax,1)
    ncm(2) = size(Ay,2)
    ncm(3) = size(Az,3)
    
    nb_traj = size(X_traj,1)
    do n=1,nb_traj
      i = int(X_traj(n)/gstep(1))+1
      j = int(Y_traj(n)/gstep(2))+1
      k = int(Z_traj(n)/gstep(3))+1
      
      
      if ((i>0) .and. (i<=ncm(1))) then
        if ((j>0) .and. (j<=ncm(2))) then
          if ((k>0) .and. (k<=ncm(3))) then
      
      	xa = X_traj(n)/gstep(1)+1 - float(i)
      	ya = Y_traj(n)/gstep(2)+1 - float(j)
      	za = Z_traj(n)/gstep(3)+1 - float(k)
      	xf = 1.0 - xa
      	yf = 1.0 - ya
      	zf = 1.0 - za
      	w1 = xa*ya*za   ! (i+1,j+1,k+1)
      	w2 = xf*ya*za   ! (i  ,j+1,k+1)
      	w3 = xa*yf*za   ! (i+1,j  ,k+1)
      	w4 = xf*yf*za   ! (i  ,j  ,k+1)
      	w5 = xa*ya*zf   ! (i+1,j+1,k  )
      	w6 = xf*ya*zf   ! (i  ,j+1,k  )
      	w7 = xa*yf*zf   ! (i+1,j  ,k  )
      	w8 = xf*yf*zf   ! (i  ,j  ,k  )
      	
      	Ax_1D(n) = ( w1*Ax(i+1,j+1,k+1) + &
      	     w2*Ax(i  ,j+1,k+1) + w3*Ax(i+1,j  ,k+1) + &
      	     w4*Ax(i  ,j  ,k+1) + w5*Ax(i+1,j+1,k  ) + &
      	     w6*Ax(i  ,j+1,k  ) + w7*Ax(i+1,j  ,k  ) + &
      	     w8*Ax(i  ,j  ,k  ))
      	Ay_1D(n) = ( w1*Ay(i+1,j+1,k+1) + &
      	     w2*Ay(i  ,j+1,k+1) + w3*Ay(i+1,j  ,k+1) + &
      	     w4*Ay(i  ,j  ,k+1) + w5*Ay(i+1,j+1,k  ) + &
      	     w6*Ay(i  ,j+1,k  ) + w7*Ay(i+1,j  ,k  ) + &
      	     w8*Ay(i  ,j  ,k  ))
      	Az_1D(n) = ( w1*Az(i+1,j+1,k+1) + &
      	     w2*Az(i  ,j+1,k+1) + w3*Az(i+1,j  ,k+1) + &
      	     w4*Az(i  ,j  ,k+1) + w5*Az(i+1,j+1,k  ) + &
      	     w6*Az(i  ,j+1,k  ) + w7*Az(i+1,j  ,k  ) + &
      	     w8*Az(i  ,j  ,k  ))
      	  endif
      	endif
     endif 	
      
    enddo
    !--Convert in MSO
    Ax_1D = Ax_1D
    Ay_1D = Ay_1D
    Az_1D =  Az_1D
    
    
    end subroutine value_along_orbit
    
    !!#############################################################
        !! m_IMPEX_diag/density_along_orbit
        !! This routine extracts electron density value from the simulation
        !! along a trajectory
        subroutine density_along_orbit(X_traj,Y_traj,Z_traj,gstep,A,A_1D)
        real(dp),dimension(:),intent(in) :: X_traj,Y_traj,Z_traj
        real(dp),dimension(:),intent(inout) :: A_1D
        real(dp),dimension(:,:,:),intent(in) :: A
        real(dp) :: gstep(3)
        !Local variables
        integer :: n,nb_traj
        integer :: i,j,k
        real(dp) :: w1,w2,w3,w4,w5,w6,w7,w8
        real(dp) :: xa,xf,ya,yf,za,zf,x,y,z
        integer::ncm(3)
        
        !-- Get array dimension
        ncm(1) = size(A,1)
        ncm(2) = size(A,2)
        ncm(3) = size(A,3)
        
        nb_traj = size(X_traj,1)
        do n=1,nb_traj
          i = int(X_traj(n)/gstep(1))
          j = int(Y_traj(n)/gstep(2))
          k = int(Z_traj(n)/gstep(3))
          
          
          if ((i>0) .and. (i<=ncm(1))) then
            if ((j>0) .and. (j<=ncm(2))) then
              if ((k>0) .and. (k<=ncm(3))) then
          
          	xa = X_traj(n)/gstep(1) - float(i)
          	ya = Y_traj(n)/gstep(2) - float(j)
          	za = Z_traj(n)/gstep(3) - float(k)
          	xf = 1.0 - xa
          	yf = 1.0 - ya
          	zf = 1.0 - za
          	w1 = xa*ya*za   ! (i+1,j+1,k+1)
          	w2 = xf*ya*za   ! (i  ,j+1,k+1)
          	w3 = xa*yf*za   ! (i+1,j  ,k+1)
          	w4 = xf*yf*za   ! (i  ,j  ,k+1)
          	w5 = xa*ya*zf   ! (i+1,j+1,k  )
          	w6 = xf*ya*zf   ! (i  ,j+1,k  )
          	w7 = xa*yf*zf   ! (i+1,j  ,k  )
          	w8 = xf*yf*zf   ! (i  ,j  ,k  )
          	
          	A_1D(n) = ( w1*A(i+1,j+1,k+1) + &
          	     w2*A(i  ,j+1,k+1) + w3*A(i+1,j  ,k+1) + &
          	     w4*A(i  ,j  ,k+1) + w5*A(i+1,j+1,k  ) + &
          	     w6*A(i  ,j+1,k  ) + w7*A(i+1,j  ,k  ) + &
          	     w8*A(i  ,j  ,k  ))
          	  endif
          	endif
         endif 	
          
        enddo
        
        
    end subroutine density_along_orbit

    !!#############################################################
    !! m_IMPEX_diag/value_species_along_orbit
    !! This routine extracts value from the simulation
    !! along a trajectory
    subroutine value_species_along_orbit(X_traj,Y_traj,Z_traj,gstep,A,Ax,Ay,Az,B,A_1D,Ax_1D,Ay_1D,Az_1D,B_1D,sgn,nrm_n,nrm_v,nrm_t)
    real(dp),dimension(:),intent(in) :: X_traj,Y_traj,Z_traj
    real(dp),dimension(:),intent(inout) :: Ax_1D,Ay_1D,Az_1D,A_1D,B_1D
    real(dp),dimension(:,:,:),intent(in) :: Ax,Ay,Az,A,B
    real(dp) :: gstep(3),nrm_n,nrm_v,nrm_t
    integer,intent(in) :: sgn
    !Local variables
    integer :: n,nb_traj
    integer :: i,j,k
    real(dp) :: w1,w2,w3,w4,w5,w6,w7,w8
    real(dp) :: xa,xf,ya,yf,za,zf,x,y,z
    integer::ncm(3)
    
    !-- Get array dimension
    ncm(1) = size(Ax,1)
    ncm(2) = size(Ay,2)
    ncm(3) = size(Az,3)
    
    nb_traj = size(X_traj,1)
    do n=1,nb_traj
      i = int(X_traj(n)/gstep(1))+1
      j = int(Y_traj(n)/gstep(2))+1
      k = int(Z_traj(n)/gstep(3))+1
      
      
      if ((i>0) .and. (i<=ncm(1))) then
        if ((j>0) .and. (j<=ncm(2))) then
          if ((k>0) .and. (k<=ncm(3))) then
      
      	xa = X_traj(n)/gstep(1)+1 - float(i)
      	ya = Y_traj(n)/gstep(2)+1 - float(j)
      	za = Z_traj(n)/gstep(3)+1 - float(k)
      	xf = 1.0 - xa
      	yf = 1.0 - ya
      	zf = 1.0 - za
      	w1 = xa*ya*za   ! (i+1,j+1,k+1)
      	w2 = xf*ya*za   ! (i  ,j+1,k+1)
      	w3 = xa*yf*za   ! (i+1,j  ,k+1)
      	w4 = xf*yf*za   ! (i  ,j  ,k+1)
      	w5 = xa*ya*zf   ! (i+1,j+1,k  )
      	w6 = xf*ya*zf   ! (i  ,j+1,k  )
      	w7 = xa*yf*zf   ! (i+1,j  ,k  )
      	w8 = xf*yf*zf   ! (i  ,j  ,k  )
      	A_1D(n) = ( w1*A(i+1,j+1,k+1) + &
      	     w2*A(i  ,j+1,k+1) + w3*A(i+1,j  ,k+1) + &
      	     w4*A(i  ,j  ,k+1) + w5*A(i+1,j+1,k  ) + &
      	     w6*A(i  ,j+1,k  ) + w7*A(i+1,j  ,k  ) + &
      	     w8*A(i  ,j  ,k  ))      	
      	Ax_1D(n) = ( w1*Ax(i+1,j+1,k+1) + &
      	     w2*Ax(i  ,j+1,k+1) + w3*Ax(i+1,j  ,k+1) + &
      	     w4*Ax(i  ,j  ,k+1) + w5*Ax(i+1,j+1,k  ) + &
      	     w6*Ax(i  ,j+1,k  ) + w7*Ax(i+1,j  ,k  ) + &
      	     w8*Ax(i  ,j  ,k  ))
      	Ay_1D(n) = ( w1*Ay(i+1,j+1,k+1) + &
      	     w2*Ay(i  ,j+1,k+1) + w3*Ay(i+1,j  ,k+1) + &
      	     w4*Ay(i  ,j  ,k+1) + w5*Ay(i+1,j+1,k  ) + &
      	     w6*Ay(i  ,j+1,k  ) + w7*Ay(i+1,j  ,k  ) + &
      	     w8*Ay(i  ,j  ,k  ))
      	Az_1D(n) = ( w1*Az(i+1,j+1,k+1) + &
      	     w2*Az(i  ,j+1,k+1) + w3*Az(i+1,j  ,k+1) + &
      	     w4*Az(i  ,j  ,k+1) + w5*Az(i+1,j+1,k  ) + &
      	     w6*Az(i  ,j+1,k  ) + w7*Az(i+1,j  ,k  ) + &
      	     w8*Az(i  ,j  ,k  ))
	B_1D(n) = ( w1*B(i+1,j+1,k+1) + &
      	     w2*B(i  ,j+1,k+1) + w3*B(i+1,j  ,k+1) + &
      	     w4*B(i  ,j  ,k+1) + w5*B(i+1,j+1,k  ) + &
      	     w6*B(i  ,j+1,k  ) + w7*B(i+1,j  ,k  ) + &
      	     w8*B(i  ,j  ,k  ))      	
      	      	     
      	  endif
      	endif
     endif 	
      
    enddo
    !--Convert in good coordinate
    Ax_1D = Ax_1D
    Ay_1D = Ay_1D
    Az_1D =  Az_1D
    A_1D = A_1D
    B_1D = B_1D
    
    
    end subroutine value_species_along_orbit    
    
    !!###################################################
    !! m_IMPEX_diag/save_value_orbit
    !! This routines dumps into a file
    !! the extracted value along a trajectory
    subroutine save_value_orbit(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss, &
   		A_1D,Ax_1D,Ay_1D,Az_1D,planetname,coord,angle)
    character(len=*),intent(in) ::traj_name,prefix,run_name,planetname,coord
    integer,dimension(:),intent(in) :: YYYY,MM,DD,hh,minut
    real(dp),dimension(:),intent(in) :: ss,A_1D,Ax_1D,Ay_1D,Az_1D
    real(dp),intent(in) :: angle
    ! Local variables
    character(len=50) :: write_name
    character(len=23)::votime
    integer :: iunit,i,nb
    character(len=5) ::msg_angle

    if (angle /= -1) then
    write(msg_angle,'(f5.1)')angle
    else
    msg_angle =""
    endif

    
    !--Create file name for orbital value 
    write(write_name,'(a3,a1,a,a1,a,a,a,a)')trim(prefix),"_",trim(run_name(1:len_trim(run_name)-7)),"_",trim(traj_name(1:len_trim(traj_name)-4)),"_",trim(adjustl(msg_angle)),".xml"
    write(*,*) 'Saving file'  ,write_name    
    
    !-- Get length of the array
    nb = size(A_1D,1)
    ! in VOTale format
    iunit = 1
     open(UNIT = iunit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
    	ACTION = 'WRITE')
    call header_orbit_VOTABLE(iunit,prefix,planetname,coord)	
    
     do i=1,nb
          write(iunit,'(a)') '<TR>' 
    call get_VOTIME(YYYY(i),MM(i),DD(i),hh(i),minut(i),ss(i),votime)
          write(iunit,'(2a,4(a,e12.3),a)') '<TD>',votime,'</TD> <TD>',A_1D(i),'</TD> <TD>',Ax_1D(i),'</TD> <TD>',Ay_1D(i),'</TD> <TD>',Az_1D(i),'</TD>'
          write(iunit,'(a)') '</TR>'	
    enddo
  ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESOURCE>'
  write(iunit,'(a)') '</VOTABLE>'      
    close(iunit)
    
    ! in ascii format
!    write(write_name,'(a3,a1,a,a1,a)')trim(prefix),"_",trim(run_name),"_",trim(traj_name)
!    iunit = 1
!    open(UNIT = iunit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
!    	ACTION = 'WRITE')
!    	
!        select case(trim(prefix))
!    case("Magw")
!      write(iunit,'(a)') 'YYYY|MM|DD|hh|mm|ss|  Btot[nT]  |  Bx[nT]  |  By[nT]  |  Bz[nT]' 	
!    case("Elew")
!      write(iunit,'(a)') 'YYYY|MM|DD|hh|mm|ss|  Etot[mV/km]  |  Ex[mV/km]  |  Ey[mV/km]  |  Ez[mV/km]'
!    case("Velw")
!      write(iunit,'(a)') 'YYYY|MM|DD|hh|mm|ss|  Z [RM] |  Utot[km/s]  |  Ux[km/s]  |  Uy[km/s]  |  Uz[km/s]'
!    end select
     	    	
!    do i=1,nb
!      write(iunit,'(i4,4(a1,i2),5(a1,f8.3))') YYYY(i)," ",MM(i)," ",DD(i)," ",hh(i)," ",minut(i), &
!      	" ",ss(i)," ",A_1D(i)," ",Ax_1D(i)," ",Ay_1D(i)," ",Az_1D(i)
!    enddo	
!    close(iunit)
    end subroutine save_value_orbit
    
    !!###################################################
    !! m_IMPEX_diag/save_density_orbit
    !! This routines dumps into a file
    !! the extracted value along a trajectory
    subroutine save_density_orbit(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss, &
   		A_1D)
    character(len=*),intent(in) ::traj_name,prefix,run_name
    integer,dimension(:),intent(in) :: YYYY,MM,DD,hh,minut
    real(dp),dimension(:),intent(in) :: ss,A_1D
    ! Local variables
    character(len=50) :: write_name
    integer :: iunit,i,nb
    
    !--Create file name for orbital value 
    write(write_name,'(a3,a1,a,a1,a)')trim(prefix),"_",trim(run_name),"_",trim(traj_name)
    write(*,*) 'Saving file'  ,write_name     
    
    !-- Get length of the array
    nb = size(A_1D,1)
    
    iunit = 1
    open(UNIT = iunit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
    	ACTION = 'WRITE')
    do i=1,nb
      write(iunit,'(i4,4(a1,i2),2(a1,f8.3))') YYYY(i)," ",MM(i)," ",DD(i)," ",hh(i)," ",minut(i), &
      	" ",ss(i)," ",A_1D(i)
    enddo	
    close(iunit)
    end subroutine save_density_orbit  
    
    !!###################################################
    !! m_IMPEX_diag/save_value_species_orbit
    !! This routines dumps into a file
    !! the extracted value along a trajectory
    subroutine save_value_species_orbit(traj_name,prefix,run_name,YYYY,MM,DD,hh,minut,ss, &
   		A_1D,Ax_1D,Ay_1D,Az_1D,B_1D,planetname,coord,angle)
    character(len=*),intent(in) ::traj_name,prefix,run_name,planetname,coord
    integer,dimension(:),intent(in) :: YYYY,MM,DD,hh,minut
    real(dp),dimension(:),intent(in) :: ss,A_1D,Ax_1D,Ay_1D,Az_1D,B_1D
    real(dp),intent(in) :: angle
    ! Local variables
    character(len=50) :: write_name
    character(len=23)::votime
    integer :: iunit,i,nb
    character(len=5) ::msg_angle
    real(dp) :: vect_norm

    if (angle /= -1) then
    write(msg_angle,'(f5.1)')angle
    else
    msg_angle =""
    endif
   
    !-- Create xml file
    !--Create file name for orbital value 
    write(write_name,'(a3,a1,a,a1,a,a,a,a)')trim(prefix),"_",trim(run_name(1:len_trim(run_name)-7)),"_",trim(traj_name(1:len_trim(traj_name)-4)),"_",trim(adjustl(msg_angle)),".xml"
    write(*,*) 'Saving file'  ,write_name
    
    !-- Get length of the array
    nb = size(A_1D,1)
    ! in VOTale format
    iunit = 1
     open(UNIT = iunit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
    	ACTION = 'WRITE')
    call header_species_orbit_VOTABLE(iunit,prefix,planetname,coord)	
    
     do i=1,nb
          write(iunit,'(a)') '<TR>' 
    call get_VOTIME(YYYY(i),MM(i),DD(i),hh(i),minut(i),ss(i),votime)
    vect_norm=sqrt(Ax_1D(i)**2+Ax_1D(i)**2+Ax_1D(i)**2)
          write(iunit,'(2a,6(a,e12.3),a)') '<TD>',votime,'</TD> <TD>',A_1D(i),'</TD><TD>',vect_norm,'</TD><TD>',Ax_1D(i),'</TD> <TD>',Ay_1D(i),'</TD> <TD>',Az_1D(i),'</TD> <TD>',B_1D(i),'</TD>'
          write(iunit,'(a)') '</TR>'	
    enddo
  ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESOURCE>'
  write(iunit,'(a)') '</VOTABLE>'      
    close(iunit)    
    
    
    !--Create file name for orbital value 
!    write(write_name,'(a3,a1,a,a1,a)')trim(prefix),"_",trim(run_name),"_",trim(traj_name)
!    write(*,*) 'Saving file'  ,write_name     
!    
!    !-- Get length of the array
!    nb = size(A_1D,1)
!    
!    iunit = 1
!    open(UNIT = iunit,FILE = write_name, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
!    	ACTION = 'WRITE')
!     write(iunit,'(a)') 'YYYY|MM|DD|hh|mm|ss| Dn[cm-3]| Ux[km/s]|  Uy[km/s] | Uz[km/s] | T[eV]'      	
!   do i=1,nb
!      write(iunit,'(i4,4(a1,i2),6(a1,f10.3))') YYYY(i)," ",MM(i)," ",DD(i)," ",hh(i)," ",minut(i), &
!      	" ",ss(i)," ",A_1D(i)," ",Ax_1D(i)," ",Ay_1D(i)," ",Az_1D(i)," ",B_1D(i)
!    enddo	
!    close(iunit)
    end subroutine save_value_species_orbit  
    
    !!######################################################
    !! m_IMPEX_diag/header_orbit_VOTABLE
    !! This routine write the header of the 
    !! VOTABLE file
    subroutine header_orbit_VOTABLE(iunit,prefix,planetname,coord)
    character(len=*),intent(in) :: prefix,planetname,coord
    integer,intent(in) :: iunit
  write(iunit,'(a)') '<?xml version="1.0"?>'
  write(iunit,'(a)') '<VOTABLE version="1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
  write(iunit,'(a)')  'xmlns="http://www.ivoa.net/xml/VOTable/v1.2" xmlns:stc="http://www.ivoa.net/xml/STC/v1.30" >'  
  write(iunit,'(a)')  '<RESOURCE name="IMPEx fly-through">'
  write(iunit,'(a)')  '<TABLE name="results">'
  
  select case (trim(planetname))
  case("mars","mercury")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="MSO">'
   write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="3393." unit="km" ucd="phys.size"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col1"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col3"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'
  case("ganymede")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="GPHIO">'
   write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="2634." unit="km" ucd="phys.size"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col1"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col3"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'  
  case("titan")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="MSO">'
   write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="2575." unit="km" ucd="phys.size"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col1"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col3"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'   
  end select 
  
  
  select case (trim(prefix))
  case("Magw")
    write(iunit,'(3a)') '<DESCRIPTION> Magnetic field components </DESCRIPTION>'
  case("Elew")
    write(iunit,'(3a)') '<DESCRIPTION> Electric field components </DESCRIPTION>'
  case("Velw")    
    write(iunit,'(a)') '<DESCRIPTION> Velocity components </DESCRIPTION>'
  end select

    write(iunit,'(3a)') '<FIELD name="Time" ID="col1" ucd="time.epoch" ', &
    	'  xtype="dateTime" utype="" datatype="char" arraysize="*" />' 

  select case (trim(prefix))
  case("Magw")
    write(iunit,'(3a)') '<FIELD name="Btot" ID="col2" ucd="phys.magField" ',&
        '  utype="" datatype="float" width="10" unit="nT" />'        
    write(iunit,'(3a)') '<FIELD name="Bx" ID="col3" ucd="phys.magField" ',&
    	'  utype="" datatype="float" width="10" unit="nT" />'
    write(iunit,'(3a)') '<FIELD name="By" ID="col4" ucd="phys.magField" ',&
        '  utype="" datatype="float" width="10" unit="nT" />'
    write(iunit,'(3a)') '<FIELD name="Bz" ID="col5" ucd="phys.magField" ',&
        '  utype="" datatype="float" width="10" unit="nT" />'
  case("Elew")
    write(iunit,'(3a)') '<FIELD name="Etot" ID="col2" ucd="phys.elecField" ',&
        '  utype="" datatype="float" width="10" unit="mV.m-1" />' 
    write(iunit,'(3a)') '<FIELD name="Ex" ID="col3" ucd="phys.elecField" ',&
    	'  utype="" datatype="float" width="10" unit="mV.m-1" />'
    write(iunit,'(3a)') '<FIELD name="Ey" ID="col4" ucd="phys.elecField" ',&
        '  utype="" datatype="float" width="10" unit="mV.m-1" />'
    write(iunit,'(3a)') '<FIELD name="Ez" ID="col5" ucd="phys.elecField" ',&
        '  utype="" datatype="float" width="10" unit="mV.m-1" />'
  case("Velw")    
    write(iunit,'(3a)') '<FIELD name="Utot" ID="col2" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="10" unit="km.s-1" />' 
    write(iunit,'(3a)') '<FIELD name="Ux" ID="col3" ucd="phys.veloc" ',&
    	'  utype="" datatype="float" width="10" unit="km.s-1" />'
    write(iunit,'(3a)') '<FIELD name="Uy" ID="col4" ucd="phys.veloc" ',&
        ' " utype="" datatype="float" width="10" unit="km.s-1" />'
    write(iunit,'(3a)') '<FIELD name="Uz" ID="col5" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="10" unit="km.s-1" />'
  end select
  
  
  write(iunit,'(a)') '<DATA>'
  write(iunit,'(a)') '<TABLEDATA>'  
    end subroutine header_orbit_VOTABLE
  !!######################################################
    !! m_IMPEX_diag/header_species_orbit_VOTABLE
    !! This routine write the header of the 
    !! VOTABLE file
    subroutine header_species_orbit_VOTABLE(iunit,prefix,planetname,coord)
    character(len=*),intent(in) :: prefix,planetname,coord
    integer,intent(in) :: iunit

 
  write(iunit,'(a)') '<?xml version="1.0"?>'
  write(iunit,'(a)') '<VOTABLE version="1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
  write(iunit,'(a)')  'xmlns="http://www.ivoa.net/xml/VOTable/v1.2" xmlns:stc="http://www.ivoa.net/xml/STC/v1.30" >'  
  write(iunit,'(a)')  '<RESOURCE name="IMPExTimeTable">'
  write(iunit,'(a)')  '<TABLE name="results">'
  
  select case (trim(planetname))
  case("mars","mercury")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="MSO">'
   write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="3393." unit="km" ucd="phys.size"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col1"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col3"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'
  case("ganymede")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="GPHIO">'
   write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="2634." unit="km" ucd="phys.size"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col1"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col3"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'  
  case("titan")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="MSO">'
   write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="2575." unit="km" ucd="phys.size"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col1"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col3"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'   
  end select   
  
   write(iunit,'(3a)') '<DESCRIPTION> ',trim(prefix),' Density, Velocity, Temperature </DESCRIPTION>'
    write(iunit,'(3a)') '<FIELD name="Time" ID="col1" ucd="time.epoch" ', &
    	'  xtype="dateTime" utype="" datatype="char" arraysize="*" />' 
    write(iunit,'(3a)') '<FIELD name="Density" ID="col2" ucd="phys.density" ',&
    	'  utype="" datatype="float" width="12" unit="cm-3"/>'
    write(iunit,'(3a)') '<FIELD name="Utot" ID="col6" ucd="phys.veloc" ',&
        '  utype="" dataype="float" width="12" unit="km.s-1"/>'         
    write(iunit,'(3a)') '<FIELD name="Ux" ID="col3" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1"/>'
    write(iunit,'(3a)') '<FIELD name="Uy" ID="col4" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1"/>'
    write(iunit,'(3a)') '<FIELD name="Uz" ID="col5" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1"/>'  
    write(iunit,'(3a)') '<FIELD name="Temperature" ID="co6" ucd="phys.temperature" ',&
        '  utype="" datatype="float" width="12" unit="eV"/>'        
  
  write(iunit,'(a)') '<DATA>'
  write(iunit,'(a)') '<TABLEDATA>'    
  
    end subroutine header_species_orbit_VOTABLE
    
    
  #endif  
 end module m_IMPEX_moment_orbit
