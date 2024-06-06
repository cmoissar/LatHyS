!=============================================================
!=============================================================
module m_IMPEX_cut

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
      read_field_cdf,	  &
      read_species_cdf,	  &
      extract_XY,	  &
      extract_XZ,	  &
      extract_YZ,	   &
      extract_species_XY,  &
      extract_species_XZ,  &
      header_2Dcut_VOTABLE, &
      header_2Dcut_J_VOTABLE, &
      header_2Dcut_species_VOTABLE
      !create_file_diag,   &
      !read_moment_species_cdf,&
      !read_field_cdf_1_array

 public ::                &
      extract_2D
contains
 !********************************************************************
 ! Auteur				:	 RModolo, SHess
 ! Date					:	 04/04/12
 ! Institution				:	LATMOS/CNRS/IPSL
 ! Derniere modification		:	05/04/12	
 ! Resume	
 ! 
 !********************************************************************
 
 
   !!############################################################
   !!=============================================================
   !!subroutine: m_IMPEX_cut/Extract_2D
   !! This routine reads merged files and 
   !! create 3 2D files according to XY, XZ and YZ planes
   
    subroutine extract_2D(run_name)
   
   character(len=*),intent(in) :: run_name
   !Local Variable
   character(len=64) :: var_name1,var_name2,var_name3,var_name4,var_name5
   character(len=4)  :: prefix
   character(len=10),dimension(:),allocatable:: ion_label
   character(len=50),dimension(:),allocatable:: diag_name
   integer :: nb_spe,ii,iunit
   character(len=50) :: write_name   
   
   prefix = "Magw"
   var_name1 = "Bx"
   var_name2 = "By"
   var_name3 = "Bz"
   
   call read_field_cdf(run_name,prefix, var_name1,var_name2, var_name3)


   prefix = "Elew"
   var_name1 = "Ex"
   var_name2 = "Ey"
   var_name3 = "Ez"
   
   call read_field_cdf(run_name,prefix, var_name1,var_name2, var_name3)  
   
!   prefix = "Velw"
!   var_name1 = "Vbulk_x"
!   var_name2 = "Vbulk_y"
!   var_name3 = "Vbulk_z"
   
!   call read_field_cdf(run_name,prefix, var_name1,var_name2, var_name3)   

   prefix = "Thew"
  !-- Read Ion species files
   var_name1 = "Density"
   var_name2 = "Ux"
   var_name3 = "Uy"
   var_name4 = "Uz"
   var_name5 = "Temperature" 
    write(write_name,'(a4,a1,2a)')trim(prefix),"_",trim(run_name),".nc"
    call read_species_cdf(write_name,run_name,'The',var_name1,var_name2,var_name3,var_name4,var_name5)



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
     call read_species_cdf(diag_name(ii),run_name,ion_label(ii),var_name1,var_name2,var_name3,var_name4,var_name5)
   enddo
    
    end subroutine extract_2D
    
    !!#############################################################
    !!=============================================================
    !! subroutine: m_IMPEX/read_field_cdf
    !! This routine reads the cdf global files
    
    subroutine read_field_cdf(run_name,prefix,varname1,varname2,varname3)
        
    character(len=*),intent(in) :: run_name
    character(len=*),intent(in) :: prefix
    character(len=*),intent(in) :: varname1,varname2,varname3
    !--Reading variables
    real(dp),dimension(:,:,:),allocatable :: Ax,Ay,Az,Amod,Jx,Jy,Jz,Jtot
    character(len=8) :: planet
    real(dp) :: sgn
    real(dp) :: B0,n0,V0,x0
        
    !--Others
    character(len=64) :: varname    
    integer  :: stId,ncid,var_id
    integer  :: varid(47),dimid(11),ns
    integer  :: ncm_tot(3),nLens,int_centr(3)
    real(dp) :: rtmp(3),radius,centr(3),gstep(3),dt,nrm,gstep_invq(3),nrm_J,box_size(6)
    logical :: file_e
    character(len=5) :: coord
    character(len=50) :: write_name
    character(len=64) :: filename
    character(len=100) :: test_name
    character(len=500) :: msg
    integer :: varfieldId,numDims,numatts,ival
    integer,dimension(nf90_max_var_dims) :: DimfieldId
   character(len=600)::numdatid

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
            
   !--Get Number of Species
   call get_simple_variable_cdf(ncid,"ns",ns)
      
   !--Allocate species and get informations
   ! call alloc_species(ns,Spe)
   ! call species_get_var_cdf(Spe,ncid)
   call get_simple_variable_cdf(ncid,"r_planet",radius)
   call get_simple_variable_cdf(ncid,"s_centr",centr)
   
   !--Get Planetname
   stId = nf90_inq_varid(ncid,"planetname",var_id)
   stId = nf90_get_var(ncid,var_id,planet)
   write(*,*) 'Planet :',trim(planet)
   !--Get normalisation values
   call get_simple_variable_cdf(ncid,"phys_mag",B0)
   call get_simple_variable_cdf(ncid,"phys_density",n0)
   call get_simple_variable_cdf(ncid,"phys_speed",V0)
   call get_simple_variable_cdf(ncid,"phys_length",x0)
   
   box_size=(/ -centr(1),ncm_tot(1)*gstep(1)-centr(1),-centr(2),ncm_tot(2)*gstep(2)-centr(2),-centr(3),ncm_tot(3)*gstep(3)-centr(3) /)
   box_size=box_size*x0

   !-- good normalization values
   select case (trim(prefix))
   case("Magw")
     nrm = 1! Mag is already normalized in nT
     nrm_J = 1./(x0*1.e3)/(four_pi*1.e-7)
   case("Velw")
     nrm = 1.
   case("Elew")
     nrm=1.! electric field is already in mV/m
   end select
    
    print *,'centr',centr
    int_centr = centr(:)/gstep(:)
    
        
    !--allocation des tableaux de lectures du fichier diag de champ
     allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
     allocate(Ay(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
     allocate(Az(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
     allocate(Amod(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
     Ax = zero ; Ay = zero ; Az = zero
     Amod = zero

     call get_simple_variable_cdf(ncid,varname1 ,Ax(:,:,:) )
     call get_simple_variable_cdf(ncid,varname2 ,Ay(:,:,:) )
     call get_simple_variable_cdf(ncid,varname3 ,Az(:,:,:) )  
     
     Amod = sqrt(Ax**2+Ay**2+Az**2)
     
    
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',filename
    
    ! coordinate transformation for each object
   ! planet="mars"
    
    ! all trasnformation coordinate have been done in the 3D cubes
    select case(trim(planet))
    case("mars")
      sgn = -1_dp
      coord='MSO'
    case("mercure")
      sgn = -1_dp
      coord='MSO'
    case ("ganymede")
      sgn=1_dp
      coord = 'GPhiO'
    case("titan")
      sgn =1_dp
      coord='TIIS'
    case default
     planet="mars    "
     sgn = -1_dp
      coord='MSO'
    end select
    
    print *,'sgn, coord',sgn,coord,planet
    
    varname = varname1(1:1)//'field'
    ival = int(int_centr(3))
    call extract_XY(Amod,Ax,Ay,Az,ival,centr,gstep,ncm_tot,varname,run_name,radius,prefix,planet,nrm,sgn,coord)
    ival = int(int_centr(2))
    call extract_XZ(Amod,Ax,Ay,Az,ival,centr,gstep,ncm_tot,varname,run_name,radius,prefix,planet,nrm,sgn,coord)   
    ival = int(int_centr(1)+2.5*radius/gstep(1))
    call extract_YZ(Amod,Ax,Ay,Az,ival,centr,gstep,ncm_tot,varname,run_name,radius,prefix,planet,nrm,sgn,coord)     

    
   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.false.,.true.,.false.,.false.,'XY',0.,&
	ncm_tot,ncm_tot,ncm_tot,ncm_tot,ncm_tot,gstep,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),(/ 0.,0.,1. /))
   call write_granule_2D_XML(varname,'XY_'//run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.false.,.true.,.false.,.false.,'YZ',0.,&
	ncm_tot,ncm_tot,ncm_tot,ncm_tot,ncm_tot,gstep,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),(/ 1.,0.,0. /))
   call write_granule_2D_XML(varname,'YZ_'//run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.false.,.true.,.false.,.false.,'XZ',0.,&
	ncm_tot,ncm_tot,ncm_tot,ncm_tot,ncm_tot,gstep,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),(/ 0.,1.,0. /))
   call write_granule_2D_XML(varname,'XZ_'//run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
    
    ! We calculate the total current
    if (trim(prefix) == "Magw") then
      gstep_invq = (quarter)/(gstep)
      allocate(Jx(ncm_tot(1)-1,ncm_tot(2)-1,ncm_tot(3)-1));	Jx(:,:,:) = zero
      allocate(Jy(ncm_tot(1)-1,ncm_tot(2)-1,ncm_tot(3)-1));	Jy(:,:,:) = zero
      allocate(Jz(ncm_tot(1)-1,ncm_tot(2)-1,ncm_tot(3)-1));	Jz(:,:,:) = zero
      allocate(Jtot(ncm_tot(1)-1,ncm_tot(2)-1,ncm_tot(3)-1));	Jtot(:,:,:) = zero
      
      Jx(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) = gstep_invq(2)*(&
             &   - Az(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   - Az(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   + Az(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   + Az(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   - Az(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   - Az(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   + Az(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
             &   + Az(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))&
             -gstep_invq(3)*(&                  
             &   - Ay(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   - Ay(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   - Ay(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   - Ay(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   + Ay(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   + Ay(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   + Ay(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
             &   + Ay(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))
             
       Jy(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) = gstep_invq(3)*(&
             &   - Ax(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   - Ax(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   - Ax(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   - Ax(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   + Ax(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   + Ax(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   + Ax(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
             &   + Ax(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  )) &
             - gstep_invq(1)*(&                     
             &   - Az(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   + Az(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   - Az(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   + Az(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   - Az(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   + Az(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   - Az(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
             &   + Az(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))  
             
        Jz(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) = gstep_invq(1)*(&
             &   - Ay(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   + Ay(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   - Ay(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   + Ay(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   - Ay(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   + Ay(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   - Ay(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
             &   + Ay(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))&
             - gstep_invq(2)*(&                       
             &   - Ax(1:ncm_tot(1)-1,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   - Ax(2:ncm_tot(1)  ,1:ncm_tot(2)-1,1:ncm_tot(3)-1) &
             &   + Ax(1:ncm_tot(1)-1,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   + Ax(2:ncm_tot(1)  ,2:ncm_tot(2)  ,1:ncm_tot(3)-1) &
             &   - Ax(1:ncm_tot(1)-1,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   - Ax(2:ncm_tot(1)  ,1:ncm_tot(2)-1,2:ncm_tot(3)  ) &
             &   + Ax(1:ncm_tot(1)-1,2:ncm_tot(2)  ,2:ncm_tot(3)  ) &
             &   + Ax(2:ncm_tot(1)  ,2:ncm_tot(2)  ,2:ncm_tot(3)  ))     
          
        Jtot = sqrt(Jx**2+Jy**2+Jz**2)
        
        varname = "Jcurre_"
        ival = int(int_centr(3))
        call extract_XY(Jtot,Jx,Jy,Jz,ival,centr,gstep,ncm_tot-1,varname,run_name,radius,prefix,planet,nrm_J,sgn,coord)
        ival = int(int_centr(2))
        call extract_XZ(Jtot,Jx,Jy,Jz,ival,centr,gstep,ncm_tot-1,varname,run_name,radius,prefix,planet,nrm_J,sgn,coord)   
        ival = int(int_centr(1)+2.5*radius/gstep(1))
        call extract_YZ(Jtot,Jx,Jy,Jz,ival,centr,gstep,ncm_tot-1,varname,run_name,radius,prefix,planet,nrm_J,sgn,coord)        
       deallocate(Jx,Jy,Jz,Jtot)
   call write_numdat_XML(fildat,planet,varname(1:3),.false.,.false.,.true.,.false.,.false.,'XY',0.,&
	ncm_tot,ncm_tot,ncm_tot,ncm_tot,ncm_tot,gstep,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),(/ 0.,0.,1. /))
   call write_granule_2D_XML(varname,'XY_'//run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
   call write_numdat_XML(fildat,planet,varname(1:3),.false.,.false.,.true.,.false.,.false.,'YZ',0.,&
	ncm_tot,ncm_tot,ncm_tot,ncm_tot,ncm_tot,gstep,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),(/ 1.,0.,0. /))
   call write_granule_2D_XML(varname,'YZ_'//run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
   call write_numdat_XML(fildat,planet,varname(1:3),.false.,.false.,.true.,.false.,.false.,'XZ',0.,&
	ncm_tot,ncm_tot,ncm_tot,ncm_tot,ncm_tot,gstep,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),(/ 0.,1.,0. /))
   call write_granule_2D_XML(varname,'XZ_'//run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
    endif
    deallocate(Ax,Ay,Az,Amod)
    
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
    subroutine read_species_cdf(diag_name,run_name,ion_label,varname1,varname2,varname3,varname4,varname5)
        
    character(len=*),intent(in) :: run_name,diag_name,ion_label
    character(len=*),intent(in) :: varname1,varname2,varname3,varname4,varname5
    
        
    !--Local variables
    real(dp),dimension(:,:,:),allocatable :: Ax,Ay,Az,A,B
    character(len=8) :: planet
    character(len=64) :: varname    
    integer  :: stId,ncid,var_id,ncid2
    integer :: ns,ncm_tot(3),int_centr(3)
    logical :: file_e
    character(len=50) :: write_name
    character(len=64) :: filename
    character(len=5) :: coord
    integer,dimension(nf90_max_var_dims) :: DimfieldId
    integer :: varfieldId,ival
    real(dp) :: gstep(3),dt,rtmp(3),radius,centr(3),box_size(6)
    real(dp) :: B0,V0,n0,x0,nrm_n,nrm_U,nrm_T,sgn,be
    real(dp),dimension(:),allocatable ::qms_d,rmds_d
   character(len=600)::numdatid
    
    !--Create file name for Field 
    !write(filename,'(a4,a1,2a)')trim(prefix),"_",trim(run_name),".nc"
    write(*,*) 'Open file'  ,diag_name  
  
    !--Inquire if the file exists
    inquire( file=trim(diag_name), exist=file_e )
    !--No file: return
    if(.not.(file_e)) then
     call wrtout(6,"File: "//trim(diag_name)//" does not exits","PERS")
     return
    endif
  
    !--Open NetCDF file
    stId = nf90_open(trim(diag_name), nf90_nowrite, ncid)
    call test_cdf(stId)

      
    !--Get number of point for fields in X,Y,Z
    StId = nf90_inq_varId(ncid,varname2,varfieldId)
    StId = nf90_Inquire_Variable(ncid,varfieldId,dimids = DimfieldId)
    StId = nf90_Inquire_Dimension(ncid,DimfieldId(1), len=ncm_tot(1))
    StId = nf90_Inquire_Dimension(ncid,DimfieldId(2), len=ncm_tot(2))
    StId = nf90_Inquire_Dimension(ncid,DimfieldId(3), len=ncm_tot(3))
    
    
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
            
   !--Get Number of Species
   call get_simple_variable_cdf(ncid,"ns",ns)
   allocate(qms_d(ns),rmds_d(ns))
   call get_simple_variable_cdf(ncid,"qms",qms_d(:))
   call get_simple_variable_cdf(ncid,"rmds",rmds_d(:))
   
   !--Allocate species and get informations
   ! call alloc_species(ns,Spe)
   ! call species_get_var_cdf(Spe,ncid)
   call get_simple_variable_cdf(ncid,"r_planet",radius)
   call get_simple_variable_cdf(ncid,"s_centr",centr)
   
   !--Get Planetname
   stId = nf90_inq_varid(ncid,"planetname",var_id)
   stId = nf90_get_var(ncid,var_id,planet)
   write(*,*) 'Planet :',planet
   !--Get normalisation values
   call get_simple_variable_cdf(ncid,"phys_mag",B0)
   call get_simple_variable_cdf(ncid,"phys_density",n0)
   call get_simple_variable_cdf(ncid,"phys_speed",V0)    
   call get_simple_variable_cdf(ncid,"phys_length",x0)
   call get_simple_variable_cdf(ncid,"betae",be)
   
   box_size=(/ -centr(1),ncm_tot(1)*gstep(1)-centr(1),-centr(2),ncm_tot(2)*gstep(2)-centr(2),-centr(3),ncm_tot(3)*gstep(3)-centr(3) /)
   box_size=box_size*x0
   !-- Allocation des tableaux de lecture
   allocate(A(ncm_tot(1),ncm_tot(2),ncm_tot(3)));	A = zero
   allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3)));	Ax = zero
   allocate(Ay(ncm_tot(1),ncm_tot(2),ncm_tot(3)));	Ay = zero
   allocate(Az(ncm_tot(1),ncm_tot(2),ncm_tot(3)));	Az = zero
   allocate(B(ncm_tot(1),ncm_tot(2),ncm_tot(3)));	B = zero
    if (varname1.ne."Density") write(*,*) "ici"
    write(*,*) varname1
         
    !if (varname1.ne."Density") call get_simple_variable_cdf(ncid,varname1 ,A(:,:,:) )
    call get_simple_variable_cdf(ncid,varname1 ,A(:,:,:) )
    call get_simple_variable_cdf(ncid,varname2 ,Ax(:,:,:) )
    call get_simple_variable_cdf(ncid,varname3 ,Ay(:,:,:) )
    call get_simple_variable_cdf(ncid,varname4 ,Az(:,:,:) )  
    !if (varname1.ne."Density") call get_simple_variable_cdf(ncid,varname5 ,B(:,:,:) )
    call get_simple_variable_cdf(ncid,varname5 ,B(:,:,:) )
    
    if (varname1.eq."Dn_tot") then 
    write(filename,'(a5,2a)')"Denw_",trim(run_name),".nc"
    stId = nf90_open(trim(filename), nf90_nowrite, ncid2)
    call test_cdf(stId)
    call get_simple_variable_cdf(ncid2,varname1 ,A(:,:,:) )
    stId = nf90_close(ncid2);  call test_cdf(stId)
    write(*,*) varname5
    write(filename,'(a5,2a)')"Prew_",trim(run_name),".nc"
    stId = nf90_open(trim(filename), nf90_nowrite, ncid2)
    call test_cdf(stId)
    call get_simple_variable_cdf(ncid2,varname5 ,B(:,:,:) )
    stId = nf90_close(ncid2);  call test_cdf(stId)
    endif
    !B=B/A
   int_centr = centr(:)/gstep(:)
     
     
    nrm_n = 1.!normalization done in the cube, in cm-3
    nrm_U = 1.     ! normalization done in the cube, in km/s
    nrm_T = 1. ! normalization done in the cube, in eV
    !if (varname1.eq."Dn_tot") nrm_T =2._dp*sum(qms_d*rmds_d)*b0**2/(2.*4.*3.141592*1e-7*1.6E-19)/n0

        ! coordinate transformation for each object
   ! coordinate system transformation is already made in the cube    
    select case (trim(planet))
     case ("mars","mercure")
        sgn = -1_dp
        coord='MSO'
     case ("ganymede")
        sgn=1_dp
        coord = 'GPhiO'
        !nrm_T = nrm_T*16._dp
     case("titan")
        sgn = 1_dp
        coord='TIIS'
        !nrm_T = nrm_T*16._dp
    case default
     planet="mars    "
     sgn = -1_dp
      coord='MSO'        
    end select
    varname = diag_name(1:4)
    if (varname.eq."Velw") varname="Thew"
     
    stId = nf90_close(ncid);  call test_cdf(stId)
    write(*,*) 'close file ',diag_name
    
    ival = int(int_centr(3))
    call extract_species_XY(A,Ax,Ay,Az,B,ival,centr,gstep,ncm_tot,varname,run_name,ion_label,radius,planet,nrm_n,nrm_U,nrm_T,sgn,coord)
    ival = int(int_centr(2))
    call extract_species_XZ(A,Ax,Ay,Az,B,ival,centr,gstep,ncm_tot,varname,run_name,ion_label,radius,planet,nrm_n,nrm_U,nrm_T,sgn,coord)   
    ival = int(int_centr(1)+2.5*radius/gstep(1))
    call extract_species_YZ(A,Ax,Ay,Az,B,ival,centr,gstep,ncm_tot,varname,run_name,ion_label,radius,planet,nrm_n,nrm_U,nrm_T,sgn,coord)     

   call write_numdat_XML(fildat,planet,varname(1:3),.false.,.false.,.true.,.false.,.false.,'XY',0.,&
	ncm_tot,ncm_tot,ncm_tot,ncm_tot,ncm_tot,gstep,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),(/ 0.,0.,1. /))
   call write_granule_2D_XML(varname,'XY_'//run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
   call write_numdat_XML(fildat,planet,varname(1:3),.false.,.false.,.true.,.false.,.false.,'YZ',0.,&
	ncm_tot,ncm_tot,ncm_tot,ncm_tot,ncm_tot,gstep,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),(/ 1.,0.,0. /))
   call write_granule_2D_XML(varname,'YZ_'//run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
   call write_numdat_XML(fildat,planet,varname(1:3),.false.,.false.,.true.,.false.,.false.,'XZ',0.,&
	ncm_tot,ncm_tot,ncm_tot,ncm_tot,ncm_tot,gstep,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),(/ 0.,1.,0. /))
   call write_granule_2D_XML(varname,'XZ_'//run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)

    
    deallocate(A,Ax,Ay,Az,B)
    deallocate(qms_d,rmds_d)
    
    end subroutine read_species_cdf     
    
  !!##############################################################################
  !! m_IMPEX_cut/extract_XY
  !! This routine extracts data in XY plane
  !! and creates a new ascii file 
  subroutine extract_XY(A,Ax,Ay,Az, cut_val,centr,gstep,ncm_tot,varname,run_name,radius,prefix,planetname,nrm,sgn,coord)
  character(len=*),intent(in) :: run_name,planetname
  character(len=*),intent(in) :: varname
  character(len=*),intent(in) :: prefix,coord  
  real(dp),dimension(:,:,:),intent(in) :: A,Ax,Ay,Az
  real(dp),intent(in) :: radius,centr(3),gstep(3),nrm,sgn
  integer,intent(in)  :: ncm_tot(3),cut_val
  ! Local variable
  real(dp),dimension(:,:),allocatable :: A_XY,Ax_XY,Ay_XY,Az_XY
  character(len=50) :: write_name
  character(len=2) :: plane
  integer :: iunit,i,j
  character(len=20) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
  
  allocate(A_XY(ncm_tot(1),ncm_tot(2)));	A_XY(:,:) = zero
  allocate(Ax_XY(ncm_tot(1),ncm_tot(2)));	Ax_XY(:,:) = zero
  allocate(Ay_XY(ncm_tot(1),ncm_tot(2)));	Ay_XY(:,:) = zero
  allocate(Az_XY(ncm_tot(1),ncm_tot(2)));	Az_XY(:,:) = zero
  ! extract information for Ax
  A_XY(:,:) =  A(:,:,cut_val)
  Ax_XY(:,:) =  Ax(:,:,cut_val)
  Ay_XY(:,:) =  Ay(:,:,cut_val)
  Az_XY(:,:) =  Az(:,:,cut_val)
  
  !-- filename of XY extracted dat
  write_name = trim(varname)//'_'//'XY'//'_'//trim(run_name(1:len_trim(run_name)-7))//".xml"
  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
  plane = 'XY'
  
  print *,'sgn,radius',sgn,radius
  
  iunit = 1
  open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
   	ACTION = 'write',STATUS = 'UNKNOWN')
   	
  if (trim(varname) /= "Jcurre") then 	
   	  call header_2Dcut_VOTABLE(iunit,prefix,plane,cut_val,planetname,coord) 	
  else
   	  call header_2Dcut_J_VOTABLE(iunit,prefix,plane,cut_val,planetname,coord) 	
  endif 	  
  do i=1,ncm_tot(1)
    do j = 1,ncm_tot(2)
         write(iunit,'(a)') '<TR>'
 call getstrf(sgn*(float(i)*gstep(1)-centr(1))/radius,tmp1,'(f9.3)')
 call getstrf(sgn*(float(j)*gstep(2)-centr(2))/radius,tmp2,'(f9.3)')
 call getstrf((float(cut_val)*gstep(3)-centr(3))/radius,tmp3,'(f9.3)')
 call  getstrf(A_XY(i,j),tmp4,'(e12.3)') ; call getstrf(Ax_XY(i,j),tmp5,'(e12.3)')
 call getstrf(Ay_XY(i,j),tmp6,'(e12.3)') ; call getstrf(Az_XY(i,j),tmp7,'(e12.3)')
         write(iunit,'(15a)') &
           '<TD>',trim(adjustl(tmp1)),'</TD><TD>',trim(adjustl(tmp2)),'</TD><TD>',trim(adjustl(tmp3)), &
	   '</TD><TD>',trim(adjustl(tmp4)),'</TD><TD>',trim(adjustl(tmp5)),'</TD><TD>',trim(adjustl(tmp6)), &
	   '</TD><TD>',trim(adjustl(tmp7)),'</TD>'
!         write(iunit,'(a,f9.3,a,f9.3,a,f9.3,a,4(f12.3,a))') &
!           '<TD>',sgn*(float(i)*gstep(1)-centr(1))/radius,'</TD><TD>',sgn*(float(j)*gstep(2)-centr(2))/radius, &
!           '</TD><TD>',(float(cut_val)*gstep(3)-centr(3))/radius,'</TD><TD>', A_XY(i,j),'</TD><TD>',&
!           Ax_XY(i,j),'</TD><TD>',Ay_XY(i,j),'</TD><TD>',Az_XY(i,j),'</TD>'
           write(iunit,'(a)') '</TR>'
    enddo
  enddo  
  ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESSOURCE>'
  write(iunit,'(a)') '</VOTABLE>'
  
  close(iunit)
  
  
!   write_name = trim(varname)//'_'//'XY'//'_'//trim(run_name)//".dat"
!    call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
!    plane = 'XY'
    
!    iunit = 1
!    open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
!     	ACTION = 'write',STATUS = 'UNKNOWN')
!    	
!    select case(trim(prefix))
!    case("Magw")
!      if (trim(varname)/= "Jcurre") then
!        write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Btot[nT]  |  Bx[nT]  |  By[nT]  |  Bz[nT]' 	
!      else
!        write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Jtot[mA/km2]  |  Jx[mA/km2]  |  Jy[mA/km2]  |  Jz[mA/km2]'
!      endif
!    case("Elew")
!      write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Etot[mV/km]  |  Ex[mV/km]  |  Ey[mV/km]  |  Ez[mV/km]'
!    case("Velw")
!      write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Utot[km/s]  |  Ux[km/s]  |  Uy[km/s]  |  Uz[km/s]'
!    end select
!   
!    do i=1,ncm_tot(1)
!      do j = 1,ncm_tot(2)
!           write(iunit,'(f9.3,f9.3,f9.3,4(f12.3))') &
!             sgn*(float(i)*gstep(1)-centr(1))/radius,sgn*(float(j)*gstep(2)-centr(2))/radius, &
!             (float(cut_val)*gstep(3)-centr(3))/radius, A_XY(i,j),&
!             Ax_XY(i,j),Ay_XY(i,j),Az_XY(i,j)
!      enddo
!    enddo  
   
!  close(iunit)
  
  deallocate(A_XY,Ax_XY,Ay_XY,Az_XY)
  end subroutine extract_XY
  
  
  !!##############################################################################
  !! m_IMPEX_cut/extract_species_XY
  !! This routine extracts data for a given ion speices in XY plane
  !! and creates a new ascii file 
  subroutine extract_species_XY(A,Ax,Ay,Az,B,cut_val,centr,gstep,ncm_tot,varname,run_name,ion_label,radius,planetname,nrm_n,nrm_U,nrm_T,sgn,coord)
  character(len=*),intent(in) :: run_name,planetname,ion_label
  character(len=*),intent(in) :: varname
  character(len=*),intent(in) :: coord  
  real(dp),dimension(:,:,:),intent(in) :: A,Ax,Ay,Az,B
  real(dp),intent(in) :: radius,centr(3),gstep(3),nrm_n,nrm_U,nrm_T,sgn
  integer,intent(in)  :: ncm_tot(3),cut_val
  ! Local variable
  real(dp),dimension(:,:),allocatable :: A_XY,Ax_XY,Ay_XY,Az_XY,Amod_XY,B_XY
  character(len=50) :: write_name
  character(len=20) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
  character(len=2) :: plane
  integer :: iunit,i,j
  
  
  allocate(A_XY(ncm_tot(1),ncm_tot(2)));	A_XY(:,:) = zero
  allocate(Ax_XY(ncm_tot(1),ncm_tot(2)));	Ax_XY(:,:) = zero
  allocate(Ay_XY(ncm_tot(1),ncm_tot(2)));	Ay_XY(:,:) = zero
  allocate(Az_XY(ncm_tot(1),ncm_tot(2)));	Az_XY(:,:) = zero
  allocate(Amod_XY(ncm_tot(1),ncm_tot(2)));	Amod_XY(:,:) = zero
  allocate(B_XY(ncm_tot(1),ncm_tot(2)));	B_XY(:,:) = zero
  ! extract information 
  A_XY(:,:) =  A(:,:,cut_val)
  Ax_XY(:,:) =  Ax(:,:,cut_val)
  Ay_XY(:,:) =  Ay(:,:,cut_val)
  Az_XY(:,:) =  Az(:,:,cut_val)
  Amod_XY = sqrt(Ax_XY**2+Ay_XY**2+Az_XY**2)
  B_XY(:,:) =  B(:,:,cut_val)
  
  !-- filename of XY extracted dat
  write_name = trim(varname)//'_'//'XY'//'_'//trim(run_name(1:len_trim(run_name)-7))//".xml"
  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
  plane = 'XY'
  
  iunit = 1
  open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
   	ACTION = 'write',STATUS = 'UNKNOWN')
   	
 call header_2Dcut_species_VOTABLE(iunit,ion_label,plane,cut_val,planetname,coord) 
  do i=1,ncm_tot(1)
    do j = 1,ncm_tot(2)
         write(iunit,'(a)') '<TR>'
!         write(iunit,'(a,f9.3,a,f9.3,a,f9.3,a,6(f12.3,a))') &
!           '<TD>',sgn*(float(i)*gstep(1)-centr(1))/radius,'</TD><TD>',sgn*(float(j)*gstep(2)-centr(2))/radius, &
!           '</TD><TD>',(float(cut_val)*gstep(3)-centr(3))/radius,'</TD><TD>', A_XY(i,j),'</TD><TD>',&
!           Ax_XY(i,j),'</TD><TD>',Ay_XY(i,j),'</TD><TD>',Az_XY(i,j),'</TD><TD>',Amod_XY(i,j), &
!           '</TD><TD>',B_XY(i,j),'</TD>'
 call getstrf(sgn*(float(i)*gstep(1)-centr(1))/radius,tmp1,'(f9.3)')
 call getstrf(sgn*(float(j)*gstep(2)-centr(2))/radius,tmp2,'(f9.3)')
 call getstrf((float(cut_val)*gstep(3)-centr(3))/radius,tmp3,'(f9.3)')
 call  getstrf(A_XY(i,j),tmp4,'(e12.3)') ; call getstrf(Ax_XY(i,j),tmp5,'(e12.3)')
 call getstrf(Ay_XY(i,j),tmp6,'(e12.3)') ; call getstrf(Az_XY(i,j),tmp7,'(e12.3)')
 call getstrf(Amod_XY(i,j),tmp8,'(e12.3)') ; call getstrf(B_XY(i,j),tmp9,'(e12.3)')
         write(iunit,'(19a)') &
           '<TD>',trim(adjustl(tmp1)),'</TD><TD>',trim(adjustl(tmp2)),'</TD><TD>',trim(adjustl(tmp3)), &
	   '</TD><TD>',trim(adjustl(tmp4)),'</TD><TD>',trim(adjustl(tmp5)),'</TD><TD>',trim(adjustl(tmp6)), &
	   '</TD><TD>',trim(adjustl(tmp7)),'</TD><TD>',trim(adjustl(tmp8)),'</TD><TD>',trim(adjustl(tmp9)),'</TD>'
           write(iunit,'(a)') '</TR>'
    enddo
  enddo  
  ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESSOURCE>'
  write(iunit,'(a)') '</VOTABLE>'
  
  close(iunit)
  
!  write_name = trim(varname)//'_'//'XY'//'_'//trim(run_name)//".dat"
!  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
!  plane = 'XY'
!  
!  iunit = 1
!  open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
!   	ACTION = 'write',STATUS = 'UNKNOWN')
!   	
!  write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] | n[cm-3]  |  Ux[km/s]  |  Uy[km/s]  |  Uz[km/s] | Utot[km/s] | T[eV]'
!  do i=1,ncm_tot(1)
!    do j = 1,ncm_tot(2)
!         write(iunit,'(f9.3,f9.3,f9.3,6(f12.3))') &
!           sgn*(float(i)*gstep(1)-centr(1))/radius,sgn*(float(j)*gstep(2)-centr(2))/radius, &
!           (float(cut_val)*gstep(3)-centr(3))/radius, A_XY(i,j),&
!           Ax_XY(i,j),Ay_XY(i,j),Az_XY(i,j),Amod_XY(i,j), &
!           B_XY(i,j)
!    enddo
!  enddo   
!  close(iunit)  
  
  
  deallocate(A_XY,Ax_XY,Ay_XY,Az_XY,B_XY,Amod_XY)
  end subroutine extract_species_XY  
  
  
  !!##############################################################################
  !! m_IMPEX_cut/extract_XZ
  !! This routine extracts data in XZ plane
  !! and creates a new ascii file 
  subroutine extract_XZ(A,Ax,Ay,Az, cut_val,centr,gstep,ncm_tot,varname,run_name,radius,prefix,planetname,nrm,sgn,coord)
  character(len=*),intent(in) :: run_name,planetname
  character(len=*),intent(in) :: varname
  character(len=*),intent(in) :: prefix,coord  
  real(dp),dimension(:,:,:),intent(in) :: A,Ax,Ay,Az
  real(dp),intent(in) :: radius,centr(3),gstep(3),nrm,sgn
  integer,intent(in)  :: ncm_tot(3),cut_val
  ! Local variable
  real(dp),dimension(:,:),allocatable :: A_XZ,Ax_XZ,Ay_XZ,Az_XZ
  character(len=50) :: write_name
  character(len=2) :: plane
  integer :: iunit,i,k
  character(len=20) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
  
  allocate(A_XZ(ncm_tot(1),ncm_tot(3)));	A_XZ(:,:) = zero
  allocate(Ax_XZ(ncm_tot(1),ncm_tot(3)));	Ax_XZ(:,:) = zero
  allocate(Ay_XZ(ncm_tot(1),ncm_tot(3)));	Ay_XZ(:,:) = zero
  allocate(Az_XZ(ncm_tot(1),ncm_tot(3)));	Az_XZ(:,:) = zero  
  ! extract information for Ax
  A_XZ(:,:) =  A(:,cut_val,:)
  Ax_XZ(:,:) =  Ax(:,cut_val,:)
  Ay_XZ(:,:) =  Ay(:,cut_val,:)
  Az_XZ(:,:) =  Az(:,cut_val,:)
  
  !-- filename of XY extracted dat
  write_name = trim(varname)//'_'//'XZ'//'_'//trim(run_name(1:len_trim(run_name)-7))//".xml"
  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
  plane = 'XZ'
  

  iunit = 1
  open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
   	ACTION = 'write',STATUS = 'UNKNOWN')
  if (trim(varname) /= "Jcurre") then 	
   	  call header_2Dcut_VOTABLE(iunit,prefix,plane,cut_val,planetname,coord) 	
  else
   	  call header_2Dcut_J_VOTABLE(iunit,prefix,plane,cut_val,planetname,coord) 	
  endif 	 	
  do i=1,ncm_tot(1)
    do k = 1,ncm_tot(3)
         write(iunit,'(a)') '<TR>'
 call getstrf(sgn*(float(i)*gstep(1)-centr(1))/radius,tmp1,'(f9.3)')
 call getstrf(sgn*(float(cut_val)*gstep(2)-centr(2))/radius,tmp2,'(f9.3)')
 call getstrf((float(k)*gstep(3)-centr(3))/radius,tmp3,'(f9.3)')
 call  getstrf(A_XZ(i,k),tmp4,'(e12.3)') ; call getstrf(Ax_XZ(i,k),tmp5,'(e12.3)')
 call getstrf(Ay_XZ(i,k),tmp6,'(e12.3)') ; call getstrf(Az_XZ(i,k),tmp7,'(e12.3)')
         write(iunit,'(15a)') &
           '<TD>',trim(adjustl(tmp1)),'</TD><TD>',trim(adjustl(tmp2)),'</TD><TD>',trim(adjustl(tmp3)), &
	   '</TD><TD>',trim(adjustl(tmp4)),'</TD><TD>',trim(adjustl(tmp5)),'</TD><TD>',trim(adjustl(tmp6)), &
	   '</TD><TD>',trim(adjustl(tmp7)),'</TD>'
!         write(iunit,'(a,f9.3,a,f9.3,a,f9.3,a,4(f12.3,a))') &
!           '<TD>',sgn*(float(i)*gstep(1)-centr(1))/radius,'</TD> <TD>', sgn*(float(cut_val)*gstep(2)-centr(2))/radius, &           
!           '</TD> <TD>',(float(k)*gstep(3)-centr(3))/radius, &
!           '</TD> <TD>',A_XZ(i,k),'</TD> <TD>', &
!           Ax_XZ(i,k),'</TD> <TD>',Ay_XZ(i,k),'</TD> <TD>',Az_XZ(i,k),'</TD>'
           write(iunit,'(a)') '</TR>'
    enddo
  enddo  
  ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESSOURCE>'
  write(iunit,'(a)') '</VOTABLE>'  
  close(iunit)
  
!-- filename of XY extracted dat
!  write_name = trim(varname)//'_'//'XZ'//'_'//trim(run_name)//".dat"
!  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
!  plane = 'XZ'
  

!  iunit = 1
!  open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
!   	ACTION = 'write',STATUS = 'UNKNOWN')
!    select case(trim(prefix))
!    case("Magw")
!      if (trim(varname)/= "Jcurre") then
!        write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Btot[nT]  |  Bx[nT]  |  By[nT]  |  Bz[nT]' 	
!      else
!        write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Jtot[mA/km2]  |  Jx[mA/km2]  |  Jy[mA/km2]  |  Jz[mA/km2]'
!      endif
!    case("Elew")
!      write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Etot[mV/km]  |  Ex[mV/km]  |  Ey[mV/km]  |  Ez[mV/km]'
!    case("Velw")
!      write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Utot[km/s]  |  Ux[km/s]  |  Uy[km/s]  |  Uz[km/s]'
!    end select  
	 	
!  do i=1,ncm_tot(1)
!    do k = 1,ncm_tot(3)
!         write(iunit,'(f9.3,f9.3,f9.3,4(f12.3))') &
!           sgn*(float(i)*gstep(1)-centr(1))/radius, sgn*(float(cut_val)*gstep(2)-centr(2))/radius, &           
!           (float(k)*gstep(3)-centr(3))/radius, &
!           A_XZ(i,k), &
!           Ax_XZ(i,k),Ay_XZ(i,k),Az_XZ(i,k)
!
!    enddo
!  enddo  
 
!  close(iunit)
    
  
  
  deallocate(A_XZ,Ax_XZ,Ay_XZ,Az_XZ)
  end subroutine extract_XZ  
  
!!##############################################################################
  !! m_IMPEX_cut/extract_species_XZ
  !! This routine extracts data in XZ plane
  !! and creates a new ascii file 
  subroutine extract_species_XZ(A,Ax,Ay,Az,B,cut_val,centr,gstep,ncm_tot,varname,run_name,ion_label,radius,planetname,nrm_n,nrm_U,nrm_T,sgn,coord)
  character(len=*),intent(in) :: run_name,planetname,ion_label
  character(len=*),intent(in) :: varname
  character(len=*),intent(in) :: coord  
  real(dp),dimension(:,:,:),intent(in) :: A,Ax,Ay,Az,B
  real(dp),intent(in) :: radius,centr(3),gstep(3),nrm_n,nrm_T,nrm_U,sgn
  integer,intent(in)  :: ncm_tot(3),cut_val
  ! Local variable
  real(dp),dimension(:,:),allocatable :: A_XZ,Ax_XZ,Ay_XZ,Az_XZ,Amod_XZ,B_XZ
  character(len=50) :: write_name
  character(len=2) :: plane
  character(len=20) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
  integer :: iunit,i,k
  
  allocate(A_XZ(ncm_tot(1),ncm_tot(3)));	A_XZ(:,:) = zero
  allocate(Ax_XZ(ncm_tot(1),ncm_tot(3)));	Ax_XZ(:,:) = zero
  allocate(Ay_XZ(ncm_tot(1),ncm_tot(3)));	Ay_XZ(:,:) = zero
  allocate(Az_XZ(ncm_tot(1),ncm_tot(3)));	Az_XZ(:,:) = zero
  allocate(Amod_XZ(ncm_tot(1),ncm_tot(3)));	Amod_XZ(:,:) = zero
  allocate(B_XZ(ncm_tot(1),ncm_tot(3)));	B_XZ(:,:) = zero
  ! extract information for Ax
  A_XZ(:,:) =  A(:,cut_val,:)
  Ax_XZ(:,:) =  Ax(:,cut_val,:)
  Ay_XZ(:,:) =  Ay(:,cut_val,:)
  Az_XZ(:,:) =  Az(:,cut_val,:)
  Amod_XZ = sqrt(Ax_XZ**2+Ay_XZ**2+Az_XZ**2)
  B_XZ(:,:) =  B(:,cut_val,:)
  
  !-- filename of XZ extracted dat
  write_name = trim(varname)//'_'//'XZ'//'_'//trim(run_name(1:len_trim(run_name)-7))//".xml"
  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
  plane = 'XZ'
  

  iunit = 1
  open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
   	ACTION = 'write',STATUS = 'UNKNOWN')
call header_2Dcut_species_VOTABLE(iunit,ion_label,plane,cut_val,planetname,coord) 

  do i=1,ncm_tot(1)
    do k = 1,ncm_tot(3)
         write(iunit,'(a)') '<TR>'
 call getstrf(sgn*(float(i)*gstep(1)-centr(1))/radius,tmp1,'(f9.3)')
 call getstrf(sgn*(float(cut_val)*gstep(2)-centr(2))/radius,tmp2,'(f9.3)')
 call getstrf((float(k)*gstep(3)-centr(3))/radius,tmp3,'(f9.3)')
 call  getstrf(A_XZ(i,k),tmp4,'(e12.3)') ; call getstrf(Ax_XZ(i,k),tmp5,'(e12.3)')
 call getstrf(Ay_XZ(i,k),tmp6,'(e12.3)') ; call getstrf(Az_XZ(i,k),tmp7,'(e12.3)')
 call getstrf(Amod_XZ(i,k),tmp8,'(e12.3)') ; call getstrf(B_XZ(i,k),tmp9,'(e12.3)')
         write(iunit,'(19a)') &
           '<TD>',trim(adjustl(tmp1)),'</TD><TD>',trim(adjustl(tmp2)),'</TD><TD>',trim(adjustl(tmp3)), &
	   '</TD><TD>',trim(adjustl(tmp4)),'</TD><TD>',trim(adjustl(tmp5)),'</TD><TD>',trim(adjustl(tmp6)), &
	   '</TD><TD>',trim(adjustl(tmp7)),'</TD><TD>',trim(adjustl(tmp8)),'</TD><TD>',trim(adjustl(tmp9)),'</TD>'
!         write(iunit,'(a,f9.3,a,f9.3,a,f9.3,a,6(f12.3,a))') &
!           '<TD>',sgn*(float(i)*gstep(1)-centr(1))/radius,'</TD> <TD>', sgn*(float(cut_val)*gstep(2)-centr(2))/radius, &           
!           '</TD> <TD>',(float(k)*gstep(3)-centr(3))/radius, &
!           '</TD> <TD>',A_XZ(i,k),'</TD> <TD>', &
!           Ax_XZ(i,k),'</TD> <TD>',Ay_XZ(i,k),'</TD> <TD>',Az_XZ(i,k),'</TD> <TD>', &
!           Amod_XZ(i,k),'</TD> <TD>',B_XZ(i,k),'</TD>'
           write(iunit,'(a)') '</TR>'
    enddo
  enddo  
  ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESSOURCE>'
  write(iunit,'(a)') '</VOTABLE>'  
  close(iunit)
  
  
  !-- filename of XZ extracted dat
!  write_name = trim(varname)//'_'//'XZ'//'_'//trim(run_name)//".dat"
!  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
!  plane = 'XZ'
!  iunit = 1
!  open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
!   	ACTION = 'write',STATUS = 'UNKNOWN')

!write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] | n[cm-3]  |  Ux[km/s]  |  Uy[km/s]  |  Uz[km/s] | Utot[km/s] | T[eV]'
!  do i=1,ncm_tot(1)
!    do k = 1,ncm_tot(3)
!         write(iunit,'(f9.3,f9.3,f9.3,6(f12.3))') &
!           sgn*(float(i)*gstep(1)-centr(1))/radius, sgn*(float(cut_val)*gstep(2)-centr(2))/radius, &           
!           (float(k)*gstep(3)-centr(3))/radius, &
!           A_XZ(i,k), &
!           Ax_XZ(i,k),Ay_XZ(i,k),Az_XZ(i,k), &
!           Amod_XZ(i,k),B_XZ(i,k)
!    enddo
!  enddo  
!  close(iunit)  
  
  deallocate(A_XZ,Ax_XZ,Ay_XZ,Az_XZ,Amod_XZ,B_XZ)
  end subroutine extract_species_XZ    

 !!##############################################################################
  !! m_IMPEX_cut/extract_YZ
  !! This routine extracts data in YZ plane
  !! and creates a new ascii file 
  subroutine extract_YZ(A,Ax,Ay,Az, cut_val,centr,gstep,ncm_tot,varname,run_name,radius,prefix,planetname,nrm,sgn,coord)
  character(len=*),intent(in) :: run_name,planetname
  character(len=*),intent(in) :: varname
  character(len=*),intent(in) :: prefix,coord  
  real(dp),dimension(:,:,:),intent(in) :: A,Ax,Ay,Az
  real(dp),intent(in) :: radius,centr(3),gstep(3),nrm,sgn
  integer,intent(in)  :: ncm_tot(3),cut_val
  ! Local variable
  real(dp),dimension(:,:),allocatable :: A_YZ,Ax_YZ,Ay_YZ,Az_YZ
  character(len=50) :: write_name
  character(len=2) :: plane
  integer :: iunit,j,k
  character(len=20) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
  
  allocate(A_YZ(ncm_tot(2),ncm_tot(3)));	A_YZ(:,:) = zero
  allocate(Ax_YZ(ncm_tot(2),ncm_tot(3)));	Ax_YZ(:,:) = zero
  allocate(Ay_YZ(ncm_tot(2),ncm_tot(3)));	Ay_YZ(:,:) = zero
  allocate(Az_YZ(ncm_tot(2),ncm_tot(3)));	Az_YZ(:,:) = zero  
  ! extract information for Ax
  A_YZ(:,:) =  A(cut_val,:,:)
  Ax_YZ(:,:) =  Ax(cut_val,:,:)
  Ay_YZ(:,:) =  Ay(cut_val,:,:)
  Az_YZ(:,:) =  Az(cut_val,:,:)
  
  !-- filename of XY extracted dat
  write_name = trim(varname)//'_'//'YZ'//'_'//trim(run_name(1:len_trim(run_name)-7))//".xml"
  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
  plane = 'YZ'
  
  
  iunit = 1
  open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
   	ACTION = 'write',STATUS = 'UNKNOWN')
    if (trim(varname) /= "Jcurre") then 	
     	  call header_2Dcut_VOTABLE(iunit,prefix,plane,cut_val,planetname,coord) 	
    else
     	  call header_2Dcut_J_VOTABLE(iunit,prefix,plane,cut_val,planetname,coord) 	
  endif 
  	
  do j=1,ncm_tot(2)
    do k = 1,ncm_tot(3)
         write(iunit,'(a)') '<TR>'    
 call getstrf(sgn*(float(cut_val)*gstep(1)-centr(1))/radius,tmp1,'(f9.3)')
 call getstrf(sgn*(float(j)*gstep(2)-centr(2))/radius,tmp2,'(f9.3)')
 call getstrf((float(k)*gstep(3)-centr(3))/radius,tmp3,'(f9.3)')
 call getstrf(A_YZ(j,k),tmp4,'(e12.3)') ; call getstrf(Ax_YZ(j,k),tmp5,'(e12.3)')
 call getstrf(Ay_YZ(j,k),tmp6,'(e12.3)') ; call getstrf(Az_YZ(j,k),tmp7,'(e12.3)')
         write(iunit,'(15a)') &
           '<TD>',trim(adjustl(tmp1)),'</TD><TD>',trim(adjustl(tmp2)),'</TD><TD>',trim(adjustl(tmp3)), &
	   '</TD><TD>',trim(adjustl(tmp4)),'</TD><TD>',trim(adjustl(tmp5)),'</TD><TD>',trim(adjustl(tmp6)), &
	   '</TD><TD>',trim(adjustl(tmp7)),'</TD>'
!         write(iunit,'(a,f9.3,a,f9.3,a,f9.3,a,4(f12.3,a))') &
!           '<TD>',sgn*(float(cut_val)*gstep(1)-centr(1))/radius,'</TD> <TD>',sgn*(float(j)*gstep(2)-centr(2))/radius, &
!           '</TD> <TD>',(float(k)*gstep(3)-centr(3))/radius, &
!           '</TD> <TD>',A_YZ(j,k),'</TD> <TD>', &
!           Ax_YZ(j,k),'</TD> <TD>',Ay_YZ(j,k),'</TD> <TD>',Az_YZ(j,k),'</TD>'
         write(iunit,'(a)') '</TR>'           
    enddo
  enddo  
  ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESSOURCE>'
  write(iunit,'(a)') '</VOTABLE>'  
  close(iunit)
  
   
  !-- filename of XY extracted dat
!  write_name = trim(varname)//'_'//'YZ'//'_'//trim(run_name)//".dat"
!  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
!  plane = 'YZ' 
!  iunit = 1
!  open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
!   	ACTION = 'write',STATUS = 'UNKNOWN')
!    select case(trim(prefix))
!    case("Magw")
!      if (trim(varname)/= "Jcurre") then
!        write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Btot[nT]  |  Bx[nT]  |  By[nT]  |  Bz[nT]' 	
!      else
!        write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Jtot[mA/km2]  |  Jx[mA/km2]  |  Jy[mA/km2]  |  Jz[mA/km2]'
!      endif
!    case("Elew")
!      write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Etot[mV/km]  |  Ex[mV/km]  |  Ey[mV/km]  |  Ez[mV/km]'
!    case("Velw")
!      write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] |  Utot[km/s]  |  Ux[km/s]  |  Uy[km/s]  |  Uz[km/s]'
!    end select
   	
!  do j=1,ncm_tot(2)
!    do k = 1,ncm_tot(3)
!         write(iunit,'(f9.3,f9.3,f9.3,4(f12.3))') &
!           sgn*(float(cut_val)*gstep(1)-centr(1))/radius,sgn*(float(j)*gstep(2)-centr(2))/radius, &
!           (float(k)*gstep(3)-centr(3))/radius, &
!           A_YZ(j,k), &
!           Ax_YZ(j,k),Ay_YZ(j,k),Az_YZ(j,k)   
!    enddo
!  enddo  
!  close(iunit) 
  
  deallocate(A_YZ,Ax_YZ,Ay_YZ,Az_YZ)
  end subroutine extract_YZ  

 !!##############################################################################
  !! m_IMPEX_cut/extract_species_YZ
  !! This routine extracts data in YZ plane
  !! and creates a new ascii file 
  subroutine extract_species_YZ(A,Ax,Ay,Az,B,cut_val,centr,gstep,ncm_tot,varname,run_name,ion_label,radius,planetname,nrm_n,nrm_U,nrm_T,sgn,coord)
  character(len=*),intent(in) :: run_name,planetname,ion_label
  character(len=*),intent(in) :: varname
  character(len=*),intent(in) :: coord  
  real(dp),dimension(:,:,:),intent(in) :: A,Ax,Ay,Az,B
  real(dp),intent(in) :: radius,centr(3),gstep(3),nrm_U,nrm_n,nrm_T,sgn
  integer,intent(in)  :: ncm_tot(3),cut_val
  ! Local variable
  real(dp),dimension(:,:),allocatable :: A_YZ,Ax_YZ,Ay_YZ,Az_YZ,Amod_YZ,B_YZ
  character(len=50) :: write_name
  character(len=2) :: plane
  integer :: iunit,j,k
  character(len=20) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
  
  allocate(A_YZ(ncm_tot(2),ncm_tot(3)));	A_YZ(:,:) = zero
  allocate(Ax_YZ(ncm_tot(2),ncm_tot(3)));	Ax_YZ(:,:) = zero
  allocate(Ay_YZ(ncm_tot(2),ncm_tot(3)));	Ay_YZ(:,:) = zero
  allocate(Az_YZ(ncm_tot(2),ncm_tot(3)));	Az_YZ(:,:) = zero  
  allocate(Amod_YZ(ncm_tot(2),ncm_tot(3)));	Amod_YZ(:,:) = zero
  allocate(B_YZ(ncm_tot(2),ncm_tot(3)));	B_YZ(:,:) = zero
  ! extract information for Ax
  A_YZ(:,:) =  A(cut_val,:,:)
  Ax_YZ(:,:) =  Ax(cut_val,:,:)
  Ay_YZ(:,:) =  Ay(cut_val,:,:)
  Az_YZ(:,:) =  Az(cut_val,:,:)
  Amod_YZ = sqrt(Ax_YZ**2+Ay_YZ**2+Az_YZ**2)
  B_YZ(:,:) = B(cut_val,:,:)
  
  !-- filename of XY extracted dat
  write_name = trim(varname)//'_'//'YZ'//'_'//trim(run_name(1:len_trim(run_name)-7))//".xml"
  call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
  plane = 'YZ'
  
  
  iunit = 1
  open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
   	ACTION = 'write',STATUS = 'UNKNOWN')
  call header_2Dcut_species_VOTABLE(iunit,ion_label,plane,cut_val,planetname,coord) 
  	
  do j=1,ncm_tot(2)
    do k = 1,ncm_tot(3)
         write(iunit,'(a)') '<TR>'    
 call getstrf(sgn*(float(cut_val)*gstep(1)-centr(1))/radius,tmp1,'(f9.3)')
 call getstrf(sgn*(float(j)*gstep(2)-centr(2))/radius,tmp2,'(f9.3)')
 call getstrf((float(k)*gstep(3)-centr(3))/radius,tmp3,'(f9.3)')
 call getstrf(A_YZ(j,k),tmp4,'(e12.3)') ; call getstrf(Ax_YZ(j,k),tmp5,'(e12.3)')
 call getstrf(Ay_YZ(j,k),tmp6,'(e12.3)') ; call getstrf(Az_YZ(j,k),tmp7,'(e12.3)')
 call getstrf(Amod_YZ(j,k),tmp8,'(e12.3)') ; call getstrf(B_YZ(j,k),tmp9,'(e12.3)')
         write(iunit,'(19a)') &
           '<TD>',trim(adjustl(tmp1)),'</TD><TD>',trim(adjustl(tmp2)),'</TD><TD>',trim(adjustl(tmp3)), &
	   '</TD><TD>',trim(adjustl(tmp4)),'</TD><TD>',trim(adjustl(tmp5)),'</TD><TD>',trim(adjustl(tmp6)), &
	   '</TD><TD>',trim(adjustl(tmp7)),'</TD><TD>',trim(adjustl(tmp8)),'</TD><TD>',trim(adjustl(tmp9)),'</TD>'
!         write(iunit,'(a,f9.3,a,f9.3,a,f9.3,a,6(f12.3,a))') &
!           '<TD>',sgn*(float(cut_val)*gstep(1)-centr(1))/radius,'</TD> <TD>',sgn*(float(j)*gstep(2)-centr(2))/radius, &
!           '</TD> <TD>',(float(k)*gstep(3)-centr(3))/radius, &
!           '</TD> <TD>',A_YZ(j,k),'</TD> <TD>', &
!           Ax_YZ(j,k),'</TD> <TD>',Ay_YZ(j,k),'</TD> <TD>',Az_YZ(j,k),'</TD> <TD>', &
!           Amod_YZ(j,k),'</TD> <TD>',B_YZ(j,k),'</TD>'
         write(iunit,'(a)') '</TR>'           
    enddo
  enddo  
  ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESSOURCE>'
  write(iunit,'(a)') '</VOTABLE>'  
  close(iunit)
    !-- filename of XY extracted dat
!    write_name = trim(varname)//'_'//'YZ'//'_'//trim(run_name)//".dat"
!    call wrtout(6," ======= Creation of file : "//trim(write_name),"PERS")
!    plane = 'YZ'
!   
!    iunit = 1
!    open(UNIT = iunit, FILE = write_name, FORM = 'FORMATTED',  &
!     	ACTION = 'write',STATUS = 'UNKNOWN')
! write(iunit,'(a)') 'X [RM]  |  Y [RM]  |  Z [RM] | n[cm-3]  |  Ux[km/s]  |  Uy[km/s]  |  Uz[km/s] | Utot[km/s] | T[eV]'  
!    	
!    do j=1,ncm_tot(2)
!      do k = 1,ncm_tot(3)  
!           write(iunit,'(f9.3,f9.3,f9.3,6(f12.3))') &
!             sgn*(float(cut_val)*gstep(1)-centr(1))/radius,sgn*(float(j)*gstep(2)-centr(2))/radius, &
!             (float(k)*gstep(3)-centr(3))/radius, &
!             A_YZ(j,k), &
!             Ax_YZ(j,k),Ay_YZ(j,k),Az_YZ(j,k), &
!             Amod_YZ(j,k),B_YZ(j,k)           
!      enddo
!    enddo  
!  close(iunit)
  
  deallocate(A_YZ,Ax_YZ,Ay_YZ,Az_YZ,Amod_YZ,B_YZ)
  end subroutine extract_species_YZ  
  
!!##########################################################
!! m_IMPEX_cut/header_2Dcut_IMPEX
!! This routine writes the header of the xml
!! file for the 2D cut IMPEX diag
subroutine header_2Dcut_VOTABLE(iunit,prefix,plane,cutval,planetname,coord)
integer, intent(in) :: iunit,cutval
character(len=*),intent(in) :: prefix,plane,planetname,coord
 
  write(iunit,'(a)') '<?xml version="1.0"?>'
  write(iunit,'(a)') '<VOTABLE version="1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
  write(iunit,'(a)')  'xmlns="http://www.ivoa.net/xml/VOTable/v1.2" xmlns:stc="http://www.ivoa.net/xml/STC/v1.30" >'  
  write(iunit,'(a)')  '<RESSOURCE name="IMPEx 2D cut">'
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
    write(iunit,'(3a)') '<FIELD name="X" ID="col1" ucd="pos.cartesian.x" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C1" datatype="float" width="9" />' 
    write(iunit,'(3a)') '<FIELD name="Y" ID="col2" ucd="pos.cartesian.y" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C2" datatype="float" width="9" />' 
    write(iunit,'(3a)') '<FIELD name="Z" ID="col3" ucd="pos.cartesian.z" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C3" datatype="float" width="9" />'
    write(iunit,'(3a)') '<FIELD name="Btot" ID="col4" ucd="phys.magField" ',&
        '  utype="" datatype="float" width="12" unit="nT" />'    	
    write(iunit,'(3a)') '<FIELD name="Bx" ID="col5" ucd="phys.magField" ',&
    	'  utype="" datatype="float" width="12" unit="nT" />'
    write(iunit,'(3a)') '<FIELD name="By" ID="col6" ucd="phys.magField" ',&
        '  utype="" datatype="float" width="12" unit="nT" />'
    write(iunit,'(3a)') '<FIELD name="Bz" ID="col7" ucd="phys.magField" ',&
        '  utype="" datatype="float" width="12" unit="nT" />'
        
    
  case("Elew")
    write(iunit,'(a)') '<DESCRIPTION> Electric field components </DESCRIPTION>'    
    write(iunit,'(3a)') '<FIELD name="X" ID="col1" ucd="pos.cartesian.x" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C1" datatype="float" width="9" />' 
    write(iunit,'(3a)') '<FIELD name="Y" ID="col2" ucd="pos.cartesian.y" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C2" datatype="float" width="9" />' 
    write(iunit,'(3a)') '<FIELD name="Z" ID="col3" ucd="pos.cartesian.z" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C3" datatype="float" width="9" />'
    write(iunit,'(3a)') '<FIELD name="Etot" ID="col4" ucd="phys.electField" ',&
        '  utype="" datatype="float" width="12" unit="mV.m-1" />'    	
    write(iunit,'(3a)') '<FIELD name="Ex" ID="col5" ucd="phys.electField" ',&
        '  utype="" datatype="float" width="12" unit="mV.m-1" />'
    write(iunit,'(3a)') '<FIELD name="Ey" ID="col6" ucd="phys.electField" ',&
        '  utype="" datatype="float" width="12" unit="mV.m-1" />'
    write(iunit,'(3a)') '<FIELD name="Ez" ID="col7" ucd="phys.electField" ',&
        '  utype="" datatype="float" width="12" unit="mV.m-1" />'
  case("Velw")
    write(iunit,'(a)') '<DESCRIPTION> Velocity components </DESCRIPTION>'  
    write(iunit,'(3a)') '<FIELD name="X" ID="col1" ucd="pos.cartesian.x" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C1" datatype="float" width="9" />' 
    write(iunit,'(3a)') '<FIELD name="Y" ID="col2" ucd="pos.cartesian.y" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C2" datatype="float" width="9" />' 
    write(iunit,'(3a)') '<FIELD name="Z" ID="col3" ucd="pos.cartesian.z" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C3" datatype="float" width="9" />'
    write(iunit,'(3a)') '<FIELD name="Utot" ID="col4" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1" />'
    write(iunit,'(3a)') '<FIELD name="Ux" ID="col5" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1" />'
    write(iunit,'(3a)') '<FIELD name="Uy" ID="col6" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1" />'
    write(iunit,'(3a)') '<FIELD name="Uz" ID="col7" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1" />'    
  end select
  
  
  write(iunit,'(a)') '<DATA>'
  write(iunit,'(a)') '<TABLEDATA>'  
  

end subroutine header_2Dcut_VOTABLE 

  
!!##########################################################
!! m_IMPEX_cut/header_2Dcut_J_IMPEX
!! This routine writes the header of the xml
!! file for the 2D cut IMPEX diag
subroutine header_2Dcut_J_VOTABLE(iunit,prefix,plane,cutval,planetname,coord)
integer, intent(in) :: iunit,cutval
character(len=*),intent(in) :: prefix,plane,planetname,coord
 
  write(iunit,'(a)') '<?xml version="1.0"?>'
  write(iunit,'(a)') '<VOTABLE version="1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
  write(iunit,'(a)')  'xmlns="http://www.ivoa.net/xml/VOTable/v1.2" xmlns:stc="http://www.ivoa.net/xml/STC/v1.30" >'  
  write(iunit,'(a)')  '<RESSOURCE name="IMPEx 2D cut">'
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
  
  

    write(iunit,'(3a)') '<DESCRIPTION> Magnetic field components </DESCRIPTION>'
    write(iunit,'(3a)') '<FIELD name="X" ID="col1" ucd="pos.cartesian.x" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C1" datatype="float" width="9"/>' 
    write(iunit,'(3a)') '<FIELD name="Y" ID="col2" ucd="pos.cartesian.y" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C2" datatype="float" width="9"/>' 
    write(iunit,'(3a)') '<FIELD name="Z" ID="col3" ucd="pos.cartesian.z" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C3" datatype="float" width="9"/>'
    write(iunit,'(3a)') '<FIELD name="Jtot" ID="col4" ucd="phys.current" ',&
        '  utype="" datatype="float" width="12" unit="mA.km-2"/>'
    write(iunit,'(3a)') '<FIELD name="Jx" ID="col5" ucd="phys.current" ',&
    	'  utype="" datatype="float" width="12" unit="mA.km-2"/>'
    write(iunit,'(3a)') '<FIELD name="Jy" ID="col6" ucd="phys.current" ',&
        '  utype="" datatype="float" width="12" unit="mA.km-2"/>'
    write(iunit,'(3a)') '<FIELD name="Jz" ID="col7" ucd="phys.current"',&
        '  utype="" datatype="float" width="12" unit="mA.km-2"/>'
 
 
 write(iunit,'(a)') '<DATA>'
  write(iunit,'(a)') '<TABLEDATA>'  
  

end subroutine header_2Dcut_J_VOTABLE         
        

!!##########################################################
!! m_IMPEX_cut/header_2Dcut_species_IMPEX
!! This routine writes the header of the xml
!! file for the 2D cut ion species IMPEX diag
subroutine header_2Dcut_species_VOTABLE(iunit,ion_label,plane,cutval,planetname,coord)
integer, intent(in) :: iunit,cutval
character(len=*),intent(in) :: ion_label,plane,planetname,coord
 
  write(iunit,'(a)') '<?xml version="1.0"?>'
  write(iunit,'(a)') '<VOTABLE version="1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
  write(iunit,'(a)')  'xmlns="http://www.ivoa.net/xml/VOTable/v1.2" xmlns:stc="http://www.ivoa.net/xml/STC/v1.30" >'  
  write(iunit,'(a)')  '<RESSOURCE name="IMPEx 2D cut">'
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

!  

    if (ion_label.ne.'The') write(iunit,'(3a)') '<DESCRIPTION> ',trim(ion_label),' Density, Velocity, Temperature </DESCRIPTION>'
    if (ion_label.eq.'The') write(iunit,'(a)') '<DESCRIPTION> Plasma Density, Velocity, Temperature </DESCRIPTION>'
    write(iunit,'(3a)') '<FIELD name="X" ID="col1" ucd="pos.cartesian.x" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C1" datatype="float" width="9"/>' 
    write(iunit,'(3a)') '<FIELD name="Y" ID="col2" ucd="pos.cartesian.y" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C2" datatype="float" width="9"/>' 
    write(iunit,'(3a)') '<FIELD name="Z" ID="col3" ucd="pos.cartesian.z" ', &
    	'  utype="stc:AstroCoords.Position3D.Value3.C3" datatype="float" width="9"/>'
    write(iunit,'(3a)') '<FIELD name="Density" ID="col4" ucd="phys.density" ',&
    	'  utype="" datatype="float" width="12" unit="cm-3"/>'
    write(iunit,'(3a)') '<FIELD name="Ux" ID="col5" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1"/>'
    write(iunit,'(3a)') '<FIELD name="Uy" ID="col6" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1"/>'
    write(iunit,'(3a)') '<FIELD name="Uz" ID="col7" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1"/>'  
    write(iunit,'(3a)') '<FIELD name="Utot" ID="col8" ucd="phys.veloc" ',&
        '  utype="" datatype="float" width="12" unit="km.s-1"/>'         
    write(iunit,'(3a)') '<FIELD name="Temperature" ID="col9" ucd="phys.temperature" ',&
        '  utype="" datatype="float" width="12" unit="eV"/>'        
    
  
  write(iunit,'(a)') '<DATA>'
  write(iunit,'(a)') '<TABLEDATA>'  
  

end subroutine header_2Dcut_species_VOTABLE 
 #endif    

end module m_IMPEX_cut
