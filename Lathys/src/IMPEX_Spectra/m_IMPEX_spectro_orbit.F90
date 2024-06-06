!!===================================================
!!===================================================
module m_IMPEX_spectro_orbit


 use defs_basis
 use defs_variable
 use defs_grid
 use defs_arr3Dtype
 use defs_particletype
 use defs_species
 use defs_parametre,only : fildat
 use m_writeout
 !use m_VO
 !use ImpexTreeXML_generator
#ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
#endif
#include "q-p_common.h"

 implicit none
 private 

#ifdef HAVE_NETCDF
 public ::		  &
      read_input_file, &
      read_file_size, &
      create_file_name,	  &
      save_spectro_orbit,  &
      compute_spectro_orbit, &
      save_angular_distribution
      !extract_XY,	  &
      !extract_XZ
      !create_file_diag,   &
      !read_moment_species_cdf,&
      !read_field_cdf_1_array


contains
 !********************************************************************
 ! Auteur				:	 RModolo, SHess
 ! Date					:	 05/04/12
 ! Institution				:	LATMOS/CNRS/IPSL
 ! Derniere modification		:	10/04/12	
 ! Resume	
 ! 
 !********************************************************************
!!=====================================================  
!! IMPEX_WS/IMPEX_WS_InterpoleFields
!! This routine reads the input file (Time and position)
!! information. Outputs are S/C position in 
!! simulation coordinate system
   subroutine read_input_file(filename,n_col,startcol,X_MSO,Y_MSO,Z_MSO,In_array,unit_traj)
   character(len=50),intent(in) :: filename
   real(dp),dimension(:),intent(inout) :: X_MSO,Y_MSO,Z_MSO
   integer,intent(in) :: n_col,startcol
   real(dp),intent(in) :: unit_traj
   character(len=*),dimension(:,:),intent(inout) :: In_array
   ! Local variables
   integer :: iunit,i, n_line,j
   integer :: Reason
   character(len=32) :: read_time,msg,a,b,c
   real(dp) :: read_x,read_y,read_z
      logical :: file_e
   
__WRT_DEBUG_IN("read_input_file")
   
   iunit = 1
   n_line = size(X_MSO,1)

   !--Inquire if the file exists
   inquire( file=trim(filename), exist=file_e )
   !--No file: return
   if(.not.(file_e)) then
     call wrtout(6,"File: "//trim(filename)//" does not exits","PERS")
      return
   endif   
   
   open(UNIT = iunit, FILE = filename, STATUS = 'OLD', ACTION = 'READ', &
        FORM = 'FORMATTED')
  print *,'n_col : ',n_col

   do i=1,n_line
!   if (n_col >3) then
!     read(iunit,*) read_time,read_x,read_y,read_z
!     Time(i) = read_time
!   else
!       read(iunit,*) read_x,read_y,read_z
!   endif
     read(iunit,*) In_array(i,1:n_col)
   !  print *,'In_array ',In_array(i,1:n_col)     
     
     !In_array(i,1) =trim(a)
     !In_array(i,2) =trim(b)
     !In_array(i,3) =trim(c)
     read(In_array(i,startcol),*) read_x
     read(In_array(i,startcol+1),*) read_y
     read(In_array(i,startcol+2),*) read_z
   !  print *,'x,y,z',read_x,read_y,read_z
     X_MSO(i) = read_x*unit_traj
     Y_MSO(i) = read_y*unit_traj
     Z_MSO(i) = read_z*unit_traj

   enddo
   close(iunit)


__WRT_DEBUG_OUT("read_input_file")
   
   end subroutine read_input_file 
   
!!#############################################################
   !!m_IMPEX_moment_orbit/read_file_size
   !! This routine reads the trajectory
   !! file and send the number of lines
   subroutine read_file_size(filename,n_line)
   character(len=50),intent(in) :: filename
   integer,intent(out) :: n_line
   !integer,intent(in) :: ncol
   !Local Variables
   integer :: Reason,j,iunit
   real(dp) :: read_x,read_y,read_z
   character(len=200) :: read_time
         logical :: file_e
__WRT_DEBUG_IN("read_file_size")

   !--Inquire if the file exists
   inquire( file=trim(filename), exist=file_e )
   !--No file: return
   if(.not.(file_e)) then
     call wrtout(6,"File: "//trim(filename)//" does not exits","PERS")
      return
   endif        
      iunit = 1
      open(UNIT = iunit, FILE = filename, STATUS = 'OLD', ACTION = 'READ', &
           FORM = 'FORMATTED')
   
      j = 0
      Reason = 0
      do while (Reason==0)
        read(iunit,*,IOSTAT=Reason) read_time
        !  print *, read_line,read_time,read_d,read_x,read_y,read_z
        j = j+1
      enddo
      n_line = j-1   
   close(iunit)
__WRT_DEBUG_OUT("read_file_size")   
   end subroutine read_file_size   
   
 
  !!=============================================================
  !!subroutine: m_IMPEX_spectro_orbit/create_file_name
  !! FUNCTION 
  !!  Create the name of the file containing particles information
  !! INPUT
  !!  filwrt=suffix containgt data
  !!  me=processus identity
  !! OUTPUT
  !!  name_file=name of the file where particles will be recorded
  subroutine create_file_name(name_file,filwrt,me)
 
   integer,intent(in) :: me
   character(len=40),intent(out) :: name_file 
   character(len=*),intent(in) :: filwrt 
   
   write(name_file,'(a3,i3.3,a1,a)')"Distrib_",me,'_',trim(filwrt)
 #ifdef HAVE_NETCDF
   name_file = trim(name_file)//".nc"
 #endif
 
 end subroutine create_file_name
 
 !!#####################################################
 !! IMPEX_Spectra/save_angular_distribution
 !! This routine extracts angular information for a point in space
 !!  and save a netcdf file
 
 subroutine save_angular_distribution(cube_name,RunID,X,Y,Z)
 character(len=*),intent(in) :: cube_name,RunID
 real(dp),intent(in) :: X,Y,Z
    character(len=200) ::name_file=""
   character(len=5) :: coord
   integer :: sgn,nproc,proc,varid(100)
   real(dp)  :: centr(3),radius,gstep(3),c_wpi,smin(3),smax(3),X_tmp,Y_tmp,Z_tmp
   integer :: ii, i,j,k,n_E,n_T,n_P,stId,ncid
   character(len=8) :: planet
   real(dp),allocatable,dimension(:) :: Energy,Thetaval,Phival
    integer :: varfieldId, numDims,var_id
    real(dp),dimension(:,:,:),allocatable :: distrib_tmp
   integer,dimension(nf90_max_var_dims) :: DimfieldId  
      integer :: dim_distrib(4)
 
 __WRT_DEBUG_IN("save_angular_resolution") 
   ! coordinate transformation for each object
   select case (trim(planet))
   case ("mars","mercure")
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
    
!-- Read the Ditrib file of proc 000 to extract information about Energy level,  theta and phi angles as well as information about the 
! processus domain box, planet position and so on.
  ! write(name_file,'(a8,i3.3,a1,a,a)')"Distrib_",0,'_',trim(cube_name),".nc"
  write(name_file,'(a,i3.3,a1,a,a)') trim(cube_name),0,'_',trim(RunID),".nc" 
   !--Open NetCDF file
   stId = nf90_open(trim(name_file), nf90_nowrite, ncid)
   call test_cdf(stId) 
   call get_simple_variable_cdf(ncid,"nEnergy" ,n_E )
   allocate(Energy(n_E));	Energy(:) = zero
   call get_simple_variable_cdf(ncid,"Energy" ,Energy )
   call test_cdf(stId) 
   call get_simple_variable_cdf(ncid,"nTheta" ,n_T )
   call test_cdf(stId)
   allocate(Thetaval(n_T));	Thetaval(:) = zero
   call get_simple_variable_cdf(ncid,"Theta" ,Thetaval )
   call test_cdf(stId)  
   call get_simple_variable_cdf(ncid,"nPhi" ,n_P )
   call test_cdf(stId)
   allocate(Phival(n_P));	Phival(:) = zero
   call get_simple_variable_cdf(ncid,"Phi" ,Phival )
   call test_cdf(stId)
   call get_simple_variable_cdf(ncid,"r_planet" ,radius )
   call test_cdf(stId)
   call get_simple_variable_cdf(ncid,"s_centr" ,centr )
   call test_cdf(stId)  
  ! call get_simple_variable_cdf(ncid,"gstep" ,gstep )
  ! call test_cdf(stId)  
    gstep(:) = 0.705
   call get_simple_variable_cdf(ncid,"phys_length" ,c_wpi )
   call test_cdf(stId)  
   call get_simple_variable_cdf(ncid,"nproc" ,nproc )
   call test_cdf(stId)     
   stId = nf90_close(ncid)   
  
     X_tmp = (sgn*X/c_wpi + centr(1))  ! Traj in c_wpi unit
     Y_tmp = (sgn*Y/c_wpi + centr(2))	! Traj in c_wpi unit
     Z_tmp =  (Z/c_wpi + centr(3))	! Traj in c_wpi unit
!-- alocating array for distribution function
  allocate(distrib_tmp(n_E,n_T,n_P));	distrib_tmp(:,:,:) = zero  
 print *,' ========= Reading Process file / Angular Distribution =============' 
    do proc=1,nproc
      print *,' Processing file # ',proc
    !do proc=1,1
      write(name_file,'(a8,i3.3,a1,a,a)')"Distrib_",proc-1,'_',trim(cube_name),".nc"
      !--Open NetCDF file
      stId = nf90_open(trim(name_file), nf90_nowrite, ncid)
      call test_cdf(stId)      
      call get_simple_variable_cdf(ncid,"s_min_loc" ,smin )
      call test_cdf(stId)  
      call get_simple_variable_cdf(ncid,"s_max_loc" ,smax )
      call test_cdf(stId)  
    !--Get number of point for fields in X,Y,Z
    stId = nf90_inq_varId(ncid,"Distribution_function",varfieldId)
    stId = nf90_Inquire_Variable(ncid,varfieldId,dimids = DimfieldId)
    stId = nf90_Inquire_Dimension(ncid,DimfieldId(1), len=ncm(1))
    stId = nf90_Inquire_Dimension(ncid,DimfieldId(2), len=ncm(2))
    stId = nf90_Inquire_Dimension(ncid,DimfieldId(3), len=ncm(3))  
 
      !-- check if part of the trajectory is in the subdomain
      !-- check if one point of X traj is in the box
 
      if ((X_tmp >= smin(1)).and.((X_tmp <= smax(1)))) then
      if ((Y_tmp >= smin(2)).and.((Y_tmp <= smax(2)))) then
      if ((Z_tmp >= smin(3)).and.((Z_tmp <= smax(3)))) then
 !-- find indexes in the trajectory which are in the subdomain
    !  print *,'Part of the traj inside the subdomain'
    !  print *,'Dimension ',ncm(1),ncm(2),ncm(3),n_E,n_T,n_P
        i = int((X_tmp-smin(1))/(2.*gstep(1)))+1
        j = int((Y_tmp-smin(2))/(2.*gstep(2)))+1
        k = int((Z_tmp-smin(3))/(2.*gstep(3)))+1
       ! print *,'ijk ', i,j,k
         if ((i>0) .and. (i<ncm(1))) then
           if ((j>0) .and. (j<ncm(2))) then
             if ((k>0) .and. (k<ncm(3))) then  
                 print *, 'Distribution found #proc :',proc
                distrib_tmp(:,:,:) = zero
                stId = nf90_inq_varid(ncid, "Distribution_function", var_id)
                call test_cdf(stId)
                stId = nf90_get_var(ncid, var_id, distrib_tmp(:, :, :), &
                        start = (/ i, j, k, 1, 1, 1 /),     &
                       count = (/ 1, 1, 1, n_E, n_T, n_P /))
 	     endif
 	   endif
 	 endif
 	 endif
 	 endif
 	 endif
    enddo
   
   !-- Saving the information in a netcdf file
    !--Create the file name
   write(name_file,'(a,a,a)')"Angular_Distribution_",trim(cube_name),".nc" 

   stId = nf90_create(name_file , nf90_clobber, ncid)
   call test_cdf(stId)   
  
   stId = nf90_def_dim(ncid, "size_Energy"          , n_E    , dim_distrib(1))
   call test_cdf(stId)
   stId = nf90_def_dim(ncid, "size_Theta"          , n_T    , dim_distrib(2))
   call test_cdf(stId)
   stId = nf90_def_dim(ncid, "size_Phi"          , n_P    , dim_distrib(3))
   call test_cdf(stId)   
   stId = nf90_def_dim(ncid, "dim_scalar"          , 1    , dim_distrib(4))
   call test_cdf(stId)    
   
      ii = 1

   stId = nf90_def_var(ncid, "EnergyDim",    nf90_int,dim_distrib(4),varid(ii))
   call test_cdf(stId); ii = ii+1
   stId = nf90_def_var(ncid, "ThetDim",    nf90_int,dim_distrib(4),varid(ii))
   call test_cdf(stId); ii = ii+1
   stId = nf90_def_var(ncid, "PhiDim",    nf90_int,dim_distrib(4),varid(ii))
   call test_cdf(stId); ii = ii+1
   stId = nf90_def_var(ncid, "Energy",    QP_NF90_DP,dim_distrib(1),varid(ii))
   call test_cdf(stId); ii = ii   +1
   stId = nf90_def_var(ncid, "Theta",    QP_NF90_DP,dim_distrib(2),varid(ii))
   call test_cdf(stId); ii = ii +1
   stId = nf90_def_var(ncid, "Phi",    QP_NF90_DP,dim_distrib(3),varid(ii))
   call test_cdf(stId); ii = ii  +1    
   stId = nf90_def_var(ncid, "Distribution",    QP_NF90_DP,dim_distrib(1:3),varid(ii))
   call test_cdf(stId); ii = ii+1
   
    !--Switch to write mode
      stId = nf90_enddef(ncid); call test_cdf(stId)
    
   ii = 1 
   stId = nf90_put_var(ncid, varid(ii), n_E)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), n_T)
   call test_cdf(stId); ii = ii+1 
   stId = nf90_put_var(ncid, varid(ii), n_P)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), Energy)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), Thetaval)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), Phival)
   call test_cdf(stId); ii = ii+1    
   stId = nf90_put_var(ncid, varid(ii), distrib_tmp)
   call test_cdf(stId); ii = ii+1    
     !--Close the file
   stId = nf90_close(ncid); call test_cdf(stId) 
   
   
   deallocate(Energy,Thetaval,Phival,distrib_tmp)
   
   
 
 __WRT_DEBUG_OUT("save_angular_resolution") 
 end subroutine save_angular_distribution
 
 !!#####################################################
 !! IMPEX_Spectra/m_IMPEX_spectro_orbit.F90
 !! This routine extracts the spectra from the Distribution files
 !! along the S/C trajectory
 
 subroutine compute_spectro_orbit(cube_name,RunID,X_traj,Y_traj,Z_traj,Out_array)
 character(len=*),intent(in) :: cube_name,RunID
  real(dp),dimension(:),intent(in) ::X_traj,Y_traj,Z_traj
   ! real(dp),intent(in) :: clockangle
 ! integer,intent(in) :: nb_var
  real(dp),dimension(:,:),intent(inout) :: Out_array
  ! Local Variables
   character(len=200) ::name_file=""
  character(len=32) :: var_name=""
  character(len=5) :: coord
  integer :: ncm(3),var,line,n_line,sgn,nproc,proc
  real(dp)  :: centr(3),radius,gstep(3),c_wpi,smin(3),smax(3) ,sum1
  real(dp) :: xa,ya,za,xf,yf,zf,w1,w2,w3,w4,w5,w6,w7,w8
  integer :: i,j,k,n_E,n_T,n_P,stId,ncid,nb_point,l,ie
  character(len=8) :: planet
  real(dp),allocatable,dimension(:) :: Energy,Thetaval,Phival
  real(dp),allocatable,dimension(:) :: X_traj_tmp,Y_traj_tmp,Z_traj_tmp
  real(dp),allocatable,dimension(:,:,:,:,:,:) :: distrib
   integer :: varfieldId, numDims,var_id
   real(dp),dimension(:,:,:),allocatable :: distrib_tmp
   integer,dimension(nf90_max_var_dims) :: DimfieldId  
  
__WRT_DEBUG_IN("compute_spectro_orbit")  

print *,' Max & Min location of the trajectory after rotation'
 !  print *, 'angle',angle
   print *, 'X ',maxval(X_traj),minval(X_traj)
   print *, 'Y ',maxval(Y_traj),minval(Y_traj)
   print *, 'Z ',maxval(Z_traj),minval(Z_traj)
       
   ! coordinate transformation for each object
   select case (trim(planet))
   case ("mars","mercure")
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
    
!-- Read the Ditrib file of proc 000 to extract information about Energy level,  theta and phi angles as well as information about the 
! processus domain box, planet position and so on.
!   write(name_file,'(a8,i3.3,a1,a,a)')"Distrib_",0,'_',trim(cube_name),".nc"
write(name_file,'(a,i3.3,a1,a,a)') trim(cube_name),0,'_',trim(RunID),".nc"
   !--Open NetCDF file
   stId = nf90_open(trim(name_file), nf90_nowrite, ncid)
   call test_cdf(stId) 
   call get_simple_variable_cdf(ncid,"nEnergy" ,n_E )
   allocate(Energy(n_E));	Energy(:) = zero
   call get_simple_variable_cdf(ncid,"Energy" ,Energy )
   call test_cdf(stId) 
   call get_simple_variable_cdf(ncid,"nTheta" ,n_T )
   call test_cdf(stId)
   allocate(Thetaval(n_T));	Thetaval(:) = zero
   call get_simple_variable_cdf(ncid,"Theta" ,Thetaval )
   call test_cdf(stId)  
   call get_simple_variable_cdf(ncid,"nPhi" ,n_P )
   call test_cdf(stId)
   allocate(Phival(n_P));	Phival(:) = zero
   call get_simple_variable_cdf(ncid,"Phi" ,Phival )
   call test_cdf(stId)
   call get_simple_variable_cdf(ncid,"r_planet" ,radius )
   call test_cdf(stId)
   call get_simple_variable_cdf(ncid,"s_centr" ,centr )
   call test_cdf(stId)  
  ! call get_simple_variable_cdf(ncid,"gstep" ,gstep )
  ! call test_cdf(stId)  
    gstep(:) = 0.705
   call get_simple_variable_cdf(ncid,"phys_length" ,c_wpi )
   call test_cdf(stId)  
   call get_simple_variable_cdf(ncid,"nproc" ,nproc )
   call test_cdf(stId)     
   stId = nf90_close(ncid)   
  
!-- alocating array for distribution function
  allocate(distrib_tmp(n_E,n_T,n_P));	distrib_tmp(:,:,:) = zero   
   

   
   !--Convert MSO trajectory in simulation unit and coordinate
   nb_point = size(X_traj,1)
   allocate(X_traj_tmp(nb_point));	X_traj_tmp(:)=zero
   allocate(Y_traj_tmp(nb_point));	Y_traj_tmp(:)=zero
   allocate(Z_traj_tmp(nb_point));	Z_traj_tmp(:)=zero
   X_traj_tmp = (sgn*X_traj/c_wpi + centr(1))  ! Traj in c_wpi unit
   Y_traj_tmp = (sgn*Y_traj/c_wpi + centr(2))	! Traj in c_wpi unit
   Z_traj_tmp =  (Z_traj/c_wpi + centr(3))	! Traj in c_wpi unit
   print *,'centr ',centr
   print *,'c_wpi ',c_wpi
   print *,'sgn ',sgn
   print *,'X_traj , X_traj_tmp ',X_traj,X_traj_tmp
   print *,'Y_traj , Y_traj_tmp ',Y_traj,Y_traj_tmp
   print *,'Z_traj , Z_traj_tmp ',Z_traj,Z_traj_tmp
   
   print *,' ========= Reading Process file =============' 
   do proc=1,nproc
     print *,' Processing file # ',proc
   !do proc=1,1
   !  write(name_file,'(a8,i3.3,a1,a,a)')"Distrib_",proc-1,'_',trim(cube_name),".nc"
   write(name_file,'(a,i3.3,a1,a,a)') trim(cube_name),proc-1,'_',trim(RunID),".nc"
     !--Open NetCDF file
     stId = nf90_open(trim(name_file), nf90_nowrite, ncid)
     call test_cdf(stId)      
     call get_simple_variable_cdf(ncid,"s_min_loc" ,smin )
     call test_cdf(stId)  
     call get_simple_variable_cdf(ncid,"s_max_loc" ,smax )
     call test_cdf(stId)  
   !--Get number of point for fields in X,Y,Z
   stId = nf90_inq_varId(ncid,"Distribution_function",varfieldId)
   stId = nf90_Inquire_Variable(ncid,varfieldId,dimids = DimfieldId)
   stId = nf90_Inquire_Dimension(ncid,DimfieldId(1), len=ncm(1))
   stId = nf90_Inquire_Dimension(ncid,DimfieldId(2), len=ncm(2))
   stId = nf90_Inquire_Dimension(ncid,DimfieldId(3), len=ncm(3))  

   !print *,' smin ',smin
   !print *,' smax ',smax
     
     !call get_simple_variable_cdf(ncid,"ncm" ,ncm )
     !call test_cdf(stId)
     !-- check if part of the trajectory is in the subdomain
     !-- check if one point of X traj is in the box

     if (any(X_traj_tmp >= smin(1)).and.(any(X_traj_tmp <= smax(1)))) then
     if (any(Y_traj_tmp >= smin(2)).and.(any(Y_traj_tmp <= smax(2)))) then
     if (any(Z_traj_tmp >= smin(3)).and.(any(Z_traj_tmp <= smax(3)))) then
!-- find indexes in the trajectory which are in the subdomain
   !  print *,'Part of the traj inside the subdomain'
   !  print *,'Dimension ',ncm(1),ncm(2),ncm(3),n_E,n_T,n_P
   do l=1,nb_point
       i = int((X_traj_tmp(l)-smin(1))/(2.*gstep(1)))+1
       j = int((Y_traj_tmp(l)-smin(2))/(2.*gstep(2)))+1
       k = int((Z_traj_tmp(l)-smin(3))/(2.*gstep(3)))+1
      ! print *,'ijk ', i,j,k
        if ((i>0) .and. (i<ncm(1))) then
          if ((j>0) .and. (j<ncm(2))) then
            if ((k>0) .and. (k<ncm(3))) then      
               distrib_tmp(:,:,:) = zero
               stId = nf90_inq_varid(ncid, "Distribution_function", var_id)
               call test_cdf(stId)
               stId = nf90_get_var(ncid, var_id, distrib_tmp(:, :, :), &
                       start = (/ i, j, k, 1, 1, 1 /),     &
                       count = (/ 1, 1, 1, n_E, n_T, n_P /))
               do ie = 1,n_E
                 sum1 = sum(distrib_tmp(ie,:,:))
                 Out_array(l,ie) = Out_array(l,ie) + sum1
                enddo                 
            endif
          endif
        endif  
   enddo
 
!     !-- then we can allocate the  distribution function
!     allocate(distrib(ncm(1),ncm(2),ncm(3),n_E,n_T,n_P));	distrib(:,:,:,:,:,:) = zero
!     call get_simple_variable_cdf(ncid,"Distribution_function" ,distrib )
!     call test_cdf(stId)    
!     do l=1,nb_point
!   !    print *,' Treating line ',l
!       i = int((X_traj_tmp(l)-smin(1))/(2.*gstep(1)))+1
!       j = int((Y_traj_tmp(2)-smin(2))/(2.*gstep(2)))+1
!       k = int((Z_traj_tmp(3)-smin(3))/(2.*gstep(3)))+1
!        if ((i>0) .and. (i<ncm(1))) then
!          if ((j>0) .and. (j<ncm(2))) then
!            if ((k>0) .and. (k<ncm(3))) then
!        
!        	xa = ((X_traj_tmp(1)-smin(1))/(2.*gstep(1)))+1 - float(i)
!        	ya = ((Y_traj_tmp(2)-smin(2))/(2.*gstep(2)))+1 - float(j)
!        	za = ((Z_traj_tmp(3)-smin(3))/(2.*gstep(3)))+1 - float(k)
!        	xf = 1.0 - xa
!        	yf = 1.0 - ya
!        	zf = 1.0 - za
!        	w1 = xa*ya*za   ! (i+1,j+1,k+1)
!        	w2 = xf*ya*za   ! (i  ,j+1,k+1)
!        	w3 = xa*yf*za   ! (i+1,j  ,k+1)
!        	w4 = xf*yf*za   ! (i  ,j  ,k+1)
!        	w5 = xa*ya*zf   ! (i+1,j+1,k  )
!        	w6 = xf*ya*zf   ! (i  ,j+1,k  )
!        	w7 = xa*yf*zf   ! (i+1,j  ,k  )
!        	w8 = xf*yf*zf   ! (i  ,j  ,k  )
!                do ie = 1,n_E
!                  sum1 = sum(distrib(i,j,k,ie,:,:))
!                  Out_array(l,ie) = Out_array(l,ie) + sum1
!                enddo
!                endif
!              endif
!            endif
!          enddo
!          deallocate(distrib)
         endif
         endif
         endif
        ! close the file
        stId = nf90_close(ncid)
       			       
    
   enddo

deallocate(distrib_tmp)
deallocate(Energy,Thetaval,Phival)
deallocate(X_traj_tmp,Y_traj_tmp,Z_traj_tmp)
__WRT_DEBUG_OUT("compute_spectro_orbit")  
 end subroutine compute_spectro_orbit
 
 
 
    !!###################################################
    !! m_IMPEX_diag/save_spectro_orbit
    !! This routines dumps into a file
    !! the extracted spectrogram along a trajectory
    subroutine save_spectro_orbit(output_filename,X_MSO,Y_MSO,Z_MSO,Out_array,In_array,n_ene,EnergyRange,ncol,startcol)
    character(len=*),intent(in) ::output_filename
    real(dp),dimension(:),intent(in) :: X_MSO,Y_MSO,Z_MSO
    real(dp),dimension(:,:),intent(in) :: Out_array
    character(len=*),dimension(:,:),intent(in) :: In_array
    integer,intent(in) :: n_ene,ncol,startcol
    real(dp),dimension(:),intent(in) :: EnergyRange
    character(len=32) :: val_ucd,unit,varname
    integer :: io,n,l
    character(len=300) :: readfile
    character(len=2000) :: msg

    ! Local variables
    integer :: iunit,i,nb
    
    __WRT_DEBUG_IN("save_spectro_orbit") 
        
        !-- Get length of the array
    nb = size(X_MSO,1)

 
    
    ! in VOTale format
    iunit = 1
     open(UNIT = iunit,FILE = output_filename, FORM = 'FORMATTED', STATUS = 'UNKNOWN', &
    	ACTION = 'READWRITE')
    io=1
    do while (io >=0)
      read(iunit,*,IOSTAT=io) readfile
   !   print *,io
    enddo
 
    ! completing the Header
    unit="Ions.m-2.s-1.str-1.eV-1"
    val_ucd = "phys.flux"
    varname = "ParticleFlux"

    if (ncol+1 <10) then
      write(iunit,'(3a,i1,5a,i2,a)') '<FIELD name="',trim(varname),'" ID="col',ncol+1,'" ucd="',trim(val_ucd),&
        '" utype="" datatype="float"  unit="',trim(unit),'" arraysize="',n_ene,'" />'
     else
       write(iunit,'(3a,i2,5a,i2,a)') '<FIELD name="',trim(varname),'" ID="col',ncol+1,'" ucd="',trim(val_ucd),&
        '" utype="" datatype="float"  unit="',trim(unit),'" arraysize="',n_ene,'" />'
     endif   
  write(msg,'(a,i2,a)') '<PARAM name="EnergyRange" unit="eV" ucd="instr.param" datatype="float" arraysize="',n_ene,'" value="   0.'
  do n=1,n_ene
    write(msg,'(a,f8.1)') trim(msg)//',',EnergyRange(n)
  enddo
  write(msg,'(a)') trim(msg)//'"/>'
  write(iunit,'(a)') trim(msg)
  
  write(iunit,'(a)') '<DATA>'
  write(iunit,'(a)') '<TABLEDATA>'  
  

   
     do i=1,nb
         write(iunit,'(a)') '<TR>' 
         msg =''
         do n=1,ncol     
           write(msg,'(a)') trim(msg)//'<TD> '//trim(In_array(i,n))//' </TD> '
         enddo  
        !   print *,' i ',i,trim(msg)
           write(msg,'(a)') trim(msg)//'<TD> '
         do n=1,n_ene
           write(msg,'(a,e10.3)') trim(msg),Out_array(i,n)
         enddo
          write(msg,'(a)') trim(msg)//'</TD> '
         
         write(iunit,'(a)') trim(msg)
                 
         write(iunit,'(a)') '</TR>'	
            
    enddo
  ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESOURCE>'
  write(iunit,'(a)') '</VOTABLE>'      
    close(iunit)

   __WRT_DEBUG_OUT("save_spectro_orbit_value") 

    end subroutine save_spectro_orbit 
     !!######################################################
    !! m_IMPEX_diag/header_orbit_VOTABLE
    !! This routine write the header of the 
    !! VOTABLE file
    subroutine header_spectro_VOTABLE(iunit,prefix,planetname,coord)
    character(len=*),intent(in) :: prefix,planetname,coord
    integer,intent(in) :: iunit
  write(iunit,'(a)') '<?xml version="1.0"?>'
  write(iunit,'(a)') '<VOTABLE version="1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
  write(iunit,'(a)')  'xmlns="http://www.ivoa.net/xml/VOTable/v1.2" xmlns:stc="http://www.ivoa.net/xml/STC/v1.30" >'  
  write(iunit,'(a)')  '<RESSOURCE name="IMPEx fly-through">'
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
  
 
    write(iunit,'(3a)') '<DESCRIPTION> Ion Spectrogram </DESCRIPTION>'
    write(iunit,'(3a)') '<FIELD name="Time" ID="col1" ucd="time" ref="UTC" ', &
    	'  utype="stc:AstroCoords.Time.TimeInstant.ISOTime" datatype="integer" width="4" />' 
    write(iunit,'(3a)') '<FIELD name="Counts" ID="col2" ucd="" ref=" ',trim(coord),&
        ' " utype="" datatype="float" width="10" unit="" arraysize="94"/>'        

  
  write(iunit,'(a)') '<DATA>'
  write(iunit,'(a)') '<TABLEDATA>'  
    end subroutine header_spectro_VOTABLE

 #endif
 
 end module m_IMPEX_spectro_orbit
