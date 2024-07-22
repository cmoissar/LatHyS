program mpi_IMPEX_diag

 use defs_basis
 use defs_arr3Dtype
 use defs_variable
 use m_writeout
 use m_logo,only        : print_date
 use defs_tregister
 use defs_parametre,only : fildat
#ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
 use m_IMPEX_cut
 use m_IMPEX_moment_orbit
 use m_IMPEX_spectro_orbit
 use IMPEXTreeXML_generator
#endif
 

 implicit none

 !********************************************************************
 ! Auteur				:	 RModolo SHess
 ! Date					:	 03/04/12
 ! Instituion				:	LATMOS/UVSQ/IPSL
 !						Quartier des Garennes
 ! 						11 bd d'alembert
 !						78280 guyancourt
 ! Dernière modification		:	03/04/12	
 ! Résumé	
 ! Ce programme lit des fichiers reconstitué (Magw, Dnw, ...)
 ! et permet d'établir des fichiers ascii ou netcdf dans 3 plans de coupes
 ! (XY), (XZ) et (YZ)
 !********************************************************************


 integer :: ireg,ncommand,i,factor
 real(dp) :: dx
 character(len=20) :: run_name,time
 character(len=50) :: traj_name=''
 character(len=50) :: spacecraft_name=''
 character(len=5)  :: time_read=''
 character(len=32) :: arg
 integer,dimension(:),allocatable :: YYYY,MM,DD,hh,minut
 real(dp),dimension(:),allocatable :: ss
 integer :: cut2D
 character(len=*), parameter :: version = '1.0'
 real(dp),dimension(:),allocatable :: X_traj,Y_traj,Z_traj 
 integer :: nbline
 real(dp) :: t1,t2

 !-- Start time
 call CPU_TIME(t1)

 !--Insert the date for which the diagnostic has to be reassmbled
 fildat = print_date()

 ncommand = command_argument_count()
 
 ! default values
 cut2D = 0

 i = 1
 do while(i<=ncommand)
  call get_command_argument(i, arg)
  select case (arg)
  case ('-d', '--date')
   i = i+1
   call get_command_argument(i, arg)
   print *,"DATE=",trim(arg)
   fildat = trim(arg)
  case ("-t","--time")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) time_read
  case ("-traj","--trajectory")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) traj_name   
  case ("-spacecraft","--spacecraft_name")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) spacecraft_name   
  case ("-2D","--2D_cut")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) cut2D      
  case ('-v', '--version')
   print '(2a)', 'diag version ', version
   stop
  case ('-h', '--help')
   call print_help()
   stop
  case default
   print '(a,a,/)', 'Unrecognized command-line option: ', arg
   call print_help()
   stop
  end select

  i = i+1
 end do

 !--Allocate and set variables for time diagnostics
 call set_tregister(r_tm,0)

  !--Factor to put any file name >1
  factor = 1
  do while(int(factor*r_tm%time(1))<1) 
   factor = factor *10
  end do

  if (time_read =='') then  
  !--Time of Diagnostic 
    write(time,'(a,i5.5)')"t",int(factor*r_tm%time(r_tm%nreg))
  else   
      write(time,'(a,a)')"t",time_read
  endif    

  !--Create the base for the Names of file
  run_name = trim(fildat)//"_"//time

  call wrtout(6,ch10//&
       " ==================== Treating file run  "//&
       &trim(run_name)//&
       " ====================",'PERS')
  !call extract_3D(run_name)

  if (cut2D == 1) then
    call extract_2D(run_name)
  endif
  
  if (traj_name /= '') then
    call read_file_size(traj_name,nbline)
    allocate(X_traj(nbline),Y_traj(nbline),Z_traj(nbline))
    X_traj = zero;	Y_traj = zero;	Z_traj = zero
    
    allocate(YYYY(nbline),MM(nbline),DD(nbline),hh(nbline))
    allocate(minut(nbline),ss(nbline))
    YYYY = zero;	MM = zero;	DD = zero
    hh = zero;		minut = zero;	ss = zero
    
    
    call read_traj_file(traj_name,X_traj,Y_traj,Z_traj, &
    		YYYY,MM,DD,hh,minut,ss)
    	
 ! if (trim(spacecraft_name) == 'MEX') then
  !do i=0,0  		
  !  print *,'Specacraft call :',spacecraft_name
  !  call extract_1D(traj_name,run_name,X_traj,Y_traj,Z_traj, &
  !  		YYYY,MM,DD,hh,minut,ss,i*90.,spacecraft_name)
  !enddo		
 ! else
 !   print *,'Specacraft call :',spacecraft_name
 !   call extract_1D(traj_name,run_name,X_traj,Y_traj,Z_traj, &
 !   		YYYY,MM,DD,hh,minut,ss,-1.,spacecraft_name)
 ! endif  		
  
 
 
 call read_traj_file(traj_name,X_traj,Y_traj,Z_traj, &
    		YYYY,MM,DD,hh,minut,ss)
    		
  do i=0,0  		
    call extract_spectro(traj_name,run_name,X_traj,Y_traj,Z_traj, &
    		YYYY,MM,DD,hh,minut,ss,0.,spacecraft_name)
  enddo		
  
    deallocate(X_traj,Y_traj,Z_traj)
    deallocate(YYYY,MM,DD,hh,minut,ss)
  endif  
  
  

 !--Clean variables for time diagnostics
 call clean_tregister(r_tm)
 
 !-- Finish time
 call CPU_TIME(t2)
  print '("CPU Time = ",f12.3," seconds.")',t2-t1


contains
 !!=============================================================
 !!routine: mpi_reassemble_file/print_help
 !!
 !! FUNCTION
 !!  Print help menu
 !!   
 subroutine print_help()
  print '(a)', 'usage: diag [OPTIONS]'
  print '(a)', ''
  print '(a)', 'Without further options, IMPEX diag create ascii files in XY, XZ and YZ planes.'
  print '(a)', ''
  print '(a)', 'diag options:'
  print '(a)', ''
  print '(a)', "  -v, --version       print version information and exit"
  print '(a)', "  -h, --help          print usage information and exit"
  print '(a)', "  -d, --date          use the date in the format dd_mm_yy"
  print '(a)', "  -t, --time          set diagnostic time in the format XXXXX"
  print '(a)', "  -traj,--trajectory  set trajectory name (less than 50 char)"
  print '(a)', "  -spacecraft,--spacecraft_name  set spacecraft name (less than 50 char)"
  print '(a)', "  -2D, --2D_cut       followed by 1 : extract2D map in files"
 end subroutine print_help
 
 subroutine extract_3D(run_name)
   character(len=20),intent(in) :: run_name
   ! Local varaibales
    character(len=50) :: write_name
   real(dp),dimension(:,:,:),allocatable :: A,Ax,Ay,Az,B
   real(dp),dimension(:),allocatable :: A_1D,Ax_1D,Ay_1D,Az_1D,B_1D
   character(len=64) :: var_name1,var_name2,var_name3,var_name4,var_name5
   character(len=4)  :: prefix
   integer :: ncm_tot(3)
   real(dp)  :: centr(3),radius,gstep(3),box_size(6),orient(3),c_wpi
   character(len=8) :: planet
   integer :: nb_point,nb_spe,ii,iunit,sgn
    character(len=5) :: coord 
    character(len=10),dimension(:),allocatable:: ion_label
    character(len=50) :: diag_name_the
    character(len=50),dimension(:),allocatable:: diag_name
   character(len=600)::numdatid
 
   
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
   
   
   
   !--Write tree.xml metadata
   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.true.,.false.,.false.,.false.,'',0.,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
   call write_granule_3D_XML(prefix,run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
   

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
   
    !-- Read Elew file
       prefix = "Elew"
       var_name1 = "Ex"
       var_name2 = "Ey"
       var_name3 = "Ez"
       
      !--Write tree.xml metadata
      call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.true.,.false.,.false.,.false.,'',0.,&
   	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
      call write_granule_3D_XML(prefix,run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
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
   prefix = trim(ion_label(ii))
   !--Write tree.xml metadata
   call write_numdat_XML(fildat,planet,prefix(1:3),.false.,.true.,.false.,.false.,.false.,'',0.,&
	YYYY,MM,DD,hh,minut,ss,numdatid,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),orient)
!	minval(X_traj),maxval(X_traj),minval(Y_traj),maxval(Y_traj),minval(Z_traj),minval(Z_traj))   
   call write_granule_3D_XML(prefix,run_name,box_size(1),box_size(2),box_size(3),box_size(4),box_size(5),box_size(6),numdatid)
   enddo   
end subroutine extract_3D

end program mpi_IMPEX_diag




