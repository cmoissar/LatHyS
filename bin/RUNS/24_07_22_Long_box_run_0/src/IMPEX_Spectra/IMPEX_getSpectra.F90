program IMPEX_getSpectra

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
 use m_IMPEX_spectro_orbit
#endif
 
#include "q-p_common.h" 

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

 character(len=*), parameter :: version = '0.0'
 character(len=200) :: arg,read_variable="",input_filename="",output_filename="",name_file = ""
 character(len=300) :: cube_name="",resource_name=""
 character(len=8) :: RunID
 integer :: ang_dist= -1
 !character(len=32),dimension(:),allocatable :: tab_var
 character(len=100),dimension(:,:),allocatable :: In_array
 real(dp),dimension(:),allocatable :: X_MSO,Y_MSO,Z_MSO
 real(dp),dimension(:,:),allocatable :: Out_array
 integer :: ncommand,i,nb_var,ind,n_line,n_col=-1,startcol=-1,l
 integer :: stId,ncid,n_ene,len_file
 real(dp),dimension(:),allocatable :: EnergyRange
 real(dp) :: t1,t2,clockangle=2000,unit_traj=1.
 logical :: file_e
  __WRT_DEBUG_IN("IMPEX_getSpectra")
 print *,' I start'
! call sleep(5)

 !-- Start time
 call CPU_TIME(t1)

 
  ncommand = command_argument_count()
  
  i = 1
  do while(i<=ncommand)
   call get_command_argument(i, arg)
   select case (arg)
   case ('-c', '--cube')
    i = i+1
    call get_command_argument(i, arg)
 !   cube_name = "/data/modolo/WWW-IMPEX/"//trim(ADJUSTL(arg))
    resource_name = trim(ADJUSTL(arg))
    print *,"Resource name:",resource_name
   case ('-i', '--input')
    i = i+1
    call get_command_argument(i, arg)
 !   input_filename = "/data/modolo/WWW-IMPEX/"//trim(adjustl(arg))//".txt" 
   input_filename = trim(adjustl(arg))//".txt"    
   case ('-o', '--output')
    i = i+1
    call get_command_argument(i, arg)
 !   output_filename = "/data/modolo/WWW-IMPEX/"//trim(adjustl(arg))   
    output_filename = trim(adjustl(arg))   
   case ('-v', '--var')
    i = i+1
    call get_command_argument(i, arg)
    read_variable = trim(ADJUSTL(arg))
   case ('-n', '--ncol')
    i = i+1
    call get_command_argument(i, arg)
     read(arg,*) n_col  
   case ('-x', '--xstart')
    i = i+1
    call get_command_argument(i, arg)
     read(arg,*) startcol     
   case ('-r', '--rangle')
    i = i+1
    call get_command_argument(i, arg)
     read(arg,*) clockangle     
   case ('-V', '--version')
    print '(2a)', 'diag version ', version
    stop
   case ('-u', '--unit')
    i = i+1
    call get_command_argument(i, arg)
     read(arg,*) unit_traj  
   case ('-a', '--ang_dist')
    i = i+1
    call get_command_argument(i, arg)
     read(arg,*) ang_dist      
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
  
  ! check argument
  if (resource_name=="") then
     print *,"Error argument -c is missing"
     stop
  endif   
 ! if (read_variable=="") then
 !    print *,"Error argument -v is missing"
 !    stop
 ! endif   
  if (input_filename=="") then
     print *,"Error argument -i is missing"
     stop
  endif   
  if (output_filename=="") then
     print *,"Error argument -o is missing"
     stop
  endif   
  if (n_col<0) then
     print *,"Error argument -n is missing"
     stop
   endif  
  if (startcol<0) then
     print *,"Error argument -x is missing"
     stop
 endif   
!! From the resource name extract the run name (it will be used to reconstruct the file name of the 
!! distribution function 
len_file = len(trim(ADJUSTL(resource_name)))
print *,'length file ',len_file
cube_name = resource_name(1:len_file-15)
RunID = resource_name(len_file-10:len_file-2)
print *,'RunID :',RunID
!cube_name = resource_name(len_file-17:len_file-10)

  call wrtout(6,ch10//&
       " ==================== Treating file run  "//&
       &trim(RunID)//&
       " ====================",'PERS')

!-- Read the nb of Energy level to dimension the output array
! write(name_file,'(a,i3.3,a1,a,a)')"http://impex.latmos.ipsl.fr/Hybrid/Mars_14_03_14/Distrib_",0,'_',trim(cube_name),".nc"
write(name_file,'(a,i3.3,a1,a,a)') trim(cube_name),0,'_',trim(RunID),".nc"
!name_file = trim(cube_name)
 !--Open NetCDF file
 print *,'Filename ',trim(name_file)
    !--Inquire if the file exists
    inquire( file=trim(name_file), exist=file_e )
    !--No file: return
    if(.not.(file_e)) then
      call wrtout(6,"File: "//trim(name_file)//" does not exits","PERS")
       return
   endif
 
 stId = nf90_open(trim(name_file), nf90_nowrite, ncid)
 call test_cdf(stId) 
 call get_simple_variable_cdf(ncid,"nEnergy" ,n_ene )
 allocate(EnergyRange(n_ene));	EnergyRange(:) = zero
 call get_simple_variable_cdf(ncid,"Energy" ,EnergyRange )
 
 stId = nf90_close(ncid);  call test_cdf(stId)
print *,'Closing Filename ',trim(name_file)
 
 
 ! Read input file size
 call read_file_size(input_filename,n_line)
 
 ! creation of array of position and time
 allocate(X_MSO(n_line), Y_MSO(n_line), Z_MSO(n_line))
 X_MSO(:) = 0.;	Y_MSO(:) = 0.;	Z_MSO(:) = 0.
 
 allocate(In_array(n_line,n_col));	In_array(:,:) = ""
 allocate(Out_array(n_line,n_ene));	Out_array(:,:) = zero
 
 call read_input_file(input_filename,n_col,startcol,X_MSO,Y_MSO,Z_MSO,In_array,unit_traj)
 
  ! do l=1,n_col
  ! print *,trim(In_array(n_line,l))
  ! enddo

 call compute_spectro_orbit(cube_name,RunID,X_MSO,Y_MSO,Z_MSO,Out_array)
 
 call save_spectro_orbit(output_filename,X_MSO,Y_MSO,Z_MSO,Out_array,In_array,n_ene,EnergyRange,n_col,startcol)    
    
 if (ang_dist /=-1) call save_angular_distribution(cube_name, RunID, X_MSO(1),Y_MSO(1),Z_MSO(1)) 

 deallocate(EnergyRange)
 deallocate(X_MSO,Y_MSO,Z_MSO,In_array,Out_array)

 !-- Finish time
 call CPU_TIME(t2)
 
 
 print '("CPU Time = ",f12.3," seconds.")',t2-t1
 __WRT_DEBUG_OUT("IMPEX_getSpectra")

contains
 !!=============================================================
 !!routine: mpi_reassemble_file/print_help
 !!
 !! FUNCTION
 !!  Print help menu
 !!   
 subroutine print_help()
  print '(a)', 'usage: interpole_spectra [OPTIONS]'
  print '(a)', ''
  print '(a)', ''
  print '(a)', 'interpole spectra arguments/options:'
  print '(a)', ''
  print '(a)', "  -V, --version       print version information and exit"
  print '(a)', "  -h, --help          print usage information and exit"
  print '(a)', "  -c, --cube          resource_ID of the run [required]"
  print '(a)', "  -v, --var           Variable name requested [required]"
  print '(a)', "  -i, --input         Input file with position (and Time) [required]"
  print '(a)', "  -o, --output        Output file with position (and Time) and Value requested [required]"
  print '(a)', "  -n, --ncol          Nb of column of input file [required]"
  print '(a)', "  -x, --xtart         Index of X poistion column [required]"
  print '(a)', "  -r, --rangle        Clock Angle incident plasma [required]"
  print '(a)', "  -u, --unit          unit value of S/C position (km/Rp)[required]"
  print '(a)', "  -a, --ang_distrib   1 if angular distribution requested [optional]"

 end subroutine print_help
 
 

end program IMPEX_getSpectra




