program IMPEX_Interpole_Fields

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
#endif
 
#include "q-p_common.h"

 implicit none

 !********************************************************************
 ! Auteur				:	 RModolo SHess
 ! Date					:	 10/04/13
 ! Instituion				:	LATMOS/UVSQ/IPSL
 !						Quartier des Garennes
 ! 						11 bd d'alembert
 !						78280 guyancourt
 ! Dernière modification		:	03/04/12	
 ! Résumé	
 ! Ce programme lit des fichiers reconstitué (Magw, Dnw, ...)
 ! et permet d'établir des fichiers VOTables avec les quantités demandées 
 ! aux points d'espace souhaité.
 !********************************************************************

 character(len=*), parameter :: version = '0.0'
 character(len=200) :: arg,read_variable="",input_filename="",output_filename=""
 character(len=300) :: cube_name=""
 character(len=32),dimension(:),allocatable :: tab_var
 character(len=200),dimension(:,:),allocatable :: In_array
 real(dp),dimension(:),allocatable :: X_MSO,Y_MSO,Z_MSO
 real(dp),dimension(:,:),allocatable :: Out_array
 integer :: ncommand,i,nb_var,ind,n_line,n_col=-1,startcol=-1,l
 real(dp) :: t1,t2,clockangle=2000,unit_traj=1.
 __WRT_DEBUG_IN("IMPEX_interpole_field") 
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
   cube_name = trim(ADJUSTL(arg))
  ! print *,"Cube ",cube_name
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
 if (cube_name=="") then
    print *,"Error argument -c is missing"
    stop
 endif   
 if (read_variable=="") then
    print *,"Error argument -v is missing"
    stop
 endif   
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
 

!  call wrtout(6,ch10//&
!       " ==================== Treating file run  "//&
!       &trim(run_name)//&
!       " ====================",'PERS')

 nb_var = count_nb_var(read_variable)
 !print *, 'nb variable ', nb_var
 
 ! allocation of variable array
 allocate(tab_var(nb_var)) ; tab_var(:) =""
 
 if (nb_var > 1) then
   do i=1,nb_var-1
     ind = INDEX(read_variable,",")
     tab_var(i) = read_variable(1:ind-1)
     read_variable = read_variable(ind+1:len(read_variable))
   enddo
 endif  
 tab_var(nb_var) = read_variable
 
 
 ! Read input file size
 call read_file_size(input_filename,n_line)
 !print *,'Le fichier comporte ',n_line,' lignes'
 
 ! creation of array of position and time
 allocate(X_MSO(n_line), Y_MSO(n_line), Z_MSO(n_line))
 X_MSO(:) = 0.;	Y_MSO(:) = 0.;	Z_MSO(:) = 0.
 
 allocate(In_array(n_line,n_col));	In_array(:,:) = ""
 
 allocate(Out_array(n_line,nb_var));	Out_array(:,:) = 0.
 
! print *, 'Requested variables are :'
! do i=1,nb_var
!   print *,tab_var(i)
! enddo
 
 
  call read_input_file(input_filename,n_col,startcol,X_MSO,Y_MSO,Z_MSO,In_array,unit_traj)
  do l=1,n_col
  print *,trim(In_array(n_line,l))
   enddo
 ! extract values

 call extract_field_value(cube_name,X_MSO,Y_MSO,Z_MSO,Out_array,nb_var,tab_var,clockangle)

 
 call save_field_value(output_filename,X_MSO,Y_MSO,Z_MSO,Out_array,In_array,nb_var,tab_var,n_col,startcol)


 deallocate(tab_var)
 deallocate(X_MSO,Y_MSO,Z_MSO,In_array,Out_array)
 !-- Finish time
 call CPU_TIME(t2)
 
__WRT_DEBUG_OUT("IMPEX_interpole_field")
  print '("CPU Time = ",f12.3," seconds.")',t2-t1
 
contains
 !!=============================================================
 !!routine: mpi_reassemble_file/print_help
 !!
 !! FUNCTION
 !!  Print help menu
 !!   
 subroutine print_help()
  print '(a)', 'usage: interpole_fields [OPTIONS]'
  print '(a)', ''
  print '(a)', ''
  print '(a)', 'interpole fields arguments/options:'
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

 end subroutine print_help
 
 !!===============================================================
 !! routine : extract_value
 !! FUNCTION
 !! extract field information from simulation files
 
 subroutine extract_field_value(cube_name,X_MSO,Y_MSO,Z_MSO,Out_array,nb_var,tab_var,clockangle)
 implicit none
 character(len=*),intent(in) :: cube_name
 real(dp),dimension(:),intent(inout) ::X_MSO,Y_MSO,Z_MSO
 character(len=*), dimension(:),intent(in) :: tab_var
 real(dp),intent(in) :: clockangle
 integer,intent(in) :: nb_var
 real(dp),dimension(:,:),intent(inout) :: Out_array
 ! Local Variables
 character(len=32) :: var_name1="",var_name2="",var_name3="",var_name4="",var_name5=""
 integer :: ncm_tot(3),ind_mag,ind_vel,ind_ele,ind_spec,var,line,n_line,ind_Jcur,ind_GPma
 real(dp)  :: centr(3),radius,gstep(3),c_wpi
 real(dp) :: xa,ya,za,xf,yf,zf,w1,w2,w3,w4,w5,w6,w7,w8
 integer :: i,j,k
 character(len=8) :: planet
 real(dp),dimension(:,:,:),allocatable :: Ax, Ay, Az,Amod,A,B
 real(dp),dimension(:),allocatable :: Ax_1D,Ay_1D,Az_1D,A_1D,B_1D
 real(dp),dimension(:),allocatable :: X0_traj,Y0_traj,Z0_traj
 
__WRT_DEBUG_IN("extract_field_value")  
 
 ind_mag = INDEX(cube_name,"Magw")
 if (ind_mag /= 0) then
    var_name1 = "Bx"; var_name2 = "By"; var_name3 = "Bz"
 endif
 
 ind_ele = INDEX(cube_name,"Elew")
  if (ind_ele /= 0) then
     var_name1 = "Ex"; var_name2 = "Ey"; var_name3 = "Ez"
  endif
  
  ind_vel = INDEX(cube_name,"Velw")
   if (ind_vel /= 0) then
      var_name1 = "Ux"; var_name2 = "Uy"; var_name3 = "Uz"
   endif
 ind_Jcur = INDEX(cube_name,"Jcur")
  if (ind_Jcur /= 0) then
     var_name1 = "Jx"; var_name2 = "Jy"; var_name3 = "Jz"
 endif

  ind_GPma = INDEX(cube_name,"GPma")
   if (ind_GPma /= 0) then
      var_name1 = "GPmagx"; var_name2 = "GPmagy"; var_name3 = "GPmagz"
 endif

 if ((ind_mag ==0).and.(ind_vel==0).and.(ind_ele==0).and.(ind_Jcur==0).and.(ind_GPma==0)) then
     var_name1 = "Density"; var_name2 = "Ux";   var_name3 = "Uy";   var_name4 = "Uz"
     var_name5 = "Temperature"
 endif    
  
 print *,'Variable',var_name1 
 
 call read_dim_field_cdf(cube_name,var_name1,ncm_tot,radius,centr,gstep,planet,c_wpi)

 allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 allocate(Ay(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 allocate(Az(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 allocate(Amod(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 Ax = zero;	Ay = zero;	Az = zero;	Amod = zero

 if ((ind_mag ==0).and.(ind_vel==0).and.(ind_ele==0)) then
    allocate(A(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
    allocate(B(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
    A = zero;	B=zero
 endif

 if ((ind_mag /=0).or.(ind_vel/=0).or.(ind_ele/=0).or.(ind_Jcur/=0).or.(ind_GPma/=0)) then
   call read_field_cdf(cube_name,var_name1,var_name2,var_name3, Ax,Ay,Az)
 else
   call read_species_cdf(cube_name,var_name1,var_name2,var_name3,var_name4,var_name5, A,Ax,Ay,Az,B)
 endif
 
 Amod = sqrt(Ax**2+Ay**2+Az**2)
 
 ! extracting value at the position
 ! 3D-linear interpolation
 n_line = size(X_MSO,1)
 allocate(Ax_1D(n_line),Ay_1D(n_line),Az_1D(n_line))
 Ax_1D(:) = zero;	Ay_1D(:) = zero; Az_1D(:) = zero
 
 
 ! rotating the trajectory if needed
 !===================================
 
 print *,'Colck angle value',clockangle
 if (clockangle /= 2000.) then
 allocate(X0_traj(n_line),Y0_traj(n_line),Z0_traj(n_line));
 X0_traj(:) = zero;	Y0_traj(:) = zero;	Z0_traj(:)=zero
 X0_traj = X_MSO; Y0_traj = Y_MSO; Z0_traj = Z_MSO
 X_MSO=X0_traj
 Y_MSO=Y0_traj*cos((clockangle+90)/180.*pi)-Z0_traj*sin((clockangle+90)/180.*pi)
 Z_MSO=Z0_traj*cos((clockangle+90)/180.*pi)+Y0_traj*sin((clockangle+90)/180.*pi)
 deallocate(X0_traj,Y0_traj,Z0_traj)
 endif


if ((ind_mag ==0).and.(ind_vel==0).and.(ind_ele==0).and.(ind_Jcur==0)) then
  allocate(A_1D(n_line),B_1D(n_line))
  A_1D = zero;	B_1D = 0
endif  
  
 
 do line=1,n_line
 
        i = int((-X_MSO(line)/c_wpi+centr(1))/gstep(1))+1
        j = int((-Y_MSO(line)/c_wpi+centr(2))/gstep(2))+1
        k = int(( Z_MSO(line)/c_wpi+centr(3))/gstep(3))+1
      
        
        if ((i>0) .and. (i<=ncm_tot(1))) then
          if ((j>0) .and. (j<=ncm_tot(2))) then
            if ((k>0) .and. (k<=ncm_tot(3))) then
        
        	xa = ((-X_MSO(line)/c_wpi+centr(1))/gstep(1))+1 - float(i)
        	ya = ((-Y_MSO(line)/c_wpi+centr(2))/gstep(2))+1 - float(j)
        	za = (( Z_MSO(line)/c_wpi+centr(3))/gstep(3))+1 - float(k)
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
          	
        	
        	Ax_1D(line) = ( w1*Ax(i+1,j+1,k+1) + &
        	     w2*Ax(i  ,j+1,k+1) + w3*Ax(i+1,j  ,k+1) + &
        	     w4*Ax(i  ,j  ,k+1) + w5*Ax(i+1,j+1,k  ) + &
        	     w6*Ax(i  ,j+1,k  ) + w7*Ax(i+1,j  ,k  ) + &
        	     w8*Ax(i  ,j  ,k  ))
 
 		Ay_1D(line) = ( w1*Ay(i+1,j+1,k+1) + &
        	     w2*Ay(i  ,j+1,k+1) + w3*Ay(i+1,j  ,k+1) + &
        	     w4*Ay(i  ,j  ,k+1) + w5*Ay(i+1,j+1,k  ) + &
        	     w6*Ay(i  ,j+1,k  ) + w7*Ay(i+1,j  ,k  ) + &
        	     w8*Ay(i  ,j  ,k  ))

  	Az_1D(line) = ( w1*Az(i+1,j+1,k+1) + &
        	     w2*Az(i  ,j+1,k+1) + w3*Az(i+1,j  ,k+1) + &
        	     w4*Az(i  ,j  ,k+1) + w5*Az(i+1,j+1,k  ) + &
        	     w6*Az(i  ,j+1,k  ) + w7*Az(i+1,j  ,k  ) + &
        	     w8*Az(i  ,j  ,k  ))

           if ((ind_mag ==0).and.(ind_vel==0).and.(ind_ele==0).and.(ind_Jcur==0)) then
                  A_1D(line) = ( w1*A(i+1,j+1,k+1) + &
		        	     w2*A(i  ,j+1,k+1) + w3*A(i+1,j  ,k+1) + &
		        	     w4*A(i  ,j  ,k+1) + w5*A(i+1,j+1,k  ) + &
		        	     w6*A(i  ,j+1,k  ) + w7*A(i+1,j  ,k  ) + &
		        	     w8*A(i  ,j  ,k  ))
        	  B_1D(line) = ( w1*B(i+1,j+1,k+1) + &
		  		     w2*B(i  ,j+1,k+1) + w3*B(i+1,j  ,k+1) + &
		  		     w4*B(i  ,j  ,k+1) + w5*B(i+1,j+1,k  ) + &
		  		     w6*B(i  ,j+1,k  ) + w7*B(i+1,j  ,k  ) + &
		  		     w8*B(i  ,j  ,k  ))
        	
        	endif
        	  endif
        	endif
       endif 	
 
 enddo
 
 do var =1, nb_var 
   print *,'nb var ', var
   if ((tab_var(var) == "Btot").or.(tab_var(var) == "Etot").or.(tab_var(var) == "Utot").or.(tab_var(var) == "Jtot")&
           .or.(tab_var(var) == "GPmagtot")) then
     Out_array(:,var) = sqrt(Ax_1D**2+Ay_1D**2+Az_1D**2)
   endif
   if ((tab_var(var) == "Bx").or.(tab_var(var) == "Ex").or.(tab_var(var) == "Ux").or.(tab_var(var) == "Jx") &
   .or.(tab_var(var) == "GPmagx")) then
     Out_array(:,var) = Ax_1D
   endif
   if ((tab_var(var) == "By").or.(tab_var(var) == "Ey").or.(tab_var(var) == "Uy").or.(tab_var(var)=="Jy") &
   .or.(tab_var(var)=="GPmagy")) then
     Out_array(:,var) = Ay_1D
   endif
   if ((tab_var(var) == "Bz").or.(tab_var(var) == "Ez").or.(tab_var(var) == "Uz").or.(tab_var(var) == "Jz")&
   .or.(tab_var(var) == "GPmagz")) then
     Out_array(:,var) = Az_1D
   endif    
   if (tab_var(var) == "Density") Out_array(:,var) = A_1D
   if (tab_var(var) == "Temperature") Out_array(:,var) = B_1D

 enddo
 
 
 
 ! Setting the position in simulation unit (c/wpi)
 ! and in the simulation frame
 !  X_MSO(:) = -X_MSO(:)/c_wpi + centr(1)
 !  Y_MSO(:) = -Y_MSO(:)/c_wpi + centr(2)
 !  Z_MSO(:) =  Z_MSO(:)/c_wpi + centr(3)
 print *, 'Max & Min X_MSO',maxval(X_MSO),minval(X_MSO)
 print *, 'Max & Min Y_MSO',maxval(Y_MSO),minval(Y_MSO)
 print *, 'Max & Min Z_MSO',maxval(Z_MSO),minval(Z_MSO)
   
   
! deallocate(Ax,Ay,Az)
! deallocate(Ax_1D,Ay_1D,Az_1D)
! if ((ind_mag ==0).and.(ind_vel==0).and.(ind_ele==0)) then
! deallocate(A,B,A_1D,B_1D)
!endif 
  if (allocated(Ax)) deallocate(Ax)
  if (allocated(Ay)) deallocate(Ay)
  if (allocated(Az)) deallocate(Az)
  if (allocated(Ax_1D)) deallocate(Ax_1D)
  if (allocated(Ay_1D)) deallocate(Ay_1D)
  if (allocated(Az_1D)) deallocate(Az_1D)
  if (allocated(A)) deallocate(A)
  if (allocated(B)) deallocate(B)
  if (allocated(A_1D)) deallocate(A_1D)
  if (allocated(B_1D)) deallocate(B_1D)  

 
__WRT_DEBUG_OUT("extract_field_value")  
 end subroutine extract_field_value
 
  !!##############################################################
   !! IMPEX_WS/read_dim_field_cdf
   !! This routines reads a simulation file
   !! and extract the dimension of a given variabe
   !! associated to a tag
   subroutine read_dim_field_cdf(cube_name,var_name,ncm_tot,radius,centr,gstep,planet,c_wpi)
        
   character(len=*),intent(in) :: cube_name
   character(len=*),intent(in) :: var_name
   integer,intent(inout)         :: ncm_tot(3)
   real(dp),intent(out)        :: radius,centr(3),gstep(3),c_wpi
   character(len=8),intent(out) :: planet
   
   !--Others
   integer  :: stId,ncid
   logical :: file_e
   character(len=50) :: write_name
   integer :: varfieldId, numDims,var_id
   integer,dimension(nf90_max_var_dims) :: DimfieldId
       
__WRT_DEBUG_IN("read_dim_field_cdf")        
     
   !--Inquire if the file exists
   inquire( file=trim(cube_name), exist=file_e )
   !--No file: return
   if(.not.(file_e)) then
     call wrtout(6,"File: "//trim(cube_name)//" does not exits","PERS")
      return
   endif
     
   !--Open NetCDF file
   stId = nf90_open(trim(cube_name), nf90_nowrite, ncid)
   call test_cdf(stId)
                
         
   !--Get number of point for fields in X,Y,Z
   StId = nf90_inq_varId(ncid,var_name,varfieldId)
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
  ! write(*,*) 'Planet :',planet
   

    stId = nf90_close(ncid);  call test_cdf(stId)
  !  write(*,*) 'close file ',cube_name   
   
__WRT_DEBUG_OUT("read_dim_field_cdf")    
   end subroutine read_dim_field_cdf
   
 !!=====================================================
 !! count_nb_var
 !! FUNCTION : count number of variable
 !! asked
 integer function count_nb_var(read_variable)
 implicit none
 character(len=*) :: read_variable
 character(len=200) :: test_var
 integer :: ind
 logical :: end_string
 
 __WRT_DEBUG_IN("count_nb_var") 
 count_nb_var = 1
 end_string = .false.
 test_var = read_variable
 do while (end_string.eqv..false.)
   ind = INDEX(test_var,",")
   if (ind /=0) then
     test_var = test_var(ind+1:len(test_var))
     count_nb_var = count_nb_var+1
   else
     end_string=.true.
   endif
 enddo
 

 __WRT_DEBUG_OUT("count_nb_var")
 end function count_nb_var
 
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


   do i=1,n_line
!   if (n_col >3) then
!     read(iunit,*) read_time,read_x,read_y,read_z
!     Time(i) = read_time
!   else
!       read(iunit,*) read_x,read_y,read_z
!   endif
     read(iunit,*) In_array(i,1:n_col)
     
     !In_array(i,1) =trim(a)
     !In_array(i,2) =trim(b)
     !In_array(i,3) =trim(c)
     read(In_array(i,startcol),*) read_x
     read(In_array(i,startcol+1),*) read_y
     read(In_array(i,startcol+2),*) read_z
     !print *,'x,y,z',read_x,read_y,read_z
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
 !!############################################################
    !! IMPEW_WS/read_field_cdf
    !! This routines reads the global output
    !! simulation files. Concern files are :
    !!		- Magw
    !!		- Elew
    !!		- Velw    
     subroutine read_field_cdf(run_name,varname1,varname2,varname3,Ax,Ay,Az)
         
     character(len=*),intent(in) :: run_name
     character(len=*),intent(in) :: varname1,varname2,varname3
     real(dp),dimension(:,:,:),intent(inout) :: Ax,Ay,Az
              
     !--Local variables
     character(len=64) :: varname    
     integer  :: stId,ncid
     logical :: file_e
     character(len=50) :: write_name
     character(len=64) :: filename
 
__WRT_DEBUG_IN("read_field_cdf")     
     
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
              
     
     call get_simple_variable_cdf(ncid,varname1 ,Ax(:,:,:) )
     call get_simple_variable_cdf(ncid,varname2 ,Ay(:,:,:) )
     call get_simple_variable_cdf(ncid,varname3 ,Az(:,:,:) )  
 
      
     
     stId = nf90_close(ncid);  call test_cdf(stId)
     write(*,*) 'close file ',run_name
__WRT_DEBUG_OUT("read_field_cdf")      
    end subroutine read_field_cdf
    
 !!############################################################
    !! IMPEW_WS/read_species_cdf
    !! This routines reads the global output
    !! simulation files. Concern files are :
    !!		- all ion species files
  
     subroutine read_species_cdf(run_name,varname1,varname2,varname3,varname4,varname5,A,Ax,Ay,Az,B)
         
     character(len=*),intent(in) :: run_name
     character(len=*),intent(in) :: varname1,varname2,varname3,varname4,varname5
     real(dp),dimension(:,:,:),intent(inout) :: A,Ax,Ay,Az,B
              
     !--Local variables
     character(len=64) :: varname    
     integer  :: stId,ncid
     logical :: file_e
     character(len=50) :: write_name
     character(len=64) :: filename
 
__WRT_DEBUG_IN("read_species_cdf")      
     
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
     
__WRT_DEBUG_OUT("read_species_cdf")      
    end subroutine read_species_cdf
    
    
    !!###################################################
    !! m_IMPEX_diag/save_value_orbit
    !! This routines dumps into a file
    !! the extracted value along a trajectory
    subroutine save_field_value(output_filename,X_MSO,Y_MSO,Z_MSO,Out_array,In_array,nb_var,tab_var,ncol,startcol)
    character(len=*),intent(in) ::output_filename
    real(dp),dimension(:),intent(in) :: X_MSO,Y_MSO,Z_MSO
    real(dp),dimension(:,:),intent(in) :: Out_array
    character(len=*),dimension(:,:),intent(in) :: In_array
    integer,intent(in) :: nb_var,ncol,startcol
    character(len=*),dimension(:),intent(in) :: tab_var
    character(len=32) :: val_ucd,unit
    integer :: io,n,l
    character(len=300) :: readfile,msg

    ! Local variables
    integer :: iunit,i,nb
__WRT_DEBUG_IN("save_field_value") 
    
    !-- Get length of the array
    nb = size(X_MSO,1)
    
    do l=1,ncol
      print *,trim(In_array(1,l))
   enddo
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
    do i=1,nb_var
     if ((tab_var(i) == "Bx").or.(tab_var(i) == "By").or.(tab_var(i) == "Bz").or.(tab_var(i) == "Btot")) then
       unit = "nT";	val_ucd = "phys.magField"
     endif
     if ((tab_var(i) == "Ex").or.(tab_var(i) == "Ey").or.(tab_var(i) == "Ez").or.(tab_var(i) == "Etot")) then
       unit = "mV.m-1";	val_ucd = "phys.elecField"
     endif    
     if ((tab_var(i) == "Ux").or.(tab_var(i) == "Uy").or.(tab_var(i) == "Uz").or.(tab_var(i) == "Utot")) then
       unit = "km.s-1";	val_ucd = "phys.veloc"
     endif
     if ((tab_var(i) == "Jx").or.(tab_var(i) == "Jy").or.(tab_var(i) == "Jz").or.(tab_var(i) == "Jtot")) then
       unit = "nA.m-2";	val_ucd = "phys.electField"
     endif
     if ((tab_var(i) == "GPmagx").or.(tab_var(i) == "GPmagy").or.(tab_var(i) == "GPmagz").or.(tab_var(i) == "GPmagtot")) then
       unit = "N.m-3";	val_ucd = "phys.electField"
     endif
     if (tab_var(i) =="Density") then
       unit = "cm-3"; val_ucd = "phys.density"
     endif
     if (tab_var(i) == "Temperature") then
       unit ="eV"; val_ucd="phys.temperature"
     endif  
     if (ncol+i <10) then
      write(iunit,'(3a,i1,5a)') '<FIELD name="',trim(tab_var(i)),'" ID="col',ncol+i,'" ucd="',trim(val_ucd),&
        '" utype="" datatype="float"  unit="',trim(unit),'" />'
     else
       write(iunit,'(3a,i2,5a)') '<FIELD name="',trim(tab_var(i)),'" ID="col',ncol+i,'" ucd="',trim(val_ucd),&
        '" utype="" datatype="float"  unit="',trim(unit),'" />'
     endif       
    enddo
    	
  write(iunit,'(a)') '<DATA>'
  write(iunit,'(a)') '<TABLEDATA>'    	
    !call header_orbit_VOTABLE(iunit,prefix,planetname,coord)	
    
     do i=1,nb
 !    write(iunit,'(5f10.3)') X_MSO(i),Y_MSO(i),Z_MSO(i),Out_array(i,1),Out_array(i,2)
         write(iunit,'(a)') '<TR>' 
    !call get_VOTIME(YYYY(i),MM(i),DD(i),hh(i),minut(i),ss(i),votime)
         msg =''
         do n=1,ncol     
           write(msg,'(a)') trim(msg)//'<TD> '//trim(In_array(i,n))//' </TD> '
 !          msg = trim(msg)
 !          print *,'msg=',msg
         enddo  
 !        print *,'OLD MSG',trim(msg)         
 !        print *,'msg=',trim(adjustl(In_array(i,n)))      
         do n=1,nb_var
           write(msg,'(a,f10.3,a)') trim(msg)//'<TD> ',Out_array(i,n),' </TD> '
         enddo
         !print *,'write in file:',trim(msg)
         write(iunit,'(a)') trim(msg)
    
!      if (nb_var  == 1) then
!            write(iunit,'(4(a,f13.3),a)') '<TD>',X_MSO(i),'</TD> <TD>',Y_MSO(i),'</TD> <TD>',Z_MSO(i), &
!            & '</TD> <TD>',Out_array(i,1),'</TD>'
!          endif
!          if (nb_var  == 2) then
!            write(iunit,'(5(a,f13.3),a)') '<TD>',X_MSO(i),'</TD> <TD>',Y_MSO(i),'</TD> <TD>',Z_MSO(i), &
!            & '</TD> <TD>',Out_array(i,1),'</TD> <TD>',Out_array(i,2),'</TD>'
!          endif
!          if (nb_var  == 3) then
!            write(iunit,'(6(a,f13.3),a)') '<TD>',X_MSO(i),'</TD> <TD>',Y_MSO(i),'</TD> <TD>',Z_MSO(i), &
!            & '</TD> <TD>',Out_array(i,1),'</TD> <TD>',Out_array(i,2),'</TD> <TD>',Out_array(i,3),'</TD>'
!          endif
!           if (nb_var  == 4) then
!            write(iunit,'(7(a,f13.3),a)') '<TD>',X_MSO(i),'</TD> <TD>',Y_MSO(i),'</TD> <TD>',Z_MSO(i), &
!            & '</TD> <TD>',Out_array(i,1),'</TD> <TD>',Out_array(i,2),'</TD> <TD>',Out_array(i,3), &
!            & '</TD> <TD>',Out_array(i,4),'</TD>'
!          endif   
!           if (nb_var  == 5) then
!            write(iunit,'(8(a,f13.3),a)') '<TD>',X_MSO(i),'</TD> <TD>',Y_MSO(i),'</TD> <TD>',Z_MSO(i), &
!            & '</TD> <TD>',Out_array(i,1),'</TD> <TD>',Out_array(i,2),'</TD> <TD>',Out_array(i,3), &
!            & '</TD> <TD>',Out_array(i,4),'</TD> <TD>',Out_array(i,5),'</TD>'
!          endif                     
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
__WRT_DEBUG_OUT("save_field_value") 
    end subroutine save_field_value
 

end program IMPEX_Interpole_Fields




