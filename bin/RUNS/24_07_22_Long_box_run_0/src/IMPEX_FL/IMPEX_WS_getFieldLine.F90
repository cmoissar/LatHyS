program IMPEX_getFieldLine

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

 character(len=*), parameter :: version = '0.1'
 character(len=200) :: arg,read_variable="",input_filename="",output_filename=""
 character(len=300) :: cube_name=""
  character(len=8) :: planetname=""
 character(len=200),dimension(:,:),allocatable :: In_array
 real(dp),dimension(:),allocatable :: X_MSO,Y_MSO,Z_MSO
 real(dp),dimension(:,:),allocatable :: Out_array
 integer :: ncommand,i,length_tab,ind,n_line,n_col=-1,startcol=-1,l
 integer :: direction=0
 integer :: lmax=1000 ! maximum of point for 1 field line
 real(dp) :: stepsize=-1.
 real(dp) :: t1,t2,clockangle=2000,unit_traj=1.
 __WRT_DEBUG_IN("IMPEX_getFieldLine") 
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
  case ('-p', '--direction')
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) direction  
  case ('-s', '--stepsize')
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) stepsize    
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
 if (((direction/=0).and.(direction/=1)).and.(direction/=-1)) then
   print *,"Warning argument -p is missing, assigned default value"
   direction = 0
 endif
 if (stepsize <0.) then
   print *,"Warning argument -s is missing, assigned default value"
   stepsize = 1.
 endif  
 

  call wrtout(6,ch10//&
       " ==================== Treating file run  "//&
       &trim(cube_name)//&
       " ====================",'PERS')

 
 
 
 ! Read input file size
 call read_file_size(input_filename,n_line)
 ! si le calcul des lignes de champ s'effectue beckward et forward alors le nombre de ligne de 
 ! champ est multiplié par 2 ( 1 ligne backward + 1 ligne forward)

 
 ! creation of array of position and time
 allocate(X_MSO(n_line), Y_MSO(n_line), Z_MSO(n_line))
 X_MSO(:) = 0.; Y_MSO(:) = 0.; Z_MSO(:) = 0.
 
 allocate(In_array(n_line,n_col)); In_array(:,:) = ""
 

 
 
  call read_input_file(input_filename,n_col,startcol,X_MSO,Y_MSO,Z_MSO,In_array,unit_traj)
 ! do l=1,n_col
 ! print *,trim(In_array(n_line,l))
 !  enddo
 ! extract values
 
  if (direction ==0) n_line=n_line*2
  allocate(Out_array(n_line*lmax,8)); Out_array(:,:) = 0.

 call calculate_field_line(cube_name,X_MSO,Y_MSO,Z_MSO,Out_array,length_tab,clockangle,direction,stepsize,unit_traj,planetname)
 print *,'Field line calculation ------done'
 print *,'Total point Field Line : ',length_tab
 print *,'Out_array ',Out_array(length_tab-1,1),Out_array(length_tab-1,2),Out_array(length_tab-1,3), &
 		Out_array(length_tab-1,4),Out_array(length_tab-1,5),Out_array(length_tab-1,6),Out_array(length_tab-1,7), &
 		Out_array(length_tab-1,8)

 
 call save_field_line(output_filename,X_MSO,Y_MSO,Z_MSO,Out_array,length_tab,cube_name,planetname)


 !deallocate(tab_var)
 deallocate(X_MSO,Y_MSO,Z_MSO,In_array,Out_array)

 !-- Finish time
 call CPU_TIME(t2)
 
 
 print '("CPU Time = ",f12.3," seconds.")',t2-t1
 __WRT_DEBUG_OUT("IMPEX_getFieldLine")
  
 
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
  print '(a)', "  -p, --direction     direction (both=0,forward=+1,backward=-1) [required]"
  print '(a)', "  -s, --stepsize      stepsize (in c/wpi) [optional]"

 end subroutine print_help
 
 !!===============================================================
 !! routine : extract_value
 !! FUNCTION
 !! extract field information from simulation files
 
 subroutine calculate_field_line(cube_name,X_MSO,Y_MSO,Z_MSO,Out_array,lcount,clockangle,direction,stepsize,unit_traj,planet)
 implicit none
 character(len=*),intent(in) :: cube_name
 real(dp),dimension(:),intent(inout) ::X_MSO,Y_MSO,Z_MSO
 real(dp),intent(in) :: clockangle,unit_traj
 real(dp),intent(inout) :: stepsize
 integer,intent(in) :: direction
 integer,intent(inout) :: lcount
 real(dp),dimension(:,:),intent(inout) :: Out_array
  character(len=8),intent(out) :: planet
 ! Local Variables
 character(len=32) :: var_name1="",var_name2="",var_name3="",var_name4="",var_name5=""
 integer :: ncm_tot(3),ind_mag,ind_vel,ind_ele,ind_spec,var,line,n_line,count_line_local
 real(dp)  :: centr(3),radius,gstep(3),c_wpi ,pos(3)
 real(dp) :: xa,ya,za,xf,yf,zf,w1,w2,w3,w4,w5,w6,w7,w8
 integer :: i,j,k,sgn_direction,sgn_coord
 logical :: continue_FL=.true.
 real(dp),dimension(:,:,:),allocatable :: Ax, Ay, Az
 real(dp) :: Ax_loc,Ay_loc,Az_loc,A_loc,dist
 real(dp),dimension(:),allocatable :: X_FL,Y_FL,Z_FL
 real(dp),dimension(:),allocatable :: X0_traj,Y0_traj,Z0_traj
 
__WRT_DEBUG_IN("calculate_field_line")  

 

 
 ind_mag = INDEX(cube_name,"Magw")
 if (ind_mag /= 0) then
    var_name1 = "Bx"; var_name2 = "By"; var_name3 = "Bz"
 endif
 
 ind_ele = INDEX(cube_name,"Elew")
  if (ind_ele /= 0) then
     var_name1 = "Ex"; var_name2 = "Ey"; var_name3 = "Ez"
  endif
  
  ind_vel = INDEX(cube_name,"Thew")
   if (ind_vel /= 0) then
      var_name1 = "Ux"; var_name2 = "Uy"; var_name3 = "Uz"
   endif
   
  if (((ind_vel ==0 ) .and. (ind_mag ==0)).and.(ind_ele == 0)) then
      var_name1 = "Ux"; var_name2 = "Uy"; var_name3 = "Uz"
   endif
  
! if ((ind_mag ==0).and.(ind_vel==0).and.(ind_ele==0)) then
!     var_name1 = "Density"; var_name2 = "Ux";	var_name3 = "Uy";	var_name4 = "Uz"
!     var_name5 = "Temperature"
! endif    
  
 print *,'Variable',var_name1 
 
 call read_dim_field_cdf(cube_name,var_name1,ncm_tot,radius,centr,gstep,planet,c_wpi)

 ! coordinate transformation for each object
  select case (trim(planet))
  case ("mars","mercury")
      sgn_coord = -1
    case ("ganymede")
      sgn_coord=+1
   case default
      sgn_coord=-1
   end select
  
 if (stepsize == 1.) stepsize=gstep(1)/4.
  
  
 allocate(Ax(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 allocate(Ay(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 allocate(Az(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 !allocate(Amod(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
 Ax = zero;	Ay = zero;	Az = zero;!	Amod = zero

! if ((ind_mag ==0).and.(ind_vel==0).and.(ind_ele==0)) then
!    allocate(A(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
!    allocate(B(ncm_tot(1),ncm_tot(2),ncm_tot(3)))
!    A = zero;	B=zero
! endif  
 
! if ((ind_mag /=0).or.(ind_vel/=0).or.(ind_ele/=0)) then
   call read_field_cdf(cube_name,var_name1,var_name2,var_name3, Ax,Ay,Az)
! else
!   call read_species_cdf(cube_name,var_name1,var_name2,var_name3,var_name4,var_name5, A,Ax,Ay,Az,B)
! endif
 
 !Amod = sqrt(Ax**2+Ay**2+Az**2)
 
 ! extracting value at the position
 ! 3D-linear interpolation
 n_line = int(size(Out_array,1)/lmax)
 print *,'n_line=',n_line

 
 
 ! rotating the trajectory if needed
 !===================================
 
 !print *,'Colck angle value',clockangle
 !if (clockangle /= 2000.) then
 !allocate(X0_traj(n_line),Y0_traj(n_line),Z0_traj(n_line));
 !X0_traj(:) = zero;	Y0_traj(:) = zero;	Z0_traj(:)=zero
 !X0_traj = X_MSO; Y0_traj = Y_MSO; Z0_traj = Z_MSO
 !X_MSO=X0_traj
 !Y_MSO=Y0_traj*cos((clockangle+90)/180.*pi)-Z0_traj*sin((clockangle+90)/180.*pi)
 !Z_MSO=Z0_traj*cos((clockangle+90)/180.*pi)+Y0_traj*sin((clockangle+90)/180.*pi) 
 !deallocate(X0_traj,Y0_traj,Z0_traj)
 !endif
 
   
 lcount = 0
 do line=1,n_line
   sgn_direction = direction
   if (direction ==0) then 
   sgn_direction=+1
   pos(1) = (sgn_coord*X_MSO(int((line+1)/2.))/c_wpi+centr(1))
   pos(2) = (sgn_coord*Y_MSO(int((line+1)/2.))/c_wpi+centr(2))
   pos(3) = ( Z_MSO(int((line+1)/2.))/c_wpi+centr(3))
   !print *,line,n_line
   else
   pos(1) = (sgn_coord*X_MSO(line)/c_wpi+centr(1))
   pos(2) = (sgn_coord*Y_MSO(line)/c_wpi+centr(2))
   pos(3) = ( Z_MSO(line)/c_wpi+centr(3))
   endif
   
   continue_FL=.true.
   count_line_local=1
   dist = sqrt(dot_product(pos-centr,pos-centr))
   if (dist < radius) continue_FL=.false.
   if ((int(pos(1)) < 0).or. (int(pos(1))>(ncm_tot(1)-2)*gstep(1))) continue_FL=.false.
   if ((int(pos(2)) < 0).or. (int(pos(2))>(ncm_tot(2)-2)*gstep(2))) continue_FL=.false.
   if ((int(pos(3)) < 0).or. (int(pos(3))>(ncm_tot(3)-2)*gstep(3))) continue_FL=.false.
   !do while ((continue_FL==.true.).and.(lcount<n_line*lmax))
   do while ((continue_FL==.true.).and.(count_line_local<lmax))
       count_line_local = count_line_local +1
      ! print *,count_line_local,lmax
       lcount = lcount+1
     ! determine the field at the position
       i = int(pos(1)/gstep(1))+1
       j = int(pos(2)/gstep(2))+1
       k = int(pos(3)/gstep(3))+1
       xa = pos(1)/gstep(1)-float(i)
       ya = pos(2)/gstep(2)-float(j)
       za = pos(3)/gstep(3)-float(k)
       
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
          	
       Ax_loc = sgn_coord*( w1*Ax(i+1,j+1,k+1) + &
        	     w2*Ax(i  ,j+1,k+1) + w3*Ax(i+1,j  ,k+1) + &
        	     w4*Ax(i  ,j  ,k+1) + w5*Ax(i+1,j+1,k  ) + &
        	     w6*Ax(i  ,j+1,k  ) + w7*Ax(i+1,j  ,k  ) + &
        	     w8*Ax(i  ,j  ,k  ))
 
       Ay_loc = sgn_coord*( w1*Ay(i+1,j+1,k+1) + &
        	     w2*Ay(i  ,j+1,k+1) + w3*Ay(i+1,j  ,k+1) + &
        	     w4*Ay(i  ,j  ,k+1) + w5*Ay(i+1,j+1,k  ) + &
        	     w6*Ay(i  ,j+1,k  ) + w7*Ay(i+1,j  ,k  ) + &
        	     w8*Ay(i  ,j  ,k  ))

       Az_loc = ( w1*Az(i+1,j+1,k+1) + &
        	     w2*Az(i  ,j+1,k+1) + w3*Az(i+1,j  ,k+1) + &
        	     w4*Az(i  ,j  ,k+1) + w5*Az(i+1,j+1,k  ) + &
        	     w6*Az(i  ,j+1,k  ) + w7*Az(i+1,j  ,k  ) + &
        	     w8*Az(i  ,j  ,k  ))  
       A_loc = sqrt(Ax_loc**2+Ay_loc**2+Az_loc**2)
       
       Out_array(lcount,1) = float(line)
       Out_array(lcount,2) = sgn_coord*(pos(1)-centr(1))*c_wpi
       Out_array(lcount,3) = sgn_coord*(pos(2)-centr(2))*c_wpi
       Out_array(lcount,4) =  (pos(3)-centr(3))*c_wpi
       Out_array(lcount,5) =  sgn_coord*Ax_loc
       Out_array(lcount,6) =  sgn_coord*Ay_loc
       Out_array(lcount,7) =  Az_loc
       Out_array(lcount,8) =  A_loc
       !print *, 'Position :',Out_array(lcount,1),Out_array(lcount,2),Out_array(lcount,3),Out_array(lcount,4)
     print *,'i,j,k, Values ',i,j,k,Ax(i,j,k),Ay(i,j,k),Az(i,j,k)

       pos(1) = pos(1) + sgn_direction*stepsize*Ax_loc/A_loc
       pos(2) = pos(2) + sgn_direction*stepsize*Ay_loc/A_loc
       pos(3) = pos(3) + sgn_direction*stepsize*Az_loc/A_loc
       ! evaluation of loop criteria
       dist = sqrt(dot_product(pos-centr,pos-centr))
       if (dist < radius) continue_FL=.false.
       if ((int(pos(1)) < 0).or. (int(pos(1))>(ncm_tot(1)-2)*gstep(1))) continue_FL=.false.
       if ((int(pos(2)) < 0).or. (int(pos(2))>(ncm_tot(2)-2)*gstep(2))) continue_FL=.false.
       if ((int(pos(3)) < 0).or. (int(pos(3))>(ncm_tot(3)-2)*gstep(3))) continue_FL=.false.
   enddo ! loop for one field line
   
   
  if (direction ==0) then
   sgn_direction=-1
   pos(1) = (sgn_coord*X_MSO(int((line+1)/2.))/c_wpi+centr(1))
   pos(2) = (sgn_coord*Y_MSO(int((line+1)/2.))/c_wpi+centr(2))
   pos(3) = ( Z_MSO(int((line+1)/2.))/c_wpi+centr(3))
   
   continue_FL=.true.
   count_line_local=1
   dist = sqrt(dot_product(pos-centr,pos-centr))
   if (dist < radius) continue_FL=.false.
   if ((int(pos(1)) < 0).or. (int(pos(1))>(ncm_tot(1)-2)*gstep(1))) continue_FL=.false.
   if ((int(pos(2)) < 0).or. (int(pos(2))>(ncm_tot(2)-2)*gstep(2))) continue_FL=.false.
   if ((int(pos(3)) < 0).or. (int(pos(3))>(ncm_tot(3)-2)*gstep(3))) continue_FL=.false.
   !do while ((continue_FL==.true.).and.(lcount<n_line*lmax))
   do while ((continue_FL==.true.).and.(count_line_local<lmax))
       count_line_local = count_line_local +1
       lcount = lcount+1
     ! determine the field at the position
       i = int(pos(1)/gstep(1))+1
       j = int(pos(2)/gstep(2))+1
       k = int(pos(3)/gstep(3))+1
       xa = pos(1)/gstep(1)-float(i)
       ya = pos(2)/gstep(2)-float(j)
       za = pos(3)/gstep(3)-float(k)
       
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
          	
       Ax_loc = sgn_coord*( w1*Ax(i+1,j+1,k+1) + &
        	     w2*Ax(i  ,j+1,k+1) + w3*Ax(i+1,j  ,k+1) + &
        	     w4*Ax(i  ,j  ,k+1) + w5*Ax(i+1,j+1,k  ) + &
        	     w6*Ax(i  ,j+1,k  ) + w7*Ax(i+1,j  ,k  ) + &
        	     w8*Ax(i  ,j  ,k  ))
 
       Ay_loc = sgn_coord*( w1*Ay(i+1,j+1,k+1) + &
        	     w2*Ay(i  ,j+1,k+1) + w3*Ay(i+1,j  ,k+1) + &
        	     w4*Ay(i  ,j  ,k+1) + w5*Ay(i+1,j+1,k  ) + &
        	     w6*Ay(i  ,j+1,k  ) + w7*Ay(i+1,j  ,k  ) + &
        	     w8*Ay(i  ,j  ,k  ))

       Az_loc = ( w1*Az(i+1,j+1,k+1) + &
        	     w2*Az(i  ,j+1,k+1) + w3*Az(i+1,j  ,k+1) + &
        	     w4*Az(i  ,j  ,k+1) + w5*Az(i+1,j+1,k  ) + &
        	     w6*Az(i  ,j+1,k  ) + w7*Az(i+1,j  ,k  ) + &
        	     w8*Az(i  ,j  ,k  ))  
       A_loc = sqrt(Ax_loc**2+Ay_loc**2+Az_loc**2)
       
       Out_array(lcount,1) = float(line)
       Out_array(lcount,2) = sgn_coord*(pos(1)-centr(1))*c_wpi
       Out_array(lcount,3) = sgn_coord*(pos(2)-centr(2))*c_wpi
       Out_array(lcount,4) =  (pos(3)-centr(3))*c_wpi
       Out_array(lcount,5) =  sgn_coord*Ax_loc
       Out_array(lcount,6) =  sgn_coord*Ay_loc
       Out_array(lcount,7) =  Az_loc
       Out_array(lcount,8) =  A_loc
       !print *, 'Position :',Out_array(lcount,1),Out_array(lcount,2),Out_array(lcount,3),Out_array(lcount,4)

       pos(1) = pos(1) + sgn_direction*stepsize*Ax_loc/A_loc
       pos(2) = pos(2) + sgn_direction*stepsize*Ay_loc/A_loc
       pos(3) = pos(3) + sgn_direction*stepsize*Az_loc/A_loc
       ! evaluation of loop criteria
       dist = sqrt(dot_product(pos-centr,pos-centr))
       if (dist < radius) continue_FL=.false.
       if ((int(pos(1)) < 0).or. (int(pos(1))>(ncm_tot(1)-2)*gstep(1))) continue_FL=.false.
       if ((int(pos(2)) < 0).or. (int(pos(2))>(ncm_tot(2)-2)*gstep(2))) continue_FL=.false.
       if ((int(pos(3)) < 0).or. (int(pos(3))>(ncm_tot(3)-2)*gstep(3))) continue_FL=.false.
   enddo ! loop for one field line   
   endif ! case both direction integration (backward and forward)
   
 enddo ! loop over all footprint
  
 
 
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
 
__WRT_DEBUG_OUT("calculate_field_line")  
 end subroutine calculate_field_line
 
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
    
 
    !!######################################################
    !! m_IMPEX_diag/header_field_VOTABLE
    !! This routine write the header of the 
    !! VOTABLE file
    subroutine header_fieldline_VOTABLE(iunit,cube_name,planetname)
    character(len=*),intent(in) :: cube_name,planetname
    integer,intent(in) :: iunit
    character(len=32) :: varname1="",varname2="",varname3="",varname4="" 
    character(len=32) :: val_ucd,unit
    integer :: ind_mag,ind_ele,ind_vel
__WRT_DEBUG_IN("header_fieldline_VOTABLE")         
 
    
    write(iunit,'(a)') '<?xml version="1.0"?>'
    write(iunit,'(a)') '<VOTABLE version="1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
    write(iunit,'(a)') 'xmlns="http://www.ivoa.net/xml/VOTable/v1.2" xmlns:stc="http://www.ivoa.net/xml/STC/v1.30" >'
    write(iunit,'(a)') '<RESOURCE name="IMPEx Field/Flow line">'
    write(iunit,'(a)') '<TABLE name="Lines">'
    write(iunit,'(a)') '<DESCRIPTION> File containing the position of each field/flow lines.'
    write(iunit,'(a)') '</DESCRIPTION>'
  

 print *,planetname 
  select case (trim(planetname))
  case("mars","mercury")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="MSO">'
   if (trim(planetname) == "mercury") then
     write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="2440." unit="km" ucd="phys.size"/>'
   else
     write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="3393." unit="km" ucd="phys.size"/>'
   endif  
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col3"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col4"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'
  case("ganymede")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="GPHIO">'
   write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="2634." unit="km" ucd="phys.size"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col3"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col4"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="GPHIO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'  
  case("titan")
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="TIIS">'
   write(iunit,'(a)') '<PARAM name="Radius" datatype="float" arraysize="*" value="2575." unit="km" ucd="phys.size"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col3"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col4"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="TIIS">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>' 
   case default
   write(iunit,'(a)')  '<GROUP ID="PosFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col2"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col3"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col4"/>'
   write(iunit,'(a)') '</GROUP>'
   write(iunit,'(a)') '<GROUP ID="FieldFrame" ref="MSO">'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C1" ref="col5"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C2" ref="col6"/>'
   write(iunit,'(a)') '<FIELDref utype="stc:AstroCoords.Position3D.Value3.C3" ref="col7"/>'
   write(iunit,'(a)') '</GROUP>'   
   
  end select 
  

  write(iunit,'(a)') '<FIELD name="Line" ID="col1" ucd="meta.number" utype="" datatype="int" width="3" />'
   ind_mag = INDEX(cube_name,"Magw")
   if (ind_mag /= 0) then
      unit = "nT"; val_ucd = "phys.magField"
      varname1 = "Bx"; varname2 = "By"; varname3 = "Bz"; varname4 = "Btot"
   endif
 
 ind_ele = INDEX(cube_name,"Elew")
  if (ind_ele /= 0) then
      unit = "mV.m-1";	val_ucd = "phys.elecField"
      varname1 = "Ex"; varname2 = "Ey"; varname3 = "Ez"; varname4 = "Etot"
  endif
  
  ind_vel = INDEX(cube_name,"Thew")
   if (ind_vel /= 0) then
      unit = "km.s-1";	val_ucd = "phys.veloc"
      varname1 = "Ux"; varname2 = "Uy"; varname3 = "Uz"; varname4 = "Utot"
   endif  
   
    if (((ind_vel ==0 ) .and. (ind_mag ==0)).and.(ind_ele == 0)) then
    unit = "km.s-1";	val_ucd = "phys.veloc"
            varname1 = "Ux"; varname2 = "Uy"; varname3 = "Uz"; varname4 = "Utot"
   endif
   
  write(iunit,'(2a)') '<FIELD name="X" ID="col2" ucd="pos.cartesian.x"',&
        ' utype="stc:AstroCoords.Position3D.Value3.C1" datatype="float"  unit="km" />'
  write(iunit,'(2a)') '<FIELD name="Y" ID="col3" ucd="pos.cartesian.y"',&
        ' utype="stc:AstroCoords.Position3D.Value3.C2" datatype="float"  unit="km" />'  
  write(iunit,'(2a)') '<FIELD name="Z" ID="col3" ucd="pos.cartesian.z"',&
        ' utype="stc:AstroCoords.Position3D.Value3.C3" datatype="float"  unit="km" />'    
  write(iunit,'(7a)') '<FIELD name="',trim(varname1),'" ID="col5" ucd="',trim(val_ucd),&
        '" utype="" datatype="float"  unit="',trim(unit),'" />'        
  write(iunit,'(7a)') '<FIELD name="',trim(varname2),'" ID="col6" ucd="',trim(val_ucd),&
        '" utype="" datatype="float"  unit="',trim(unit),'" />'
  write(iunit,'(7a)') '<FIELD name="',trim(varname3),'" ID="col7" ucd="',trim(val_ucd),&
        '" utype="" datatype="float"  unit="',trim(unit),'" />'
  write(iunit,'(7a)') '<FIELD name="',trim(varname4),'" ID="col8" ucd="',trim(val_ucd),&
        '" utype="" datatype="float"  unit="',trim(unit),'" />'          
  
  write(iunit,'(a)') '<DATA>'
  write(iunit,'(a)') '<TABLEDATA>'  
__WRT_DEBUG_OUT("header_fieldline_VOTABLE")           
end subroutine header_fieldline_VOTABLE

    
    
    !!###################################################
    !! m_IMPEX_diag/save_value_orbit
    !! This routines dumps into a file
    !! the extracted value along a trajectory
    subroutine save_field_line(output_filename,X_MSO,Y_MSO,Z_MSO,Out_array,lcount,cube_name,planetname)
    character(len=*),intent(in) ::output_filename
    character(len=*),intent(in) :: cube_name,planetname    
    real(dp),dimension(:),intent(in) :: X_MSO,Y_MSO,Z_MSO
    real(dp),dimension(:,:),intent(in) :: Out_array
    integer,intent(in) :: lcount
    character(len=32) :: val_ucd,unit
    integer :: io,n,l
    character(len=300) :: readfile,msg

    ! Local variables
    integer :: iunit,i,nb
__WRT_DEBUG_IN("save_field_value") 
    
    !-- Get length of the array
  !  nb = size(X_MSO,1)
    
!    do l=1,ncol
!      print *,trim(In_array(1,l))
!   enddo
  print *, 'writing header ',trim(output_filename)
!    ! in VOTale format
    iunit = 1
   open(UNIT = iunit, FILE = output_filename, STATUS = 'UNKNOWN', ACTION = 'WRITE', &
        FORM = 'FORMATTED',IOSTAT=io)
   !print *,io   	
! Header
    call header_fieldline_VOTABLE(iunit,cube_name,planetname)    
   
print *,'header called successfully'

 !   
     do i=1,lcount
         write(iunit,'(a)') '<TR>' 
         msg =''   
         write(msg,'(a,i3,a)') trim(msg)//'<TD> ',int(Out_array(i,1)),' </TD> '
         do n=2,8
           write(msg,'(a,f10.3,a)') trim(msg)//'<TD> ',Out_array(i,n),' </TD> '
         enddo
         write(iunit,'(a)') trim(msg)
                        
          write(iunit,'(a)') '</TR>'	
            
    enddo
 ! ! finalize the file
  write(iunit,'(a)') '</TABLEDATA>'
  write(iunit,'(a)') '</DATA>'
  write(iunit,'(a)') '</TABLE>'
  write(iunit,'(a)') '</RESOURCE>'
  write(iunit,'(a)') '</VOTABLE>'      
  close(iunit)
    

__WRT_DEBUG_OUT("save_field_line") 
    end subroutine save_field_line
 

end program IMPEX_getFieldLine




