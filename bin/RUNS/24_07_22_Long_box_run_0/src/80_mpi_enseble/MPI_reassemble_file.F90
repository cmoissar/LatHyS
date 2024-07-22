program mpi_reassemble_file

 use defs_basis
 use defs_arr3Dtype
 use defs_variable
 use m_writeout
 use defs_parametre,only : fildat,dt,ncx0,ncy0,ncz0,nhm,gstep
 use m_logo,only        : print_date
 use defs_tregister
#ifdef HAVE_NETCDF 
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
 use m_split_cdf
 use m_merge_global_cdf
#else
 use m_reassemble_unform
#endif


 implicit none

 !********************************************************************
 ! Auteur				:	 RModolo MMancini
 ! Date					:	 23/06/03
 ! Instituion				:	CETP/CNRS/IPSL
 !					10-12 avenue de l'Europe
 !					78140 VELIZY
 ! Dernière modification		:	19/07/10	
 ! Résumé	
 ! Ce programme lit des fichiers fortran non formatté à accès séquentiel.
 ! Chaque  fichier correspond à l'enregistrement des données d'un processus
 ! le programme traite ensuite les informations liées à chaque procesus 
 ! pour regrouper les informations 
 ! et sortir les tableaux de résultats telle que le fait un mono-processeur.
 ! On écrit ensuite le résultat dans un fichier séquentiel
 !********************************************************************

#ifdef HAVE_NETCDF
 integer :: cdf2bin=0,bin2cdf=0
#endif
 integer :: ireg,ncommand,i,factor
 real(dp) :: dx
 character(len=20) :: run_name,time
 character(len=32) :: arg
 character(len=5)  :: time_read='' 
 character(len=*), parameter :: version = '1.0'

 !--Insert the date for which the diagnostic has to be reassmbled
 fildat = print_date()

 ncommand = command_argument_count()

 i = 1
 do while(i<=ncommand)
  call get_command_argument(i, arg)
  select case (arg)
#ifdef HAVE_NETCDF
  case ('-cdf2bin')
   cdf2bin = 1
   stop "option not yet implemented"
  case ('-bin2cdf')
   bin2cdf = 1
   stop "option not yet implemented"
#endif
  case ('-d', '--date')
   i = i+1
   call get_command_argument(i, arg)
   print *,"DATE=",trim(arg)
   fildat = trim(arg)
  case ("-nhm", "--maxiter")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) nhm 
  case ("-ncxyz", "--cell_size")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) ncx0 
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) ncy0 
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) ncz0 
  case ("-ncx")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) ncx0 
  case ("-ncy")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) ncy0 
  case ("-ncz")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) ncz0 
  case ("-dt","--time_step")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) dt 
  case ("-dx","--gstep")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) dx
   gstep = dx
  case ("-t","--time")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) time_read    
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
 !call set_tregister(r_tm,0)

  !--Factor to put any file name >1
 ! factor = 1
 ! do while(int(factor*r_tm%time(1))<1) 
 !  factor = factor *10
 ! end do

  !print *, 'time_read ', time_read 
   if (time_read =='') then  
    !--Time of Diagnostic 
      call set_tregister(r_tm,0)
      write(time,'(a,i5.5)')"t",int(factor*r_tm%time(r_tm%nreg))
    else   
       if (time_read/='all')  write(time,'(a,a)')"t",time_read
       print *,'time ',time
  endif   

  !--Create the base for the Names of file
  run_name = trim(fildat)//"_"//time

  call wrtout(6,ch10//&
       " ==================== Joining file run the "//&
       &trim(run_name)//&
       " ====================",'PERS')


#ifdef HAVE_NETCDF 
  if(cdf2bin==1) then
   !   call reassemble_cdf2bin(run_name)
  else if(bin2cdf==1) then
   !   call reassemble_bin2cdf(run_name)
  else
  ! call reassemble_cdf(run_name)
  call split_cdf(run_name)
  call merge_cdf(run_name)
  endif
#else 
  call reassemble_unform(run_name)
#endif


   if (time_read =='all') then  
    !--Time of Diagnostic 
    do i=1,r_tm%nreg
       write(time,'(a,i5.5)')"t",int(factor*r_tm%time(i))   
    !--Create the base for the Names of file
     run_name = trim(fildat)//"_"//time

     call wrtout(6,ch10//&
       " ==================== Joining file run the "//&
       &trim(run_name)//&
       " ====================",'PERS')


#ifdef HAVE_NETCDF 
       if(cdf2bin==1) then
       !   call reassemble_cdf2bin(run_name)
       else if(bin2cdf==1) then
       !   call reassemble_bin2cdf(run_name)
       else
       ! call reassemble_cdf(run_name)
       call split_cdf(run_name)
       call merge_cdf(run_name)
       endif
#else 
      call reassemble_unform(run_name)
#endif
     enddo !-- iteration on diangostic files
     endif !-- il time_read='all'  


 !--Clean variables for time diagnostics
 call clean_tregister(r_tm)

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
  print '(a)', 'Without further options, diag rassebles files with the current data.'
  print '(a)', ''
  print '(a)', 'diag options:'
  print '(a)', ''
  print '(a)', '  -v, --version     print version information and exit'
  print '(a)', '  -h, --help        print usage information and exit'
  print '(a)', '  -d, --date        use the date in the format dd_mm_yy'
  print '(a)', '  -cdf2bin          convert from cdf to binary'
  print '(a)', '                         (only if netcdf is available)'
  print '(a)', '  -bin2cdf          convert from binary to cdf'
  print '(a)', '                         (only if netcdf is available)'
  print '(a)', "  -nhm, --maxiter       set the max number of iteration,"
  print '(a)', "                        followed by 1 integer."
  print '(a)', "  -ncxyz, --cell_size   set the cell number ncx,ncy,ncz at once"
  print '(a)', "                        followed by 3 integers."
  print '(a)', "  -ncx,                 set the cell number in X (ncx)"
  print '(a)', "                        followed by 1 integer."
  print '(a)', "  -ncy,                 set the cell number in Y (ncy)"
  print '(a)', "                        followed by 1 integer."
  print '(a)', "  -ncz,                 set the cell number in Z (ncz)"
  print '(a)', "                        followed by 1 integer."
  print '(a)', "  -dt, --time_step      set time step (dt)"
  print '(a)', "                        followed by 1 real."
  print '(a)', "  -t, --time            set diagnostic time"
  print '(a)', "                        in format XXXXX  with X an integer."
  print '(a)', "  -dx, --gstep          set the spatial step (gstep)"
  print '(a)', "                        followed by 1 real (all the "
  print '(a)', "                         composent of gstep have the same value)."
 end subroutine print_help



end program mpi_reassemble_file




