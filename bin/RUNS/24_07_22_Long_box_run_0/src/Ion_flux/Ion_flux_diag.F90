program Ion_flux_diag

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
 use m_flux_computation
#endif

#include "q-p_common.h"

 implicit none
!********************************************************************
 ! Auteur				:	 RModolo 
 ! Date					:	 18/01/14
 ! Instituion				:	LATMOS/UVSQ/IPSL
 !						Quartier des Garennes
 ! 						11 bd d'alembert
 !						78280 guyancourt
 ! Dernière modification		:	17/01/14	
 ! Résumé	
 ! Ce programme lit des fichiers p3 (particules), et construit des
 ! tableaux de flux (n*|v|) pour chaque espece planétaires
 ! Pour les O+ on sépare les cartes de flux en deux parties 
 ! E<30eV et E> 30eV
 ! On calcule ensuite les flux d'échappement (total, par face, à differentes distances
 ! dans le sillage)
 !********************************************************************
 integer :: ireg,ncommand,i,factor
 character(len=20) :: run_name,time
 character(len=32) :: arg
 character(len=5)  :: time_read=''
 real(dp) :: Esep=30.
 character(len=*), parameter :: version = '1.0'

 !--Insert the date for which the diagnostic has to be reassmbled
 fildat = print_date()

 ncommand = command_argument_count()
 
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
 case ("-E","--Energy")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) Esep        
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
       
   call extract_fluxes(run_name,Esep) 
   
 
 

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
  print '(a)', "  -E, --Energy        set diagnostic Energy separation [in eV]"
 end subroutine print_help



end program Ion_flux_diag