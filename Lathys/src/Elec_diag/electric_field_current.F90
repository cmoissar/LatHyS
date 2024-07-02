program electric_field_current

 use defs_basis
 use defs_arr3Dtype
 use defs_variable
 use m_writeout
 use m_logo,only        : print_date
 use defs_tregister
 use defs_parametre,only : fildat
 use netcdf
 use defs_basic_cdf
 use diag_wrt_common_cdf
 use m_current_computation

 implicit none

  !********************************************************************
  ! Auteur				:	 RModolo
  ! Date					:	 26/02/19
  ! Instituion				:	LATMOS/UVSQ/IPSL
  !						Quartier des Garennes
  ! 						11 bd d'alembert
  !						78280 guyancourt
  ! Dernière modification		:	26/02/19
  ! Résumé
  ! This program reads the reconstructed fiel (Magw, Dnw, ...)
  ! and create Jcurrent fiel as well an electric field file in which the contribution
  ! from each term is indicated
  !********************************************************************

 integer :: ireg,ncommand,i,factor
 real(dp) :: dx
 character(len=20) :: run_name,time
 character(len=5)  :: time_read=''
 character(len=32) :: arg
 integer :: Jcurr
 character(len=*), parameter :: version = '1.0'
 real(dp) :: t1,t2

 !-- Start time
 call CPU_TIME(t1)

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
  case ("-Jcurr","--Jcurrent")
   i = i+1
   call get_command_argument(i, arg)
   read(arg,*) Jcurr
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
!! call set_tregister(r_tm,0)
!!
!!  !--Factor to put any file name >1
!!  factor = 1
!!  do while(int(factor*r_tm%time(1))<1)
!!   factor = factor *10
!!  end do

 !! if (time_read =='') then
 !! !--Time of Diagnostic
 !!   write(time,'(a,i5.5)')"t",int(factor*r_tm%time(r_tm%nreg))
 !! else
      write(time,'(a,a)')"t",time_read
 !! endif

  !--Create the base for the Names of file
  run_name = trim(fildat)//"_"//time

  call wrtout(6,ch10//&
       " ==================== Treating file run  "//&
       &trim(run_name)//&
       " ====================",'PERS')

 ! if (Jcurr == 1) then
    print *, 'Ok pour diag Jcurr'
    call compute_current(run_name)
 ! endif

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
  print '(a)', ''
  print '(a)', 'diag options:'
  print '(a)', ''
  print '(a)', "  -v, --version       print version information and exit"
  print '(a)', "  -h, --help          print usage information and exit"
  print '(a)', "  -d, --date          use the date in the format dd_mm_yy"
  print '(a)', "  -t, --time          set diagnostic time in the format XXXXX"
  print '(a)', "  -Jcurr,--Jcurrent   set diagnostic for Electric current"
 end subroutine print_help


end program electric_field_current