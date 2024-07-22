!!=============================================================
!!=============================================================
!!module: m_cmdline
!! NAME
!!  m_cmdline (MMancini)
!!
!! FUNCTION
!!  Contains command lines to set some variables help ..  
!!
module m_cmdline
 
 use defs_basis
 use m_logo,only        : print_date

 implicit none
 private

  public  ::        &
       cmd_line

  private ::        &
       print_help

contains
 !!#####################################################################

 !!=============================================================
 !!routine: m_cmdline/cmd_line
 !!
 !! FUNCTION
 !!  Command line
 !!   
 subroutine cmd_line(re_st)
  use defs_parametre

  integer,intent(out) :: re_st

  integer :: i,ncommand,stop_opt
  real(dp) :: dx
  character(len=*), parameter :: version = "1.0"
  character(len=40) :: arg
  character(len=500) ::cmd

  stop_opt = 0

  call get_command(cmd)

  ncommand = command_argument_count()
  
  i = 1
  do while(i<=ncommand)
   call get_command_argument(i, arg)
   
   select case(arg)
   case("-pn", "--planetname")
    i = i+1
    call get_command_argument(i, arg)
    planetname = trim(arg)
   case("-s", "--qp_out")
    i = i+1
    call get_command_argument(i, arg)
    qp_out_name = trim(arg)
   case ('-d', '--date')
       i = i+1
       call get_command_argument(i, arg)
       print *,"DATE=",trim(arg)
   fildat = trim(arg)    
   case("-r", "--restart")
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) restart 
    if(restart<0 .or. restart>4) call print_help()
   case("-rf", "--restart_file")
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) re_st
    if(re_st/=0 .and. re_st/=1) call print_help()
   case("-nhm", "--maxiter")
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) nhm 
   case("-ncxyz", "--cell_size")
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) ncx0 
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) ncy0 
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) ncz0 
   case("-ncx")
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) ncx0 
   case("-ncy")
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) ncy0 
   case("-ncz")
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) ncz0 
   case("-dt","--time_step")
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) dt 
   case("-dx","--gstep")
    i = i+1
    call get_command_argument(i, arg)
    read(arg,*) dx
    gstep = dx
   case("-ug","--uniform_grid")
     i = i+1
     call get_command_argument(i, arg)
     ROI%excess = 0
     ROI%excess_part = 0
   case("-le","--load_exosphere")
     i = i+1
    call get_command_argument(i, arg)
    exospherename = trim(arg)
   case("-la","--load_atmosphere")
     i = i+1
    call get_command_argument(i, arg) 
    atmospherename = trim(arg)
   case("-li","--load_ionosphere")
     i = i+1
    call get_command_argument(i, arg) 
    ionospherename = trim(arg)    
   case("-v", "--version")
    print "(2a)", "quiet_plasma version ", version
    stop_opt = 1
   case("-h", "--help")
    call print_help();
    stop_opt = 1
   case("-nm","--no_mag")
     i = i+1
     call get_command_argument(i, arg)
     no_env_mag=1  !if 1 then no planetary magnetic field
   case("-ssl","--planet_ssl")
     i = i+1
     call get_command_argument(i, arg)
     read(arg,*) planet_ssl
   case("-sslat","--planet_sslat")
     i = i+1
     call get_command_argument(i, arg)
     read(arg,*) planet_sslat
   case default
    print "(a,a,/)", "Unrecognized command-line option: ", arg
    call print_help(); 
   end select

   i = i+1
  end do

  if(stop_opt==1) stop

 end subroutine cmd_line


 !!=============================================================
 !!routine: m_cmdline/print_help
 !!
 !! FUNCTION
 !!  Print help menu
 !!   
 subroutine print_help()

  use defs_mpitype,only    : mpiinfo
  use m_logo

  call logo(6,date_ok=0)
  if(mpiinfo%me==0)then
   print '(a)', ""
   print '(a)', " quiet_plasma options:"
   print '(a)', ""
   print '(a)', "  -v, --version         print version information and exit"
   print '(a)', "  -h, --help            print usage information and exit"
   print '(a)', "  -t, --time            print time"
   print '(a)', "  -s, --qp_out          set the ouput file name"
   print '(a)', "  -pn, --planetname     set the planet environment,"
   print '(a)', "                        followed by a string:"
   print '(a)', "                               mars"
   print '(a)', "                               mercure"
   print '(a)', "                               moon"
   print '(a)', "                               ganymede"
   print '(a)', "                               mars3try"
   print '(a)', "                               titan"
   print '(a)', "                               shock"
   print '(a)', "                               venus"
   print '(a)', "  -r, --restart         set the restart,"
   print '(a)', "                        followed by:"
   print '(a)', "                            0, restart is off"
   print '(a)', "                            1, simple restart"
   print '(a)', "                            2, restart adding cells in Y"
   print '(a)', "                            3, restart adding cells in X"
   print '(a)', "                        Note: when restarting be careful to impose the same"
   print '(a)', "                              number of cells (ncx,ncy,ncz) and  planet name"
   print '(a)', "                              used in the previuos run."
   print '(a)', "                              dt,dx are selected automatically."
   print '(a)', "  -rf, --restart_file   set the restart file,"
   print '(a)', "                        followed by 0 or 1" 
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
   print '(a)', "  -dx, --gstep          set the spatial step (gstep)"
   print '(a)', "                        followed by 1 real (all the "
   print '(a)', "                         composent of gstep will have the same value)."
   print '(a)', "  -ug, --uniform_grid   set a uniform domain decomposition"
   print '(a)', "  -le, --load_exosphere load exospheric neutral density from a file"
   print '(a)', "                        the file name corresponds to input name"
   print '(a)', "  -la, --load_atmosphere load atmospheric neutral density from a file"
   print '(a)', "                        the file name corresponds to input name"   
   print '(a)', "  -li, --load_ionosphere load ionospheric density from a file"
   print '(a)', "                        the file name corresponds to input name"    
   print '(a)', "  -nm, --no_mag        set to 1 for unmagnetized Mars object"
   print '(a)', "  -ssl, --planet_ss    fix the solar longitude of the main crustal source"
   print '(a)', "  -sslat, --planet_sslat       fix the solar latitude of the main crustal source"
  end if
  stop
 end subroutine print_help

end module m_cmdline
