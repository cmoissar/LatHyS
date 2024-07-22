!!=============================================================
!!=============================================================
!!module: hybrid_model/m_writeout
!! NAME
!!  m_writeout (MMancini)
!!
!! FUNCTION
!!  Some function to print 
module m_writeout

 use defs_basis
 use defs_mpitype,only    : nproc,rang

 implicit none
 private

 integer,public,save  :: wrtdisk = 0   !--0=wrt on disk; 1=no
 integer,public,save  :: wrtscreen = 0 !--0=wrt on screen; 1=no
 integer,private,save :: dbg_adv = 1   !--Step shown for any function

 private::          &
      wrtout_myproc    !--To write only on a proc

 public::           &
      wrtout,       &  !--To write in general
      wrt_double,   &  !--To write on screen and disk
      wrt_debug        !--To write on the debug in-out from routines

contains
 !!############################################################


 !=============================================================
 !!subroutine: hybrid_model/wrtout_myproc
 !! FUNCTION (From ABINIT)
 !!  Do the output for one proc. For parallel or sequential output use wrtout()
 !!  instead. Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
 !!  
 !!
 !! INPUTS
 !!  unit=unit number for writing
 !!  message=(character(len=*)) message to be written
 !!
 !! OUTPUT
 !!  (only writing)
 subroutine wrtout_myproc(unit,message)

  integer,intent(in) :: unit
  character(len=*),intent(in) :: message
  integer,save :: iexit=0,ncomment=0,nwarning=0
  integer :: lenmessage,rtnpos
  character(len=len(message)) :: messtmp
  !******************************************************************
  if(message/=' ') then
   messtmp=message
   lenmessage=len(message)
   !  Here, split the message, according to the char(10)
   !  characters (carriage return). This technique is portable accross different OS.
   rtnpos=index(messtmp,ch10)
   !  do while(rtnpos/=0)
   do while(rtnpos > 1)
    write(unit, '(a)' ) trim(messtmp(1:rtnpos-1))
    messtmp=messtmp(rtnpos+1:lenmessage)
    lenmessage=lenmessage-rtnpos
    rtnpos=index(messtmp,ch10)
   end do
   write(unit, '(a)' ) trim(messtmp)
  else
   write(unit,*)
  end if

  if( index(trim(message),'BUG') /= 0 )then
   write(unit, '(a)' ) '  Action : contact hybrid_model group.'
   write(unit,*)
  end if

  if( index(trim(message),'BUG') /= 0   .or. &
       & index(trim(message),'Calculation completed') /= 0 )then
   if(nwarning<10000 .and. ncomment<1000)then
    write(unit, '(a,i5,a,i4,a)' ) &
         &     '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
   else
    write(unit, '(a,i6,a,i6,a)' ) &
         &     '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
   end if
   if(iexit/=0)then
    write(unit, '(a)' ) ' Note : exit requested by the user.'
   end if
  end if

  if( index(trim(message),'Exit') /= 0 )  iexit=1

  !Count the number of warnings and comments. Only take into
  !account unit 6, in order not to duplicate these numbers.
  if( index(trim(message),'WARNING') /= 0 .and. unit==6 )then
   nwarning=nwarning+1
  end if
  if( index(trim(message),'COMMENT') /= 0 .and. unit==6 )then
   ncomment=ncomment+1
  end if

 end subroutine wrtout_myproc
 


 !=============================================================
 !!subroutine: hybrid_model/wrtout
 !! NAME
 !!  wrtout
 !!
 !! FUNCTION
 !!  Organizes the sequential or parallel version of the write intrinsic
 !!  Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
 !!
 !!
 !! INPUTS
 !!  msg=(character(len=*)) message to be written
 !!  unit=unit number for writing
 !!  mode_paral=
 !!   'COLL' if all procs are calling the routine with the same message to be written once only 
 !!   'PERS' if the procs are calling the routine with different messages each to be written, 
 !!          or if one proc is calling the routine
 !!
 !! OUTPUT
 !!  (only writing)
 subroutine wrtout(unit,msg,mode_paral)

  integer,intent(in) :: unit
  character(len=4),intent(in) :: mode_paral
  character(len=*),intent(in) :: msg
  integer :: rtnpos
  character(len=7) :: tag
  character(len=len(msg)) :: my_msg
  character(len=len(msg)+50) :: string
  !Variables introduced for MPI version
  integer,save :: master=0
  integer :: me

  !******************************************************************

  !Be careful with the coding  of the parallel case ...

  !MG: Be careful**2, me and nproc are defined in MPI_COMM_WORLD.
  !One should pass the MPI communicator 

  !Determine who I am
  me = rang

  !msg is not changed therefore we can pass literal strings as well.
  my_msg = msg

  if( (mode_paral=='COLL') .or. (nproc==1) ) then
   if(me==master)  call wrtout_myproc(unit, my_msg)     
  elseif(mode_paral=='PERS') then

   if(me<10) then
    write(tag,'("-P-000",i1)') me
   elseif(me<100) then
    write(tag,'("-P-00",i2)') me
   elseif(me<1000) then
    write(tag,'("-P-0",i3)') me
   elseif(me<10000) then
    write(tag,'("-P-",i4)') me
   else
    tag=' ######'
   end if

   rtnpos=index(my_msg,ch10)
   do while(rtnpos/=0)
    write(string,'(3a)') tag, ' ', my_msg(1:rtnpos-1)
    write(unit,'(A)') trim(string)
    my_msg=my_msg(rtnpos+1:len(my_msg))
    rtnpos=index(my_msg,ch10)
   end do
   write(string, "(3a)") tag, ' ', my_msg
   write(unit,'(A)') trim(string)

  elseif(mode_paral=='INIT') then
   master=unit
  else
   write(string,'(7a)')ch10,&
        &   '  wrtout: ERROR -',ch10,&
        &   '  Unknown write mode: ',mode_paral,ch10,&
        &   '  Continuing anyway ...'
   write(unit, '(A)' ) trim(string)

  end if
 end subroutine wrtout
 

 !=============================================================
 !!subroutine: hybrid_model/wrt_double
 !! NAME
 !!  wrt_double
 !!
 !! FUNCTION
 !!  Double writing function screen and disk
 !!
 !!
 !! INPUTS
 !!  msg=(character(len=*)) message to be written
 !!  unit=unit number for writing
 !!  wrt_screen= 0 write (1 no write) on screen
 !!  wrt_disk= 0 write (1 no write) on disk
 !!
 !!
 !! OUTPUT
 !!  (only writing)

 subroutine wrt_double(unit,msg,wrt_screen,wrt_disk)

  integer,intent(in) :: unit
  integer,intent(in) :: wrt_screen,wrt_disk
  character(len=*),intent(in) :: msg

  !--writing on screen
  if (wrt_screen == 0) call wrtout(6,msg,'COLL') 

  !--writing on disk
  if (wrt_disk == 0 .and. unit/=6) call wrtout(unit,msg,'COLL') 

 end subroutine wrt_double


 !=============================================================
 !!subroutine: hybrid_model/wrt_debug
 !! NAME
 !!  wrt_debug
 !!
 !! FUNCTION
 !!  Double writing function screen and disk  
 !!  the debug in and out from subroutine
 !!
 !! INPUTS
 !!  msg=(character(len=*)) message to be written
 !!  unit=unit number for writing
 !!  wrt_screen= 0 write (1 no write) on screen
 !!  wrt_disk= 0 write (1 no write) on disk
 !!  inout_opt=0 in (1 out)
 !!
 !!
 !! OUTPUT
 !!  (only writing)
 subroutine wrt_debug(msg,inout_opt)
  integer,intent(in) :: inout_opt
  character(len=*),intent(in) :: msg
  integer :: ii
  character(len=500) :: msg_dbg
 
  if(inout_opt==0)    dbg_adv = dbg_adv + 1

  write(msg_dbg,'(40a)')('--',ii=1,dbg_adv),('  ',ii=1,20-dbg_adv)
  msg_dbg = trim(msg_dbg) // msg 

  if(inout_opt==0) then
   msg_dbg = trim(msg_dbg) // " (IN)"
  else
   msg_dbg = trim(msg_dbg) // " (OUT)"
  endif

  call wrt_double(6,msg_dbg,0,1)

  if(inout_opt==1)  dbg_adv = dbg_adv - 1

 end subroutine wrt_debug

end module m_writeout


