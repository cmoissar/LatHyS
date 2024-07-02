!!=============================================================
!!=============================================================
!!module: hybrid_model/m_timing
!! NAME
!!  m_timing (MMancini)
!!
!! FUNCTION
!!  Module for the timing of quiet plasma
module m_timing

 use defs_basis
 use m_writeout
 use mpi
#include "q-p_common.h"

 implicit none
 private
 !!=============================================================
 !!--Type for count_exit_type 

 integer,parameter :: timer_key_length=100
 integer,parameter :: timer_arr_length=__TIMER_ARR_LENGTH

 !!=============================================================
 !!--Type for time-analysis of a single variable
 type,private :: entry_type

  !--The name of the entry
  character(len=timer_key_length) :: key=""

  !--The number of time the section of code has been executed.
  integer :: ncalls=0

  integer :: status=0
  ! 0 --> entry is not defined.
  ! 1 --> timing analysis is in execution.
  ! 2 --> timing analysis is completed.

  !--The CPU time of the entry
  real(dpo) :: cpu_time_tot=zero

  !--The wall time of the entry
  real(dpo) :: wall_time_tot=zero

  !--Temporary array used to store the reference CPU and WALL time.
  real(dpo) :: tzero(2)

 end type entry_type

 type(entry_type),pointer,dimension(:) :: time_arr => NULL()
 private::              &
      time_init_names,  &
      time_add,         &
      time_sprint,      &
      time_separator,   &
      time_partition    

 public::               &
      time_init,        &
      time_clean,       &
      time_get

contains
 !!############################################################

 !=============================================================
 !=============================================================
 !!module: m_timing/time_init
 !! NAME
 !!  time_init (MMancini)
 !!
 !! FUNCTION
 !!  Initialize the time variables
 subroutine time_init
  if(.not.(associated(time_arr))) allocate(time_arr(timer_arr_length))
  call time_init_names
 end subroutine time_init

 !=============================================================
 !=============================================================
 !!module: m_timing/time_init_names
 !! NAME
 !!  time_init_names (MMancini)
 !!
 !! FUNCTION
 !!  Initialize the time variables names
 subroutine time_init_names()
  ! if(.not.(associated(time_arr))) return
  time_arr(1)%key  = "Total  time" 
#ifdef HAVE_TIMING
  time_arr(2)%key  = " Allocation"
  time_arr(3)%key  = "     h3init"
  time_arr(4)%key  = "   re_start"
  time_arr(5)%key  = "      first"
  time_arr(6)%key  = "       cam3"
  time_arr(7)%key  = "     Dmtsp3"
  time_arr(8)%key  = "wrt_fields"
  time_arr(9)%key  = "wrt_particles"
  time_arr(10)%key = "wrt_tm_results"   
  time_arr(11)%key = "wrt_restart"   
  time_arr(12)%key = "     Emtsp3"   
  time_arr(13)%key = "       last"   
  time_arr(14)%key = ""   
  time_arr(15)%key = ""   
  time_arr(16)%key = ""   
  time_arr(17)%key = ""   
  time_arr(18)%key = ""   
  time_arr(19)%key = ""   
  !--in move    
  time_arr(20)%key = "move"   
  time_arr(21)%key = "pecalc"   
  time_arr(22)%key = "ecalc3"   
  time_arr(23)%key = "mvsp3r"
  time_arr(24)%key = "smth_func" 
  time_arr(25)%key = ""      
  time_arr(26)%key = ""   
  !--in mvsp3r
  time_arr(27)%key = "vcalc3"    
  time_arr(28)%key = "Amtsp3"    
  time_arr(29)%key = "xcalc3"   
  time_arr(30)%key = "pack_com"   
  time_arr(31)%key = "new_particles"   
  time_arr(32)%key = "Bmtsp3"   
  time_arr(33)%key = ""   
  !--pack_com
  time_arr(34)%key = "pack_com"   
  time_arr(35)%key = "pack_part"   
  time_arr(36)%key = "pre_communication"   
  time_arr(37)%key = "communication"   
  time_arr(38)%key = "rangement"   
  time_arr(39)%key = ""   
  !--cam3
  time_arr(41)%key = "calc_field"   
  time_arr(42)%key = "momad"   
  time_arr(43)%key = "energy_proc"   
  time_arr(44)%key = ""
  time_arr(45)%key = ""   
  !--calc_field
  time_arr(46)%key = "ecalc3"   
  time_arr(47)%key = "bcalc3"   
  time_arr(48)%key = ""   
  time_arr(49)%key = ""   
  !--test
  time_arr(50)%key = "test"   
  !--ecalc3
  time_arr(51)%key = "ecalc3"
  time_arr(52)%key = "dni_loop"   
  time_arr(53)%key = "b_eight"
  time_arr(54)%key = "JxB"       
  time_arr(55)%key = "gred_PE"   
  time_arr(56)%key = "curl_B"       
  time_arr(57)%key = "e_local"
  time_arr(58)%key = "cond_limit"   
  time_arr(59)%key = ""   
  time_arr(60)%key = ""   
  time_arr(61)%key = ""   
  time_arr(62)%key = ""  
  !--PBC+SMTH   
  time_arr(63)%key = "PBC"   
  time_arr(64)%key = "SMTH"   
  time_arr(66)%key = ""   
  time_arr(66)%key = ""   
  !--moments
  time_arr(67)%key = "momad3r"   
  time_arr(68)%key = "Amtsp3"   
  time_arr(69)%key = "Bmtsp3" 
  time_arr(70)%key = "Cmtsp3" 
  time_arr(71)%key = ""
  !--Exosphere and Photoproduction
  time_arr(72)%key = "exosphere"
  time_arr(73)%key = "photoproduction"
  time_arr(74)%key = "charge_exchange"
  time_arr(75)%key = "ionization"
  time_arr(76)%key = "ionosphere"
  time_arr(77)%key = "split"
  time_arr(78)%key = "feed ionosphere"
  time_arr(79)%key = ""  

#endif            
 end subroutine time_init_names

 !=============================================================
 !=============================================================
 !!module: m_timing/time_get
 !! NAME
 !!  time_get (MMancini)
 !!
 !! FUNCTION
 !!  Get the time for a entry
 !!
 !! IN 
 !!  timenume=entry number of time_arr
 !!  verb=what to do: 
 !!                  1 start timer
 !!                  2 stop  timer
 !!                  3 restart timer
 !!
 subroutine time_get(timenume,verb)
  
 ! use mpi
  
  integer,intent(in) :: timenume,verb
  real(dpo) :: loctime
  character(len=500) :: msg
  type(entry_type),pointer :: locentry =>NULL()

#ifndef HAVE_TIMING 
  if(timenume>1) return
#endif
  if(timenume>timer_arr_length) then
   write(msg,'(5a)')" ERROR: the dimension",ch10,&
        &        " time array is too small",ch10,&
        &        " Leaving now..."
   call wrtout(6,msg,'COLL')
   stop
  endif

  locentry => time_arr(timenume)

  loctime = MPI_WTIME()

  select case(verb)
  case(1)!--Start timer
   locentry%status = 1
   locentry%tzero(1) = loctime
   locentry%ncalls = locentry%ncalls+1
  case(2)!--Stop timer
   locentry%tzero(1) =  loctime - locentry%tzero(1)
   locentry%cpu_time_tot = locentry%cpu_time_tot + locentry%tzero(1)
   locentry%status = 2
   locentry%tzero(1) = zero    
  case(3)!--Restart timer
   locentry%status = 1
   locentry%ncalls = 1
   locentry%cpu_time_tot = zero
   locentry%tzero(1) = loctime

  case default
   write(msg,'(5a)')" ERROR: no good option",ch10,&
        &        " in subroutine time_get ",ch10,&
        &        " Leaving now..."
   call wrtout(6,msg,'COLL')
   stop
  end select
  nullify(locentry)
 end subroutine time_get

 !=============================================================
 !=============================================================
 !!module: m_timing/time_add
 !! NAME
 !!  time_add (MMancini)
 !!
 !! FUNCTION
 !!  Sum time on all proc  
 !!
 !! IN 
 !!  timenume=entry number of time_arr
 !! OUT 
 !!  locentry=overall time on all procs
 !!  std_time=standard sqrt(<time_proc**2>-<time_proc>**2)
 !! 
 subroutine time_add(timenume,locentry,std_time)

  use defs_mpitype,only     : mpiinfo
  !use mpi

  integer,intent(in)  :: timenume
  real(dpo),intent(out)    :: std_time
  type(entry_type),intent(out) :: locentry 

  integer :: code
  real(dpo) :: meantime
  character(len=500) :: msg

  if(timenume > timer_arr_length) then
   write(msg,'(5a)')" ERROR: the dimension",ch10,&
        &        " time array is too small",ch10,&
        &        " Leaving now..."
   call wrtout(6,msg,'COLL')
   stop
  endif

  !--Initialization
  locentry = time_arr(timenume)
  std_time = zero
  call MPI_REDUCE(time_arr(timenume)%cpu_time_tot,locentry%cpu_time_tot,&
       &          1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpiinfo%comm,code)
  call MPI_REDUCE((time_arr(timenume)%cpu_time_tot)**two,std_time,&
       &          1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpiinfo%comm,code)
  call MPI_REDUCE(time_arr(timenume)%ncalls,locentry%ncalls,&
       &          1,MPI_INTEGER,MPI_SUM,0,mpiinfo%comm,code)
  meantime = locentry%cpu_time_tot/real(mpiinfo%nproc,dpo)  
  if(locentry%cpu_time_tot>tol6) then  
   std_time = sqrt(std_time/real(mpiinfo%nproc,dpo)-meantime**two)
   std_time = std_time/meantime
  endif
 end subroutine time_add


 !=============================================================
 !=============================================================
 !!module: m_timing/time_separator
 !! NAME
 !!  time_separator (MMancini)
 !!
 !! FUNCTION
 !! Print separator between timing partition
 !!
 subroutine time_separator(namepar)

  character(len=*),intent(in) :: namepar
  integer :: ii
  character(len=150) :: msg

  write(msg,'(92a,4a)')&
       & (' ',ii=1,10),('-',ii=1,82),ch10,&
       & ' -',namepar,':'
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 end subroutine time_separator
 
 !=============================================================
 !=============================================================
 !!module: m_timing/time_partition
 !! NAME
 !!  time_partition (MMancini)
 !!
 !! FUNCTION
 !!  Print the time analysis for a partitiones routine
 !! IN
 !!  namepar(character) = name of the partition analysed
 !!  pos = entry of the analysed routine in time_arr    
 !!  posarr(:) = position of the subroutines which composed the
 !!              analysed one
 !!
 subroutine time_partition(namepar,pos,posarr)

  integer,intent(in) :: pos
  character(len=*),intent(in) :: namepar
  integer,intent(in),dimension(:) :: posarr
  integer :: ii,jj
  real(dpo) :: std_time,accum,percent
  character(len=500) :: msg
  type(entry_type) :: locentry,totentry

  call time_separator('Partition of '//namepar)
  call time_add(pos,totentry,std_time)
  accum = zero
  !--Print time analysis for move
  if(size(posarr)/=1 .or. posarr(1)/=0) then
   do jj=1,size(posarr)
    ii = posarr(jj)
    if(time_arr(ii)%status == 0) cycle
    call time_add(ii,locentry,std_time)
    percent = totentry%cpu_time_tot
    if(percent<tol10) percent = one
    write(msg,'(a32,f18.3,f16.2,i15,f10.2)')&
         & trim(locentry%key),&
         & locentry%cpu_time_tot,&
         & locentry%cpu_time_tot/percent*hundred,&
         & locentry%ncalls,&
         & std_time
    call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
    accum = accum+locentry%cpu_time_tot/percent
   enddo
  endif
  write(msg,'(a32,f18.3,f16.2,i15)')&
        & '++total of '//namepar,&
        & totentry%cpu_time_tot,&
        & accum*hundred,&
        & totentry%ncalls
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 end subroutine time_partition

 !=============================================================
 !=============================================================
 !!module: m_timing/time_sprint
 !! NAME
 !!  time_sprint (MMancini)
 !!
 !! FUNCTION
 !!  Print the Time analysis
 !!
 subroutine time_sprint

  integer  :: ii
  real(dpo) :: std_time
  character(len=500) :: msg
  type(entry_type) :: locentry

  write(msg,'(92a,a,36a,20a,36a,a,92a)')&
       & ('=',ii=1,92),ch10,&
       & ('=',ii=1,36),"   Time analysis    ",&
       & ('=',ii=1,36),ch10,&
       & ('=',ii=1,92)
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

#ifdef HAVE_TIMING
  write(msg,'(a32,a18,a16,a15,a10)')&
       & 'Sub Program',&
       & 'time(sec)',&
       & '%','calls','sdt'
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)

  call time_partition('quiet_plasma',1,(/(ii,ii=2,19)/))

  call time_partition('initialization',3,(/72,73,76/))

  call time_partition('cam3',6,(/20,41,42,43/))

  call time_partition('calc_field',41,(/21,46,47,48/))

  !call time_partition('ecalc3',51,(/52,53,54,55,56,57,58/))

  call time_partition('move',20,(/21,22,23,24/))

  call time_partition('mvsp3r',23,(/(ii,ii=27,32)/))

  call time_partition('vcalc3',27,(/74,77,78/))

  call time_partition('pack_com',34,(/36,35,37,38/))

  call time_partition('Moment on Q_P',1,(/67,68,69,70,7,12/))

  call time_partition('Comm on Q_P',1,(/36,37,63,64/))

  !call time_partition('test',50,(/0/))

#endif
  !--Print Total time 
  call time_add(1,locentry,std_time)
  write(msg,'(92a,a,a13,f34.3,a,a13,f34.3,f30.2,a,80a)')&
       & (' ',ii=1,10),('-',ii=1,82),ch10,&
       & 'Indiv. time',&
       & time_arr(1)%cpu_time_tot,ch10,&
       & trim(locentry%key),&
       & locentry%cpu_time_tot,&
       & std_time,ch10,&
       & ('=',ii=1,80)
  call wrt_double(qp_out,msg,wrtscreen,wrtdisk)
 end subroutine time_sprint

 !=============================================================
 !=============================================================
 !!module: m_timing/time_clean
 !! NAME
 !!  time_clean (MMancini)
 !!
 !! FUNCTION
 !!  Clean the time variables
 subroutine time_clean

  call time_sprint
  if(associated(time_arr))   deallocate(time_arr)
  nullify(time_arr)
 end subroutine time_clean

end module m_timing
