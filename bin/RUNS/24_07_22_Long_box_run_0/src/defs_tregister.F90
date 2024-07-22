!!=============================================================
!!=============================================================
!!module: hybrid_model/m_tregister
!! NAME
!!  m_tregister (MMancini)
!!
!! FUNCTION
!!  Module for Time Registers
module defs_tregister

 use defs_basis
 use defs_parametre
#include "q-p_common.h"

 implicit none
 private

 !!=============================================================
 !!--Type for diagnostic time. 
 type,public :: treg_type 
  integer  :: nreg  !--N of diagnostics
  real(dp),pointer :: time(:) !--Times of output files
  character(len = 60),pointer :: file(:) !--Names of output files
 end type treg_type

 public::              &
      set_tregister,   &!--Set the registers for diagnostic
      clean_tregister   !--Clean the registers for diagnostic
contains
 !!############################################################

 !!=============================================================
 !!routine: m_tregister/set_tregister
 !!
 !! FUNCTION
 !!  Intialise registers containing the time where to print the
 !!  diagnostic 
 !!
 !! IN
 !!  ifout=Optional if present and different from 0 then write also on 
 !!  qp_out file
 !!
 !! OUT
 !!  reg_tm%nreg=number of diagnostic
 !!  reg_tm%time(nr)=Diagnostic times
 !!  reg_tm%file(nr)=diagnostic files
 subroutine set_tregister(reg_tm,ifout)

  use m_writeout

  integer,optional,intent(in) :: ifout
  integer :: jreg
  type(treg_type),intent(inout) :: reg_tm

  integer :: ireg,unit,factor
  character(len=200) :: msg

  integer:: tmax

  !--When doing a restart, what's the time of the last output file you have?
  !--IF NOT DOING A RESTART, KEEP AT 10
  integer::t_last_dump=10

  integer::div1=5
  integer::until1=10
  integer::div2=5
  integer::until2=20
  integer::div3=5
  integer::until3=30
  integer::div4=5
  integer::until4=40
  integer::div5=5
  integer::until5

  integer::last_until
  integer::nreg

  tmax = nhm*dt
  until5 = tmax

  if(present(ifout)) then 
   unit = 6
  else
   unit = qp_out
  endif
  

  !--Number of diagnostics to do
  reg_tm%nreg = 1  !reg_tm%time(1)=t_last_dump

  last_until = t_last_dump
  if (until1 > last_until) then
    reg_tm%nreg = reg_tm%nreg + (until1 - last_until) / div1
    if (mod(until1 - last_until, div1) .ne. 0) reg_tm%nreg = reg_tm%nreg + 1
  endif
  last_until = max(until1, last_until)
  if (until2 > last_until) then
    reg_tm%nreg = reg_tm%nreg + (until2 - last_until)/div2
    if (mod(until2 - last_until, div2) .ne. 0) reg_tm%nreg = reg_tm%nreg + 1
  endif
  last_until = max(until2, last_until)
  if (until3 > last_until) then
    reg_tm%nreg = reg_tm%nreg + (until3 - last_until)/div3
    if (mod(until3 - last_until, div3) .ne. 0) reg_tm%nreg = reg_tm%nreg + 1
  endif
  last_until = max(until3, last_until)
  if (until4 > last_until) then
    reg_tm%nreg = reg_tm%nreg + (until4 - last_until)/div4
    if (mod(until4 - last_until, div4) .ne. 0) reg_tm%nreg = reg_tm%nreg + 1
  endif
  last_until = max(until4, last_until)
  if (until5 > last_until) then
    reg_tm%nreg = reg_tm%nreg + (until5 - last_until)/div5
    if (mod(until5 - last_until, div5) .ne. 0) reg_tm%nreg = reg_tm%nreg + 1
  endif
  
  nullify(reg_tm%time);  allocate(reg_tm%time(0:reg_tm%nreg))
  nullify(reg_tm%file);  allocate(reg_tm%file(0:reg_tm%nreg))

  !--Temps des diagnostiques
  reg_tm%time( 0) = 0._dp !--DO NOT CHANGE  0.
  reg_tm%time( 1) = t_last_dump
 
  do jreg=1, reg_tm%nreg-1

    if (reg_tm%time(jreg) + div1 <= until1) then
      reg_tm%time(jreg+1) = reg_tm%time(jreg) + div1
    endif
    if (reg_tm%time(jreg) < until1 .and. until1 < reg_tm%time(jreg)+div1) then
      reg_tm%time(jreg+1) = until1
    endif

    if (until1 >= until2) cycle
    if (until1 <= reg_tm%time(jreg) .and. reg_tm%time(jreg) + div2 <= until2) then
      reg_tm%time(jreg+1) = reg_tm%time(jreg) + div2
    endif
    if (reg_tm%time(jreg) < until2 .and. until2 < reg_tm%time(jreg)+div2) then
      reg_tm%time(jreg+1) = until2
    endif
   
    if (until2 >= until3) cycle
    if (until2 <= reg_tm%time(jreg) .and. reg_tm%time(jreg) + div3 <= until3) then
      reg_tm%time(jreg+1) = reg_tm%time(jreg) + div3
    endif
    if (reg_tm%time(jreg) < until3 .and. until3 < reg_tm%time(jreg)+div3) then
      reg_tm%time(jreg+1) = until3
    endif

    if (until3 >= until4) cycle
    if (until3 <= reg_tm%time(jreg) .and. reg_tm%time(jreg) + div4 <= until4) then
      reg_tm%time(jreg+1) = reg_tm%time(jreg) + div4
    endif
    if (reg_tm%time(jreg) < until4 .and. until4 < reg_tm%time(jreg)+div4) then
      reg_tm%time(jreg+1) = until4
    endif

    if (until4 >= until5) cycle
    if (until4 <= reg_tm%time(jreg) .and. reg_tm%time(jreg) + div5 <= until5) then
      reg_tm%time(jreg+1) = reg_tm%time(jreg) + div5
    endif
    if (reg_tm%time(jreg) < until5 .and. until5 < reg_tm%time(jreg)+div5) then
      reg_tm%time(jreg+1) = until5
    endif

  enddo

  !--Factor to put any file name >1
  factor = 1
  do while(int(factor*reg_tm%time(1))<1) 
   factor = factor *10
  end do

  if(maxval(reg_tm%time)>(real(nhm,dp)+1)*dt) then
   print *," ERROR : diagnostic time > max time",nhm,dt
   print *," To Solve: see in defs_tregister.F90"
   stop
  endif

  write(msg,'(a18,a15,a12,a10)')" Files:         ","    Output no. "," Iteration ","t  "
  call wrt_double(unit,msg,wrtscreen,wrtdisk)

  do ireg=0,reg_tm%nreg
   !--Noms des fichiers de diagnostiques
   write(reg_tm%file(ireg),'(a8,a2,i5.5)')trim(fildat),'_t',int(factor*reg_tm%time(ireg))
   
   !--Write files name, iteration and time 
   write(msg,'(a20,i12,i12,f10.1)') &
        & trim(reg_tm%file(ireg)), ireg,&
        & int(reg_tm%time(ireg)/dt),reg_tm%time(ireg)
   call wrt_double(unit,msg,wrtscreen,wrtdisk)
  enddo

 end subroutine set_tregister

 !!=============================================================
 !!routine: m_tregister/clean_tregister
 !!
 !! FUNCTION
 !!  Clean variables concerning time diagnostic
 !!         
 !! OUT
 !!  reg_tm%nreg=number of diagnostic set to zero
 subroutine clean_tregister(reg_tm)

  type(treg_type),intent(inout) :: reg_tm

  if(associated(reg_tm%time)) then
   deallocate(reg_tm%time)
   nullify(reg_tm%time)
  endif

  if(associated(reg_tm%file)) then
   deallocate(reg_tm%file)
   nullify(reg_tm%file)
  endif
  reg_tm%nreg = 0 

 end subroutine clean_tregister

end module defs_tregister
