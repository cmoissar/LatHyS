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
  type(treg_type),intent(inout) :: reg_tm

  integer :: ireg,unit,factor
  character(len=200) :: msg
  
  if(present(ifout)) then 
   unit = 6
  else
   unit = qp_out
  endif
   

  reg_tm%nreg = 13  !--Number of diagnostics to do

  nullify(reg_tm%time);  allocate(reg_tm%time(0:reg_tm%nreg))
  
  nullify(reg_tm%file);  allocate(reg_tm%file(0:reg_tm%nreg))

  !--Temps des diagnostiques
  reg_tm%time( 0) = 0._dp !--DO NOT CHANGE  0.
  
  !reg_tm%time( 1) = 610._dp
  !reg_tm%time( 2) = 620._dp
  !reg_tm%time( 3) = 630._dp
  !reg_tm%time( 4) = 645._dp
  !reg_tm%time( 5) = 655._dp
  !reg_tm%time( 6) = 665._dp
  !reg_tm%time( 7) = real(nhm,dp)*dt

  !reg_tm%time( 1) = real(nhm,dp)*dt/5.
  !reg_tm%time( 2) = real(nhm,dp)*dt/5.*2.
  !reg_tm%time( 3) = real(nhm,dp)*dt/5.*3.
  !reg_tm%time( 4) = real(nhm,dp)*dt/5.*4.
  !reg_tm%time( 5) = real(nhm,dp)*dt
  reg_tm%time( 1) = real(nhm,dp)*dt/13.
  reg_tm%time( 2) = real(nhm,dp)*dt/13.*2.
  reg_tm%time( 3) = real(nhm,dp)*dt/13.*3.
  reg_tm%time( 4) = real(nhm,dp)*dt/13.*4.
  reg_tm%time( 5) = real(nhm,dp)*dt/13.*5.
  reg_tm%time( 6) = real(nhm,dp)*dt/13.*6.
  reg_tm%time( 7) = real(nhm,dp)*dt/13.*7.
  reg_tm%time( 8) = real(nhm,dp)*dt/13.*8.
  reg_tm%time( 9) = real(nhm,dp)*dt/13.*9.
  reg_tm%time(10) = real(nhm,dp)*dt/13.*10.
  reg_tm%time(11) = real(nhm,dp)*dt/13.*11.
  reg_tm%time(12) = real(nhm,dp)*dt/13.*12.
  reg_tm%time(13) = real(nhm,dp)*dt


  !--Factor to put any file name >1
  factor = 1
  do while(int(factor*reg_tm%time(1))<1) 
   factor = factor *10
  end do

  if(maxval(reg_tm%time)>real(nhm,dp)*dt) then
   print *," ERROR : diagnostic time > max time"
   print *," To Solve: see in m_tregister.F90"
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
