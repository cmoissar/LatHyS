!!=============================================================
!!=============================================================
!!module: m_logo
!! NAME
!!  m_logo (MMancini)
!!
!! FUNCTION
!!  To print the logo of hybrid-model/quit-plasma
module m_logo

 use defs_basis
 use m_writeout
#include "q-p_common.h"

 implicit none
 private

 public ::          &
      take_date,    &!--Get data
      logo,         &!--To print the logo
      print_date     !--Print data on a string
 private ::         &
      compil_info    !--Print infos concerning compilation
contains
 !!############################################################

 !!=============================================================
 !!subroutine: m_logo/logo
 !! FUNCTION 
 !!  Only estetic write of the logo  
 !!
 !! INPUTS
 !!  unit=unit number for writing
 !!  wrt_disk=if 0 write on disk
 !!  wrt_screen=if 0 write on screen
 !!
 !! OUTPUT
 !!  (only writing)
 subroutine logo(unit,date_ok)

  integer,intent(in) :: unit,date_ok
  integer,dimension(8) :: valeur1
  character(len=500) :: msg


  !--Logo date et heure de début du run
  valeur1 = take_date()
  write(msg,'(3a)')&
       &' ____________________________________________________________',ch10,&
       &' ____________________________________________________________'
  call wrt_double(unit,msg,wrtdisk,wrtscreen)
  write(msg,'(16a)')&
       &'||        ||     __      || ',ch10,&
       &'|| \\\\  || ||   //  \\\\    || ',ch10,&
       &'||  \\\\ || ||       ||    || ',ch10,&
       &'||   \\\\|| ||       //    || ',ch10,&
       &'||\\\\  \\\|| ||\\\\  ===|   //|| ',ch10,&
       &'|| \\\\  || || \\\\    \\\\ // || ',ch10,&
       &'|| ||  || || //    || \\\\ || ',ch10,&
       &'|| ||  || ||// \\\\__//  \\\\|| '
  call wrt_double(unit,msg,wrtdisk,wrtscreen)
  write(msg,'(9a)')&
       &'                ||||  Hyb 3D                                ',ch10,&
       &'   An Hybrid-Model 3D-code for Mars/Solar Wind Interaction  ',ch10,&
       &'   MARS3D : periodic bcs in Y and Z  open bcs in X          ',ch10,&
       &' ____________________________________________________________',ch10,& 
       &' ____________________________________________________________'
  call wrt_double(unit,msg,wrtdisk,wrtscreen)
  if(date_ok==1)then
   write(msg,'(a,5(i2.2,a),i4.4)')&
        &'                   start time : '&
        &,valeur1(5),':',valeur1(6),':',valeur1(7),&
        '  date: ',valeur1(3),'/',valeur1(2),'/',valeur1(1)
   call wrt_double(unit,msg,wrtdisk,wrtscreen)
  end if

  call compil_info(unit)

  write(msg,'(2a,3(2a),a)')ch10,ch10,&
       &'   Known problems:',ch10,&
       &'   1) ',ch10,&
       &'   2) ',ch10,&
       &'      '
  call wrt_double(unit,msg,wrtdisk,wrtscreen)

 end subroutine logo

 !!=============================================================
 !!subroutine: m_logo/take_date
 !! FUNCTION 
 !!  Get the data in a integer(8) array
 !!
 !! OUTPUT
 !!  valeur1=data and time in a 8-array
 function take_date()

  integer,dimension(8) :: take_date
  character(len = 8)   :: date1
  character(len = 10)  :: temps1
  character(len = 5)   :: zone1

  call date_and_time(date1,temps1,zone1,take_date)

 end function take_date

 !!=============================================================
 !!subroutine: m_logo/print_date
 !! FUNCTION 
 !!  Get the data in a integer(8) array
 !!
 !! OUTPUT
 !!  date_str=character containig data
 function print_date()

  !  integer,intent(in) :: len_date_str
  character(len=8) :: print_date

  integer,dimension(8) :: valeur1

  character(len=4) ::msg

  valeur1 = take_date()

  write(msg,'(i4)')valeur1(1)
  write(print_date,'(i2.2,a1,i2.2,a1,a2)')&
       & valeur1(3),'_',&   !--Day
       & valeur1(2),'_',&   !--Month
       & trim(msg(3:))      !--Year

 end function print_date

 !!=============================================================
 !!subroutine: m_logo/compil_info
 !! FUNCTION 
 !!  Print information concerning Compilation, and Revision
 !!  at compiling time
 subroutine compil_info(unit)

  integer,intent(in) :: unit

  integer ::  bzrrevno
  character(len=100) :: compiler,bzr_msg,dp_msg
  character(len=100) :: wt_msg,npt_msg,cdf_msg,deb_msg
  character(len=500) :: msg

  !--Print information about Bzr revision
#ifdef BZR_REVNO
  bzrrevno = BZR_REVNO
#else
  bzrrevno = -1
#endif
  if(bzrrevno/=-1)then
   write(bzr_msg,'(a,i6)')&
        &"   Last BZR revision number:                        ",&
        & bzrrevno
  else
   write(bzr_msg,'(a)')&
        &"   No BZR revision "
  endif


  !--Test if double or single precision
  write(dp_msg,'(a13,a45)')"   Precision:",&
#ifdef HAVE_DOUBLE_PRECISION
       &                      "double"
#else
       &                      "simple"
#endif

  !--Test if have wavetest
  write(wt_msg,'(a18,a40)')"   HAVE_WAVE_TEST:",&
#ifdef HAVE_WAVE_TEST
       &                      "   Yes"
#else
       &                      "    No"
#endif

  !--Test if have no_planet
  write(npt_msg,'(a18,a40)')"   HAVE_NO_PLANET:",&
#ifdef HAVE_NO_PLANET
       &                      "   Yes"
#else
       &                      "    No"
#endif

  !--Test if have no_planet
  write(cdf_msg,'(a15,a43)')"   HAVE_NETCDF:",&
#ifdef HAVE_NETCDF
       &                      "   Yes"
#else
       &                      "    No"
#endif

  !--Test if have debug
  write(deb_msg,'(a14,a44)')"   HAVE_DEbUG:",&
#ifdef HAVE_DEBUG
       &                      "   Yes"
#else
       &                      "    No"
#endif

  !--Print information about Compiler if possible
  msg = ""
  compiler =&
       &"   Compiler:                                       Unknown  "

#ifdef __GFORTRAN__
  write(compiler,'(a,a)')&
       &"   Compiler:                         gfortran (gcc)  ",&
       & __VERSION__
#endif      

#ifdef __INTEL_COMPILER 
  write(compiler,'(a,f6.2)')&
       &"   Compiler:                         intel fortran  ",&
       & real(__INTEL_COMPILER,dp)/100.
#endif

  write(msg,'(18a)')ch10,&
       &  trim(bzr_msg),ch10,&
       &  trim(compiler),ch10,&
       &  trim(dp_msg),ch10,&
       &  trim(npt_msg),ch10,&
       &  trim(wt_msg),ch10,&
       &  trim(cdf_msg),ch10,&
       &  trim(deb_msg),ch10,&
       &"   Compiling date:                             ",__DATE__

  call wrt_double(unit,msg,wrtdisk,wrtscreen)

 end subroutine compil_info


end module m_logo
