!!===================================================
!!===================================================
module m_VO
 use defs_basis
 use defs_variable
#include "q-p_common.h"

 implicit none
 private 
 public ::  Get_VOTIME, getstri,getstrf

contains
subroutine get_VOTIME(YYYY,MM,DD,hh,minut,ss,votime)
 character(len=23),intent(out) :: votime
 integer,intent(in) :: YYYY,MM,DD,hh,minut
 real(dp),intent(in) :: ss   
 write(votime,'(i4,5(a,i2.2),a,i4.4)') YYYY,"-",MM,"-",DD,"T",hh,":",minut,":",int(ss),".",int(ss*1000.-int(ss)*1000.)
end subroutine Get_VOTIME

subroutine getstri(inp,outp,frmt)
 integer,intent(in) :: inp
 character(len=20),intent(out) :: outp
 character(len=*),intent(in) :: frmt
 write(outp,frmt) inp
end subroutine getstri

subroutine getstrf(inp,outp,frmt)
 real(dp),intent(in) :: inp
 character(len=*),intent(in) :: frmt
 character(len=20),intent(out) :: outp
 write(outp,frmt) inp
end subroutine getstrf

end module m_VO

