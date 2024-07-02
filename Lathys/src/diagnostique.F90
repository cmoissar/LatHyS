!!=============================================================
!!=============================================================
module diagnostique

 use diag_fields
 use diag_particles
 use diag_iono
 use diag_prod
 use diag_tm_results
 use diag_flux_part_imp
 use diag_moment_species
 use m_distribution_function

#include "q-p_common.h"

 implicit none
 private

 public ::                       &
      diag_all

contains
 !!############################################################

 subroutine diag_all(filwrtname)
  
  use defs_variable,only : r_tm

  character(len=*),intent(in) :: filwrtname
  character(len=5) :: time_str, treg_str
  integer :: jreg

   call wrt_fields(filwrtname)
   call wrt_iono(filwrtname)
   call wrt_prod(filwrtname)
   call wrt_particles(filwrtname)
   call wrt_tm_results(filwrtname)
   call wrt_moment_species(filwrtname)
   if (distrib_activated.eq.1) call wrt_distribution_function(filwrtname)
   !if (diag_part_imp.ne.0) call wrt_flux_part_imp(filwrtname)

   !if the run is heavy, saving particles only at specified times may help
   !if (210==iter*dt) then
   !    call wrt_particles(filwrtname)
   !endif


   if (mpiinfo%me==0) then
     write (time_str, '(i5.5)') int(iter*dt)
     call system('mkdir t' // time_str)
     call system('mv *_t' // time_str // '.* t' // time_str // '/' )
 
    
! In case some files stubbornly stayed in the ncfiles directory
! (which especially happens during restarts)
     do jreg=0,r_tm%nreg
       write(treg_str, '(i5.5)') int(nhm*dt/r_tm%nreg*jreg)
       call system('mv *3*' // treg_str // '.nc t' // treg_str // '/')
     enddo 

     call system('tar -cf t' // time_str // '/Pack' // time_str // '.tar t' // time_str // '/*.nc')
     call system('rm t' // time_str // '/*.nc')

   endif



end subroutine diag_all

end module diagnostique
