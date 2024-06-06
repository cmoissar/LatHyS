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
  
  character(len=*),intent(in) :: filwrtname

   call wrt_fields(filwrtname)
   call wrt_iono(filwrtname)
   call wrt_prod(filwrtname)
   call wrt_particles(filwrtname)
   call wrt_tm_results(filwrtname)
   call wrt_moment_species(filwrtname)
   if (distrib_activated.eq.1) call wrt_distribution_function(filwrtname)
   !if (diag_part_imp.ne.0) call wrt_flux_part_imp(filwrtname)

end subroutine diag_all

end module diagnostique
