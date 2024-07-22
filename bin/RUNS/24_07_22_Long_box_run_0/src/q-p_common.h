/*hyb_model/src/q-p_common.h*/

/*Macros for Hybrid model/quiet_plasma
Contains:
         macro for get timing
         macro for debug_printig
author  :MMancini
*/

#ifndef QP_COMMON_H
#  define QP_COMMON_H
#  include "configure.h"

/*Max size of real array which are written at time is machine
  dipendent we set it for ciclad cluster*/
#  define QP_MAXSIZEWRITE huge(1)

/*Macro for double or single precision variables*/
#  ifdef HAVE_DOUBLE_PRECISION
#    define QP_MPI_DP MPI_DOUBLE_PRECISION
#    define QP_NF90_DP nf90_double
#    ifdef HAVE_MAXSIZEWRITE
#      define QP_MAXSIZEWRITE 1500000
#    endif
#  else 
#    define QP_MPI_DP MPI_REAL
#    define QP_NF90_DP nf90_float
#    ifdef HAVE_MAXSIZEWRITE
#      define QP_MAXSIZEWRITE 3000000
#    endif
#  endif

/*Macro for time function*/
#  ifdef HAVE_TIMING
#    define __GETTIME(a,b) call time_get(a,b) 
#    define __TIMER_ARR_LENGTH 100
#  else 
#    define __GETTIME(a,b)
#    define __TIMER_ARR_LENGTH 1
#  endif

/*Macro for debug_printing function*/
#  ifdef HAVE_DEBUG
#    define __WRT_DEBUG_IN(a) call wrt_debug(a,0)
#    define __WRT_DEBUG_OUT(a) call wrt_debug(a,1)
#  else 
#    define __WRT_DEBUG_IN(a) 
#    define __WRT_DEBUG_OUT(a)
#  endif

/*Macro for no_planet function*/
#  ifdef HAVE_NO_PLANET
#    define __NO_PLANET_CYCLE return
#  else 
#    define __NO_PLANET_CYCLE
#  endif



# define WPUNTO(a)  print *,"PUNTO",a

#endif
