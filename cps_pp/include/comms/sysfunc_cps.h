#include<config.h>

/*!\file
  \brief Declarations for the communications layer.
  
  $Id: sysfunc_cps.h,v 1.3 2008-02-08 22:01:05 chulwoo Exp $
*/

#ifdef USE_QMP
#include <comms/sysfunc_qmp.h>
#elif TARGET == QCDOC || TARGET == QCDSP
#include <sysfunc.h>
#elif TARGET == cpsMPI 
#include <comms/sysfunc_mpi.h>
#elif TARGET == BGL
#include <comms/sysfunc_bgl.h>
#else
#include <comms/sysfunc_noarch.h>
#endif // TARGET




