#include<config.h>

/*!\file
  \brief Declarations for the communications layer.
  
  $Id: sysfunc.h,v 1.13 2006-11-25 19:09:47 chulwoo Exp $
*/

#ifdef USE_QMP
#include <comms/sysfunc_qmp.h>
#elif TARGET == QCDOC || TARGET == QCDSP
#include <sysfunc.h>
#elif TARGET == cpsMPI
#include <comms/sysfunc_mpi.h>
#else
#include <comms/sysfunc_noarch.h>
#endif // TARGET




