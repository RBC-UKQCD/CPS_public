#include<config.h>

/*!\file
  \brief Declarations for the communications layer.
  
  $Id: sysfunc.h,v 1.11 2004-09-03 12:34:14 zs Exp $
*/

#if TARGET == QCDOC || TARGET == QDCSP
#include <sysfunc.h>
#elif TARGET == cpsMPI
#include <comms/sysfunc_mpi.h>
#else
#include <comms/sysfunc_noarch.h>
#endif // TARGET




