#include<config.h>

/*!\file
  \brief Declarations for the communications layer.
  
  $Id: sysfunc.h,v 1.12 2004-09-03 15:49:26 zs Exp $
*/

#if TARGET == QCDOC || TARGET == QCDSP
#include <sysfunc.h>
#elif TARGET == cpsMPI
#include <comms/sysfunc_mpi.h>
#else
#include <comms/sysfunc_noarch.h>
#endif // TARGET




