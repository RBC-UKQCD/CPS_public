#include<config.h>

/*!\file
  \brief Declarations for the communications layer.
  
  $Id: sysfunc_cps.h,v 1.4 2012-03-27 21:17:55 chulwoo Exp $
*/

#ifdef USE_QMP
#include <comms/sysfunc_qmp.h>
#elif TARGET == BGL
#include <comms/sysfunc_bgl.h>
#else
#include <comms/sysfunc_noarch.h>
#endif // TARGET




