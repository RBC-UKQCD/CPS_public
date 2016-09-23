#ifndef _QLA_H
#define _QLA_H

/* The following macros specify the prevailing color and precision */
/* Defaults to SU(3) single precision */
#ifndef QLA_Nc
#define QLA_Nc 3
#endif

#ifndef QLA_Precision
#define QLA_Precision 'F'
#endif

/* allow numeric precision specification */
#if QLA_Precision == 1
#undef QLA_Precision
#define QLA_Precision 'F'
#endif
#if QLA_Precision == 2
#undef QLA_Precision
#define QLA_Precision 'D'
#endif

/* Apply rule for default color namespace, based on QLA_Nc */
#if QLA_Nc == 2
#if QLA_Colors != 'N'
#define QLA_Colors 2
#endif
#elif QLA_Nc == 3
#if QLA_Colors != 'N'
#define QLA_Colors 3
#endif
#else
#define QLA_Colors 'N'
#endif

/* Define specific types and map them to generic types */
#include <qla_types.h>

/* Random number functions */
#include <qla_random.h>

/* Prototypes for complex functions */
#include <qla_cmath.h>

/* Headers we always define */
#include <qla_int.h>
#include <qla_df.h>
#include <qla_dq.h>

/* Headers we define regardless of precision */
#if ( QLA_Colors == 3 )
#include <qla_df3.h>
#include <qla_dq3.h>
#elif ( QLA_Colors == 2 )
#include <qla_df2.h>
#include <qla_dq2.h>
#elif ( QLA_Colors == 'N' )
#include <qla_dfn.h>
#include <qla_dqn.h>
#else
NONSTANDARD_QLA_Colors
#endif

/* Headers we define regardless of color */
#if ( QLA_Precision == 'F' )
#include <qla_f.h>
#elif ( QLA_Precision == 'D' )
#include <qla_d.h>
#elif ( QLA_Precision == 'Q' )
#include <qla_q.h>
#else
NONSTANDARD_QLA_Precision
#endif

/* Headers specific to color and precision */

#if ( QLA_Precision == 'F' )

#if   ( QLA_Colors == 3 )
#include <qla_f3.h>
#elif ( QLA_Colors == 2 )
#include <qla_f2.h>
#elif ( QLA_Colors == 'N' )
#include <qla_fn.h>
#endif

#elif ( QLA_Precision == 'D' )

#if   ( QLA_Colors == 3 )
#include <qla_d3.h>
#elif ( QLA_Colors == 2 )
#include <qla_d2.h>
#elif ( QLA_Colors == 'N' )
#include <qla_dn.h>
#endif

#elif ( QLA_Precision == 'Q' )

#if   ( QLA_Colors == 3 )
#include <qla_q3.h>
#elif ( QLA_Colors == 2 )
#include <qla_q2.h>
#elif ( QLA_Colors == 'N' )
#include <qla_qn.h>
#endif

#endif /* QLA_Precision */

#endif /* _QLA_H */
