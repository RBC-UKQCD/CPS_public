#include<config.h>

// If not MPI then must be either QCDOC or QCDSP
#if TARGET == QCDOC
#include <sysfunc.h>
#elif TARGET == cpsMPI
#include <comms/sysfunc_mpi.h>
#else
#include <comms/sysfunc_noarch.h>
#endif // TARGET




