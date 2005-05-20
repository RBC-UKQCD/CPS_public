#include <config.h>
#include <util/data_types.h>
CPS_START_NAMESPACE
const int MAX_BUF = 72;
void glb_sum_internal(Float * float_p, int dir,int len);
void glb_sum_internal2(Float * float_p, int ndir);
// sum_flag = 1 for integer sum, 0 for exclusive OR
void glb_sum_internal2(unsigned int * uint_p, int ndir, int sum_flag = 1);
CPS_END_NAMESPACE
