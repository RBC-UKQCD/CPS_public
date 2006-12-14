#include<config.h>
CPS_START_NAMESPACE

//--------------------------------------------------------------------
/*
 *  glb_sum_five.C
 */
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include<comms/glb.h>
#include<util/error.h>
CPS_START_NAMESPACE

void glb_sum_five(Float * float_p)
{
    VRB.Func("glb_sum_five","glb_sum_five");
    glb_sum(float_p);
}


CPS_END_NAMESPACE
