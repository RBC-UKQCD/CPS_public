#include <config.h>
#include <util/gjp.h>


extern "C" void _mcleanup();
CPS_START_NAMESPACE
void Start(){
// DefaultSetup();
}
void End(){
  _mcleanup();
}
CPS_END_NAMESPACE
