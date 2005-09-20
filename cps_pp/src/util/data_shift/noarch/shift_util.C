#include <config.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/data_shift.h>
CPS_START_NAMESPACE
void GlobalDataShift::Shift(int i, int n_disp){
  if (n_disp==0) return;
  ERR.NotImplemented(cname, "Shift(int, int)", "");
}
CPS_END_NAMESPACE
