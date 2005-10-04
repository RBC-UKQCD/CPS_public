/*!<

  Statistics held by each integrator with regards to the conjugate
  gradient iterations.

 */

#ifndef _CG_STATS
#define _CG_STATS

#include <config.h>
CPS_START_NAMESPACE

class CgStats {
public:
  int cg_calls;

  int cg_iter_total;
  int cg_iter_min;
  int cg_iter_max;

  Float cg_iter_av;

  Float true_rsd_total;
  Float true_rsd_av;
  Float true_rsd_min;
  Float true_rsd_max;

  CgStats() {
    init();
  }

  void init() {
    cg_calls = 0;
    cg_iter_total = 0;
    cg_iter_av = 0.0;
    cg_iter_max = 0;
    cg_iter_min = 1000000;

    true_rsd_total = 0.0;
    true_rsd_av = 0.0;
    true_rsd_min = 3.4e38;
    true_rsd_max = 0.0;
  }


};

CPS_END_NAMESPACE

#endif /* !_CG_STATS */
