#include <alg/gparity_contract_arg.h>
CPS_START_NAMESPACE

GparityContractArg::GparityContractArg(): conf_start(0), conf_incr(1), conf_lessthan(1){
  fix_gauge.fix_gauge_kind = FIX_GAUGE_NONE;
  fix_gauge.stop_cond = 1e-08;
  fix_gauge.max_iter_num = 10000;
};

CPS_END_NAMESPACE
