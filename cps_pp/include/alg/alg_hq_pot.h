#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief  Definitions of the AlgHQPotential class.

  $Id: alg_hq_pot.h,v 1.2 2007-06-05 15:44:19 chulwoo Exp $
*/
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_HQPOTENTIAL_H
#define INCLUDED_ALG_HQPOTENTIAL_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/fix_gauge_arg.h>
CPS_START_NAMESPACE

class AlgHQPotential : public Alg
{
 private:
    char *cname;

    NoArg *alg_HQPotential_arg;
        // The argument structure for the AlgFixGauge algorithm
 

 public:
    AlgHQPotential(Lattice & latt, CommonArg *c_arg, NoArg *arg);

    virtual ~AlgHQPotential();

    void run(int, int, int);  
    // run(direction, number of time seps, start time seps);
};

#endif

CPS_END_NAMESPACE
