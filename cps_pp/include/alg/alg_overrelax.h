#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of the AlgOverRelax class.
  
  $Id: alg_overrelax.h,v 1.2 2004/09/21 18:07:14 chulwoo Exp $
*/
//------------------------------------------------------------------

#ifndef INCLUDED_ALG_OVERRELAX_H
#define INCLUDED_ALG_OVERRELAX_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/overrelax_arg.h>
CPS_START_NAMESPACE

//-----------------------------------------------------------------------------
//! Class implementing the gauge field global heatbath algorithm.
/*!
This is Creutz overrelaxation

  \ingroup alg 
 */
//------------------------------------------------------------------
class AlgOverRelax : public Alg
{
 private:
    char *cname;

    OverRelaxArg *alg_overrelax_arg;
        // The argument structure for the GheatBath algorithm
    Float fGamma;
    void kernel(  Float *sigma, Float *U);
    void UpdateLink(Matrix* mp, const Matrix & stap);

 public:
    AlgOverRelax(Lattice & latt, CommonArg *c_arg, OverRelaxArg *arg);

    virtual ~AlgOverRelax();

    void run(void);
   
};

#endif





CPS_END_NAMESPACE
