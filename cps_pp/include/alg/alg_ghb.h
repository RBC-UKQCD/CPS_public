#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of the AlgGheatBath class.
  
  $Id: alg_ghb.h,v 1.3 2004-08-18 11:57:35 zs Exp $
*/
//------------------------------------------------------------------

#ifndef INCLUDED_ALG_GHB_H
#define INCLUDED_ALG_GHB_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/ghb_arg.h>
CPS_START_NAMESPACE

//-----------------------------------------------------------------------------
//! Class implementing the gauge field global heatbath algorithm.
/*!
  The algorithm used is the Cabbibo-Marinari SU(N) heatbath update with
  the Kennedy-Pendleton method for updating SU(2) subgroups

  \ingroup alg 
 */
//------------------------------------------------------------------
class AlgGheatBath : public Alg
{
 private:
    char *cname;

    GhbArg *alg_ghb_arg;
        // The argument structure for the GheatBath algorithm

    void relocate();
    void preserve_seed();
    void UpdateLink(Matrix* mp, const Matrix & stap);

 public:
    AlgGheatBath(Lattice & latt, CommonArg *c_arg, GhbArg *arg);

    virtual ~AlgGheatBath();

    void run(void);
    void NoCheckerBoardRun();
    void NodeCheckerBoardRun();
};

#endif





CPS_END_NAMESPACE
