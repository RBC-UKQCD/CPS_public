#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief  Definitions of the AlgPlaq class.

  $Id: alg_plaq.h,v 1.3 2004-08-18 11:57:35 zs Exp $
*/
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_PLAQ_H
#define INCLUDED_ALG_PLAQ_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/no_arg.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
//! A class implementing calculation of the average plaquette.
/*!
  This class computes the real trace of the plaquette averaged over
  the total number of plaquettes and the number of colours (which is three).
  Also computed is the variance of this mean, and one third of the real trace
  of the plaquette at the origin in the X-Y plane.

  \ingroup alg
*/
//------------------------------------------------------------------
class AlgPlaq : public Alg
{
 private:
    char *cname;

    NoArg *alg_plaq_arg;
        // The argument structure for the plaquette
 
    Float norm_fac;       
        // normalization factor

 public:
    AlgPlaq(Lattice & latt, CommonArg *c_arg, NoArg *arg);

    virtual ~AlgPlaq();

    void run(void);
};



#endif





CPS_END_NAMESPACE
