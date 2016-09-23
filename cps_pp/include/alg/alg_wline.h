#include<config.h>
CPS_START_NAMESPACE
//------------------------------------------------------------------
/*!\file
  \brief Definition of AlgWline class.

  $Id: alg_wline.h,v 1.3 2004/09/02 16:53:10 zs Exp $
*/
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_WLINE_H
#define INCLUDED_ALG_WLINE_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/no_arg.h>
CPS_START_NAMESPACE

//! Measures the %Wilson line in all directions.
/*!
  In each direction, the %Wilson line is calculated and averaged over the
  lattice.
  
  \ingroup alg
*/

class AlgWline : public Alg
{
 private:
    char *cname;

    NoArg *alg_wline_arg;
        // The argument structure for the plaquette
 
 public:
    AlgWline(Lattice & latt, CommonArg *c_arg, NoArg *arg);

    virtual ~AlgWline();

    //! Measures the %Wilson line in all directions.
    void run(void);
};



#endif





CPS_END_NAMESPACE
