#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Definitions of the Alg class.

  $Id: alg_base.h,v 1.3 2004/08/18 11:57:35 zs Exp $
*/
//------------------------------------------------------------------
//
// alg_base.h
//
// Header file for the base alg class.
// The type of glue or fermion is given as
// an argument of type Lattice& to the constructor.
//
//------------------------------------------------------------------


#ifndef INCLUDED_ALG_BASE_H
#define INCLUDED_ALG_BASE_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/common_arg.h>
CPS_START_NAMESPACE

/*! \defgroup alg Algorithms */
//------------------------------------------------------------------
//! Alg is the base abstract class for all algorithms.
//------------------------------------------------------------------
class Alg
{
 private:
    char *cname;

    Lattice& alg_lattice;
        // Local reference to the Lattice object

 protected:
    CommonArg *common_arg;
        //!< The common argument structure for all algorithms


 public:
    Alg(Lattice& latt, CommonArg *c_arg);
    //!< Basic initialisation.

    virtual ~Alg();

    Lattice& AlgLattice();
        //!< Returns the lattice,
};


#endif





CPS_END_NAMESPACE
