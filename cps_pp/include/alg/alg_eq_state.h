/*!\file
  \brief Definition of AlgEqState class..

  $Id: alg_eq_state.h,v 1.3 2004-09-02 17:00:10 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-09-02 17:00:10 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/alg_eq_state.h,v 1.3 2004-09-02 17:00:10 zs Exp $
//  $Id: alg_eq_state.h,v 1.3 2004-09-02 17:00:10 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_eq_state.h,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/alg_eq_state.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


#include<config.h>
CPS_START_NAMESPACE


#ifndef INCLUDED_ALG_EQ_STATE_H
#define INCLUDED_ALG_EQ_STATE_H

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/eq_state_arg.h>
CPS_START_NAMESPACE

//! Plaquette measurement 
/*!
  The normalised real trace of the average plaquette is measured in the
  hyperplane containing a specified direction and in the hyperplane
  perpendicular to that direction.

  On an anisotropic lattice, the specified direction must be the the
  anisotropic direction (why?).
*/
class AlgEqState : public Alg
{
 private:
    char *cname;

    EqStateArg *alg_eq_state_arg;
        // The argument structure for AlgEqState
 
    Float norm_fac;       
        // normalization factor

 public:
    AlgEqState(Lattice & latt, CommonArg *c_arg, EqStateArg *arg);

    virtual ~AlgEqState();

    //! Do the calculation.
    void run(void);
};



#endif





CPS_END_NAMESPACE
