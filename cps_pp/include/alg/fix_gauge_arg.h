//------------------------------------------------------------------
/*!\file
  \brief Parameter container for the AlgFixGauge algorithm.

  $Id: fix_gauge_arg.h,v 1.3 2004-09-02 16:56:33 zs Exp $
*/
//------------------------------------------------------------------

#include<config.h>
CPS_START_NAMESPACE

#ifndef INCLUDED_FIX_GAUGE_ARG_H
#define INCLUDED_FIX_GAUGE_ARG_H

CPS_END_NAMESPACE
#include <util/enum.h>
#include <util/vector.h>
CPS_START_NAMESPACE

//! Container for gauge fixing parameters.
/*!  \ingroup algargs */
struct FixGaugeArg{
    
    FixGaugeType fix_gauge_kind; /*!< The kind of gauge fixing. */
    int hyperplane_start;
/*!<
  The global lattice coordinate of the first hyperplane of gauge fixing
  matrices. This is ignored for Landau gauge fixing
*/
    int hyperplane_step;
/*!<
  The coordinate step between hyperplanes on which the gauge is fixed.
  This is ignored for Landau gauge fixing.
*/
    int hyperplane_num;
/*!<
  The total number of the hyperplanes on which to fix the 
  gauge. If set to zero when Coulomb gauge is requested, then treated as
  a request to fix all hyperplanes in the global lattice.
  This is ignored for Landau gauge fixing.
*/

    Float stop_cond;               /*!< The stopping condition. */
    int max_iter_num;
/*!< Maximum number of iterations. If zero, then there is no maximum. */

};


#endif

CPS_END_NAMESPACE
