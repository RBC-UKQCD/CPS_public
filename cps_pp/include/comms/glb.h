#include<config.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Declarations of collective communications routines

  $Id: glb.h,v 1.4 2004-08-18 11:57:36 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:36 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/glb.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Id: glb.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/glb.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  glb.h
 */
#ifndef INCLUDED_GLB_H
#define INCLUDED_GLB_H

CPS_END_NAMESPACE
#include <util/vector.h>
#include <util/lattice.h>
#include <comms/nga_reg.h>
CPS_START_NAMESPACE

/*! \defgroup collectivecomms Collective communications routines
  \ingroup comms      
  @{ */  

//--------------------------------------------------------------
//! Sums a floating point number over all nodes in a 4-dimensinal grid.
//--------------------------------------------------------------
extern void glb_sum(Float * float_p);

//--------------------------------------------------------------
//! Sums a floating point number over all nodes in a 5-dimensinal grid.
// Sum over all nodes of the "virtual" 5-dimensional volume.
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes(), GJP.Snodes}
// Relevant for spread-out DWF (GJP.s_nodes not 1) only.
//--------------------------------------------------------------
extern void glb_sum_five(Float * float_p);

//--------------------------------------------------------------
//! Sums a floating point number over all nodes along a single direction.
//--------------------------------------------------------------
extern void glb_sum_dir(Float * float_p, int dir);
//! Sums a vector of floating point numbers over all nodes along a single direction.
extern void glb_sum_multi_dir(Float * float_p, int dir, int len);
//! Sums a Matrix over all nodes along a single direction.
extern void glb_sum_matrix_dir(Matrix * float_p, int dir);

//--------------------------------------------------------------
//! Finds the maximum floating point number over all nodes.
//--------------------------------------------------------------
extern void glb_max(Float * float_p);

//--------------------------------------------------------------
//! Finds the minimum floating point number over all nodes.
//--------------------------------------------------------------
extern void glb_min(Float * float_p);

//--------------------------------------------------------------
//! Sums a vector of floating point numbers over all nodes in a hyperplane
//--------------------------------------------------------------
extern void 
slice_sum(Float * float_p, int blcklength, int dir);

/*! @} */

#endif

CPS_END_NAMESPACE
