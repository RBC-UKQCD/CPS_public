#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/comms/glb.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: glb.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.6  2001/08/16 12:54:13  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.5  2001/08/16 10:50:02  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:13  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:03  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: glb.h,v $
//  $Revision: 1.1.1.1 $
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
#include<util/vector.h>
#include<util/lattice.h>
#include<comms/nga_reg.h>
CPS_START_NAMESPACE

//--------------------------------------------------------------
// Sum over all nodes 
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()}
//--------------------------------------------------------------
extern void glb_sum(Float * float_p);

//--------------------------------------------------------------
// Sum over all nodes of the "virtual" 5-dimensional volume.
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes(), GJP.Snodes}
// Relevant for spread-out DWF (GJP.s_nodes not 1) only.
//--------------------------------------------------------------
extern void glb_sum_five(Float * float_p);

//--------------------------------------------------------------
// Sum over all nodes along a direction
// (0,1,2,3,4) <-> (x,y,z,t,s)
//--------------------------------------------------------------
extern void glb_sum_dir(Float * float_p, int dir);
extern void glb_sum_multi_dir(Float * float_p, int dir, int len);
extern void glb_sum_matrix_dir(Matrix * float_p, int dir);

//--------------------------------------------------------------
// Max over all nodes 
//--------------------------------------------------------------
extern void glb_max(Float * float_p);

//--------------------------------------------------------------
// Min over all nodes 
//--------------------------------------------------------------
extern void glb_min(Float * float_p);

//--------------------------------------------------------------
//  sum over all nodes on the 3D slice which is orthogonal to
//  the direction dir:
//
//  float_p[i] = sum_over_3D_nodes(float_p[i]) 
//  with 0 <= i < blcklength.
//
//  dir = 0, 1, 2, 3 corresponds to physics x, y, z, t
//--------------------------------------------------------------
extern void 
slice_sum(Float * float_p, int blcklength, int dir);

#endif
CPS_END_NAMESPACE
