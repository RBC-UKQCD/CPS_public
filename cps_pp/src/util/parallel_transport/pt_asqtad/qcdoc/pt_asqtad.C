#include <config.h>
#include <stdio.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpStag class methods.

  $Id: pt_asqtad.C,v 1.4 2004-04-27 03:51:21 cwj Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: cwj $
//  $Date: 2004-04-27 03:51:21 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_asqtad/qcdoc/pt_asqtad.C,v 1.4 2004-04-27 03:51:21 cwj Exp $
//  $Id: pt_asqtad.C,v 1.4 2004-04-27 03:51:21 cwj Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.3.2.1  2004/04/26 02:14:12  cwj
//  *** empty log message ***
//
//  Revision 1.3  2004/01/13 23:25:24  chulwoo
//  *** empty log message ***
//
//  Revision 1.1.2.3  2003/12/27 21:05:31  cwj
//
//  (somewhat) cleaned up for QCDOC + qos-1-8-5
//
//  Revision 1.1.2.2  2003/12/11 20:22:53  cwj
//  *** empty log message ***
//
//  Revision 1.1.2.1  2003/11/06 21:02:09  cwj
//  *** empty log message ***
//
//  Revision 1.1.1.1  2003/11/04 05:05:17  chulwoo
//
//  starting again
//
//
//  Revision 1.1.1.1  2003/09/18 22:30:55  chulwoo
//  Mike's files for single node QCDOC + Parallel transport
//  I added some hacks for PARALLEL without MPI_SCU
//  PARALLEL=2 set PARALLEL without MPI_SCU
//
//
//  Revision 1.3  2003/08/29 21:00:24  mike
//  Removed specific MatMInv function as not needed.
//
//  Revision 1.2  2003/07/24 16:53:54  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.7  2002/03/11 22:27:05  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.4.2.1  2002/03/08 16:36:41  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.4  2001/08/16 10:50:19  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:43  anj
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
//  Revision 1.2  2001/05/25 06:16:05  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: pt_asqtad.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_asqtad/qcdoc/pt_asqtad.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// pt_asqtad.C
//
// ParTransAsqtad is derived from the ParTransStagTypes class.
// ParTransAsqtad is the front end for a library that contains
// all Dirac operators associated with Staggered fermions.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/vector.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/pt.h>
#include <util/error.h>
#include <util/stag.h>
#include <comms/cbuf.h>
#include <comms/glb.h>
#include <math.h>
CPS_START_NAMESPACE



//------------------------------------------------------------------
/*!
  Only one instance of this class is allowed to be in existence at
  any time.
  \param latt The lattice on which this Dirac operator is defined
  \param f_field_out A (pointer to) a spin-colour field (optionally). 
  \param f_field_in A (pointer to) a spin-colour field (optionally).
  \param arg Parameters for the solver.
  \param cnv_frm_flag Whether the lattice fields should be converted to
  to a new storage order appropriate for the type of fermion action.
  are also converted: This assumes they are initially in the same order as
  the gauge field.
 */
//------------------------------------------------------------------
ParTransAsqtad::ParTransAsqtad(Lattice & latt) :
			 ParTransStagTypes(latt)
{
  cname = "ParTransAsqtad";
  char *fname = "ParTransAsqtad(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Do the necessary conversions
  //----------------------------------------------------------------
    lat.Convert(STAG);

  //----------------------------------------------------------------
  // Set the node checkerboard size of the fermion field
  //----------------------------------------------------------------
  f_size_cb = GJP.VolNodeSites() * lat.FsiteSize() / 2;

  pt_init(lat.GaugeField());
  pt_init_g();

#if 0
  //----------------------------------------------------------------
  // Allocate memory for the temporary fermion vector frm_tmp.
  //----------------------------------------------------------------
  frm_tmp = (Vector *) smalloc(f_size_cb * sizeof(Float));
  if(frm_tmp == 0)
    ERR.Pointer(cname,fname, "frm_tmp");
  VRB.Smalloc(cname,fname, "frm_tmp", 
	      frm_tmp, f_size_cb * sizeof(Float));
#endif

}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
ParTransAsqtad::~ParTransAsqtad() {
  char *fname = "~ParTransAsqtad()";
  VRB.Func(cname,fname);

//    lat.Convert(CANONICAL);

  //----------------------------------------------------------------
  // Free memory
  //----------------------------------------------------------------
  pt_delete_g();
  pt_delete();
}

CPS_END_NAMESPACE
