#include <config.h>
#include <stdio.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of ParallelTransport class methods.
  
  $Id: pt_base.C,v 1.3 2004-04-27 03:51:21 cwj Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: cwj $
//  $Date: 2004-04-27 03:51:21 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/pt_base.C,v 1.3 2004-04-27 03:51:21 cwj Exp $
//  $Id: pt_base.C,v 1.3 2004-04-27 03:51:21 cwj Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2.2.2  2004/04/26 02:14:12  cwj
//  *** empty log message ***
//
//  Revision 1.2.2.1  2004/04/01 05:20:25  cwj
//  *** empty log message ***
//
//  Revision 1.2  2004/01/13 20:39:53  chulwoo
//  Merging with multibuild
//
//  Revision 1.1.2.1  2003/11/06 21:04:06  cwj
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
//  Revision 1.2  2003/07/24 16:53:54  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.8  2002/03/11 22:27:02  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.5.2.1  2002/03/08 16:36:36  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.5  2001/08/16 12:54:29  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.4  2001/08/16 10:50:15  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:36  anj
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
//  $RCSfile: pt_base.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/pt_base.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------


CPS_END_NAMESPACE
#include <util/pt.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/dirac_op.h>
#include <util/error.h>
#include <util/gjp.h>
#include <comms/glb.h>
#include <comms/cbuf.h>
#include <stdlib.h>	// exit()
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Initialize static variables
//------------------------------------------------------------------
int ParTrans::scope_lock = 0;

//------------------------------------------------------------------
// static variables
//------------------------------------------------------------------
static int nx[4];
static int bc[4] = {0,0,0,0};	// boundary condition on this node


static void BondCond(Lattice& lat, Matrix *u_base)
{
  for(int u = 0; u < 4; ++u) {
    if(bc[u]) {
      int x[4];
      for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
        for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
          for(x[2] = 0; x[2] < nx[2]; ++x[2]) {
            for(x[3] = 0; x[3] < nx[3]; ++x[3]) {
	      if(x[u]==nx[u]-1) {
	        Matrix *m = u_base+lat.GsiteOffset(x)+u;
	        *m *= -1;
	      }
	    }
          }
	}
      }
    }
  }
}


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
  also converted: This assumes they are initially in the same order as
  the gauge field.
 */
//------------------------------------------------------------------
ParTrans::ParTrans(Lattice & latt) :           // Lattice object
		 lat(latt)
{
  cname = "ParTrans";
  char *fname = "ParTrans(L&,CnvFrmType)";
  VRB.Clock(cname,fname,"Just entered\n");
  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // Check and set the scope_lock
  //----------------------------------------------------------------
  if(scope_lock != 0){
    ERR.General(cname,fname,
		"Only one ParTrans object is allowed to be on scope\n");
  }
  scope_lock = 1;

  //----------------------------------------------------------------
  // Set the gauge field pointer
  //----------------------------------------------------------------
  gauge_field = lat.GaugeField();

  //----------------------------------------------------------------
  // Save the circular buffer control register values
  //----------------------------------------------------------------
  saveCbufCntrlReg();

  //----------------------------------------------------------------
  // Initialize boundary condition on this node
  //----------------------------------------------------------------
  bc[0] = 0;
  bc[1] = 0;
  bc[2] = 0;
  bc[3] = 0;
  if(GJP.Xbc() == BND_CND_APRD)
  	bc[0] = GJP.XnodeCoor() == (GJP.Xnodes()-1) ? 1 : 0;
  if(GJP.Ybc() == BND_CND_APRD)
  	bc[1] = GJP.YnodeCoor() == (GJP.Ynodes()-1) ? 1 : 0;
  if(GJP.Zbc() == BND_CND_APRD)
  	bc[2] = GJP.ZnodeCoor() == (GJP.Znodes()-1) ? 1 : 0;
  if(GJP.Tbc() == BND_CND_APRD)
  	bc[3] = GJP.TnodeCoor() == (GJP.Tnodes()-1) ? 1 : 0;

  nx[0] = GJP.XnodeSites();
  nx[1] = GJP.YnodeSites();
  nx[2] = GJP.ZnodeSites();
  nx[3] = GJP.TnodeSites();


  //----------------------------------------------------------------
  // turn on the boundary condition, if not inside a dirac operator
  //----------------------------------------------------------------
  if(DiracOp::scope_lock ==0)
  BondCond(latt, gauge_field);

  //???

  // Added in by Ping for anisotropic lattices
  //------------------------------------------------------------------
  // No data manipulations if the scaling factor is 1.0
  // [built into lat.MltFloat].
  {
    Float factor = GJP.XiVXi() / GJP.XiV();
    lat.MltFloat(factor, GJP.XiDir());
  }
}


//------------------------------------------------------------------
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
ParTrans::~ParTrans() {
  char *fname = "~ParTrans()";

  VRB.Func(cname,fname);

  //----------------------------------------------------------------
  // turn off the boundary condition
  //----------------------------------------------------------------
  if(DiracOp::scope_lock ==0)
  BondCond(lat, gauge_field);

  VRB.Clock(cname,fname,"Exiting\n");

  //----------------------------------------------------------------
  // Restore the circular buffer control register values
  //----------------------------------------------------------------
  restoreCbufCntrlReg();


  //----------------------------------------------------------------
  // Release the scope_lock
  //----------------------------------------------------------------
  scope_lock = 0;

  //???

  // Added in by Ping for anisotropic lattices
  //------------------------------------------------------------------
  // No data manipulations if the scaling factor is 1.0
  // [built into lat.MltFloat].
  {
    Float factor = GJP.XiV() / GJP.XiVXi();
    lat.MltFloat(factor, GJP.XiDir());
  }
}

CPS_END_NAMESPACE
