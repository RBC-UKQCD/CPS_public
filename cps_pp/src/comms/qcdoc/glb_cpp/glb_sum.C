#include<config.h>
#include<stdio.h>
#include<qcdoc_align.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief Definition of glb_sum routine.

  $Id: glb_sum.C,v 1.2 2004-01-13 20:39:06 chulwoo Exp $ 
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:06 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_cpp/glb_sum.C,v 1.2 2004-01-13 20:39:06 chulwoo Exp $
//  $Id: glb_sum.C,v 1.2 2004-01-13 20:39:06 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.2.1.2.3  2003/12/11 20:22:52  cwj
//  *** empty log message ***
//
//  Revision 1.1.2.1.2.2  2003/11/28 19:12:47  cwj
//  *** empty log message ***
//
//  Revision 1.1.2.1.2.1  2003/11/06 20:22:20  cwj
//  *** empty log message ***
//
//  Revision 1.1.1.1  2003/11/04 05:05:03  chulwoo
//
//  starting again
//
//
//  Revision 1.2  2003/10/10 21:30:28  chulwoo
//
//  added sizeof's to compensate different convention
//
//  Revision 1.1  2003/10/08 18:38:29  chulwoo
//  start from vanilla_comms
//  start from QCDSP comms
//
//  Revision 1.1.1.1  2003/09/18 22:30:43  chulwoo
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
//  Revision 1.1.1.1  2003/06/22 13:34:47  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.9  2003/02/24 12:10:35  anj
//  Fixed a nasty hack with a slightly less nasty one in glb_sum.C, fixed getseeds
//  so that it includes the global RNG object, and added new keychain and ssh
//  docs.
//
//  Revision 1.8  2001/08/16 12:54:01  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.7  2001/08/16 10:49:48  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.5  2001/07/19 12:58:18  anj
//  B
//  C
//  My new expression for determining the size of the Double64 was not working, as Double64 actually appears to be a pointer to a double, and so is four bytes instead of 8.  Fixed up for now, until we neCVS: ----------------------------------------------------------------------
//
//  Revision 1.4  2001/07/03 17:00:48  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/29 12:04:16  anj
//
//  A few minor fixes and tests, but mostly a change in the I/O handling.
//  Off QCDSP, the I/O functions printf and fprintf are overriden by my
//  own qcdio.h library.  (This should eventually become part of the
//  general i/o spec.)  All this does is stop all processors from sending
//  out indentical output. Anj.
//
//  Revision 1.2  2001/06/19 18:11:50  anj
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
//  Revision 1.2  2001/05/25 06:16:01  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: glb_sum.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_cpp/glb_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum
//
// Sum over all nodes 
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()}
//--------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc.h>
CPS_START_NAMESPACE



static Double64 transmit_buf PEC_ALIGN LOCATE("Edramnoncache");
static Double64 receive_buf PEC_ALIGN LOCATE("Edramnoncache");
static Double64 gsum_buf PEC_ALIGN LOCATE("Edramnoncache");

static int output = 1;


//----------------------------------------------------------------------
/*!
  \param float_p The number to be summed.
  \post The number pointed to by \a float_p is summed over all nodes
  and that sum is written back to \a float_p, which is identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 
void glb_sum(Float * float_p)
{
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};
  static int counter = 0;

  if (output) printf("glb_sum %d before = %e ", counter, (double)*float_p);
  gsum_buf = (Double64)*float_p;

  for(int i = 0; i < 4; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	SCUDirArg send(&transmit_buf, gjp_scu_dir[2*i], SCU_SEND, sizeof(Double64) );
	SCUDirArg rcv(&receive_buf, gjp_scu_dir[2*i+1], SCU_REC, sizeof(Double64) );

//	printf("sent=%e received %e\n",transmit_buf,receive_buf);
	send.StartTrans();
	rcv.StartTrans();
	send.TransComplete();
	rcv.TransComplete();

        gsum_buf += receive_buf;
        transmit_buf = receive_buf;
      }
  }
  *float_p = (Float)gsum_buf;
if (output)   printf("after = %e\n", (double)*float_p);
  counter++;
}
CPS_END_NAMESPACE
