#include<config.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//-------------------------------------------------------------------
/*!\file
  \brief  Definition of glb_min and glb_max routine.

  $Id: glb_min_max.C,v 1.2 2004-03-15 03:26:46 cwj Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: cwj $
//  $Date: 2004-03-15 03:26:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_passthru/glb_min_max.C,v 1.2 2004-03-15 03:26:46 cwj Exp $
//  $Id: glb_min_max.C,v 1.2 2004-03-15 03:26:46 cwj Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.2.1  2004/03/11 17:09:40  cwj
//  *** empty log message ***
//
//  Revision 1.2  2004/01/13 20:39:06  chulwoo
//  Merging with multibuild
//
//  Revision 1.1.2.1.2.1  2003/11/06 20:22:20  cwj
//  *** empty log message ***
//
//  Revision 1.1.1.1  2003/11/04 05:05:03  chulwoo
//
//  starting again
//
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
//  Revision 1.5  2001/08/16 12:54:01  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.4  2001/08/16 10:49:47  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
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
//  $RCSfile: glb_min_max.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/glb_passthru/glb_min_max.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  glb_min_max.C (C++ version)
 */


CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include <comms/sysfunc.h>
CPS_START_NAMESPACE

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))


const SCUDir dir[] = { SCU_XP, SCU_XM, SCU_YP, SCU_YM,
                       SCU_ZP, SCU_ZM, SCU_TP, SCU_TM };



static Float *transmit_buf = NULL;
static Float *receive_buf = NULL;
static Float *gsum_buf = NULL;


//----------------------------------------------------------------------
/*!
  The number pointed to by \a float_p on all nodes is compared.
  \param float_p The number to be compared.
  \post The maximum is found and is written back to \a float_p, which is
  then identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 

void glb_max(Float * float_p)
{
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};

  if (transmit_buf == NULL) 
      transmit_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));
  if (receive_buf == NULL) 
      receive_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));
  if (gsum_buf == NULL) 
      gsum_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));
  *gsum_buf = *float_p;

  for(int i = 0; i < 4; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	SCUDirArg send(transmit_buf, dir[2*i], SCU_SEND, sizeof(Float));
	SCUDirArg rcv(receive_buf, dir[2*i+1], SCU_REC, sizeof(Float));

	send.StartTrans();
	rcv.StartTrans();
	send.TransComplete();
	rcv.TransComplete();

        *gsum_buf = max(*gsum_buf, *receive_buf) ;
        *transmit_buf = *receive_buf;
      }
  }
  *float_p = *gsum_buf;
}


//----------------------------------------------------------------------
/*!
  The number pointed to by \a float_p on all nodes is compared.
  \param float_p The number to be compared.
  \post The minimum is found and is written back to \a float_p, which is
  then identical on all nodes.

  \ingroup comms
*/
//---------------------------------------------------------------------- 
void glb_min(Float * float_p)
{
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};
  if (transmit_buf == NULL) 
      transmit_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));
  if (receive_buf == NULL) 
      receive_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));
  if (gsum_buf == NULL) 
      gsum_buf = (Float *)qalloc(QFAST|QNONCACHE,sizeof(Float));

  *gsum_buf = *float_p;

  for(int i = 0; i < 4; ++i) {

      *transmit_buf = *gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	SCUDirArg send(transmit_buf, dir[2*i], SCU_SEND, sizeof(Float));
	SCUDirArg rcv(receive_buf, dir[2*i+1], SCU_REC, sizeof(Float));

	send.StartTrans();
	rcv.StartTrans();
	send.TransComplete();
	rcv.TransComplete();

        *gsum_buf = min(*gsum_buf, *receive_buf) ;
        *transmit_buf = *receive_buf;
      }
  }
  *float_p = *gsum_buf;
}

CPS_END_NAMESPACE
