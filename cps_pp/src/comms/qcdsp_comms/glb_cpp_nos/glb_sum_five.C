#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp_nos/glb_sum_five.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: glb_sum_five.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.6  2001/08/16 12:54:10  anj
//  Some fixes follosin the float-> IFloat change, mostly of the (variable
//  anme) IFloat_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.5  2001/08/16 10:49:59  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.3  2001/07/03 17:00:53  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.2  2001/06/19 18:12:07  anj
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
//  $RCSfile: glb_sum_five.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp_nos/glb_sum_five.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
//  glb_sum_five
//
// Sum over all nodes of the "virtual" 5-dimensional volume.
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes(), GJP.Snodes}
// Relevant for spread-out DWF (GJP.s_nodes not 1) only.
//--------------------------------------------------------------


CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <sysfunc.h>
CPS_START_NAMESPACE


static Double64 transmit_buf;
static Double64 receive_buf;
static Double64 gsum_buf;
static IFloat *send_buf = (IFloat *) &transmit_buf;
static IFloat *rcv_buf = (IFloat *) &receive_buf;

static volatile unsigned* dsp_scu_base0x10 = 
        (volatile unsigned* )(DSP_SCU_BASE + 0x10);


void glb_sum_five(Float * float_p)
{
  int NP[5] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes(), GJP.Snodes()};


  // Sum over the "virtual" 5-dimensional mesh
  //------------------------------------------------------------
  gsum_buf = *float_p;

  int i;
  for(i = 0; i < 5; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
        bsm(send_buf, 2, 0, 1, gjp_scu_wire_map[2*i], TRANSMIT);
        bsm(rcv_buf, 2, 0, 1, gjp_scu_wire_map[2*i+1], RECEIVE);
        while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*i+1]) ) ;
        while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*i]) ) ;

        gsum_buf += receive_buf;
        transmit_buf = receive_buf;
      }
  }

  // Broadcast the result of node (0,0,0,0,0)
  //------------------------------------------------------------
  if(GJP.XnodeCoor() != 0 || 
     GJP.YnodeCoor() != 0 || 
     GJP.ZnodeCoor() != 0 || 
     GJP.TnodeCoor() != 0 || 
     GJP.SnodeCoor() != 0 ) {
    gsum_buf = 0;
  }

  for(i = 0; i < 5; ++i) {
    
      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
        bsm(send_buf, 2, 0, 1, gjp_scu_wire_map[2*i], TRANSMIT);
        bsm(rcv_buf, 2, 0, 1, gjp_scu_wire_map[2*i+1], RECEIVE);
        while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*i+1]) ) ;
        while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*i]) ) ;

        gsum_buf += receive_buf;
        transmit_buf = receive_buf;
      }
  }

  *float_p = gsum_buf;
}


CPS_END_NAMESPACE
