#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/utils/getseeds.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: getseeds.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.3  2003/02/24 12:10:35  anj
//  Fixed a nasty hack with a slightly less nasty one in glb_sum.C, fixed getseeds
//  so that it includes the global RNG object, and added new keychain and ssh
//  docs.
//
//  Revision 1.2  2001/08/01 15:17:30  anj
//  Minor changes to ensure painless compilarion.Anj
//
//  Revision 1.1  2001/07/31 10:12:56  anj
//  Added a directory for any testing utilities that are required.  In
//  particular, I have written a simple program to grab the seeds on
//  QCDSP, so that (near) binary agreement can be ensured across
//  platforms.  Anj.
//
//  Revision 1.6  2001/07/03 17:00:56  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.5  2001/06/29 12:04:17  anj
//
//  A few minor fixes and tests, but mostly a change in the I/O handling.
//  Off QCDSP, the I/O functions printf and fprintf are overriden by my
//  own qcdio.h library.  (This should eventually become part of the
//  general i/o spec.)  All this does is stop all processors from sending
//  out indentical output. Anj.
//
//  Revision 1.4  2001/06/22 13:06:04  anj
//  Minor alterations for testing and debugging. Anj
//
//  Revision 1.3  2001/06/21 09:20:34  anj
//  *** empty log message ***
//
//  Revision 1.2  2001/06/19 18:12:22  anj
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
//  Revision 1.2  2001/05/25 06:16:04  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: getseeds.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/utils/getseeds.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <stdio.h>
#include <stdlib.h>	// exit()
#include<config.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/alg_pbp.h>
#include<alg/do_arg.h>
#include<alg/common_arg.h>
#include<alg/pbp_arg.h>
CPS_START_NAMESPACE

#ifdef PARALLEL
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif

GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;

int main(int argc,char *argv[]){

  //----------------------------------------------------------------
  // Initializes all Global Job Parameters
  //----------------------------------------------------------------
  DoArg do_arg;

  do_arg.x_node_sites = 4;
  do_arg.y_node_sites = 4;
  do_arg.z_node_sites = 2;
  do_arg.t_node_sites = 4;
  do_arg.s_node_sites = 0;

#ifdef PARALLEL
  do_arg.x_nodes = 2;
  do_arg.y_nodes = 2;
  do_arg.z_nodes = 2;
  do_arg.t_nodes = 2;
  do_arg.s_nodes = 1;
#else
  do_arg.x_nodes = 1;
  do_arg.y_nodes = 1;
  do_arg.z_nodes = 1;
  do_arg.t_nodes = 1;
  do_arg.s_nodes = 1;
#endif

  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
  do_arg.start_seed_kind = START_SEED_FIXED;
  do_arg.colors = 3;
  do_arg.beta = 6.0;
  do_arg.dwf_height = 0.9;
  do_arg.verbose_level = -100502;

  GJP.Initialize(do_arg);


  //----------------------------------------------------------------
  // Set verbose level
  //----------------------------------------------------------------
  VRB.Level(GJP.VerboseLevel());

  //----------------------------------------------------------------
  // Get all of the seeds:
  //----------------------------------------------------------------

#ifdef PARALLEL

  printf("[%i,%i,%i,%i] %i %i %i %i\n",
	 CoorT(),CoorX(),CoorY(),CoorZ(),
	 Seed(),SeedS(),SeedT(),SeedST()
	 );

#else

  printf("[0,0,0,0] %i\n", SERIAL_SEED );

#endif

  return(0);
}




CPS_END_NAMESPACE
