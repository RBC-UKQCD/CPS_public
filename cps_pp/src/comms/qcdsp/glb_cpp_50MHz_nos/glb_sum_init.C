#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:12 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_50MHz_nos/glb_sum_init.C,v 1.2 2004-01-13 20:39:12 chulwoo Exp $
//  $Id: glb_sum_init.C,v 1.2 2004-01-13 20:39:12 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.2.1  2003/11/05 18:12:13  mike
//  Changing directory structure.
//
//  Revision 1.1.1.1  2003/06/22 13:34:47  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.2  2001/06/19 18:11:59  anj
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
//  Revision 1.2  2001/05/25 06:16:02  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: glb_sum_init.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_50MHz_nos/glb_sum_init.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//
// glb_sum_init.C
//

CPS_END_NAMESPACE
#include<comms/glb_sum_init.h>
CPS_START_NAMESPACE

//--------------------------------------------------------------------------
// Constants used by glb_sum.asm only.
//--------------------------------------------------------------------------
extern const unsigned int GLB_SUM_INFO_BUF = 0;
extern const int GLB_SUM_MAX_TRY = 0;


//--------------------------------------------------------------------------
// !!!! IMPORTANT !!!!   MUST BE THE SAME AS THAT IN glb_sum.asm.
//--------------------------------------------------------------------------
struct GsumInfo
{
  int not_master;       /* the flag indicating whether this node is the 
			   master node. The master node has the final
			   result during a bare Sum/Max. */
  int max_mode;         /* the pre-determined value for the passthrough
			   control register if Max desired. Max and Sum
			   share the same wire selection. */
  int brdcst_mode;      /* the pre-determined value for the passthrough
			   control register if broadcast desired. For the
			   master node it however holds the wires selection
			   since the register has to be 0. */
  int send_wire;        /* the special regular send wire during broadcast */
  int recv_wire;        /* the special regular recv wire during max/sum */
  int nbits;            /* number of bits reserved for overflow during
			   global sum. It equals the rounded up base 2 log
			   of the total number of nodes */
};
 

//--------------------------------------------------------------------------
// void glb_sum_init(...)
//--------------------------------------------------------------------------
void glb_sum_init(const int masterMachineCoord[], 
		  const int thisMachineCoord[], 
		  const int size[],
		  unsigned int sixWords,  int maxTRY)
{	
  //------------------------------------------------------------------------
  // The following coordinates are in MACHINE order [0/1, 2/3, 4/5, 6/7].
  // Odd wires are the positive directions.
  //------------------------------------------------------------------------
  int maxCoord[4], thisCoord[4], masterCoord[4];
  int iDim;    
  int send_wire, recv_wire; 
  struct GsumInfo *buffer = (struct GsumInfo *)sixWords;
    
  //--------------------------------------------------------------------------
  // define a union for the convenience of setting the passthrough registers
  //--------------------------------------------------------------------------
  union {	
    int _intval;
    struct     {
      unsigned int brdcst_src  :8;            // the least significant bits 
      unsigned int cmbn_src    :8;
      unsigned int final_carry :7;
      unsigned int add_mode    :1;
    } _bitval;
  } brdcst_mode, max_mode;

  //--------------------------------------------------------------------------
  // set the constants for glb_sum.asm
  //--------------------------------------------------------------------------
  { 
    unsigned int *buffer_p = (unsigned int *)&GLB_SUM_INFO_BUF;
    *buffer_p = sixWords;    

    int * maxTRY_p = (int *)&GLB_SUM_MAX_TRY;
    *maxTRY_p = (maxTRY > 1 ? maxTRY : 1);
  }

  
  //--------------------------------------------------------------------------
  // initialize sum_mode, brdcst_mode
  //--------------------------------------------------------------------------
  max_mode._intval = brdcst_mode._intval = 0;

  for( iDim = 0; iDim < 4; ++iDim)    {
    thisCoord[iDim]  = thisMachineCoord[iDim];
    masterCoord[iDim] = masterMachineCoord[iDim];
    maxCoord[iDim]   = size[iDim] - 1;    
  }
  
  //--------------------------------------------------------------------------
  // translate the coordinates by the amount of maxCoord[i]/2 - masterCoord[i]
  // so that the master is in the scu network center.
  //--------------------------------------------------------------------------
  for( iDim = 0; iDim < 4; ++iDim)    {
    thisCoord[iDim] += maxCoord[iDim] / 2 - masterCoord[iDim];        
    thisCoord[iDim]  = (thisCoord[iDim] + size[iDim]) % size[iDim];  
    // plus size[iDim] because operator % on a negative number is undefined.
    masterCoord[iDim]  = maxCoord[iDim] / 2;    
  }

  //-------------------------------------------------------------------------
  // Compare coordinates of this node with those of the master node until
  // a difference occurs from wire 6, 7 to wire 0, 1.
  // This order is chosen to take advantage that
  // number of nodes in wire 0,1 is fixed to 4 and likely is the smallest.
  // If the difference occurs say in coord23, wires 7,6,5,4 will be
  // selected to participate in passthrough mode unless this node 
  // is on the boundary in these directions.
  //-------------------------------------------------------------------------
  for( iDim = 3; iDim >= 0; --iDim )  {
    if (iDim == 0 && thisCoord[iDim] == masterCoord[iDim]) { // master node
      send_wire = recv_wire = 1;                // recv_wire will not be used
      buffer->not_master = 0;
      max_mode._bitval.cmbn_src      |= (1<<1);  
      brdcst_mode._bitval.brdcst_src |= (1<<1) | (1<<0);     
    } else {                                    // could be the master node.
      if (thisCoord[iDim] != 0) 	{
	// if this node is not on the lower boundary
	brdcst_mode._bitval.brdcst_src |= (1<<2*iDim);
	max_mode._bitval.cmbn_src      |= (1<<2*iDim);
      }
      if (thisCoord[iDim] != maxCoord[iDim]) {
	// if this node is not on the upper boundary
	brdcst_mode._bitval.brdcst_src |= (1<<2*iDim+1);
	max_mode._bitval.cmbn_src      |= (1<<2*iDim+1);
      }
      if (thisCoord[iDim] != masterCoord[iDim]) { // not the master node
	buffer->not_master = 1;
	send_wire = recv_wire = (thisCoord[iDim] > masterCoord[iDim] ? 
				 2*iDim : 2*iDim+1);
	if ((thisCoord[iDim] == 0 || thisCoord[iDim] == maxCoord[iDim])
	    && iDim == 3) {                       // a leaf node
	  max_mode._intval = brdcst_mode._intval = 0;   
	}	
	if( iDim==0 ) {// different from the master node only in 0,1 direction
	  send_wire = 1;                         
	  max_mode._bitval.cmbn_src |= (1<<0) | (1<<1);
	}  
	break;                                    // !!!!
      }
    }
  }

  // fill in the log2 (rounded up!) of numOfNodes, which is nbits to be
  // reserved during global sum.
  // nbits = 1, 2, 3 if num_nodes = [1..2], [3..4], [5..8]
  {
    int nnodes = 1;
    for (iDim = 0; iDim < 4; ++iDim)
      nnodes *= size[iDim];
    int nbits = 1;
    int powerOf2 = 2;
    while (nnodes > powerOf2)  
      {++nbits; powerOf2 *= 2;}   
    buffer->nbits = nbits;  
  }
  
  // store these scu passthrough control register value
  buffer->max_mode = max_mode._intval;
  buffer->brdcst_mode = brdcst_mode._intval;
  buffer->send_wire = send_wire;
  buffer->recv_wire= recv_wire;

}

CPS_END_NAMESPACE
