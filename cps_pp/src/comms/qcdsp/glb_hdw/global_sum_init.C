#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:16 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_hdw/global_sum_init.C,v 1.2 2004-01-13 20:39:16 chulwoo Exp $
//  $Id: global_sum_init.C,v 1.2 2004-01-13 20:39:16 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: global_sum_init.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_hdw/global_sum_init.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//
//  global_sum_init.C
//       implementating the initialization part of global_sum.h
//
// 
//  See comments at the end on how the passthrough path is setup.
//
// 
//  Author:                     Ping Chen          April 1998
//
//
//  Last modified by:           Ping Chen          March 4, 1999
//
//
//

CPS_END_NAMESPACE
#include<global_sum.h>
#include<global_sum_info.h>
#include <stdlib.h>                   // void *malloc(unsigned size)
CPS_START_NAMESPACE
                                    
//--------------------------------------------------------------------------
// Constants used by global_sum.asm only.
//--------------------------------------------------------------------------
GlobalSumInfo *  GLOBAL_SUM_INFO_BUF = 0;


//--------------------------------------------------------------------------
// local inline functions
//--------------------------------------------------------------------------
inline unsigned int plusWire(unsigned int direction) 
{ return 2*direction+1; }

inline unsigned int minusWire(unsigned int direction) 
{ return 2*direction; }

inline unsigned int bitCode(unsigned int wire) 
{ return 1 << wire; }


//--------------------------------------------------------------------------
//  union Passthrough {	
//   int _intval;
//   struct     {
//    unsigned int brdcst_src  :8;            // the least significant bits 
//    unsigned int cmbn_src    :8;            
//    unsigned int final_carry :7;
//    unsigned int add_mode    :1;
//   } _bitval;
//  }
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
// void global_sum_init(int, int) 
//--------------------------------------------------------------------------
CPS_END_NAMESPACE
#include <comms/sysfunc.h>                   // SizeT() and etc.
CPS_START_NAMESPACE
void global_sum_init(int fast_mode, int max_num_try)
{
  int masterCoord[4];
  int thisCoord[4];
  int size[4];
  int scu_mask;
  
  // remap the physics coordinates into machine coordinates.
  //-----------------------------------------------------------------------
  {
    int physSize[4] =  {SizeT(), SizeX(), SizeY(), SizeZ()};
    int masterPhysCoord[4] = {0, 0, 0, 0};
    int thisPhysCoord[4] = { CoorT(), CoorX(), CoorY(), CoorZ()};
      
    for (int physDim = 0; physDim < 4; ++physDim) {
      int machDim = SCURemap((SCUDir)(2 * physDim)) / 2; 
      size[machDim]        = physSize[physDim];    
      masterCoord[machDim] = masterPhysCoord[physDim];
      thisCoord[machDim]   = thisPhysCoord[physDim];
    }
  }

  // calculate scu_mask
  //-----------------------------------------------------------------------
  if (fast_mode) {
    scu_mask = *(int *)0x813078;
    scu_mask |= 10;             // fast_mode and idle receive
  } else {
    scu_mask = 0xaaa02;
  }
  

  // call the general initialization routine
  //-----------------------------------------------------------------------
  global_sum_init(masterCoord, thisCoord, size, max_num_try, scu_mask);
}




//--------------------------------------------------------------------------
// void global_sum_init(...)
//--------------------------------------------------------------------------
void global_sum_init(const int masterMachineCoord[], 
		     const int thisMachineCoord[], 
		     const int size[],
		     int maxTRY,
		     int scuMask)
{	
  //--------------------------------------------------------------------------
  // set the values for max_num_try and scu_reset_mask
  //--------------------------------------------------------------------------
  if (!GLOBAL_SUM_INFO_BUF) 
    GLOBAL_SUM_INFO_BUF = (GlobalSumInfo *)malloc(sizeof(GlobalSumInfo));
  GlobalSumInfo *buffer = (GlobalSumInfo *)GLOBAL_SUM_INFO_BUF;
  buffer->max_num_try = (maxTRY > 1 ? maxTRY : 1);
  buffer->scu_reset_mask = scuMask;
    
  //------------------------------------------------------------------------
  // The following coordinates are in MACHINE order [0/1, 2/3, 4/5, 6/7].
  // Odd wires are the positive directions.
  //------------------------------------------------------------------------
  int maxCoord[4], thisCoord[4], masterCoord[4];
  int dim;

  for( dim = 0; dim < 4; ++dim)    {
    thisCoord[dim]  = thisMachineCoord[dim];
    masterCoord[dim] = masterMachineCoord[dim];
    maxCoord[dim]   = size[dim] - 1;    
  }

  //--------------------------------------------------------------------------
  // translate the coordinates by the amount of maxCoord[i]/2 - masterCoord[i]
  // so that the master is in the center of serial communication network.
  // Caution:
  //   The behaviour of operator % on a negative number is undefined.
  //--------------------------------------------------------------------------
  for( dim = 0; dim < 4; ++dim)    {
    thisCoord[dim] += maxCoord[dim] / 2 - masterCoord[dim];        
    thisCoord[dim]  = (thisCoord[dim] + size[dim]) % size[dim];  
    masterCoord[dim]  = maxCoord[dim] / 2;    
  }

  //------------------------------------------------------------------------
  // ordered_dim[] arranges the four dimensions according to the number of 
  // nodes in these dimensions from largest to smallest.
  // e.g., a single motherboard gives [0, 1, 2, 3].  Use bubble sort.
  //------------------------------------------------------------------------
  int ordered_dim[4];
  int non_local_dims = 0;
  int i;  
  for (i = 0; i < 4; ++i)    {
    ordered_dim[i] = i;
    if (size[i] > 1) 
      ++non_local_dims;    
  }
  for (i = 0; i < 4; ++i) {    
    for (int j = i+1; j < 4; ++j) {
      int idim = ordered_dim[i];
      int jdim = ordered_dim[j];
      if (size[idim] < size[jdim]) {
	ordered_dim[i] = jdim;
	ordered_dim[j] = idim;
      }
    }
  }
  
  // TASK:  if [0, 2, 1, 3], error: Wire 7 recv TIMEOUT   ???????
  // Seems: the shortest one cannot be 2
  ordered_dim[0] = 0;
  ordered_dim[1] = 2;
  ordered_dim[2] = 3;
  ordered_dim[3] = 1;


  //-------------------------------------------------------------------------
  // Figure out the wires selection to participate the passthrough mode.
  // 1. We choose the identical route for broadcast, max and sum to minimize
  //    the buffer storage space.
  // 2. Algorithm:
  // Compare coordinates of this node with those of the master node until
  // a difference occurs from wire 6, 7 to wire 0, 1.
  // This order is chosen to take advantage that
  // number of nodes in wire 0,1 is fixed to 4 and likely is the smallest.
  // If the difference occurs say in coord23, wires 7,6,5,4 plus either
  // 2 or 3 will beselected to participate in passthrough mode unless this 
  // node is on the boundary in these dimensions.
  //-------------------------------------------------------------------------
  {
    unsigned int all_wires = 0;  
    int special_wire; 

    for(int i = 0; i < non_local_dims; ++i )  {
      dim = ordered_dim[i];    
      // if this node is not on the lower boundary, or for the shortest dim
      if (thisCoord[dim] != 0 || i == non_local_dims-1) 	
	all_wires |= bitCode(minusWire(dim));  
      // if this node is not on the upper boundary, or for the shortest dim
      if (thisCoord[dim] != maxCoord[dim] || i == non_local_dims-1) 	
	all_wires |= bitCode(plusWire(dim));	
      // if this node is not the master node, break out of the loop after ...
      if (thisCoord[dim] != masterCoord[dim]) { 
	buffer->is_master = 0;
	if (thisCoord[dim] < masterCoord[dim] || i == non_local_dims-1)
	  special_wire = plusWire(dim);
	else
	  special_wire = minusWire(dim);
	// a leaf node (Definition: only recv/send during brdcst/sum_or_max)
	if (i == 0 && (thisCoord[dim]==0 || thisCoord[dim]==maxCoord[dim]))
	  all_wires = 0;   
	break;                                    // END of the loop
      }
      // if this node is the master node
      else if (i == non_local_dims - 1) {
	buffer->is_master = 1;
	special_wire = plusWire(dim);
	buffer->is_master = all_wires ^ bitCode(plusWire(dim));  
	all_wires ^= bitCode(minusWire(dim));  
      } 
    }

    // store the passthrough route into persistent buffer.
    buffer->all_wires = all_wires;
    buffer->special_wire = special_wire;
  }

  //--------------------------------------------------------------------------
  // fill in the log2 (rounded up!) of numOfNodes, which is nbits to be
  // reserved during global sum.
  // nbits = 1, 2, 3 if num_nodes = [1..2], [3..4], [5..8]
  //--------------------------------------------------------------------------
  {
    int nnodes = 1;
    for (dim = 0; dim < 4; ++dim)
      nnodes *= size[dim];
    int nbits = 1;
    int powerOf2 = 2;
    while (nnodes > powerOf2)  
      {++nbits; powerOf2 *= 2;}   
    buffer->nbits = nbits;  
  }
}



/////////////////////////////////////////////////////////////////////////////
// Comments:                        Ping Chen                     10/1/98
//
//

//              The path diagram for broadcast in 2-d
/////////////////////////////////////////////////////////////////////////////
// 
//                          dimension 1
//                             wire 3
//                                |
//                                |
// 
//                    Node       Node       Node        Node
//                     ^          ^          ^           ^
//                     ^          ^          ^           ^
//                    Node       Node       Node        Node
//                     ^          ^          ^           ^
//                     ^          ^          ^           ^
//    wire 0 --< <--  Node  < <  MASTER < < Node < < <  Node -<- wire 1    
// dimension 0         v          v          v           v       dimension 0
//                     v          v          v           v
//                    Node       Node       Node        Node
//               
//                                |
//                                |
//                             wire 2
//                           dimension 1
//
// The master node starts out by sending in regular scu mode. 
// All other nodes initiate their passthrough mode by waiting on their
// pre-determined receiving wire.
 

//
//                  The path diagram for max/sum in 2-d
/////////////////////////////////////////////////////////////////////////////
// 
//                          dimension 1
//                             wire 3
//                                |
//                                |
// 
//                    Node       Node       Node        Node
//                     v          v          v           v
//                     v          v          v           v
//                    Node       Node       Node        Node
//                     v          v          v           v
//                     v          v          v           v
//    wire 0 ------>> Node  > >  MASTER > > Node > > >  Node >> ---wire 1    
// dimension 0         ^          ^          ^           ^       dimension 0
//                     ^          ^          ^           ^
//                    Node       Node       Node        Node
//               
//                                |
//                                |
//                             wire 2
//                           dimension 1
// All nodes (including the master node) start out by initiate the 
// passthrough mode by sending in regular scu mode.
// The master node ends the whole process by waiting on the final
// result in regular scu receive mode.

// 
//                    How is the path determined?
/////////////////////////////////////////////////////////////////////////////
// 
// The master node is the node who broadcasts its information to every node
// when broadcast mode is selected and gets the final answer when max/sum
// mode is selected. 
//
// A node/db knows its role by comparing its coordinates with those of the 
// master node.
//
// Let master node have coordinates (M0, M1, M2, M3) and a general
// node have coordinates (X0, X1, X2, X3), the path is determined by 
// the following procedure.
// 
// Broadcast:
// 0>   if this node is the master node, i.e., if (X0, X1, X2, X3) =
//      (M0, M1, M2, M3), it will send to all directions on which it
//      is not the boundary node.
// 1>   if this node is on the line (X0 != M0, M1, M2, M3):
//      if it's on the upper side of the master node (i.e. X0 > M0) 
//      it receives from wire 0; 
//      otherwise, it receives from wire 1 otherwise.
//      In either cases, it will send to all other directions/wires 
//      (excluding, however, those on which it is the boundary node). 
//      
// 2>   if this node is on the surface (X0 != M0, X1 != M1, M2, M3):
//      if it's on the upper side of line (X0 != M0, M1, M2, M3) (i.e. X1>M1)
//      it receives from wire 2 and sends to wires 3,4,5,6,7
//      (excluding however, those on which it is the boundary node);
//      otherwise,  it receives from wire 3 and sends to wires 2,4,5,6,7
//      (excluding, however, those on which it is the boundary node). 
// 3>   if this node is on the hypersurface (X0 != M0, X1 != M1, X2 !=M2, M3):
//      if it's on the upper side of surface(X0 != M0, X1 !=M1, M2, M3) 
//      (i.e. if X2 > M2), it receives from wire 4 and sends to wires 5,6,7
//      (excluding however, those on which it is the boundary node);
//      otherwise,  it receives from wire 5 and sends to wires 4,6,7
//      (excluding, however, those on which it is the boundary node). 
// 4>   if this node is on the hyperspace (X0 != M0, X1 != M1, X2 !=M2,
//      X3 != M3):
//      if it's on the upper side of hypersurface(X0 != M0, X1 !=M1, X2 !=M2,
//      M3)(i.e. if X2 > M2), it receives from wire 6 and sends to wire 7
//      (if this node is the boundary node on direction 7, however, it
//      won't send, furthmore special care is needed since this node is
//      a leaf on this path, the scu passthrough contrl register at address
//      0x813068 has to be set to 0 instead of choosing broadcast mode).
//      Otherwise,  it receives from wire 7 and sends to wire 6
//      (if this node is the boundary node on direction 6, however, it
//      won't send, furthmore special care is needed since this node is
//      a leaf on this path, the scu passthrough contrl register at address
//      0x813068 has to be set to 0 instead of choosing broadcast mode).
// Note: 
//      During the process of determining the path
//      and storing the path information to some system information area,
//      we will offset the coordinates of every node such that the master node
//      will be at the center of coordinates so that it is easier to talk about
//      on which side of the master node a node is and it is easier to 
//      implement a balanced broadcast.


//
//                          Special cases - leaf nodes
/////////////////////////////////////////////////////////////////////////////
//
// Leaf node:
//     a node only sending during passthrough max/sum or only receiving 
//     during passthrough broadcast.
//     In this case, the passthrough control register has to be set to 0.
//

CPS_END_NAMESPACE
