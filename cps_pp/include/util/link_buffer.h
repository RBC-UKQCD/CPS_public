#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/link_buffer.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: link_buffer.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:13:18  anj
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
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: link_buffer.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/link_buffer.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#ifndef LINK_BUFFER_H
#define LINK_BUFFER_H

CPS_END_NAMESPACE
#include<util/list.h>
CPS_START_NAMESPACE

class Matrix;
class Lattice;
struct LinkEntry;

//-----------------------------------------------------------------
//
// LinkBuffer is a class that buffers the offnode links. 
// 
//-----------------------------------------------------------------

class LinkBuffer{
  char * cname;
  int site_stride[4];
  int node_stride[4];
  int site_size[4];
  int tab_size;
    //size of the hash table
  int buf_sz;
    //the size of the buffer, used for debuging 

  list_head * hash_tab; 
  //the hash table is indexed by the local coordiates and 
  //the direction of the link. in this way, it would be easy to remove
  //all links with the same local coordinates and orientaion. 
  //this is important because when a local link is changed, the same
  //links on other nodes are most likely also changed, if we have the 
  //buffered version of thoses links we should remove them from the buffer. 
  //Note we need to clean up these links when we try to update the
  //corresponding local link, no matter whether it's accepted or not.
  //---------------------------------------------------------------

  list_head flush_list;
  //flush_list defines the order in which the links should be removed,
  //if we run out of memory, that is, if free_list becomes empty.
  //---------------------------------------------------------------
  
  list_head free_list;
  //free_list is the pool of memory for buffering links.
  //---------------------------------------------------------------

  Lattice &lat;
  //reference to the lattice
  
  LinkEntry * link_entry_buf;

  int hash(const int *x, int mu);
  int GetNodeId(const int * x );
  // int check_lists();

 public: 
  LinkBuffer(Lattice & lat, int size); 
  ~LinkBuffer();

  //----------------------------------------------------------------
  // GetBufferedLink looks into the buffer to search for the a given 
  // link if it's not there, it gets the link from the neighbours 
  //----------------------------------------------------------------

  Matrix * GetBufferedLink(const int *x, int mu);
  
  //----------------------------------------------------------------
  // ClearBufferedLink looks into the buffer to search for a given link
  // if it's there, remove it from the buffer. 
  // In the heatbath algorithm when a local link is changed we need to make
  // sure that the offnode links in the buffer that has the same
  // local coordiates (in its local coordinates) get removed from 
  // the buffer  
  //----------------------------------------------------------------
  void ClearBufferedLink(const int *x, int mu);
  void ClearAll();

};

#endif

CPS_END_NAMESPACE
