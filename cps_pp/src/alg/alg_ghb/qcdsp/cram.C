#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:45 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/cram.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Id: cram.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:11:24  anj
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
//  Revision 1.2  2001/05/25 06:15:59  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: cram.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_ghb/qcdsp/cram.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/****************************************************************
*  These routines allow a user to push the cram		*
* vectors and cram status information onto a doubly linked	*
* list (stack) so that cram environments can be used and	*
* restored.							*
*   whipe_cram() just destroys the entire contents of the	*
*  cram.   							*
****************************************************************/
CPS_END_NAMESPACE
#include	<stdlib.h>
#include	<stdio.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE

struct Cram_contents;
struct Cram_contents {
  Cram_contents*		last;
  Cram_contents*		next;
  unsigned	cram_data[0x400];
  };

static Cram_contents*	cram_stack[2] 
  = { (Cram_contents*)0x000000, (Cram_contents*)0x000000 };

static unsigned* cram_base[2]
  = { (unsigned*)0x809800, (unsigned*)0x809c00 };

void  push_cram( int block )
{
  int i;
  void* malloc( unsigned );

  /* Check if this is the first call to push_cram0().  If it is,
  ** make sure that the stack pointer points to something! */
  if( cram_stack[block] == 0x000000 ) 
  {
    cram_stack[block] 
      = (Cram_contents*)malloc(sizeof(Cram_contents));
    cram_stack[block]->last = (Cram_contents*)0x000000;
  }

  /*  make space in memory for the next item. */
  cram_stack[block]->next 
    = (Cram_contents*)malloc(sizeof(Cram_contents));

  /*  Tell the next item where this item is. */
  cram_stack[block]->next->last = cram_stack[block];

  /*  Store all of the current information. */
  for(i=0; i<0x400; i++)
    cram_stack[block]->cram_data[i] = *(cram_base[block]+i);

  /* Set the stack pointer to the correct place. */
  cram_stack[block] = cram_stack[block]->next;

  return;
}

void  pop_cram( int block )
{
  int i;

  /* Check to be sure this is not the bottom of the stack!	*/
  if( cram_stack[block]->last == 0x000000 )
  { printf("Fatal CRAM error!\n");  exit(-1);	}

  /* Decrement the stack pointer. */
  cram_stack[block] = cram_stack[block]->last;

  for(i=0; i<0x400; i++)
    *(cram_base[block]+i) = cram_stack[block]->cram_data[i];

  /* Eliminate the unused entry. */
  free( cram_stack[block]->next );

  /* Set the pointer-ahead to null. */
  cram_stack[block]->next = (Cram_contents*) 0x000000;
}

void  whipe_cram( int block )
{
  int i;

  for(i=0; i<0x400; i++)
    *(cram_base[block] + i) = 0x00000000;

  return;

}
CPS_END_NAMESPACE
