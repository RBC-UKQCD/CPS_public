#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Routines implementing a closed doubly linked list.

  $Id $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/list.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: list.h,v 1.2 2003-07-24 16:53:53 zs Exp $
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
//  $RCSfile: list.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/list.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#ifndef _LINUX_LIST_H
#define _LINUX_LIST_H                //!< Prevent multiple inclusion.

/*
 * Simple doubly linked list implementation.
 *
 * Some of the internal functions ("__xxx") are useful when
 * manipulating whole lists rather than single entries, as
 * sometimes we already know the next/prev entries and we can
 * generate better code by using them directly rather than
 * using the generic single-entry routines.
 */
//! A list item.
//extern "C" {
struct list_head {
	struct list_head *next, *prev;
};
//}

//! Declare a minimal linked list.
/*! A single link, pointing to itself */
#define LIST_HEAD(name) \
	struct list_head name = { &name, &name }

//! Make a list into a minimal linked list.
/*! Makes the list a single link, pointing to itself */
#define INIT_LIST_HEAD(ptr) do { \
	(ptr)->next = (ptr); (ptr)->prev = (ptr); \
} while (0)


//! Insert a new entry between two known consecutive entries. 
/*!
 * This is only for internal list manipulation where we know
 * the prev/next entries already!
 */
inline void __list_add(struct list_head * new_item,
	struct list_head * prev,
	struct list_head * next )
{
	next->prev = new_item;
	new_item->next = next;
	new_item->prev = prev;
	prev->next = new_item;
}

//! Insert a new entry after the specified list item
/*!
  \param new_item The list item to add
  \param head The list item after which the new list item is to be placed
 */
inline void list_add(struct list_head *new_item, struct list_head *head)
{
	__list_add(new_item, head, head->next);
}


//!Delete a list entry 
/*!
   by making the prev/next entries point to each other.
 * This is only for internal list manipulation where we know
 * the prev/next entries already!
 */
inline void __list_del(struct list_head * prev,
				  struct list_head * next)
{
	next->prev = prev;
	prev->next = next;
}

//! Delete a list entry after the specified list item
/*!
  \param entry The list item to be deleted.
 */
inline void list_del(struct list_head *entry)
{
	__list_del(entry->prev, entry->next);
}

//! Inquire whether a list is minimal (single link)
/*!
  \param head The list item to test.
  \return True if \a head is the only link in the list, false otherwise.
 */
inline int list_empty(struct list_head *head)
{
	return head->next == head;
}


//! Insert one list into another.
/*!
  \param list The list to be inserted
  \param head The item after which \a list is to be inserted
 */
inline void list_splice(struct list_head *list, struct list_head *head)
{
	struct list_head *first = list->next;

	if (first != list) {
		struct list_head *last = list->prev;
		struct list_head *at = head->next;

		first->prev = head;
		head->next = first;

		last->next = at;
		at->prev = last;
	}
}

#define list_entry(ptr, type, member) \
	((type *)((char *)(ptr)-(unsigned long)(&((type *)0)->member)))


#endif

CPS_END_NAMESPACE
