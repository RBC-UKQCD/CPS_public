#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/verbose/main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/07/03 17:00:59  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/21 15:40:16  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:12:34  anj
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
//  $RCSfile: main.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/verbose/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <stdio.h>
#include <stdlib.h>	// exit()
#include<config.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/pmalloc.h>
#include<util/smalloc.h>
CPS_START_NAMESPACE

#ifdef PARALLEL
CPS_END_NAMESPACE
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif
 
Verbose VRB;
Error ERR;


int
main(int, char **)
{
  int size = 40;
  int *ptr;
  char *cname = "test";
  char *fname = "main";
  
  
  VRB.Level(20);
  printf("verbose level = %d\n",VRB.Level()  );
  

  VRB.Func(cname,fname);
  // Use when entering a function for flow control purposes.
  // It prints out   
  // class_name::func_name : Entered
  // If func_clock_active = 1 it also prints Clock = value
  // on the same line.


  VRB.Clock(cname,fname);
  // Use to print the clock value. It prints
  // class_name::func_name : Clock = value


  VRB.Clock(cname,fname,"This is the time at the beginning of main\n");
  // Use to print the clock value. It prints
  // class_name::func_name : Clock = value
  // The rest is as in printf.


  ptr = (int *) pmalloc(size * sizeof(int));
  if(ptr == 0)
    ERR.Pointer(cname,fname, "ptr");
  VRB.Pmalloc(cname,fname,
	      "ptr",ptr, size * sizeof(int));
  // Use when initializing a pointer with pmalloc. It prints out   
  // func_name : pmalloc initialized pointer 
  //    "ptr_name" to ptr with size size


  VRB.Pfree(cname,fname, "ptr", ptr);
  pfree(ptr);
  // Use before freeing a pointer with pfree. It prints out   
  // func_name : pfree will free pointer "ptr_name" = ptr 


  VRB.Pclear(cname,fname);
  // Use after calling pclear(). It prints out   
  // func_name : pclear() called


  ptr = (int *) smalloc(size * sizeof(int));
  if(ptr == 0)
    ERR.Pointer(cname,fname, "ptr");
  VRB.Smalloc(cname,fname,
	      "ptr",ptr, size * sizeof(int));
  // Use when initializing a pointer with smalloc. It prints out   
  // func_name : smalloc initialized pointer 
  //    "ptr_name" to ptr with size size


  VRB.Sfree(cname,fname, "ptr", ptr);
  sfree(ptr);
  // Use before freeing a pointer with sfree. It prints out   
  // func_name : sfree will free pointer "ptr_name" = ptr 


  VRB.Sclear(cname,fname);
  // Use after calling sclear(). It prints out   
  // func_name : sclear() called


  VRB.Flow(cname,fname, "This program is doing something that you must know.\n");
  // Use to print info in order follow the flow 
  // inside a function. Usage is as in printf.


  int par = 1997;
  VRB.Input(cname,fname, "The input parameters are = %d\n", par);
  // Use inside a function to print any relevant input 
  // values. Usage is as in printf.
        

  par = 1998;
  VRB.Result(cname,fname, "The result is = %d\n", par);
  // Use to print any results. Usage is as in printf.
        

  VRB.Warn(cname,fname, "Something may be wrong\n");
  // Use to print any Warnings. Usage is as 
  // in printf. It appends to the output WARNING.

  
  VRB.Debug(cname,fname, "Use some debugging\n");
  // Use to print debugging info. Usage is as 
  // in printf.

  VRB.Debug("Use some debugging just like printf.\n");
  // Use to print debugging info. Usage is as 
  // in printf.


  VRB.LedOn(cname,fname);
  // Use to turn on the LED and print 
  // class_name::func_name : LED on


  VRB.LedOff(cname,fname);
  // Use to turn off the LED and print 
  // class_name::func_name : LED off


  int number = 20;
  VRB.LedFlash(cname,fname, number);
  // Use to flash the LED and print 
  // class_name::func_name : LED flashing
  // "number" number of times.
  // It leaves the LED on.
  

  VRB.Clock(cname,fname);
  // Use to print the clock value. It prints
  // class_name::func_name : Clock = value


  VRB.Clock(cname,fname, "The program is almost done know and you may want to know the time\n");
  // Use to print the clock value. It prints
  // class_name::func_name : Clock = value
  // The rest is as in printf.


  VRB.FuncEnd(cname,fname);
  // Use when exiting a function for flow control purposes.
  // It prints out   
  // class_name::func_name : Exiting
  // If func_clock_active = 1 it also prints Clock = value
  // on the same line.

}


CPS_END_NAMESPACE
