#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/verbose.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: verbose.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.3  2001/09/06 11:51:47  anj
//  Minor modifications to the test suite, e.g. standardizing the
//  verbosity and such.  Collected the output from the original and the
//  latest versions using the new test suite, and checked them. Anj
//
//  Revision 1.2  2001/06/19 18:13:19  anj
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
//  $RCSfile: verbose.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/verbose.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// verbose.h
//
// Header file for the verbose class. An object of this class
// called VRB should be created at the highest scope (outside 
// main). The header file declares VRB as external.
//
//------------------------------------------------------------------

#ifndef INCLUDED_VERBOSE_H
#define INCLUDED_VERBOSE_H

// String length for vsprintf
//------------------------------------------------------------------
#define MAX_VRB_STR_LEN 500 
// The maximum string length to be used with vsprintf.
// This is needed by all verbose functions with ... argument.                              

// Verbose base
//------------------------------------------------------------------
#define VERBOSE_BASE 100
// When level is negative VRB_BASE is used
// as a base in order to extract the various function levels from 
// the single value of level. If a funcyion level is present
// the coresponding function level active flag is set.
// For example, if VERBOSE_BASE = 10 then
// if level = -142 then it returns 1 if level_value is
// equal to 1 or to 4 or to 2.


// Verbose function levels.
//------------------------------------------------------------------
enum VerboseLevelType {
  VERBOSE_WARN_LEVEL = 1,        // Controls Warn
  VERBOSE_RESULT_LEVEL = 2,      // Controls Result
  VERBOSE_PMALLOC_LEVEL = 3,     // Controls Pmalloc, Pfree, Pclear
  VERBOSE_SMALLOC_LEVEL = 3,     // Controls Smalloc, Sfree, Sclear
  VERBOSE_FUNC_LEVEL = 4,        // Controls Func, FuncEnd
  VERBOSE_FLOW_LEVEL = 5,        // Controls Flow
  VERBOSE_INPUT_LEVEL = 6,       // Controls Input
  VERBOSE_DEBUG_LEVEL = 7,       // Controls Debug
  VERBOSE_FUNC_CLOCK_LEVEL = 8,  // Clock output in Func, FuncEnd 
  VERBOSE_FLOW_CLOCK_LEVEL = 9,  // Clock output in Flow
  VERBOSE_CLOCK_LEVEL = 10,      // Controls Clock
  VERBOSE_LED_LEVEL = 11,        // Controls LedOn, LedOff
  VERBOSE_RNGSEED_LEVEL = 12     // Controls RNG seed output
};

//------------------------------------------------------------------
//
// Verbose is the base class
//
//------------------------------------------------------------------
class Verbose
{
  private:
    char *cname;    // Class name.
    char v_string[MAX_VRB_STR_LEN];  
       // Needed by all verbose functions with ... argument.                              

    int level;
       // The verbose level controls printing. All verbose function
       // level active flags are set by comparing level with
       // the corresponding verbose function level.
 
    int func_active;
       // Func, FuncEnd  active flag.

    int func_clock_active;
       // Clock output in Func, FuncEnd  active flag.

    int pmalloc_active;
       // Pmalloc, Pfree, Pclear active flag.

    int smalloc_active;
       // Smalloc, Sfree, Sclear active flag.

    int flow_active;
       // FuncFlow active flag.

    int flow_clock_active;
       // Clock output in Flow active flag.

    int input_active;
       // FuncInput active flag.

    int result_active;
       // Result active flag.

    int warn_active;
       // Warn active flag.

    int debug_active;
       // Debug active flag.

    int led_active;
       // Led active flag.

    int clock_active;
       // Clock active flag.

    int rngseed_active;
       // RNG seed flag.

  public:
    Verbose();

    virtual ~Verbose();


    int Level(void);
       // Returns the value of level.


    void Level(int);
       // Sets the value of level and activates the appropriate 
       // function level active flags.


    int Active(VerboseLevelType level_value);
       // It returns 1 if level_value < level.
       // If level is negative it rerturns 1 if
       // level_value matches an entry in level as expanded
       // in base VERBOSE_BASE (for example if VRB_BASE = 10 then
       // if level = -142 then it returns 1 if level_value is
       // equal to 1 or to 4 or to 2).
       // Otherwise it returns 0.


    void Func(char *class_name, char *func_name);
       // Use when enetring a function for flow control purposes.
       // It prints out   
       // class_name::func_name : Entered
       // If func_clock_active = 1 it also prints Clock = value
       // on the same line.


    void FuncEnd(char *class_name, char *func_name);
       // Use when exiting a function for flow control purposes.
       // It prints out   
       // class_name::func_name : Exiting
       // If func_clock_active = 1 it also prints Clock = value
       // on the same line.


    void Pmalloc(char *class_name, char *func_name,
			  char *ptr_name,    // pointer name
			  void *ptr,         // pointer address
			  int size);         // allocation size
       // Use when initializing a pointer with pmalloc. It prints out   
       // func_name : pmalloc initialized pointer 
       //    "ptr_name" to ptr with size size


    void Pfree(char *class_name, char *func_name,
			char *ptr_name,    // pointer name
			void *ptr);        // pointer address
       // Use before freeing a pointer with pfree. It prints out   
       // func_name : pfree wil free pointer "ptr_name" = ptr 


    void Pclear(char *class_name, char *func_name);
       // Use after calling pclear(). It prints out   
       // func_name : pclear() called


    void Smalloc(char *class_name, char *func_name, 
                          char *ptr_name,    // pointer name
			  void *ptr,         // pointer address
			  int size);         // allocation size
       // Use after initializing a pointer with smalloc. It prints out   
       // func_name : smalloc initialized pointer 
       //     "ptr_name" to ptr with size size


    void Sfree(char *class_name, char *func_name,
			char *ptr_name,    // pointer name
			void *ptr);        // pointer address
       // Use before freeing a pointer with sfree. It prints out   
       // func_name : sfree wil free pointer "ptr_name" = ptr 


    void Sclear(char *class_name, char *func_name);
       // Use after calling sclear(). It prints out   
       // func_name : sclear() called


   void Flow(char *class_name, char *func_name,
	     const char *format,  // format
	     ...);                // argument list   
        // Use to print info in order follow the flow 
        // inside a function. Usage is as in printf.


   void Input(char *class_name, char *func_name,
	      const char *format,  // format
	      ...);                // argument list   
        // Use inside a function to print any relevant input 
        // values. Usage is as in printf.
        

   void Result(char *class_name, char *func_name,
			 const char *format,  // format
			 ...);                // argument list   
        // Use to print any results. Usage is as in printf.
        

   void Warn(char *class_name, char *func_name,
			 const char *format,  // format
			 ...);                // argument list   
        // Use to print any Warnings. Usage is as 
	// in printf. It appends to the output WARNING.


   void Debug(char *class_name, char *func_name,
	                 const char *format,  // format
	                 ...);                // argument list   
        // Use to print debugging info. Usage is as 
	// in printf.

   void Debug(const char *format,  // format
	                 ...);     // argument list   
        // Use to print debugging info. Usage is as 
	// in printf.


   void LedOn(char *class_name, char *func_name);
        // Use to turn on the LED and print 
	// class_name::func_name : LED on


   void LedOff(char *class_name, char *func_name);
        // Use to turn off the LED and print 
	// class_name::func_name : LED off


   void LedFlash(char *class_name, char *func_name, int number);
       // Use to flash the LED and print 
       // class_name::func_name : LED flashing
       // "number" number of times.
       // It leaves the LED on.


   void Clock(char *class_name, char *func_name,
	                 const char *format,  // format
	                 ...);                // argument list   
	// Use to print the clock value. It prints
        // class_name::func_name : Clock = value
	// The rest is as in printf.


   void Clock(char *class_name, char *func_name);
	// Use to print the clock value. It prints
        // class_name::func_name : Clock = value


   void RNGSeed(char *class_name, char *func_name,
	     const char *format,  // format
	     ...);                // argument list   
        // Use to print the RNG seeds.
        // Usage is as in printf.

};

//------------------------------------------------------------------
// External declarations.
//------------------------------------------------------------------
extern Verbose VRB;


#endif


CPS_END_NAMESPACE
