#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of Verbose class.

  $Id: verbose.h,v 1.2 2003-07-24 16:53:53 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/verbose.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: verbose.h,v 1.2 2003-07-24 16:53:53 zs Exp $
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
//  $Revision: 1.2 $
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
#define INCLUDED_VERBOSE_H       //!< Prevent multiple inclusion.

#define MAX_VRB_STR_LEN 500 
//!< The maximum length of message strings

#define VERBOSE_BASE 100
//!< Allows certain categories of message to be selected.
/*!< See the description of the Verbose class for more details. */ 


//------------------------------------------------------------------
//! Definition of the types of message, related to various aspects of code activity.
/*! See the description of the Verbose class for more details. */

enum VerboseLevelType {
  VERBOSE_WARN_LEVEL = 1,        //!< Warning messages.
  VERBOSE_RESULT_LEVEL = 2,      //!< Messages printing results
  VERBOSE_PMALLOC_LEVEL = 3,     //!< Dynamic memory management messages
  VERBOSE_SMALLOC_LEVEL = 3,     //!< Dynamic memory management messages
  VERBOSE_FUNC_LEVEL = 4,        //!< Function entry and exit messages
  VERBOSE_FLOW_LEVEL = 5,        //!< Function progress messages
  VERBOSE_INPUT_LEVEL = 6,       //!< Function argument messages
  VERBOSE_DEBUG_LEVEL = 7,       //!< Debugging messages
  VERBOSE_FUNC_CLOCK_LEVEL = 8,  //!< Clock output on function entry and exit 
  VERBOSE_FLOW_CLOCK_LEVEL = 9,  //!< Clock output during progress through function. 
  VERBOSE_CLOCK_LEVEL = 10,      //!< Clock output
  VERBOSE_LED_LEVEL = 11,        //!< LED
  VERBOSE_RNGSEED_LEVEL = 12     //!< RNG seeding information.
};

//! Class to control printing of informational messages.
/*!
  Various different types of message, relating to different types of code
  activity, are defined by the #VerboseLevelType. This class manages an
  overall verbosity level, which, depending on its numerical value, permits the
  printing of some of these types of informational messages.
  
  A message category is enabled if the verbosity level is larger than the
  number assigned to the category in the #VerboseLevelType enumeration type.
  The larger this verbosity level is, the more messages can be printed.

  If the verbosity level is negative, certain categories of message can be
  specifically enabled in a peculiar way determined by the value of
  #VERBOSE_BASE: Essentially, if #VERBOSE_BASE is defined to be 10, then
  the categories of message corresponding to the verbosity levels digits
  are enabled. Defining #VERBOSE_BASE to any other value appears to be less
  than useful, IMHO.

  \todo Most of these methods should be const with const arguments.
  \todo Why not use vprintf rather than vsprintf followed by printf?
*/

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
    //!< Gets the value of the verbosity level.

    void Level(int);
    //!< Sets the value of the verbosity level

    int Active(VerboseLevelType level_value);
    //!< Decides whether a message category should be enabled.

    void Func(char *class_name, char *func_name);
    //!< For printing a message upon function entry.

    void FuncEnd(char *class_name, char *func_name);
    //!< For printing a message upon function exit.    

    void Pmalloc(char *class_name, char *func_name,
		 char *ptr_name, void *ptr, int size);         
    //!< For printing a message upon memory allocation.

    void Pfree(char *class_name, char *func_name,
	       char *ptr_name, void *ptr);      
    //!< For printing a message when freeing allocated memory.
    
    void Pclear(char *class_name, char *func_name);
    //!< For printing a message after calling pclear(). 

    void Smalloc(char *class_name, char *func_name, 
		 char *ptr_name, void *ptr, int size);
    //!< For printing a message upon memory allocation.

    void Sfree(char *class_name, char *func_name,
	       char *ptr_name, void *ptr);
    //!< For printing a message when freeing allocated memory.
    
    void Sclear(char *class_name, char *func_name);
    //!< For printing a message after calling pclear(). 

    void Flow(char *class_name, char *func_name,
	      const char *format, ...);                
    //!< For printing a message during function execution.

    void Input(char *class_name, char *func_name,
	       const char *format, ...);
   //!< For printing function arguments.

   void Result(char *class_name, char *func_name,
	       const char *format, ...);
   //!< For printing results.
    

   void Warn(char *class_name, char *func_name,
	     const char *format, ...);                
   //!< For printing warnings.

   void Debug(char *class_name, char *func_name,
	      const char *format, ...);
   //!< For printing debugging messages.

   void Debug(const char *format, ...);     
   //!< For printing debugging messages.

   void LedOn(char *class_name, char *func_name);
   //!< For printing LED messages.

   void LedOff(char *class_name, char *func_name);
   //!< For printing LED messages.   

   void LedFlash(char *class_name, char *func_name, int number);
   //!< For printing LED messages.      

   void Clock(char *class_name, char *func_name,
	      const char *format, ...);
   //!< For printing a message with the clock value.

   void Clock(char *class_name, char *func_name);
   //!< For printing the clock value.

   void RNGSeed(char *class_name, char *func_name,
		const char *format, ...);
   //!< For printing RNG seeding information.

};

//------------------------------------------------------------------
/*! An instance of the Verbose class, named VRB, should be
  created at the highest scope (outside main). This external declaration
  allows control of and access to the message printing mechanisms.
*/
//------------------------------------------------------------------
extern Verbose VRB;


#endif



CPS_END_NAMESPACE
