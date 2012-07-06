#ifndef INCLUDED_VERBOSE_H
#define INCLUDED_VERBOSE_H       //!< Prevent multiple inclusion.

#include<config.h>
/*!\file
  \brief  Definition of Verbose class.

  $Id: verbose.h,v 1.5 2012-07-06 20:22:08 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-07-06 20:22:08 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/verbose.h,v 1.5 2012-07-06 20:22:08 chulwoo Exp $
//  $Id: verbose.h,v 1.5 2012-07-06 20:22:08 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: verbose.h,v $
//  $Revision: 1.5 $
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


CPS_START_NAMESPACE

//------------------------------------------------------------------
//! Definition of the types of message, related to various aspects of code activity.
/*! See the description of the Verbose class for more details. */

enum VerboseLevelType{
    VERBOSE_NONE_LEVEL,        //!< No messages.
    VERBOSE_WARN_LEVEL,        //!< Warning messages.
    VERBOSE_RESULT_LEVEL,      //!< Messages printing results
    VERBOSE_PMALLOC_LEVEL,     //!< Dynamic memory management messages
    VERBOSE_SMALLOC_LEVEL = VERBOSE_PMALLOC_LEVEL,     //!< Dynamic memory management messages
    VERBOSE_FUNC_LEVEL,        //!< Function entry and exit messages
    VERBOSE_FLOW_LEVEL,        //!< Function progress messages
    VERBOSE_INPUT_LEVEL,       //!< Function argument messages
    VERBOSE_DEBUG_LEVEL,       //!< Debugging messages
    VERBOSE_FUNC_CLOCK_LEVEL,  //!< Clock output on function entry and exit 
  VERBOSE_FLOW_CLOCK_LEVEL,  //!< Clock output during progress through function
    VERBOSE_CLOCK_LEVEL,      //!< Clock output
    VERBOSE_LED_LEVEL,        //!< LED
    VERBOSE_RNGSEED_LEVEL,     //!< RNG seeding information.
    N_VERBOSE_LEVELS	     //!< Number of message levels.
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
*/

class Verbose
{
  private:
    const char *cname;    // Class name.

    // Message type enabled flags
    int active[N_VERBOSE_LEVELS];

    int Active(int, int);
    // Decides whether a message category should be enabled.

  public:
    Verbose();

    virtual ~Verbose();

    int Level();
    //!< Gets the value of the verbosity level.

    void Level(int);
    //!< Sets the value of the verbosity level

    void Level(VerboseLevelType);
    //!< Sets the value of the verbosity level
    
    void Func(const char*, const char*);
    //!< For printing a message upon function entry.

    void FuncEnd(const char*, const char*);
    //!< For printing a message upon function exit.    

    void Pmalloc(const char*, const char*, const char*, const void*, int size);
    //!< For printing a message upon memory allocation.

    void Pfree(const char*, const char*, const char*, const void*);      
    //!< For printing a message when freeing allocated memory.
    
    void Pclear(const char*, const char*);
    //!< For printing a message after calling pclear(). 

    void Smalloc(const char*, const char*, const char*, const void*, int size);
    //!< For printing a message upon memory allocation.

    void Sfree(const char*, const char*, const char*, const void*);
    //!< For printing a message when freeing allocated memory.
    
    void Sclear(const char*, const char*);
    //!< For printing a message after calling pclear(). 

    void Flow(const char*, const char*, const char*, ...);                
    //!< For printing a message during function execution.

    void Input(const char*, const char*, const char*, ...);
   //!< For printing function arguments.

   void Result(const char*, const char*, const char*, ...);
   //!< For printing results.

   void Warn(const char*, const char*, const char*, ...);                
   //!< For printing warnings.

   void Debug(const char*, const char*, const char*, ...);
   //!< For printing debugging messages.

   void Debug(const char *format, ...);     
   //!< For printing debugging messages.

   void LedOn(const char*, const char*);
   //!< For printing LED messages.

   void LedOff(const char*, const char*);
   //!< For printing LED messages.   

   void LedFlash(const char*, const char*, int number);
   //!< For printing LED messages.      

   void Clock(const char*, const char*, const char*, ...);
   //!< For printing a message with the clock value.

   void Clock(const char*, const char*);
   //!< For printing the clock value.

   void RNGSeed(const char*, const char*, const char*, ...);
   //!< For printing RNG seeding information.

   //! Disable the printing of messages of a particular type.
   void DeactivateLevel(VerboseLevelType);

   //! Enable the printing of messages of a particular type.
   void ActivateLevel(VerboseLevelType);
   
   //! Disable the printing of all messages.
   void DeactivateAll();

   //! Test whether a type of message is enabled.
   int IsActivated(VerboseLevelType) const;
};

//------------------------------------------------------------------
/*! An instance of the Verbose class, named VRB, should be
  created at the highest scope (outside main). This external declaration
  allows control of and access to the message printing mechanisms.
*/
//------------------------------------------------------------------
extern Verbose VRB;


CPS_END_NAMESPACE
#endif
