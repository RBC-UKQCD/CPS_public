#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of Verbose class methods.

  $Id: verbose.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/verbose/verbose.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: verbose.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2001/09/06 11:51:48  anj
//  Minor modifications to the test suite, e.g. standardizing the
//  verbosity and such.  Collected the output from the original and the
//  latest versions using the new test suite, and checked them. Anj
//
//  Revision 1.4  2001/07/03 17:01:03  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.3  2001/06/29 12:04:17  anj
//
//  A few minor fixes and tests, but mostly a change in the I/O handling.
//  Off QCDSP, the I/O functions printf and fprintf are overriden by my
//  own qcdio.h library.  (This should eventually become part of the
//  general i/o spec.)  All this does is stop all processors from sending
//  out indentical output. Anj.
//
//  Revision 1.2  2001/06/19 18:13:41  anj
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
//  Revision 1.2  2001/05/25 06:16:11  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: verbose.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/verbose/verbose.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE

// GRF: Some compilers don't have this defined in their header files,
//      so here it is just in case.

#ifndef CLOCKS_PER_SEC
# define CLOCKS_PER_SEC 1000000
#endif

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
/*!
  The constructor sets the verbosity level to zero (no messages) by default.
*/
Verbose::Verbose() {
cname = "Verbose";

// Default verbose level value.
//------------------------------------------------------------------
level = 0;

// Default active levels.
//------------------------------------------------------------------
warn_active = 0;
result_active = 0;
pmalloc_active = 0;
smalloc_active = 0;
func_active = 0;
flow_active = 0;
input_active = 0;
rngseed_active = 0;

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Verbose::~Verbose() {}


//------------------------------------------------------------------
/*!
  \return The value of the verbosity level.
*/
//------------------------------------------------------------------
int Verbose::Level()
{
  return level;
}


//------------------------------------------------------------------
/*!
  Depending on the level chosen, the categories of message are enabled.
  \param value The value to which the verbosity level is set.
 */
//------------------------------------------------------------------
void Verbose::Level(int value)
{
  level = value;

  warn_active = Active(VERBOSE_WARN_LEVEL);
  result_active = Active(VERBOSE_RESULT_LEVEL); 
  pmalloc_active = Active(VERBOSE_PMALLOC_LEVEL); 
  smalloc_active = Active(VERBOSE_SMALLOC_LEVEL);
  func_active = Active(VERBOSE_FUNC_LEVEL);
  flow_active = Active(VERBOSE_FLOW_LEVEL);
  input_active = Active(VERBOSE_INPUT_LEVEL);
  debug_active = Active(VERBOSE_DEBUG_LEVEL);
  func_clock_active = Active(VERBOSE_FUNC_CLOCK_LEVEL);
  flow_clock_active = Active(VERBOSE_FLOW_CLOCK_LEVEL);
  clock_active = Active(VERBOSE_CLOCK_LEVEL);
  led_active = Active(VERBOSE_LED_LEVEL);
  rngseed_active = Active(VERBOSE_RNGSEED_LEVEL);

}



//------------------------------------------------------------------
/*!
  Based on the verbosity level, decides whether a particular category
  is enabled.

  Why is this method public?
  \param level_value The message category.
  \return 1 if this message category should be enabled, 0 otherwise.
*/
//------------------------------------------------------------------
int Verbose::Active(VerboseLevelType level_value)
{
  if( int(level_value) < level) return 1;
  
  if(level < 0) {
    int x = -level;
    int base = VERBOSE_BASE;
    int i;
    int y;
    for(i=0; i<base; i++){
      y = x % base;
      if( int(level_value) == y) return 1;
      x = (x - y) / base;
      if(x <= 0) break;
    }
  }

  return 0;
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : Entered</tt>\n
  is printed to \c stdout.  
  If the clock output on function entry and exit message category is enabled,
  the message
  \n<tt><class name>::<function name> : Entered  Clock (<clock frequency> MHz) = <clock></tt>\n
  is printed to \c stdout.
  \param class_name The name of the class whose method is being entered.
  \param func_name The name of the function/method being entered.
*/
//------------------------------------------------------------------
void Verbose::Func(char *class_name, char *func_name)
{
  if( func_active ){
    printf("%s::%s : Entered :", class_name, func_name);
    if(func_clock_active){
#ifdef _TARTAN
      printf("  Clock (12.5 MHz) = %d\n", clock());
#else
      int cps = CLOCKS_PER_SEC;
      printf("  Clock (%2.1f MHz) = %d\n", cps/1.e+06, clock());
#endif
    }
    else {
      printf("\n");
    }
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : Exiting</tt>\n
  is printed to \c stdout.  
  If the clock output on function entry and exit message category is enabled,
  the message
  \n<tt><class name>::<function name> : Exiting  Clock (<clock frequency> MHz) = <clock></tt>\n
  is printed to \c stdout.
  \param class_name The name of the class whose method is being entered.
  \param func_name The name of the function/method being entered.
*/
//------------------------------------------------------------------
void Verbose::FuncEnd(char *class_name, char *func_name)
{
  if( func_active ){
    printf("%s::%s : Exiting :", class_name, func_name);
    if(func_clock_active){
#ifdef _TARTAN
      printf("  Clock (12.5 MHz) = %d\n", clock());
#else
      int cps = CLOCKS_PER_SEC;
      printf("  Clock (%2.1f MHz) = %d\n", cps/1.e+06, clock());
#endif
    }
    else {
      printf("\n");
    }
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : pmalloc initialized pointer</tt>\n
  <tt>        <pointer name> to <pointer address> with size <size></tt>\n  
  is printed to \c stdout.  
  \param class_name The name of a class
  \param func_name The name of a function or method.
  \param ptr_name The variable name of the pointer
  \param ptr The pointer itself.
  \param size The amount of memory (in bytes) allocated.
*/
//------------------------------------------------------------------
void Verbose::Pmalloc(char *class_name, char *func_name,
		      char *ptr_name,    // pointer name
		      void *ptr,         // pointer address
		      int size)          // allocation size
{
  if( pmalloc_active ){
    printf("%s::%s : pmalloc initialized pointer\n\t%s to %x with size %x\n", 
	   class_name, func_name, ptr_name, int(ptr), size);
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : pfree will free pointer</tt>\n
  <tt>        <pointer name> = <pointer address></tt>\n  
  is printed to \c stdout.  
  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param ptr_name The variable name of the pointer.
  \param ptr The pointer itself.
 */
//------------------------------------------------------------------
void Verbose::Pfree(char *class_name, char *func_name,
		    char *ptr_name,    // pointer name
		    void *ptr)         // pointer address
{
  if( pmalloc_active ){
    printf("%s::%s : pfree will free pointer\n\t%s = %x\n", 
	   class_name, func_name, ptr_name, int(ptr));
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : pclear called</tt>\n
  is printed to \c stdout.  
  \param class_name The name of a class
  \param func_name The name of a function or method
 */
//------------------------------------------------------------------
void Verbose::Pclear(char *class_name, char *func_name) 
{
  if( pmalloc_active ){
    printf("%s::%s : pclear called\n", class_name, func_name);
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : smalloc initialized pointer</tt>\n
  <tt>        <pointer name> to <pointer address> with size <size></tt>\n  
  is printed to \c stdout.  
  \param class_name The name of a class
  \param func_name The name of a function or method.
  \param ptr_name The variable name of the pointer
  \param ptr The pointer itself.
  \param size The amount of memory (in bytes) allocated.
*/
//------------------------------------------------------------------
void Verbose::Smalloc(char *class_name, char *func_name,
		      char *ptr_name,    // pointer name
		      void *ptr,         // pointer address
		      int size)          // allocation size
{
  if( smalloc_active ){
    printf("%s::%s : smalloc initialized pointer\n\t%s to %x with size %x\n", 
	   class_name, func_name, ptr_name, int(ptr), size);
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : sfree will free pointer</tt>\n
  <tt>        <pointer name> = <pointer address></tt>\n  
  is printed to \c stdout.  
  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param ptr_name The variable name of the pointer.
  \param ptr The pointer itself.
 */
//------------------------------------------------------------------
void Verbose::Sfree(char *class_name, char *func_name,
		    char *ptr_name,    // pointer name
		    void *ptr)         // pointer address
{
  if( smalloc_active ){
    printf("%s::%s : sfree will free pointer\n\t%s = %x\n", 
	   class_name, func_name, ptr_name, int(ptr));
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : sclear called</tt>\n
  is printed to \c stdout.  
  \param class_name The name of a class
  \param func_name The name of a function or method
 */
//------------------------------------------------------------------
void Verbose::Sclear(char *class_name, char *func_name) 
{
  if( smalloc_active ){
    printf("%s::%s : sclear called\n", class_name, func_name);
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> :</tt>\n
  <tt><message></tt>\n  
  is printed to \c stdout.
  If the clock output on function entry and exit message category is enabled,
  the message
  \n<tt><class name>::<function name> : Clock (<clock frequency> MHz) = <clock></tt>\n
  <tt><message></tt>\n
  is printed to \c stdout.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
 */
//------------------------------------------------------------------
void Verbose::Flow(char *class_name, char *func_name,
		   const char *format,  // format
		   ...)                 // argument list   
{
  if( flow_active ){
    printf("%s::%s :", class_name, func_name);
    if(flow_clock_active){
#ifdef _TARTAN
      printf(" Clock (12.5 MHz) = %d\n\t", clock());
#else
      int cps = CLOCKS_PER_SEC;
      printf(" Clock (%2.1f MHz) = %d\n\t", cps/1.e+06, clock());
#endif
    }
    else {
      printf("\n\t");
    }
    va_list args;
    va_start(args, format);
    vsprintf(v_string, format, args);
    printf("%s",v_string);
  }
}
        

//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> :</tt>\n
  <tt><message></tt>\n
  is printed to \c stdout.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
 */
//------------------------------------------------------------------
void Verbose::Input(char *class_name, char *func_name, 
		    const char *format,  // format
		    ...)                 // argument list   
{
  if( input_active ){
    va_list args;
    va_start(args, format);
    printf("%s::%s :\n\t", class_name, func_name);
    vsprintf(v_string, format, args);
    printf("%s",v_string);
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> :</tt>\n
  <tt><message></tt>\n
  is printed to \c stdout.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
 */
//------------------------------------------------------------------
void Verbose::Result(char *class_name, char *func_name,
		      const char *format,  // format
		      ...)                 // argument list   
{
  if( result_active ){
    va_list args;
    va_start(args, format);
    printf("%s::%s :\n\t", class_name, func_name);
    vsprintf(v_string, format, args);
    printf("%s",v_string);
  }
}

        
//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt>WARNING <class name>::<function name> :</tt>\n
  <tt><message></tt>\n
  is printed to \c stdout and the file \c phys.warn.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
 */
//------------------------------------------------------------------
void Verbose::Warn(char *class_name, char *func_name,
		      const char *format,  // format
		      ...)                 // argument list   
{
  if( warn_active ){
    va_list args;
    va_start(args, format);
    printf("WARNING %s::%s :\n\t", class_name, func_name);
    vsprintf(v_string, format, args);
    printf("%s",v_string);
    
    FILE *fp;
    char *filename = "phys.warn";
    if( (fp = fopen(filename, "a")) == NULL ) {
      ERR.FileA("Verbose","Warn", filename);
    }
    fprintf(fp,"WARNING %s::%s :\n\t", class_name, func_name);
    fprintf(fp,"%s",v_string);
    fclose(fp);
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> :</tt>\n
  <tt><message></tt>\n
  is printed to \c stdout.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
*/
//------------------------------------------------------------------
void Verbose::Debug(char *class_name, char *func_name,
		    const char *format,  // format
		    ...)                 // argument list   
{
  if( debug_active ){
    va_list args;
    va_start(args, format);
    printf("%s::%s :\n\t", class_name, func_name);
    vsprintf(v_string, format, args);
    printf("%s",v_string);
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><message></tt>\n
  is printed to \c stdout.
  
  \param format A format string for the message (<em>&aacute; la</em> \c printf) 
  \param ... Optional arguments to the format string.
*/
//------------------------------------------------------------------
void Verbose::Debug(const char *format,  // format
		               ...)      // argument list   
{
  if( debug_active ){
    va_list args;
    va_start(args, format);
    vsprintf(v_string, format, args);
    printf("%s",v_string);
  }
}



//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : LED on</tt>\n
  is printed to \c stdout. The LED is turned on.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
*/
//------------------------------------------------------------------
void Verbose::LedOn(char *class_name, char *func_name)
{
  if( led_active ){
    printf("%s::%s : LED on\n", class_name, func_name);
#ifdef _TARTAN
    asm("   LDI     2h,     IOF");
#endif
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : LED off</tt>\n
  is printed to \c stdout. The LED is turned off.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
*/
//------------------------------------------------------------------
void Verbose::LedOff(char *class_name, char *func_name)
{
  if( led_active ){
    printf("%s::%s : LED off\n", class_name, func_name);
#ifdef _TARTAN
    asm("   LDI     6h,     IOF");
#endif
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : LED flashing</tt>\n
  is printed to \c stdout. The LED is flashed some number of times and left
  turned on.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
  \param number The number of time to flash.
*/

//------------------------------------------------------------------
void Verbose::LedFlash(char *class_name, char *func_name, int number)
{
  if( led_active ){
    printf("%s::%s : LED flashing\n", class_name, func_name);
#ifdef _TARTAN
    int repeat = 1500000;
    int dum = 0;
    int i, j;
    for(j = 0; j < number; j++){
      asm("   LDI     6h,     IOF");
      for(i = 0; i < repeat; i++){
	dum = dum + 1;
      }
      asm("   LDI     2h,     IOF");
      for(i = 0; i < repeat; i++){
	dum = dum + 1;
      }
    }
#endif
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : Clock (<clock frequency> MHz) = <clock></tt>\n
  <tt><message></tt>\n
  is printed to \c stdout.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
 */
//------------------------------------------------------------------
void Verbose::Clock(char *class_name, char *func_name,
		    const char *format,  // format
		    ...)                 // argument list   
{
  if( clock_active ){
    printf("%s::%s :", class_name, func_name);
#ifdef _TARTAN
    printf(" Clock (12.5 MHz cycles) = %d\n\t", clock());
#else
    int cps = CLOCKS_PER_SEC;
    printf(" Clock (%2.1f MHz) = %d\n\t", cps/1.e+06, clock());
#endif
    va_list args;
    va_start(args, format);
    vsprintf(v_string, format, args);
    printf("%s",v_string);
  }
}


//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> : Clock (<clock frequency> MHz) = <clock></tt>\n
  is printed to \c stdout.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
 */
//------------------------------------------------------------------
void Verbose::Clock(char *class_name, char *func_name)
{
  if( clock_active ){
    printf("%s::%s :", class_name, func_name);
#ifdef _TARTAN
    printf(" Clock (12.5 MHz cycles) = %d\n", clock());
#else
    int cps = CLOCKS_PER_SEC;
    printf(" Clock (%2.1f MHz) = %d\n", cps/1.e+06, clock());
#endif
  }
}



//------------------------------------------------------------------
/*!
  If the verbosity level enables this message category, the message
  \n<tt><class name>::<function name> :</tt>\n
  <tt><message></tt>\n
  is printed to \c stdout.
  
  \param class_name The name of a class
  \param func_name The name of a function or method
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
 */
//------------------------------------------------------------------
void Verbose::RNGSeed(char *class_name, char *func_name,
		   const char *format,  // format
		   ...)                 // argument list   
{
  if( rngseed_active ){
    printf("%s::%s :", class_name, func_name);
    va_list args;
    va_start(args, format);
    vsprintf(v_string, format, args);
    printf("%s",v_string);

  }
}
        






CPS_END_NAMESPACE
