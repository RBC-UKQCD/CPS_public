#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of Verbose class methods.

>>>>>>> 1.18.12.1
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-07-13 15:27:42 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/verbose/verbose.C,v 1.21 2012-07-13 15:27:42 chulwoo Exp $
//  $Id: verbose.C,v 1.21 2012-07-13 15:27:42 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: verbose.C,v $
//  $Revision: 1.21 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/verbose/verbose.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
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
const int MAX_STRING = 256;

//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
/*!
  The constructor sets the verbosity level to report results,
  function progress messages, clock output messages and RNG seed
  information by default.
*/
static void vrb_printf(const char *cname,const char *fname,const char *str){
#if TARGET == BGL || TARGET == BGP || (defined USE_QMP )
  if (!UniqueID()) printf("Node %d: %s::%s: %s",UniqueID(),cname,fname,str);
#else
  printf("%s::%s: %s",cname,fname,str);
#endif
}

Verbose VRB;
Verbose::Verbose(){
    
    cname = "Verbose";

// Default verbose level value.

    for(int i=0; i<N_VERBOSE_LEVELS; i++) active[i] = 0;
    ActivateLevel(VERBOSE_RESULT_LEVEL);
    ActivateLevel(VERBOSE_FLOW_LEVEL);
    ActivateLevel(VERBOSE_CLOCK_LEVEL);
    ActivateLevel(VERBOSE_RNGSEED_LEVEL);
   
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Verbose::~Verbose() {}


//------------------------------------------------------------------
/*!
  \return A value identifying the verbosity level.
*/
//------------------------------------------------------------------
//It is done this way for backwards compatibility.
int Verbose::Level() 
{
    int factor = 1,  level = 0;
    
    for(int i=0; i<N_VERBOSE_LEVELS; i++)
	if(active[i]){
	    level -= factor*i;
	    factor *= 100;
	}
    return level;
}


//------------------------------------------------------------------
/*!
  Depending on the level chosen, the categories of message are enabled.
  \param value The value to which the verbosity level is set.
 */
//------------------------------------------------------------------
// This is for backwards compatibilty
void Verbose::Level(int value){
    
    int i;
    if(value<0)
	for(i=0; i<N_VERBOSE_LEVELS; i++) active[i] = Active(i, value);
    else{
	if(value>=N_VERBOSE_LEVELS) value = N_VERBOSE_LEVELS-1;
	for(i=0; i<=value; i++) active[i] = 1;
	for(i=value; i<N_VERBOSE_LEVELS; i++) active[i] = 0;
    }
    
}


//------------------------------------------------------------------
/*!
  Depending on the level chosen, the categories of message are enabled.
  \param value The value to which the verbosity level is set.
 */
//------------------------------------------------------------------
void Verbose::Level(VerboseLevelType value){

    int i;
    for(i=0; i<=value; i++) active[i] = 1;
    for(i=value+1; i<N_VERBOSE_LEVELS; i++) active[i] = 0;
   
}

//------------------------------------------------------------------
/*!
  Based on the verbosity level, decides whether a particular category
  is enabled.

  \param level_value The message category.
  \return 1 if this message category should be enabled, 0 otherwise.
*/
//------------------------------------------------------------------

int Verbose::Active(int test_level, int level){

    if(test_level<level) return 1;
  
    if(level<0) {
	int x = -level;
	const int base = 100;
	int i;
	int y;
	for(i=0; i<base; i++){
	    y = x % base;
	    if( int(test_level) == y){
             if (!UniqueID())
             printf("Verbose::Acitve(%d,%d)=1\n",test_level,level);
             return 1;
             }
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
void Verbose::Func(const char *class_name, const char *func_name) {

    if(!active[VERBOSE_FUNC_LEVEL]) return;

//    printf("%s::%s : Entered :", class_name, func_name);
    vrb_printf(class_name,func_name,"Entered :");
    if(active[VERBOSE_FUNC_CLOCK_LEVEL]){
#ifdef _TARTAN
	printf("  Clock (12.5 MHz) = %d\n", (int)clock());
#else
	int cps = CLOCKS_PER_SEC;
	printf("  Clock (%2.1f MHz) = %d\n", cps/1.e+06, (int)clock());
#endif
    }
    else {
	if (!UniqueID()) printf("\n");
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
void Verbose::FuncEnd(const char *class_name, const char *func_name){
    
    if(!active[VERBOSE_FUNC_LEVEL]) return;

//    printf("%s::%s : Exiting :", class_name, func_name);
    vrb_printf(class_name,func_name,"Exiting :");
    if(active[VERBOSE_FUNC_CLOCK_LEVEL]){
#ifdef _TARTAN
	printf("  Clock (12.5 MHz) = %d\n", (int)clock());
#else
	int cps = CLOCKS_PER_SEC;
     if(!UniqueID())
	printf("  Clock (%2.1f MHz) = %d\n", cps/1.e+06, (int)clock());
#endif
    }
    else {
	if (!UniqueID()) printf("\n");
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
void Verbose::Pmalloc(const char *class_name, const char *func_name,
		      const char *ptr_name, const void *ptr, int size) {

    if(!active[VERBOSE_PMALLOC_LEVEL]) return;
    char vrb_string[MAX_STRING];
    
    sprintf(vrb_string,"pmalloc initialized pointer\t%s to %p with size 0x%x\n",
	   ptr_name, ptr, size);
    vrb_printf(class_name,func_name,vrb_string);

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
void Verbose::Pfree(const char *class_name, const char *func_name,
		    const char *ptr_name, const void *ptr){
    
    if(!active[VERBOSE_PMALLOC_LEVEL] ) return;
    char vrb_string[MAX_STRING];

    sprintf(vrb_string,"pfree will free pointer\t%s = %p\n",
	   ptr_name, ptr);
    vrb_printf(class_name,func_name,vrb_string);
    
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
void Verbose::Pclear(const char *class_name, const char *func_name) {

    if(!active[VERBOSE_PMALLOC_LEVEL]) return;

    printf("%s::%s : pclear called\n", class_name, func_name);
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
void Verbose::Smalloc(const char *class_name, const char *func_name,
		      const char *ptr_name, const void *ptr, int size){
    
    if(!active[VERBOSE_SMALLOC_LEVEL]) return;
    char vrb_string[MAX_STRING];

    sprintf(vrb_string,"smalloc initialized pointer\t%s to %p with size 0x%x\n",
	   ptr_name, ptr, size);
    vrb_printf(class_name,func_name,vrb_string);
    
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
void Verbose::Sfree(const char *class_name, const char *func_name,
		    const char *ptr_name, const void *ptr) {
    
    if(! active[VERBOSE_SMALLOC_LEVEL] ) return;
    char vrb_string[MAX_STRING];

    sprintf(vrb_string,"sfree will free pointer\t%s = %p\n",
	   ptr_name, ptr);
    vrb_printf(class_name,func_name,vrb_string);

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
void Verbose::Sclear(const char *class_name, const char *func_name) {
    
    if(! active[VERBOSE_SMALLOC_LEVEL] ) return;

    printf("%s::%s : sclear called\n", class_name, func_name);
    
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
void Verbose::Flow(const char *class_name, const char *func_name,
		   const char *format, ...){
    
    if(!active[VERBOSE_FLOW_LEVEL]) return;
    char vrb_string[MAX_STRING];
    
//    printf("%s::%s :", class_name, func_name);
    va_list args;
    va_start(args, format);
//    vprintf(format, args);
    vsnprintf(vrb_string, MAX_STRING,format, args);
    vrb_printf(class_name,func_name,vrb_string);
    if(active[VERBOSE_FLOW_CLOCK_LEVEL]){
#ifdef _TARTAN
	printf(" Clock (12.5 MHz) = %d\n\t", (int)clock());
#else
	int cps = CLOCKS_PER_SEC;
      if (!UniqueID())
	printf(" Clock (%2.1f MHz) = %d\n\t", cps/1.e+06, (int)clock());
#endif
    }
    else {
//	printf("\n\t");
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
void Verbose::Input(const char *class_name, const char *func_name, 
		    const char *format, ...) {
    
    if( !active[VERBOSE_INPUT_LEVEL] ) return;
    
    va_list args;
    va_start(args, format);
    printf("%s::%s :\n\t", class_name, func_name);
    vprintf(format, args);
    
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
void Verbose::Result(const char *class_name, const char *func_name,
		     const char *format, ...) {
    
    char vrb_string[MAX_STRING];
    if( !active[VERBOSE_RESULT_LEVEL] ) return;
    
    va_list args;
    va_start(args, format);
    vsnprintf(vrb_string, MAX_STRING,format, args);
    vrb_printf(class_name,func_name,vrb_string);
//    printf("%s::%s :\n\t", class_name, func_name);
    
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
void Verbose::Warn(const char *class_name, const char *func_name,
		   const char *format, ...) {
    
    if(!active[VERBOSE_WARN_LEVEL]) return;

    va_list args;
    va_start(args, format);
    printf("WARNING %s::%s :\n\t", class_name, func_name);
    vprintf(format, args);
    
    FILE *fp;
    char *filename = "phys.warn";
    if( (fp = Fopen(filename, "a")) == NULL ) {
	ERR.FileA("Verbose","Warn", filename);
    }
    Fprintf(fp,"WARNING %s::%s :\n\t", class_name, func_name);
    Vfprintf(fp, format, args);
    Fclose(fp);
    
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
void Verbose::Debug(const char *class_name, const char *func_name,
		    const char *format, ...) {
    
    if(!active[VERBOSE_DEBUG_LEVEL]) return;
    char vrb_string[MAX_STRING];
    
    va_list args;
    va_start(args, format);
    vsnprintf(vrb_string, MAX_STRING,format, args);
    vrb_printf(class_name,func_name,vrb_string);

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
void Verbose::Debug(const char *format, ...) {
    
    if(! active[VERBOSE_DEBUG_LEVEL]) return;
    
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    
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
void Verbose::LedOn(const char *class_name, const char *func_name) {
    
    if(!active[VERBOSE_LED_LEVEL]) return;
    
    printf("%s::%s : LED on\n", class_name, func_name);
#ifdef _TARTAN
    asm("   LDI     2h,     IOF");
#endif
    
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
void Verbose::LedOff(const char *class_name, const char *func_name) {
    
    if(!active[VERBOSE_LED_LEVEL]) return;

    printf("%s::%s : LED off\n", class_name, func_name);
#ifdef _TARTAN
    asm("   LDI     6h,     IOF");
#endif
    
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
void Verbose::LedFlash(const char *class_name, const char *func_name,
		       int number){
    
    if(!active[VERBOSE_LED_LEVEL]) return;

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
void Verbose::Clock(const char *class_name, const char *func_name,
		    const char *format, ...){
    
    if(!active[VERBOSE_CLOCK_LEVEL]) return;
    char vrb_string[MAX_STRING];

//    printf("%s::%s :", class_name, func_name);
    va_list args;
    va_start(args, format);
//    vprintf(format, args);
    vsnprintf(vrb_string, MAX_STRING,format, args);
    vrb_printf(class_name,func_name,vrb_string);
    int cps = CLOCKS_PER_SEC;
#if TARGET == BGL
  if(!UniqueID())
#endif
    printf(" Clock (%2.1f MHz) = %d\n\t", cps/1.e+06, (int)clock());
    
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
void Verbose::Clock(const char *class_name, const char *func_name){

    if(!active[VERBOSE_CLOCK_LEVEL]) return;

    printf("%s::%s :", class_name, func_name);
#ifdef _TARTAN
    printf(" Clock (12.5 MHz cycles) = %d\n", clock());
#else
    int cps = CLOCKS_PER_SEC;
    printf(" Clock (%2.1f MHz) = %d\n", cps/1.e+06, (int)clock());
#endif
    

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
void Verbose::RNGSeed(const char *class_name, const char *func_name,
		      const char *format, ...){
    
    if(!active[VERBOSE_RNGSEED_LEVEL]) return;

    printf("%s::%s :", class_name, func_name);
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    

}
        
/*!
  \param value The type of message to enable.
 */
void Verbose::ActivateLevel(VerboseLevelType value){

    active[value] = 1;

}

/*!
  \param value The type of message to enable.
 */
void Verbose::DeactivateLevel(VerboseLevelType value){

    active[value] = 0;

}

void Verbose::DeactivateAll(){

    for(int i=0; i<N_VERBOSE_LEVELS; i++) active[i] = 0;

}

/*!
  \param t The level of the type of message.
  \return 1 if messages of this type are currently enabled, zero otherwise.
*/
int Verbose::IsActivated(VerboseLevelType t) const{

    return active[t];
    
}
    



CPS_END_NAMESPACE
