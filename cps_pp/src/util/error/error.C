#include<config.h>
CPS_START_NAMESPACE
/*!\file 
  \brief   Definition of Error class methods.

  $Id: error.C,v 1.4 2004-01-13 20:39:49 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:49 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/error/error.C,v 1.4 2004-01-13 20:39:49 chulwoo Exp $
//  $Id: error.C,v 1.4 2004-01-13 20:39:49 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.3.2.1  2003/11/05 16:32:16  mike
//  Attempting to create new working branch!
//
//  Revision 1.3  2003/09/11 15:07:34  zs
//  Corrected documentation.
//
//  Revision 1.2  2003/07/24 16:53:54  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.2  2001/06/19 18:13:15  anj
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
//  Revision 1.2  2001/05/25 06:16:08  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: error.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/error/error.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>  //exit
#include <util/error.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Error messages
// Why are these here?
//------------------------------------------------------------------
char *pointer_str = "Error in %s::%s :\n\tpointer %s is not initialized.\n";
char *file_r_str = "Error in %s::%s :\n\tcan not open file %s to read.\n";
char *file_w_str = "Error in %s::%s :\n\tcan not open file %s to write.\n";
char *file_a_str = "Error in %s::%s :\n\tcan not open file %s to append.\n";
char *not_implemented_str = "Error in %s::%s :\n\tnot implemented.\n\t";
char *hardware_str = "Hardware error in %s::%s :\n\t";
char *general_str = "Error in %s::%s :\n\t";

//------------------------------------------------------------------
// Error names
//------------------------------------------------------------------
char *error_file_name = CWDPREFIX("phys.error");
char *error_class_name = "Error";

//------------------------------------------------------------------
// Constructor
/*!
   The exit codes are initialised here.
*/
//------------------------------------------------------------------
Error::Error() {

    cname = "Error";
    pointer_exit_value = -1;
    file_r_exit_value = -2;
    file_w_exit_value = -3;
    file_a_exit_value = -4;
    not_implemented_exit_value = -5;
    hardware_exit_value = -6;
    general_exit_value = -7;

}



//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
Error::~Error() {}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> : </tt>\n
  <tt>pointer <i>pointer name</i> is not initialized.</tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -1.

  \param class_name The name of a class whose method is being entered.
  \param func_name The name of a function or method.
  \param ptr_name The variable name of the pointer.
*/
//------------------------------------------------------------------
void Error::Pointer(char *class_name, char *func_name, 
		    char *ptr_name)
{
  FILE *fp;
  char *error_func_name = "Pointer";

  printf(pointer_str, class_name, func_name, ptr_name); 

  if( (fp = fopen(error_file_name, "w")) == NULL ) { 
    printf(file_w_str, error_class_name, error_func_name, error_file_name); 
    exit(file_w_exit_value); 
  } 
  fprintf(fp, pointer_str, class_name, func_name, ptr_name); 
  fclose(fp); 
  
  exit(pointer_exit_value); 
  
}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> : </tt>\n
  <tt>can not open file <i>file name</i> to read.</tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -2.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param file_name The name of the file.
*/
//------------------------------------------------------------------
void Error::FileR(char *class_name, char *func_name,
		  char *file_name)   // file name
{
  FILE *fp;
  char *error_func_name = "FileR";

  printf(file_r_str, class_name, func_name, file_name);

  if( (fp = fopen(error_file_name, "w")) == NULL ) {
    printf(file_w_str, error_class_name, error_func_name, error_file_name);
    exit(file_w_exit_value);
  }
  fprintf(fp, file_r_str, class_name, func_name, file_name);
  fclose(fp);

  exit(file_r_exit_value);
}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> : </tt>\n
  <tt>can not open file <i>file name</i> to write.</tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -3.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param file_name The name of the file.
*/
//------------------------------------------------------------------
void Error::FileW(char *class_name, char *func_name,
		  char *file_name)   // file name
{
  FILE *fp;
  char *error_func_name = "FileW";

  printf(file_w_str,class_name, func_name, file_name);

  if( (fp = fopen(error_file_name, "w")) == NULL ) {
    printf(file_w_str, error_class_name, error_func_name, error_file_name);
    exit(file_w_exit_value);
  }
  fprintf(fp, file_w_str,class_name, func_name, file_name);
  fclose(fp);

  exit(file_w_exit_value);
}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> :</tt>\n
  <tt>can not open file <i>file name</i> to append.</tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -4.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param file_name The name of the file.
*/
//------------------------------------------------------------------
void Error::FileA(char *class_name, char *func_name, 
		  char *file_name)   // file name
{
  FILE *fp;
  char *error_func_name = "FileA";

  printf(file_a_str, class_name, func_name, file_name);

  if( (fp = fopen(error_file_name, "w")) == NULL ) {
    printf(file_w_str, error_class_name, error_func_name, error_file_name);
    exit(file_w_exit_value);
  }
  fprintf(fp, file_a_str, class_name, func_name, file_name);
  fclose(fp);

  exit(file_a_exit_value);
}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> :</tt>\n
  <tt>not implemented.</tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -5.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
*/
//------------------------------------------------------------------
void Error::NotImplemented(char *class_name, char *func_name) 
{
  FILE *fp;
  char *error_func_name = "NotImplemented";

  printf(not_implemented_str, class_name, func_name);

  if( (fp = fopen(error_file_name, "w")) == NULL ) {
    printf(file_w_str, error_class_name, error_func_name, error_file_name);
    exit(file_w_exit_value);
  }
  fprintf(fp, not_implemented_str, class_name, func_name);
  fclose(fp);

  exit(not_implemented_exit_value);
}


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> :</tt>\n
  <tt>not implemented.</tt>\n
  <tt><i>message</i></tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -5.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.

  \todo Why not use vprintf rather than vsprintf followed by printf?      
*/
//------------------------------------------------------------------
void Error::NotImplemented(char *class_name, char *func_name,
			   const char *format,  // format of message
			   ...)                 // argument list of message
{
  FILE *fp;
  char *error_func_name = "NotImplemented";
  va_list args;
  va_start(args, format);

  printf(not_implemented_str, class_name, func_name);
  vsprintf(e_string, format, args);
  printf("%s",e_string);

  if( (fp = fopen(error_file_name, "w")) == NULL ) {
    printf(file_w_str, error_class_name, error_func_name, error_file_name);
    exit(file_w_exit_value);
  }
  fprintf(fp, not_implemented_str, class_name, func_name);
  fprintf(fp, "%s",e_string);
  fclose(fp);

  exit(not_implemented_exit_value);
}
        


//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Hardware error in <i>class name</i>::<i>function name</i> :</tt>\n
  <tt><i>message</i></tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -6.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
*/
//------------------------------------------------------------------
void Error::Hardware(char *class_name, char *func_name,
		    const char *format,  // format of message
		    ...)                 // argument list of message
{
  FILE *fp;
  char *error_func_name = "Hardware";
  va_list args;
  va_start(args, format);

  printf(hardware_str, class_name, func_name);
  vsprintf(e_string, format, args);
  printf("%s",e_string);

  if( (fp = fopen(error_file_name, "w")) == NULL ) {
    printf(file_w_str, error_class_name, error_func_name, error_file_name);
    exit(file_w_exit_value);
  }
  fprintf(fp, hardware_str, class_name, func_name);
  fprintf(fp, "%s",e_string);
  fclose(fp);

  exit(hardware_exit_value);
}

        

//------------------------------------------------------------------
/*!
  Prints the message
  \n<tt>Error in <i>class name</i>::<i>function name</i> :</tt>\n
  <tt><i>message</i></tt>\n
  to \c stdout and to the file \c phys.error.
  Exits with -7.

  \param class_name The name of a class.
  \param func_name The name of a function or method.
  \param format A format string for the message (<em>&aacute; la</em> \c printf)
  \param ... Optional arguments to the format string.
*/
//------------------------------------------------------------------
void Error::General(char *class_name, char *func_name,
		    const char *format,  // format of message
		    ...)                 // argument list of message
{
  FILE *fp;
  char *error_func_name = "General";
  va_list args;
  va_start(args, format);

  printf(general_str, class_name, func_name);
  vsprintf(e_string, format, args);
  printf("%s",e_string);

  if( (fp = fopen(error_file_name, "w")) == NULL ) {
    printf(file_w_str, error_class_name, error_func_name, error_file_name);
    exit(file_w_exit_value);
  }
  fprintf(fp, general_str, class_name, func_name);
  fprintf(fp, "%s",e_string);
  fclose(fp);

  exit(general_exit_value);
}
        



CPS_END_NAMESPACE
