#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/error/error.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: error.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
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
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/error/error.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// error.C
//
// Error is the base class. There are no derived classes.
// The functions in this class print a message and exit 
// with a specific value. 
//
//------------------------------------------------------------------
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>  //exit
#include<util/error.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// Error messages
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
char *error_file_name = "phys.error";
char *error_class_name = "Error";

//------------------------------------------------------------------
// Constructor
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
// Pointer: Pointer not initialized.
//------------------------------------------------------------------
void Error::Pointer(char *class_name, char *func_name,
		    char *ptr_name)    // pointer name
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
// FileR: Can not open file to read.
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
// FileW: Can not open file to write.
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
// FileA: Can not open file to append.
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
// NotImplemented: Not implemented.
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
// NotImplemented: NotImplemented error (write your own message).
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
// Hardware: Hardware error (write your own message).
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
// General: General error (write your own message).
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
