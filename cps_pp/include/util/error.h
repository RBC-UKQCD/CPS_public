#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/error.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: error.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:13:17  anj
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
//  $RCSfile: error.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/error.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// error.h
//
// Header file for the Error base class. There are no derived 
// classes. The functions in this class print a message and exit 
// with a specific value. An object of this class called
// ERR should be created at the highest scope (outside main). 
// The header file declares ERR as external.
//
//------------------------------------------------------------------

#ifndef INCLUDED_ERROR_H
#define INCLUDED_ERROR_H


#define MAX_ERR_STR_LEN 500  
// The maximum string length to be used with vsprintf.
// This is needed by all error functions with ... argument.                              



//------------------------------------------------------------------
//
// Error is the base class
//
//------------------------------------------------------------------
class Error
{
  private:
    char *cname;    // Class name.
    char e_string[MAX_ERR_STR_LEN];  
       // Needed by all error functions with ... argument.                              

    int pointer_exit_value;
       // Pointer exits with exit(pointer_exit_value)

    int file_r_exit_value;
       // FileR exits with exit(file_r_exit_value)

    int file_w_exit_value;
       // FileW exits with exit(file_w_exit_value)

    int file_a_exit_value;
       // FileA exits with exit(file_a_exit_value)

    int not_implemented_exit_value;
       // NotImplemented exits with exit(not_implemented_exit_value)

    int hardware_exit_value;
       // Hardware exits with exit(hardware_exit_value)

    int general_exit_value;
       // General exits with exit(general_exit_value)

  public:
    Error();

    virtual ~Error();

    void Error::Pointer(char *class_name, char *func_name,
			char *ptr_name);   // pointer name
       // Error message for: Pointer not initialized.


     void Error::FileR(char *class_name, char *func_name, 
		       char *file_name);  // file name
       // Error message for: Can not open file to read.


     void Error::FileW(char *class_name, char *func_name, 
		       char *file_name);  // file name
       // Error message for: Can not open file to write.


     void Error::FileA(char *class_name, char *func_name, 
		       char *file_name);  // file name
       // Error message for: Can not open file to append.


     void Error::NotImplemented(char *class_name, char *func_name); 
       // Error message for: Not implemented.


     void Error::NotImplemented(char *class_name, char *func_name, 
                         const char *format,  // format of message
                         ...);                // arg. list of message  
       // Error message for: NotImplemented (write your own message).


     void Error::Hardware(char *class_name, char *func_name, 
			  const char *format,  // format of message
			  ...);                // arg. list of message  
       // Error message for: Hardware error (write your own message).
        

     void Error::General(char *class_name, char *func_name, 
                         const char *format,  // format of message
                         ...);                // arg. list of message  
       // Error message for: General error (write your own message).
        
};


//------------------------------------------------------------------
// External declarations.
//------------------------------------------------------------------
extern Error ERR;


#endif


CPS_END_NAMESPACE
