#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration and definition of Error class.

  $Id: error.h,v 1.2 2003-07-24 16:53:53 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/error.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: error.h,v 1.2 2003-07-24 16:53:53 zs Exp $
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
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/error.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

#ifndef INCLUDED_ERROR_H
#define INCLUDED_ERROR_H  //!< Prevent multiple inclusions


#define MAX_ERR_STR_LEN 500
//!< The maximum length of error message strings

//! A class to handle the printing of error messages.
/*! Upon encountering an error from which the code cannot recover,
  one of these methods can be called to print a diagnostic message and exit
  with a predefined value.
  \todo Is it worth making things const correct in a fatal class?
*/
  
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

    void Pointer(char *class_name, char *func_name,
			char *ptr_name);  
    //!< Error message for an uninitialized pointer.

    void FileR(char *class_name, char *func_name,
		      char *file_name);
    //!< Error message for failure to open a file to read.

    void FileW(char *class_name, char *func_name,
		      char *file_name);  
    //!< Error message for failure to open a file to write.

    void FileA(char *class_name, char *func_name,
		      char *file_name); 
    //!< Error message for failure to open a file to append.

    void NotImplemented(char *class_name, char *func_name);
    //!< Error message when something is not implemented.

    void NotImplemented(char *class_name, char *func_name,
			       const char *format, ...);
    //!< Error message when something is not implemented.

    void Hardware(char *class_name, char *func_name,
			 const char *format, ...);
    //!< Error message for hardware failure.

    void General(char *class_name, char *func_name,
			const char *format, ...);
    //!< Error message for miscellaneous failure.    
        
};


//------------------------------------------------------------------
/*! An instance of the Error class, named ERR, should be
  created at the highest scope (outside main). This external declaration
  allows  access to the error message printing mechanisms.
*/
//------------------------------------------------------------------
extern Error ERR;


#endif



CPS_END_NAMESPACE
