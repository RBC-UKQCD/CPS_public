#ifndef INCLUDED_ERROR_H
#define INCLUDED_ERROR_H  //!< Prevent multiple inclusions

#include<config.h>
#ifdef HAVE_QCDOCOS_SCU_CHECKSUM_H
#include <qcdocos/scu_checksum.h>
#endif
/*!\file
  \brief  Declaration and definition of Error class.

*/


CPS_START_NAMESPACE

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


    char *error_class_name;    
    char *error_func_name;	
    char *error_file_name;

    // Different types of error.
    enum{
	pointer,
        file_r,
	file_w,
	file_a,
	not_implemented,
	hardware,
        general,
	n_error_types
    };

    int exit_value[n_error_types];
    char* error_string[n_error_types];


  public:
    
    Error();

    virtual ~Error();

    void Pointer(const char*, const char*, const char*);  
    //!< Error message for an uninitialized pointer.

    void FileR(const char*, const char*, const char*);
    //!< Error message for failure to open a file to read.

    void FileW(const char*, const char*, const char*);  
    //!< Error message for failure to open a file to write.

    void FileA(const char*, const char*, const char*); 
    //!< Error message for failure to open a file to append.

    void NotImplemented(const char*, const char*);
    //!< Error message when something is not implemented.

    void NotImplemented(const char*, const char*, const char*, ...);
    //!< Error message when something is not implemented.

    void Hardware(const char*, const char*, const char*, ...);
    //!< Error message for hardware failure.

    void General(const char*, const char*, const char*, ...);
    //!< Error message for miscellaneous failure.    

    void HdwCheck(const char *,const char *);
        
};


//------------------------------------------------------------------
/*! An instance of the Error class, named ERR, should be
  created at the highest scope (outside main). This external declaration
  allows  access to the error message printing mechanisms.
*/
//------------------------------------------------------------------
extern Error ERR;


CPS_END_NAMESPACE
#endif
