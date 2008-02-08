#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the CommonArg structure.

  $Id: common_arg.h,v 1.6 2008-02-08 22:01:05 chulwoo Exp $
*/
//--------------------------------------------------------------------------
#ifndef INCLUDED_COMMON_ARG_H
#define INCLUDED_COMMON_ARG_H  //!< Prevent multiple inclusion
CPS_END_NAMESPACE
#include <string.h>
#include <util/smalloc.h>
#include <util/error.h>
CPS_START_NAMESPACE

const int MAX_STRING_LEN = 20;


/*! \defgroup algargs Algorithm parameters
  \ingroup alg */

//! A container for parameters common to all algorithms.
/*! \ingroup algargs */
struct CommonArg {

    void *results;  /*!< This can point to a string containing the name
		      of a file to which to write results, or an object to
    		    which results might be passed.*/

    //! Optional initialisation with a filename and label.
    /*!\param l A label identifying this job.
      \param f The name of a file to which to write results.
    */
    CommonArg(const char *l=0, const char *f=0){
	filename = label = 0;
	results = 0;
	if(l) set_label(l);
	if(f) set_filename(f);
    }

    char *filename;   /*!< The name of a file to which to write results. */
    char *label;       /*!< A label which uniquely identifies this task. */ 

    //! Set the file name
    /*! Results produced by the algorithm are written to this file.
      \param s The name of the file
    */
    void set_filename(const char *s){

	if(filename) sfree(filename);
	if(! (filename = (char*)smalloc(strlen(s)+1)) )
	    ERR.Pointer("CommonArg", "set_filename", "filename");
	strcpy(filename, s);
	results = (void*)filename;

    }

    //! Assign the label name
    /*! A label can be assigned to a algorithm run.
      \param s The label
    */
    void set_label(const char *s){

	if(label) sfree(label); 
	if(! (label = (char*)smalloc(strlen(s)+1)) )
	    ERR.Pointer("CommonArg", "set_label", "label");
	strcpy(label, s);

    }

    ~CommonArg(){
    	if(filename) sfree(filename);
	if(label) sfree(label);
    }
};

#endif /* !INCLUDED_COMMON_ARG_H */

CPS_END_NAMESPACE
