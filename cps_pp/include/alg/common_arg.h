#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the CommonArg structure.

  $Id: common_arg.h,v 1.2 2003-07-24 16:53:53 zs Exp $
*/
//--------------------------------------------------------------------------
#ifndef INCLUDED_COMMON_ARG_H
#define INCLUDED_COMMON_ARG_H  //!< Prevent multiple inclusion
CPS_END_NAMESPACE
#include <string.h>
CPS_START_NAMESPACE

const int MAX_STRING_LEN = 20;  //!< Maximum length of strings

/*! \defgroup algargs Algorithm parameters
  \ingroup alg */

//! A container for parameters common to all algorithms.
/*! \ingroup algargs */
struct CommonArg {

    char label[MAX_STRING_LEN]; /*!< The label which uniquely 
				  identifies this task. */

    void *results;  /*!< The name of a file to which to write results. */

    // CTORS
    //! Optional initialisation with a filename.
    /*!\param l A label, less than MAX_STRING_LEN characters long,
      identifying this job.
      \param mem The name of a file to which to write results.
    */
  CommonArg(char* l = "", void* mem = 0) : results(mem)
  { strncpy(label, l, MAX_STRING_LEN); 
    label[MAX_STRING_LEN-1] = '\0'; }

  ~CommonArg() {}
};

#endif /* !INCLUDED_COMMON_ARG_H */

CPS_END_NAMESPACE
