#include<config.h>
CPS_START_NAMESPACE
/*  CIM Wed Jul  2 12:46:41 EDT 1997 */
#ifndef INCLUDED_COMMON_ARG_H
#define INCLUDED_COMMON_ARG_H
CPS_END_NAMESPACE
#include <string.h>
CPS_START_NAMESPACE

const int MAX_STRING_LEN = 20;


struct CommonArg {

  char label[MAX_STRING_LEN]; //  The label which uniquely
                              //  identifies this task.

  void *results;  // Pointer to location in memory where to
                  // store results.

  // CTORS
  CommonArg(char* l = "", void* mem = 0) : results(mem)
  { strncpy(label, l, MAX_STRING_LEN); 
    label[MAX_STRING_LEN-1] = '\0'; }

  ~CommonArg() {}
};

#endif /* !INCLUDED_COMMON_ARG_H */
CPS_END_NAMESPACE
