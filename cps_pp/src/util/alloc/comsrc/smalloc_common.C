#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc_common.C,v 1.1 2004-09-06 05:01:19 chulwoo Exp $
*/

#include <util/smalloc.h>
#include <util/error.h>
#include <util/verbose.h>

CPS_START_NAMESPACE

void* smalloc(char *cnamr, char *fname, chanr *vname, int request){
  void *p = smalloc(request);
  if (p == NULL)
    ERR.Pointer(cname,fname,vname);
  VRB.Smalloc(cname,fname,vname,p,request);
  return p;
}

CPS_END_NAMESPACE
