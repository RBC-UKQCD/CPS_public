#include<config.h>
/*!\file
  \brief  Implementation of dynamic memory management routines.	

  $Id: smalloc_common.C,v 1.3 2004-09-28 16:10:28 chulwoo Exp $
*/

#include <util/smalloc.h>
#include <util/error.h>
#include <util/verbose.h>

CPS_START_NAMESPACE

void* smalloc(char *cname, char *fname, char *vname, int request){
  if (request < 1)
    ERR.General(cname,fname,"requested size of %s is %d\n",vname,request);
  void *p = smalloc(request);
  if (p == 0)
    ERR.Pointer(cname,fname,vname);
  VRB.Smalloc(cname,fname,vname,p,request);
  return p;
}
void sfree(char *cname, char *fname, char *vname, void *p){
  VRB.Sfree(cname,fname,vname,p);
  if (p == 0)
    ERR.Pointer(cname,fname,vname);
  sfree(p);
}

void* fmalloc(char *cname, char *fname, char *vname, int request){
  void *p = fmalloc(request);
  if (p == 0)
    ERR.Pointer(cname,fname,vname);
  VRB.Smalloc(cname,fname,vname,p,request);
  return p;
}
void ffree(char *cname, char *fname, char *vname, void *p){
  VRB.Sfree(cname,fname,vname,p);
  if (p == 0)
    ERR.Pointer(cname,fname,vname);
  sfree(p);
}

CPS_END_NAMESPACE
