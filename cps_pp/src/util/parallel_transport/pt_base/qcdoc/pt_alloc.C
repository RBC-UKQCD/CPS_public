#include <util/pt_int.h>
void *PT::Alloc(char *cname, char *fname, char *vname, int request,unsigned
int flag ){
    if (request<0){
      printf("Alloc(): %s::%s: %s %d bytes\n",cname,fname,vname,request);
      exit(-42);
    }
    if (request==0) return NULL;
    void *p= qalloc(flag,request);
    if (!p) PointerErr(cname,fname,vname);
    return p;
}
void *PT::FastAlloc(char *cname, char *fname, char *vname, int request ){
   if (request<0){
      printf("FastAlloc(): %s::%s: %s %d bytes\n",cname,fname,vname,request);
      exit(-42);
    }
    if (request==0) return NULL;
    void *p= FastAlloc(request);
    if (!p) PointerErr(cname,fname,vname);
    return p;
}
