#include <alg/dterms.h>
#if TARGET == QCDOC
#include <comms/sysfunc_cps.h>
#endif

CPS_START_NAMESPACE

DTerms::DTerms():term(NULL)
{
  cname = "DTerms" ;
}
  
DTerms::DTerms(int N):term(NULL), num(N)
{
  cname = "DTerms" ;
  //char *fname = "DTerms(int)";
  
  //VRB.Func(cname,fname);
  Allocate() ;
}

/*!
 Purpose:
   get the Term d at specified coordinates (s+vec).
    can deal with s+vec being off node.

 Arguments:
\li   d:       term index
\li   vec:     coordinates [x,y,z,t] of the Term we want to fetch. 
               These coordinates are relative to the [0,0,0,0]
               site of the node. They  could be out-of-range, i.e., 
               located off-node.
\li   tmp:     Buffer for the Complex number if communication is needed.
\li   return:  a reference to the Complex number. 
               If off-node, it points to tmp.


**/
Complex& DTerms::GetTerm(const int d, const int *vec, Complex& tmp) const
{
  // offset out-of-range coordinates site[] into on_node_site[]
  // in order to locate the Term
  //------------------------------------------------------------------------
  int on_node_site[4],site[4];
  int on_node = 1;
  Complex *on_node_term;
  {
    for (int i = 0; i < 4; ++i) {
      site[i] = on_node_site[i] = vec[i] ;
      while (on_node_site[i] < 0) {
        on_node_site[i] += GJP.NodeSites(i) ;
      }
      on_node_site[i] %= GJP.NodeSites(i);
      if (on_node_site[i] != site[i]) {
        on_node = 0;
      }
    }
    on_node_term = term + d + num*(on_node_site[0] + GJP.XnodeSites()*(
			           on_node_site[1] + GJP.YnodeSites()*(
                                   on_node_site[2] + GJP.ZnodeSites()*
                                   on_node_site[3]))) ;
  }

// Huey Wen's fix, TARTAN -> PARALLEL!!!
#ifndef PARALLEL
//VRB.FuncEnd(cname, fname) ;
  return *on_node_term;
#endif

  // send to the destination node if the site is off-node
  //------------------------------------------------------------------------
  if (on_node) {
//  VRB.FuncEnd(cname, fname) ;
    return *on_node_term;
  } else {
    Complex send = *on_node_term ;
    Complex &recv = tmp ;
    for (int i = 0; i < 4; ++i) {
      while (site[i] != on_node_site[i]) {
        if (site[i] < 0) {
	  // the Complex has 2 number of floats 
          getMinusData((IFloat *)&recv, (IFloat *)&send, sizeof(recv)/sizeof(IFloat), i);
          on_node_site[i] -= GJP.NodeSites(i) ;
        } else {
	  // the Complex has 2 number of floats 
          getPlusData((IFloat *)&recv, (IFloat *)&send, sizeof(recv)/sizeof(IFloat), i);
          on_node_site[i] += GJP.NodeSites(i) ;
        }
        send = recv;
      }
    }
//  VRB.FuncEnd(cname, fname) ;
    return recv ;
  }
}

void DTerms::Allocate()
{
  char *fname = "Allocate()";
  VRB.Func(cname, fname);
  term=(Complex*) smalloc(GJP.VolNodeSites()*sizeof(Complex)*num);
  if (term == 0) ERR.Pointer(cname, fname, "term");
  VRB.Smalloc(cname, fname, "term", term,
	      GJP.VolNodeSites()*sizeof(Complex)*num);
}
void DTerms::Delete()
{
  char *fname = "Delete()";
  VRB.Func(cname, fname);
  if(term != NULL)
    {
      VRB.Sfree(cname, fname, "term", term);
      sfree(term);
      term = NULL;
    }
}

CPS_END_NAMESPACE
