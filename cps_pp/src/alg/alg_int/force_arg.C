#include <util/qcdio.h>
#include <alg/force_arg.h>

CPS_START_NAMESPACE

ForceArg::ForceArg() : cname("ForceArg"), L1(0.0), L2(0.0), Linf(0.0)
{
//  char *fname="ForceArg()";
//  printf("%s::%s: allocated (%p)\n",cname,fname,this);
}

ForceArg::ForceArg(Float L1, Float L2, Float Linf) : cname("ForceArg"), L1(fabs(L1)), L2(fabs(L2)), Linf(fabs(Linf))
{
//  char *fname="ForceArg()";
//  printf("%s::%s: allocated (%p)\n",cname,fname,this);
}

ForceArg::~ForceArg()
{
//  char *fname="~ForceArg()";
//  printf("%s::%s: deallocated (%p)\n",cname,fname,this);
}


void ForceArg::print(Float dt, char *label) {
  char *fname = "print(Float, char*)";
  FILE *fp;
  
#if TARGET==cpsMPI
  using MPISCU::fprintf;
#endif
  
  // Print out monitor info
  //---------------------------------------------------------------
  if( (fp = Fopen("force.dat", "a")) == NULL ) {
    ERR.FileA(cname,fname, "force.dat");
  }
  Fprintf(fp,"%s L1 = %e L2 = %e Linf = %e dt = %f\n", 
	  label, (IFloat)L1, (IFloat)L2, (IFloat)Linf, (IFloat)dt);
  Fclose(fp);
    
}

CPS_END_NAMESPACE
