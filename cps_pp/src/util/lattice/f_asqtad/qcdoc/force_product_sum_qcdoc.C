#include<config.h>
#include<qalloc.h>

/*!\file
  \brief Helper routines for the fermionic force term

  $Id: force_product_sum_qcdoc.C,v 1.2 2004-07-01 17:43:49 chulwoo Exp $
*/
//--------------------------------------------------------------------

#include <stdio.h>
#include <util/lattice.h>
#include <util/time.h>
#include <util/gjp.h>

USING_NAMESPACE_CPS

// The outer product of v and w, multiplied by the parity factor and
// coeffient, and added to the force f.
// N.B. The force term is multiplied by -1 on ODD parity sites. This the
// the MILC convention and is here to test against the MILC code.
// This should eventually be changed to EVEN for the CPS.

extern "C" void Force_cross2dag(const Vector *chi, const Vector *phi, Matrix *result, int counter, double *fac);

#undef PROFILE

void Fasqtad::force_product_sum(const Vector *v, const Vector *w,
				    IFloat coeff, Matrix *f){

  static int vol  = 0;
  char *fname = "force_product_sum(*V,*V,F,*M)";
	if (vol==0)
	vol = GJP.XnodeSites() * GJP.YnodeSites() * GJP.ZnodeSites() * GJP.TnodeSites() ;
  ForceFlops +=78*vol;
	unsigned long v2 = (unsigned long)v;
  if( qalloc_is_fast((Vector *)v) &&
  qalloc_is_fast((Vector *)w) &&
   qalloc_is_fast((Matrix *)f) )
	v2 = v2 - 0xb0000000 + 0x9c000000;
	
#ifdef PROFILE
  Float dtime = -dclock();
#endif
  IFloat coeff2 = 2.0*coeff;
  Force_cross2dag((Vector *)v2, w, f, vol/2, &coeff2);
#ifdef PROFILE
	dtime += dclock();
	print_flops(cname,fname,78*vol,dtime);
#endif
#if 0
  {IFloat *tmp = (IFloat *)f;
   printf("result[0]=");
   for(int i = 0;i<18;i++){
	printf(" %0.8e",*(tmp++));
	if(i%6==5) printf("\n");
   }
   tmp = (IFloat *)&(f[vol-1]);
   printf("result[%d]=",vol-1);
   for(int i = 0;i<18;i++){
	printf(" %0.8e",*(tmp++));
	if(i%6==5) printf("\n");
   }
  }
  exit(54);
#endif
}
