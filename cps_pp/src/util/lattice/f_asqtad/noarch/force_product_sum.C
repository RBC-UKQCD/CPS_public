#include<config.h>

/*!\file
  \brief Helper routines for the fermionic force term

  $Id: force_product_sum.C,v 1.4 2004/08/19 22:20:49 mclark Exp $
*/
//--------------------------------------------------------------------

#include <stdio.h>
#include <util/lattice.h>
#include <util/gjp.h>

USING_NAMESPACE_CPS

extern "C" void Force_cross2dag(const Vector *chi, const Vector *phi, Matrix *result, int counter, double *fac);

// f += coeff * U_dir(x) * v(x)^dagger
// Note that STAG order stores hermitian conjugate of links.

void Fasqtad::force_product_sum(const Matrix *v,  int dir,
				IFloat coeff, Matrix *f){

  Matrix m, md;
  int n;

  int s[4];
  for(s[3]=0; s[3]<node_sites[3]; s[3]++)
    for(s[2]=0; s[2]<node_sites[2]; s[2]++)
      for(s[1]=0; s[1]<node_sites[1]; s[1]++) {
	for(s[0]=0; s[0]<node_sites[0]; s[0]++){
	  n = FsiteOffset(s);
	  m.DotMEqual(v[n], *(GaugeField()+GsiteOffset(s)+dir));
	  md.Dagger(m);
	  md *= coeff;
	  f[n] += md;
	}
      }

  int vol = node_sites[0]*node_sites[1]*node_sites[2]*node_sites[3];
  ForceFlops += vol * 234;

}


// f += u(vw)^dagger

void Fasqtad::force_product_sum(const Matrix *u, const Matrix *v,
				const Matrix *w, IFloat coeff, 
				Matrix *f, Matrix *tmp){

  Matrix m, md;
  int n;

  int s[4];
  for(s[3]=0; s[3]<node_sites[3]; s[3]++)
    for(s[2]=0; s[2]<node_sites[2]; s[2]++)
      for(s[1]=0; s[1]<node_sites[1]; s[1]++) {
	for(s[0]=0; s[0]<node_sites[0]; s[0]++){
	  n = FsiteOffset(s);
	  m.DotMEqual(v[n], w[n]);
	  md.Dagger(m);
	  md *= coeff;
	  f[n].DotMPlus(u[n], md);
	}    
      }

  int vol = node_sites[0]*node_sites[1]*node_sites[2]*node_sites[3];
  ForceFlops += vol * 432;

}


// f += coeff v w^dagger

void Fasqtad::force_product_sum(const Matrix *v, const Matrix *w,
				IFloat coeff, Matrix *f){

  Matrix m, md;
  int n;

  int s[4];
  for(s[3]=0; s[3]<node_sites[3]; s[3]++)
    for(s[2]=0; s[2]<node_sites[2]; s[2]++)
      for(s[1]=0; s[1]<node_sites[1]; s[1]++) {
	for(s[0]=0; s[0]<node_sites[0]; s[0]++){
	  n = FsiteOffset(s);
	  md.Dagger(w[n]);
	  m.DotMEqual(v[n], md);
	  m *= coeff;
	  f[n] += m;
	}
      }

  int vol = node_sites[0]*node_sites[1]*node_sites[2]*node_sites[3];
  ForceFlops += vol * 234;

}

// f += coeff (v w)^dagger

void Fasqtad::force_product_d_sum(const Matrix *v, const Matrix *w,
				    IFloat coeff, Matrix *f){

  Matrix m, md;
  int n;
  
  int s[4];
  for(s[3]=0; s[3]<node_sites[3]; s[3]++)
    for(s[2]=0; s[2]<node_sites[2]; s[2]++)
      for(s[1]=0; s[1]<node_sites[1]; s[1]++) {
	for(s[0]=0; s[0]<node_sites[0]; s[0]++) {
	  n = FsiteOffset(s);
	  m.DotMEqual(v[n], w[n]);
	  md.Dagger(m);
	  md *= coeff;
	  f[n] += md;
	}
      }

  int vol = node_sites[0]*node_sites[1]*node_sites[2]*node_sites[3];
  ForceFlops += vol * 234;

}



// The outer product of v and w, multiplied by the parity factor and
// coeffient, and added to the force f.
// N.B. The force term is multiplied by -1 on ODD parity sites. This the
// the MILC convention and is here to test against the MILC code.
// This should eventually be changed to EVEN for the CPS.

void Fasqtad::force_product_sum(const Vector *v, const Vector *w,
				    IFloat coeff, Matrix *f){

  int vol = GJP.XnodeSites() * GJP.YnodeSites() * GJP.ZnodeSites() * GJP.TnodeSites() ;
  ForceFlops +=78*vol;

  Matrix m;
  int n, s[4];
  
  //    printf("force_product_sum\n");
  for(s[3]=0; s[3]<GJP.TnodeSites(); s[3]++)
    for(s[2]=0; s[2]<GJP.ZnodeSites(); s[2]++)
      for(s[1]=0; s[1]<GJP.YnodeSites(); s[1]++)
	for(s[0]=0; s[0]<GJP.XnodeSites(); s[0]++){
	  n = FsiteOffset(s);
	  m.Cross2(v[n], w[n]);
	  // Note the extra factor of 2 from Matrix::Cross2; this is corrected by
	  // the factor of 1/2 in Matrix::TrLessAntiHermMatrix used in
	  // Fasqtad::update_momenta.
#if 0
	  if(parity(s)==CHKB_ODD) m *= -coeff;    
	  else m *= coeff;
#else
	  m *= coeff;
#endif
	  f[n] += m;
	  //		printf("%d %d %d %d f{%d]=%p\n",s[0],s[1],s[2],s[3],n,&(f[n]));
	}
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
