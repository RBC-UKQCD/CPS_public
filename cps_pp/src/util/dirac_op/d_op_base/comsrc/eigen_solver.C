#include <config.h>
CPS_START_NAMESPACE
CPS_END_NAMESPACE
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <util/vector.h>
#include <util/dirac_op.h>
#include<math.h>
#include <util/smalloc.h>
#include <util/vector_template_func.h>
CPS_START_NAMESPACE

Float pythag(Float a, Float b);
void tred2(Float **a, int n, Float *d, Float *e);
int tqli(Float *d, Float *e, int n, Float **z);
int tqli(Float *d, Float *e, int n, Float **z,int iters);
int tqli(Float *d, Float *e, int n, Float **z, int numshifts, Float* shifts);
int tqri(Float *d, Float *e, int n, Float **z,int iters);
void eigsrt(Float *d, int n);
void eigsrt(Float *d, Float **v, int n);
void eigsrt(Float *d, Float **v1, int n, Float **v2, int n2);
void eigsrt(Float *d, Float *e, Float **v1, int n, Float **v2, int n2);
void eigsrt(Float *d, Float *e, Float *tte, Float **v1, int n, Float **v2, int n2);
void eigsrt_low(Float *d, Float *e, Float *tte, Float **v1, int n, Float **v2, int n2);
void eigsrt_low(Float *d, Float **v1, int n, Float **v2, int n2);
IFloat vecNorm(const Float *a, int len);
int tqri(Float *d, Float *e, Float *s, int n, Float **Q, int steps);
int tqri(Float *d, Float *e, int nk, int n, Float **z);
int tqri2(Float *d, Float *e, Float *s, int n, Float **Q, int steps);
int tqri3(Float *d, Float *e, Float *s, int n, Float **Q, int steps);
void mult_vec_by_tridiag(Float* out, Float *in, Float *d, Float *e, int n);
void eigsrt_low(Float *d, int n);
void lanczos_GramSchm(Float *psi, Float **vec, int Nvec, int f_size, Float* alpha);
void QRtrf(Float *d, Float *e,int nk, int n, Float **Q, Float dsh, int kmin, int kmax);
void mvfloattoFloat(Float* out, float* in, int f_size);
void mvFloattofloat(float* out, Float* in, int f_size);

//int eigen_solver(Float **A, int n, Float *Eval, Float *aux)
int eigen_solver(Float *A, Float *Evec, Float *Eval, int n)
//use the same call as the old version of eigen solver, so do not need to change the function that call it
//Float A[n*(n+1)/2], EV[n*n], E[n];
{
	int i,j;
	Float **a=(Float **)smalloc(n*sizeof(Float*));
	Float *aux=(Float *)smalloc(n*sizeof(Float));
	for(i=0;i<n;i++)
		a[i]=&Evec[i*n];
	for(i=0;i<n;i++)
		for(j=0;j<=i;j++)
			a[j][i]=a[i][j]=A[i*(i+1)/2+j];
	tred2(a,n,Eval,aux);
	int success = tqli(Eval,aux,n,a);
	sfree(a);
	sfree(aux);
	return success;
}
void tred2(Float **a, int n, Float *d, Float *e)
{
    int     l, k, j, i;
    Float  scale, hh, h, g, f;

    for (i = n-1; i > 0; i--)
    {
        l = i - 1;
        h = scale = 0.0;
        if (l > 0)
        {
	  for (k = 0; k < i; k++)
                scale += fabs(a[i][k]);
            if (scale == 0.0)
                e[i] = a[i][l];
            else
            {
                for (k = 0; k < i; k++)
                {
                    a[i][k] /= scale;
                    h += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = f > 0 ? -sqrt(h) : sqrt(h);
                e[i] = scale * g;
                h -= f * g;
                a[i][l] = f - g;
                f = 0.0;
                for (j = 0; j < i; j++)
                {
                    a[j][i] = a[i][j] / h;
                    g = 0.0;
                    for (k = 0; k < j+1; k++)
                        g += a[j][k] * a[i][k];
                    for (k = j + 1; k < i; k++)
                        g += a[k][j] * a[i][k];
                    e[j] = g / h;
                    f += e[j] * a[i][j];
                }
                hh = f / (h + h);
                for (j = 0; j < i; j++)
                {
                    f = a[i][j];
                    e[j] = g = e[j] - hh * f;
                    for (k = 0; k < j+1; k++)
                        a[j][k] -= (f * e[k] + g * a[i][k]);
                }
            }
        } else
            e[i] = a[i][l];
        d[i] = h;
    }
    d[0] = 0.0;
    e[0] = 0.0;
    for (i = 0; i < n; i++)
    {
        if (d[i]!=0.0)
        {
            for (j = 0; j < i; j++)
            {
                g = 0.0;
                for (k = 0; k < i; k++)
                    g += a[i][k] * a[k][j];
                for (k = 0; k < i; k++)
                    a[k][j] -= g * a[k][i];
            }
        }
        d[i] = a[i][i];
        a[i][i] = 1.0;
        for (j = 0; j < i; j++)
            a[j][i] = a[i][j] = 0.0;
    }
}

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
int tqli(Float *d, Float *e, int n, Float **z)
{
    int     m, l, iter, i, k;
    Float  s, r, p, g, f, dd, c, b;

    for (i = 1; i < n; i++) e[i - 1] = e[i];
    e[n-1] = 0.0;
    for (l = 0; l < n; l++)
    {
        iter = 0;
        do
        {
            for (m = l; m < n - 1; m++)
            {
                dd = fabs(d[m]) + fabs(d[m + 1]);
                if (fabs(e[m]) + dd == dd)
                    break;
            }
            if (m != l)
            {
                if (iter++ == 30)
                    return (0);
                g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                r = sqrt((g * g) + 1.0);
                g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                s = c = 1.0;
                p = 0.0;
                for (i = m - 1; i >= l; i--)
                {
                    f = s * e[i];
                    b = c * e[i];
                    if (fabs(f) >= fabs(g))
                    {
                        c = g / f;
                        r = sqrt((c * c) + 1.0);
                        e[i + 1] = f * r;
                        c *= (s = 1.0 / r);
                    } else
                    {
                        s = f / g;
                        r = sqrt((s * s) + 1.0);
                        e[i + 1] = g * r;
                        s *= (c = 1.0 / r);
                    }
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;
                    //for (k = 0; k < n; k++)
                    //{
		    //  f = z[k][i + 1];
		    //  z[k][i + 1] = s * z[k][i] + c * f;
		    //  z[k][i] = c * z[k][i] - s * f;
                    //}
                }
                d[l] = d[l] - p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }

    return (1);
}

Float pythag(Float a, Float b){
  Float absa=fabs(a);
  Float absb=fabs(b);
  return( absa > absb ? absa*sqrt(1.0+pow(absb/absa,2)) :
	  (absb == 0.0 ? 0.0 : absb*sqrt(1.0+pow(absa/absb,2))));
}

int tqri(Float *d, Float *e, int nk, int n, Float **z)
{
  int j;
  int kmin = 0;
  int kmax = nk-1;
  int iter=0;
  int maxiter=n*100;

  while(iter < maxiter){
    
    // diagonalize leading submatrix to compute shift
    Float dsub = d[kmax]-d[kmax-1];
    //Float dd = sqrt(dsub*dsub+4.0*e[kmax-1]*e[kmax-1]);
    Float dd = pythag(dsub,2.0*e[kmax-1]);
    Float dsh = 0.5*(d[kmax-1]+d[kmax]+SIGN(dd,dsub));

    // do the rotation and chase bulge
    QRtrf(d,e,nk,n,z,dsh,kmin,kmax);

    for(j=kmax;j>=kmin+1;j--){
      dd = fabs(d[j-1]) + fabs(d[j]);
      if (fabs(e[j-1]) + dd > dd){
	kmax=j;
	break;
      }
    }

    if(j==kmin)return(iter);

    for(j=0;j<kmax;j++){
      dd = fabs(d[j]) + fabs(d[j + 1]);
      if (fabs(e[j]) + dd > dd){
	kmin=j;
	break;
      }
    }
    iter++;
  }
  
  return (iter);
}

void QRtrf(Float *d, Float *e,int nk, int n, Float **Q, Float dsh, int kmin, int kmax){


  // printf("enetring QRTrf shift=%.16f\n",dsh);

  // for(int i=0;i<n;++i)
  // printf("alpha beta %d %.16f %.16f\n", i, d[i], e[i]);
  // printf("\n");
  
  int k=kmin;
  //Float den = 1.0/sqrt((d[k]-dsh)*(d[k]-dsh)+e[k]*e[k]);
  Float den = 1.0/pythag(d[k]-dsh,e[k]);
  Float c = (d[k]-dsh)*den;
  Float s = -e[k]*den;
  //printf("%d D,e,dsh,c,s =%.16f %.16f %.16f, %.16f %.16f\n",k+1, d[k]-dsh,e[k],dsh,c,s);
  Float tmp1 = d[k];
  Float tmp2 = d[k+1];
  Float tmpb = e[k];
  d[k] = c*c*tmp1 + s*s*tmp2 - 2.0*c*s*tmpb;
  d[k+1] = s*s*tmp1 + c*c*tmp2 + 2.0*c*s*tmpb;
  e[k] = c*s*(tmp1-tmp2) + (c*c-s*s)*tmpb;
  Float x = -s*e[k+1];
  e[k+1] *= c;
  for(int i=0;i<nk;i++){
    Float Qtmp1 = Q[i][k];
    Float Qtmp2 = Q[i][k+1];
    Q[i][k] = c*Qtmp1 - s*Qtmp2;
    Q[i][k+1] = s*Qtmp1 + c*Qtmp2; 
  }
  // givens rotations
  for(int k = kmin+1; k < kmax; k++){
    //den = 1.0/sqrt(x*x + e[k-1]*e[k-1]);
    den = 1.0/pythag(x,e[k-1]);
    c = e[k-1]*den;
    s = - x*den;
    //printf("%d d,e,dsh,c,s =%.16f  %.16f %.16f, %.16f %.16f\n",k+1,x,e[k-1],dsh,c,s);
    //printf("d[nk-2] %.16f\n",d[nk-2]);
    tmp1 = d[k];
    tmp2 = d[k+1];
    tmpb = e[k];
    d[k] = c*c*tmp1 + s*s*tmp2 - 2.0*c*s*tmpb;
    d[k+1] = s*s*tmp1 + c*c*tmp2 + 2.0*c*s*tmpb;
    e[k] = c*s*(tmp1-tmp2) + (c*c-s*s)*tmpb;
    e[k-1] = c*e[k-1] - s*x;
    if(k!=kmax-1){
      // printf("if %d,%d, yes\n",k,kmax-1);
      x = -s*e[k+1];
      e[k+1] *= c;
    }
    for(int i = 0;i<nk;i++){
      Float Qtmp1 = Q[i][k];
      Float Qtmp2 = Q[i][k+1];
      Q[i][k] = c*Qtmp1 - s*Qtmp2;
      Q[i][k+1] = s*Qtmp1 + c*Qtmp2;
    }
  }

}

void QRtrf_debug(Float *d, Float *e,int nk, int n, Float **Q, Float dsh, int kmin, int kmax){

#if 1
  n=5;
  nk=n;
  kmin=0;
  kmax=n-1;
  for(int i=1;i<=n;++i){
    d[i-1]=i*i+1.23; e[i-1]=i*i*i+5.432;
  }
  dsh=0;
#endif

  //debug
  Float* d0; Float* e0;
  {
    d0 = (Float*)malloc( sizeof(Float)* n);
    e0 = (Float*)malloc( sizeof(Float)* n);
    for(int i=0;i<n;++i){ d0[i]=d[i]; e0[i]=e[i];}
  }

  int k=kmin;
  Float den = 1.0/sqrt((d[k]-dsh)*(d[k]-dsh)+e[k]*e[k]);
  Float c = (d[k]-dsh)*den;
  Float s = -e[k]*den;

  Float tmp1 = d[k];
  Float tmp2 = d[k+1];
  Float tmpb = e[k];

  printf("debug k=%d %f %f %f %f %f %f\n",k,den, c,s,tmp1,tmp2,tmpb);
  
  //       TDa(k)   = c*c*tmpa1 + s*s*tmpa2 - 2.D0*c*s*tmpb
  d[k] = c*c*tmp1 + s*s*tmp2 - 2.0*c*s*tmpb;
  //       TDa(k+1) = s*s*tmpa1 + c*c*tmpa2 + 2.D0*c*s*tmpb
  d[k+1] = s*s*tmp1 + c*c*tmp2 + 2.0*c*s*tmpb;
  //       TDb(k)   = c*s*(tmpa1-tmpa2) + (c*c-s*s)*tmpb
  e[k] = c*s*(tmp1-tmp2) + (c*c-s*s)*tmpb;
  //       x = -s*TDb(k+1)
  Float x = -s*e[k+1];
  //       TDb(k+1) = c*TDb(k+1)
  e[k+1] *= c;
  for(int i=0;i<nk;i++){
    Float Qtmp1 = Q[i][k];
    Float Qtmp2 = Q[i][k+1];
    Q[i][k] = c*Qtmp1 - s*Qtmp2;
    Q[i][k+1] = s*Qtmp1 + c*Qtmp2; 
  }

  for(int i=0;i<n;++i)
    printf("debug k=%d %d %f %f\n",k,i,d[i],e[i]);
  

  //debug
  {
    //Float *QA= (Float*)malloc(sizeof(Float)*n*n);
    Float tmp;
    printf("after QA %d\n",k);
    for(int i=0;i<n;++i){
      printf("%d ",i);
      for(int j=0;j<n;++j){

	//   Qdag[i][k]   (A-dsh)[k][j]
	// = + Q[j-1][i]  (A-dsh)[j-1][j]
	// = + Q[j  ][i]  (A-dsh)[j  ][j]
	// = + Q[j+1][i]  (A-dsh)[j+1][j]
	//
	//   (A-dsh)[j-1][j]  = e0[j-1]
	//   (A-dsh)[j  ][j]  = d0[j]-dsh 
	//   (A-dsh)[j+1][j]  = e0[j]
#if 1
	tmp = Q[j][i]* (d0[j]-dsh);
	if( j>1)
	  tmp += Q[j-1][i]* e0[j-1];
	if(j<n-1)
	  tmp += Q[j+1][i]* e0[j];
	//QA[n*i+j]=tmp;
	printf("%e ",tmp);
#else
	tmp = Q[i][j]* (d0[j]-dsh);
	if( j>1)
	  tmp += Q[i][j-1]* e0[j-1];
	if(j<n-1)
	  tmp += Q[i][j+1]* e0[j];
	QA[n*i+j]=tmp;
	printf("%e ",tmp);
#endif
      }
      printf("\n");
    }
  }//debug


  // givens rotations
  for(int k = kmin+1; k < kmax; k++){
    //       Fden = 1.D0/DSQRT( x**2 + TDb(k-1)**2 )
    den = 1.0/sqrt(x*x + e[k-1]*e[k-1]);
    //       c = TDb(k-1)*Fden
    c = e[k-1]*den;
    //       s = - x*Fden
    s = - x*den;

    //        tmpa1 = TDa(k)
    tmp1 = d[k];
    //        tmpa2 = TDa(k+1)
    tmp2 = d[k+1];
    //        tmpb  = TDb(k)
    tmpb = e[k];
    //       TDa(k)   = c*c*tmpa1 + s*s*tmpa2 - 2.D0*c*s*tmpb
    d[k] = c*c*tmp1 + s*s*tmp2 - 2.0*c*s*tmpb;
    //       TDa(k+1) = s*s*tmpa1 + c*c*tmpa2 + 2.D0*c*s*tmpb
    d[k+1] = s*s*tmp1 + c*c*tmp2 + 2.0*c*s*tmpb;
    //       TDb(k)   = c*s*(tmpa1-tmpa2) + (c*c-s*s)*tmpb
    e[k] = c*s*(tmp1-tmp2) + (c*c-s*s)*tmpb;
    //       TDb(k-1) = c*TDb(k-1) - s*x
    e[k-1] = c*e[k-1] - s*x;
    
    if(k!=kmax-1){
      //x = -s*TDb(k+1)
      x = -s*e[k+1];
      // TDb(k+1) = c*TDb(k+1)
      e[k+1] *= c;
    }
    for(int i = 0;i<nk;i++){
      Float Qtmp1 = Q[i][k];
      Float Qtmp2 = Q[i][k+1];
      Q[i][k] = c*Qtmp1 - s*Qtmp2;
      Q[i][k+1] = s*Qtmp1 + c*Qtmp2;
    }

    for(int i=0;i<n;++i)
    printf("debug k=%d %d %f %f\n",k,i,d[i],e[i]);

    
  //debug
  {
    //Float *QA= (Float*)malloc(sizeof(Float)*n*n);
    Float tmp;
    printf("after QA %d\n",k);
    for(int i=0;i<n;++i){
      printf("%d ",i);
      for(int j=0;j<n;++j){

	//   Qdag[i][k]   (A-dsh)[k][j]
	// = + Q[j-1][i]  (A-dsh)[j-1][j]
	// = + Q[j  ][i]  (A-dsh)[j  ][j]
	// = + Q[j+1][i]  (A-dsh)[j+1][j]
	//
	//   (A-dsh)[j-1][j]  = e0[j-1]
	//   (A-dsh)[j  ][j]  = d0[j]-dsh 
	//   (A-dsh)[j+1][j]  = e0[j]
#if 1
	tmp = Q[j][i]* (d0[j]-dsh);
	if( j>1)
	  tmp += Q[j-1][i]* e0[j-1];
	if(j<n-1)
	  tmp += Q[j+1][i]* e0[j];
	//QA[n*i+j]=tmp;
	printf("%e ",tmp);
#else
	tmp = Q[i][j]* (d0[j]-dsh);
	if( j>1)
	  tmp += Q[i][j-1]* e0[j-1];
	if(j<n-1)
	  tmp += Q[i][j+1]* e0[j];
	QA[n*i+j]=tmp;
	printf("%f ",tmp);
#endif
      }
      printf("\n");
    }
  }//debug


  printf("debug Q at k=%d\n",k);
  for(int i=0;i<n;++i){
    for(int j=0;j<n;++j) printf("%f ",Q[i][j]);
    printf("\n");
  }


  
  }//loop(k)


  
  //debug
  {

    Float *QA= (Float*)malloc(sizeof(Float)*n*n);
    Float tmp;
    printf("final Q^T A\n");
    for(int i=0;i<n;++i){
      printf("%d ",i);
      for(int j=0;j<n;++j){

#if 1
	//   Qdag[i][k]   (A-dsh)[k][j]
	// = + Q[j-1][i]  (A-dsh)[j-1][j]
	// = + Q[j  ][i]  (A-dsh)[j  ][j]
	// = + Q[j+1][i]  (A-dsh)[j+1][j]
	//
	//   (A-dsh)[j-1][j]  = e0[j-1]
	//   (A-dsh)[j  ][j]  = d0[j]-dsh 
	//   (A-dsh)[j+1][j]  = e0[j]
#if 1
	tmp = Q[j][i]* (d0[j]-dsh);
	if( j>0)
	  tmp += Q[j-1][i]* e0[j-1];
	if(j<n-1)
	  tmp += Q[j+1][i]* e0[j];
	QA[n*i+j]=tmp;
	printf("%e ",tmp);
#else
	tmp = Q[i][j]* (d0[j]-dsh);
	if( j>1)
	  tmp += Q[i][j-1]* e0[j-1];
	if(j<n-1)
	  tmp += Q[i][j+1]* e0[j];
	QA[n*i+j]=tmp;
	printf("%e ",tmp);
#endif
#else
	//    (A-dsh)[i][k] Qdag[k][j]
	// = + (A-dsh)[i][i-1] Q[j][i-1] 
	// = + (A-dsh)[i][i]   Q[j  ][i]  
	// = + (A-dsh)[i][i+1] Q[j][i+1]  
	//
	//   (A-dsh)[i][i-1]  = e0[i-1]
	//   (A-dsh)[i  ][i]  = d0[i]-dsh 
	//   (A-dsh)[i][i+1]  = e0[i]

#if 1
	tmp = (d0[i]-dsh)*Q[i][j]; 
	if( i>1)
	  tmp += e0[i-1]*Q[j][i-1];
	if(i<n-1)
	  tmp += e0[i]*Q[j][i+1];
	QA[n*i+j]=tmp;
	printf("%e ",tmp);
#else
	tmp = Q[i][j]* (d0[j]-dsh);
	if( j>1)
	  tmp += Q[i][j-1]* e0[j-1];
	if(j<n-1)
	  tmp += Q[i][j+1]* e0[j];
	QA[n*i+j]=tmp;
	printf("%e ",tmp);
#endif
#endif
      }
      printf("\n");   
    }
    

    printf("QA Q\n");
    for(int i=0;i<n;++i){
      printf("%d ",i);
      for(int j=0;j<n;++j){
	//   Qdag[i][k]   (A-dsh)[k][j]
	// = + Q[j-1][i]  (A-dsh)[j-1][j]
	// = + Q[j  ][i]  (A-dsh)[j  ][j]
	// = + Q[j+1][i]  (A-dsh)[j+1][j]
	//
	//   (A-dsh)[j-1][j]  = e0[j-1]
	//   (A-dsh)[j  ][j]  = d0[k]-dsh 
	//   (A-dsh)[j+1][j]  = e0[j]

#if 1
	tmp=0;
	for(int k=0;k<n;++k)
	  tmp += QA[n*i+k]*Q[k][j];
	printf("%e ",tmp);
#else
	tmp=0;
	for(int k=0;k<n;++k)
	  tmp += QA[n*i+k]*Q[j][k];
	printf("%e ",tmp);
#endif
      }
      printf("\n");   
    }

    
    free(d0); free(e0);
    free(QA);
  }
  
}



int tqli(Float *d, Float *e, int n, Float **z, int numshifts)
{
    int     m, l, iter, i, k;
    Float  s, r, p, g, f, dd, c, b;

    for (i = 1; i < n; i++) e[i - 1] = e[i];
    e[n-1] = 0.0;
    for (l = 0; l < n; l++)
      {
        iter = 0;
        do{
            for (m = l; m < n - 1; m++)
	      {
                dd = fabs(d[m]) + fabs(d[m + 1]);
                if (fabs(e[m]) + dd == dd)
		  break;
	      }
            if (m != l)
	      {
                if (iter++ == numshifts)
		  return (0);
                g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                r = sqrt((g * g) + 1.0);
                g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                s = c = 1.0;
                p = 0.0;
                for (i = m - 1; i >= l; i--)
		  {
                    f = s * e[i];
                    b = c * e[i];
                    if (fabs(f) >= fabs(g))
		      {
                        c = g / f;
                        r = sqrt((c * c) + 1.0);
                        e[i + 1] = f * r;
                        c *= (s = 1.0 / r);
		      } else
		      {
                        s = f / g;
                        r = sqrt((s * s) + 1.0);
                        e[i + 1] = g * r;
                        s *= (c = 1.0 / r);
		      }
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;
                    for (k = 0; k < n; k++)
		      {
                        f = z[k][i + 1];
                        z[k][i + 1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
		      }
		  }
                d[l] = d[l] - p;
                e[l] = g;
                e[m] = 0.0;
	      }
	    } while (m != l);
      }
    
    return (1);
}
int tqli(Float *d, Float *e, int n, Float **z, int numshifts, Float* shifts)
{
    int     m, l, iter, i, k;
    Float  s, r, p, g, f, dd, c, b;

    //sort shifts into reverse order, smallest to largest
    eigsrt_low(shifts,n);
    //eigsrt(shifts,n);

    for (i = 1; i < n; i++) e[i - 1] = e[i];
    e[n-1] = 0.0;
    for (l = 0; l < numshifts; l++)
      {
        //iter = 0;
        //do{
	  // deflate
	  for (m = l; m < n - 1; m++)
	    {
	      dd = fabs(d[m]) + fabs(d[m + 1]);
	      if (fabs(e[m]) + dd == dd)
		break;
	    }
	  //m=n-1;
            if (m != l)
	      {
                //if (iter++ == numshifts)
		//return (0);
                //g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                //r = sqrt((g * g) + 1.0);
                //g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                g = d[m] - shifts[l];
                s = c = 1.0;
                p = 0.0;
                for (i = m - 1; i >= l; i--)
		  {
                    f = s * e[i];
                    b = c * e[i];
                    if (fabs(f) >= fabs(g))
		      {
                        c = g / f;
                        r = sqrt((c * c) + 1.0);
                        e[i + 1] = f * r;
                        c *= (s = 1.0 / r);
		      } else
		      {
                        s = f / g;
                        r = sqrt((s * s) + 1.0);
                        e[i + 1] = g * r;
                        s *= (c = 1.0 / r);
		      }
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;
                    for (k = 0; k < n; k++)
		      {
                        f = z[k][i + 1];
                        z[k][i + 1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
		      }
		  }
                d[l] = d[l] - p;
                e[l] = g;
                e[m] = 0.0;
	      }
	    //} while (m != l);
      }
    
    return (1);
}





// implicitly shifted QR iteration of degree "number of shifts"
// for a symmetric tridiagonal matrix
// steps = number of shifts
// d diagonal
// e sub/super-diagonal (e[n-1] arbitrary)
// s shifts (there are "steps" of them).
// Q unitary matrix initialized to 1 (Q[i] is the ith column of Q)
int tqri(Float *d, Float *e, Float *s, int n, Float **Q, int steps)
{
  Float *temp = (Float *) smalloc("eigen_solver.C","tqri", "temp", n * sizeof(Float));
  Float *alpha = (Float *) smalloc("eigen_solver.C","tqri", "alpha", n * sizeof(Float));
  Float *beta = (Float *) smalloc("eigen_solver.C","tqri", "beta", n * sizeof(Float));
  
  LRG.SetSigma(0.5);
  //Float EPS = DBL_EPSILON;
  Float EPS = 1e-8;

  //sort shifts into reverse order, smallest to largest
  eigsrt_low(s,n);
  //for(int j=0;j<n;j++) printf("Sorted Shifts: %d %.16e\n",j, s[j]);
  //debug
  //s[0]=-3.2;
  //s[1]=2.07;
  // form the vector p(T) e_1 = (T-shift_0)...(T-shift_(numshifts-1)) (1,0,0...)^T
  // Q0 = (T-shift_(numshifts-1)) (1,0,0...)^T
  if(steps>0){
    //*((Float*)Q[0]) = d[0]-s[0];
    //*((Float*)Q[0]+1) = e[0];
    Q[0][0] = d[0]-s[0];
    Q[0][1] = e[0];
    for(int i=1;i<steps;i++){
      // temp = A.Q[0]
      mult_vec_by_tridiag(temp, Q[0], d, e, i+2);
      // Q[0] = -shift*Q[0]+temp
      fTimesV1PlusV2(Q[0], -s[i], Q[0], temp, i+2);    
    }
    //normalize
    Float norm = 1.0/sqrt(vecNorm(Q[0],steps+1));
    vecTimesEquFloat(Q[0],norm,steps+1);
  }else{
    Q[0][0] = 1.0;
  }
  // now form unitary Q (Q_0 is first column). Successive columns j are T^j x
  // following Watkins (2004) (same as in Numerical Recipes for QL algorithm)
  // start with first vector (column) of Q_0 = x, recursion to get Q_1, Q_2, ..., Q_m 
  // build new tridiagonal matrix along the way
  // column 1 (2nd column)
  mult_vec_by_tridiag(Q[1], Q[0], d, e, n);
  for(int i=0;i<n;i++)printf("Q1 %d %.16e \n",i,Q[1][i]);
  alpha[0] = dotProduct(Q[0],Q[1],n);
  fTimesV1PlusV2(Q[1], -alpha[0],  Q[0], Q[1], n);
  for(int i=0;i<n;i++)printf("Q1 %d %.16e \n",i,Q[1][i]);
  beta[0] = sqrt(vecNorm(Q[1],n));
  if(beta[0]<EPS){
    for(int i=0;i<n;i++){
      Q[1][i] = LRG.Grand();
    }
    lanczos_GramSchm(Q[1], Q, 1, n, 0); 
    Float norm = 1.0/sqrt(vecNorm(Q[1],n));
    vecTimesEquFloat(Q[1],norm,n);   
  }else{
    vecTimesEquFloat(Q[1],1.0/beta[0],n);
  }
  //for(int i=0;i<n;i++)printf("Q0 Q1 %d %e %e \n",i,Q[0][i],Q[1][i]);
  // now the rest (chase the bulge)
  for(int j=2;j<n;j++){
    mult_vec_by_tridiag(Q[j], Q[j-1], d, e, n);
    alpha[j-1] = dotProduct(Q[j-1],Q[j],n);
    fTimesV1PlusV2(Q[j], -beta[j-2], Q[j-2], Q[j], n);
    fTimesV1PlusV2(Q[j], -alpha[j-1], Q[j-1], Q[j], n);
    for(int i=0;i<n;i++)printf("Q[%d] %d %.16e \n",j,i,Q[j][i]);
    beta[j-1] = sqrt(vecNorm(Q[j],n));
    if(beta[j-1]<EPS){
      for(int i=0;i<n;i++){
	Q[j][i] = LRG.Grand();
      }
      lanczos_GramSchm(Q[j], Q, j, n, 0); 
      Float norm = 1.0/sqrt(vecNorm(Q[j],n));
      vecTimesEquFloat(Q[j],norm,n);   
    }else{
      vecTimesEquFloat(Q[j],1.0/beta[j-1],n);
    }
  }
  // finish up with last diagonal element
  if(n>1){
    mult_vec_by_tridiag(temp, Q[n-1], d, e, n);    
    alpha[n-1] = dotProduct(Q[n-1],temp,n);
  }
  for(int j=0;j<n;j++){
    d[j]=alpha[j];
  }
  for(int j=0;j<n-1;j++){
    e[j]=beta[j];
  }
  sfree(beta);
  sfree(alpha);
  sfree(temp);

}


// "Watkins" version. Orthogonal transform with the first column (row) =  p(A)e_1
// then use tred2 from NR to "chase the bulge"
//  
int tqri2(Float *d, Float *e, Float *s, int n, Float **Q, int steps)
{
  Float *temp = (Float *) smalloc("eigen_solver.C","tqri2", "temp", n * sizeof(Float));
  Float *temp2 = (Float *) smalloc("eigen_solver.C","tqri2", "temp2", n * sizeof(Float));
  Float **A = (Float **) smalloc("eigen_solver.C","tqri2", "A", (n)* sizeof(Float*));
  for(int i=0;i<n;i++)
    A[i] = (Float *) smalloc("eigen_solver.C","tqri2", "A", (n)* sizeof(Float*));
  Float **scratch = (Float **) smalloc("eigen_solver.C","tqri2", "scratch", (n)* sizeof(Float*));
  for(int i=0;i<n;i++)
    scratch[i] = (Float *) smalloc("eigen_solver.C","tqri2", "scratch", (n)* sizeof(Float*));

  LRG.SetSigma(0.5);
  Float EPS = DBL_EPSILON;

  //sort shifts into reverse order, smallest to largest
  eigsrt_low(s,n);

  if(steps>0){
    Q[0][0] = d[0]-s[0];
    Q[0][1] = Q[1][0]=e[0];
    for(int i=1;i<steps;i++){
      // temp = A.Q[0]
      mult_vec_by_tridiag(temp, Q[0], d, e, i+2);
      // Q[0] = -shift*Q[0]+temp
      fTimesV1PlusV2(Q[0], -s[i], Q[0], temp, i+2);    
    }
    //normalize
    Float norm = 1.0/sqrt(vecNorm(Q[0],steps+1));
    vecTimesEquFloat(Q[0],norm,steps+1);
  }else{
    return(0);//Q[0][0] = 1.for;
  }
  for(int j=0;j<n;j++){
    printf("X= %d %e \n",j,Q[0][j]);
  }

  // reflector Q0 = 1 - 2 u u*, u=(e_1-x)/2u0, |u|=1
  Float u0= sqrt(0.5*(1.0-Q[0][0]));
  temp[0]=0.5*(1.0-Q[0][0])/u0;
  for(int i=1;i<n;i++){
    temp[i]=-0.5*(Q[0][i])/u0;
  }
  // A Q_0 = A - 2 A u u*
  // A u:
  mult_vec_by_tridiag(temp2, temp, d, e, n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      scratch[i][j] = -2.0*temp2[i]*temp[j];
      if(i==j)scratch[i][j]+=d[i];
      if(i==j-1)scratch[i][j]+=e[i];
      if(i==j+1)scratch[i][j]+=e[j];
    }
  }
  // Q_0^T A Q_0 = A Q_0 - 2 u u* A Q_0
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      Float sum = 0.0;
      for(int k=0;k<n;k++){
	sum += temp[k]*scratch[k][j];
      }
      A[i][j] = scratch[i][j]-2.0* sum * temp[i];
      //printf("A %d %d %e\n",i,j,A[i][j]);
    }
  }
  // now reduce A to tridiag again
  tred2(A, n, d, e);
  // tred2 takes e[0] to be arbitrary
  for(int j=1;j<n;j++){
    e[j-1]=e[j];
  }
  /*
  printf("Q house\n");
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      printf(" %e ",A[i][j]);
    }
    printf("\n");
  }
  */

  // Q= Q0 Q1...Qn-2 = (1-2 uu*)Q
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      Float sum = 0.0;
      for(int k=0;k<n;k++){
	sum += temp[k]*A[j][k];
	//sum += temp[k]*A[k][j];
      }
      scratch[i][j] = A[j][i]-2.0*sum*temp[i];
    }
  }

  // my Q is stored in transpose order
  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){
      Q[j][i] = scratch[i][j];
      //printf("Q %d%d %e %e \n",j,i,Q[j][i], scratch[i][j]);
    }
  }
  /*
  printf("Q total\n");
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      printf(" %e ",Q[i][j]);
    }
    printf("\n");
  }
  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){
      scratch[i][j]=0.0;
      for(int k=0;k<n;k++){
	scratch[i][j] += Q[k][i]*Q[k][j];
      }
      printf("QQ* %d%d %e\n",i,j,scratch[i][j]);
    }
  }
  */
  for(int ii=n-1;ii>=0;ii--){
    sfree(scratch[ii]);
  }
  sfree(scratch);
  for(int ii=n-1;ii>=0;ii--){
    sfree(A[ii]);
  }
  sfree(A);
  sfree(temp2);
  sfree(temp);
}


// "Watkins" version. Orthogonal transform with the first column (row) =  p(A)e_1
// then use tred2 from NR to "chase the bulge"
// form Q from krylov subspace
int tqri3(Float *d, Float *e, Float *s, int n, Float **Q, int steps)
{
  Float *temp = (Float *) smalloc("eigen_solver.C","tqri2", "temp", n * sizeof(Float));
  Float *temp2 = (Float *) smalloc("eigen_solver.C","tqri2", "temp2", n * sizeof(Float));
  Float **A = (Float **) smalloc("eigen_solver.C","tqri2", "A", (n)* sizeof(Float*));
  for(int i=0;i<n;i++)
    A[i] = (Float *) smalloc("eigen_solver.C","tqri2", "A", (n)* sizeof(Float*));
  Float **scratch = (Float **) smalloc("eigen_solver.C","tqri2", "scratch", (n)* sizeof(Float*));
  for(int i=0;i<n;i++)
    scratch[i] = (Float *) smalloc("eigen_solver.C","tqri2", "scratch", (n)* sizeof(Float*));

  LRG.SetSigma(0.5);
  Float EPS = DBL_EPSILON;

  //sort shifts into reverse order, smallest to largest
  eigsrt_low(s,n);
  if(steps>0){
    Q[0][0] = d[0]-s[0];
    Q[0][1] = Q[1][0]=e[0];
    for(int i=1;i<steps;i++){
      // temp = A.Q[0]
      mult_vec_by_tridiag(temp, Q[0], d, e, i+2);
      // Q[0] = -shift*Q[0]+temp
      fTimesV1PlusV2(Q[0], -s[i], Q[0], temp, i+2);    
    }
    //normalize
    //Float norm = 1.0/sqrt(vecNorm(Q[0],steps+1));
    //vecTimesEquFloat(Q[0],norm,steps+1);
  }else{
    return(0);//Q[0][0] = 1.0;
  }

  // reflector Q0=q0,q1...qn-1, qi=A^i-1 qi-1
  for(int i=1;i<n;i++){
    mult_vec_by_tridiag(Q[i], Q[i-1], d, e, n);
    //Float norm = 1.0/sqrt(vecNorm(Q[i],n));
    //vecTimesEquFloat(Q[i],norm,n);
    //GramSchm(Q[i], Q, i, n);
  }
  for(int i=0;i<n;i++){
    Float norm = 1.0/sqrt(vecNorm(Q[i],n));
    vecTimesEquFloat(Q[i],norm,n);
  }
  printf("Q0\n");
  for(int i=0;i<n;i++){
    Float norm = sqrt(vecNorm(Q[i],n));
    printf("column %d (norm %e)  \n",i,norm );
    for(int j=0;j<n;j++){
      printf(" %e ",Q[i][j]);
    }
    printf("\n");
  }
  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){
      scratch[i][j]=0.0;
      for(int k=0;k<n;k++){
	scratch[i][j] += Q[k][i]*Q[k][j];
      }
      printf("QQ* %d%d %e\n",i,j,scratch[i][j]);
    }
  }
  exit(0);
  // A Q_0
  // Q0 is stored in transpose order
  for(int i=0;i<n;i++){
    mult_vec_by_tridiag(scratch[i], Q[i], d, e, n);
  }
  // Q_0^T A Q_0
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      Float sum = 0.0;
      for(int k=0;k<n;k++){
	sum += Q[i][k]*scratch[j][k];
      }
      A[i][j] = sum;
      //printf("A %d %d %e\n",i,j,A[i][j]);
    }
  }
  // now reduce A to tridiag again
  tred2(A, n, d, e);
  // tred2 takes e[0] to be arbitrary
  for(int j=1;j<n;j++){
    e[j-1]=e[j];
  }
  printf("Q house\n");
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      printf(" %e ",A[i][j]);
    }
    printf("\n");
  }
  //exit(0);
  // Q= Q0 Q1...Qn-2 = (1-2 uu*)Q
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      scratch[i][j]= 0.0;
      for(int k=0;k<n;k++){
	scratch[i][j] += Q[k][i]*A[k][j];
      }
    }
  }

  // my Q is stored in transpose order
  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){
      Q[j][i] = scratch[i][j];
      //printf("Q %d%d %e %e \n",j,i,Q[j][i], scratch[i][j]);
    }
  }
  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){
      scratch[i][j]=0.0;
      for(int k=0;k<n;k++){
	scratch[i][j] += Q[k][i]*Q[j][k];
      }
      printf("QQ* %d%d %e\n",i,j,scratch[i][j]);
    }
  }

  for(int ii=n-1;ii>=0;ii--){
    sfree(scratch[ii]);
  }
  sfree(scratch);
  for(int ii=n-1;ii>=0;ii--){
    sfree(A[ii]);
  }
  sfree(A);
  sfree(temp2);
  sfree(temp);
}





IFloat vecNorm(const Float *a, int len){
  IFloat sum = 0.0;
  Float dum;
  for(int i = 0; i < len; ++i) {
    dum = *a++;
    sum += dum*dum;
  }
  return sum;
}

void mult_vec_by_tridiag(Float* out, Float *in, Float *d, Float *e, int n){

  out[0] = d[0]*in[0]+e[0]*in[1];
  for(int j=1;j<n-1;j++){
    //superdiagonal
    out[j]=e[j]*in[j+1];
    //diagonal
    out[j]+=d[j]*in[j];
    //subdiagonal
    out[j]+=e[j-1]*in[j-1];
  }
  if(n>1)
    out[n-1]=e[n-2]*in[n-2]+d[n-1]*in[n-1];
}

// orthogonalize psi w/r to vec's only one vector
void GramSchm_save(Float *psi, Float **vec, int Nvec, int f_size) 
{
  Float xp;

  for(int i = 0; i < Nvec; ++i)
    {
      xp = dotProduct(psi, vec[i], f_size);
      /* psi = psi - <vec[i],psi> vec[i] */
      fTimesV1PlusV2(psi, -xp, vec[i], psi, f_size);
    }
}

//#define USE_BLAS
#ifndef USE_BLAS
#define MOVE_FLOAT( pa, pb, n )  moveFloat(pa, pb, n)
#define VEC_TIMESEQU_FLOAT(py, fact, n ) vecTimesEquFloat( py, fact, n)
#define AXPY(n, fact, px, py)  fTimesV1PlusV2(py, fact, px, py, n)
#define glb_DDOT(n, px, py, p_dot) { *(p_dot)=dotProduct(px,py,n); glb_sum((p_dot)); }
#define ZDOT(n,px,py,p_dot) compDotProduct( p_dot, p_dot+1, px, py, n)

#define BASE float
#define CONJ_SIGN (-1.0)
#define INDEX int
#define CONST_REAL(a,i) (((const BASE *) a)[2*(i)])
#define CONST_IMAG(a,i) (((const BASE *) a)[2*(i)+1])
#define REAL(a,i) (((BASE *) a)[2*(i)])
#define IMAG(a,i) (((BASE *) a)[2*(i)+1])
//#define ZDOTfloat(n,px,py,p_dot) cblas_cdotc_sub( n/2,px,1, py,1, p_dot)
#if 1
inline void
ZDOTfloat (const int N, const float *X,  const float *Y,
             float *result)
{
IFloat c_r,c_i;
compDotProduct<float,float>(&c_r,&c_i,X,Y,N);
*result = c_r;
*(result+1) = c_i;
}
#else
inline void
ZDOTfloat (const int N, const float *X,  const float *Y,
             float *result)
{
  BASE r_real = 0.0;
  BASE  r_imag = 0.0;
  INDEX i;
  INDEX ix = 0;
  INDEX iy = 0;
  for (i = 0; i < (N/2); i++) {
    const BASE x_real = CONST_REAL(X, ix);
    const BASE x_imag = CONST_IMAG(X, ix);
    const BASE y_real = CONST_REAL(Y, iy);
    const BASE y_imag = CONST_IMAG(Y, iy);
    r_real += x_real * y_real - CONJ_SIGN * x_imag * y_imag;
    r_imag += x_real * y_imag + CONJ_SIGN * x_imag * y_real;
    ix += 1;
    iy += 1;
  }
  result[0] = (float) r_real;
  result[1] = (float) r_imag;
}
#endif
//#define ZAXPYfloat(n, fact, px, py)  cblas_caxpy(n/2, fact, px,1,py,1)
inline void
ZAXPYfloat (const int N, const float *alpha, const float *X,  float *Y)
{
  INDEX i;
  INDEX ix = 0;
  INDEX iy = 0;

  const BASE alpha_real = alpha[0];
  const BASE alpha_imag = alpha[1];

  if (fabs(alpha_real) == 0 && fabs(alpha_imag) == 0) {
    return;
  }

  for (i = 0; i < (N/2); i++) {
    const BASE x_real = CONST_REAL(X, ix);
    const BASE x_imag = CONST_IMAG(X, ix);
    REAL(Y, iy) += (alpha_real * x_real - alpha_imag * x_imag);
    IMAG(Y, iy) += (alpha_real * x_imag + alpha_imag * x_real);
    ix += 1;
    iy += 1;
  }
}
#undef CONJ_SIGN
#undef BASE
#undef INDEX
#undef CONST_REAL
#undef CONST_IMAG

#define ZAXPY(n, fact, px, py)   cTimesV1PlusV2(py, fact[0], fact[1], px, py, n);
#else
#include <util/qblas_extend.h>
#define MOVE_FLOAT( pa, pb, n )  cblas_dcopy(n, pb, 1, pa, 1)
#define VEC_TIMESEQU_FLOAT(py, fact, n ) cblas_dscal( n,  fact, py,1 )
#define AXPY(n, fact, px, py)  cblas_daxpy(n, fact, px,1,py,1)
#define glb_DDOT(n, px, py, p_dot) { *(p_dot) = cblas_ddot(n,px,py); glb_sum((p_dot)); }
#define ZDOT(n,px,py,p_dot) cblas_zdotc_sub( n/2, px,py, p_dot)
#define ZAXPY(n, fact, px, py)  cblas_zaxpy(n/2, fact, px,1,py,1)
#define ZDOTfloat(n,px,py,p_dot) cblas_cdotc_sub( n/2,px,1, py,1, p_dot)
//      ZDOTfloat(f_size, (float*)(vec[i]), (float*)vtmp, xp);
//      cblas_cdotc_sub(f_size/2, vec[i], 1, (float*)vtmp, 1, xp);
#define ZAXPYfloat(n, fact, px, py)  cblas_caxpy(n/2, fact, px,1,py,1)
#endif

// orthogonalize psi w/r to vec's only one vector
// treating the imaginary part of dot product exactly
void lanczos_GramSchm(Float *psi, Float **vec, int Nvec, int f_size, Float *alpha) 
{
  Float xp[2]; //xp_r, xp_i;
  
  //if(alpha)  {
  //compDotProduct(&xp_r, &xp_i, (Float*)(vec[Nvec-1]), (Float*)psi ,  f_size);
  //   *alpha = xp_r;   //  Re ( vec[Nvec-1],  psi ) needed for Lanczos' alpha
  // }

  for(int i = 0; i<Nvec; ++i)
    //for(int i = Nvec-1; i>=0; --i)
    {
      //compDotProduct(xp, xp+1, (Float*)(vec[i]), (Float*)psi ,  f_size);
      ZDOT(f_size, (Float*)(vec[i]), (Float*)psi, xp);
      //glb_sum_five((Float *)&xp_r);
      //glb_sum_five((Float *)&xp_i);
      //slice_sum(Float * float_p, int blcklength, int dir);
      slice_sum(xp, 2, 1970);
      //if(!UniqueID())printf("[%d][%d] %e %e\n", UniqueID(),i, xp[0],xp[1]);

      /* psi = psi - <vec[i],psi> vec[i] */
      //cTimesV1PlusV2(psi, -xp[0], -xp[1], vec[i], psi, f_size);
      xp[0] =-xp[0]; xp[1] = -xp[1];
      ZAXPY(f_size, xp, vec[i], psi);
      // printf("GS %d %.16e %16e\n", i+1, xp_r, xp_i);

      if(i==Nvec-1 && alpha)  *alpha = xp[0];   //  Re ( vec[Nvec-1],  psi ) needed for Lanczos' alpha
    }
  
}


void lanczos_GramSchm_test(Float *psi, Float **vec, int Nvec, int f_size, Float *alpha) 
{
  Float xp[2]; //xp_r, xp_i;
  
  //if(alpha)  {
  //compDotProduct(&xp_r, &xp_i, (Float*)(vec[Nvec-1]), (Float*)psi ,  f_size);
  //   *alpha = xp_r;   //  Re ( vec[Nvec-1],  psi ) needed for Lanczos' alpha
  // }

  for(int i = 0; i<Nvec; ++i)
    //for(int i = Nvec-1; i>=0; --i)
    {
      //compDotProduct(xp, xp+1, (Float*)(vec[i]), (Float*)psi ,  f_size);
      ZDOT(f_size, (Float*)(vec[i]), (Float*)psi, xp);
      //glb_sum_five((Float *)&xp_r);
      //glb_sum_five((Float *)&xp_i);
      //slice_sum(Float * float_p, int blcklength, int dir);
      slice_sum(xp, 2, 1970);
      //if(!UniqueID())printf("[%d][%d] %e %e\n", UniqueID(),i, xp[0],xp[1]);

      if(!UniqueID() && fabs(xp[1])> 1e-13)
	printf("[%d][%d] %e %e\n", UniqueID(),i, xp[0],xp[1]);
      
      /* psi = psi - <vec[i],psi> vec[i] */
      //cTimesV1PlusV2(psi, -xp[0], -xp[1], vec[i], psi, f_size);
      xp[0] =-xp[0]; xp[1] = -xp[1];
      ZAXPY(f_size, xp, vec[i], psi);
      // printf("GS %d %.16e %16e\n", i+1, xp_r, xp_i);

      if(i==Nvec-1 && alpha)  *alpha = xp[0];   //  Re ( vec[Nvec-1],  psi ) needed for Lanczos' alpha
    }
  
}


// float version
#if 0
void lanczos_GramSchm_test(Float *psi, float **vec, int Nvec, int f_size, Float *alpha) 
{
  Float xp[2]; //xp_r, xp_i;
  Vector *vtmp = (Vector *) smalloc("","lanczos_GramSchm_test", "r", f_size * sizeof(Float));
  
  for(int i = 0; i<Nvec; ++i)
    {
      mvfloattoFloat((Float*)vtmp, (float*)vec[i], f_size );
      //ZDOT(f_size, (Float*)(vec[i]), (Float*)psi, xp);
      ZDOT(f_size, (Float*)vtmp, (Float*)psi, xp);
      slice_sum(xp, 2, 1970);

      if(!UniqueID() && fabs(xp[1])> 1e-13)
	printf("[%d][%d] %e %e\n", UniqueID(),i, xp[0],xp[1]);
      
      /* psi = psi - <vec[i],psi> vec[i] */
      xp[0] =-xp[0]; xp[1] = -xp[1];
      //ZAXPY(f_size, xp, vec[i], psi);
      ZAXPY(f_size, xp, (Float*)vtmp, psi);

      if(i==Nvec-1 && alpha)  *alpha = xp[0];   //  Re ( vec[Nvec-1],  psi ) needed for Lanczos' alpha
    }
  sfree(vtmp);
}
#endif
// another float version (dot prods. are floats, not Floats)
void lanczos_GramSchm_test(Float *psi, float **vec, int Nvec, int f_size, Float *alpha) 
{

  float xp[2]; //xp_r, xp_i;
  Float xpd[2]; //xp_r, xp_i; glb sum in double
  for(int i=0;i<2;i++) xp[i]=xpd[i]=0.;
  Vector *vtmp = (Vector *) smalloc("","lanczos_GramSchm_test", "r", f_size * sizeof(float));

  mvFloattofloat((float*)vtmp, psi, f_size );
  
  for(int i = 0; i<Nvec; ++i)
    {
	float * ftmp1 = (float*)(vec[i]);
	float *ftmp2 = (float*)vtmp;
	VRB.Debug("","lanczos_GramSchm_test(before)","vec[%d]=%0.14e %0.14e vtmp=%0.14e %0.14e\n",i,
	*ftmp1,*(ftmp1+1), *ftmp2,*(ftmp2+1));

      ZDOTfloat(f_size, (float*)(vec[i]), (float*)vtmp, xp);
//      cblas_cdotc_sub(f_size/2, vec[i], 1, (float*)vtmp, 1, xp);

      if(fabs(xp[1])> 1e-13)
	VRB.Debug("","lanczos_GramSchm_test()","[%d][%d] %0.14e %0.14e\n", 0, i, xp[0],xp[1]);

      xpd[0] = (Float)xp[0];
      xpd[1] = (Float)xp[1];
      slice_sum((Float*)xpd, 2, 1970);
      xp[0] = (float)xpd[0];
      xp[1] = (float)xpd[1];
      
      /* psi = psi - <vec[i],psi> vec[i] */
      xp[0] =-xp[0]; xp[1] = -xp[1];
      //ZAXPY(f_size, xp, vec[i], psi);
      ZAXPYfloat(f_size, xp, vec[i], (float*)vtmp);
	VRB.Debug("","lanczos_GramSchm_test(after)","vec[%d]=%0.14e %0.14e vtmp=%0.14e %0.14e\n",i,
	*ftmp1,*(ftmp1+1), *ftmp2,*(ftmp2+1));


     if (VRB.IsActivated(VERBOSE_DEBUG_LEVEL))
     {
      ZDOTfloat(f_size, (float*)(vec[i]), (float*)vtmp, xp);
	VRB.Debug("","lanczos_GramSchm_test(after)","[%d][%d] %0.14e %0.14e\n", 0, i, xp[0],xp[1]);
     }

      if(i==Nvec-1 && alpha)  *alpha = xp[0];   //  Re ( vec[Nvec-1],  psi ) needed for Lanczos' alpha
    }

  mvfloattoFloat(psi, (float*)vtmp, f_size );

  sfree(vtmp);
}

void mvfloattoFloat(Float* out, float* in, int f_size)
{
#if 1
  float flt;
  for(int i=0;i<f_size;i++){
    flt = in[i];
    out[i] = (Float)flt;
  }
#endif
};
void mvFloattofloat(float* out, Float* in, int f_size)
{
#if 1
  float flt;
  for(int i=0;i<f_size;i++){
    flt = (float)in[i];
    out[i] = flt;
  }
#endif
};




// orthogonalize psi w/r to vec's only one vector
// assuming the imaginary part of dot product zero
void lanczos_GramSchm_real(Float *psi, Float **vec, int Nvec, int f_size, Float *alpha) 
{
  Float xp;

  for(int i = 0; i<Nvec; ++i)
    {
      glb_DDOT(f_size, (Float*)(vec[i]), (Float*)psi, &xp);

      /* psi = psi - <vec[i],psi> vec[i] */
      //cTimesV1PlusV2(psi, -xp[0], -xp[1], vec[i], psi, f_size);
      //xp[0] =-xp[0]; xp[1] = -xp[1];
      AXPY(f_size, -xp, vec[i], psi);
      // printf("GS %d %.16e %16e\n", i+1, xp_r, xp_i);

      if(i==Nvec-1 && alpha)  *alpha = xp;   //  Re ( vec[Nvec-1],  psi ) needed for Lanczos' alpha
    }


  
}


// d is eigenvalues, to be sorted
// v is corresponding vectors
void eigsrt(Float *d, Float **v, int n)
{
  int k,j,i;
  Float p;
  
  for (i=0;i<n-1;i++) {
    p=d[i];k=i;
    for (j=i+1;j<n;j++)
      if (d[j] > p){ p=d[j];k=j;}
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      // swap vectors
      for (j=0;j<n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}

// d is eigenvalues, to be sorted
// v is corresponding vectors
// sort smallest to biggest
void eigsrt_low(Float *d, int n)
{
  int k,j,i;
  Float p;
  
  for (i=0;i<n-1;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++)
      if (d[j] <= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
    }
  }
}

void eigsrt(Float *d, int n)
{
  int k,j,i;
  Float p;
  
  for (i=0;i<n-1;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
    }
  }
}

void eigsrt(Float *d, Float **v, int n, Vector **v2, int n2)
{
  int k,j,i;
  Float p;
  Vector *p2;
  Float *p3;
  
  for (i=0;i<n-1;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++)
      if (d[j] >= p){
	p=d[k=j];
      }
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      // swap vectors
      p2=v2[k];
      v2[k]=v2[i];
      v2[i]=p2;
      // swap vectors
      p3=v[k];
      v[k]=v[i];
      v[i]=p3;
      // swap vectors
      /*for (j=0;j<n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }      
      
      for (j=0;j<n2;j++) {
	p=*((Float*)v2[i]+j);
	*((Float*)v2[i]+j)=*((Float*)v2[k]+j);
	*((Float*)v2[k]+j)=p;
      }
      */
    }
  }
}

void eigsrt(Float *d, Float *e, Float **v, int n, Vector **v2, int n2)
{
  int k,j,i;
  Float p;
  Vector *p2;
  Float p3;
 
  for (i=0;i<n-1;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++)
      if (d[j] <= p){
	p=d[k=j];
	p3=e[j];
      }
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      e[k]=e[i];
      e[i]=p3;
      // swap vectors
      p2=v2[k];
      v2[k]=v2[i];
      v2[i]=p2;
      for (j=0;j<n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}
void eigsrt(Float *d, Float *e, Float *tte, Float **v, int n, Vector **v2, int n2)
{
  int k,j,i;
  Float p;
  Vector *p2;
  Float p3,p4;
 
  for (i=0;i<n-1;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++)
      if (d[j] <= p){
	p=d[k=j];
	p3=e[j];
	p4=tte[j];
      }
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      e[k]=e[i];
      e[i]=p3;
      tte[k]=tte[i];
      tte[i]=p4;
      // swap vectors
      p2=v2[k];
      v2[k]=v2[i];
      v2[i]=p2;
      for (j=0;j<n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}

void eigsrt_low(Float *d, Float **v, int n, Vector **v2, int n2)
{
  int k,j,i;
  Float p;
  Vector *p2;
  Float p3,p4;
 
  for (i=0;i<n-1;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++)
      if (d[j] <= p){
	p=d[k=j];
      }
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      // swap vectors
      p2=v2[k];
      v2[k]=v2[i];
      v2[i]=p2;
      for (j=0;j<n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}

/* 
Finds all eigenvalues of an upper Hessenberg matrix a[0..n-1][0..n-1]. 
On input a can be exactly as output from elmhes and eltra 11.6. 
On output, wri[0..n-1] contains the eigen- values of a, while zz[0..n-1][0..n-1] is a matrix 
whose columns contain the corresponding eigenvectors. The eigenvalues are not sorted, 
except that complex conjugate pairs appear consecutively with the eigenvalue having 
the positive imaginary part first. For a complex eigenvalue, only the eigenvector corresponding 
to the eigenvalue with positive imaginary part is stored, with real part in zz[0..n-1][i] 
and imaginary part in h.zz[0..n-1][i+1]. The eigenvectors are not normalized.
*/
/*
void hqr2()
{
  int nn,m,l,k,j,its,i,mmin,na; Doub z,y,x,w,v,u,t,s,r,q,p,anorm=0.0,ra,sa,vr,vi;
  const Doub EPS=numeric_limits<Doub>::epsilon();
  for (i=0;i<n;i++)            //Compute matrix norm for possible use in lo- 
    for (j=MAX(i-1,0);j<n;j++) //cating single small subdiagonal element.
      anorm += abs(a[i][j]);
  nn=n-1; 
  t=0.0;                       //Gets changed only by an exceptional shift. 
  while (nn >= 0) {            //Begin search for next eigenvalue.
    its=0; 
    do {
      for (l=nn;l>0;l--) {                 //Begin iteration: look for single small subdi-
	s=abs(a[l-1][l-1])+abs(a[l][l]);   //agonal element. 
	if (s == 0.0) s=anorm;
	if (abs(a[l][l-1]) <= EPS*s) {
	  a[l][l-1] = 0.0; 
	  break;
	}
      } 
      x=a[nn][nn]; 
      if (l == nn) {	//One root found.
	wri[nn]=a[nn][nn]=x+t;
	nn--; 
      } else {
	y=a[nn-1][nn-1]; 
	w=a[nn][nn-1]*a[nn-1][nn]; 
	if (l == nn-1) {	//Two roots found...
	  p=0.5*(y-x); 
	  q=p*p+w;
	  z=sqrt(abs(q)); 
	  x += t; 
	  a[nn][nn]=x; 
	  a[nn-1][nn-1]=y+t; 
	  if (q >= 0.0) { //...a real pair.
	    z=p+SIGN(z,p); 
	    wri[nn-1]=wri[nn]=x+z; 
	    if (z != 0.0) wri[nn]=x-w/z; 
	    x=a[nn][nn-1]; 
	    s=abs(x)+abs(z); 
	    p=x/s; q=z/s; 
	    r=sqrt(p*p+q*q); 
	    p /= r; 
	    q /= r; 
	    for (j=nn-1;j<n;j++) {	//Row modification.
	      z=a[nn-1][j]; 
	      a[nn-1][j]=q*z+p*a[nn][j]; 
	      a[nn][j]=q*a[nn][j]-p*z;
	    } 
	    for (i=0;i<=nn;i++) {	//Column modification.
	      z=a[i][nn-1]; 
	      a[i][nn-1]=q*z+p*a[i][nn]; 
	      a[i][nn]=q*a[i][nn]-p*z;
	    } 
	    for (i=0;i<n;i++) {	//Accumulate transformations.
	      z=zz[i][nn-1]; 
	      zz[i][nn-1]=q*z+p*zz[i][nn]; 
	      zz[i][nn]=q*zz[i][nn]-p*z;
	    } 
	  } else {	.//..a complex pair.
	    wri[nn]=Complex(x+p,-z); 
	    wri[nn-1]=conj(wri[nn]);
	  }
	  nn -= 2; 
	} else {	//No roots found. Continue iteration.
	  if (its == 30) throw("Too many iterations in hqr"); 
	  if (its == 10 || its == 20) {	//Form exceptional shift.
	    t += x; 
	    for (i=0;i<nn+1;i++) a[i][i] -= x; 
	    s=abs(a[nn][nn-1])+abs(a[nn-1][nn-2]); 
	    y=x=0.75*s; w = -0.4375*s*s;
	  } 
	  ++its; 
	  for (m=nn-2;m>=l;m--) { //Form shift and then look for 
	    z=a[m][m];            //2 consecutive small sub- diagonal elements.
	    r=x-z; 
	    s=y-z; 
	    p=(r*s-w)/a[m+1][m]+a[m][m+1]; //Equation (Webnote 16.21).
	    q=a[m+1][m+1]-z-r-s; 
	    r=a[m+2][m+1]; 
	    s=abs(p)+abs(q)+abs(r); //Scale to prevent overflow or underflow.
	    p /= s;
	    q /= s; 
	    r /= s; 
	    if (m == l) break; 
	    u=abs(a[m][m-1])*(abs(q)+abs(r)); 
	    v=abs(p)*(abs(a[m-1][m-1])+abs(z)+abs(a[m+1][m+1]));
	    if (u <= EPS*v) break;  //Equation (Webnote 16.24).
	  } 
	  for (i=m;i<nn-1;i++) {
	    a[i+2][i]=0.0; 
	    if (i != m) a[i+2][i-1]=0.0;
	  } 
	  for (k=m;k<nn;k++) { 
	    //Double QR step on rows l to nn and columns m to nn.
	    if (k != m) { 
	      p=a[k][k-1];	// Begin setup of Householder 
	      q=a[k+1][k-1];	// vector. 
	      r=0.0;
	      if (k+1 != nn) r=a[k+2][k-1]; 
	      if ((x=abs(p)+abs(q)+abs(r)) != 0.0) {
		p /= x;	//Scale to prevent overflow or 
		q /= x;	//underflow. 
		r /= x;
	      }
	    } 
	    if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
	      if (k == m) { 
		if (l != m)
		  a[k][k-1] = -a[k][k-1]; 
	      } else
		a[k][k-1] = -s*x; 
	      p += s;		//Equations(Webnote16.22).
	      x=p/s; 
	      y=q/s; 
	      z=r/s; 
	      q /= p; 
	      r /= p; 
	      for (j=k;j<n;j++) { //Row modification.
		p=a[k][j]+q*a[k+1][j]; 
		if (k+1 != nn) {
		  p += r*a[k+2][j]; a[k+2][j] -= p*z;
		} 
		a[k+1][j] -= p*y; 
		a[k][j] -= p*x;
	      }
	      mmin = nn < k+3 ? nn : k+3; 
	      for (i=0;i<mmin+1;i++) { //Column modification.
		p=x*a[i][k]+y*a[i][k+1]; 
		if (k+1 != nn) {
		  p += z*a[i][k+2]; 
		  a[i][k+2] -= p*r;
		}
		a[i][k+1] -= p*q; 
		a[i][k] -= p;
	      }
	      for (i=0; i<n; i++) { //Accumulate transformations.
		p=x*zz[i][k]+y*zz[i][k+1]; 
		if (k+1 != nn) {
		  p += z*zz[i][k+2]; 
		  zz[i][k+2] -= p*r;
		}
		zz[i][k+1] -= p*q; 
		zz[i][k] -= p;
	      }
	    }
	  }
	}
      }
    } while (l+1 < nn);
  }
  //All roots found. Backsubstitute to find vectors of upper triangular form.
  if (anorm != 0.0) { 
    for (nn=n-1;nn>=0;nn--) {
      p=real(wri[nn]);
      q=imag(wri[nn]); 
      na=nn-1; 
      if (q == 0.0) {	//Real vector.
	m=nn; 
	a[nn][nn]=1.0; 
	for (i=nn-1;i>=0;i--) {
	  w=a[i][i]-p; 
	  r=0.0; 
	  for (j=m;j<=nn;j++)
	    r += a[i][j]*a[j][nn]; 
	  if (imag(wri[i]) < 0.0) {
	    z=w;
	    s=r; 
	  } else { 
	    m=i;
	  }
	  if (imag(wri[i]) == 0.0) { 
	    t=w;
	    if (t == 0.0) t=EPS*anorm;
	    a[i][nn]=-r/t; 
	  } else {	// Solve real equations.
	    x=a[i][i+1]; 
	    y=a[i+1][i]; 
	    q=SQR(real(wri[i])-p)+SQR(imag(wri[i])); 
	    t=(x*s-z*r)/q; 
	    a[i][nn]=t; 
	    if (abs(x) > abs(z))
	      a[i+1][nn]=(-r-w*t)/x; 
	    else
	      a[i+1][nn]=(-s-y*t)/z;
	  } 
	  t=abs(a[i][nn]); // Overflow control.
	  if (EPS*t*t > 1)
	    for (j=i;j<=nn;j++) a[j][nn] /= t;
	}
      } else if (q < 0.0) {                      //Complex vector, only do one case.
	m=na;                                    //Last vector component chosen imag-
	if (abs(a[nn][na]) > abs(a[na][nn])) {	 //inary so that eigenvec-
	    a[na][na]=q/a[nn][na];	         //tor matrix is triangular.
	  a[na][nn]=-(a[nn][nn]-p)/a[nn][na];
	} else {
	  Complex temp=Complex(0.0,-a[na][nn])/Complex(a[na][na]-p,q); 
	  a[na][na]=real(temp); 
	  a[na][nn]=imag(temp);
	}
	a[nn][na]=0.0; 
	a[nn][nn]=1.0; 
	for (i=nn-2;i>=0;i--) {
	  w=a[i][i]-p; 
	  ra=sa=0.0; 
	  for (j=m;j<=nn;j++) {
	    ra += a[i][j]*a[j][na]; 
	    sa += a[i][j]*a[j][nn];
	  } if (imag(wri[i]) < 0.0) {
	    z=w; 
	    r=ra; 
	    s=sa;
	  } else { 
	    m=i;
	    if (imag(wri[i]) == 0.0) {
	      Complex temp = Complex(-ra,-sa)/Complex(w,q); 
	      a[i][na]=real(temp); 
	      a[i][nn]=imag(temp);
	    } else {	//Solve complex equations. 
	      x=a[i][i+1];
	      y=a[i+1][i]; 
	      vr=SQR(real(wri[i])-p)+SQR(imag(wri[i]))-q*q; 
	      vi=2.0*q*(real(wri[i])-p); 
	      if (vr == 0.0 && vi == 0.0)
		vr=EPS*anorm*(abs(w)+abs(q)+abs(x)+abs(y)+abs(z)); 
	      Complex temp=Complex(x*r-z*ra+q*sa,x*s-z*sa-q*ra)/
		Complex(vr,vi); 
	      a[i][na]=real(temp); 
	      a[i][nn]=imag(temp); 
	      if (abs(x) > abs(z)+abs(q)) {
		a[i+1][na]=(-ra-w*a[i][na]+q*a[i][nn])/x;
		a[i+1][nn]=(-sa-w*a[i][nn]-q*a[i][na])/x; 
	      } else {
		Complex temp=Complex(-r-y*a[i][na],-s-y*a[i][nn])/ Complex(z,q);
		a[i+1][na]=real(temp); 
		a[i+1][nn]=imag(temp);
	      }
	    }
	  }
	  t=MAX(abs(a[i][na]),abs(a[i][nn])); 
	  if (EPS*t*t > 1)
	    for (j=i;j<=nn;j++) { 
	      a[j][na] /= t; 
	      a[j][nn] /= t;
	    }
	}
      }
    }
    for (j=n-1;j>=0;j--) //Multiply by transformation matrix to 
      for (i=0;i<n;i++) { //give vectors of original full matrix. 
	z=0.0;
	for (k=0;k<=j;k++) 
	  z += zz[i][k]*a[k][j];
	zz[i][j]=z;
      }
  }
}

void sort(){
  //Sort the eigenvalues in descending order of their real parts using straight insertion.
  int i; 
  for (int j=1;j<n;j++) {
    Complex x=wri[j]; 
    for (i=j-1;i>=0;i--) {
      if (real(wri[i]) >= real(x)) break; 
      wri[i+1]=wri[i];
    } 
    wri[i+1]=x;
  }
}

void sortvecs()
//Sort the eigenvalues in descending order of their real parts using straight insertion, and simul- 
//taneously rearrange the eigenvectors.
{
  int i; 
  VecDoub temp(n); 
  for (int j=1;j<n;j++) {
    Complex x=wri[j]; 
    for (int k=0;k<n;k++) 
      temp[k]=zz[k][j];
    for (i=j-1;i>=0;i--) { 
      if (real(wri[i]) >= real(x)) break; 
      wri[i+1]=wri[i]; 
      for (int k=0;k<n;k++)
	zz[k][i+1]=zz[k][i];
    } 
    wri[i+1]=x; 
    for (Int k=0;k<n;k++)
      zz[k][i+1]=temp[k];
  }
}

*/



CPS_END_NAMESPACE
