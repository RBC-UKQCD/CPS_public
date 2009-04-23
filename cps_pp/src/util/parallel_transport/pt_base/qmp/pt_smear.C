#ifdef USE_QMP
#include <string.h>
#include "asq_data_types.h"
#include "pt_int.h"

extern "C" {
void vaxpy3(vector *res,Float *scale,vector *mult,vector *add, int ncvec);
  inline void vaxpy3_m(matrix *res,Float *scale,matrix *mult,matrix *add, 
  int ncvec){
    vaxpy3((vector *)res, scale, (vector *)mult,(vector *)add,ncvec);
  }
}

enum {NUM_DIR=8,POS_DIR=4};

static int Rotate (int mu, int i){
        int mu_p = (mu+i)%4;
        if( mu >= 4)
                mu_p += 4;
//      printf("Rotate(%d, %d)=%d)\n",mu,i,mu_p);
        return mu_p;
}
static int NotParallel( int mu, int nu){
        int dif = mu-nu;
        if (dif==0 ||dif==-4||dif==4)return 0;
        else return 1;
}

void PT::asqtad_fat(AsqDArg *asq_arg, matrix *fatlink){
  char *fname = "Smear()";

//--------------------------------------------------------------------
// c1 -> one link; c2 -> 3-link; c3 -> 3-link staple; c5 -> 5-link staple;
// c7 -> 7-link staple; c6 -> 5-link "straight" staple
//--------------------------------------------------------------------

  IFloat c1 = asq_arg->c1;
  IFloat c3 = asq_arg->c3;
  IFloat c5 = asq_arg->c5;
  IFloat c7 = asq_arg->c7;
  IFloat c6 = asq_arg->c6;
#if 0
  IFloat c1 = GJP.KS_coeff();
  IFloat c2 = GJP.Naik_coeff();
  IFloat c3 = -GJP.staple3_coeff();
  IFloat c5 = GJP.staple5_coeff();
  IFloat c7 = -GJP.staple7_coeff();
  IFloat c6 = GJP.Lepage_coeff();
#endif

  int i, j; 
  const int N = 4;
//  int vol = GJP.VolNodeSites();
//  ParTransAsqtad pt(*this);
  matrix *result[NUM_DIR];
  matrix *Unit;
  matrix *P3[NUM_DIR];
  matrix *P3mu[NUM_DIR];
  matrix *P5[NUM_DIR];
  matrix *P5mu[NUM_DIR];
  matrix *P7[NUM_DIR];
  matrix *P7mu[NUM_DIR];
  matrix *P6[NUM_DIR];
  matrix *P6mu[NUM_DIR];
  matrix *Pmumu[NUM_DIR];
  matrix *Pmumumu[NUM_DIR];
//VRB.Flow(cname,fname,"vol=%d\n",vol);
  Unit = (matrix *)Alloc(cname,fname,"Unit",sizeof(matrix)*vol);
  for(i = 0;i<N;i++)
    Pmumumu[i] = P7mu[i] = (matrix *)Alloc(cname,fname,"Pmumumu[i]",sizeof(matrix)*vol);
  for(i = 0;i<N;i++)
    Pmumu[i] = P7[i] = (matrix *)Alloc(cname,fname,"Pmumu[i]",sizeof(matrix)*vol);
  for(i = 0;i<N;i++)
    P6mu[i] = P5mu[i] = (matrix *)Alloc(cname,fname,"P6mu[i]",sizeof(matrix)*vol);
  for(i = 0;i<N;i++)
    P6[i] = P5[i] = (matrix *)Alloc(cname,fname,"P6[i]",sizeof(matrix)*vol);
  for(i = 0;i<N;i++)
    P3mu[i] = (matrix *)Alloc(cname,fname,"P3mu[i]",sizeof(matrix)*vol);
  for(i = 0;i<N;i++)
    P3[i] = (matrix *)Alloc(cname,fname,"P3[i]",sizeof(matrix)*vol);
  for(j = 0;j<POS_DIR;j++){
//     result[j] = fields[0] + vol*j;
     result[j] = fatlink + vol*j;
//    result[j+POS_DIR] = fields[1] + vol*j;
//     VRB.Flow(cname,fname,"result[%d]=%p\n",j,result[j]);
  }
//     result[j] = (matrix *)fmalloc(sizeof(matrix)*vol);

  for(j = 0;j<vol;j++) Unit[j] = c1;
  for(j = 0;j<POS_DIR;j++)
#if 1
    bzero((char *)result[j],vol*18*sizeof(Float));
#else
     for(int k = 0;k<vol;k++) (result[j]+k)->Zeromatrix();
#endif
//  Float dtime = -dclock();
  int nflops = 0;
//  ParTrans::PTflops=0;
  int dir[] = {6,0,2,4,7,1,3,5},dirs[8]; //mapping between ParTrans and DiracOpAsqtad
  matrix *min[NUM_DIR],*mout[NUM_DIR];
  if (NUM_DIR%N !=0) fprintf(stderr,"%s:NUM_DIR(%d)is not divisible by N(%d)\n",fname,NUM_DIR,N);
  for(int mu = 0;mu<POS_DIR;mu += N){
	int mu_p,nu_p,rho_p,sigma_p;
    for(i  = 0;i<N;i++){
      mu_p = Rotate(mu,i);
      min[i] = Unit;
      mout[i] = result[mu_p];
      dirs[i] = dir[mu_p];
    }
    mat(N,mout,min,dirs);
    for(int nu = 0;nu<NUM_DIR;nu++)
    if(NotParallel(mu,nu)){
      for(i  = 0;i<N;i++){
	nu_p = Rotate(nu,i);
        min[i] = Unit;
        mout[i] = P3[i];
        dirs[i] = dir[nu_p];
      }
      mat(N,mout,min,dirs);
      for(i  = 0;i<N;i++){
		mu_p = Rotate(mu,i);
        min[i] = P3[i];
        mout[i] = P3mu[i];
        dirs[i] = dir[mu_p];
      }
      mat(N,mout,min,dirs);
      for(int rho = 0;rho<NUM_DIR;rho++)
      if(NotParallel(mu,rho) && NotParallel(nu,rho)){
        for(i  = 0;i<N;i++){
  	   	  rho_p = Rotate(rho,i);
          min[i] = P3[i];
          mout[i] = P5[i];
          dirs[i] = dir[rho_p];
        }
        mat(N,mout,min,dirs);
        for(i  = 0;i<N;i++){
	  	mu_p = Rotate(mu,i);
          min[i] = P5[i];
          mout[i] = P5mu[i];
          dirs[i] = dir[mu_p];
        }
        mat(N,mout,min,dirs);
        for(int sigma = 0;sigma<NUM_DIR;sigma++)
        if(NotParallel(mu,sigma) && NotParallel(nu,sigma)&&NotParallel(rho,sigma)){
          for(i  = 0;i<N;i++){
    		sigma_p = Rotate(sigma,i);
            min[i] = P5[i];
            mout[i] = P7[i];
            dirs[i] = dir[sigma_p];
          }
          mat(N,mout,min,dirs);
          for(i  = 0;i<N;i++){
    		mu_p = Rotate(mu,i);
            min[i] = P7[i];
            mout[i] = P7mu[i];
            dirs[i] = dir[mu_p];
          }
          mat(N,mout,min,dirs);
          for(i  = 0;i<N;i++){
    		sigma_p = Rotate(sigma,i);
            int sig_n = (sigma_p+4)%8;
            min[i] = P7mu[i];
            mout[i] = P7[i];
            dirs[i] = dir[sig_n];
          }
          mat(N,mout,min,dirs);
		  Float c75 = c7/c5;
          for(i  = 0;i<N;i++)
            vaxpy3_m(P5mu[i],&c75,P7[i],P5mu[i],vol*3);
	  nflops +=vol*18*2*N;
        }
        for(i  = 0;i<N;i++){
            rho_p = Rotate(rho,i);
            int rho_n = (rho_p+4)%8;
            min[i] = P5mu[i];
            mout[i] = P5[i];
            dirs[i] = dir[rho_n];
        }
        mat(N,mout,min,dirs);
		Float c53 = c5/c3;
        for(i  = 0;i<N;i++)
          vaxpy3_m(P3mu[i],&c53,P5[i],P3mu[i],vol*3);
		nflops +=vol*18*2*N;
      }
      for(i  = 0;i<N;i++){
        nu_p = Rotate(nu,i);
        min[i] = P3[i];
        mout[i] = P6[i];
        dirs[i] = dir[nu_p];
      }
      mat(N,mout,min,dirs);
      for(i  = 0;i<N;i++){
        mu_p = Rotate(mu,i);
        min[i] = P6[i];
        mout[i] = P6mu[i];
        dirs[i] = dir[mu_p];
      }
      mat(N,mout,min,dirs);
      for(i  = 0;i<N;i++){
        nu_p = Rotate(nu,i);
        int nu_n = (nu_p+4)%8;
        min[i] = P6mu[i];
        mout[i] = P6[i];
        dirs[i] = dir[nu_n];
      }
      mat(N,mout,min,dirs);
	  Float c63 = c6/c3;
      for(i  = 0;i<N;i++)
        vaxpy3_m(P3mu[i],&c63,P6[i],P3mu[i],vol*3);
	  nflops +=vol*18*2*N;
      for(i  = 0;i<N;i++){
        nu_p = Rotate(nu,i);
        int nu_n = (nu_p+4)%8;
        min[i] = P3mu[i];
        mout[i] = P3[i];
        dirs[i] = dir[nu_n];
      }
      mat(N,mout,min,dirs);

	  Float c31 = c3/c1;
	  for(i  = 0;i<N;i++){
	    mu_p =Rotate(mu,i);
        vaxpy3_m(result[mu_p],&c31,P3[i],result[mu_p],vol*3);
      }
	  nflops +=vol*18*2*N;
    }
  }

  Free( Unit);
  for(j = 0;j<N;j++){
    Free( P3[j]);
    Free( P3mu[j]);
    Free( P5[j]);
    Free( P5mu[j]);
    Free( P7[j]);
    Free( P7mu[j]);
  }
}
  
void PT::asqtad_long(AsqDArg *asq_arg, matrix *longlink, matrix *longlink_m){

  IFloat c2 = asq_arg->c2;
  int i, j; 
  const int N = 4;

  matrix *result[NUM_DIR];
  matrix *Unit;
  matrix *Pmumu[NUM_DIR];
  int dir[] = {6,0,2,4,7,1,3,5},dirs[8]; //mapping between ParTrans and DiracOpAsqtad
  matrix *min[NUM_DIR],*mout[NUM_DIR];

  for(j = 0;j<POS_DIR;j++){
     result[j] = longlink + vol*j;
//     VRB.Flow(cname,fname,"result[%d]=%p\n",j,result[j]);
  }
  if (longlink_m != NULL)
  for(j = 0;j<POS_DIR;j++){
     result[j+POS_DIR] = longlink_m + vol*j;
  }
  for(j = 0;j<vol;j++) Unit[j] = c2;
  int n_dirs = POS_DIR;
  if (longlink_m != NULL) n_dirs = NUM_DIR;
  for(int mu = 0;mu<n_dirs;mu += N){
    if (mu == 4) for(j = 0;j<vol;j++) Unit[j] = -c2;
    int mu_p;
    for(i  = 0;i<N;i++){
      mu_p = Rotate(mu,i);
      min[i] =  Unit;
      mout[i] = result[mu_p];
      dirs[i] = dir[mu_p];
    }
    mat(N,mout,min,dirs);
    for(i  = 0;i<N;i++){
      mu_p = Rotate(mu,i);
      min[i] =  result[mu_p];
      mout[i] = Pmumu[i];
      dirs[i] = dir[mu_p];
    }
    mat(N,mout,min,dirs);
    for(i  = 0;i<N;i++){
      mu_p = Rotate(mu,i);
      min[i] = Pmumu[i];
      mout[i] =  result[mu_p];
      dirs[i] = dir[mu_p];
    }
    mat(N,mout,min,dirs);
  }

  Free( Unit);
  for(j = 0;j<N;j++){
    Free( Pmumu[j]);
  }

//  dtime +=dclock();
//  nflops += ParTrans::PTflops;
//  printf("%s:%s:",cname,fname);
//  print_flops(nflops,dtime);

#if 0
  printf("%s:%s:",cname,fname);
  for(i = 0;i<3;i++){
  for(j = 0;j<POS_DIR;j++){
	IFloat *tmp = (IFloat *)(fields[i]+vol*j);
  printf("%e ",*tmp);
  }
  printf("\n");
  }
#endif

}
#endif
