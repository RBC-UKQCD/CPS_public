/* -*- mode:c++; c-basic-offset:4 -*- */
#ifndef BFM_HDCG_H
#define BFM_HDCG_H

/* EigCG code for Mpc^dagger Mpc solver 2012 by Qi
 * copied from the CG code (2012 ?version in bfm) and made modifications and 
 * adding the eigCG part
 */

#include <bfm.h>
#include <bfm_qdp.h>

#include <util/lattice/bfm_evo.h>
#include <util/lattice/hdcg_controller.h>

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <cstring>
#include <cassert>

#ifdef LIST_ENGINE
#include <pprefetch.h>
#endif

#include<qmp.h>

//template<class Float> void matrix_dgemm (const int M,const int N, const int K, Float **A, const double *B, double *C);
//void min_eig_index(int *INDEX, int nev,double *EIG, int n);
//void invert_H_matrix(complex<double> *data, int n); //if n is large enough, must be parallerized!!!
//template<class Float> void eigcg_vec_mult(Float* V, const int m, double *QZ, const int n, const int f_size_cb, const int nthread, const int me);

template<class Float, class Float_h>
int bfm_evo<Float>::HD_CGNE_M(BfmMultiGrid<Float_h> &hdcg, Fermion_t solution[2], Fermion_t source[2])
{
    int me = this->thread_barrier();
    Fermion_t src = this->threadedAllocFermion(); 
    Fermion_t tmp = this->threadedAllocFermion(); 
    Fermion_t Mtmp= this->threadedAllocFermion(); 

    // src_o = Mdag * (source_o - Moe MeeInv source_e)
    this->MooeeInv(source[Even],tmp,DaggerNo);
    this->Meo(tmp,src,Odd,DaggerNo);
    this->axpy(tmp,src,source[Odd],-1.0);
    this->Mprec(tmp,src,Mtmp,DaggerYes);  
  
    if (isBoss() && !me) printf("hdcg.Pcg(solution[%d](%p),src(%p),tmp(%p)\n");
    int iter = hdcg.Pcg(solution[Odd],src,tmp,1.0e-6,5.0e-3);
//    int iter = this->HD_CGNE_prec(solution[Odd], src);

    // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
    this->Meo(solution[Odd],tmp,Even,DaggerNo);
    this->axpy(src,tmp,source[Even],-1.0);
    this->MooeeInv(src,solution[Even],DaggerNo);
  
    this->threadedFreeFermion(tmp);
    this->threadedFreeFermion(src);
    this->threadedFreeFermion(Mtmp);

    return iter;
}

#if 0
template <class Float>
void myaxpy(Float *r, Float *x, Float *y, double a, int len)
{
  double aa_b[4+4]; // 64 bytes
  double *aa = ALIGNIT(aa_b);
  for(int i=0;i<4;i++) aa[i] = a;

  if ( sizeof(Float) == sizeof(double) ) {
    vmx_vaxpy((double *)r,(double *)aa,(double *)x,(double *)y,len);
  } else { 
    vmx_vaxpy_s((float *)r,(double *)aa,(float *)x,(float *)y,len);
  }

  // thread_barrier();

  return;
}

template<class Float>
void eigcg_vec_mult3(Float* V, const int m, double *QZ, const int n, const int f_size_cb,
                     const int nthread, const int me,
                     bfm_evo<Float> &bfmobj)
//QZ is saved in column major format. 
//perform V = V*QZ;
{
    std::vector<Float *> ret(n, NULL);
    for(int i = 0; i < n; ++i) {
        ret[i] = (Float *)bfmobj.threadedAllocFermion();
        bfmobj.set_zero(ret[i]);
    }

    assert(nthread >= n);
    int mywork, myoff;
    bfmobj.thread_work_partial_nobarrier(f_size_cb / 24, me / n, nthread / n, mywork, myoff);
    myoff *= 24;

    if(me < nthread / n * n) {
        int i = me % n;
        for(int j = 0; j < m; ++j) {
            myaxpy(ret[i] + myoff, // &r
                   V + j * f_size_cb + myoff, // &x
                   ret[i] + myoff, // &y
                   QZ[i * m + j], // a
                   mywork); // length
        }
    }

    for(int i = 0; i < n; ++i) {
        bfmobj.copy((Fermion_t)(V + i * f_size_cb), ret[i]);
        bfmobj.threadedFreeFermion(ret[i]);
    }
}
#endif

#endif
