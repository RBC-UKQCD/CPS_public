/* -*- mode:c++; c-basic-offset:2 -*- */
/****************************************************************************/
/* Aug 2012                                                                 */
/* Hantao Yin                                                               */
/*                                                                          */
/* bfm_evo.h                                                                */
/* functions added by Hantao (mainly DWF-like fermion evolution related).   */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_BFM_EVO_HT_H
#define INCLUDED_BFM_EVO_HT_H

#include <stdio.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <omp.h>
#include <math.h>
#include <vector>
#include "bfm_evo_aux.h"

// FIXME: it inherits from bfm_qdp for the sole reason of using its
// importGauge() function. I'm too lazy to do any manual
// shifts/conjugates ......
template <class Float>
class bfm_evo : public bfm_qdp<Float> {
public:
  // BFM has this
  // enum {Even = 0, Odd};

  integer cps_idx_cb(int x[4], int s, int reim, int i, int i_size);

  // s outer most
  integer cps_idx(int x[4], int s, int reim, int i, int i_size);

  // s inner most (but outside color and spin)
  integer cps_idx_s(int x[4], int s, int reim, int i, int i_size);

  // compute the vector pair (v1, v2) needed to calculate fermion force.
  void calcMDForceVecs(Fermion_t v1[2], Fermion_t v2[2],
                       Fermion_t phi1, Fermion_t phi2);

  void Booee(Fermion_t psi, Fermion_t chi, int dag);

  // Ritz method used to compute the maximum/minimum eigenvalue of M^\dag M.
  // Use algorithm presented in arXiv: hep-lat/9507023.
  //
  // If compute_min == true then we compute the minmum eigenvalue of
  // M^\dag M, otherwise we compute the maximum eigenvalue, i.e. the
  // negative of the minimum eigenvalue of -M^\dag M.
  double ritz(Fermion_t x, int compute_min);

  // solve a propagator, for HtCayleyTanh this is just unpreconditioned CG
  // for HmCayleyTanh this is D^{-1} Dminus acting on in[2].
  // int prop_solve(Fermion_t out[2], Fermion_t in[2]);

  // ======================================================================
  // these functions need to be rewritten to fit into bfm style.
  // currently they contain both bfm and CPS style gauge/fermion fields.
private:
  // auxiliary functions used by compute_force().
  void fforce_site(Float *mom, Float *gauge,
                   Float *v1, Float *v1p,
                   Float *v2, Float *v2p, int mu, Float coef);

  void fforce_internal(Float *mom, Float *gauge,
                       Float *v1, Float *v2, // internal data
                       Float coef, int mu,
                       int me, int nthreads);

  void fforce_surface(Float *mom, Float *gauge,
                      Float *v1, Float *v2, // internal data
                      Float *v1_s, Float *v2_s, // surface data
                      Float coef, int mu);

  void copySendFrmData(Float v3d[], Float v4d[], int mu, bool send_neg);

  // complex version of axpy()
  void axpy_c(Fermion_t r, Fermion_t x, Fermion_t y, std::complex<double> a, Fermion_t tmp) {
    this->zaxpy(r, x, y, a);
  }
public:
  void thread_work_partial_nobarrier(int nwork, int me, int nthreads,
                                     int &mywork, int &myoff)
  {
    int basework = nwork / nthreads;
    int backfill = nthreads - (nwork % nthreads);
    mywork = (nwork + me) / nthreads;
    myoff  = basework * me;
    if ( me > backfill ) 
        myoff += (me-backfill);
  }

  // compute fermion force:
  //
  // mom += coef * (phiL^\dag e_i(M) \phiR + \phiR^\dag e_i(M^\dag) \phiL)
  //
  // For BFM M is M = M_oo - M_oe M^{-1}_ee M_eo
  void compute_force(Float *mom, Float *gauge, Fermion_t phiL, Fermion_t phiR, double coef);
  
  // psi assumes the following order: (color, spin, s, x, y, z, t),
  // mainly used to import/export the "v1" and "v2" vectors in evolution.
  template<typename FloatEXT>
  void thread_impexFermion_s(FloatEXT *psi, Fermion_t handle[2], int doimport);

  Float *threadedAllocFloat(size_t size, int mem_type=mem_slow);
  void threadedFreeFloat(Float *);

  // bicg_M: Biconjugate gradient method on preconditioned Dirac
  // operator (It never converges).
  //
  // FIXME: test code only, don't use it unless you know what you are
  // doing.
  int bicg_M(Fermion_t sol, Fermion_t src);

  // bicgstab_M: Biconjugate gradient stabilized method on
  // preconditioned Dirac operator.
  //
  // FIXME: test code only, don't use it unless you know what you are
  // doing.
  int bicgstab_M(Fermion_t sol, Fermion_t src);

  // GCR, solves M x = b
  int gcr_M(Fermion_t sol, Fermion_t src);

  // GMRES(m) solves M x = b.
  //
  // Restarts after m iterations.
  int gmres_M(Fermion_t sol, Fermion_t src, const int m);
public:
  //======================================================================
  // the following member functions are single-threaded functions:
  // ======================================================================

  // psi assumes 5D even/odd preconditioned order: (color, spin, x, y, z, t, s)/2
  template<typename FloatEXT>
  void cps_impexcbFermion(FloatEXT *psi, Fermion_t handle, int doimport, int cb);

  // psi assumes regular canonical order: (color, spin, x, y, z, t, s)
  template<typename FloatEXT>
  void cps_impexFermion(FloatEXT *psi, Fermion_t handle[2], int doimport);

  // psi assumes the following order: (color, spin, s, x, y, z, t),
  // mainly used to import/export the "v1" and "v2" vectors in evolution.
  template<typename FloatEXT>
  void cps_impexFermion_s(FloatEXT *psi, Fermion_t handle[2], int doimport);

  template<typename FloatEXT>
  void cps_importGauge(FloatEXT *importme);

  //EigCG
  Fermion_t allocCompactFermion   (int mem_type=mem_slow);
  Fermion_t threadedAllocCompactFermion   (int mem_type=mem_slow);
  void* threaded_alloc(int length, int mem_type=mem_slow);
  void threaded_free(void *handle);
  int EIG_CGNE_M(Fermion_t solution[2], Fermion_t source[2]);
  int Eig_CGNE_prec(Fermion_t psi, Fermion_t src);

  // copied from Jianglei's bfm
  double CompactMprec(Fermion_t compact_psi,
                      Fermion_t compact_chi,
                      Fermion_t psi,
                      Fermion_t chi,
                      Fermion_t tmp,
                      int dag,int donrm=0) ;

  // copied from Jianglei's bfm
  void CompactMunprec(Fermion_t compact_psi[2],
                      Fermion_t compact_chi[2],
                      Fermion_t psi[2],
                      Fermion_t chi[2],
                      Fermion_t tmp,
                      int dag);

  // do deflation using eigenvectors/eigenvalues from Rudy's Lanczos code.
  void deflate(Fermion_t out, Fermion_t in,
               const multi1d<Fermion_t [2]> *evec,
               const multi1d<Float> *eval, int N);
};

template<class Float>
integer bfm_evo<Float>::cps_idx_cb(int x[4], int s, int reim, int i, int i_size)
{
  // int cb   = ( x[0]+x[1]+x[2]+x[3] )&0x1;

  int csite
    =x[0] +this->node_latt[0]
    *(x[1] + this->node_latt[1]
      *(x[2] +this->node_latt[2]
        *(x[3] +s*this->node_latt[3])));
  csite /= 2;

  // int cbvol = (this->node_latt[0]*
  //              this->node_latt[1]*
  //              this->node_latt[2]*
  //              this->node_latt[3]*
  //              this->cbLs)/2;

  // return (cb*cbvol+csite)*i_size*2 + i*2 + reim;
  return csite*i_size*2 + i*2 + reim;
}

template<class Float>
integer bfm_evo<Float>::cps_idx(int x[4], int s, int reim, int i, int i_size)
{
  int csite =
    x[0] + this->node_latt[0]
    *(x[1] + this->node_latt[1]
      *(x[2] +this->node_latt[2]
        *(x[3] +s*this->node_latt[3])));

  return (csite*i_size + i)*2 + reim;
}

template<class Float>
integer bfm_evo<Float>::cps_idx_s(int x[4], int s, int reim, int i, int i_size)
{
  int csite =
    s + this->Ls
    *(x[0] + this->node_latt[0]
      *(x[1] + this->node_latt[1]
        *(x[2] +this->node_latt[2]
          *x[3])));

  return (csite*i_size + i)*2 + reim;
}

template <class Float> template<typename FloatEXT>
void bfm_evo<Float>::cps_impexcbFermion(FloatEXT *psi, Fermion_t handle, int doimport, int cb)
{
  int Nspinco=12;
  int i_inc = this->simd() * 2;
  int vol5d =
    this->node_latt[0] *
    this->node_latt[1] *
    this->node_latt[2] *
    this->node_latt[3] *
    this->Ls;

  Float *bagel = (Float *)handle;
  omp_set_num_threads(this->nthread);

#pragma omp parallel for 
  for (int site = 0; site < vol5d; site++) {
    
    int x[4], s;
    int si=site;
    x[0]=si%this->node_latt[0];    si=si/this->node_latt[0];
    x[1]=si%this->node_latt[1];    si=si/this->node_latt[1];
    x[2]=si%this->node_latt[2];    si=si/this->node_latt[2];
    x[3]=si%this->node_latt[3];
    s   =si/this->node_latt[3];
    
    int sp = this->precon_5d ? s : 0;
    if ( (x[0]+x[1]+x[2]+x[3] + sp &0x1) == cb ) {

      int bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1);
      int cidx_base = this->cps_idx_cb(x, s, 0, 0, Nspinco);

      for ( int co=0;co<Nspinco;co++ ) { 
        for ( int reim=0;reim<2;reim++ ) {
          // int bidx = bagel_idx(x, reim, co + Nspinco * (s / 2), Nspinco * this->cbLs, 1);
          // int bidx = this->bagel_idx5d(x, s, reim, co, Nspinco, 1);
          // int cidx = cps_idx_cb(x, s, reim, co, Nspinco);
          int bidx = bidx_base + reim + co * i_inc;
          int cidx = cidx_base + reim + co * 2;

          if ( doimport ) bagel[bidx] = psi[cidx];
          else psi[cidx] = bagel[bidx] ;
        }}//co,reim
    }//cb
  }//xyzts
}

template <class Float> template<typename FloatEXT>
void bfm_evo<Float>::cps_impexFermion(FloatEXT *psi, Fermion_t handle[2], int doimport)
{
  int Nspinco=12;
  int i_inc = this->simd() * 2;
  int vol5d =
    this->node_latt[0] *
    this->node_latt[1] *
    this->node_latt[2] *
    this->node_latt[3] *
    this->Ls;

  omp_set_num_threads(this->nthread);
  Float *bagel[2] = { (Float *)handle[0], (Float *)handle[1] };

#pragma omp parallel for 
  for (int site = 0; site < vol5d; site++) {
    int x[4], s;
    int si=site;
    x[0]=si%this->node_latt[0];    si=si/this->node_latt[0];
    x[1]=si%this->node_latt[1];    si=si/this->node_latt[1];
    x[2]=si%this->node_latt[2];    si=si/this->node_latt[2];
    x[3]=si%this->node_latt[3];
    s   =si/this->node_latt[3];
    
    int sp = this->precon_5d ? s : 0;
    int cb = x[0]+x[1]+x[2]+x[3]+sp &0x1;

    int bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1);
    int cidx_base = this->cps_idx(x, s, 0, 0, Nspinco);

    for ( int co=0;co<Nspinco;co++ ) { 
      for ( int reim=0;reim<2;reim++ ) { 
        // int bidx = bagel_idx(x, reim, co + Nspinco * (s / 2), Nspinco * this->cbLs, 1);
        // int bidx = this->bagel_idx5d(x, s, reim, co, Nspinco, 1);
        // int cidx = cps_idx(x, s, reim, co, Nspinco);
        int bidx = bidx_base + reim + co * i_inc;
        int cidx = cidx_base + reim + co * 2;

        if ( doimport ) bagel[cb][bidx] = psi[cidx];
        else psi[cidx] = bagel[cb][bidx];
      }}//co, reim
  }//xyzts
}

template <class Float> template<typename FloatEXT>
void bfm_evo<Float>::cps_impexFermion_s(FloatEXT *psi, Fermion_t handle[2], int doimport)
{
  int Nspinco=12;
  int i_inc = this->simd() * 2;
  int vol5d =
    this->node_latt[0] *
    this->node_latt[1] *
    this->node_latt[2] *
    this->node_latt[3] *
    this->Ls;

  omp_set_num_threads(this->nthread);
  Float *bagel[2] = { (Float *)handle[0], (Float *)handle[1] };

#pragma omp parallel for 
  for (int site = 0; site < vol5d; site++) {
    int x[4], s;
    int si=site;
    s   =si%this->Ls;              si=si/this->Ls;
    x[0]=si%this->node_latt[0];    si=si/this->node_latt[0];
    x[1]=si%this->node_latt[1];    si=si/this->node_latt[1];
    x[2]=si%this->node_latt[2];
    x[3]=si/this->node_latt[2];

    int sp = this->precon_5d ? s : 0;
    int cb = x[0]+x[1]+x[2]+x[3]+sp & 0x1;

    int bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1);
    int cidx_base = this->cps_idx_s(x, s, 0, 0, Nspinco);

    for ( int co=0;co<Nspinco;co++ ) {
      for ( int reim=0;reim<2;reim++ ) {
        // int bidx = this->bagel_idx5d(x, s, reim, co, Nspinco, 1);
        // int cidx = cps_idx_s(x, s, reim, co, Nspinco);
        int bidx = bidx_base + reim + co * i_inc;
        int cidx = cidx_base + reim + co * 2;

        if ( doimport ) bagel[cb][bidx] = psi[cidx];
        else psi[cidx] = bagel[cb][bidx];
      }}//co, reim
  }//xyzts
}

template <class Float> template<typename FloatEXT>
void bfm_evo<Float>::thread_impexFermion_s(FloatEXT *psi, Fermion_t handle[2], int doimport)
{
  int Nspinco=12;
  int i_inc = this->simd() * 2;
  int vol5d =
    this->node_latt[0] *
    this->node_latt[1] *
    this->node_latt[2] *
    this->node_latt[3] *
    this->Ls;

  int me, thrlen, throff;
  this->thread_work(vol5d, me, thrlen, throff);

  Float *bagel[2] = { (Float *)handle[0], (Float *)handle[1] };

  for (int site = 0; site < thrlen; ++site) {
    int x[4], s;
    int si=site + throff;
    s   =si%this->Ls;              si=si/this->Ls;
    x[0]=si%this->node_latt[0];    si=si/this->node_latt[0];
    x[1]=si%this->node_latt[1];    si=si/this->node_latt[1];
    x[2]=si%this->node_latt[2];
    x[3]=si/this->node_latt[2];

    int sp = this->precon_5d ? s : 0;
    int cb = x[0]+x[1]+x[2]+x[3]+sp & 0x1;

    int bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1);
    int cidx_base = this->cps_idx_s(x, s, 0, 0, Nspinco);

    for ( int co=0;co<Nspinco;co++ ) {
      for ( int reim=0;reim<2;reim++ ) {
        // int bidx = this->bagel_idx5d(x, s, reim, co, Nspinco, 1);
        // int cidx = cps_idx_s(x, s, reim, co, Nspinco);
        int bidx = bidx_base + reim + co * i_inc;
        int cidx = cidx_base + reim + co * 2;

        if ( doimport ) bagel[cb][bidx] = psi[cidx];
        else psi[cidx] = bagel[cb][bidx];
      }}//co, reim
  }//xyzts
}

template <class Float>
Float * bfm_evo<Float>::threadedAllocFloat(size_t size, int mem_type)
{
  int me = this->thread_barrier();

  void *ret;
  if ( me == 0 ) {
    ret = bfm_alloc(size * sizeof(Float), mem_type);
  }
  ret = this->thread_bcast(me, ret);
  this->thread_barrier();
  return (Float *)ret;
}

template <class Float>
void bfm_evo<Float>::threadedFreeFloat(Float *f)
{
  int me = this->thread_barrier();
  if ( me == 0 ) { 
    bfm_free(f);
  }
  this->thread_barrier();
}

static inline int idx_4d(const int x[4], const int lx[4]) {
  int ret = 0;
  for(int i = 3; i >= 0; --i) {
    ret = ret * lx[i] + x[i];
  }
  return ret;
}

static inline int idx_5d(const int x[5], const int lx[5]) {
  int ret = 0;
  for(int i = 4; i >= 0; --i) {
    ret = ret * lx[i] + x[i];
  }
  return ret;
}

static inline int idx_4d_surf(const int x[4], const int lx[4], int mu) {
  int ret = 0;
  for(int i = 3; i >= 0; --i) {
    if(i == mu) continue;
    ret = ret * lx[i] + x[i];
  }
  return ret;
}

static inline int idx_5d_surf(const int x[5], const int lx[5], int mu) {
  int ret = 0;
  for(int i = 4; i >= 0; --i) {
    if(i == mu) continue;
    ret = ret * lx[i] + x[i];
  }
  return ret;
}

template <class Float> template<typename FloatEXT>
void bfm_evo<Float>::cps_importGauge(FloatEXT *importme)
{
  multi1d<LatticeColorMatrix> U(Nd);

  omp_set_num_threads(this->nthread);

  int Ndircoco=72;
  int Ncoco = 9;
  QDPdouble *U_p;

  int vol4d =
    this->node_latt[0] *
    this->node_latt[1] *
    this->node_latt[2] *
    this->node_latt[3];

  for (int mu=0;mu<Nd;mu++) {
    U_p = (QDPdouble *)&(U[mu].elem(0).elem());

#pragma omp parallel for 
    for (int site=0;site<vol4d;site++ ) {
      int x[4];
      int s=site;
      x[0]=s%this->node_latt[0];    s=s/this->node_latt[0];
      x[1]=s%this->node_latt[1];    s=s/this->node_latt[1];
      x[2]=s%this->node_latt[2];    s=s/this->node_latt[2];
      x[3]=s%this->node_latt[3];
      
      int qidx_base = this->chroma_idx(x, 0, 0, Ncoco);

      for(int coco = 0; coco < Ncoco; ++coco) {
        for ( int reim = 0; reim < 2; ++reim) {
          // int qidx = this->chroma_idx(x,reim,coco,Ncoco);
          int qidx = qidx_base + reim + coco * 2;

          int siteoff = mu + Nd * site;
          int cidx = reim + 2 * (coco + Ncoco * siteoff);
          U_p[qidx] = importme[cidx];
        }} // reim,coco
    } // x
  }//mu

  // to bfm
  this->importGauge(U);
}

template <class Float>
void bfm_evo<Float>::calcMDForceVecs(Fermion_t v1[2], Fermion_t v2[2],
                       Fermion_t phi1, Fermion_t phi2)
{
  // Meo is Wilson D times a matrix (see page 27 in Peter's draft).
  // Moe/Meo: check bfmbase<Float>::G5D_Meo() in bfmdperp.C.
  // Mee/Moo: check bfmbase<Float>::G5D_Mooee().
  // Mee/Moo inverse: check bfmbase<Float>::G5D_MooeeInv().
  
  // v2e
  this->Meo(phi2, v1[Odd], Even, DaggerNo);
  this->MooeeInv(v1[Odd], v1[Even], DaggerNo);
  this->Booee(v1[Even], v2[Even], DaggerNo);

  // v2o
  this->Booee(phi2, v2[Odd], DaggerNo);

  // v1e
  this->Meo(phi1, v1[Odd], Even, DaggerYes);
  this->MooeeInv(v1[Odd], v1[Even], DaggerYes);

  // v1o
  this->copy(v1[Odd], phi1);
}

template <class Float>
void bfm_evo<Float>::Booee(Fermion_t psi, Fermion_t chi, int dag)
{
  int Pminus=-1;
  int Pplus=1;

  // just copied the relevant part in G5D_Meo() over.
  if ( (this->solver == HmCayleyTanh)
       || (this->solver == HtCayleyTanh)
       || (this->solver == HwCayleyTanh)
       || (this->solver == HwCayleyZolo)
       || (this->solver == HtCayleyZolo)
       ) {

    if ( dag ) { 

      // Assemble the 5d matrix
      for(int s=0;s<this->Ls;s++){
        if ( s==0 ) {
          axpby_ssp_proj(chi,this->beo[s],psi,   -this->ceo[s+1]  ,psi,s,s+1,Pplus);
          axpby_ssp_proj(chi,   1.0,chi,this->mass*this->ceo[this->Ls-1],psi,s,this->Ls-1,Pminus);
        } else if ( s==(this->Ls-1)) { 
          axpby_ssp_proj(chi,this->beo[s],psi,this->mass*this->ceo[0],psi,s,0,Pplus);
          axpby_ssp_proj(chi,1.0,chi,-this->ceo[s-1],psi,s,s-1,Pminus);
        } else {
          axpby_ssp_proj(chi,this->beo[s],psi,-this->ceo[s+1],psi,s,s+1,Pplus);
          axpby_ssp_proj(chi,1.0   ,chi,-this->ceo[s-1],psi,s,s-1,Pminus);
        }
      }

    } else { 

      // Assemble the 5d matrix
      for(int s=0;s<this->Ls;s++){
        if ( s==0 ) {
          //	chi = bs psi[s] + cs[s] psi[s+1}
          //    chi += -mass*cs[s] psi[s+1}
          axpby_ssp_proj(chi,this->beo[s],psi,-this->ceo[s],psi ,s, s+1,Pminus);
          axpby_ssp_proj(chi,1.0,chi,this->mass*this->ceo[s],psi,s,this->Ls-1,Pplus);
        } else if ( s==(this->Ls-1)) { 
          axpby_ssp_proj(chi,this->beo[s],psi,this->mass*this->ceo[s],psi,s,0,Pminus);
          axpby_ssp_proj(chi,1.0,chi,-this->ceo[s],psi,s,s-1,Pplus);
        } else {
          axpby_ssp_proj(chi,this->beo[s],psi,-this->ceo[s],psi,s,s+1,Pminus);
          axpby_ssp_proj(chi,1.0,chi,-this->ceo[s],psi,s,s-1,Pplus);
        }
      }
    }
  } else if(this->solver == DWF && this->precon_5d == 1) {
    // Booee is the identity matrix in this case.
    this->copy(chi, psi);
    return;
  } else {
    if ( this->isBoss() ) {
      printf("Booee: preconditioned unimplemented \n");
    }
    exit(-1);
  }
}

static inline double quad_solve(double *ct, double *st,
                                double a, double b, double c,
                                double d, double e, double f)
{
  double p = b * (d - f) + e * (c - a);
  double q = b * (d + f) - e * (c + a);
  double r = 2 * (c * d - a * f);

  // solve p + q * cos(2t) + r * sin(2t) = 0
  double den = sqrt(q * q + r * r);
  double ca = q / den;

  double ci = sqrt(0.5 * (1 + ca));
  double si = sqrt(0.5 * (1 - ca));
  if(r < 0) si = -si;

  double cb = -p / den;
  if(fabs(cb) > 1.) {
    printf("Panic: cos(psi) > 1\n");
    exit(-1);
  }
  double cj = sqrt(0.5 * (1 + cb));
  double sj = sqrt(0.5 * (1 - cb));

  double ct1 = ci * cj + si * sj;
  double st1 = si * cj - ci * sj;
  double v1 =
    (a * ct1 * ct1 + b * st1 * ct1 + c * st1 * st1)
    / (d * ct1 * ct1 + e * st1 * ct1 + f * st1 * st1);

  double ct2 = ci * cj - si * sj;
  double st2 = si * cj + ci * sj;
  double v2 =
    (a * ct2 * ct2 + b * st2 * ct2 + c * st2 * st2)
    / (d * ct2 * ct2 + e * st2 * ct2 + f * st2 * st2);

  if(v1 < v2) {
    *ct = ct1; *st = st1;
    return v1;
  } else {
    *ct = ct2; *st = st2;
    return v2;
  }
}

// Ritz method used to compute the maximum/minimum eigenvalue of M^\dag M.
// Use algorithm presented in arXiv: hep-lat/9507023.
template <class Float>
double bfm_evo<Float>::ritz(Fermion_t x, int compute_min)
{
  int me = this->thread_barrier();

  double stop_rsd = this->residual * this->residual;

  Fermion_t y = this->threadedAllocFermion();
  Fermion_t p = this->threadedAllocFermion();
  Fermion_t z = this->threadedAllocFermion();
  Fermion_t t = this->threadedAllocFermion();
  Fermion_t u = this->threadedAllocFermion();

  double mu, pnorm, gnorm2;

  // normalize x
  double fact = this->norm(x);
  fact = sqrt(1./ fact);
  this->scale(x, fact);

  if(this->isBoss() && !me) {
    printf("bfm_evo::ritz <x, x> = %17.10e\n", 1. / (fact * fact));
  }

  // y = A x, A = MdagM or -MdagM
  mu = this->Mprec(x, t, y, 0, 1);
  this->Mprec(t, y, u, 1);
  if(! compute_min) {
    this->scale(y, -1.);
    mu = -mu;
  }

  gnorm2 = this->axpy_norm(p, x, y, -mu);
  pnorm = sqrt(gnorm2);

  int i;
  for(i = 0; i < this->max_iter; ++i) {
    if(gnorm2 < stop_rsd) break;

    // if(i % 100 == 0 && this->isBoss() && !me) {
    //   printf("bfm_evo::ritz iter = %6d gnorm2 = %17.10e, mu = %17.10e\n", i, gnorm2, mu); 
    // }

    // z = A p
    double pap = this->Mprec(p, t, z, 0, 1);
    this->Mprec(t, z, u, 1);
    if(! compute_min) {
      this->scale(z, -1.);
      pap = -pap;
    }

    // minimize x cos(theta) + p / pnorm * sin(theta) via theta
    double d = this->norm(x);
    double e = 2. * this->inner_real(x, p) / pnorm;
    double f = 1.;
    // double a = this->inner_real(x, y);
    double a = mu * d;
    double b = 2. * this->inner_real(x, z) / pnorm;
    double c = pap / (pnorm * pnorm);

    double ct,st;
    mu = quad_solve(&ct, &st, a, b, c, d, e, f);

    this->axpby(x, x, p, ct, st / pnorm);
    this->axpby(y, y, z, ct, st / pnorm);

    double gnew = this->axpy_norm(t, x, y, -mu);
    double beta = ct * gnew / gnorm2;
    gnorm2 = gnew;

    // this->axpy(u, x, p, -st * pnorm); // ! not stable
    double xpp = this->inner_real(x, p);
    this->axpy(u, x, p, -xpp);

    pnorm = sqrt(this->axpy_norm(p, u, t, beta));
  }

  if(! compute_min) mu = -mu;

  // check eigenvalue again
  double xnorm = this->norm(x);
  double mux = this->Mprec(x, y, t, 0, 1);
  this->Mprec(y, t, u, 1);
  double mu_sq = this->norm(t);
  
  if(this->isBoss() && !me) {
    if(i < this->max_iter) {
      printf("bfm_evo::ritz converged at iteration %d.\n", i);
    } else {
      printf("bfm_evo::ritz maximum iteration number reached!\n");
    }
    printf("bfm_evo::ritz ||x|| = %17.10e\n", sqrt(xnorm));

    printf("bfm_evo::ritz three ways of computing the eigenvalue should agree.\n");
    printf("bfm_evo::ritz  eig1 = %17.10e\n", mu / xnorm);
    printf("bfm_evo::ritz  eig2 = %17.10e\n", mux / xnorm);
    printf("bfm_evo::ritz  eig3 = %17.10e\n", sqrt(mu_sq / xnorm));
  }

  this->threadedFreeFermion(y);
  this->threadedFreeFermion(p);
  this->threadedFreeFermion(z);
  this->threadedFreeFermion(t);
  this->threadedFreeFermion(u);
  return mu / xnorm;
}

// FIXME: I'll need to replace getPlusData by something else.
// For now it works.
#include <comms/scu.h>

template <class Float>
void bfm_evo<Float>::copySendFrmData(Float v3d[], Float v4d[], int mu, bool send_neg)
{
  int lclx[5] = {this->node_latt[0],
                 this->node_latt[1],
                 this->node_latt[2],
                 this->node_latt[3],
                 this->Ls};

  int low[4] = { 0, 0, 0, 0 };
  int high[4] = {lclx[0], lclx[1], lclx[2], lclx[3] };
  low[mu] = send_neg ? 0 : lclx[mu] - 1;
  high[mu] = low[mu] + 1;

  int block_size = 24 * lclx[4]; // s inner most

  const int hl[4] = {high[0] - low[0],
                     high[1] - low[1],
                     high[2] - low[2],
                     high[3] - low[3] };
  const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

  int me, thrlen, throff;
  this->thread_work(hl_sites, me, thrlen, throff);

  for(int i = 0; i < thrlen; ++i) {
    int x[4];
    int tmp = i + throff;
    x[0] = tmp % hl[0] + low[0]; tmp /= hl[0];
    x[1] = tmp % hl[1] + low[1]; tmp /= hl[1];
    x[2] = tmp % hl[2] + low[2]; tmp /= hl[2];
    x[3] = tmp % hl[3] + low[3];

    int off_4d = idx_4d(x, lclx);
    int off_3d = idx_4d_surf(x, lclx, mu);
        
    memcpy(v3d + off_3d * block_size,
           v4d + off_4d * block_size,
           sizeof(Float) * block_size);
  }
}

// Calculate fermion force on a specific site, also do the
// summation over s direction.
//
// FIXME: need to add a line sum in s direction to support splitting
// in s direction.
template<class Float>
void bfm_evo<Float>::fforce_site(Float *mom, Float *gauge,
                                 Float *v1, Float *v1p,
                                 Float *v2, Float *v2p, int mu, Float coef)
{
    Float t1[18], t2[18];

    switch(mu) {
    case 0:
      bfm_evo_aux::sprojTrXm(t1, v1p, v2, this->Ls, 0, 0);
      bfm_evo_aux::sprojTrXp(t2, v2p, v1, this->Ls, 0, 0);
      break;
    case 1:
      bfm_evo_aux::sprojTrYm(t1, v1p, v2, this->Ls, 0, 0);
      bfm_evo_aux::sprojTrYp(t2, v2p, v1, this->Ls, 0, 0);
      break;
    case 2:
      bfm_evo_aux::sprojTrZm(t1, v1p, v2, this->Ls, 0, 0);
      bfm_evo_aux::sprojTrZp(t2, v2p, v1, this->Ls, 0, 0);
      break;
    default:
      bfm_evo_aux::sprojTrTm(t1, v1p, v2, this->Ls, 0, 0);
      bfm_evo_aux::sprojTrTp(t2, v2p, v1, this->Ls, 0, 0);
    }

    bfm_evo_aux::su3_add(t1, t2);
    bfm_evo_aux::mDotMEqual(t2, gauge, t1);
    bfm_evo_aux::trless_am(t2, -coef);
    bfm_evo_aux::su3_add(mom, t2);
}

template<class Float>
void bfm_evo<Float>::fforce_internal(Float *mom, Float *gauge,
                                     Float *v1, Float *v2, // internal data
                                     Float coef, int mu,
                                     int me, int nthreads)
{
  int lclx[5] = {this->node_latt[0],
                 this->node_latt[1],
                 this->node_latt[2],
                 this->node_latt[3],
                 this->Ls};
  int low[4] = { 0, 0, 0, 0 };
  int high[4] = { lclx[0], lclx[1], lclx[2], lclx[3] };
  --high[mu];

  int block_size = 24 * lclx[4];

  const int hl[4] = {high[0] - low[0],
                     high[1] - low[1],
                     high[2] - low[2],
                     high[3] - low[3] };
  const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

  // note: some of the threads are dedicated to communication. There
  // must be exactly *nthreads* threads executing this function, the
  // variable *me* must range from 0 to nthreads - 1, inclusive.
  int thrlen, throff;
  this->thread_work_partial_nobarrier(hl_sites, me, nthreads,
                                      thrlen, throff);

  for(int i = 0; i < thrlen; ++i) {
    int x[4];
    int tmp = i + throff;
    x[0] = tmp % hl[0] + low[0]; tmp /= hl[0];
    x[1] = tmp % hl[1] + low[1]; tmp /= hl[1];
    x[2] = tmp % hl[2] + low[2]; tmp /= hl[2];
    x[3] = tmp % hl[3] + low[3];

    int off_4d = idx_4d(x, lclx);
    int gid = mu + 4 * off_4d;
    int fid = block_size * off_4d;

    int y[4] = {x[0], x[1], x[2], x[3]};
    ++y[mu];
    int fidp = block_size * idx_4d(y, lclx);

    this->fforce_site(mom + 18 * gid, gauge + 18 * gid,
                      v2 + fid, v2 + fidp,
                      v1 + fid, v1 + fidp, mu, coef);
  }
}

template<class Float>
void bfm_evo<Float>::fforce_surface(Float *mom, Float *gauge,
                                    Float *v1, Float *v2, // internal data
                                    Float *v1_s, Float *v2_s, // surface data
                                    Float coef, int mu)
{
  int lclx[5] = {this->node_latt[0],
                 this->node_latt[1],
                 this->node_latt[2],
                 this->node_latt[3],
                 this->Ls};
  int low[4] = { 0, 0, 0, 0 };
  int high[4] = { lclx[0], lclx[1], lclx[2], lclx[3] };
  low[mu] = lclx[mu] - 1;

  int block_size = 24 * lclx[4];

  const int hl[4] = {high[0] - low[0],
                     high[1] - low[1],
                     high[2] - low[2],
                     high[3] - low[3] };
  int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

  int me, thrlen, throff;
  this->thread_work(hl_sites, me, thrlen, throff);

  for(int i = 0; i < thrlen; ++i) {
    int x[4];
    int tmp = i + throff;
    x[0] = tmp % hl[0] + low[0]; tmp /= hl[0];
    x[1] = tmp % hl[1] + low[1]; tmp /= hl[1];
    x[2] = tmp % hl[2] + low[2]; tmp /= hl[2];
    x[3] = tmp % hl[3] + low[3];

    int off_4d = idx_4d(x, lclx);
    int gid = mu + 4 * off_4d;
    int fid = block_size * off_4d;
    int fid_s = block_size * idx_4d_surf(x, lclx, mu);

    this->fforce_site(mom + 18 * gid, gauge + 18 * gid,
                      v2 + fid, v2_s + fid_s,
                      v1 + fid, v1_s + fid_s, mu, coef);
  }
}

// compute fermion force for Mobius class fermions:
// This is the threaded equivalent of fbfm::EvolveMemFforceBase() in CPS.
//
// mom += coef * (phiL^\dag e_i(M) \phiR + \phiR^\dag e_i(M^\dag) \phiL)
// M = M_oo - M_oe M^{-1}_ee M_eo
//
// IMPORTANT: at least 5 threads are needed for this function to work
// correctly since we want to interleave communication and the
// evaluation of internal forces.
template<class Float>
void bfm_evo<Float>::compute_force(Float *mom, Float *gauge, Fermion_t phiL, Fermion_t phiR, double coef)
{
  int me = this->thread_barrier();

  Fermion_t v1[2] = {this->threadedAllocFermion(), this->threadedAllocFermion()};
  Fermion_t v2[2] = {this->threadedAllocFermion(), this->threadedAllocFermion()};

  this->calcMDForceVecs(v1, v2, phiL, phiR);

  // compute various sizes
  int lclx[5] = {this->node_latt[0], this->node_latt[1], this->node_latt[2], this->node_latt[3], this->Ls};
  int vol_5d = 24 * lclx[0] * lclx[1] * lclx[2] * lclx[3] * lclx[4];

  int surf_size[4];
  int surf_size_all = 0;
  for(int i = 0; i < 4; ++i) {
    surf_size[i] = vol_5d / lclx[i];
    surf_size_all += 2 * surf_size[i];
  }
  
  // calculate offset of surface vectors v1 and v2
  int surf_v1[4], surf_v2[4];
  surf_v1[0] = 0;
  surf_v2[0] = surf_size[0];
  for(int i = 1; i < 4; ++i) {
    surf_v1[i] = surf_v1[i-1] + surf_size[i-1] * 2;
    surf_v2[i] = surf_v1[i] + surf_size[i];
  }
  
  Float *v1f = this->threadedAllocFloat(vol_5d);
  Float *v2f = this->threadedAllocFloat(vol_5d);
  Float *sndbuf = this->threadedAllocFloat(surf_size_all);
  Float *rcvbuf = this->threadedAllocFloat(surf_size_all);

  this->thread_impexFermion_s(v1f, v1, 0);
  this->thread_impexFermion_s(v2f, v2, 0);

  for(int i = 0; i < 4; ++i) {
    this->copySendFrmData(sndbuf + surf_v1[i], v1f, i, true);
    this->copySendFrmData(sndbuf + surf_v2[i], v2f, i, true);
  }

  // Fused comm/internal force.
  //
  // The last 4 threads (typically 60-63) are used for
  // communication. All other threads (typically 0-59) are used to
  // calculate internal forces.
  if(this->nthread <= 4) {
    if(!me) {
      printf("compute_force: Oops, at least 5 threads are needed.\n");
    }
    exit(-1);
  }

  // parallelize comm/internal force calculation
  if(me >= this->nthread - 4) {
    int dir = this->nthread - me - 1;
    cps::getPlusData(rcvbuf + surf_v1[dir], sndbuf + surf_v1[dir],
                     surf_size[dir] * 2, dir);
  } else {
    for(int i = 0; i < 4; ++i) {
      fforce_internal(mom, gauge, v1f, v2f, coef, i, me, this->nthread - 4);
    }
  }

  this->thread_barrier();

  for(int i = 0; i < 4; ++i) {
    fforce_surface(mom, gauge, v1f, v2f,
                   rcvbuf + surf_v1[i], // v1 surface
                   rcvbuf + surf_v2[i], // v2 surface
                   coef, i);
  }

  this->threadedFreeFloat(v1f);
  this->threadedFreeFloat(v2f);
  this->threadedFreeFloat(sndbuf);
  this->threadedFreeFloat(rcvbuf);

  this->threadedFreeFermion(v1[0]);
  this->threadedFreeFermion(v1[1]);
  this->threadedFreeFermion(v2[0]);
  this->threadedFreeFermion(v2[1]);
}

// complex version of axpy()
// template<class Float>
// void bfm_evo<Float>::axpy_c(Fermion_t r, Fermion_t x, Fermion_t y, std::complex<double> a, Fermion_t tmp)
// {
//   this->copy(tmp, x);
//   this->scale(tmp, std::real(a), std::imag(a));
//   this->axpy(r, tmp, y, 1.0);
// }

// bicg_M: Biconjugate gradient method on preconditioned Dirac
// operator (It never converges).
//
// FIXME: test code only, don't use it unless you know what you are
// doing.
template<class Float>
int bfm_evo<Float>::bicg_M(Fermion_t sol, Fermion_t src)
{
  int me = this->thread_barrier();

  Fermion_t r   = this->threadedAllocFermion();
  Fermion_t rd  = this->threadedAllocFermion();
  Fermion_t p   = this->threadedAllocFermion();
  Fermion_t pd  = this->threadedAllocFermion();
  Fermion_t mp  = this->threadedAllocFermion();
  Fermion_t mdpd= this->threadedAllocFermion();
  Fermion_t x   = sol;
  Fermion_t xd  = this->threadedAllocFermion();
  this->copy(xd, x);

  Fermion_t tv1 = this->threadedAllocFermion();
  Fermion_t tv2 = this->threadedAllocFermion();

  const double src_norm = this->norm(src);
  const double stop = src_norm * this->residual * this->residual;

  this->Mprec(x , r , tv1, 0, 0);
  this->Mprec(xd, rd, tv1, 1, 0);
  double rnorm  = this->axpy_norm(r , r , src, -1.0); // r0 <- b-M*x0
  double rdnorm = this->axpy_norm(rd, rd, src, -1.0);// r0d <- b-Md*x0

  if ( this->isBoss() && !me ) {
    printf("iter = %5d rsd = %17.10e true rsd = %17.10e\n", 0, rnorm, rnorm);
  }

  this->copy(p, r);
  this->copy(pd, rd);

  std::complex<double> rddr = this->inner(rd, r);

  int k = 1;
  for(; k <= this->max_iter; ++k) {
    this->Mprec(p, mp, tv1, 0, 0);
    this->Mprec(pd, mdpd, tv1, 1, 0);

    std::complex<double> pddmp = this->inner(pd, mp);
    std::complex<double> alpha = rddr / pddmp;

    this->axpy_c(x , p , x , alpha, tv1); // x <- x + alpha * p
    this->axpy_c(xd, pd, xd, alpha, tv1); // xd <- xd + alpha * pd

    this->axpy_c(r , mp  , r , -alpha, tv1); // r <- r - alpha * Mp
    this->axpy_c(rd, mdpd, rd, -alpha, tv1); // rd <- rd - alpha * Mdpd

    rnorm = this->norm(r);
    rdnorm = this->norm(rd);

    // check stopping condition
    if(rnorm < stop) {
      // compute true residual
      this->Mprec(x, tv2, tv1, 0, 0);
      double true_rsd = this->axpy_norm(tv1, tv2, src, -1.0);

      if(this->isBoss() && !me) {
        printf("bicg_M: converged in %d iterations.\n", k);
        printf("bicg_M: acc_rsd = %9.3e %9.3e true_rsd = %9.3e\n",
               sqrt(rnorm/src_norm), sqrt(rdnorm/src_norm), sqrt(true_rsd/src_norm));
      }
      break;
    }

    std::complex<double> tmp = this->inner(rd, r);
    std::complex<double> beta = tmp / rddr;
    rddr = tmp;

    this->axpy_c(p , p , r , beta, tv1); // p <- r + beta * p
    this->axpy_c(pd, pd, rd, beta, tv1); // pd <- rd + beta * pd

    // ======================================================================
    // compare rsd and true rsd
    this->Mprec(x, tv2, tv1, 0, 0);
    double true_rsd = this->axpy_norm(tv2, tv2, src, -1.0);

    if ( this->isBoss() && !me ) {
        printf("iter = %5d rsd = %9.3e true rsd = %9.3e a = (%9.3e %9.3e) b = (%9.3e %9.3e)\n",
               k, rnorm, true_rsd, real(alpha), imag(alpha), real(beta), imag(beta));
    }
    // ======================================================================
  }

  if(k > this->max_iter) {
    if(this->isBoss() && !me) {
      printf("bicg_M: not converged in %d iterations.\n", k);
    }
  }

  this->threadedFreeFermion(r);
  this->threadedFreeFermion(rd);
  this->threadedFreeFermion(p);
  this->threadedFreeFermion(pd);
  this->threadedFreeFermion(mp);
  this->threadedFreeFermion(mdpd);
  this->threadedFreeFermion(xd);
  this->threadedFreeFermion(tv1);
  this->threadedFreeFermion(tv2);
  
  return k;
}

// bicgstab_M: Biconjugate gradient stabilized method on
// preconditioned Dirac operator.
//
// FIXME: test code only, don't use it unless you know what you are
// doing.
template<class Float>
int bfm_evo<Float>::bicgstab_M(Fermion_t sol, Fermion_t src)
{
  int me = this->thread_barrier();

  Fermion_t r0  = this->threadedAllocFermion();
  Fermion_t r   = this->threadedAllocFermion();
  Fermion_t p   = this->threadedAllocFermion();
  Fermion_t v   = this->threadedAllocFermion();
  Fermion_t s   = this->threadedAllocFermion();
  Fermion_t t   = this->threadedAllocFermion();
  Fermion_t x   = sol;
  Fermion_t tv1 = this->threadedAllocFermion();
  Fermion_t tv2 = this->threadedAllocFermion();

  const double src_norm = this->norm(src);
  const double stop = src_norm * this->residual * this->residual;

  this->Mprec(x, r0, tv1, 0, 0);
  double r0n = this->axpy_norm(r0, r0, src, -1.0); // r0 <- b-M*x0, r0^hat = r0
  this->copy(r, r0);

  if ( this->isBoss() && !me ) {
    printf("iter = %5d rsd = %17.10e true rsd = %17.10e\n", 0, r0n, r0n);
  }

  std::complex<double> rho(1, 0);
  std::complex<double> alpha(1, 0);
  std::complex<double> omega(1, 0);

  this->set_zero(v);
  this->set_zero(p);

  int k = 1;
  for(; k <= this->max_iter; ++k) {
    std::complex<double> rho_k = this->inner(r0, r);
    std::complex<double> beta = rho_k / rho * alpha / omega;
    rho = rho_k;

    this->axpy_c(tv1, v, p, -omega, tv2);
    this->axpy_c(p, tv1, r, beta, tv2);

    this->Mprec(p, v, tv1, 0, 0);

    alpha = rho / this->inner(r0, v);

    this->axpy_c(s, v, r, -alpha, tv1);
    this->Mprec(s, t, tv1, 0, 0);

    omega = this->inner(t, s) / this->norm(t);
    
    this->axpy_c(x, p, x, alpha, tv1);
    this->axpy_c(x, s, x, omega, tv1);
    this->axpy_c(r, t, s, -omega, tv1);

    // compute true residual
    this->Mprec(x, tv2, tv1, 0, 0);
    double true_rsd = this->axpy_norm(tv1, tv2, src, -1.0);

    // check stopping condition
    if(true_rsd < stop) {
      if(this->isBoss() && !me) {
        printf("bicgstab_M: converged in %d iterations.\n", k);
        printf("bicgstab_M: true_rsd = %10.3e\n", sqrt(true_rsd/src_norm));
      }
      break;
    }

    // ======================================================================
    // debug information
    if ( this->isBoss() && !me ) {
        printf("iter = %5d true rsd = %10.3e "
               "rho = (%10.3e %10.3e) alpha = (%10.3e %10.3e) omega = (%10.3e %10.3e)\n",
               k, true_rsd,
               real(rho), imag(rho),
               real(alpha), imag(alpha),
               real(omega), imag(omega));
    }
    // ======================================================================
  }

  if(k > this->max_iter) {
    if(this->isBoss() && !me) {
      printf("bicgstab_M: not converged in %d iterations.\n", k);
    }
  }

  this->threadedFreeFermion(r0);
  this->threadedFreeFermion(r);
  this->threadedFreeFermion(p);
  this->threadedFreeFermion(v);
  this->threadedFreeFermion(s);
  this->threadedFreeFermion(t);
  this->threadedFreeFermion(tv1);
  this->threadedFreeFermion(tv2);
  
  return k;
}

// copied from Jianglei's bfm
template<typename Float>
double bfm_evo<Float>::CompactMprec(Fermion_t compact_psi,
                                    Fermion_t compact_chi,
                                    Fermion_t psi,
                                    Fermion_t chi,
                                    Fermion_t tmp,
                                    int dag,int donrm)
{
  this->copy(psi, compact_psi);
  double result = this->Mprec(psi, chi, tmp, dag, donrm);
  this->copy(compact_chi, chi);
  return result;
}

// copied from Jianglei's bfm
template<typename Float>
void bfm_evo<Float>::CompactMunprec(Fermion_t compact_psi[2],
                                    Fermion_t compact_chi[2],
                                    Fermion_t psi[2],
                                    Fermion_t chi[2],
                                    Fermion_t tmp,
                                    int dag)
{
  this->copy(psi[0], compact_psi[0]);
  this->copy(psi[1], compact_psi[1]);
  this->Munprec(psi, chi, tmp, dag);
  this->copy(compact_chi[0], chi[0]);
  this->copy(compact_chi[1], chi[1]);
}

template<typename Float>
void bfm_evo<Float>::deflate(Fermion_t out, Fermion_t in,
                             const multi1d<Fermion_t [2]> *evec,
                             const multi1d<Float> *eval,
                             int N)
{
  if(N == 0 || evec == NULL || eval == NULL) {
    if(this->isBoss()) {
      printf("bfm_evo::deflate() must provide eigenvectors.\n");
    }
    exit(-1);
  }

  this->set_zero(out);
  for(int i = 0; i < N; ++i) {
    std::complex<double> dot = this->inner((*evec)[i][1], in);
    this->zaxpy(out, (*evec)[i][1], out, dot / double((*eval)[i]));
  }
}

// GCR, the matrix is preconditioned M.
template<class Float>
int bfm_evo<Float>::gcr_M(Fermion_t sol, Fermion_t src)
{
  int me = this->thread_barrier();

  Fermion_t r   = this->threadedAllocFermion();
  Fermion_t gr  = this->threadedAllocFermion();
  Fermion_t agr = this->threadedAllocFermion();
  Fermion_t p   = this->threadedAllocFermion();
  Fermion_t ap  = this->threadedAllocFermion();
  Fermion_t x   = sol;
  Fermion_t tv1 = this->threadedAllocFermion();
  Fermion_t tv2 = this->threadedAllocFermion();

  const double src_norm = this->norm(src);
  const double stop = src_norm * this->residual * this->residual;

  this->Mprec(x, r, tv2, 0, 0);
  double rnorm = this->axpy_norm(r, r, src, -1.0); // r <- b - M x

  if ( this->isBoss() && !me ) {
    std::printf("gcr_M: iter = %5d rsd = %10.3e true rsd = %10.3e\n",
                0, std::sqrt(rnorm / src_norm),
                std::sqrt(rnorm / src_norm));
  }

  this->g5r5(gr, r);
  this->Mprec(gr, agr, tv1, 0, 0);
  this->copy(p, gr);
  this->copy(ap, agr);

  std::complex<double> ragr = this->inner(r, agr);

  int k = 1;
  for(; k <= this->max_iter; ++k) {
    double pdmmp = this->norm(ap);

    std::complex<double> alpha = ragr / pdmmp;
    this->zaxpy(x, p, x, alpha);
    this->zaxpy(r, ap, r, -alpha);
    rnorm = this->norm(r);

    if(rnorm < stop) {
      if(this->isBoss() && !me) {
        std::printf("gcr_M: converged in %d iterations.\n", k);
        std::printf("gcr_M: rsd = %10.3e\n", std::sqrt(rnorm/src_norm));
      }
      break;
    }

    this->g5r5(gr, r);
    this->Mprec(gr, agr, tv2, 0, 0);

    std::complex<double> ragrn = this->inner(r, agr);
    std::complex<double> beta = ragrn / ragr;
    ragr = ragrn;

    this->zaxpy(p, p, gr, beta);
    this->zaxpy(ap, ap, agr, beta);

    // ======================================================================
    // Computing true residual and other information, the
    // following can be removed without any effect on convergence.
    this->Mprec(x, tv1, tv2, 0, 0);
    double true_rsd = this->axpy_norm(tv1, tv1, src, -1.0);

    if ( this->isBoss() && !me ) {
      std::printf("gcr_M: iter = %5d rsd = %10.3e true_rsd = %10.3e\n",
                  k,
                  std::sqrt(rnorm / src_norm),
                  std::sqrt(true_rsd / src_norm));
    }
    // ======================================================================
  }

  if(k > this->max_iter) {
    if(this->isBoss() && !me) {
      std::printf("gcr_M: not converged in %d iterations.\n", k);
    }
  }

  this->Mprec(x, tv1, tv2, 0, 0);
  double true_rsd = this->axpy_norm(tv1, tv1, src, -1.0);

  if(this->isBoss() && !me) {
    std::printf("gcr_M: true_rsd = %10.3e\n",
                std::sqrt(true_rsd/src_norm));
  }

  this->threadedFreeFermion(r);
  this->threadedFreeFermion(gr);
  this->threadedFreeFermion(agr);
  this->threadedFreeFermion(p);
  this->threadedFreeFermion(ap);
  this->threadedFreeFermion(tv1);
  this->threadedFreeFermion(tv2);
  
  return k;
}

// GMRES(m), we restart after m iterations.
template<class Float>
int bfm_evo<Float>::gmres_M(Fermion_t sol, Fermion_t src, const int m)
{
  using namespace std;
  typedef complex<double> cmplx;

  int me = this->thread_barrier();

  Fermion_t r   = this->threadedAllocFermion();
  Fermion_t w   = this->threadedAllocFermion();
  Fermion_t tv1 = this->threadedAllocFermion();

  // the history of search directions
  vector<Fermion_t> v(m + 1, NULL);
  for(int i = 0; i <= m; ++i) {
    v[i] = this->threadedAllocFermion();
  }

  vector<cmplx> H((m + 1) * m, 0);
  vector<cmplx> R((m + 1) * m, 0);
  vector<cmplx> B(m, 0);

  vector<cmplx> C(m, 0);
  vector<cmplx> S(m, 0);
  vector<cmplx> Y(m, 0);

  const double len = sqrt(this->norm(src));
  const double stop = len * this->residual;

  this->Mprec(sol, r, tv1, 0, 0);
  double rsq = this->axpy_norm(r, r, src, -1.0); // r <- b - M x

  int j = 0;
  for(; j < this->max_iter / m; ++j) {
    double beta = sqrt(rsq);
    this->axpy(v[0], r, r, 1/beta - 1); // v[0] <- r / beta

    B.assign(m, 0);
    B[0] = beta;

    int nr = m;
    double rho = len;

    for(int i = 0; i < m; ++i) {
      this->Mprec(v[i], w, tv1, 0, 0);

      // Arnoldi iteration
      for(int k = 0; k <= i; ++k) {
        H[k*m+i] = this->inner(v[k], w);
        this->zaxpy(w, v[k], w, -H[k*m+i]);
      }
      double w2 = sqrt(this->norm(w));

      H[(i+1)*m+i] = w2;
      this->axpy(v[i+1], w, w, 1/w2 - 1);
            
      R[0*m+i] = H[0*m+i];

      // Givens transformation
      for(int k = 1; k <= i; ++k) {
        cmplx gamma = C[k-1] * R[(k-1)*m+i] + conj(S[k-1]) * H[k*m+i];
        R[k*m+i] = -S[k-1] * R[(k-1)*m+i] + C[k-1] * H[k*m+i];
        R[(k-1)*m+i] = gamma;
      }

      double rii = norm(R[i*m+i]);
      double hii = norm(H[(i+1)*m+i]);
      double delta = sqrt(rii + hii);

      cmplx mu, tau;
      if(rii < hii) {
        mu = R[i*m+i] / H[(i+1)*m+i];
        tau = conj(mu) / abs(mu);
      } else {
        mu = H[(i+1)*m+i] / R[i*m+i];
        tau = mu / abs(mu);
      }

      C[i] = sqrt(rii) / delta;
      S[i] = sqrt(hii) * tau / delta;

      R[i*m+i] = C[i] * R[i*m+i] + conj(S[i]) * H[(i+1)*m+i];
      B[i+1] = -S[i] * B[i];
      B[i] *= C[i];

      rho = abs(B[i+1]);

      if(this->isBoss() && !me) {
        std::printf("gmres: (j i) = %4d %4d rsd = %10.3e\n",
                    j, i, rho / len);
      }

      if(rho < stop) {
        nr = i;
        break;
      }
    }

    for(int k = nr - 1; k >= 0; --k) {
      Y[k] = B[k];
      for(int i = k + 1; i < nr; ++i) {
        Y[k] -= R[k*m+i] * Y[i];
      }
      Y[k] /= R[k*m+k];

      this->zaxpy(sol, v[k], sol, Y[k]);
    }

    this->Mprec(sol, r, tv1, 0, 0);
    rsq = this->axpy_norm(r, r, src, -1.0);
    if(rho < stop) break;
  }

  if(j >= this->max_iter / m) {
    if(this->isBoss() && !me) {
      std::printf("gmres: not converged in %d iterations.\n", j);
    }
  }

  if(this->isBoss() && !me) {
    std::printf("gmres: true_rsd = %10.3e\n",
                std::sqrt(rsq) / len);
  }

  this->threadedFreeFermion(r);
  this->threadedFreeFermion(w);
  this->threadedFreeFermion(tv1);
  
  for(int i = 0; i <= m; ++i) {
    this->threadedFreeFermion(v[i]);
  }

  return j;
}

#endif
