/* -*- mode:c++; c-basic-offset:2 -*- */
/****************************************************************************/
/* Apr 2012                                                                 */
/* Hantao Yin                                                               */
/*                                                                          */
/* bfm_evo.h                                                                */
/* functions added by Hantao (mainly mobius evolution related).             */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_BFM_EVO_HT_H
#define INCLUDED_BFM_EVO_HT_H

#include <stdio.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <omp.h>
#include <math.h>
#include "bfm_evo_aux.h"

struct rquo_approx {
  double bsn_a0;
  double frm_a0;
  vector<double> bsn_alpha;
  vector<double> bsn_beta;
  vector<double> frm_alpha;
  vector<double> frm_beta;
};

// currently it uses 4D preconditioning only.
//
// FIXME: it inherits from bfm_qdp for the sole reason of using its
// importGauge() function. I'm too lazy to do any manual
// shifts/conjugates ......
template <class Float>
class bfm_evo : public bfm_qdp<Float> {
public:
  enum {Even = 0, Odd};

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
  int prop_solve(Fermion_t out[2], Fermion_t in[2]);

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
                              Float coef, int mu);

  void fforce_surface(Float *mom, Float *gauge,
                      Float *v1, Float *v2, // internal data
                      Float *v1_s, Float *v2_s, // surface data
                      Float coef, int mu);

  void copySendFrmData(Float v3d[], Float v4d[], int mu, bool send_neg);
public:
  // compute fermion force:
  //
  // mom += coef * (phiL^\dag e_i(M) \phiR + \phiR^\dag e_i(M^\dag) \phiL)
  //
  // For BFM M is M = M_oo - M_oe M^{-1}_ee M_eo
  void compute_force(Float *mom, Float *gauge, Fermion_t phiL, Fermion_t phiR, double coef);
  
  // compute quotient and rational quotient fermion forces
  void compute_force_quotient(Float *mom, Float *gauge, Fermion_t phi, double mb, double mf);
  // not implemented yet.

  void compute_force_rquotient(Float *mom, Float *gauge, Fermion_t phi, double mb, double mf,
                               const rquo_approx &rq_param);
  // not implemented yet.

  // psi assumes the following order: (color, spin, s, x, y, z, t),
  // mainly used to import/export the "v1" and "v2" vectors in evolution.
  template<typename FloatEXT>
  void thread_impexFermion_s(FloatEXT *psi, Fermion_t handle[2], int doimport);

  Float *threadedAllocFloat(size_t size, int mem_type=mem_slow);
  void threadedFreeFloat(Float *);

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

#pragma omp parallel
  {
#pragma omp for 
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
  }//omp
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

#pragma omp parallel
  {
#pragma omp for 
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
  }//omp
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

#pragma omp parallel
  {
#pragma omp for 
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
  }//omp
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

#pragma omp parallel
    {

#pragma omp for 
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
    }
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
  this->MooeeInv(v1[Odd], v1[Even], DaggerYes); // v1e in place

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

template <class Float>
int bfm_evo<Float>::prop_solve(Fermion_t out[2], Fermion_t in[2])
{
  int me = this->thread_barrier();

  // DWF
  if(this->solver == DWF && this->precon_5d == 1 ||
     this->solver == HtCayleyTanh && this->precon_5d == 0) {
    return this->CGNE_M(out, in);
  }

  // Mobius
  if(this->solver == HmCayleyTanh && this->precon_5d == 0) {
    Fermion_t tmp[2] = {
      this->threadedAllocFermion(),
      this->threadedAllocFermion()};

    this->set_zero(tmp[0]);
    this->set_zero(tmp[1]);

    this->G5D_Dminus(in, tmp, 0);
    int i = this->CGNE_M(out, tmp);

    this->threadedFreeFermion(tmp[0]);
    this->threadedFreeFermion(tmp[1]);

    return i;
  }

  if(this->isBoss() && !me) {
    printf("prop_solve not implemented\n");
  }
  exit(-1);
  return 0;
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
      sprojTrXm(t1, v1p, v2, this->Ls, 0, 0);
      sprojTrXp(t2, v2p, v1, this->Ls, 0, 0);
      break;
    case 1:
      sprojTrYm(t1, v1p, v2, this->Ls, 0, 0);
      sprojTrYp(t2, v2p, v1, this->Ls, 0, 0);
      break;
    case 2:
      sprojTrZm(t1, v1p, v2, this->Ls, 0, 0);
      sprojTrZp(t2, v2p, v1, this->Ls, 0, 0);
      break;
    default:
      sprojTrTm(t1, v1p, v2, this->Ls, 0, 0);
      sprojTrTp(t2, v2p, v1, this->Ls, 0, 0);
    }

    su3_add(t1, t2);
    mDotMEqual(t2, gauge, t1);
    trless_am(t2, -coef);
    su3_add(mom, t2);
}

template<class Float>
void bfm_evo<Float>::fforce_internal(Float *mom, Float *gauge,
                                     Float *v1, Float *v2, // internal data
                                     Float coef, int mu)
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

// compute fermion force:
//
// mom += coef * (phiL^\dag e_i(M) \phiR + \phiR^\dag e_i(M^\dag) \phiL)
//
// For BFM M is M = M_oo - M_oe M^{-1}_ee M_eo
//
// Be careful about the boundary conditions.
// Currently only for DWF like fermions.
// This is the threaded equivalent of fbfm::EvolveMemFforceBase() in CPS.
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

  // FIXME: can I have thread 0 do the communication only while the
  // rest of the threads compute internal forces? If this is possible
  // then we can turn blocking communications into non-blocking ones
  // (it blocks only 1 thread).
  if(!me) {
    for(int i = 0; i < 4; ++i) {
      cps::getPlusData(rcvbuf + surf_v1[i], sndbuf + surf_v1[i],
                       surf_size[i] * 2, i); 
    }
  }

  this->thread_barrier();
  for(int i = 0; i < 4; ++i) {
    fforce_internal(mom, gauge, v1f, v2f, coef, i);
  }

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

#endif
