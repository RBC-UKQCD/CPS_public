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
#include <bagel_int.h>
#include <bfm.h>
#include <bfm_qdp.h>
#include <omp.h>
#include <math.h>
#include <vector>
#include <util/gjp.h>
#include "bfm_evo_aux.h"

enum
{ Export = 0, Import = 1 };
// FIXME: it inherits from bfm_qdp for the sole reason of using its
// importGauge() function. I'm too lazy to do any manual
// shifts/conjugates ......
template <class Float>
class bfm_evo : public bfm_qdp<Float> {
public:
  // BFM has this
  // enum {Even = 0, Odd};

  integer cps_idx_cb(int x[4], int s, int reim, int i, int i_size);
  integer cps_idx_cb_gparity(int x[4], int s, int reim, int i, int i_size, int flav);

  // s outer most
  integer cps_idx(int x[4], int s, int reim, int i, int i_size);
  integer cps_idx_gparity(int x[4], int s, int reim, int i, int i_size, int flav);

  // s inner most (but outside color and spin)
  integer cps_idx_s(int x[4], int s, int reim, int i, int i_size);
    // index for 4d fermion field
  integer cps_idx_4d (int x[4], int reim, int i, int i_size);
  integer cps_idx_s_gparity(int x[4], int s, int reim, int i, int i_size, int flav);

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
#if 0 //testing
private:
#endif
  // auxiliary functions used by compute_force().
  //CK: gpf1_offset_p is the offset to reach the second G-parity flavour in the vectors v1p and v2p.
  //    Its value depends on whether v1p/v2p are internal vectors (24*5dvol) or in the buffer send 
  //    from the next node (24*Ls*3dsurfvol  where 3dsurfvol is the 3d surface volume in the comms direction)
  void fforce_site(Float *mom, Float *gauge,
                   Float *v1, Float *v1p,
                   Float *v2, Float *v2p, int mu, Float coef, int gpf1_offset_p = 0);

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
    printf("void axpy_c temporarily disabled\n");
    exit(-1);
#if 0
    this->zaxpy(r, x, y, a);
#endif
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

#if 0
  //CHECK CODE
  template<typename FloatEXT>
  void thread_impexFermion_s_test(FloatEXT *psi, Fermion_t handle[2], int doimport);
#endif  

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

//  template<typename FloatEXT>
//  void cps_importGauge(FloatEXT *importme);
  // Imports a 4D CPS fermion to a 5D BFM fermion, putting the left-handed
  // part at s=0 and the right-handed part at s=Ls-1. (Or does the inverse,
  // exporting a 5D BFM fermion to a 4D CPS fermion).
  // psi assumes regular canonical order: (color, spin, x, y, z, t)
  template < typename FloatEXT >
    void cps_impexFermion_4d (FloatEXT * psi, Fermion_t handle[2],
			      int doimport, bool prezero = true)
// Imports a 4D CPS fermion to a 5d BFM fermion, putting the left-handed
// part at s=0 and the right-handed part at s=Ls-1. (Or does the inverse,
// exporting a 5D BFM fermion to a 4D CPS fermion).
//template < class Float > template < typename FloatEXT >
//  void bfm_evo < Float >::cps_impexFermion_4d (FloatEXT * psi,
//					       Fermion_t handle[2],
//					       int doimport, bool prezero)
{
  if (doimport && prezero)
    {
#pragma omp parallel
      {
	// zero out 5d bulk since we only import to the walls
	this->set_zero (handle[Even]);
	this->set_zero (handle[Odd]);
      }
    }

  int Nspinco = 12;
  int i_inc = this->simd () * 2;
  int vol4d =
    this->node_latt[0] *
    this->node_latt[1] * this->node_latt[2] * this->node_latt[3];

  omp_set_num_threads (this->nthread);
  Float *bagel[2] = { (Float *) handle[0], (Float *) handle[1] };

#pragma omp parallel for
  for (int site = 0; site < vol4d; site++)
    {
      int x[4];
      int si = site;
      x[0] = si % this->node_latt[0];
      si = si / this->node_latt[0];
      x[1] = si % this->node_latt[1];
      si = si / this->node_latt[1];
      x[2] = si % this->node_latt[2];
      si = si / this->node_latt[2];
      x[3] = si % this->node_latt[3];

      int bidx_base_left = this->bagel_idx5d (x, 0, 0, 0, Nspinco, 1);
      int bidx_base_right =
	this->bagel_idx5d (x, this->Ls - 1, 0, 0, Nspinco, 1);
      int cidx_base = this->cps_idx_4d (x, 0, 0, Nspinco);

      for (int co = 0; co < Nspinco; co++)
	{
	  // right-handed components are first six spin-color components
	  // left-handed components are last six spin-color components
	  int bidx_base;
	  int s;
	  if (co < 6)
	    {
	      bidx_base = bidx_base_right;
	      s = this->Ls - 1;
	    }
	  else
	    {
	      bidx_base = bidx_base_left;
	      s = 0;
	    }
	  int sp = this->precon_5d ? s : 0;
	  int cb = (x[0] + x[1] + x[2] + x[3] + sp) & 0x1;

	  for (int reim = 0; reim < 2; reim++)
	    {
	      int bidx = bidx_base + reim + co * i_inc;
	      int cidx = cidx_base + reim + co * 2;

	      if (doimport)
		bagel[cb][bidx] = psi[cidx];
	      else
		psi[cidx] = bagel[cb][bidx];
	    }
	}			//co, reim
    }				//xyzts
}

  template < typename FloatEXT > void cps_importGauge (FloatEXT * importme);
#if 0
//CK: Appears to assume 'importme' is in canonical ordering
//template <class Float> template<typename FloatEXT>
//void bfm_evo<Float>::cps_importGauge(FloatEXT *importme)
{
  int u_sz = Nd;
  if(cps::GJP.Gparity()) u_sz *= 2; //U* fields are stacked on second set of Nd LatticeColorMatrix objects in the array
  
  multi1d<LatticeColorMatrix> U(u_sz);
  omp_set_num_threads(this->nthread);

  int Ndircoco = 72;
  int Ncoco = 9;
  QDPdouble *U_p;

  int vol4d =
    this->node_latt[0] *
    this->node_latt[1] * this->node_latt[2] * this->node_latt[3];
  assert (vol4d>0 );

  for (int muu=0;muu<u_sz;muu++) {
    U_p = (QDPdouble *)&(U[muu].elem(0).elem());
    int flav = muu / Nd; int mu = muu % Nd;
#pragma omp parallel for 
    for (int site=0;site<vol4d;site++ ) {
      int x[4];
      int s=site;
      x[0]=s%this->node_latt[0];    s/=this->node_latt[0];
      x[1]=s%this->node_latt[1];    s/=this->node_latt[1];
      x[2]=s%this->node_latt[2];    s/=this->node_latt[2];
      x[3]=s%this->node_latt[3];
      int qidx_base = this->chroma_idx(x, 0, 0, Ncoco);

      for(int coco = 0; coco < Ncoco; ++coco) {
         for ( int reim = 0; reim < 2; ++reim) {
      
		int qidx = qidx_base + reim + coco * 2;
          int siteoff = mu + Nd * site + flav*Nd*vol4d; //Second G-parity flavour offset by Nd*vol4d
          int cidx = reim + 2 * (coco + Ncoco * siteoff);
          U_p[qidx] = importme[cidx];
        }} // reim,coco
    } // x
  }//mu
//  if(this->isBoss()) printf("before importGauge\n");
  // to bfm
  this->importGauge (U);
//  if(this->isBoss()) printf("after importGauge\n");
}
#endif

  //EigCG
#if 0 //THESE ARE IN BFM
  Fermion_t allocCompactFermion   (int mem_type=mem_slow);
  Fermion_t threadedAllocCompactFermion   (int mem_type=mem_slow);
  void* threaded_alloc(int length, int mem_type=mem_slow);
  void threaded_free(void *handle);
#endif

  int EIG_CGNE_M(Fermion_t solution[2], Fermion_t source[2]);
  int Eig_CGNE_prec(Fermion_t psi, Fermion_t src);

#if 0 //CK: leaving them in BFM
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
#endif

  // do deflation using eigenvectors/eigenvalues from Rudy's Lanczos code.
  void deflate(Fermion_t out, Fermion_t in,
               const multi1d<Fermion_t [2]> *evec,
               const multi1d<Float> *eval, int N);
  void set_mass (double mass);

//#ifdef USE_NEW_BFM_GPARITY
#if 1
  inline void axpby_ssp_proj(Fermion_t out, std::complex<double> a,Fermion_t x, std::complex<double> b,Fermion_t y,int sxo,int sy,int psign){
    this->axpby_ssp_proj_complex(out,a.real(),a.imag(),x,b.real(),b.imag(),y,sxo,sy,psign);
  }
#endif
};

  // Simple utility function to set the mass and reinit if necessary.

template < class Float > void bfm_evo < Float >::set_mass (double mass)
{
  if (this->mass != mass)
    {
      this->mass = mass;
      this->GeneralisedFiveDimEnd ();
      this->GeneralisedFiveDimInit ();
    }
}

//CK: this function gives the offset within a checkerboarded vector
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

//For G-parity the WILSON layout is 5d preconditioned
//|       s=0       |        s=1        | ......... |        s = 0      |.....
//| odd f0 | odd f1 | even f0 | even f1 | ......... | even f0 | even f1 |.....
//where the blocks on the lowest line have their *4d* parity indicated. (5d parity) = [(4d parity) + s] % 2
//hence the first half of the full WILSON vector had 5d parity odd, and the second half 5d parity even

template<class Float>
integer bfm_evo<Float>::cps_idx_cb_gparity(int x[4], int s, int reim, int i, int i_size, int flav)
{
  int s_off = this->node_latt[0] * this->node_latt[1] * this->node_latt[2] * this->node_latt[3]; //2 4D half-volumes, one for each flavour
  int f_off = s_off/2;

  int csite
    =x[0] +this->node_latt[0]
    *(x[1] + this->node_latt[1]
      *(x[2] +this->node_latt[2] * x[3]));
  csite /= 2;
  
  csite += flav * f_off + s*s_off;
      
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
integer bfm_evo<Float>::cps_idx_gparity(int x[4], int s, int reim, int i, int i_size, int flav)
{
  //For G-parity we have 2 flavours on each s-slice
  int s_off = 2*this->node_latt[0] * this->node_latt[1] * this->node_latt[2] * this->node_latt[3]; //2 4D volumes, one for each flavour
  int f_off = s_off/2;

  int csite =
    x[0] + this->node_latt[0]
    *(x[1] + this->node_latt[1]
      *(x[2] +this->node_latt[2]*x[3]));
  csite += s*s_off + flav * f_off;

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

template < class Float >
  integer bfm_evo < Float >::cps_idx_4d (int x[4], int reim, int i,
                                         int i_size)
{
  int csite =
    x[0] + this->node_latt[0]
    * (x[1] + this->node_latt[1] * (x[2] + this->node_latt[2] * (x[3])));

  return (csite * i_size + i) * 2 + reim;
}

template<class Float>
integer bfm_evo<Float>::cps_idx_s_gparity(int x[4], int s, int reim, int i, int i_size, int flav)
{
  //This s-inner mapping is new here. Offset the second flavour by 1 5D volume, just like in bfm
  int f_off = this->Ls * this->node_latt[0] * this->node_latt[1] * this->node_latt[2] * this->node_latt[3];

  int csite =
    s + this->Ls
    *(x[0] + this->node_latt[0]
      *(x[1] + this->node_latt[1]
        *(x[2] +this->node_latt[2]
          *x[3])));
  csite += flav * f_off;

  return (csite*i_size + i)*2 + reim;
}

//CK: Note if the BFM preconditioning is 4D then the 4D checkerboard of the imported field will be the opposite of the 5D checkerboard of the CPS field! cb is the output checkerboard.
//The set of all sites with  x+y+z+t+s odd is the same as the set of sites with x+y+z+t even, and vice versa.
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

  int work = vol5d;
  if(cps::GJP.Gparity()) work*=2;

#pragma omp parallel for 
  for (int sf = 0; sf < work; sf++) {
    int flav = sf;
    int site = flav % vol5d; flav /= vol5d;
    
    int x[4], s;
    int si=site;
    x[0]=si%this->node_latt[0];    si=si/this->node_latt[0];
    x[1]=si%this->node_latt[1];    si=si/this->node_latt[1];
    x[2]=si%this->node_latt[2];    si=si/this->node_latt[2];
    x[3]=si%this->node_latt[3];
    s   =si/this->node_latt[3];
    
    int sp = this->precon_5d ? s : 0;
    if ( (x[0]+x[1]+x[2]+x[3] + (sp &0x1)) == cb ) {

      int bidx_base;
      int cidx_base;
#ifdef BFM_GPARITY
      bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1,flav);      
      cidx_base = cps::GJP.Gparity() ? this->cps_idx_cb_gparity(x, s, 0, 0, Nspinco, flav) : this->cps_idx_cb(x, s, 0, 0, Nspinco);
#else
	bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1);
	cidx_base = this->cps_idx_cb(x, s, 0, 0, Nspinco);
#endif
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


//Convert a bfm-style Fermion_t pair to or from a CANONICAL format CPS-style fermion
//if doimport == 0 psi is the output and handle the input
//if doimport == 1 handle is the output and psi the input
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

  int work = vol5d;
  if(cps::GJP.Gparity()) work*=2;

#pragma omp parallel for 
  for (int sf = 0; sf < work; sf++) {
    int flav = sf;
    int site = flav % vol5d; flav /= vol5d;

    int x[4], s;
    int si=site;
    x[0]=si%this->node_latt[0];    si=si/this->node_latt[0];
    x[1]=si%this->node_latt[1];    si=si/this->node_latt[1];
    x[2]=si%this->node_latt[2];    si=si/this->node_latt[2];
    x[3]=si%this->node_latt[3];
    s   =si/this->node_latt[3];
    
    int sp = this->precon_5d ? s : 0;
    int cb = x[0]+x[1]+x[2]+x[3]+sp &0x1;

    int bidx_base;
    int cidx_base;

#ifdef BFM_GPARITY

    bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1, flav);
    cidx_base = cps::GJP.Gparity() ? this->cps_idx_gparity(x, s, 0, 0, Nspinco, flav) : this->cps_idx(x, s, 0, 0, Nspinco);
#else
      bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1);
      cidx_base = this->cps_idx(x, s, 0, 0, Nspinco);
#endif
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

//Convert a bfm style Fermion_t pair (left,right) to a 's-ordered' fermion
//if doimport == 0 the input is handle and the output psi
//if doimport == 1 the input is psi and the output handle
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

  int work = vol5d;
  if(cps::GJP.Gparity()) work*=2;

#pragma omp parallel for 
  for (int sf = 0; sf < work; sf++) {
    int flav = sf;
    int site = flav % vol5d; flav /= vol5d;

    int x[4], s;
    int si=site;
    s   =si%this->Ls;              si=si/this->Ls;
    x[0]=si%this->node_latt[0];    si=si/this->node_latt[0];
    x[1]=si%this->node_latt[1];    si=si/this->node_latt[1];
    x[2]=si%this->node_latt[2];
    x[3]=si/this->node_latt[2];

    int sp = this->precon_5d ? s : 0;
    int cb = (x[0]+x[1]+x[2]+x[3]+sp) & 0x1;


    int bidx_base;
    int cidx_base;
#ifdef BFM_GPARITY
    bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1, flav);
    cidx_base = cps::GJP.Gparity() ? this->cps_idx_s_gparity(x, s, 0, 0, Nspinco, flav) : this->cps_idx_s(x, s, 0, 0, Nspinco);
#else
      bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1);
      cidx_base = this->cps_idx_s(x, s, 0, 0, Nspinco);
#endif
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


#if 0
//This is check code

template <class Float> template<typename FloatEXT>
void bfm_evo<Float>::thread_impexFermion_s_test(FloatEXT *psi, Fermion_t handle[2], int doimport)
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

  int work = vol5d;
  if(cps::GJP.Gparity()) work*=2;
  this->thread_work(work, me, thrlen, throff);

  Float *bagel[2] = { (Float *)handle[0], (Float *)handle[1] };

  for (int site = 0; site < thrlen; ++site) {
    int flav = site + throff;
    int site = flav % vol5d; flav /= vol5d;

    int x[4], s;
    int si=site;
    s   =si%this->Ls;              si=si/this->Ls;
    x[0]=si%this->node_latt[0];    si=si/this->node_latt[0];
    x[1]=si%this->node_latt[1];    si=si/this->node_latt[1];
    x[2]=si%this->node_latt[2];
    x[3]=si/this->node_latt[2];

    int sp = this->precon_5d ? s : 0;
    int cb = x[0]+x[1]+x[2]+x[3]+sp & 0x1;

    int bidx_base;
    int cidx_base;
    if(cps::GJP.Gparity()){
      bidx_base = this->bagel_gparity_idx5d(x, s, 0, 0, Nspinco, 1, flav);
      cidx_base = this->cps_idx_s_gparity(x, s, 0, 0, Nspinco, flav);
    }else{
      bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1);
      cidx_base = this->cps_idx_s(x, s, 0, 0, Nspinco);
    }

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

#endif








//Convert a bfm style Fermion_t pair (left,right) to a 's-ordered' fermion
//if doimport == 0 the input is handle and the output psi
//if doimport == 1 the input is psi and the output handle
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

  int work = vol5d;
  if(cps::GJP.Gparity()) work*=2;
  this->thread_work(work, me, thrlen, throff);

  Float *bagel[2] = { (Float *)handle[0], (Float *)handle[1] };

  for (int sf = 0; sf < thrlen; ++sf) {
    int flav = sf + throff;
    int site = flav % vol5d; flav /= vol5d;

    int x[4], s;
    int si=site;
    s   =si%this->Ls;              si=si/this->Ls;
    x[0]=si%this->node_latt[0];    si=si/this->node_latt[0];
    x[1]=si%this->node_latt[1];    si=si/this->node_latt[1];
    x[2]=si%this->node_latt[2];
    x[3]=si/this->node_latt[2];

    int sp = this->precon_5d ? s : 0;
    int cb = x[0]+x[1]+x[2]+x[3]+sp & 0x1;

    int bidx_base;
    int cidx_base;
#ifdef BFM_GPARITY
    bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1, flav);
    cidx_base = cps::GJP.Gparity() ? this->cps_idx_s_gparity(x, s, 0, 0, Nspinco, flav) : this->cps_idx_s(x, s, 0, 0, Nspinco);
#else
      bidx_base = this->bagel_idx5d(x, s, 0, 0, Nspinco, 1);
      cidx_base = this->cps_idx_s(x, s, 0, 0, Nspinco);
#endif
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

//CK: Appears to assume 'importme' is in canonical ordering
template <class Float> template<typename FloatEXT>
void bfm_evo<Float>::cps_importGauge(FloatEXT *importme)
{
  int u_sz = Nd;
  if(cps::GJP.Gparity()) u_sz *= 2; //U* fields are stacked on second set of Nd LatticeColorMatrix objects in the array
  
  multi1d<LatticeColorMatrix> U(u_sz);
  omp_set_num_threads(this->nthread);

  int Ndircoco=72;
  int Ncoco = 9;
  QDPdouble *U_p;

  int vol4d =
    this->node_latt[0] *
    this->node_latt[1] *
    this->node_latt[2] *
    this->node_latt[3];

  for (int muu=0;muu<u_sz;muu++) {
    U_p = (QDPdouble *)&(U[muu].elem(0).elem());
    int flav = muu / Nd; int mu = muu % Nd;

#pragma omp parallel for 
    for (int site=0;site<vol4d;site++ ) {
      int x[4];
      int s=site;
      x[0]=s%this->node_latt[0];    s/=this->node_latt[0];
      x[1]=s%this->node_latt[1];    s/=this->node_latt[1];
      x[2]=s%this->node_latt[2];    s/=this->node_latt[2];
      x[3]=s%this->node_latt[3];
      
      int qidx_base = this->chroma_idx(x, 0, 0, Ncoco);

      for(int coco = 0; coco < Ncoco; ++coco) {
        for ( int reim = 0; reim < 2; ++reim) {
          // int qidx = this->chroma_idx(x,reim,coco,Ncoco);
          int qidx = qidx_base + reim + coco * 2;

          int siteoff = mu + Nd * site + flav*Nd*vol4d; //Second G-parity flavour offset by Nd*vol4d
          int cidx = reim + 2 * (coco + Ncoco * siteoff);
          U_p[qidx] = importme[cidx];
        }} // reim,coco
    } // x
  }//mu

  // to bfm
  this->importGauge(U);
}


//CK: phi1 = Mprec phi2 for fermionic vectors.

//Calculates: (odd,even)
//v2 = (Boo phi2, Bee MeeInv Meo phi2)
//v1 = (phi1,  MeeInv^dag Meo^dag phi1)

template <class Float>
void bfm_evo<Float>::calcMDForceVecs(Fermion_t v1[2], Fermion_t v2[2],
                       Fermion_t phi1, Fermion_t phi2)
{
  // Meo is Wilson D times a matrix (see page 27 in Peter's draft).
  // Moe/Meo: check bfmbase<Float>::G5D_Meo() in bfmdperp.C.
  // Mee/Moo: check bfmbase<Float>::G5D_Mooee().
  // Mee/Moo inverse: check bfmbase<Float>::G5D_MooeeInv().
  
  //2kappa = 1/(5-M5)
  
  // v2e      =  Bee * 2kappa * Meo phi2
  this->Meo(phi2, v1[Odd], Even, DaggerNo); //Uses v1[Odd] as temp storage
  this->MooeeInv(v1[Odd], v1[Even], DaggerNo); 
  this->Booee(v1[Even], v2[Even], DaggerNo);

  // v2o = Boo phi2
  this->Booee(phi2, v2[Odd], DaggerNo);

  // v1e  =  2kappa Meo^dag phi1
  this->Meo(phi1, v1[Odd], Even, DaggerYes);
  this->MooeeInv(v1[Odd], v1[Even], DaggerYes);

  //CK: For WilsonTM, comparison to CPS version
  //MooeeInv = 2 kappa g5theta(ctheta,-stheta)
  //kappa = 1/[2 sqrt( (m+4)^2 + eps^2 )]
  //ctheta = 2 (m+4) kappa
  //stheta = 2 eps kappa
  //g5theta(ctheta,stheta) = ctheta + i stheta g5


  // v1o = 1oo phi1
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
          this->axpby_ssp_proj(chi,this->beo[s],psi,   -this->ceo[s+1]  ,psi,s,s+1,Pplus);
          this->axpby_ssp_proj(chi,   1.0,chi,this->mass*this->ceo[this->Ls-1],psi,s,this->Ls-1,Pminus);
        } else if ( s==(this->Ls-1)) { 
          this->axpby_ssp_proj(chi,this->beo[s],psi,this->mass*this->ceo[0],psi,s,0,Pplus);
          this->axpby_ssp_proj(chi,1.0,chi,-this->ceo[s-1],psi,s,s-1,Pminus);
        } else {
          this->axpby_ssp_proj(chi,this->beo[s],psi,-this->ceo[s+1],psi,s,s+1,Pplus);
          this->axpby_ssp_proj(chi,1.0   ,chi,-this->ceo[s-1],psi,s,s-1,Pminus);
        }
      }
    } else { 

      // Assemble the 5d matrix
      for(int s=0;s<this->Ls;s++){
        if ( s==0 ) {
          //	chi = bs psi[s] + cs[s] psi[s+1}
          //    chi += -mass*cs[s] psi[s+1}
          this->axpby_ssp_proj(chi,this->beo[s],psi,-this->ceo[s],psi ,s, s+1,Pminus);
          this->axpby_ssp_proj(chi,1.0,chi,this->mass*this->ceo[s],psi,s,this->Ls-1,Pplus);
        } else if ( s==(this->Ls-1)) { 
          this->axpby_ssp_proj(chi,this->beo[s],psi,this->mass*this->ceo[s],psi,s,0,Pminus);
          this->axpby_ssp_proj(chi,1.0,chi,-this->ceo[s],psi,s,s-1,Pplus);
        } else {
          this->axpby_ssp_proj(chi,this->beo[s],psi,-this->ceo[s],psi,s,s+1,Pminus);
          this->axpby_ssp_proj(chi,1.0,chi,-this->ceo[s],psi,s,s-1,Pplus);
        }
      }
    }
  } else if(this->solver == DWF && this->precon_5d == 1) {
    // Booee is the identity matrix in this case.
    this->copy(chi, psi);
    return;
  } else if(this->solver == WilsonTM && this->precon_5d ==0){ //CK: I hope this is correct
    this->copy(chi, psi);
    return;
  } else {
    if ( this->isBoss() ) {
      printf("Booee: method not implemented for this fermion type / preconditioning type\n");
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
  mu = this->Mprec(x, t, y, 0, 1); // t = Mpc x  (y temp)
  this->Mprec(t, y, u, 1); //y = Mpc^dag t (u temp)
  if(! compute_min) {
    this->scale(y, -1.); //y=-y
    mu = -mu;
  }

  gnorm2 = this->axpy_norm(p, x, y, -mu); //p = -mu * x + y
  pnorm = sqrt(gnorm2);

  int i;
  for(i = 0; i < this->max_iter; ++i) {
    if(this->isBoss() && !me && i%100==0) {
      printf("bfm_evo::ritz iter = %6d gnorm2 = %17.10e, targ gnorm2 = %17.10e,  mu = %17.10e\n", i, gnorm2, stop_rsd, mu); 
    }

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
  low[mu] = send_neg ? 0 : lclx[mu] - 1; //pick out the slice at the boundary in the send direction
  high[mu] = low[mu] + 1;

  int block_size = 24 * lclx[4]; // s inner most

  const int hl[4] = {high[0] - low[0],
                     high[1] - low[1],
                     high[2] - low[2],
                     high[3] - low[3] };
  const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3]; //3-volume on surface (hl[mu]=1) [in units of blocks of size 24*Ls]
  const int vol4d = this->node_latt[0] * this->node_latt[1] * this->node_latt[2] * this->node_latt[3];

  int me, thrlen, throff;
  int work = hl_sites; if(cps::GJP.Gparity()) work *=2;
  this->thread_work(work, me, thrlen, throff);

  for(int i = 0; i < thrlen; ++i) {
    int x[4], flav;
    int tmp = i + throff;

    //For G-parity, fermion data blocks [each of size 24*Ls] increment in x,y,z,t,flav
    //Use similar mapping for surface volume, with flav changing slowest (offset 1 * surface 3-volume blocks)
    flav = tmp/hl_sites; tmp = tmp % hl_sites;

    x[0] = tmp % hl[0] + low[0]; tmp /= hl[0];
    x[1] = tmp % hl[1] + low[1]; tmp /= hl[1];
    x[2] = tmp % hl[2] + low[2]; tmp /= hl[2];
    x[3] = tmp % hl[3] + low[3];

    int off_4d = idx_4d(x, lclx);
    int off_3d = idx_4d_surf(x, lclx, mu);
        
    if(cps::GJP.Gparity()){
      //Implement G-parity flavour twist where appropriate. Note that the boundary sign on the boundary between C \bar{u}^T and d fields is implemented on the gauge links
      //here so we do not need to explicitly apply it to the communicated data.
      if(cps::GJP.Bc(mu) == cps::BND_CND_GPARITY && (send_neg && cps::GJP.NodeCoor(mu) == 0) || (!send_neg && cps::GJP.NodeCoor(mu) == cps::GJP.Nodes(mu)-1)   ){
	if(flav==0) memcpy(v3d + off_3d * block_size + hl_sites * block_size,  v4d + off_4d * block_size,  sizeof(Float) * block_size); //d -> CubarT buf
	else memcpy(v3d + off_3d * block_size,  v4d + off_4d * block_size + vol4d * block_size,  sizeof(Float) * block_size); //CubarT -> d buf
      }else{ //copy both flavours to their respective buffers
	memcpy(v3d + off_3d * block_size,  v4d + off_4d * block_size,  sizeof(Float) * block_size); //d -> d
	memcpy(v3d + off_3d * block_size + hl_sites * block_size,  v4d + off_4d * block_size + vol4d * block_size,  sizeof(Float) * block_size); //CubarT -> CubarT
      }
    }else{
      memcpy(v3d + off_3d * block_size,
	     v4d + off_4d * block_size,
	     sizeof(Float) * block_size);
    }
  }
}

// Calculate fermion force on a specific site, also do the
// summation over s direction.
//
// FIXME: need to add a line sum in s direction to support splitting
// in s direction.

//CK: v1p = v1[x+mu]
//    fermion vectors appear to be in CANONICAL ordering
template<class Float>
void bfm_evo<Float>::fforce_site(Float *mom, Float *gauge,
                                 Float *v1, Float *v1p,
                                 Float *v2, Float *v2p, int mu, Float coef,int gpf1_offset_p)
{
    Float t1[18], t2[18];

    if(cps::GJP.Gparity()) printf("flav 0\n");

    printf("v1: ");
    printf("%f %f ... %f",v1[0],v1[1],v1[24*this->Ls-1]);
    printf("\n");

    printf("v2: ");
    printf("%f %f ... %f",v2[0],v2[1],v2[24*this->Ls-1]);
    printf("\n");

    printf("v1p: ");
    printf("%f %f ... %f",v1p[0],v1p[1],v1p[24*this->Ls-1]);
    printf("\n");

    printf("v2p: ");
    printf("%f %f ... %f",v2p[0],v2p[1],v2p[24*this->Ls-1]);
    printf("\n");
    
    printf("gauge: %f %f ...%f\n",gauge[0],gauge[1],gauge[17]); 
    printf("mom: %f %f ...%f\n",mom[0],mom[1],mom[17]); 

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
    printf("Minus proj contrib: %f %f ... %f\n",t1[0],t1[1],t1[17]);
    printf("Plus proj contrib: %f %f ... %f\n",t2[0],t2[1],t2[17]);

    bfm_evo_aux::su3_add(t1, t2); //t1 -> t1 + t2
    bfm_evo_aux::mDotMEqual(t2, gauge, t1); //t2 -> gauge * t1
    bfm_evo_aux::trless_am(t2, -coef);

    printf("Traceless AHmat contrib: %f %f ... %f\n",t2[0],t2[1],t2[17]);

    if(cps::GJP.Gparity1fX()) for(int i=0;i<18;i++) t2[i]*=2.0;  //double latt testing, not production code

    bfm_evo_aux::su3_add(mom, t2);

    if(cps::GJP.Gparity()){
      //add force from second flavour
      Float t1_f1[18], t2_f1[18];
      const int vol4d = this->node_latt[0] * this->node_latt[1] * this->node_latt[2] * this->node_latt[3];
      const int f1_off = 24*this->Ls * vol4d; //f1 offset by 5d volume in this ordering scheme

      v1+=f1_off; v2+=f1_off; 
      v1p+=gpf1_offset_p; v2p+=gpf1_offset_p; //offset for 'plus' site depends on whether the data is stored in the buffer or the on-node vector

      printf("flav 1\n");

    printf("v1: ");
    printf("%f %f ... %f",v1[0],v1[1],v1[24*this->Ls-1]);
    printf("\n");

    printf("v2: ");
    printf("%f %f ... %f",v2[0],v2[1],v2[24*this->Ls-1]);
    printf("\n");

    printf("v1p: ");
    printf("%f %f ... %f",v1p[0],v1p[1],v1p[24*this->Ls-1]);
    printf("\n");

    printf("v2p: ");
    printf("%f %f ... %f",v2p[0],v2p[1],v2p[24*this->Ls-1]);
    printf("\n");

    Float *gauge_f1 = gauge + vol4d*18*4;
    Float *mom_f1 = mom + vol4d*18*4;
    printf("gauge: %f %f ...%f\n",gauge_f1[0],gauge_f1[1],gauge_f1[17]); 
    printf("mom: %f %f ...%f\n",mom_f1[0],mom_f1[1],mom_f1[17]); 

      switch(mu) {
      case 0:
	bfm_evo_aux::sprojTrXp(t1_f1, v1, v2p, this->Ls, 0, 0);   
	bfm_evo_aux::sprojTrXm(t2_f1, v2, v1p, this->Ls, 0, 0);
	break;
      case 1:
	bfm_evo_aux::sprojTrYp(t1_f1, v1, v2p, this->Ls, 0, 0);
	bfm_evo_aux::sprojTrYm(t2_f1, v2, v1p, this->Ls, 0, 0);
	break;
      case 2:
	bfm_evo_aux::sprojTrZp(t1_f1, v1, v2p, this->Ls, 0, 0);
	bfm_evo_aux::sprojTrZm(t2_f1, v2, v1p, this->Ls, 0, 0);
	break;
      default:
	bfm_evo_aux::sprojTrTp(t1_f1, v1, v2p, this->Ls, 0, 0);
	bfm_evo_aux::sprojTrTm(t2_f1, v2, v1p, this->Ls, 0, 0);
      }

      {
	cps::Matrix a; a.Trans(t1_f1); Float *aa = (Float*) &a[0];
	cps::Matrix b; b.Trans(t2_f1); Float *bb = (Float*) &b[0];
	printf("Minus proj contrib: %f %f ... %f\n",aa[0],aa[1],aa[17]);
	printf("Plus proj contrib: %f %f ... %f\n",bb[0],bb[1],bb[17]);
      }

      bfm_evo_aux::su3_add(t1_f1, t2_f1); //t1_f1 -> t1_f1 + t2_f1
      
      //set it up to use the f1 gauge field (sign*U*), such that the boundary sign comes free
      //this will need to be complex conjugated

      bfm_evo_aux::mStarDotMTransEqual(t2_f1, gauge_f1, t1_f1); // do  (U*)* t^T 
      bfm_evo_aux::trless_am(t2_f1, -coef);

      printf("Traceless AHmat contrib: %f %f ... %f\n",t2_f1[0],t2_f1[1],t2_f1[17]);

      bfm_evo_aux::su3_add(mom, t2_f1);

      //setup momentum for the second flavour
      for(int i=1;i<18;i+=2){ t2[i]*=-1; t2_f1[i]*=-1; } //mom[f1] is mom[f0]*
      bfm_evo_aux::su3_add(mom_f1, t2);
      bfm_evo_aux::su3_add(mom_f1, t2_f1);
    }

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
  --high[mu]; //exclude the site on the boundary

  int block_size = 24 * lclx[4];

  const int hl[4] = {high[0] - low[0],
                     high[1] - low[1],
                     high[2] - low[2],
                     high[3] - low[3] };
  const int hl_sites = hl[0] * hl[1] * hl[2] * hl[3];

  const int gparity_vp_off = block_size * lclx[0] * lclx[1] * lclx[2] * lclx[3]; //offset of second flavour (not used when G-parity is off)
  
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

    //testing
    int gx[4] = { x[0] + cps::GJP.XnodeSites()*cps::GJP.XnodeCoor(), 
		  x[1] + cps::GJP.YnodeSites()*cps::GJP.YnodeCoor(),
		  x[2] + cps::GJP.ZnodeSites()*cps::GJP.ZnodeCoor(),
		  x[3] + cps::GJP.TnodeSites()*cps::GJP.TnodeCoor() };

    if(cps::GJP.Gparity1fX()){
      int flav = 0;    
      if( gx[0] >= cps::GJP.XnodeSites()*cps::GJP.Xnodes()/2 ){ gx[0] -= cps::GJP.XnodeSites()*cps::GJP.Xnodes()/2; flav = 1; }
      printf("1f GP coord (%d %d %d %d) flav %d\n",gx[0],gx[1],gx[2],gx[3],flav);      
    }else if(cps::GJP.Gparity()){
      printf("2f GP coord (%d %d %d %d)\n",gx[0],gx[1],gx[2],gx[3]);
    }


    //Note fforce_site computes the force on this site from both flavours in the case of G-parity BCs
    this->fforce_site(mom + 18 * gid, gauge + 18 * gid,
                      v2 + fid, v2 + fidp,
                      v1 + fid, v1 + fidp, mu, coef, gparity_vp_off);
  }

  //GPARITY TESTING: COMPARE 1F AND 2F METHODS (NOT USED IN PRODUCTION CODE)
  if(cps::GJP.Gparity1fX() && me==0){ //use only first thread for this (does not need to be fast as it is only testing)
    this->thread_barrier();

    printf("Patching up 1f G-parity force\n");

    //want  p_0' = p_0 + delta p_0 + cconj(delta p_1)
    //      p_1' = p_1 + delta p_1 + cconj(delta p_0)
    //we did p_i' = p_i + 2 * delta p_i
    //and we know p_1 = cconj(p_0)
    //so we now do  p_0' = 0.5* p_0' + 0.5* cconj(p_1')
    //so we now do  p_1' = 0.5* p_1' + 0.5* cconj(p_0')
    //to fix this
    int momsz = 4*18*cps::GJP.VolNodeSites();
    Float *buf = (Float *)bfm_alloc(momsz * sizeof(Float) );
    for(int ii=0;ii<momsz;ii++) buf[ii] = 0.0;

    //Communicate \delta p from first half onto second half and vice versa
    Float *data_buf = mom;
    Float *send_buf = data_buf;
    Float *recv_buf = buf;

    if(cps::GJP.Xnodes()>1){
      //pass between nodes
      for(int i=0;i<cps::GJP.Xnodes()/2;i++){
	cps::getMinusData((Float *)recv_buf, (Float *)send_buf, momsz , 0);
	data_buf = recv_buf;
	recv_buf = send_buf;
	send_buf = data_buf;
      }
    }else{
      //shift mom[mu] field by xsites/2
      for(long i=0;i<cps::GJP.VolNodeSites();i++){
	//i = (x + Lx*(y+Ly*(z+Lz*t) ) )
	int x = i % cps::GJP.XnodeSites();
	int pos_rem = i/cps::GJP.XnodeSites(); //(y+Ly*(z+Lz*t)

	int x_from = (x + cps::GJP.XnodeSites()/2) % cps::GJP.XnodeSites();

	int i_from = 18*mu + 18*4*(x_from + cps::GJP.XnodeSites()*pos_rem);
	int i_to = 18*mu + 18*4*i;

	for(int j=0;j<18;j++) buf[i_to+j] = mom[i_from+j];
      }
      data_buf = buf;
    }
    for(int i=0;i<cps::GJP.VolNodeSites();i++){ //do fixup step
      int mat_off = 18*mu + 18*4*i;
      
      for(int j=0;j<18;j++){
	if(j%2==0) mom[mat_off+j] = mom[mat_off+j]/2.0 + data_buf[mat_off+j]/2.0;
	else mom[mat_off+j] = mom[mat_off+j]/2.0 - data_buf[mat_off+j]/2.0;
      }
    }
    bfm_free(buf);
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

  const int gparity_vp_off = block_size * hl_sites; //offset of second flavour (not used when G-parity is off)

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
                      v1 + fid, v1_s + fid_s, mu, coef,gparity_vp_off);
  }


  //GPARITY TESTING: COMPARE 1F AND 2F METHODS (NOT USED IN PRODUCTION CODE)
  if(cps::GJP.Gparity1fX() && me==0){ //use only first thread for this (does not need to be fast as it is only testing)
    this->thread_barrier();

    printf("Patching up 1f G-parity force\n");

    //want  p_0' = p_0 + delta p_0 + cconj(delta p_1)
    //      p_1' = p_1 + delta p_1 + cconj(delta p_0)
    //we did p_i' = p_i + 2 * delta p_i
    //and we know p_1 = cconj(p_0)
    //so we now do  p_0' = 0.5* p_0' + 0.5* cconj(p_1')
    //so we now do  p_1' = 0.5* p_1' + 0.5* cconj(p_0')
    //to fix this
    int momsz = 4*18*cps::GJP.VolNodeSites();
    Float *buf = (Float *)bfm_alloc(momsz * sizeof(Float) );
    for(int ii=0;ii<momsz;ii++) buf[ii] = 0.0;

    //Communicate \delta p from first half onto second half and vice versa
    Float *data_buf = mom;
    Float *send_buf = data_buf;
    Float *recv_buf = buf;

    if(cps::GJP.Xnodes()>1){
      //pass between nodes
      for(int i=0;i<cps::GJP.Xnodes()/2;i++){
	cps::getMinusData((Float *)recv_buf, (Float *)send_buf, momsz , 0);
	data_buf = recv_buf;
	recv_buf = send_buf;
	send_buf = data_buf;
      }
    }else{
      //shift mom[mu] field by xsites/2
      for(long i=0;i<cps::GJP.VolNodeSites();i++){
	//i = (x + Lx*(y+Ly*(z+Lz*t) ) )
	int x = i % cps::GJP.XnodeSites();
	int pos_rem = i/cps::GJP.XnodeSites(); //(y+Ly*(z+Lz*t)

	int x_from = (x + cps::GJP.XnodeSites()/2) % cps::GJP.XnodeSites();

	int i_from = 18*mu + 18*4*(x_from + cps::GJP.XnodeSites()*pos_rem);
	int i_to = 18*mu + 18*4*i;

	for(int j=0;j<18;j++) buf[i_to+j] = mom[i_from+j];
      }
      data_buf = buf;
    }
    for(int i=0;i<cps::GJP.VolNodeSites();i++){ //do fixup step
      int mat_off = 18*mu + 18*4*i;
      
      for(int j=0;j<18;j++){
	if(j%2==0) mom[mat_off+j] = mom[mat_off+j]/2.0 + data_buf[mat_off+j]/2.0;
	else mom[mat_off+j] = mom[mat_off+j]/2.0 - data_buf[mat_off+j]/2.0;
      }
    }
    bfm_free(buf);
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
    if(cps::GJP.Gparity()) surf_size[i] *=2; //2 flavours

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
  
  int fsize = vol_5d; if(cps::GJP.Gparity()) fsize*=2; 

  Float *v1f = this->threadedAllocFloat(fsize);
  Float *v2f = this->threadedAllocFloat(fsize);
  Float *sndbuf = this->threadedAllocFloat(surf_size_all);
  Float *rcvbuf = this->threadedAllocFloat(surf_size_all);

  this->thread_impexFermion_s(v1f, v1, 0);
  this->thread_impexFermion_s(v2f, v2, 0);

  for(int i = 0; i < 4; ++i) {
    this->copySendFrmData(sndbuf + surf_v1[i], v1f, i, true);
    this->copySendFrmData(sndbuf + surf_v2[i], v2f, i, true);
  }

  if(this->nthread <= 4) {
//#define DROPOUT_LT5THREADS
#ifdef DROPOUT_LT5THREADS
    if(!me) {
      printf("compute_force: Oops, at least 5 threads are needed.\n");
    }
    exit(-1);
#else
    //CK: We can do it with less than 5 threads, but less efficiently (so this will work on a cluster/laptop)
    if(me==0){ //Do comms on single thread
      for(int dir=0; dir<4; dir++)
	cps::getPlusData(rcvbuf + surf_v1[dir], sndbuf + surf_v1[dir],
			 surf_size[dir] * 2, dir);
    }    
    for(int i = 0; i < 4; ++i) {
      fforce_internal(mom, gauge, v1f, v2f, coef, i, me, this->nthread); //run over however many threads we have
    }
#endif
  }else{
    // Fused comm/internal force.
    //
    // The last 4 threads (typically 60-63) are used for
    // communication. All other threads (typically 0-59) are used to
    // calculate internal forces.
    
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

#if 0 //CK: in BFM, leaving them there!
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
#endif

template<typename Float>
void bfm_evo<Float>::deflate(Fermion_t out, Fermion_t in,
                             const multi1d<Fermion_t [2]> *evec,
                             const multi1d<Float> *eval,
                             int N)
{
  
  //CK: Why was this code disabled?? I have re-enabled it!
  //printf("void bfm_evo<Float>::deflate temporarily disabled\n");
  //exit(-1);

  if(N == 0 || evec == NULL || eval == NULL) {
    if(this->isBoss()) {
      printf("bfm_evo::deflate() must provide eigenvectors.\n");
    }
    exit(-1);
  }
  this->axpby(out, in, in, 0., 0.);
  //this->set_zero(out);
  for(int i = 0; i < N; ++i) {
    std::complex<double> dot = this->inner((*evec)[i][1], in);
//#ifdef BFM_GPARITY
#if 1
    this->caxpy(out, (*evec)[i][1], out, dot.real() / double((*eval)[i]),  dot.imag() / double((*eval)[i]) );
#else
    this->zaxpy(out, (*evec)[i][1], out, dot / double((*eval)[i]));
#endif
  }
}

// GCR, the matrix is preconditioned M.
template<class Float>
int bfm_evo<Float>::gcr_M(Fermion_t sol, Fermion_t src)
{
  printf("int bfm_evo<Float>::gcr_M temporarily disabled");
  exit(-1);
  
#if 0
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
#endif
}

// GMRES(m), we restart after m iterations.
template<class Float>
int bfm_evo<Float>::gmres_M(Fermion_t sol, Fermion_t src, const int m)
{
  printf("int bfm_evo<Float>::gmres_M temporarily disabled\n");
  exit(-1);

#if 0
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
#endif
}

#endif
