// Is STL available ?
#include <vector>
#include <config.h>

CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOp class Low Mode Approximation 

*/
  CPS_END_NAMESPACE
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/checksum.h>
#include <comms/glb.h>
#include <math.h>
#include <stdio.h>
#include <comms/sysfunc_cps.h>
#include <util/eigen_container.h>
#include<util/time_cps.h>
#include<util/conversion.h>
#ifdef USE_BLAS
#include <util/qblas_extend.h>
#endif
#ifdef USE_QUDA
#include "cublas_v2.h"
#include <cuComplex.h>
#endif
  CPS_START_NAMESPACE
//void cDot (std::vector < Float * >)
#if 0
template < FLOAT1, FLOAT2 >
  void compDotProduct (std::vector < Float > &result,
                       const std::vector < FLOATA * >a,
                       const std::vector < FLOATB * >b, int len)
{

  int a_size = a.size ();
  int b_size = b.size ();
  result.resize (2 * a_size * b_size, 0.);
//#pragma omp parallel for reduction(+:re,im)
  for (int i = 0; i < a_size; i += a_step)
    for (int j = 0; j < b_size; j += b_step)
      for (int k = 0; k < len; j += c_step) {
#pragma omp parallel for
        for (int kk = 0; (kk < c_step) && ((k + kk) < c_size); kk++) {
          std::vector < Float > re (a_step * b_step, 0);
          std::vector < Float > in (a_step * b_step, 0);
          size_t ind = 2 * (k + kk);
          for (int ii = 0; (ii < a_step) && ((i + ii) < a_size); ii++) {
            FLOATA *a_p = a[i + ii];
            for (int jj = 0; (jj < b_step) && ((j + jj) < b_size); jj++) {
              FLOATB *b_p = b[j + jj];
              re[ii + a_step * jj] +=
                *(a_p + ind) * *(b_p + ind) + *(a_p + ind + 1) * *(b_p + ind +
                                                                   1);
              im[ii + a_step * jj] +=
                *(a_p + ind) * *(b_p + ind + 1) - *(a_p + ind + 1) * *(b_p +
                                                                       ind);
//       im+= a[i] * b[i+1] - a[i+1] * b[i];       // imag part
            }
          }
          for (int ii = 0; (ii < a_step) && ((i + ii) < a_size); ii++)
            for (int jj = 0; (jj < b_step) && ((j + jj) < b_size); jj++)
#pragma omp critical
            {
              int ind = (i + ii) + a_size * (j + jj);
              result[2 * ind] += re[ii + a_step * jj];
              result[2 * ind + 1] += im[ii + a_step * jj];
            }
        }
      }

}
#endif

#undef USE_CUDA_MM
#ifdef USE_CUDA_MM
#define CUDA_MALLOC cudaMallocManaged
#else
#define CUDA_MALLOC cudaMalloc
#endif
//#define CUDA_COMPLEX cuDoubleComplex
//#define ZAXPY cublasZaxpy
//#define ZDOTC cublasZdotc
#define CUDA_COMPLEX cuComplex
#define MAKE_COMPLEX make_cuFloatComplex
#define ZAXPY cublasCaxpy
#define ZDOTC cublasCdotc
static int src_max = 16;
static int evec_max = 64;
static std::vector < CUDA_COMPLEX * >h_src (src_max, NULL), h_dest (src_max,
                                                                    NULL),
h_evec (evec_max, NULL);
static std::vector < CUDA_COMPLEX * >d_src (src_max, NULL), d_dest (src_max,
                                                                    NULL),
d_evec (evec_max, NULL);
void cublasGS (std::vector < Vector * >&dest, std::vector < Vector * >&src,
               size_t evec_len, Float ** evec, size_t f_len, Float * eval)
{
  const char *fname = "cublasGS()";
  size_t numElem = f_len / 2;
  size_t size = numElem * sizeof (CUDA_COMPLEX);
  assert (dest.size () == src.size ());
  int src_len = dest.size ();
  assert (src_len <= src_max);
  assert (evec_len <= evec_max);
  CUDA_COMPLEX *z_cublas, *z_cublas2;
  Float *z;

  Float dtime = time_elapse ();
#ifndef USE_CUDA_MM
  for (int i = 0; i < src_len; i++)
    if (!h_src[i])
      cudaMallocHost ((CUDA_COMPLEX **) & (h_src[i]), size);
  for (int i = 0; i < src_len; i++)
    if (!h_dest[i])
      cudaMallocHost ((CUDA_COMPLEX **) & (h_dest[i]), size);
  for (int i = 0; i < evec_len; i++)
    if (!h_evec[i])
      cudaMallocHost ((CUDA_COMPLEX **) & (h_evec[i]), size);
#endif

  cudaMallocHost ((CUDA_COMPLEX **) & z_cublas,
                  src_len * evec_len * sizeof (CUDA_COMPLEX));
  cudaMallocHost ((CUDA_COMPLEX **) & z_cublas2,
                  src_len * evec_len * sizeof (CUDA_COMPLEX));
  z = (Float *) smalloc (2 * src_len * evec_len * sizeof (Float));

  dtime = time_elapse ();
  VRB.Result ("", fname, "cudaMallocHost(%d %d %d ) %e seconds\n",
              src_len, evec_len, f_len, dtime);

  for (int i = 0; i < src_len; i++)
    if (!d_src[i])
      cudaMallocManaged ((void **) &(d_src[i]), size);
  for (int i = 0; i < src_len; i++)
    if (!d_dest[i])
      cudaMallocManaged ((void **) &(d_dest[i]), size);
  for (int i = 0; i < evec_len; i++)
    if (!d_evec[i])
      cudaMallocManaged ((void **) &(d_evec[i]), size);

  dtime = time_elapse ();
  VRB.Result ("", fname, "cudaMalloc(%d %d %d ) %e seconds\n",
              src_len, evec_len, f_len, dtime);

  for (int j = 0; j < src_len; j++) {
    Float *src_p = (Float *) src[j];
    Float *dest_p = (Float *) dest[j];
#pragma omp parallel for
    for (size_t i = 0; i < numElem; i++) {
#ifdef USE_CUDA_MM
      *(d_src[j] + i) = MAKE_COMPLEX (src_p[i * 2], src_p[i * 2 + 1]);
      *(d_dest[j] + i) = MAKE_COMPLEX (dest_p[i * 2], dest_p[i * 2 + 1]);
#else
      *(h_src[j] + i) = MAKE_COMPLEX (src_p[i * 2], src_p[i * 2 + 1]);
      *(h_dest[j] + i) = MAKE_COMPLEX (dest_p[i * 2], dest_p[i * 2 + 1]);
#endif
    }
  }

  for (int j = 0; j < evec_len; j++) {
    float *evec_p = (float *) evec[j];
#pragma omp parallel for
    for (size_t i = 0; i < numElem; i++) {
#ifdef USE_CUDA_MM
      *(d_evec[j] + i) = MAKE_COMPLEX (evec_p[i * 2], evec_p[i * 2 + 1]);
#else
      *(h_evec[j] + i) = MAKE_COMPLEX (evec_p[i * 2], evec_p[i * 2 + 1]);
#endif
    }
  }
  dtime = time_elapse ();
  VRB.Result ("", fname, "make_cuComplex(%d %d %d ) %e seconds\n",
              src_len, evec_len, f_len, dtime);


  cublasHandle_t handle;
  cublasCreate (&handle);
  cublasStatus_t stat;

#ifndef USE_CUDA_MM
  for (int i = 0; i < src_len; i++)
    cudaMemcpy (d_src[i], h_src[i], size, cudaMemcpyHostToDevice);
  for (int i = 0; i < src_len; i++)
    cudaMemcpy (d_dest[i], h_dest[i], size, cudaMemcpyHostToDevice);
  for (int i = 0; i < evec_len; i++)
    cudaMemcpy (d_evec[i], h_evec[i], size, cudaMemcpyHostToDevice);
#endif
  dtime = time_elapse ();
  VRB.Result ("", fname, "cudaMemcpy(%d %d %d ) %e seconds\n", src_len,
              evec_len, f_len, dtime);
  for (int gs = 0; gs < 1; gs++) {
    dtime = time_elapse ();
    for (int j = 0; j < evec_len; j++) {
      for (int i = 0; i < src_len; i++) {
        stat =
          ZDOTC (handle, numElem, d_evec[j], 1, d_src[i], 1,
                 z_cublas + (i * evec_len + j));
      }
    }
    for (int i = 0; i < evec_len * src_len; i++) {
      z[i * 2] = cuCrealf (z_cublas[i]);
      z[i * 2 + 1] = cuCimagf (z_cublas[i]);
    }
    glb_sum (z, src_len * evec_len * 2);
    for (int j = 0; j < evec_len; j++) {
      Float fac = 1. / eval[j];
      VRB.Debug ("", fname, "fac[%d]=%e\n", j, fac);
      for (int i = 0; i < src_len; i++) {
        z_cublas2[i * evec_len + j] =
//          MAKE_COMPLEX (fac * cuCrealf (z_cublas[i * evec_len + j]),
//                                fac * cuCimagf (z_cublas[i * evec_len + j]));
          MAKE_COMPLEX (fac * z[2 * (i * evec_len + j)],
                        fac * z[1 + 2 * (i * evec_len + j)]);
        VRB.Debug ("", fname, "coef[%d][%d]=%e %e\n", i, j,
                    fac * z[2 * (i * evec_len + j)],
                    fac * z[1 + 2 * (i * evec_len + j)]);
        stat =
          ZAXPY (handle, numElem, z_cublas2 + (i * evec_len + j),
                 d_evec[j], 1, d_dest[i], 1);
      }
    }
    dtime = time_elapse ();
    VRB.Result ("", fname, "cublasZdotc+Zaxpy(%d %d %d ) %e seconds\n", src_len,
                evec_len, f_len, dtime);
    double *z_p = (double *) z_cublas;
    if (0)
      for (int i = 0; i < src_len; i++)
        for (int j = 0; j < evec_len; j++) {
          VRB.Debug ("", fname, "cublasZdotc[%d][%d]  = (%e %e) \n", i, j,
                     cuCrealf (z_cublas[i * evec_len + j]),
                     cuCimagf (z_cublas[i * evec_len + j]));
//     z_cublas[i*evec_len+j]=make_cuDoubleComplex(-cuCreal(z_cublas[i*evec_len+j]),-cuCimagf(z_cublas[i*evec_len+j]));
        }
  }

  if (0) {
    for (int j = 0; j < evec_len; j++)
      for (int i = 0; i < src_len; i++)
        stat =
          ZDOTC (handle, numElem, d_evec[j], 1, d_src[i], 1,
                 z_cublas + (i * evec_len + j));

    for (int i = 0; i < src_len; i++)
      for (int j = 0; j < evec_len; j++) {
        VRB.Debug ("", fname, "cublasZdotc[%d][%d] (after) = (%e %e) \n", i, j,
                   cuCrealf (z_cublas[i * evec_len + j]),
                   cuCimagf (z_cublas[i * evec_len + j]));
      }
  }
//  dtime = time_elapse ();
#ifndef USE_CUDA_MM
  for (int i = 0; i < src_len; i++)
    cudaMemcpy (h_dest[i], d_dest[i], size, cudaMemcpyDeviceToHost);
#endif

  for (int j = 0; j < src_len; j++) {
    Float *dest_p = (Float *) dest[j];
#pragma omp parallel for
    for (size_t i = 0; i < numElem; i++) {
//    *(h_src[j]+i)=make_cuDoubleComplex(src_p[i*2],src_p[i*2+1]);
#ifdef USE_CUDA_MM
      dest_p[i * 2] = cuCrealf (*(d_dest[j] + i));
      dest_p[i * 2 + 1] = cuCimagf (*(d_dest[j] + i));
#else
      dest_p[i * 2] = cuCrealf (*(h_dest[j] + i));
      dest_p[i * 2 + 1] = cuCimagf (*(h_dest[j] + i));
#endif
    }
  }
  dtime = time_elapse ();
  VRB.Result ("", fname, "cudaMemcpy+copy(%d %d %d ) %e seconds\n", src_len,
              evec_len, f_len, dtime);

  for (int j = 0; j < src_len; j++) {
//    cudaFree (d_src[j]);
//    cudaFree (d_dest[j]);
  }
  for (int j = 0; j < evec_len; j++) {
//    cudaFree (d_evec[j]);
  }
  dtime = time_elapse ();
  VRB.Result ("", fname, "cudaFree(%d %d %d ) %e seconds\n", src_len,
              evec_len, f_len, dtime);

#ifndef USE_CUDA_MM
  for (int j = 0; j < src_len; j++) {
//    cudaFreeHost (h_src[j]);
//    cudaFreeHost (h_dest[j]);
  }
  for (int j = 0; j < evec_len; j++) {
//    cudaFreeHost (h_evec[j]);
  }
#endif
  cudaFreeHost (z_cublas);
  cudaFreeHost (z_cublas2);
  sfree (z);
  dtime = time_elapse ();
  VRB.Result ("", fname, "cudaFreeHost(%d %d %d ) %e seconds\n", src_len,
              evec_len, f_len, dtime);
}

#undef CUDA_MALLOC
#undef CUDA_COMPLEX
#undef MAKE_COMPLEX
#undef ZAXPY

void movefloattoFloat (Float * out, float *in, size_t f_size);
int DiracOp::InvLowModeApprox (std::vector < Vector * >&out,
                               std::vector < Vector * >&in,
                               char *fname_eig_root, int neig)
{


  time_elapse ();
  size_t f_size_cb;             // Node checkerboard size of the fermion field
//  int itr;                      // Current number of CG iterations
  int max_itr;                  // Max number of CG iterations
  Float stp_cnd;                // Stop if residual^2 <= stp_cnd
  Float res_norm_sq_prv;        // The previous step |residual|^2
  Float res_norm_sq_cur;        // The current step |residual|^2 
  Float a;
  Float b;
  Float d;
//  int i, j;
  assert (out.size () == in.size ());
  int len = out.size ();
  const char *fname = "InvLowModeApprox(V*,V*,C,I,F*)";
  VRB.Result (cname, fname, "len=%d\n", len);
  static Timer timer (cname, fname);
  timer.start ();
//  ERR.General(cname,fname,"Broken without BLAS\n");
  Float mdagm_time = 0.;
  Float gsum_time = 0.;
  Float linalg_time = 0.;
//------------------------------------------------------------------
// Initializations
// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------
  if (lat.Fclass () == F_CLASS_CLOVER) {
    f_size_cb = (size_t) GJP.VolNodeSites () * lat.FsiteSize () / 2;
  } else {
    f_size_cb =
      (size_t) GJP.VolNodeSites () * lat.FsiteSize () / (lat.FchkbEvl () + 1);
  }
  assert (f_size_cb == lat.half_size);
  const int n_fields = GJP.SnodeSites ();       //   *nk ; 
  const size_t f_size_per_site =
    lat.FsiteSize () / GJP.SnodeSites () / (lat.FchkbEvl () + 1);
  char fname_eig_root_bc[1024];
  snprintf (fname_eig_root_bc, 1024, "%s", fname_eig_root);
//         GJP.Bc(0),GJP.Bc(1),GJP.Bc(2),GJP.Bc(3));
  VRB.Debug (cname, fname, "fname=%s\n", fname_eig_root_bc);
  // search for eigen cache
  EigenCache *ecache;
  if ((ecache = EigenCacheListSearch (fname_eig_root_bc, neig)) == 0) {

    ERR.General (cname, fname,
                 "Eigenvector cache does not exist: neig %d name %s \n",
                 neig, fname_eig_root_bc);
//    printf("cache_name=%s neig=%d ecache=%p\n",fname_eig_root_bc,neig,ecache);
    //ecache = new EigenCache();  
    //EigenCacheList. push_back( ecache );
  } else
    VRB.Result (cname, fname, "returned %s %d %d\n", ecache->Name (),
                ecache->Neig (), ecache->NeigCoarse ());
//  exit(-42);
  Float mass = dirac_arg->mass;
  Float *eval = ecache->evals.data ();
  //
  // make eigen values squared if needed
  //
  switch (dirac_arg->RitzMatOper) {
  case MAT_HERM:
  case MATPC_HERM:
    for (int j = 0; j < neig; ++j)
      eval[j] = eval[j] * eval[j];
    break;
  case MATPCDAG_MATPC:
  case NEG_MATPCDAG_MATPC:
  case MATDAG_MAT:
  case NEG_MATDAG_MAT:
  case MATDAG_MAT_NORM:
  case NEG_MATDAG_MAT_NORM:
  case MATPCDAG_MATPC_SHIFT:
    break;
  default:
    ERR.General (cname, fname, "RiztMatType %d not implemented\n",
                 dirac_arg->RitzMatOper);
  }

//------------------------------------------------------------------
// Set the source vector pointer
//------------------------------------------------------------------
// Set the solution vector pointer
//------------------------------------------------------------------
  time_elapse ();
  if (0) {
    std::vector < Float > apb (2 * len * neig, 0.);
    Float *apb_p = apb.data ();
    std::vector < float *>evecs;
    std::vector < Float * >sol_p;
    int num = neig;
    if (num > len)
      num = len;
    if (num) {
      for (int iev = 0; iev < num; ++iev) {
        evecs.push_back ((float *) ecache->evecs[iev]);
        VRB.Flow (cname, fname, "evecs[%d]=%p\n", iev, evecs[iev]);
      }
      for (int iev = 0; iev < len; ++iev) {
        sol_p.push_back ((Float *) in[iev]);
        VRB.Flow (cname, fname, "sol=%p\n", in[iev]);
      }

      time_elapse ();
      compDotProduct (apb, evecs, sol_p, f_size_cb / 2);        // 2 for complex
      print_flops (fname, "compDotProduct",
                   len * (8. / 2.) * f_size_cb * num, time_elapse ());
      glb_sum (apb_p, apb.size ());
    }
  }

  if(  ecache->Neig() <neig) {
    std::vector < Float * >sol_f;
    std::vector < Float * >src_f;
    src_f.resize (len);
    sol_f.resize (len);
    for (int i = 0; i < len; i++) {
      sol_f[i] = (Float *) out[i];
      src_f[i] = (Float *) in[i];
    }
    std::vector < Complex > coef;
// printf 
//    VRB.Result(cname,fname,"ecache=%p neig=%d neig_b=%d\n",ecache,ecache->Neig(),ecache->NeigCoarse());
// old implementation
    
    if (0) {
      for (int i = 0; i < len; i++) {
        out[i]->VecZero (f_size_cb);
      }
      ecache->blockProj (sol_f, src_f, 0, neig, coef);
      for (int i = 0; i < len; i++) {
        VRB.Result (cname, fname, "%d: src=%e norm(sol,ecache_proj) =%e\n", i,
                    in[i]->NormSqGlbSum (f_size_cb),
                    out[i]->NormSqGlbSum (f_size_cb));
      }
      for (int i = 0; i < neig; i++)
        for (int j = 0; j < len; j++) {
          VRB.Flow (cname, fname, "coef[%d][%d] = %e %e\n", i, j,
                    coef[j + len * i].real (), coef[j + len * i].imag ());
        }
    } else {
      for (int i = 0; i < len; i++) {
        out[i]->VecZero (f_size_cb);
      }
      ecache->blockProj2 (sol_f, src_f, 0, neig, coef);
      for (int i = 0; i < len; i++) {
        VRB.Result (cname, fname, "%d: src=%e norm(sol,ecache_proj2) =%e\n", i,
                    in[i]->NormSqGlbSum (f_size_cb),
                    out[i]->NormSqGlbSum (f_size_cb));
      }
      for (int i = 0; i < neig; i++)
        for (int j = 0; j < len; j++) {
          VRB.Flow (cname, fname, "coef[%d][%d] = %e %e\n", i, j,
                    coef[j + len * i].real (), coef[j + len * i].imag ());
        }
    }

  } else {
#ifdef USE_QUDA
    Float dtime;
    if (1) {
      dtime = -dclock ();
      for (int i = 0; i < len; i++) {
        out[i]->VecZero (f_size_cb);
      }
      int index = 0;
      while (index < neig) {
        int evec_len = evec_max;
        if ((neig - index) < evec_len)
          evec_len = neig - index;
        cublasGS (out, in, evec_len, &(ecache->evecs[index]), f_size_cb,
                  &(eval[index]));
        index += evec_len;
      }
      dtime += dclock ();
      print_flops (fname, "cublas", len * f_size_cb * 7 * neig, dtime);
      for (int i = 0; i < len; i++) {
        VRB.Result (cname, fname, "norm(sol,cublas) %d =%e\n", i,
                    out[i]->NormSqGlbSum (f_size_cb));
      }
    } else
#endif
    {
      time_elapse ();
      dtime = -dclock ();
      for (int i = 0; i < len; i++) {
        out[i]->VecZero (f_size_cb);
      }

      std::vector < Complex > z (len);
      Float *evecFloat = (Float *) smalloc (f_size_cb * sizeof (Float));
      Vector *evecVec = (Vector *) evecFloat;
      for (int iev = 0; iev < neig; ++iev) {

        movefloattoFloat (evecFloat, (float *) ecache->evecs[iev], f_size_cb);
        for (int i = 0; i < len; i++) {
          IFloat *src_tmp = (IFloat *) in[i];
#ifndef USE_BLAS
          z[i] = evecVec->CompDotProductGlbSum (in[i], f_size_cb);
#else
          Complex z;
          cblas_zdotc_sub (f_size_cb / 2,
                           (double *) evecFloat, 1, (double *) in[i], 1,
                           (double *) &z);
          glb_sum ((double *) &z);
          glb_sum ((double *) &z + 1);
#endif
          VRB.Debug (cname, fname,
                     "evec %d: %e %e src: %e %e z: %e %e eval=%e\n", iev,
                     *(evecFloat), *(evecFloat + 1), *(src_tmp),
                     *(src_tmp + 1), z[i].real (), z[i].imag (), eval[iev]);
          z[i] /= eval[iev];
        }
        for (int i = 0; i < len; i++) {
#ifndef USE_BLAS
          out[i]->CTimesV1PlusV2 (z[i], evecVec, out[i], f_size_cb);
#else
          cblas_zaxpy (f_size_cb / 2, (double *) (z.data () + i),
                       (Float *) evecFloat, 1, (Float *) out[i], 1);
#endif
//    print_flops("inv_lowmode_approx","ctimesv1pv2", f_size_cb*4, time_elapse());
        }
      }
      dtime += dclock ();
      print_flops (fname, "cpu", len * f_size_cb * 7 * neig, dtime);
      for (int i = 0; i < len; i++) {
        VRB.Result (cname, fname, "norm(sol) %d =%e\n", i,
                    out[i]->NormSqGlbSum (f_size_cb));
      }

      sfree (evecFloat);
//  *true_res = 0.0;              // could compute the residue norm here, but saving a dirac operator for now
    }
  }
  timer.stop (true);
//  exit (-43);
  return neig;
}


void DiracOp::InvLowModeProj (Vector * in, char *fname_eig_root, int neig)
{

  time_elapse ();
  size_t f_size_cb;             // Node checkerboard size of the fermion field
  int itr;                      // Current number of CG iterations
  int max_itr;                  // Max number of CG iterations
  Float stp_cnd;                // Stop if residual^2 <= stp_cnd
  Float a;
  Float b;
  Float d;
  int i, j;
  char *fname = "InvLowModeApprox(V*,V*,C,I,F*)";
//  ERR.General(cname,fname,"Broken without BLAS\n");
//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------
// Set the source vector pointer
//------------------------------------------------------------------
  Vector *src = in;
  IFloat *src_tmp = (IFloat *) src;
// Set the node checkerboard size of the fermion field
//------------------------------------------------------------------
  if (lat.Fclass () == F_CLASS_CLOVER) {
    f_size_cb = GJP.VolNodeSites () * lat.FsiteSize () / 2;
  } else {
    f_size_cb = GJP.VolNodeSites () * lat.FsiteSize () / (lat.FchkbEvl () + 1);
  }
  assert (f_size_cb == lat.half_size);
  const int n_fields = GJP.SnodeSites ();       //   *nk ; 
  const size_t f_size_per_site =
    lat.FsiteSize () / GJP.SnodeSites () / (lat.FchkbEvl () + 1);
  char fname_eig_root_bc[1024];
  snprintf (fname_eig_root_bc, 1024, "%s.bc%d%d%d%d",
            fname_eig_root, GJP.Bc (0), GJP.Bc (1), GJP.Bc (2), GJP.Bc (3));
  // search for eigen cache
  EigenCache *ecache;
  if ((ecache = EigenCacheListSearch (fname_eig_root_bc, neig)) == 0) {
    ecache = new EigenCache ();
    EigenCacheList.push_back (ecache);
  }

  Float mass = dirac_arg->mass;
#if 0
  EigenContainer eigcon (lat, fname_eig_root_bc, neig,
                         f_size_per_site, n_fields, ecache);
  Float *eval = eigcon.load_eval ();
#else
  Float *eval = ecache->evals.data ();
#endif
  //
  // make eigen values squared if needed
  //
  switch (dirac_arg->RitzMatOper) {
  case MAT_HERM:
  case MATPC_HERM:
    for (int i = 0; i < neig; ++i)
      eval[i] = eval[i] * eval[i];
    break;
  case MATPCDAG_MATPC:
  case NEG_MATPCDAG_MATPC:
  case MATDAG_MAT:
  case NEG_MATDAG_MAT:
  case MATDAG_MAT_NORM:
  case NEG_MATDAG_MAT_NORM:
  case MATPCDAG_MATPC_SHIFT:
    break;
  default:
    ERR.General (cname, fname, "RiztMatType %d not implemented\n",
                 dirac_arg->RitzMatOper);
  }
  Float *evecFloat = (Float *) smalloc (f_size_cb * sizeof (Float));
  Vector *evecVec = (Vector *) evecFloat;
  FILE *fp = Fopen ("LowmodeSrcProj.dat", "o");
  for (int iev = 0; iev < neig; ++iev) {
    movefloattoFloat (evecFloat, (float *) ecache->evecs[iev], f_size_cb);
    //time_elapse();
//    Vector *evec = eigcon.nev_load (iev);
    //print_time("inv_lowmode_approx","loading", time_elapse());
#ifndef USE_BLAS
    Complex z = evecVec->CompDotProductGlbSum (src, f_size_cb);
#else
    Complex z;
    cblas_zdotc_sub (f_size_cb / 2,
                     (double *) evecVec, 1, (double *) src, 1, (double *) &z);
    glb_sum ((double *) &z);
    glb_sum ((double *) &z + 1);
#endif
    z /= eval[iev];
    if (!UniqueID ())
      fprintf (fp, "%d %.14e %.14e\n", iev, z.real (), z.imag ());
  }
  fclose (fp);
  return;
}

CPS_END_NAMESPACE
