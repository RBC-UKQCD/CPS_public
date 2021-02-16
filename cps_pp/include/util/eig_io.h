#ifndef INCLUDED_EIG_IO_H
#define INCLUDED_EIG_IO_H

#define _FILE_OFFSET_BITS 64
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include <complex>
#include <map>
#include <vector>
#include <memory.h>
#include <iostream>

#include <sys/stat.h>
#include <comms/sysfunc_cps.h>
#include <comms/glb.h>
#include <util/sumarray.h>
#include <util/time_cps.h>
//#include <util/eigen_container.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <unistd.h>


#include <zlib.h>
static uint32_t crc32_loop (uint32_t previousCrc32, unsigned char *data,
                            size_t len)
{

  // CL: crc32 of zlib was incorrect for very large sizes, so do it block-wise
  uint32_t crc = previousCrc32;
  size_t blk = 0;
  size_t step = (size_t) 1024 * 1024 * 1024;
  while (len > step) {
    crc = crc32 (crc, &data[blk], step);
    blk += step;
    len -= step;
  }

  if (len > 0)
    crc = crc32 (crc, &data[blk], len);
  return crc;

}





CPS_START_NAMESPACE
#if 1
  struct _evc_meta_
{
  int s[5];
  int b[5];
  int nkeep;
  int nkeep_single;

  // derived
  int nb[5];
  int blocks;

  int neig;

  int index;

    std::vector < uint32_t > crc32_header;

  int FP16_COEF_EXP_SHARE_FLOATS;
//    bool big_endian;
};

#if 1
struct evec_write
{
//  int s[5];
  int b[5];
  int nkeep;
  int nkeep_single;

  int findex;
//  int filesperdir;
  int n_dir;
  int bigendian;

  int vrb_nkeep_res;
  int vrb_evec_res;

  int concur;

  // derived
//  int nb[5];
//  int blocks;
//  int prec;
};
#endif

class EigenCache;               //forward declaration


class EvecReader
{
private:
  const char *cname;
  int nb_per_node;
  size_t f_size;
//    size_t f_size_block;

public:
  const int n_cycle = 32;
  bool machine_is_little_endian;
  bool bigendian;               // read/write data back in big endian
  _evc_meta_ args;
  bool crc32_checked;
  uint32_t nprocessors;
  uint32_t nfile;
//      std::vector < int >nn;
    EvecReader ():cname ("EvecReader"), bigendian (false), crc32_checked (false), nprocessors (1), nfile (1)    // fixed to little endian
  {
    machine_is_little_endian = machine_endian ();
    VRB.Result (cname, "EvecReader()", "machine_is_little_endian=%d\n",
                machine_is_little_endian);
  };
  ~EvecReader () {
  };
  typedef float OPT;
  int nthreads;

  static const char *header;
  std::vector < double >evals;

  size_t vol4d, vol5d;
  size_t f_size_block;



  std::vector < std::vector < OPT > >block_data;
  std::vector < std::vector < OPT > >block_data_ortho;
  std::vector < std::vector < OPT > >block_coef;

  int machine_endian ()
  {
    char endian_c[4] = { 1, 0, 0, 0 };
    uint32_t *endian_i = (uint32_t *) endian_c;
    if (*endian_i == 0x1)
//        printf("little endian\n");
      return 1;
    else
//        printf("big endian\n");
      return 0;

  }



  template < class M > int sumArray (M * vs, const long n_elem)
  {
    // M can be double or long
    int status = 0;
    status = glb_sum (vs, n_elem);
    return status;
  }

  void fix_short_endian (unsigned short *dest, int nshorts)
  {
    if (nshorts > 128)
      VRB.Result (cname, "fix_short_endian", "nshorts=%d\n", nshorts);
    if ((bigendian && machine_is_little_endian) ||      // written for readability
        (!bigendian && !machine_is_little_endian))
      for (int i = 0; i < nshorts; i++) {
        char *c1 = (char *) &dest[i];
        char tmp = c1[0];
        c1[0] = c1[1];
        c1[1] = tmp;
      }
  }

  void fix_float_endian (float *dest, int nfloats)
  {
//      int n_endian_test = 1;
    //bool machine_is_little_endian = *(char *)&n_endian_test == 1;
    if (nfloats > 128)
      VRB.Debug (cname, "fix_float_endian", "nfloats=%d\n", nfloats);

    if ((bigendian && machine_is_little_endian) ||      // written for readability
        (!bigendian && !machine_is_little_endian)) {
#pragma omp parallel for
      for (int i = 0; i < nfloats; i++) {
        float before = dest[i];
        char *c = (char *) &dest[i];
        char tmp;
        int j;
        for (j = 0; j < 2; j++) {
          tmp = c[j];
          c[j] = c[3 - j];
          c[3 - j] = tmp;
        }
        float after = dest[i];
        if (!i)
          VRB.Result (cname, "fix_float_endian",
                      "fix_float_endian: %g ->%g\n", before, after);
      }
    }

  }


  int get_bfm_index (int *pos, int co, int *s_l)
  {

    int ls = s_l[4];
//      int vol_4d_oo = vol4d / 2;
//      int vol_5d = vol_4d_oo * ls;

    int NtHalf = s_l[3] / 2;
    int simd_coor = pos[3] / NtHalf;
    int regu_coor =
      (pos[0] +
       s_l[0] * (pos[1] + s_l[1] * (pos[2] + s_l[2] * (pos[3] % NtHalf)))) / 2;
//      int regu_vol = vol_4d_oo / 2;

    return +regu_coor * ls * 48 + pos[4] * 48 + co * 4 + simd_coor * 2;
  }

  size_t get_cps_index (int *pos, int co, int *s_l)
  {

#ifdef USE_BFM
    return get_bfm_index (pos, co, s_l);
#else
//      int ls = s_l[4];
//      int vol_4d_oo = vol4d / 2;
//      int vol_5d = vol_4d_oo * ls;

    int regu_coor = (pos[0] + s_l[0] *
                     (pos[1] + s_l[1] *
                      (pos[2] + s_l[2] * (pos[3] + s_l[3] * pos[4])))) / 2;

    return ((regu_coor) * 12 + co) * 2;
#endif
  }

  inline void index_to_pos (int i, int *pos, int *latt)
  {
    int d;
    for (d = 0; d < 5; d++) {
      pos[d] = i % latt[d];
      i /= latt[d];
    }
  }

  inline int pos_to_index (int *pos, int *latt)
  {
    return pos[0] + latt[0] * (pos[1] +
                               latt[1] * (pos[2] +
                                          latt[2] * (pos[3] +
                                                     latt[3] * pos[4])));
  }


  inline void pos_to_blocked_pos (int *pos, int *pos_in_block, int *block_coor)
  {
    int d;
    for (d = 0; d < 5; d++) {
      block_coor[d] = pos[d] / args.b[d];
      pos_in_block[d] = pos[d] - block_coor[d] * args.b[d];
    }
  }

  template < class T >
    void caxpy_single (T * res, std::complex < T > ca, T * x, T * y,
                       size_t f_size)
  {
    std::complex < T > *cx = (std::complex < T > *)x;
    std::complex < T > *cy = (std::complex < T > *)y;
    std::complex < T > *cres = (std::complex < T > *)res;
    size_t c_size = f_size / 2;

    for (size_t i = 0; i < c_size; i++)
      cres[i] = ca * cx[i] + cy[i];
  }

  template < class T >
    void caxpy_threaded (T * res, std::complex < T > ca, T * x, T * y,
                         size_t f_size)
  {
    std::complex < T > *cx = (std::complex < T > *)x;
    std::complex < T > *cy = (std::complex < T > *)y;
    std::complex < T > *cres = (std::complex < T > *)res;
    size_t c_size = f_size / 2;

#pragma omp parallel for
    for (size_t i = 0; i < c_size; i++)
      cres[i] = ca * cx[i] + cy[i];
  }

  template < class T > void scale_single (T * res, T s, size_t f_size)
  {
#pragma omp parallel for
    for (size_t i = 0; i < f_size; i++)
      res[i] *= s;
  }

  template < class T >
    void caxpy (T * res, std::complex < T > ca, T * x, T * y, size_t f_size)
  {
    std::complex < T > *cx = (std::complex < T > *)x;
    std::complex < T > *cy = (std::complex < T > *)y;
    std::complex < T > *cres = (std::complex < T > *)res;
    size_t c_size = f_size / 2;

#pragma omp parallel for
    for (size_t i = 0; i < c_size; i++)
      cres[i] = ca * cx[i] + cy[i];
  }

  template < class T > std::complex < T > sp_single (T * a, T * b,
                                                     size_t f_size) {
    std::complex < T > *ca = (std::complex < T > *)a;
    std::complex < T > *cb = (std::complex < T > *)b;
    size_t c_size = f_size / 2;

    std::complex < T > ret = 0.0;
#pragma omp parallel for
    for (size_t i = 0; i < c_size; i++)
      ret += conj (ca[i]) * cb[i];

    return ret;
  }

  template < class T > std::complex < T > sp (T * a, T * b, size_t f_size) {
    std::complex < T > *ca = (std::complex < T > *)a;
    std::complex < T > *cb = (std::complex < T > *)b;
    size_t c_size = f_size / 2;

    std::complex < T > res = 0.0;
#pragma omp parallel shared(res)
    {
      std::complex < T > resl = 0.0;
#pragma omp for
      for (size_t i = 0; i < c_size; i++)
        resl += conj (ca[i]) * cb[i];

#pragma omp critical
      {
        res += resl;
      }
    }
    return res;
  }

  template < class T > T norm_of_evec (std::vector < std::vector < T > >&v,
                                       int j) {
    T gg = 0.0;
#pragma omp parallel shared(gg)
    {
      T ggl = 0.0;
#pragma omp for
      for (int nb = 0; nb < args.blocks; nb++) {
        T *res = &v[nb][(int64_t) f_size_block * j];
        ggl += sp_single (res, res, f_size_block).real ();
      }

#pragma omp critical
      {
        gg += ggl;
      }
    }
    return gg;
  }

    void write_bytes (void *buf, int64_t s, FILE * f, uint32_t & crc)
    {
      static double data_counter = 0.0;
      const char *fname="write_bytes()";
      if( s==0) return;
      if(s<0) 
        ERR.General(cname,fname, "size(%d) negative!\n", s);

//      printf("write_bytes(%p %ld %p %d\n",buf,s,f,&crc);

    // checksum
//      crc = crc32_fast (buf, s, crc);
    crc = crc32_loop (crc, (Bytef *) buf, s);

      double t0 = dclock ();
      int ret;
      if ((ret = fwrite (buf, s, 1, f)) != 1) {
        ERR.General(cname,fname, "Write failed! %d %d\n", s, ret);
      }
      double t1 = dclock ();

    data_counter += (double) s;
    if (data_counter > 1024. * 1024. * 256) {
      VRB.Result (cname, "write_bytes", "Writing at %g GB/s\n",
                  (double) s / 1024. / 1024. / 1024. / (t1 - t0));
      data_counter = 0.0;
    }
  }

  void write_floats (FILE * f, uint32_t & crc, OPT * in, int64_t n)
  {
    const char *fname = "write_floats()";
    float *buf = (float *) smalloc (cname, fname, "buf", sizeof (float) * n);
    if (!buf) {
      ERR.General (cname, fname, "Out of mem\n");
    }
    // convert to float if needed
#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
      buf[i] = in[i];

    fix_float_endian (buf, n);

    write_bytes (buf, n * sizeof (float), f, crc);

    free (buf);
  }

  void read_floats (char *&ptr, OPT * out, int64_t n)
  {
    float *in = (float *) ptr;
    ptr += 4 * n;

#pragma omp parallel for
    for (size_t i = 0; i < n; i++)
      out[i] = in[i];
    fix_float_endian (out, n);
//      std::cout << "out[0]= " << out[0] << std::endl;
    if (std::isnan (out[0]))
      std::cout << "read_floats out[0]= " << out[0] << std::endl;
  }

  int fp_map (float in, float min, float max, int N)
  {
    // Idea:
    //
    // min=-6
    // max=6
    //
    // N=1
    // [-6,0] -> 0, [0,6] -> 1;  reconstruct 0 -> -3, 1-> 3
    //
    // N=2
    // [-6,-2] -> 0, [-2,2] -> 1, [2,6] -> 2;  reconstruct 0 -> -4, 1->0, 2->4
    int ret = (int) ((float) (N + 1) * ((in - min) / (max - min)));
    if (ret == N + 1) {
      ret = N;
    }
    return ret;
  }

  float fp_unmap (unsigned short val, float min, float max, int N)
  {
    unsigned short tmp = val;
    fix_short_endian (&val, 1);
    if ((float) ((int) val + 0.5) != (float) (val + 0.5))
      std::cout << tmp << "after fix \t" << val << std::endl;
    return min + (float) ((int) val + 0.5) * (max - min) / (float) (N + 1);
  }

#define SHRT_UMAX 65535
#define BASE 1.4142135623730950488

  float unmap_fp16_exp (unsigned short e)
  {
    float de = (float) ((int) e - SHRT_UMAX / 2);
    return pow (BASE, de);
  }

  void read_floats_fp16 (char *&ptr, OPT * out, int64_t n, int nsc)
  {

    const char *fname = "read_floats_fp16()";

    int64_t nsites = n / nsc;
    if (n % nsc) {
      ERR.General (cname, fname, "Invalid size in read_floats_fp16\n");
    }

    unsigned short *in = (unsigned short *) ptr;
    ptr += 2 * (n + nsites);

//#define assert(exp)  { if ( !(exp) ) { fprintf(stderr,"Assert " #exp " failed\n"); exit(84); } }

    // do for each site
#pragma omp parallel for
    for (int64_t site = 0; site < nsites; site++) {

      OPT *ev = &out[site * nsc];

      unsigned short *bptr = &in[site * (nsc + 1)];

      unsigned short exp = *bptr++;
      fix_short_endian (&exp, 1);
      OPT max = unmap_fp16_exp (exp);
      OPT min = -max;

      for (int i = 0; i < nsc; i++) {
        ev[i] = fp_unmap (*bptr++, min, max, SHRT_UMAX);
      }
      if (std::isnan (ev[0]))
        std::cout << "read_floats_fp16 ev[0]= " << ev[0] << std::endl;

    }

  }

  int read_metadata (const char *root, std::vector < Float > &_evals)
  {
    const char *fname = "read_metadata()";
    std::string path (root);

    char buf[1024];
    sprintf (buf, "%s/metadata.txt", root);
    FILE *f = NULL;
    uint32_t status = 0;

    if (UniqueID () == 0) {
      VRB.Result (cname, fname, "node 0, before fopen %s\n", buf);
      f = fopen (buf, "r");
      status = f ? 1 : 0;
      VRB.Result (cname, fname, "node 0, fopen %s returns %d\n", buf, status);
    }
    sumArray (&status, 1);
    printf ("node %d: status %d\n", UniqueID (), status);
    if (!status) {
      ERR.General (cname, fname, "failed to open %s \n", buf);
      // should terminate here
      return false;
    }
    for (int i = 0; i < 5; i++) {
      args.s[i] = 0;
      args.b[i] = 0;
      args.nb[i] = 0;
    }
    args.neig = 0;
    args.nkeep = 0;
    args.nkeep_single = 0;
    args.blocks = 0;
    args.FP16_COEF_EXP_SHARE_FLOATS = 0;

#define _IRL_READ_INT(buf,p) if (f) { assert(fscanf(f,buf,p)==1); } else { *(p) = 0; }
    if (UniqueID () == 0) {
      VRB.Debug (cname, fname, "node 0, before reading metadata\n");
      for (int i = 0; i < 5; i++) {
        sprintf (buf, "s[%d] = %%d\n", i);
        _IRL_READ_INT (buf, &args.s[i]);
      }
      for (int i = 0; i < 5; i++) {
        sprintf (buf, "b[%d] = %%d\n", i);
        _IRL_READ_INT (buf, &args.b[i]);
      }
      for (int i = 0; i < 5; i++) {
        sprintf (buf, "nb[%d] = %%d\n", i);
        _IRL_READ_INT (buf, &args.nb[i]);
      }
      _IRL_READ_INT ("neig = %d\n", &args.neig);
      _IRL_READ_INT ("nkeep = %d\n", &args.nkeep);
      _IRL_READ_INT ("nkeep_single = %d\n", &args.nkeep_single);
      _IRL_READ_INT ("blocks = %d\n", &args.blocks);
//      _IRL_READ_INT ("big_endian = %d\n", &args.big_endian);
      _IRL_READ_INT ("FP16_COEF_EXP_SHARE_FLOATS = %d\n",
                     &args.FP16_COEF_EXP_SHARE_FLOATS);
      VRB.Result (cname, fname, "node 0, after reading metadata \n");
    }
    //     bigendian = args.big_endian; //currently fixed to be little endian
    sumArray (args.s, 5);
    sumArray (args.b, 5);
    sumArray (args.nb, 5);
    sumArray (&args.neig, 1);
    sumArray (&args.nkeep, 1);
    sumArray (&args.nkeep_single, 1);
    sumArray (&args.blocks, 1);
    sumArray (&args.FP16_COEF_EXP_SHARE_FLOATS, 1);
//    exit (-42);

//we do not divide the fifth dimension
    std::vector < int >nn (5, 1);
    nn[4] = 1;
    for (int i = 0; i < 4; i++) {
      nn[i] = GJP.Nodes (i);
      nprocessors *= nn[i];
      nfile *= (GJP.Sites (i) / args.s[i]);
    }
    double barrier = 0;
    sumArray (&barrier, 1);

    if (UniqueID () == 0) {
      VRB.Debug (cname, fname, "node-layout %d %d %d %d %d nprocessors %d\n",
                 nn[0], nn[1], nn[2], nn[3], nn[4], nprocessors);
      std::cout << "nprocessor= " << nprocessors << " nfile= " << nfile <<
        std::endl;
    }
    args.crc32_header.resize (nfile);
    VRB.Debug (cname, fname, "node 0, before reading crc32\n");
    if (UniqueID () == 0) {
      for (uint32_t i = 0; i < nfile; i++) {
        sprintf (buf, "crc32[%d] = %%X\n", i);
        _IRL_READ_INT (buf, &args.crc32_header[i]);
        VRB.Debug (cname, fname, "crc32[%d] = %X\n", i, args.crc32_header[i]);
      }
      VRB.Debug (cname, fname, "node 0, after reading crc32\n");
    }
    //printf("node %d, before sumarray crc32\n", UniqueID());
    sumArray (&args.crc32_header[0], nfile);

#undef _IRL_READ_INT
    {

//first debug for this
#if 1
//      int ngroup = (nprocessors - 1) / nfile + 1;       //number of nodes reading the same file, maybe overkill? NOT WORKING. TO BE DEBUGGED
      int ngroup = 1;

      int nodeID = UniqueID ();
      int slot = nodeID / ngroup;
      int first = 0;
      int nperdir = nfile / n_cycle;
      if (nfile % n_cycle)
        nperdir += 1;
      char *raw_in = NULL;
      off_t raw_size = 0;
      std::vector < uint32_t > crc32_part (nprocessors);
      while (slot < nfile) {
        int dir = slot / nperdir;

        sprintf (buf, "%s/%2.2d/%10.10d.compressed", root, dir, slot);
        FILE *f2 = fopen (buf, "rb");
        if (!f2) {
          fprintf (stderr, "Could not open %s\n", buf);
          //return 3;
          sleep (2);
          f2 = fopen (buf, "rb");
          if (!f2) {
            ERR.General (cname, fname, "Could not open %s again.\n", buf);
          }
        }
        fseeko (f2, 0, SEEK_END);
        off_t size0 = ftello (f2);
//          std::cout << "size0= " << size0 << std::endl;
        off_t read_size = size0 / ngroup;
//          std::cout << "read_size= " << read_size << std::endl;
        off_t offset = read_size * (nodeID % ngroup);
        if (((nodeID % ngroup) == (ngroup - 1))
            || (nodeID == (nprocessors - 1))) {
          read_size = size0 - offset;
        }
//          std::cout << "offset= " << offset << std::endl;

        if (raw_size < read_size) {
          VRB.Result (cname, fname, "offset=%d read_size=%d\n", offset,
                      read_size);
          if (raw_in)
            sfree (raw_in);
          raw_in = (char *) smalloc (cname, fname, "raw_in", read_size);
          raw_size = read_size;
        }
        if (0) {
          off_t half = read_size / 2;
          fseeko (f2, 0, SEEK_SET);
          fread (raw_in, 1, half, f2);
          uint32_t first = crc32_loop (0, (Bytef *) raw_in, half);
          fseeko (f2, half, SEEK_SET);
          fread (raw_in, 1, half, f2);
          uint32_t second = crc32_loop (0, (Bytef *) raw_in, half);
//          printf ("half first second = %d %x %x\n", half, first, second);
//            std::cout << "half " << half << " first " << first << " second " << second << std::endl;
        }
        fseeko (f2, offset, SEEK_SET);
        if (!raw_in) {
          fprintf (stderr, "Out of mem\n");
          return 5;
        }

        off_t t_pos = ftello (f2);
        if (fread (raw_in, 1, read_size, f2) != read_size) {
          fprintf (stderr, "Invalid fread\n");
          return 6;
        }
//        off_t t_pos2 = ftello (f2);
        VRB.Debug (cname, fname, "%d read: %d %d\n", nodeID, t_pos, t_pos);

//        std::vector < uint32_t > crc32_part (nprocessors);
//        uint32_t crc32_part[nprocessors];
        for (uint32_t i = 0; i < nprocessors; i++)
          crc32_part[i] = 0;
        crc32_part[nodeID] = crc32_loop (0, (Bytef *) raw_in, read_size);
        VRB.Debug (cname, fname,
                   "%d %d: ngroup read_size offset crc32 : %d %d %d %x\n",
                   nprocessors, nodeID, ngroup, read_size, offset,
                   crc32_part[nodeID]);
//        printf ("sumArray: node %d first %d slot %d\n", nodeID, first, slot);
        if (ngroup > 1)
          sumArray (crc32_part.data (), nprocessors);
        if (nodeID % ngroup == 0) {
          uint32_t crc32_all = crc32_part[nodeID];
          VRB.Debug (cname, fname, "%d: crc32_all: %x\n", nodeID, crc32_all);
          for (int i = 1; i < ngroup; i++) {
            VRB.Debug (cname, fname, "%d: crc32_part[%d]: %x\n", nodeID,
                       nodeID + i, crc32_part[nodeID + i], crc32_all);
            if (i < (ngroup - 1))
              crc32_all =
                crc32_combine (crc32_all, crc32_part[nodeID + i], read_size);
            else
              crc32_all =
                crc32_combine (crc32_all, crc32_part[nodeID + i],
                               size0 - (read_size * (ngroup - 1)));
          }
          if (VRB.IsActivated (VERBOSE_DEBUG_LEVEL)) 
          printf ("%d: crc32_all: %x crc32_header %x\n", slot, crc32_all,
                  args.crc32_header[slot]);
//          assert (crc32_all == args.crc32_header[slot]);
          if (crc32_all != args.crc32_header[slot])
	  ERR.General(cname,fname,"slot %d: crc32(%x) does not agree with the header(%x)\n",
		 slot, crc32_all, args.crc32_header[slot]);

        }
        if (f2)
          fclose (f2);
        slot += nprocessors / ngroup;   // in case nfile > nprocessors
        first += nprocessors / ngroup;
      }
// in case there were idle nodes at the end
      if ((first < nfile) && (slot >= nfile)) {
//        printf ("idle one needs sumArray: node %d first %d slot %d\n", nodeID, first, slot);
        for (uint32_t i = 0; i < nprocessors; i++)
          crc32_part[i] = 0;
        if (ngroup > 1)
          sumArray (crc32_part.data (), nprocessors);
      }
      if (raw_in)
        free (raw_in);
      raw_in = NULL;
      raw_size = 0;
#endif
      crc32_checked = true;

    }

    if (f)
      fclose (f);

//    double vals[nvec];
    long nvec = 0;
    evals.resize (args.neig);
    _evals.resize (args.neig);
    double *vals = evals.data ();
//    memset (vals, 0, sizeof (vals));
    if (!UniqueID ()) {
      std::cout << "Reading eigenvalues \n";
      std::string filename = path + "/eigen-values.txt.smoothed";
      FILE *file = fopen (filename.c_str (), "r");
      if (!file) {
        filename.clear ();
        filename = path + "/eigen-values.txt";
        VRB.Result (cname, fname,
                    "smoothed eignevalues not available. Trying %s\n",
                    filename.c_str ());
        file = fopen (filename.c_str (), "r");
      }
      fscanf (file, "%ld\n", &nvec);
      assert (nvec <= args.neig);
      for (int i = 0; i < args.neig; i++) {
        fscanf (file, "%lE\n", vals + i);
//        std::cout << (evals[i]) << std::endl;
      }
      fclose (file);
    }
    sumArray (vals, args.neig);
    if (!UniqueID ())
      std::cout << "End Reading eigenvalues, read_metadata done\n";
    for (int i = 0; i < args.neig; i++) {
      _evals[i] = evals[i];
      VRB.Result (cname, fname, "eval[%d]=%0.14e\n", i, _evals[i]);
    }

    return 1;
  }

  int globalToLocalCanonicalBlock (int slot,
                                   const std::vector < int >&src_nodes, int nb);
  void get_read_geometry (const std::vector < int >&cnodes, std::map < int,
                          std::vector < int > >&slots,
                          std::vector < int >&slot_lvol,
                          std::vector < int >&lvol, int64_t & slot_lsites,
                          int &ntotal);

//    int decompress (const char *root_, std::vector < OPT * >&dest_all);
  int read_compressed_vectors (const char *root, const char *checksum_dir,
                               std::vector < OPT * >&dest_all,
                               int start = -1, int end = -1,
                               float time_out = 100000);

  int read_compressed_blocks (const char *root, const char *checksum_dir,
                              float time_out = 100000);
  int build_evecs (std::vector < OPT * >&dest_all, int start = -1, int end =
                   -1);
  int dot (std::vector < OPT * >&_rhs, int start, int end,
           std::vector < Complex > &coef);
  template < class FLOAT >
    int proj (std::vector < FLOAT * >&_sol, std::vector < FLOAT * >&_rhs,
              int start, int end, std::vector < Complex > &_coef)
  {

    const char *fname = "proj()";
    static Timer timer (cname, fname);
    timer.start ();

    if (start < 0)
      start = 0;

    int nrhs = _rhs.size ();

    VRB.Result (cname, fname, "rhs.size() %d  args.neig %d \n",
                nrhs, args.neig);

    if (end < 0)
      end = args.neig;
    if (end > args.neig)
      end = args.neig;
    VRB.Result (cname, fname, "start %d end %d \n", start, end);

    int nvec = end - start;
    _coef.resize (nrhs * nvec);

    size_t f_size_block = f_size / nb_per_node;
    VRB.Result (cname, fname, "nb_per_node %d f_size_block %d \n",
                nb_per_node, f_size_block);

//      std::vector < std::vector < EvecReader::OPT > >block_data;
    block_data.resize (nb_per_node);
    for (int i = 0; i < nb_per_node; i++) {
      block_data[i].resize (f_size_block);
      memset (&block_data[i][0], 0, sizeof (float) * f_size_block);
    }

//  for (int i = start; i < end; i++)
//    memset (dest_all[i - start], 0, f_size * sizeof (OPT));

    double t0 = dclock ();
    int s_l[5], nb_l[5];

    s_l[0] = GJP.NodeSites (0);
    s_l[1] = GJP.NodeSites (1);
    s_l[2] = GJP.NodeSites (2);
    s_l[3] = GJP.NodeSites (3);
    s_l[4] = GJP.NodeSites (4);

    nb_l[0] = s_l[0] / args.b[0];
    nb_l[1] = s_l[1] / args.b[1];
    nb_l[2] = s_l[2] / args.b[2];
    nb_l[3] = s_l[3] / args.b[3];
    nb_l[4] = s_l[4] / args.b[4];

//#pragma omp parallel
    for (int j = start; j < end; j++) {

      double ta, tb;
      int tid = omp_get_thread_num ();

//      if (!tid)
//        ta = dclock ();
      time_elapse ();


#pragma omp parallel for
      for (int nb = 0; nb < nb_per_node; nb++) {

        float *dest_block = &block_data[nb][0];

        {
          // do reconstruction of this block
          memset (dest_block, 0, sizeof (float) * f_size_block);
          for (int i = 0; i < args.nkeep; i++) {
            OPT *ev_i = &block_data_ortho[nb][(int64_t) f_size_block * i];
            OPT *coef = &block_coef[nb][2 * (i + args.nkeep * j)];
            caxpy_single (dest_block, *(std::complex < float >*) coef, ev_i,
                          dest_block, f_size_block);
          }
        }
      }

//      if (!tid) {
//        tb = dclock ();
//      if (j % 100 == 0)
      VRB.Result (cname, fname, "1 evec reconstruction: %d - %g seconds\n", j,
                  time_elapse ());
//      }
      //int[5] loc_s = {GJP.NodeSites(0), GJP.NodeSites(1), GJP.NodeSites(2), GJP.NodeSites(3),GJP.NodeSites(4)};

      Float factor = 1. / evals[j];
      int n_print = 0;
      for (int rhs = 0; rhs < nrhs; rhs++) {
//      OPT *dest = dest_all[j - start];
        FLOAT *rhs_p = _rhs[rhs];
        Float re = 0., im = 0.;
#pragma omp parallel for reduction(+:re,im)
        for (int idx = 0; idx < vol4d; idx++) {
          int pos[5], pos_in_block[5], block_coor[5];
          index_to_pos (idx, pos, s_l);

          int parity = (pos[0] + pos[1] + pos[2] + pos[3]) % 2;
          if (parity == 1) {

            for (pos[4] = 0; pos[4] < s_l[4]; pos[4]++) {
              pos_to_blocked_pos (pos, pos_in_block, block_coor);

              int bid = pos_to_index (block_coor, nb_l);
              int ii = pos_to_index (pos_in_block, args.b) / 2;
              float *vec = &block_data[bid][ii * 24];
              size_t cps_index = get_cps_index (pos, 0, s_l);

              int co;
              for (co = 0; co < 12; co++) {
//              OPT *src = &dest[get_cps_index (pos, co, s_l)];
//                FLOAT *src = rhs_p + get_cps_index (pos, co, s_l);
                FLOAT *src = rhs_p + cps_index + 2 * co;
                re += src[0] * vec[2 * co + 0] + src[1] * vec[2 * co + 1];
                im += -src[0] * vec[2 * co + 1] + src[1] * vec[2 * co + 0];
                if (0)
                  if (bid > 1532 && (!co) && n_print < 100 && (!rhs)
                      && (!pos[4])) {
                    VRB.Result (cname, fname,
                                "bid=%d ii=%d pos=%d %d %d %d %d co=%d  vec=%e %e cps_index=%p %p %d src=%e %e\n",
                                bid, ii, pos[0], pos[1], pos[2], pos[3], pos[4],
                                co, vec[2 * co + 0], vec[2 * co + 1], src,
                                rhs_p, (src - rhs_p), src[0], src[1]);
                    n_print++;
                  }
//              out[0] = dst[2 * co + 0];
//              out[1] = dst[2 * co + 1];
              }
            }
          }
        }
        glb_sum (&re), glb_sum (&im);
        _coef[rhs + nrhs * (j - start)] += Complex (re * factor, im * factor);
      }
      //     if (j % 100 == 0)
      VRB.Result (cname, fname, "Inner product : %d - %g seconds\n", j,
                  time_elapse ());

#pragma omp parallel for
      for (int idx = 0; idx < vol4d; idx++) {
        int pos[5], pos_in_block[5], block_coor[5];
        index_to_pos (idx, pos, s_l);

        int parity = (pos[0] + pos[1] + pos[2] + pos[3]) % 2;
        if (parity == 1) {

          for (pos[4] = 0; pos[4] < s_l[4]; pos[4]++) {
            pos_to_blocked_pos (pos, pos_in_block, block_coor);

            int bid = pos_to_index (block_coor, nb_l);
            int ii = pos_to_index (pos_in_block, args.b) / 2;
            float *vec = &block_data[bid][ii * 24];
            size_t cps_index = get_cps_index (pos, 0, s_l);
            for (int rhs = 0; rhs < nrhs; rhs++) {
              FLOAT *sol_p = _sol[rhs];
              int co;
              for (co = 0; co < 12; co++) {
//              OPT *src = &dest[get_cps_index (pos, co, s_l)];
                FLOAT *sol = sol_p + cps_index + 2 * co;
                double re = _coef[rhs + nrhs * (j - start)].real ();
                double im = _coef[rhs + nrhs * (j - start)].imag ();
                sol[0] += re * vec[2 * co + 0] - im * vec[2 * co + 1];
                sol[1] += im * vec[2 * co + 0] + re * vec[2 * co + 1];
//              out[0] = dst[2 * co + 0];
//              out[1] = dst[2 * co + 1];
              }
            }
          }
        }
      }
//      if (j % 100 == 0)
      VRB.Result (cname, fname, "Projection : %d - %g seconds\n", j,
                  time_elapse ());
    }

    double t1 = dclock ();

    cps::sync ();
    VRB.Flow (cname, fname,
              "Reconstruct and inner product  of eigenvectors in %g seconds\n",
              t1 - t0);

    cps::sync ();
    timer.stop (true);
#undef FP_16_SIZE
    return true;
  }

  template < class FLOAT >
    int proj2 (std::vector < FLOAT * >&_sol, std::vector < FLOAT * >&_rhs,
               int start, int end, std::vector < Complex > &_coef)
  {

    const char *fname = "proj2()";
    static Timer timer (cname, fname);
    timer.start ();
    time_elapse ();

    if (start < 0)
      start = 0;

    int nrhs = _rhs.size ();

    VRB.Result (cname, fname, "rhs.size() %d  args.neig %d \n",
                nrhs, args.neig);

    if (end < 0)
      end = args.neig;
    if (end > args.neig)
      end = args.neig;
    VRB.Result (cname, fname, "start %d end %d \n", start, end);

    int nvec = end - start;
    _coef.resize (nrhs * nvec);

    size_t f_size_block = f_size / nb_per_node;
    VRB.Result (cname, fname, "nb_per_node %d f_size_block %d \n",
                nb_per_node, f_size_block);

    double t0 = dclock ();
    int s_l[5], nb_l[5];

    s_l[0] = GJP.NodeSites (0);
    s_l[1] = GJP.NodeSites (1);
    s_l[2] = GJP.NodeSites (2);
    s_l[3] = GJP.NodeSites (3);
    s_l[4] = GJP.NodeSites (4);

    nb_l[0] = s_l[0] / args.b[0];
    nb_l[1] = s_l[1] / args.b[1];
    nb_l[2] = s_l[2] / args.b[2];
    nb_l[3] = s_l[3] / args.b[3];
    nb_l[4] = s_l[4] / args.b[4];

    std::vector < std::vector < std::vector < Float > > >coef3;
    coef3.resize (nrhs);
#ifdef PROJ2_DEBUG
    std::vector < std::vector < std::vector < Float > > >coef4;
    std::vector < std::vector < std::vector < Float > > >coef5;
    std::vector < std::vector < std::vector < Float > > >coef6;
    coef4.resize (nrhs);
    coef5.resize (nrhs);
    coef6.resize (nrhs);
#endif
    for (int rhs = 0; rhs < nrhs; rhs++) {
      coef3[rhs].resize (nb_per_node);
#ifdef PROJ2_DEBUG
      coef4[rhs].resize (nb_per_node);
      coef5[rhs].resize (nb_per_node);
      coef6[rhs].resize (nb_per_node);
#endif
      for (int nb = 0; nb < nb_per_node; nb++) {
        coef3[rhs][nb].resize (args.nkeep * 2);
#ifdef PROJ2_DEBUG
        coef4[rhs][nb].resize (args.nkeep * 2);
        coef5[rhs][nb].resize (args.nkeep);
        coef6[rhs][nb].resize (args.nkeep);
#endif
        for (int i = 0; i < args.nkeep; i++) {
          coef3[rhs][nb][i * 2] = 0.;
          coef3[rhs][nb][i * 2 + 1] = 0.;
#ifdef PROJ2_DEBUG
          coef4[rhs][nb][i * 2] = 0.;
          coef4[rhs][nb][i * 2 + 1] = 0.;
          coef5[rhs][nb][i] = 0.;
          coef6[rhs][nb][i] = 0.;
#endif
    }}}

    std::vector < Float > all_src (nrhs, 0.);
    int n_print = 0;
    VRB.Result (cname, fname, "Setup time=%e sec\n", time_elapse ());
    for (size_t idx = 0; idx < vol4d; idx++) {
//      if (idx % 1000 == 0) VRB.Result (cname, fname, "idx=%d time=%e sec\n", idx, time_elapse ());
      int pos[5], pos_in_block[5], block_coor[5];
      index_to_pos (idx, pos, s_l);

      int parity = (pos[0] + pos[1] + pos[2] + pos[3]) % 2;
      if (parity == 1) {

        for (pos[4] = 0; pos[4] < s_l[4]; pos[4]++) {
          pos_to_blocked_pos (pos, pos_in_block, block_coor);

          int bid = pos_to_index (block_coor, nb_l);
          int ii = pos_to_index (pos_in_block, args.b) / 2;
          size_t cps_index = get_cps_index (pos, 0, s_l);
#pragma omp parallel for
          for (int i = 0; i < args.nkeep; i++) {
            float *ev_i =
              &block_data_ortho[bid][ii * 24 + (int64_t) f_size_block * i];

            for (int rhs = 0; rhs < nrhs; rhs++) {
              FLOAT *rhs_p = _rhs[rhs];
              for (int co = 0; co < 12; co++) {
                FLOAT *src = rhs_p + get_cps_index (pos, co, s_l);
                if (i == 0)
                  all_src[rhs] += src[0] * src[0] + src[1] * src[1];
                if (0)
                  if (bid > 1532 && (!co) && n_print < 100 && (!i) && (!rhs)
                      && (!pos[4])) {
                    VRB.Result (cname, fname,
                                "bid=%d ii=%d pos=%d %d %d %d %d co=%d i=%d ev_i=%e %e cps_index=%p %p %d src=%e %e\n",
                                bid, ii, pos[0], pos[1], pos[2], pos[3], pos[4],
                                co, i, ev_i[2 * co + 0], ev_i[2 * co + 1], src,
                                rhs_p, (src - rhs_p), src[0], src[1]);
                    n_print++;
                  }
                Float re =
                  src[0] * ev_i[2 * co + 0] + src[1] * ev_i[2 * co + 1];
                Float im =
                  -src[0] * ev_i[2 * co + 1] + src[1] * ev_i[2 * co + 0];
                coef3[rhs][bid][i * 2] += re;
                coef3[rhs][bid][i * 2 + 1] += im;
#ifdef PROJ2_DEBUG
                coef4[rhs][bid][i * 2] += 1.;
                coef4[rhs][bid][i * 2 + 1] += 0.;
                coef5[rhs][bid][i] +=
                  ev_i[2 * co + 0] * ev_i[2 * co + 0] + ev_i[2 * co +
                                                             1] * ev_i[2 * co +
                                                                       1];
                coef6[rhs][bid][i] += src[0] * src[0] + src[1] * src[1];
#endif
              }
            }
          }
        }
      }
    }
    VRB.Result (cname, fname, "Initial projection time=%e sec\n",
                time_elapse ());
#ifdef PROJ2_DEBUG
    glb_sum (all_src.data (), nrhs);
    for (int rhs = 0; rhs < nrhs; rhs++)
      VRB.Result (cname, fname, "all_src[%d]=%e\n", rhs, all_src[rhs]);
    for (int rhs = 0; rhs < nrhs; rhs++) {
      Float all_rhs = 0.;
      Float all_block = 0.;
      for (int nb = 0; nb < nb_per_node; nb++) {
        for (int i = 0; i < args.nkeep; i++) {
          all_block += coef5[rhs][nb][i];
          all_rhs += coef6[rhs][nb][i];
#if 0
          VRB.Result (cname, fname,
                      "coef3[%d][%d][%d]=%0.4e %0.4e coef4[%d][%d][%d]=%0.4e %0.4e\n",
                      rhs, nb, i, coef3[rhs][nb][2 * i],
                      coef3[rhs][nb][2 * i + 1], rhs, nb, i,
                      coef4[rhs][nb][2 * i], coef4[rhs][nb][2 * i + 1]);
#endif
      }}
      glb_sum (&all_rhs);
      glb_sum (&all_block);
      VRB.Result (cname, fname, "block[%d]=%e rhs[%d]=%e \n", rhs, all_block,
                  rhs, all_rhs);
    }
#endif
//`exit(-54);

#ifdef PROJ2_DEBUG
    Float coef7[end - start][nb_per_node];
#pragma omp parallel for
    for (int j = start; j < end; j++) {
      for (int bid = 0; bid < nb_per_node; bid++) {
        coef7[j - start][bid] = 0.;
    }}
    VRB.Result (cname, fname, "Debug info time=%e sec\n", time_elapse ());
#endif

    Float coef8[nrhs][2 * (end - start)];
#pragma omp parallel for
    for (int j = start; j < end; j++) {
      for (int rhs = 0; rhs < nrhs; rhs++) {
        coef8[rhs][2 * (j - start)] = 0.;
        coef8[rhs][2 * (j - start) + 1] = 0.;
    }}

#pragma omp parallel for
    for (int rhs = 0; rhs < nrhs; rhs++) {
      for (int bid = 0; bid < nb_per_node; bid++) {
        for (int i = 0; i < args.nkeep; i++) {
          Float re = coef3[rhs][bid][i * 2];
          Float im = coef3[rhs][bid][i * 2 + 1];
          Float re2 = 0., im2 = 0.;
          for (int j = start; j < end; j++) {
            Float factor = 1. / evals[j];
            OPT *coef = &block_coef[bid][2 * (i + args.nkeep * j)];
            coef8[rhs][2 * (j - start)] +=
              factor * (re * coef[0] + im * coef[1]);
            coef8[rhs][2 * (j - start) + 1] +=
              factor * (im * coef[0] - re * coef[1]);
#ifdef PROJ2_DEBUG
            if (rhs == 0)
              coef7[j - start][bid] += coef[0] * coef[0] + coef[1] * coef[1];
#endif
          }
        }
      }
    }
    glb_sum (&coef8[0][0], 2 * (end - start) * nrhs);
#ifdef PROJ2_DEBUG
    glb_sum (&coef7[0][0], (end - start) * nb_per_node);
#endif
    for (int j = start; j < end; j++) {
      for (int rhs = 0; rhs < nrhs; rhs++) {
        VRB.Flow (cname, fname, "coef[%d][%d]=%e %e\n", j - start, rhs,
                  coef8[rhs][2 * (j - start)], coef8[rhs][2 * (j - start) + 1]);
        _coef[rhs + nrhs * (j - start)] =
          Complex (coef8[rhs][2 * (j - start)],
                   coef8[rhs][2 * (j - start) + 1]);
    }}
#ifdef PROJ2_DEBUG
    for (int j = start; j < end; j++) {
      Float all = 0.;
      for (int bid = 0; bid < nb_per_node; bid++) {
        all += coef7[j - start][bid];
      }
      VRB.Result (cname, fname, "all[%d]= %e\n", j, all);
    }
    VRB.Result (cname, fname, "Debug time=%e sec\n", time_elapse ());
#endif

    Float coef9[nrhs][nb_per_node][2 * args.nkeep];
#pragma omp parallel for
    for (int bid = 0; bid < nb_per_node; bid++) {
      for (int i = 0; i < args.nkeep; i++) {
        for (int rhs = 0; rhs < nrhs; rhs++) {
          Float re = 0., im = 0.;
          for (int j = start; j < end; j++) {
            Float *tmp1 = &coef8[rhs][2 * (j - start)];
            OPT *tmp2 = &block_coef[bid][2 * (i + args.nkeep * j)];
            re += tmp1[0] * tmp2[0] - tmp1[1] * tmp2[1];
            im += tmp1[0] * tmp2[1] + tmp1[1] * tmp2[0];
          }
          coef9[rhs][bid][2 * i] = re;
          coef9[rhs][bid][2 * i + 1] = im;
        }
      }
    }
    VRB.Result (cname, fname, "Coef calc time=%e sec\n", time_elapse ());

#pragma omp parallel for
    for (int idx = 0; idx < vol4d; idx++) {
      int pos[5], pos_in_block[5], block_coor[5];
      index_to_pos (idx, pos, s_l);
      int parity = (pos[0] + pos[1] + pos[2] + pos[3]) % 2;
      if (parity == 1) {
        for (pos[4] = 0; pos[4] < s_l[4]; pos[4]++) {
          pos_to_blocked_pos (pos, pos_in_block, block_coor);
          int bid = pos_to_index (block_coor, nb_l);
          int ii = pos_to_index (pos_in_block, args.b) / 2;
          int cps_index = get_cps_index (pos, 0, s_l);
          for (int i = 0; i < args.nkeep; i++) {
            OPT *vec =
              &block_data_ortho[bid][ii * 24 + (int64_t) f_size_block * i];
            for (int rhs = 0; rhs < nrhs; rhs++) {
              FLOAT *sol_p = _sol[rhs];
              Float re = 0., im = 0.;
#if 0
              for (int j = start; j < end; j++) {
                Float *tmp1 = &coef8[rhs][2 * (j - start)];
                OPT *tmp2 = &block_coef[bid][2 * (i + args.nkeep * j)];
                re += tmp1[0] * tmp2[0] - tmp1[1] * tmp2[1];
                im += tmp1[0] * tmp2[1] + tmp1[1] * tmp2[0];
              }
#else
              re = coef9[rhs][bid][2 * i];
              im = coef9[rhs][bid][2 * i + 1];
#endif
              for (int co = 0; co < 12; co++) {
//                  FLOAT *sol = sol_p + get_cps_index (pos, co, s_l);
                FLOAT *sol = sol_p + cps_index + 2 * co;
                sol[0] += re * vec[2 * co + 0] - im * vec[2 * co + 1];
                sol[1] += im * vec[2 * co + 0] + re * vec[2 * co + 1];
              }
            }
          }
        }
      }
    }
    VRB.Result (cname, fname, "Construction time=%e sec\n", time_elapse ());
    timer.stop (true);
#undef FP_16_SIZE
    return true;
  }


};

class EvecWriter:public EvecReader
{


  const char *cname;
public:
    EvecWriter ():EvecReader (), cname ("EvecWriter")
  {
//      machine_is_little_endian = machine_endian ();
  };
  ~EvecWriter () {
  }
  int writeCompressedVector (const
                             char
                             *dir,
                             OPT * V,
                             struct
                             evec_write &warg, std::vector < Float > &evals);
  template < class T > std::complex <
    T > sp_single (T * a, T * b, size_t f_size) {
    std::complex < T > *ca = (std::complex < T > *)a;
    std::complex < T > *cb = (std::complex < T > *)b;
    size_t c_size = f_size / 2;
    std::complex < T > ret = 0.0;
    for (size_t i = 0; i < c_size; i++)
      ret += conj (ca[i]) * cb[i];
    return ret;
  }

#if 1
  void get_coef (int nb, int i, int j)
  {
    if (f_size_block > 128) VRB.Flow ("", "get_coef()", "nb i j f_size_block=%d %d %d %d\n", nb, i, j, f_size_block);
    OPT *res = &block_data[nb][(size_t)
                               f_size_block * j];
    OPT *ev_i = &block_data_ortho[nb][(size_t)
                                      f_size_block * i];
    std::complex < OPT > c = sp_single (ev_i, res, f_size_block);
    OPT *cptr = &block_coef[nb][2 * (i + args.nkeep * j)];
    cptr[0] = c.real ();
    cptr[1] = c.imag ();
    if (f_size_block > 128)
      caxpy_threaded (res, -c, ev_i, res, f_size_block);
    else
      caxpy_single (res, -c, ev_i, res, f_size_block);
  }
#endif


// can assume that v >=0 and need to guarantee that unmap_fp16_exp(map_fp16_exp(v)) >= v
  unsigned short map_fp16_exp (float v)
  {
    const char *fname = "map_fp16_exp()";
    // float has exponents 10^{-44.85} .. 10^{38.53}
    int exp = (int) ceil (log (v) / log (BASE)) + SHRT_UMAX / 2;
    if (exp < 0)
      return 0;
    if (exp < 0 || exp > SHRT_UMAX) {
      ERR.General (cname, fname, "Error in map_fp16_exp(%0.14e,%d)\n", v, exp);
    }

    return (unsigned short) exp;
  }

  void write_floats_fp16 (FILE * f,
                          uint32_t & crc, OPT * in, int64_t n, int nsc)
  {

    const char *fname = "write_fp16_exp()";
    int64_t nsites = n / nsc;
    if (n % nsc) {
      ERR.General (cname, fname,
                   "Invalid size in write_floats_fp16 %ld %d\n", n, nsc);
    }

    unsigned short *buf = (unsigned short *) smalloc (cname,
                                                      fname,
                                                      "buf",
                                                      sizeof (short)
                                                      * (n + nsites));
    if (!buf) {
      ERR.General (cname, fname, "Out of mem\n");
    }
    // do for each site
#pragma omp parallel for
    for (size_t site = 0; site < nsites; site++) {

      OPT *ev = &in[site * nsc];
      unsigned short *bptr = &buf[site * (nsc + 1)];
      OPT max = fabs (ev[0]);
      OPT min;
      for (int i = 0; i < nsc; i++) {
        if (fabs (ev[i]) > max)
          max = fabs (ev[i]);
      }

      unsigned short exp = map_fp16_exp (max);
      max = unmap_fp16_exp (exp);
      min = -max;
      *bptr++ = exp;
      for (int i = 0; i < nsc; i++) {
        int val = fp_map (ev[i], min, max,
                          SHRT_UMAX);
        if (val < 0)
          val = 0;
        if (val < 0 || val > SHRT_UMAX) {
          ERR.General (cname, fname,
                       "Node %d: Assert failed: val = %d (%d), ev[i] = %.15g, min  = %0.15g max = %.15g, exp = %d\n",
                       UniqueID (), val, SHRT_UMAX, ev[i], min, max, (int) exp);
        }
        *bptr++ = (unsigned short) val;
      }

    }

    write_bytes (buf, sizeof (short) * (n + nsites), f, crc);
    free (buf);
  }

};
#endif
class Lexicographic
{
public:

  static inline void
    CoorFromIndex (std::vector <
                   int >&coor, int index, std::vector < int >&dims)
  {
    int nd = dims.size ();
      coor.resize (nd);
    for (int d = 0; d < nd; d++)
    {
      coor[d] = index % dims[d];
      index = index / dims[d];
    }
  }

  static inline void
    IndexFromCoor (std::vector <
                   int >&coor, int &index, std::vector < int >&dims)
  {
    int nd = dims.size ();
    int stride = 1;
    index = 0;
    for (int d = 0; d < nd; d++) {
      index = index + stride * coor[d];
      stride = stride * dims[d];
    }
  }
};

void alcf_evecs_save (char *dest, EigenCache * ec, int nkeep);
void movefloattoFloat (Float * out, float *in, size_t f_size);
CPS_END_NAMESPACE
#endif
