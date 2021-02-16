#define _FILE_OFFSET_BITS 64
#include <omp.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include <complex>
#include <vector>
#include <memory.h>
#include <iostream>

#include <sys/stat.h>
#include <comms/sysfunc_cps.h>
//#include "sumarray.h"
#include <util/time_cps.h>
#include <util/eigen_container.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/eig_io.h>
#include <unistd.h>

using namespace std;

USING_NAMESPACE_CPS CPS_START_NAMESPACE

size_t fread_check (void *ptr, off_t offset, size_t read_size, FILE * f)
{
  size_t total = 0, i_try = 0;
  while (total < read_size) {
    fseeko (f, offset, SEEK_SET);
    total = fread (ptr, 1, read_size, f);
    if (total < read_size)
      usleep (100);
    i_try++;
    if ((i_try % 100) == 0) {
      printf ("Node %d: fread failed %d times read_size=%d %d\n", UniqueID (),
              i_try, read_size, total);
//            total=1;//make it pass for now.
    }
    if ((i_try % 1000) == 0) {
      printf
        ("Node %d: fread failed %d times read_size=%d %d, letting it pass for now\n",
         UniqueID (), i_try, read_size, total);
      total = read_size;        //make it pass for now.
    }
  }
}

int EvecReader::globalToLocalCanonicalBlock (int slot,
                                             const std::vector < int >&src_nodes, int nb)
{

  // processor coordinate
  int _nd = (int) src_nodes.size ();
  std::vector < int >_src_nodes = src_nodes;
  std::vector < int >pco (_nd);
  Lexicographic::CoorFromIndex (pco, slot, _src_nodes);
  std::vector < int >cpco (pco);

  // get local block

  std::vector < int >_nb;
  _nb.resize (5);
  _nb[4] = GJP.SnodeSites () / args.b[4];       // b is the block size
  _nb[0] = GJP.XnodeSites () / args.b[0];
  _nb[1] = GJP.YnodeSites () / args.b[1];
  _nb[2] = GJP.ZnodeSites () / args.b[2];
  _nb[3] = GJP.TnodeSites () / args.b[3];
  std::vector < int >_nbc (_nb);
  assert (_nd == 5);
  std::vector < int >c_src_local_blocks (_nd);
  for (int i = 0; i < _nd; i++) {
    assert (GJP.Sites (i) % (src_nodes[i] * args.b[i]) == 0);
    c_src_local_blocks[i] = GJP.Sites (i) / src_nodes[i] / args.b[i];
  }
  std::vector < int >cbcoor (_nd);      // coordinate of block in slot in canonical form
  Lexicographic::CoorFromIndex (cbcoor, nb, c_src_local_blocks);

  // cpco, cbcoor
  std::vector < int >clbcoor (_nd);
  for (int i = 0; i < _nd; i++) {
    int cgcoor = cpco[i] * c_src_local_blocks[i] + cbcoor[i];   // global block coordinate
    int pcoor = cgcoor / _nbc[i];       // processor coordinate in my Grid
    int tpcoor = GJP.NodeCoor (i);
    if (pcoor != tpcoor)
      return -1;
    clbcoor[i] = cgcoor - tpcoor * _nbc[i];     // canonical local block coordinate for canonical dimension i
  }

  int lnb;
  Lexicographic::IndexFromCoor (clbcoor, lnb, _nbc);
  //std::cout << "Mapped slot = " << slot << " nb = " << nb << " to " << lnb << std::endl;
  return lnb;
}


void EvecReader::get_read_geometry (const std::vector < int >&cnodes,
                                    std::map < int,
                                    std::vector < int > >&slots,
                                    std::vector < int >&slot_lvol,
                                    std::vector < int >&lvol,
                                    int64_t & slot_lsites, int &ntotal)
{

  int _nd = (int) cnodes.size ();
  std::vector < int >nodes = cnodes;

  int glb_i[5];
  glb_i[0] = GJP.XnodeSites () * GJP.Xnodes ();
  glb_i[1] = GJP.YnodeSites () * GJP.Ynodes ();
  glb_i[2] = GJP.ZnodeSites () * GJP.Znodes ();
  glb_i[3] = GJP.TnodeSites () * GJP.Tnodes ();
  glb_i[4] = GJP.SnodeSites () * GJP.Snodes ();

  slots.clear ();
  slot_lvol.clear ();
  lvol.clear ();

  int i;
  ntotal = 1;
  int64_t lsites = 1;
  slot_lsites = 1;
  for (i = 0; i < _nd; i++) {
    assert (glb_i[i] % nodes[i] == 0);
    slot_lvol.push_back (glb_i[i] / nodes[i]);
    lvol.push_back (glb_i[i] / GJP.Nodes (i));
    lsites *= lvol.back ();
    slot_lsites *= slot_lvol.back ();
    ntotal *= nodes[i];
  }

  std::vector < int >lcoor, gcoor, scoor;
  lcoor.resize (_nd);
  gcoor.resize (_nd);
  scoor.resize (_nd);

  // create mapping of indices to slots
  for (int lidx = 0; lidx < lsites; lidx++) {
    Lexicographic::CoorFromIndex (lcoor, lidx, lvol);
    for (int i = 0; i < _nd; i++) {
      gcoor[i] = lcoor[i] + GJP.NodeCoor (i) * lvol[i];
      scoor[i] = gcoor[i] / slot_lvol[i];
    }
    int slot;
//      printf("%d: slot %d lidx %d\n",UniqueID(),slot,lidx);
    Lexicographic::IndexFromCoor (scoor, slot, nodes);
    std::map < int, std::vector < int >>::iterator sl = slots.find (slot);
    if (sl == slots.end ()) {
//        printf("%d: slot %d created\n",UniqueID(),slot);
      slots[slot] = std::vector < int >();
    }
    slots[slot].push_back (lidx);
  }

}




int EvecReader::read_compressed_vectors (const char *root,
                                         const char *cdir,
                                         std::vector < OPT * >&dest_all,
                                         int start, int end, float time_out)
{
  std::string path (root);
  const char *fname = "read_compressed_vectors()";
  static Timer timer (cname, fname);
  timer.start ();

  int dest_total = dest_all.size ();
  if (start < 0)
    start = 0;
  if (end < 0)
    end = dest_total;
  VRB.Result (cname, fname, "start=%d end=%d\n", start, end);


//    vector < vector < OPT > >block_data;
//    vector < vector < OPT > >block_data_ortho;
//    vector < vector < OPT > >block_coef;

  int glb_i[5];
  glb_i[0] = GJP.XnodeSites () * GJP.Xnodes ();
  glb_i[1] = GJP.YnodeSites () * GJP.Ynodes ();
  glb_i[2] = GJP.ZnodeSites () * GJP.Znodes ();
  glb_i[3] = GJP.TnodeSites () * GJP.Tnodes ();
  glb_i[4] = GJP.SnodeSites () * GJP.Snodes ();

  vol4d =
    GJP.NodeSites (0) * GJP.NodeSites (1) * GJP.NodeSites (2) *
    GJP.NodeSites (3);
  vol5d = vol4d * GJP.SnodeSites ();

  //.........................Reading eigenvalues ...........................
  time_elapse();
  long nvec = 0;
  if (!UniqueID ()) {
    const std::string filename = path + "/eigen-values.txt";
    FILE *file = fopen (filename.c_str (), "r");
    fscanf (file, "%ld\n", &nvec);
    fclose (file);
  }
  sumArray (&nvec, 1);
  VRB.Result (cname, fname, "dest_all.size() %d  args.neig %d nvec %d \n",
              (int) dest_all.size (), args.neig, nvec);
//    assert (dest_all.size () <= args.neig);
  assert (end <= dest_all.size ());
  assert (nvec >= args.neig);
//    assert (end <= args.neig );
  if (end > args.neig)
    end = args.neig;
  VRB.Result (cname, fname, "start %d end %d \n", start, end);


  //.......................Reading metadata........................

  char hostname[1024];
  sprintf (hostname, "%d", UniqueID ());

  uint32_t nprocessors_real = 1;
//    uint32_t status = 0;
  for (int i = 0; i < 5; i++) {
    nprocessors_real *= GJP.Nodes (i);;
  }

  vector < int >nn (5, 1);
  for (int i = 0; i < 5; i++) {
//      assert (glb_i[i] % args.s[i] == 0);
    if ((glb_i[i] % args.s[i]) != 0)
      ERR.General (cname, fname,
                   "glb_i[%d](%d) should be divided by args.s[%d])(%d)\n", i,
                   glb_i[i], i, args.s[i]);
    nn[i] = glb_i[i] / args.s[i];
  }
  print_time(fname,"reading metadata",time_elapse());

  cps::sync ();

  if (nn[4] != 1) {
    std::cout << "nn[4] != 1. forget Grid -> GJP conversion? \n";
  }
  assert (nn[4] == 1);


  VRB.Result (cname, fname,
              "Reading data that was generated on node-layout: %d %d %d %d %d\n ",
              nn[0], nn[1], nn[2], nn[3], nn[4]);

  int nb_per_node = 1;          //# of blocks per node
  for (int dir = 0; dir < 5; dir++) {
    nb_per_node *= GJP.NodeSites (dir) / args.b[dir];
  }
  cps::sync ();

  size_t f_size_coef_block = args.neig * 2 * args.nkeep;
  //int vol5d = GJP.XnodeSites() * GJP.YnodeSites() *GJP.ZnodeSites() *GJP.TnodeSites() *GJP.SnodeSites();
  f_size = vol5d / 2 * 24;
  size_t f_size_block = f_size / nb_per_node;
  VRB.Result(cname,fname,"nb_per_node=%d f_size_block=%d args.nkeep=%d\n",nb_per_node,f_size_block, args.nkeep);

  block_coef.resize (nb_per_node);
  for (int i = 0; i < nb_per_node; i++) {
    block_coef[i].resize (f_size_coef_block);
    memset (&block_coef[i][0], 0, f_size_coef_block * sizeof (float));
  }
  block_data_ortho.resize (nb_per_node);
  for (int i = 0; i < nb_per_node; i++) {
    block_data_ortho[i].resize (f_size_block * args.nkeep);
    memset (&block_data_ortho[i][0], 0,
            f_size_block * args.nkeep * sizeof (float));
  }
  block_data.resize (nb_per_node);
  for (int i = 0; i < nb_per_node; i++) {
    block_data[i].resize (f_size_block);
    memset (&block_data[i][0], 0, sizeof (float) * f_size_block);
  }


  // now get read geometry
  std::map < int, std::vector < int >>slots;
  std::vector < int >slot_lvol, lvol;
  int64_t slot_lsites;
  int ntotal;
  std::vector < int >_nn (nn.begin (), nn.end ());
  get_read_geometry (_nn, slots, slot_lvol, lvol, slot_lsites, ntotal);
  if (UniqueID () == 0)
//    std::cout << "After get_read_geometry()" << endl;
  print_time(fname,"After get_read_geometry()",time_elapse());
  cps::sync ();
//    int _nd = (int) lvol.size ();

  // types
  //typedef typename Field::scalar_type Coeff_t;
  //typedef typename CoarseField::scalar_type CoeffCoarse_t;
  typedef double RealD;

  // slot layout
  int num_dir = 32;
  int nperdir = ntotal / num_dir;
  if (ntotal % num_dir)
    nperdir += 1;
//    if (nperdir < 1) nperdir = 1;
  FILE *f;

  // load all necessary slots and store them appropriately
  for (std::map < int, std::vector < int > >::iterator sl = slots.begin ();
       sl != slots.end (); sl++) {
    std::vector < int >&idx = sl->second;
    int slot = sl->first;
    VRB.Debug (cname, fname, "%d: slot %d len %d \n", UniqueID (), slot,
               idx.size ());
    std::vector < float >rdata;
    vector < uint32_t > slot_checksums (args.blocks * 2 + 1);

    if ((!crc32_checked) && cdir) {
      char checksum_filename[1024];
      sprintf (checksum_filename, "%s/%d_checksum.txt", cdir, slot);
      f = fopen (checksum_filename, "r");
      if (!f) {
        ERR.General (cname, fname, "Reading %s failed.\n", checksum_filename);
//        assert (f);
      }
      for (int line = 0; line < 2 * args.blocks + 1; line++) {
        fscanf (f, "%X\n", &slot_checksums[line]);
      }
      fclose (f);
    }

    char buf[4096];

    // load one slot vector
    sprintf (buf, "%s/%2.2d/%10.10d.compressed", path.c_str (),
             slot / nperdir, slot);
    f = fopen (buf, "rb");
    if (!f) {
      fprintf (stderr, "Node %s cannot read %s\n", hostname, buf);
      fflush (stderr);
      return false;
    }

    uint32_t crc = 0x0;
    off_t size;

    //GridStopWatch gsw;
    //gsw.Start();
    //double t0 = -dclock();

//      fseeko (f, 0, SEEK_END);
//      size = ftello (f);
    fseeko (f, 0, SEEK_SET);


    double t1 = -dclock ();
    //{
    int nsingleCap = args.nkeep_single;

    int64_t _cf_block_size = slot_lsites * 12 / 2 / args.blocks;

#define FP_16_SIZE(a,b)  (( (a) + (a/b) )*2)

    // first read single precision basis vectors
    for (int nb = 0; nb < args.blocks; nb++) {

      int mnb = globalToLocalCanonicalBlock (slot, _nn, nb);
      if (mnb != -1) {
        //read now
        size_t read_size = (size_t) _cf_block_size * 2 * nsingleCap * 4;
//      if(!UniqueID()) std::cout <<"read_size= "<< read_size <<std::endl;
        fseeko (f, read_size * nb, SEEK_SET);
        std::vector < char >raw_in (read_size);
        assert (fread (&raw_in[0], read_size, 1, f) == 1);
        uint32_t crc_comp = crc32_loop (0, (Bytef *) & raw_in[0], read_size);
        if ((!crc32_checked) && cdir)   // turns off if checksum directory is not specified
          if (crc_comp != slot_checksums[nb]) {
            ERR.General (cname, fname,
                         "nb = %d, crc_compute = %X, crc_read[nb] = %X\n", nb,
                         crc_comp, slot_checksums[nb]);
//            assert (crc_comp == slot_checksums[nb]);
          }
        char *ptr = &raw_in[0];

#pragma omp parallel
        {
          std::vector < float >buff (_cf_block_size * 2, 0);
#pragma omp for
          for (int i = 0; i < nsingleCap; i++) {
            char *lptr = ptr + buff.size () * (i) * 4;
            read_floats (lptr, &buff[0], buff.size ());
            memcpy (&block_data_ortho[mnb][i * buff.size ()], &buff[0],
                    buff.size () * sizeof (float));
          }
        }
      }
    }

//#pragma omp barrier
//#pragma omp single
    cps::sync ();
  print_time(fname,"Finished reading block_data_ortho 0-nsingleCap",time_elapse());
//      if (UniqueID () == 0) std::cout << "Finished reading block_data_ortho 0-nsingleCap" << std::endl;

    //ptr = ptr + _cf_block_size * 2*nsingleCap*args.blocks*4;



    // TODO: at this point I should add a checksum test for block_sp(nb,v,v) for all blocks, then I would know that the mapping
    // to blocks is OK at this point; after that ...

    // then read fixed precision basis vectors
//#pragma omp parallel

    for (int nb = 0; nb < args.blocks; nb++) {
      int mnb = globalToLocalCanonicalBlock (slot, _nn, nb);
      if (mnb != -1) {
        size_t read_size = FP_16_SIZE (2 * _cf_block_size,
                                       24) * (size_t) (args.nkeep - nsingleCap);
        off_t seek_size =
          (off_t) _cf_block_size * 2 * nsingleCap * args.blocks * 4 +
          (args.nkeep - nsingleCap) * nb * FP_16_SIZE (2 * _cf_block_size,
                                                       24);
//      if(!UniqueID()) std::cout <<"read_size= "<< read_size <<std::endl;
        assert (seek_size >= 0);
        fseeko (f, seek_size, SEEK_SET);

        std::vector < char >raw_in (read_size);
//        assert (fread (&raw_in[0], read_size, 1, f) == 1);
        fread_check (&raw_in[0], seek_size, read_size, f);
        uint32_t crc_comp = crc32_loop (0, (Bytef *) & raw_in[0], read_size);
        if (cdir)
          if (crc_comp != slot_checksums[nb + args.blocks]) {
            ERR.General (cname, fname,
                         "nb = %d, crc_compute = %X, crc_read[nb+blocks] = %X\n",
                         nb, crc_comp, slot_checksums[nb + args.blocks]);
//            assert (crc_comp == slot_checksums[nb + args.blocks]);
          }
        char *ptr = &raw_in[0];

#pragma omp parallel
        {
          std::vector < float >buff (_cf_block_size * 2, 0);
#pragma omp for
          for (int i = nsingleCap; i < args.nkeep; i++) {
            char *lptr = ptr + FP_16_SIZE (buff.size (), 24) * (i - nsingleCap);
            read_floats_fp16 (lptr, &buff[0], buff.size (), 24);
            memcpy (&block_data_ortho[mnb][i * buff.size ()], &buff[0],
                    buff.size () * sizeof (float));
          }
        }
      }
    }


    //      ptr = ptr + FP_16_SIZE( _cf_block_size*2*(args.nkeep - nsingleCap)*args.blocks, 24 );
    print_time(fname, "Finished reading block_data_ortho nsingleCap-nkeep\n" ,time_elapse());
    cps::sync ();
    VRB.Debug (cname, fname, "sizeof(off_t)=%d sizeof(size_t)=%d\n",
               sizeof (off_t), sizeof (size_t));

    off_t seek_size =
      (off_t) _cf_block_size * 2 * nsingleCap * args.blocks * 4 +
      FP_16_SIZE (_cf_block_size * 2 * (args.nkeep - nsingleCap) *
                  args.blocks, 24);

//faster, but needs more memory
#if TARGET != BGQ

    fseeko (f, seek_size, SEEK_SET);

    size_t read_size =
      (4 * (size_t) args.nkeep_single * 2 +
       FP_16_SIZE ((args.nkeep - args.nkeep_single) * 2,
                   args.FP16_COEF_EXP_SHARE_FLOATS))
      * (args.neig * args.blocks);
    VRB.Flow (cname, fname, "read_size= %d\n", read_size);


    std::vector < char >raw_in (read_size);
    assert (fread (&raw_in[0], read_size, 1, f) == 1);
    uint32_t crc_comp = crc32_loop (0, (Bytef *) & raw_in[0], read_size);
    if (cdir)
      if (crc_comp != slot_checksums[2 * args.blocks]) {
        ERR.General (cname, fname,
                     "Reading block_coef crc_compute = %X, crc_read[2*blocks] = %X\n",
                     crc_comp, slot_checksums[2 * args.blocks]);
//        assert (crc_comp == slot_checksums[2 * args.blocks]);
      }

    char *ptr = &raw_in[0];

#pragma omp parallel
    {
      std::vector < float >buf1 (args.nkeep_single * 2);
      std::vector < float >buf2 ((args.nkeep - args.nkeep_single) * 2);

#pragma omp for
      for (int j = 0; j < args.neig; j++)
        for (int nb = 0; nb < args.blocks; nb++) {
          int ii, oi;
          int mnb = globalToLocalCanonicalBlock (slot, _nn, nb);
          //int mnb = (nb < nb_per_node) ? nb : -1;
          if (mnb != -1) {

            char *lptr = ptr + (4 * buf1.size () + FP_16_SIZE (buf2.size (),
                                                               args.FP16_COEF_EXP_SHARE_FLOATS))
              * (nb + j * args.blocks);
            int l;
            VRB.Debug (cname, fname, "lptr=%p %d %d\n", lptr - ptr,
                       (int) *lptr, (int) *(lptr + 1));
            read_floats (lptr, &buf1[0], buf1.size ());
            //automatically increase lptr
            memcpy (&block_coef[mnb][j * (buf1.size () + buf2.size ())],
                    &buf1[0], buf1.size () * sizeof (float));
            VRB.Debug (cname, fname, "lptr=%p %d %d\n", lptr - ptr,
                       (int) *lptr, (int) *(lptr + 1));
            read_floats_fp16 (lptr, &buf2[0], buf2.size (),
                              args.FP16_COEF_EXP_SHARE_FLOATS);
            memcpy (&block_coef[mnb]
                    [j * (buf1.size () + buf2.size ()) + buf1.size ()],
                    &buf2[0], buf2.size () * sizeof (float));

          }
        }

    }
#else


//#pragma omp parallel
    {
      std::vector < float >buf1 (args.nkeep_single * 2);
      std::vector < float >buf2 ((args.nkeep - args.nkeep_single) * 2);

//#pragma omp for
      for (int j = 0; j < args.neig; j++)
        for (int nb = 0; nb < args.blocks; nb++) {
          int ii, oi;
          int mnb = globalToLocalCanonicalBlock (slot, _nn, nb);
          //int mnb = (nb < nb_per_node) ? nb : -1;
          if (mnb != -1) {

            off_t offset =
              seek_size + (4 * buf1.size () + FP_16_SIZE (buf2.size (),
                                                          args.FP16_COEF_EXP_SHARE_FLOATS))
              * (nb + j * args.blocks);
            std::vector <
              char >raw_tmp (sizeof (float) * (buf1.size () + buf2.size ()));
            fread_check (&raw_tmp[0], offset, raw_tmp.size (), f);
            int l;
            char *lptr = raw_tmp.data ();
            int *iptr = (int *) raw_tmp.data ();
            VRB.Debug (cname, fname, "lptr=%p %d %d\n",
                       (char *) (offset - seek_size), *iptr, *(iptr + 1));
            read_floats (lptr, &buf1[0], buf1.size ());
            //automatically increase lptr
            memcpy (&block_coef[mnb][j * (buf1.size () + buf2.size ())],
                    &buf1[0], buf1.size () * sizeof (float));

            VRB.Debug (cname, fname, "lptr=%p %d %d\n",
                       (char *) (offset - seek_size), *iptr, *(iptr + 1));

            read_floats_fp16 (lptr, &buf2[0], buf2.size (),
                              args.FP16_COEF_EXP_SHARE_FLOATS);
            memcpy (&block_coef[mnb]
                    [j * (buf1.size () + buf2.size ()) + buf1.size ()],
                    &buf2[0], buf2.size () * sizeof (float));

          }
        }

    }
#endif


    fclose (f);
    t1 += dclock ();
    //std::cout << "Processed " << totalGB << " GB of compressed data at " << totalGB/t1<< " GB/s" << std::endl;
//      if (UniqueID () == 0)
    VRB.Result (cname, fname, "Processed compressed data in %g sec\n", t1);
    //}
  }
  cps::sync ();
  for (int i = start; i < end; i++)
    memset (dest_all[i - start], 0, f_size * sizeof (OPT));

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

#pragma omp parallel
  {
    for (int j = start; j < end; j++) {
//      if (UniqueID() == 0) std::cout << "j = " << j << std::endl;

//                      float* dest = &dest_all[ (int64_t)f_size * j ];
      OPT *dest = dest_all[j - start];

      double ta, tb;
      int tid = omp_get_thread_num ();

      //if (!tid)
      ta = dclock ();


#pragma omp for
      for (int nb = 0; nb < nb_per_node; nb++) {

        float *dest_block = &block_data[nb][0];

        {
          // do reconstruction of this block
          memset (dest_block, 0, sizeof (float) * f_size_block);
          for (int i = 0; i < args.nkeep; i++) {
            float *ev_i = &block_data_ortho[nb][(int64_t) f_size_block * i];
            float *coef = &block_coef[nb][2 * (i + args.nkeep * j)];
            caxpy_single (dest_block, *(complex < float >*) coef, ev_i,
                          dest_block, f_size_block);
          }
        }
      }

      if (!tid) {
        tb = dclock ();
        if (j % 100 == 0)
          VRB.Flow (cname, fname, "%d - %g seconds\n", j, tb - ta);
      }
      //int[5] loc_s = {GJP.NodeSites(0), GJP.NodeSites(1), GJP.NodeSites(2), GJP.NodeSites(3),GJP.NodeSites(4)};

#pragma omp for
      for (int idx = 0; idx < vol4d; idx++) {
        int pos[5], pos_in_block[5], block_coor[5];
        index_to_pos (idx, pos, s_l);

        int parity = (pos[0] + pos[1] + pos[2] + pos[3]) % 2;
        if (parity == 1) {

          for (pos[4] = 0; pos[4] < s_l[4]; pos[4]++) {
            pos_to_blocked_pos (pos, pos_in_block, block_coor);

            int bid = pos_to_index (block_coor, nb_l);
            int ii = pos_to_index (pos_in_block, args.b) / 2;
            float *dst = &block_data[bid][ii * 24];

            int co;
            for (co = 0; co < 12; co++) {
              float *out = &dest[get_cps_index (pos, co, s_l)];
              out[0] = dst[2 * co + 0];
              out[1] = dst[2 * co + 1];
            }
          }
        }
      }

    }
  }

  double t1 = dclock ();

  cps::sync ();
  VRB.Flow (cname, fname, "Reconstruct eigenvectors in %g seconds\n", t1 - t0);

  if (0) {
    for (int i = 0; i < 10; i++) {
      std::cout << dest_all[i] << std::endl;
    }
    std::cout << "end decompressing" << endl;
  }
  cps::sync ();
  timer.stop (true);
#undef FP_16_SIZE
  return true;
}

int EvecReader::read_compressed_blocks (const char *root,
                                        const char *cdir, float time_out)
{
  std::string path (root);
  const char *fname = "read_compressed_blocks()";
  static Timer timer (cname, fname);
  timer.start ();


  char *cname = "read_compressed_vectors";

  int glb_i[5];
  glb_i[0] = GJP.XnodeSites () * GJP.Xnodes ();
  glb_i[1] = GJP.YnodeSites () * GJP.Ynodes ();
  glb_i[2] = GJP.ZnodeSites () * GJP.Znodes ();
  glb_i[3] = GJP.TnodeSites () * GJP.Tnodes ();
  glb_i[4] = GJP.SnodeSites () * GJP.Snodes ();

  vol4d =
    GJP.NodeSites (0) * GJP.NodeSites (1) * GJP.NodeSites (2) *
    GJP.NodeSites (3);
  vol5d = vol4d * GJP.SnodeSites ();

  //.........................Reading eigenvalues ...........................
  long nvec = 0;
  if (!UniqueID ()) {
    const std::string filename = path + "/eigen-values.txt";
    FILE *file = fopen (filename.c_str (), "r");
    fscanf (file, "%ld\n", &nvec);
    fclose (file);
  }
  sumArray (&nvec, 1);


  //.......................Reading metadata........................

  char hostname[1024];
  sprintf (hostname, "%d", UniqueID ());

  uint32_t nprocessors_real = 1;
//    uint32_t status = 0;
  for (int i = 0; i < 5; i++) {
    nprocessors_real *= GJP.Nodes (i);;
  }

  vector < int >nn (5, 1);
  for (int i = 0; i < 5; i++) {
//      assert (glb_i[i] % args.s[i] == 0);
    if ((glb_i[i] % args.s[i]) != 0)
      ERR.General (cname, fname,
                   "glb_i[%d](%d) should be divided by args.s[%d])(%d)\n", i,
                   glb_i[i], i, args.s[i]);
    nn[i] = glb_i[i] / args.s[i];
  }

//  cps::sync ();

  if (nn[4] != 1) {
    std::cout << "nn[4] != 1. forget Grid -> GJP conversion? \n";
  }
  assert (nn[4] == 1);


  VRB.Result (cname, fname,
              "Reading data that was generated on node-layout: %d %d %d %d %d\n ",
              nn[0], nn[1], nn[2], nn[3], nn[4]);

  nb_per_node = 1;          //# of blocks per node
  for (int dir = 0; dir < 5; dir++) {
    nb_per_node *= GJP.NodeSites (dir) / args.b[dir];
  }
//  cps::sync ();

  size_t f_size_coef_block = args.neig * 2 * args.nkeep;
  //int vol5d = GJP.XnodeSites() * GJP.YnodeSites() *GJP.ZnodeSites() *GJP.TnodeSites() *GJP.SnodeSites();
  f_size = vol5d / 2 * 24;
  size_t f_size_block = f_size / nb_per_node;

  block_coef.resize (nb_per_node);
  for (int i = 0; i < nb_per_node; i++) {
    block_coef[i].resize (f_size_coef_block);
    memset (&block_coef[i][0], 0, f_size_coef_block * sizeof (float));
  }
  block_data_ortho.resize (nb_per_node);
  for (int i = 0; i < nb_per_node; i++) {
    block_data_ortho[i].resize (f_size_block * args.nkeep);
    memset (&block_data_ortho[i][0], 0,
            f_size_block * args.nkeep * sizeof (float));
  }


  // now get read geometry
  std::map < int, std::vector < int >>slots;
  std::vector < int >slot_lvol, lvol;
  int64_t slot_lsites;
  int ntotal;
  std::vector < int >_nn (nn.begin (), nn.end ());
  get_read_geometry (_nn, slots, slot_lvol, lvol, slot_lsites, ntotal);
  if (UniqueID () == 0)
    std::cout << "After get_read_geometry()" << endl;
  cps::sync ();
//    int _nd = (int) lvol.size ();

  // types
  //typedef typename Field::scalar_type Coeff_t;
  //typedef typename CoarseField::scalar_type CoeffCoarse_t;
  typedef double RealD;

  // slot layout
  int num_dir = 32;
  int nperdir = ntotal / num_dir;
  if (ntotal % num_dir)
    nperdir += 1;
//    if (nperdir < 1) nperdir = 1;
  FILE *f;

  // load all necessary slots and store them appropriately
  for (std::map < int, std::vector < int > >::iterator sl = slots.begin ();
       sl != slots.end (); sl++) {
    std::vector < int >&idx = sl->second;
    int slot = sl->first;
    VRB.Debug (cname, fname, "%d: slot %d len %d \n", UniqueID (), slot,
               idx.size ());
    std::vector < float >rdata;
    vector < uint32_t > slot_checksums (args.blocks * 2 + 1);

    if ((!crc32_checked) && cdir) {
      char checksum_filename[1024];
      sprintf (checksum_filename, "%s/%d_checksum.txt", cdir, slot);
      f = fopen (checksum_filename, "r");
      if (!f) {
        ERR.General (cname, fname, "Reading %s failed.\n", checksum_filename);
//        assert (f);
      }
      for (int line = 0; line < 2 * args.blocks + 1; line++) {
        fscanf (f, "%X\n", &slot_checksums[line]);
      }
      fclose (f);
    }

    char buf[4096];

    // load one slot vector
    sprintf (buf, "%s/%2.2d/%10.10d.compressed", path.c_str (),
             slot / nperdir, slot);
    f = fopen (buf, "rb");
    if (!f) {
      fprintf (stderr, "Node %s cannot read %s\n", hostname, buf);
      fflush (stderr);
      return false;
    }

    uint32_t crc = 0x0;
    off_t size;

    //GridStopWatch gsw;
    //gsw.Start();
    //double t0 = -dclock();

//      fseeko (f, 0, SEEK_END);
//      size = ftello (f);
    fseeko (f, 0, SEEK_SET);


    double t1 = -dclock ();
    //{
    int nsingleCap = args.nkeep_single;

    int64_t _cf_block_size = slot_lsites * 12 / 2 / args.blocks;

#define FP_16_SIZE(a,b)  (( (a) + (a/b) )*2)

    // first read single precision basis vectors
    for (int nb = 0; nb < args.blocks; nb++) {

      int mnb = globalToLocalCanonicalBlock (slot, _nn, nb);
      if (mnb != -1) {
        //read now
        size_t read_size = (size_t) _cf_block_size * 2 * nsingleCap * 4;
//      if(!UniqueID()) std::cout <<"read_size= "<< read_size <<std::endl;
        fseeko (f, read_size * nb, SEEK_SET);
        std::vector < char >raw_in (read_size);
        assert (fread (&raw_in[0], read_size, 1, f) == 1);
        uint32_t crc_comp = crc32_loop (0, (Bytef *) & raw_in[0], read_size);
        if ((!crc32_checked) && cdir)   // turns off if checksum directory is not specified
          if (crc_comp != slot_checksums[nb]) {
            ERR.General (cname, fname,
                         "nb = %d, crc_compute = %X, crc_read[nb] = %X\n", nb,
                         crc_comp, slot_checksums[nb]);
//            assert (crc_comp == slot_checksums[nb]);
          }
        char *ptr = &raw_in[0];

#pragma omp parallel
        {
          std::vector < float >buff (_cf_block_size * 2, 0);
#pragma omp for
          for (int i = 0; i < nsingleCap; i++) {
            char *lptr = ptr + buff.size () * (i) * 4;
            read_floats (lptr, &buff[0], buff.size ());
            memcpy (&block_data_ortho[mnb][i * buff.size ()], &buff[0],
                    buff.size () * sizeof (float));
          }
        }
      }
    }

//#pragma omp barrier
//#pragma omp single
//    cps::sync ();
//      if (UniqueID () == 0) std::cout << "Finished reading block_data_ortho 0-nsingleCap" << std::endl;

    //ptr = ptr + _cf_block_size * 2*nsingleCap*args.blocks*4;



    // TODO: at this point I should add a checksum test for block_sp(nb,v,v) for all blocks, then I would know that the mapping
    // to blocks is OK at this point; after that ...

    // then read fixed precision basis vectors
//#pragma omp parallel

    for (int nb = 0; nb < args.blocks; nb++) {
      int mnb = globalToLocalCanonicalBlock (slot, _nn, nb);
      if (mnb != -1) {
        size_t read_size = FP_16_SIZE (2 * _cf_block_size,
                                       24) * (size_t) (args.nkeep - nsingleCap);
        off_t seek_size =
          (off_t) _cf_block_size * 2 * nsingleCap * args.blocks * 4 +
          (args.nkeep - nsingleCap) * nb * FP_16_SIZE (2 * _cf_block_size,
                                                       24);
//      if(!UniqueID()) std::cout <<"read_size= "<< read_size <<std::endl;
        assert (seek_size >= 0);
        fseeko (f, seek_size, SEEK_SET);

        std::vector < char >raw_in (read_size);
//        assert (fread (&raw_in[0], read_size, 1, f) == 1);
        fread_check (&raw_in[0], seek_size, read_size, f);
        uint32_t crc_comp = crc32_loop (0, (Bytef *) & raw_in[0], read_size);
        if (cdir)
          if (crc_comp != slot_checksums[nb + args.blocks]) {
            ERR.General (cname, fname,
                         "nb = %d, crc_compute = %X, crc_read[nb+blocks] = %X\n",
                         nb, crc_comp, slot_checksums[nb + args.blocks]);
//            assert (crc_comp == slot_checksums[nb + args.blocks]);
          }
        char *ptr = &raw_in[0];

#pragma omp parallel
        {
          std::vector < float >buff (_cf_block_size * 2, 0);
#pragma omp for
          for (int i = nsingleCap; i < args.nkeep; i++) {
            char *lptr = ptr + FP_16_SIZE (buff.size (), 24) * (i - nsingleCap);
            read_floats_fp16 (lptr, &buff[0], buff.size (), 24);
            memcpy (&block_data_ortho[mnb][i * buff.size ()], &buff[0],
                    buff.size () * sizeof (float));
          }
        }
      }
    }


    //      ptr = ptr + FP_16_SIZE( _cf_block_size*2*(args.nkeep - nsingleCap)*args.blocks, 24 );
    VRB.Flow (cname, fname,
              "Finished reading block_data_ortho nsingleCap-nkeep\n");
//    cps::sync ();
    VRB.Debug (cname, fname, "sizeof(off_t)=%d sizeof(size_t)=%d\n",
               sizeof (off_t), sizeof (size_t));

    off_t seek_size =
      (off_t) _cf_block_size * 2 * nsingleCap * args.blocks * 4 +
      FP_16_SIZE (_cf_block_size * 2 * (args.nkeep - nsingleCap) *
                  args.blocks, 24);

//faster, but needs more memory
#if TARGET != BGQ

    fseeko (f, seek_size, SEEK_SET);

    size_t read_size =
      (4 * (size_t) args.nkeep_single * 2 +
       FP_16_SIZE ((args.nkeep - args.nkeep_single) * 2,
                   args.FP16_COEF_EXP_SHARE_FLOATS))
      * (args.neig * args.blocks);
    VRB.Flow (cname, fname, "read_size= %d\n", read_size);


    std::vector < char >raw_in (read_size);
    assert (fread (&raw_in[0], read_size, 1, f) == 1);
    uint32_t crc_comp = crc32_loop (0, (Bytef *) & raw_in[0], read_size);
    if (cdir)
      if (crc_comp != slot_checksums[2 * args.blocks]) {
        ERR.General (cname, fname,
                     "Reading block_coef crc_compute = %X, crc_read[2*blocks] = %X\n",
                     crc_comp, slot_checksums[2 * args.blocks]);
//        assert (crc_comp == slot_checksums[2 * args.blocks]);
      }

    char *ptr = &raw_in[0];

#pragma omp parallel
    {
      std::vector < float >buf1 (args.nkeep_single * 2);
      std::vector < float >buf2 ((args.nkeep - args.nkeep_single) * 2);

#pragma omp for
      for (int j = 0; j < args.neig; j++)
        for (int nb = 0; nb < args.blocks; nb++) {
          int ii, oi;
          int mnb = globalToLocalCanonicalBlock (slot, _nn, nb);
          //int mnb = (nb < nb_per_node) ? nb : -1;
          if (mnb != -1) {

            char *lptr = ptr + (4 * buf1.size () + FP_16_SIZE (buf2.size (),
                                                               args.FP16_COEF_EXP_SHARE_FLOATS))
              * (nb + j * args.blocks);
            int l;
            VRB.Debug (cname, fname, "lptr=%p %d %d\n", lptr - ptr,
                       (int) *lptr, (int) *(lptr + 1));
            read_floats (lptr, &buf1[0], buf1.size ());
            //automatically increase lptr
            memcpy (&block_coef[mnb][j * (buf1.size () + buf2.size ())],
                    &buf1[0], buf1.size () * sizeof (float));
            VRB.Debug (cname, fname, "lptr=%p %d %d\n", lptr - ptr,
                       (int) *lptr, (int) *(lptr + 1));
            read_floats_fp16 (lptr, &buf2[0], buf2.size (),
                              args.FP16_COEF_EXP_SHARE_FLOATS);
            memcpy (&block_coef[mnb]
                    [j * (buf1.size () + buf2.size ()) + buf1.size ()],
                    &buf2[0], buf2.size () * sizeof (float));

          }
        }

    }
#else


//#pragma omp parallel
    {
      std::vector < float >buf1 (args.nkeep_single * 2);
      std::vector < float >buf2 ((args.nkeep - args.nkeep_single) * 2);

//#pragma omp for
      for (int j = 0; j < args.neig; j++)
        for (int nb = 0; nb < args.blocks; nb++) {
          int ii, oi;
          int mnb = globalToLocalCanonicalBlock (slot, _nn, nb);
          //int mnb = (nb < nb_per_node) ? nb : -1;
          if (mnb != -1) {

            off_t offset =
              seek_size + (4 * buf1.size () + FP_16_SIZE (buf2.size (),
                                                          args.FP16_COEF_EXP_SHARE_FLOATS))
              * (nb + j * args.blocks);
            std::vector <
              char >raw_tmp (sizeof (float) * (buf1.size () + buf2.size ()));
            fread_check (&raw_tmp[0], offset, raw_tmp.size (), f);
            int l;
            char *lptr = raw_tmp.data ();
            int *iptr = (int *) raw_tmp.data ();
            VRB.Debug (cname, fname, "lptr=%p %d %d\n",
                       (char *) (offset - seek_size), *iptr, *(iptr + 1));
            read_floats (lptr, &buf1[0], buf1.size ());
            //automatically increase lptr
            memcpy (&block_coef[mnb][j * (buf1.size () + buf2.size ())],
                    &buf1[0], buf1.size () * sizeof (float));

            VRB.Debug (cname, fname, "lptr=%p %d %d\n",
                       (char *) (offset - seek_size), *iptr, *(iptr + 1));

            read_floats_fp16 (lptr, &buf2[0], buf2.size (),
                              args.FP16_COEF_EXP_SHARE_FLOATS);
            memcpy (&block_coef[mnb]
                    [j * (buf1.size () + buf2.size ()) + buf1.size ()],
                    &buf2[0], buf2.size () * sizeof (float));

          }
        }

    }
#endif


    fclose (f);
    t1 += dclock ();
    //std::cout << "Processed " << totalGB << " GB of compressed data at " << totalGB/t1<< " GB/s" << std::endl;
//      if (UniqueID () == 0)
    VRB.Result (cname, fname, "Processed compressed data in %g sec\n", t1);
    //}
  }
  timer.stop(true);
}

int EvecReader::build_evecs (std::vector < OPT * >&dest_all, int start, int end)
{

  const char *fname = "build_evecs()";
  static Timer timer (cname, fname);
  timer.start ();

  int dest_total = dest_all.size ();
  if (start < 0)
    start = 0;
  if (end < 0)
    end = dest_total;
  VRB.Result (cname, fname, "start=%d end=%d\n", start, end);

  VRB.Result (cname, fname, "dest_all.size() %d  args.neig %d \n",
              (int) dest_all.size (), args.neig);
  assert (end <= dest_all.size ());
//  assert (nvec >= args.neig);
  if (end > args.neig)
    end = args.neig;
  VRB.Result (cname, fname, "start %d end %d \n", start, end);
  size_t f_size_block = f_size / nb_per_node;
  VRB.Result (cname, fname, "nb_per_node %d f_size_block %d \n", nb_per_node, f_size_block);

  vector < vector < OPT > >block_data;
  block_data.resize (nb_per_node);
  for (int i = 0; i < nb_per_node; i++) {
    block_data[i].resize (f_size_block);
    memset (&block_data[i][0], 0, sizeof (float) * f_size_block);
  }

  for (int i = start; i < end; i++)
    memset (dest_all[i - start], 0, f_size * sizeof (OPT));

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
  {
    for (int j = start; j < end; j++) {
//      if (UniqueID() == 0) std::cout << "j = " << j << std::endl;

//                      float* dest = &dest_all[ (int64_t)f_size * j ];
      OPT *dest = dest_all[j - start];

      double ta, tb;
      int tid = omp_get_thread_num ();

      //if (!tid)
      ta = dclock ();


#pragma omp parallel for
      for (int nb = 0; nb < nb_per_node; nb++) {

        float *dest_block = &block_data[nb][0];

        {
          // do reconstruction of this block
          memset (dest_block, 0, sizeof (float) * f_size_block);
          for (int i = 0; i < args.nkeep; i++) {
            float *ev_i = &block_data_ortho[nb][(int64_t) f_size_block * i];
            float *coef = &block_coef[nb][2 * (i + args.nkeep * j)];
            caxpy_single (dest_block, *(complex < float >*) coef, ev_i,
                          dest_block, f_size_block);
          }
        }
      }

      if (!tid) {
        tb = dclock ();
        if (j % 100 == 0)
          VRB.Flow (cname, fname, "%d - %g seconds\n", j, tb - ta);
      }
      //int[5] loc_s = {GJP.NodeSites(0), GJP.NodeSites(1), GJP.NodeSites(2), GJP.NodeSites(3),GJP.NodeSites(4)};

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
            float *dst = &block_data[bid][ii * 24];

            int co;
            for (co = 0; co < 12; co++) {
              float *out = &dest[get_cps_index (pos, co, s_l)];
              out[0] = dst[2 * co + 0];
              out[1] = dst[2 * co + 1];
            }
          }
        }
      }

    }
  }

  double t1 = dclock ();

//  cps::sync ();
  VRB.Flow (cname, fname, "Reconstructed eigenvectors in %g seconds\n", t1 - t0);

  timer.stop (true);
#undef FP_16_SIZE
  return true;
}

int EvecReader::dot (std::vector < OPT * >&_rhs, int start, int end, std::vector< Complex > &coef)
{

  const char *fname = "dot()";
  static Timer timer (cname, fname);
  timer.start ();

  if (start < 0)
    start = 0;

  int nrhs = _rhs.size();
  VRB.Result (cname, fname, "rhs.size() %d  args.neig %d \n", nrhs , args.neig);

  if (end <0) 
    end = args.neig;
  if (end > args.neig)
    end = args.neig;
  VRB.Result (cname, fname, "start %d end %d \n", start, end);
  
  int nvec = end-start;
  coef.resize(nrhs*nvec);

  size_t f_size_block = f_size / nb_per_node;
  VRB.Result (cname, fname, "nb_per_node %d f_size_block %d \n", nb_per_node, f_size_block);

  vector < vector < OPT > >block_data;
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

      //if (!tid)
      ta = dclock ();


#pragma omp parallel for
      for (int nb = 0; nb < nb_per_node; nb++) {

        float *dest_block = &block_data[nb][0];

        {
          // do reconstruction of this block
          memset (dest_block, 0, sizeof (float) * f_size_block);
          for (int i = 0; i < args.nkeep; i++) {
            float *ev_i = &block_data_ortho[nb][(int64_t) f_size_block * i];
            float *coef = &block_coef[nb][2 * (i + args.nkeep * j)];
            caxpy_single (dest_block, *(complex < float >*) coef, ev_i,
                          dest_block, f_size_block);
          }
        }
      }

      if (!tid) {
        tb = dclock ();
        if (j % 100 == 0)
          VRB.Flow (cname, fname, "%d - %g seconds\n", j, tb - ta);
      }
      //int[5] loc_s = {GJP.NodeSites(0), GJP.NodeSites(1), GJP.NodeSites(2), GJP.NodeSites(3),GJP.NodeSites(4)};

#pragma omp parallel for
  for (int rhs = 0; rhs < nrhs; rhs++) {
//      OPT *dest = dest_all[j - start];
      OPT *dest = _rhs[rhs];
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

            int co;
            for (co = 0; co < 12; co++) {
              OPT *src = &dest[get_cps_index (pos, co, s_l)];
              double re = src[0] * vec[2 * co + 0] + src[1] * vec[2 * co + 1];
              double im = -src[0] * vec[2 * co + 1] + src[1] * vec[2 * co + 0];
              coef[j-start+nvec*rhs] += Complex(re,im);
//              out[0] = dst[2 * co + 0];
//              out[1] = dst[2 * co + 1];
            }
          }
        }
      }

    }
  }

  double t1 = dclock ();

  cps::sync ();
  VRB.Flow (cname, fname, "Reconstruct and inner product  of eigenvectors in %g seconds\n", t1 - t0);

  cps::sync ();
  timer.stop (true);
#undef FP_16_SIZE
  return true;
}



CPS_END_NAMESPACE
