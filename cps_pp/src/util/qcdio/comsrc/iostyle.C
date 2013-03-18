// -*- mode:c++; c-basic-offset:4 -*-

// a slight modification on ReadLattice class (if not inheritance)
// to enable parallel reading/writing of "Gauge Connection Format" lattice data

// Read the format of Gauge Connection Format
// from QCDSP {load,unload}_lattice format

#include <config.h>

#include <math.h>
#include <util/iostyle.h>
#include <util/qcdio.h>
#include <util/fpconv.h>
#include <util/time_cps.h>
#include <comms/scu.h>

#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <cassert>

using namespace std;
CPS_START_NAMESPACE

// Compute global offset from local offset and node ID.
unsigned long long lcl2glb(unsigned long long lclid,
                           const unsigned glb[5],
                           const unsigned lcl[5],
                           const unsigned node_coor[5],
                           unsigned dim)
{
    unsigned x[5];
    for(int i = 0; i < dim; ++i) {
        x[i] = lclid % lcl[i] + lcl[i] * node_coor[i];
        lclid /= lcl[i];
    }

    unsigned long long glbid = 0;
    for(int i = dim - 1; i >= 0; --i) {
        glbid = glbid * glb[i] + x[i];
    }
    return glbid;
}

// Compute local data offset and node ID from global offset.
void glb2lcl(unsigned long long *node_id,
             unsigned long long *site_id,
             unsigned long long glbid,
             const unsigned glb[5],
             const unsigned lcl[5],
             const unsigned node[5],
             unsigned dim)
{
    unsigned x[5];
    for(int i = 0; i < dim; ++i) {
        x[i] = glbid % glb[i];
        glbid /= glb[i];
    }

    *node_id = 0;
    *site_id = 0;

    for(int i = dim - 1; i >= 0; --i) {
        *node_id = *node_id * node[i] + x[i] / lcl[i];
        *site_id = *site_id * lcl[i]  + x[i] % lcl[i];
    }
}

// Shift local data.
//
// Before calling this function, node 0 must have node (0 + i)'s data.
// After calling this function node 0 has node (0 + i + 1)'s data.
//
// For other nodes it's a bit more complicated, work it out.
//
// So a series of calls to this function with i = 0, 1, 2, ...  let
// each node loop over the entire lattice. Node 0 will go over the
// global lattice in a natural order, other nodes go over the global
// lattice in a different (but still simple) order.
void shift_data(char **o, char **n, unsigned long i,
                size_t size, 
                const unsigned node[5])
{
    // FIXME: Rewrite getPlusData() to remove this restriction.
    assert(size % sizeof(IFloat) == 0);

    for(int j = 0; j < 5; ++j) {
        getPlusData((IFloat *)*n, (IFloat *)*o, size / sizeof(IFloat), j);
        swap(*o, *n);
        if(i % node[j] < node[j] - 1) break;
        i /= node[j];
    }
}

// remap: Do a remapping so the data layout is suitable for I/O.
//
// If from_lcl == true then we are mapping from local storage order to
// global storage order (and thus will store the data to file
// subsequently), otherwise we are mapping from global to local (and
// thus we will recover the data from file subsequently).
void remap(char *out, char *in, char *tmp,
           const unsigned glb[5],
           const unsigned lcl[5],
           const unsigned node[5],
           const unsigned node_coor[5],
           size_t site_size, // bytes per site
           const int dims,
           bool from_lcl)
{
    const char *cname = "";
    const char *fname = "remap()";

    VRB.Result(cname, fname, "Start remapping, dims = %d.\n", dims);

    const unsigned long long lcl_vol =
        (unsigned long long)lcl[0] * lcl[1] * lcl[2] * lcl[3] * lcl[4];
    // x direction is contiguous
    const size_t block_size = site_size * lcl[0];
    const size_t node_size = block_size * lcl[1] * lcl[2] * lcl[3] * lcl[4];

    // How many shifts (i.e., calls to getPlusData()) we need to do.
    // Note: There is no need to shift in the direction of the last
    // dimension.
    unsigned long shifts = 1;
    for(int i = 0; i < dims - 1; ++i) {
        shifts *= node[i];
    }

    unsigned long long mynodeid = 0;
    for(int i = dims - 1; i >= 0; --i) {
        mynodeid = mynodeid * node[i] + node_coor[i];
    }

    char *o = in;
    char *n = tmp;

    for(unsigned long i = 0; i < shifts; ++i) {
        // compute which node's data we have now
        unsigned node_x[5];

        unsigned long long id = i;
        for(int k = 0; k < 5; ++k) {
            node_x[k] = (id + node_coor[k]) % node[k]; id /= node[k];
        }

        // this will only be used if from_lcl == false
        unsigned long long file_node = 0;
        for(int k = dims - 1; k >= 0; --k) {
            file_node = file_node * node[k] + node_x[k];
        }

        // loop over current local data and pick up blocks
        // FIXME: put an OpenMP directive here?
        for(unsigned long long k = 0; k < lcl_vol; k += lcl[0]) {
            if(from_lcl) { // local ==> global
                unsigned long long glbid = lcl2glb(k, glb, lcl, node_x, dims);
                
                if(glbid / lcl_vol == mynodeid) {
                    memcpy(out + glbid % lcl_vol * site_size,
                           o + k * site_size,
                           block_size);
                }
            } else { // global ==> local
                unsigned long long glbid = file_node * lcl_vol + k;
                
                unsigned long long node_id, site_id;
                glb2lcl(&node_id, &site_id, glbid, glb, lcl, node, dims);
                
                if(node_id == mynodeid) {
                    memcpy(out + site_id * site_size,
                           o + k * site_size,
                           block_size);
                }
            }
        } // end for local sites

        shift_data(&o, &n, i, node_size, node);
    } // end for shifts

    if(n == in) {
        memcpy(n, o, node_size);
    }
    VRB.Result(cname, fname, "End remapping.\n");
}

// convert to file format, csum and pdcsum are computed outside.
//
// we have data in msite, need to fill fsite.
void convert2file(char *fsite, char *msite,
                  size_t fsize, size_t msize, int data_per_site,
                  const LatHeaderBase &hd, const DataConversion &dconv,
                  Float *RandSum, Float *Rand2Sum,
                  TempBufAlloc &rng)
{
    if(hd.headerType() == LatHeaderBase::LATTICE_HEADER) { // Gauge
        for(int mu = 0; mu < 4; ++mu) {
            dconv.host2file(fsite + fsize / 4 * mu,
                            msite + msize / 4 * mu,
                            data_per_site / 4);
        }
    } else { // rng
        UGrandomGenerator *ugran = (UGrandomGenerator*)msite;
        ugran->store(rng.IntPtr());
        dconv.host2file(fsite, rng, data_per_site);
        // next rand
        Float rn = ugran->Grand();
        *RandSum += rn;
        *Rand2Sum += rn*rn;
        // recover
        ugran->load(rng.IntPtr());
    }
}


// convert to memory format, csum and pdcsum are computed outside.
//
// Note: now we have data in fsite, need to fill msite.
void convert2mem(char *fsite, char *msite,
                 size_t fsize, size_t msize, int data_per_site,
                 const LatHeaderBase &hd, const DataConversion &dconv,
                 Float *RandSum, Float *Rand2Sum,
                 TempBufAlloc &rng)
{
    if(hd.headerType() == LatHeaderBase::LATTICE_HEADER) { // Gauge
        for(int mu = 0; mu < 4; ++mu) {
            dconv.file2host(msite + msize / 4 * mu,
                            fsite + fsize / 4 * mu,
                            data_per_site / 4);
        }
    } else { // rng
        dconv.file2host(rng, fsite, data_per_site);
        UGrandomGenerator *ugran = (UGrandomGenerator*)msite;
        ugran->load(rng.IntPtr());
        // next rand
        Float rn = ugran->Grand();
        *RandSum += rn;
        *Rand2Sum += rn*rn;
        // recover
        ugran->load(rng.IntPtr());
    }
}

// convert_data: convert data suitable for I/O.
//
// fsize: size per site in file
// msize: size per site in memory
//
// If mem2file == true then we are converting to file format,
// otherwise converting to memory format.
//
// Whether we have data in fdata or mdata depends on the flag
// mem2file.
void convert_data(char *fdata, char *mdata, size_t fsize, size_t msize,
                  int data_per_site, unsigned long long lcl_vol,
                  const LatHeaderBase &hd, const DataConversion &dconv,
                  int dimension, QioArg &qio_arg, 
                  unsigned *csum, unsigned *pdcsum, Float *RandSum, Float *Rand2Sum,
                  bool mem2file)
{
    const char *cname = "";
    const char *fname = "convert_data()";

    VRB.Result(cname, fname, "Start converting data.\n");
    TempBufAlloc rng(data_per_site * dconv.hostDataSize());

    const unsigned glb[5] = {
        qio_arg.XnodeSites() * qio_arg.Xnodes(),
        qio_arg.YnodeSites() * qio_arg.Ynodes(),
        qio_arg.ZnodeSites() * qio_arg.Znodes(),
        qio_arg.TnodeSites() * qio_arg.Tnodes(),
        dimension == 4 ? 1 : qio_arg.SnodeSites() * qio_arg.Snodes(),
    };
    const unsigned lcl[5] = {
        qio_arg.XnodeSites(),
        qio_arg.YnodeSites(),
        qio_arg.ZnodeSites(),
        qio_arg.TnodeSites(),
        dimension == 4 ? 1 : qio_arg.SnodeSites(),
    };
    const unsigned node_coor[5] = {
        qio_arg.Xcoor(),
        qio_arg.Ycoor(),
        qio_arg.Zcoor(),
        qio_arg.Tcoor(),
        qio_arg.Scoor(),
    };

    for(unsigned long long xst = 0; xst < lcl_vol; ++xst) {
        char *fsite = fdata + xst * fsize;
        char *msite = mdata + xst * msize;

        // it's int for backward compatibility
        int global_id = lcl2glb(xst, glb, lcl, node_coor, dimension);

        if(mem2file) { // memory ==> file
            convert2file(fsite, msite, fsize, msize, data_per_site,
                         hd, dconv, RandSum, Rand2Sum, rng);

            *csum += dconv.checksum(fsite, data_per_site);
            *pdcsum += dconv.posDepCsum(fsite, data_per_site, dimension, qio_arg,
                                        -1, global_id);
        } else { // file ==> memory
            // the order of function calls is different!
            *csum += dconv.checksum(fsite, data_per_site);
            *pdcsum += dconv.posDepCsum(fsite, data_per_site, dimension, qio_arg,
                                        -1, global_id);

            convert2mem (fsite, msite, fsize, msize, data_per_site,
                         hd, dconv, RandSum, Rand2Sum, rng);
        }

    } //xst
    VRB.Result(cname, fname, "End converting data.\n");
}

/*********************************************************************/
/* ParallelIO functions ***********************************************/
/*********************************************************************/

// the last three pointers used to return information when loading Lattice Random Generators
int ParallelIO::load(char *data, const int data_per_site, const int site_mem,
		     const LatHeaderBase &hd, const DataConversion &dconv, 
		     const int dimension /* 4 or 5 */,
		     unsigned int *ptrcsum, unsigned int *ptrpdcsum,
		     Float *rand_sum, Float *rand_2_sum)
{
    const char *fname = "load()";

    const size_t chars_per_site = data_per_site * dconv.fileDataSize();

    const unsigned glb[5] = {
        qio_arg.XnodeSites() * qio_arg.Xnodes(),
        qio_arg.YnodeSites() * qio_arg.Ynodes(),
        qio_arg.ZnodeSites() * qio_arg.Znodes(),
        qio_arg.TnodeSites() * qio_arg.Tnodes(),
        dimension == 4 ? 1 : qio_arg.SnodeSites() * qio_arg.Snodes(),
    };
    const unsigned lcl[5] = {
        qio_arg.XnodeSites(),
        qio_arg.YnodeSites(),
        qio_arg.ZnodeSites(),
        qio_arg.TnodeSites(),
        dimension == 4 ? 1 : qio_arg.SnodeSites(),
    };
    const unsigned node[5] = {
        qio_arg.Xnodes(),
        qio_arg.Ynodes(),
        qio_arg.Znodes(),
        qio_arg.Tnodes(),
        qio_arg.Snodes(),
    };
    const unsigned node_coor[5] = {
        qio_arg.Xcoor(),
        qio_arg.Ycoor(),
        qio_arg.Zcoor(),
        qio_arg.Tcoor(),
        qio_arg.Scoor(),
    };

    const unsigned long long lcl_vol = (unsigned long long)lcl[0] * lcl[1] * lcl[2] * lcl[3] * lcl[4];

    char *fdata = new char[lcl_vol * chars_per_site];
    char *rdata = new char[lcl_vol * chars_per_site];
    char * temp = new char[lcl_vol * chars_per_site];

    int error = 0;
    //////////////////////////////////////////////////////////////////////
    // step 1: load data from file
    //
    // Note: we don't need sSpread() because all nodes will
    // participate in loading.
    VRB.Result(cname, fname, "Parallel loading starting\n");
    setConcurIONumber(qio_arg.ConcurIONumber);

    getIOTimeSlot();

    unsigned long long mynodeid = 0;
    for(int i = dimension - 1; i >= 0; --i) {
        mynodeid = mynodeid * node[i] + node_coor[i];
    }

    FILE *fp = fopen(qio_arg.FileName, "rb");
    if(fp == NULL) {
        error = 1;
        ERR.FileR(cname, fname, qio_arg.FileName);
    }
        
    fseek(fp, hd.dataStart() + (long)(mynodeid * lcl_vol * chars_per_site), SEEK_SET);
        
    if(fread(fdata, chars_per_site, lcl_vol, fp) != lcl_vol) {
        error = 1;
        ERR.FileR(cname, fname, qio_arg.FileName);
    }
        
    fclose(fp);
        
    VRB.Result(cname, fname, "Parallel loading finishing\n");

    finishIOTimeSlot();
  
    if(synchronize(error) != 0) {
        ERR.FileR(cname, fname, qio_arg.FileName);
    }

    //////////////////////////////////////////////////////////////////////
    // step 2: do remapping
    remap(rdata, fdata, temp, glb, lcl, node, node_coor, chars_per_site, dimension, false);

    //////////////////////////////////////////////////////////////////////
    // step 3: convert data from file format to memory format
    unsigned int csum = 0;
    unsigned int pdcsum = 0;
    Float RandSum = 0;
    Float Rand2Sum = 0;

    convert_data(rdata, data, chars_per_site, site_mem,
                 data_per_site, lcl_vol,
                 hd, dconv,
                 dimension, qio_arg,
                 &csum, &pdcsum, &RandSum, &Rand2Sum, false);

    //////////////////////////////////////////////////////////////////////
    // step 4: for 4D data we don't need duplicated checksums.
    if(dimension == 4 && qio_arg.Scoor() != 0) {
        csum = pdcsum = 0;
        RandSum = Rand2Sum = 0;
    }

    VRB.Result(cname, fname, "Parallel Loading done!\n");

    if(ptrcsum) *ptrcsum = csum;
    if(ptrpdcsum) *ptrpdcsum = pdcsum;
    if(rand_sum) *rand_sum = RandSum;
    if(rand_2_sum) *rand_2_sum = Rand2Sum;

    delete[] fdata;
    delete[] rdata;
    delete[] temp;

    return 1;
}

int ParallelIO::store(iostream & output,
		      char * data, const int data_per_site, const int site_mem,
		      LatHeaderBase & hd, const DataConversion & dconv,
		      const int dimension /* 4 or 5 */,
		      unsigned int * ptrcsum, unsigned int * ptrpdcsum,
		      Float * rand_sum, Float * rand_2_sum)
{
    const char *fname = "store()";

    const size_t chars_per_site = data_per_site * dconv.fileDataSize();

    const unsigned glb[5] = {
        qio_arg.XnodeSites() * qio_arg.Xnodes(),
        qio_arg.YnodeSites() * qio_arg.Ynodes(),
        qio_arg.ZnodeSites() * qio_arg.Znodes(),
        qio_arg.TnodeSites() * qio_arg.Tnodes(),
        dimension == 4 ? 1 : qio_arg.SnodeSites() * qio_arg.Snodes(),
    };
    const unsigned lcl[5] = {
        qio_arg.XnodeSites(),
        qio_arg.YnodeSites(),
        qio_arg.ZnodeSites(),
        qio_arg.TnodeSites(),
        dimension == 4 ? 1 : qio_arg.SnodeSites(),
    };
    const unsigned node[5] = {
        qio_arg.Xnodes(),
        qio_arg.Ynodes(),
        qio_arg.Znodes(),
        qio_arg.Tnodes(),
        qio_arg.Snodes(),
    };
    const unsigned node_coor[5] = {
        qio_arg.Xcoor(),
        qio_arg.Ycoor(),
        qio_arg.Zcoor(),
        qio_arg.Tcoor(),
        qio_arg.Scoor(),
    };

    const unsigned long long lcl_vol = (unsigned long long)lcl[0] * lcl[1] * lcl[2] * lcl[3] * lcl[4];

    char *fdata = new char[lcl_vol * chars_per_site];

    unsigned int csum = 0, pdcsum = 0;
    Float RandSum = 0, Rand2Sum = 0;
    ////////////////////////////////////////////////////////////////////////
    // converting data
    convert_data(fdata, data, chars_per_site, site_mem,
                 data_per_site, lcl_vol,
                 hd, dconv,
                 dimension, qio_arg,
                 &csum, &pdcsum, &RandSum, &Rand2Sum, true);

    // simple hack to prevent duplicated checksum values.
    if(dimension == 4 && qio_arg.Scoor() != 0) {
        csum = pdcsum = 0;
        RandSum = Rand2Sum = 0;
    }

    char *rdata = new char[lcl_vol * chars_per_site];
    char * temp = new char[lcl_vol * chars_per_site];
    remap(rdata, fdata, temp, glb, lcl, node, node_coor, chars_per_site, dimension, true);

    //////////////////////////////////////////////////////////////////////
    // start parallel writing
    VRB.Result(cname, fname, "Parallel unloading starting\n");
    setConcurIONumber(qio_arg.ConcurIONumber);
    VRB.Result(cname, fname, "ConcurIONumber = %d\n", qio_arg.ConcurIONumber);

    int error = 0;

    getIOTimeSlot();

    // this condition differs from read()
    if(dimension == 5 || qio_arg.Scoor() == 0) {
        unsigned long long mynodeid = 0;
        for(int i = dimension - 1; i >= 0; --i) {
            mynodeid = mynodeid * node[i] + node_coor[i];
        }

        int retry = 0;
        do {
            FILE *fp = fopen(qio_arg.FileName, "r+b");
            if(fp == NULL) {
                error = 1;
                goto sync_error;
            }

            fseek(fp, hd.dataStart() + (long)(mynodeid * lcl_vol * chars_per_site), SEEK_SET);

            if(fwrite(rdata, chars_per_site, lcl_vol, fp) != lcl_vol) {
                error = 1;
                goto sync_error;
            }

            fclose(fp);

        sync_error:
            if(error) {
                printf("Node %d: Write error\n", UniqueID());
                ++retry;
            }
        } while (error != 0 && retry < 20);
    }

    finishIOTimeSlot();

    if(synchronize(error) > 0) {
        ERR.FileW(cname, fname, qio_arg.FileName);
    }

    VRB.Result(cname,fname,"Parallel Unloading done!\n");

    if(ptrcsum) *ptrcsum = csum;
    if(ptrpdcsum) *ptrpdcsum = pdcsum;
    if(rand_sum) *rand_sum = RandSum;
    if(rand_2_sum) *rand_2_sum = Rand2Sum;

    delete[] fdata;
    delete[] rdata;
    delete[] temp;

    return 1;
}

#if 1

/*********************************************************************/
/* SerialIO functions ***********************************************/
/*********************************************************************/

int SerialIO::load(char *data, const int data_per_site, const int site_mem,
		   const LatHeaderBase &hd, const DataConversion &dconv, 
		   const int dimension /* 4 or 5 */,
		   unsigned int *ptrcsum, unsigned int *ptrpdcsum,
		   Float *rand_sum, Float *rand_2_sum)
{
    const char *fname = "load()";

    const size_t chars_per_site = data_per_site * dconv.fileDataSize();

    const unsigned glb[5] = {
        qio_arg.XnodeSites() * qio_arg.Xnodes(),
        qio_arg.YnodeSites() * qio_arg.Ynodes(),
        qio_arg.ZnodeSites() * qio_arg.Znodes(),
        qio_arg.TnodeSites() * qio_arg.Tnodes(),
        dimension == 4 ? 1 : qio_arg.SnodeSites() * qio_arg.Snodes(),
    };
    const unsigned lcl[5] = {
        qio_arg.XnodeSites(),
        qio_arg.YnodeSites(),
        qio_arg.ZnodeSites(),
        qio_arg.TnodeSites(),
        dimension == 4 ? 1 : qio_arg.SnodeSites(),
    };
    const unsigned node[5] = {
        qio_arg.Xnodes(),
        qio_arg.Ynodes(),
        qio_arg.Znodes(),
        qio_arg.Tnodes(),
        qio_arg.Snodes(),
    };
    const unsigned node_coor[5] = {
        qio_arg.Xcoor(),
        qio_arg.Ycoor(),
        qio_arg.Zcoor(),
        qio_arg.Tcoor(),
        qio_arg.Scoor(),
    };

    const unsigned long long lcl_vol = (unsigned long long)lcl[0] * lcl[1] * lcl[2] * lcl[3] * lcl[4];

    char *fdata = new char[lcl_vol * chars_per_site];
    char *rdata = new char[lcl_vol * chars_per_site];
    char * temp = new char[lcl_vol * chars_per_site];

    //////////////////////////////////////////////////////////////////////
    // step 1: read from the file
    unsigned long shifts = 1;
    for(int i = 0; i < dimension; ++i) {
        shifts *= node[i];
    }

    FILE *fp = Fopen(qio_arg.FileName, "r");
    if(fp == NULL) {
        ERR.FileW(cname, fname, qio_arg.FileName);
    }
    if(UniqueID() == 0) {
        fseek(fp, hd.dataStart(), SEEK_SET);
    }

    VRB.Result(cname, fname, "Serial loading <thru node 0> starting\n");

    char *o = fdata;
    char *n = temp;

    int error = 0;
    for(unsigned long i = 0; i < shifts; ++i) {
        if(UniqueID() == 0) {
            size_t ret = fread(o, chars_per_site, lcl_vol, fp);
            if(ret != lcl_vol) {
                error = 1;
                ERR.FileR(cname, fname, qio_arg.FileName);
            }
        }

        shift_data(&o, &n, i, lcl_vol * chars_per_site, node);
    }
    Fclose(fp);

    VRB.Result(cname, fname, "Serial loading <thru node 0> finishing\n");

    if(synchronize(error) > 0) {
        ERR.FileR(cname, fname, qio_arg.FileName);
    }

    //////////////////////////////////////////////////////////////////////
    // step 2: spread (clone) lattice data along s-dim
    if(dimension == 4) {
        sSpread(o, chars_per_site * qio_arg.VolNodeSites());
    }
  
    //////////////////////////////////////////////////////////////////////
    // step 3: do remapping
    remap(rdata, o, n, glb, lcl, node, node_coor, chars_per_site, dimension, false);

    //////////////////////////////////////////////////////////////////////
    // step 4: convert data from file format to memory format
    unsigned int csum = 0;
    unsigned int pdcsum = 0;
    Float RandSum = 0;
    Float Rand2Sum = 0;

    convert_data(rdata, data, chars_per_site, site_mem,
                 data_per_site, lcl_vol,
                 hd, dconv,
                 dimension, qio_arg,
                 &csum, &pdcsum, &RandSum, &Rand2Sum, false);

    //////////////////////////////////////////////////////////////////////
    // step 5: for 4D data we don't need duplicated checksums.
    if(dimension == 4 && qio_arg.Scoor() != 0) {
        csum = pdcsum = 0;
        RandSum = Rand2Sum = 0;
    }

    VRB.Result(cname, fname, "Serial Loading done!\n");

    if(ptrcsum) *ptrcsum = csum;
    if(ptrpdcsum) *ptrpdcsum = pdcsum;
    if(rand_sum) *rand_sum = RandSum;
    if(rand_2_sum) *rand_2_sum = Rand2Sum;

    delete[] fdata;
    delete[] rdata;
    delete[] temp;

  VRB.FuncEnd(cname,fname);
  return 1;
}

// data_per_site : how many numbers (floating point numbers for gauge
// field or unsigned integers for RNG) per site we will store in the
// file.
//
// site_mem : The size of per site data in bytes in memory (e.g. for gauge
// field this is always sizeof(Matrix) * 4).
int SerialIO::store(iostream &output, char *data,
                    const int data_per_site, const int site_mem,
		    LatHeaderBase &hd, const DataConversion &dconv,
		    const int dimension /* 4 or 5 */,
		    unsigned int *ptrcsum, unsigned int *ptrpdcsum,
		    Float *rand_sum, Float *rand_2_sum)
{ 
    const char *fname = "store()";
  
    const size_t chars_per_site = data_per_site * dconv.fileDataSize();

    const unsigned glb[5] = {
        qio_arg.XnodeSites() * qio_arg.Xnodes(),
        qio_arg.YnodeSites() * qio_arg.Ynodes(),
        qio_arg.ZnodeSites() * qio_arg.Znodes(),
        qio_arg.TnodeSites() * qio_arg.Tnodes(),
        dimension == 4 ? 1 : qio_arg.SnodeSites() * qio_arg.Snodes(),
    };
    const unsigned lcl[5] = {
        qio_arg.XnodeSites(),
        qio_arg.YnodeSites(),
        qio_arg.ZnodeSites(),
        qio_arg.TnodeSites(),
        dimension == 4 ? 1 : qio_arg.SnodeSites(),
    };
    const unsigned node[5] = {
        qio_arg.Xnodes(),
        qio_arg.Ynodes(),
        qio_arg.Znodes(),
        qio_arg.Tnodes(),
        qio_arg.Snodes(),
    };
    const unsigned node_coor[5] = {
        qio_arg.Xcoor(),
        qio_arg.Ycoor(),
        qio_arg.Zcoor(),
        qio_arg.Tcoor(),
        qio_arg.Scoor(),
    };

    const unsigned long long lcl_vol = (unsigned long long)lcl[0] * lcl[1] * lcl[2] * lcl[3] * lcl[4];

    char *fdata = new char[lcl_vol * chars_per_site];

    unsigned csum = 0, pdcsum = 0;
    Float RandSum = 0, Rand2Sum = 0;
    //////////////////////////////////////////////////////////////////////
    // start converting data
    convert_data(fdata, data, chars_per_site, site_mem,
                 data_per_site, lcl_vol,
                 hd, dconv,
                 dimension, qio_arg,
                 &csum, &pdcsum, &RandSum, &Rand2Sum, true);

    // simple hack to prevent duplicated checksum values.
    if(dimension == 4 && qio_arg.Scoor() != 0) {
        csum = pdcsum = 0;
        RandSum = Rand2Sum = 0;
    }

    char *rdata = new char[lcl_vol * chars_per_site];
    char * temp = new char[lcl_vol * chars_per_site];
    remap(rdata, fdata, temp, glb, lcl, node, node_coor, chars_per_site, dimension, true);

    VRB.Result(cname, fname, "Serial unloading <thru node 0> starting\n");
    
    FILE *fp = Fopen(qio_arg.FileName, "ab");

    unsigned long shifts = 1;
    for(int i = 0; i < dimension; ++i) {
        shifts *= node[i];
    }

    char *o = rdata;
    char *n = temp;

    int error = 0;
    for(unsigned long i = 0; i < shifts; ++i) {
        size_t ret = Fwrite(o, chars_per_site, lcl_vol, fp);
        if(ret != lcl_vol) {
            error = 1;
            goto sync_error;
        }
        Fflush(fp);

        shift_data(&o, &n, i, lcl_vol * chars_per_site, node);
    }

    Fclose(fp);

    // if(isNode0())     output.close();
    // if(!input.good()) error = 1;

    VRB.Result(cname, fname, "Serial unloading <thru node 0> finishing\n");

 sync_error:
    if(synchronize(error)>0) 
        ERR.FileW(cname,fname,qio_arg.FileName);

    VRB.Result(cname, fname, "Serial Unloading done!\n");

    if(ptrcsum) *ptrcsum = csum;
    if(ptrpdcsum) *ptrpdcsum = pdcsum;
    if(rand_sum) *rand_sum = RandSum;
    if(rand_2_sum) *rand_2_sum = Rand2Sum;


    delete[] fdata;
    delete[] rdata;
    delete[] temp;

    return 1;
}

// NOTE: !!!
// the x-shift is a little different from y-, z-, & t-shift.
// x-shift shift 1 NODE, while others shift a SITE
void SerialIO::xShiftNode(char * data, const int xblk, const int dir) const {

  if(qio_arg.Xnodes() <= 1) return;

  if(isRow0 ()) {  // x rotation only apply to nodes (0,0,0,x)
    //    const SCUDir pos_dir[] = { SCU_XP, SCU_YP, SCU_ZP, SCU_TP };
    //    const SCUDir neg_dir[] = { SCU_XM, SCU_YM, SCU_ZM, SCU_TM };
    
    char * sendbuf = data;
    VRB.Func(cname,"xShift");

    int fsize = xblk/sizeof(IFloat);
    if (xblk%sizeof(IFloat)>0) fsize++;
//    int  *tmp_p = (int *)sendbuf;
    TempBufAlloc  recvbuf (fsize*sizeof(IFloat));
    if(dir>0) 
      getMinusData(recvbuf.FPtr(),(IFloat *)sendbuf,fsize,0);
    else
      getPlusData(recvbuf.FPtr(),(IFloat *)sendbuf,fsize,0);
//    tmp_p = (int *)recvbuf.FPtr();
    
    memcpy(sendbuf, recvbuf.FPtr(), xblk);
  }
}

void SerialIO::yShift(char * data, const int xblk, const int dir) const {
  int useSCU = 1;
  if(qio_arg.Ynodes() <= 1) useSCU = 0;
    VRB.Func(cname,"yShift");

  if(isFace0()) {
    int fsize = xblk/sizeof(IFloat);
    if (xblk%sizeof(IFloat)>0) fsize++;
    TempBufAlloc  sendbuf (fsize*sizeof(IFloat));
    TempBufAlloc  recvbuf (fsize*sizeof(IFloat));
    if(dir>0){ 
      memcpy(sendbuf, data + xblk * (qio_arg.YnodeSites() - 1), xblk);
      getMinusData(recvbuf.FPtr(),sendbuf.FPtr(),fsize,1);
    } else {
      memcpy(sendbuf, data, xblk);
      getPlusData(recvbuf.FPtr(),sendbuf.FPtr(),fsize,1);
    }


#if 0
    SCUDirArg send, recv;

    if(dir>0) {
      if(useSCU) {
	memcpy(sendbuf, data + xblk * (qio_arg.YnodeSites() - 1), xblk);
	send.Init(sendbuf, SCU_YP, SCU_SEND, xblk);
	recv.Init(recvbuf, SCU_YM, SCU_REC,  xblk);
      }	
      else {
	memcpy(recvbuf, data + xblk * (qio_arg.YnodeSites() - 1), xblk);
      }
    }
    else {
      if(useSCU) {
	memcpy(sendbuf, data, xblk);
	send.Init(sendbuf, SCU_YM, SCU_SEND, xblk);
	recv.Init(recvbuf, SCU_YP, SCU_REC,  xblk);
      }
      else {
	memcpy(recvbuf, data, xblk);
      }
    }

    if(useSCU) {
      SCUTrans(&send);
      SCUTrans(&recv);
    }
#endif

    // doing memory move at the same time
    if(dir>0) {
      for(int i=qio_arg.YnodeSites()-1;i>0;i--) 
	memcpy(data + i*xblk, data + (i-1)*xblk, xblk);
      //      memmove(data + xblk, data, xblk * (qio_arg.YnodeSites()-1));
    }
    else{
      for(int i=0;i<qio_arg.YnodeSites()-1;i++)
	memcpy(data+i*xblk, data+(i+1)*xblk, xblk);
      //      memmove(data, data+xblk, xblk * (qio_arg.YnodeSites()-1));
    }

//    if(useSCU) {
//      SCUTransComplete();
//    }

    if(dir>0) {
      memcpy(data, recvbuf, xblk);
    }
    else{
      memcpy(data + xblk*(qio_arg.YnodeSites()-1), recvbuf, xblk);
    }
  }
}


void SerialIO::zShift(char * data, const int xblk, const int dir) const {
  int useSCU = 1;
  if(qio_arg.Znodes() <= 1) useSCU = 0;
    VRB.Func(cname,"zShift");

  if(isCube0()) {
    int yblk = xblk * qio_arg.YnodeSites();

    int fsize = xblk/sizeof(IFloat);
    if (xblk%sizeof(IFloat)>0) fsize++;
    TempBufAlloc  sendbuf (fsize*sizeof(IFloat));
    TempBufAlloc  recvbuf (fsize*sizeof(IFloat));

#if 0
    TempBufAlloc sendbuf(xblk);
    TempBufAlloc recvbuf(xblk);

    SCUDirArg send, recv;
#endif

    char * loface = data;
    char * hiface = data + yblk * (qio_arg.ZnodeSites()-1);

    for(int yc=0;yc<qio_arg.YnodeSites();yc++) {
      if(dir>0){ 
	memcpy(sendbuf, hiface + yc*xblk, xblk);
        getMinusData(recvbuf.FPtr(),sendbuf.FPtr(),fsize,2);
      } else {
	memcpy(sendbuf, loface + yc*xblk, xblk);
        getPlusData(recvbuf.FPtr(),sendbuf.FPtr(),fsize,2);
      }

#if 0
      if(dir>0) {
	if(useSCU) {
	  memcpy(sendbuf, hiface + yc*xblk, xblk);
	  send.Init(sendbuf, SCU_ZP, SCU_SEND, xblk);
	  recv.Init(recvbuf, SCU_ZM, SCU_REC,  xblk);
	}
	else {
	  memcpy(recvbuf, hiface + yc*xblk, xblk);
	}
      }
      else {
	if(useSCU) {
	  memcpy(sendbuf, loface + yc*xblk, xblk);
	  send.Init(sendbuf, SCU_ZM, SCU_SEND, xblk);
	  recv.Init(recvbuf, SCU_ZP, SCU_REC,  xblk);
	}
	else {
	  memcpy(recvbuf, loface + yc*xblk, xblk);
	}
      }

      if(useSCU) {
	SCUTrans(&send);
	SCUTrans(&recv);
      }
#endif

      if(dir > 0) {
	for(int zc=qio_arg.ZnodeSites()-1;zc>0;zc--) {
	  memcpy(data + zc*yblk + yc*xblk, data + (zc-1)*yblk + yc*xblk, xblk);
	}
      }
      else {
	for(int zc=0;zc<qio_arg.ZnodeSites()-1;zc++) {
	  memcpy(data + zc*yblk + yc*xblk, data + (zc+1)*yblk + yc*xblk, xblk);
	}
      }
	  
//      if(useSCU) {
//	SCUTransComplete();
//      }

      if(dir>0) {
	memcpy(loface + yc*xblk, recvbuf, xblk);
      }
      else{
	memcpy(hiface + yc*xblk, recvbuf, xblk);
      }
    }
  }
}

void SerialIO::tShift(char * data, const int xblk, const int dir) const {
  int useSCU = 1;
  if(qio_arg.Tnodes() <= 1) useSCU = 0;
    VRB.Func(cname,"tShift");

  if(isSdim0()) {
    int yblk = xblk * qio_arg.YnodeSites();
    int zblk = yblk * qio_arg.ZnodeSites();

    int fsize = xblk/sizeof(IFloat);
    if (xblk%sizeof(IFloat)>0) fsize++;
    TempBufAlloc  sendbuf (fsize*sizeof(IFloat));
    TempBufAlloc  recvbuf (fsize*sizeof(IFloat));

#if 0
    TempBufAlloc sendbuf(xblk);
    TempBufAlloc recvbuf(xblk);

    SCUDirArg send, recv;
#endif

    char * locube = data;
    char * hicube = data + zblk * (qio_arg.TnodeSites()-1);

    for(int zc=0;zc<qio_arg.ZnodeSites();zc++) {
      char * loface = locube + zc*yblk;
      char * hiface = hicube + zc*yblk;
      for(int yc=0;yc<qio_arg.YnodeSites();yc++) {

        if(dir>0){ 
	  memcpy(sendbuf, hiface + yc*xblk, xblk);
          getMinusData(recvbuf.FPtr(),sendbuf.FPtr(),fsize,3);
        } else {
	  memcpy(sendbuf, loface + yc*xblk, xblk);
          getPlusData(recvbuf.FPtr(),sendbuf.FPtr(),fsize,3);
        }

#if 0
	if(dir>0) {
	  if(useSCU) {
	    memcpy(sendbuf, hiface + yc*xblk, xblk);
	    send.Init(sendbuf, SCU_TP, SCU_SEND, xblk);
	    recv.Init(recvbuf, SCU_TM, SCU_REC,  xblk);
	  }
	  else {
	    memcpy(recvbuf, hiface + yc*xblk, xblk);
	  }
	}
	else {
	  if(useSCU) {
	    memcpy(sendbuf, loface + yc*xblk, xblk);
	    send.Init(sendbuf, SCU_TM, SCU_SEND, xblk);
	    recv.Init(recvbuf, SCU_TP, SCU_REC,  xblk);
	  }
	  else {
	    memcpy(recvbuf, loface + yc*xblk, xblk);
	  }
	}

	if(useSCU) {
	  SCUTrans(&send);
	  SCUTrans(&recv);
	}
#endif

	if(dir > 0) {
	  for(int tc=qio_arg.TnodeSites()-1;tc>0;tc--) {
	    memcpy(data + tc*zblk + zc*yblk + yc*xblk, data + (tc-1)*zblk + zc*yblk + yc*xblk, xblk);
	  }
	}
	else {
	  for(int tc=0;tc<qio_arg.TnodeSites()-1;tc++) {
	    memcpy(data + tc*zblk + zc*yblk + yc*xblk, data + (tc+1)*zblk + zc*yblk + yc*xblk, xblk);
	  }
	}
	 
//	if(useSCU) {
//	  SCUTransComplete();
//	}

	if(dir>0) {
	  memcpy(loface + yc*xblk, recvbuf, xblk);
	}
	else{
	  memcpy(hiface + yc*xblk, recvbuf, xblk);
	}
      }
    }
  }
}

void SerialIO::sShift(char * data, const int xblk, const int dir) const {
  if(qio_arg.Snodes() * qio_arg.SnodeSites() == 1) return;

  int useSCU = 1;
  if(qio_arg.Snodes() <= 1) useSCU = 0;
    VRB.Func(cname,"sShift");

  int yblk = xblk * qio_arg.YnodeSites();
  int zblk = yblk * qio_arg.ZnodeSites();
  int tblk = zblk * qio_arg.TnodeSites();

    int fsize = xblk/sizeof(IFloat);
    if (xblk%sizeof(IFloat)>0) fsize++;
    TempBufAlloc  sendbuf (fsize*sizeof(IFloat));
    TempBufAlloc  recvbuf (fsize*sizeof(IFloat));

#if 0  
  TempBufAlloc sendbuf(xblk);
  TempBufAlloc recvbuf(xblk);
  
  SCUDirArg send, recv;
#endif
  
  char * lohypcb = data;
  char * hihypcb = data + tblk * (qio_arg.SnodeSites()-1);
  
  for(int tc=0;tc<qio_arg.TnodeSites();tc++) {
    char * locube = lohypcb + tc*zblk;
    char * hicube = hihypcb + tc*zblk;
    
    for(int zc=0;zc<qio_arg.ZnodeSites();zc++) {
      char * loface = locube + zc*yblk;
      char * hiface = hicube + zc*yblk;
      for(int yc=0;yc<qio_arg.YnodeSites();yc++) {

        if(dir>0){ 
	  memcpy(sendbuf, hiface + yc*xblk, xblk);
          getMinusData(recvbuf.FPtr(),sendbuf.FPtr(),fsize,4);
        } else {
	  memcpy(sendbuf, loface + yc*xblk, xblk);
          getPlusData(recvbuf.FPtr(),sendbuf.FPtr(),fsize,4);
        }
#if 0
	if(dir>0) {
	  if(useSCU) {
	    memcpy(sendbuf, hiface + yc*xblk, xblk);
	    send.Init(sendbuf, SCU_SP, SCU_SEND, xblk);
	    recv.Init(recvbuf, SCU_SM, SCU_REC,  xblk);
	  }
	  else {
	    memcpy(recvbuf, hiface + yc*xblk, xblk);
	  }
	}
	else {
	  if(useSCU) {
	    memcpy(sendbuf, loface + yc*xblk, xblk);
	    send.Init(sendbuf, SCU_SM, SCU_SEND, xblk);
	    recv.Init(recvbuf, SCU_SP, SCU_REC,  xblk);
	  }
	  else {
	    memcpy(recvbuf, loface + yc*xblk, xblk);
	  }	    
	}

	if(useSCU) {
	  SCUTrans(&send);
	  SCUTrans(&recv);
	}
#endif
	
	if(dir > 0) {
	  for(int sc=qio_arg.SnodeSites()-1;sc>0;sc--) {
	    memcpy(data + sc*tblk     + tc*zblk + zc*yblk + yc*xblk, 
		   data + (sc-1)*tblk + tc*zblk + zc*yblk + yc*xblk,
		   xblk);
	  }
	}
	else {
	  for(int sc=0;sc<qio_arg.SnodeSites()-1;sc++) {
	    memcpy(data + sc*tblk     + tc*zblk + zc*yblk + yc*xblk, 
		   data + (sc+1)*tblk + tc*zblk + zc*yblk + yc*xblk, 
		   xblk);
	  }
	}
	
//	if(useSCU) {
//	  SCUTransComplete();
//	}
		
	if(dir>0) {
	  memcpy(loface + yc*xblk, recvbuf, xblk);
	}
	else{
	  memcpy(hiface + yc*xblk, recvbuf, xblk);
	}
      }
    }
  }
}


#if TARGET == QCDOC
void SerialIO::sSpread(char * data, const int datablk) const {  

  const char *fname = "sSpread()";

  // clone the lattice from s==0 nodes to all s>0 nodes
  // only used in loading process
  if(qio_arg.Snodes() <= 1)  return;

  // eg. 8 nodes
  // step:  i   ii       iii       iv
  //        0   -->  1   -->   2   -->   3 
  //        |     
  //        \/    
  //        7   -->  6   -->   5   -->   4

  SCUDirArg  socket;

  VRB.Flow(cname, fname, "Spread on S dimension:: 0 ==> %d\n", qio_arg.Snodes()-1);

  if(qio_arg.Scoor() == 0) {
    socket.Init(data, SCU_SM, SCU_SEND, datablk);
    SCUTrans(&socket);
    SCUTransComplete();
  }
  else if(qio_arg.Scoor() == qio_arg.Snodes()-1) {
    socket.Init(data, SCU_SP, SCU_REC, datablk);
    SCUTrans(&socket);
    SCUTransComplete();
  }

  // spread simultaneously along +S and -S directions
  int sender[2] = { 0, qio_arg.Snodes()-1};
  int receiver[2] = { 1, sender[1]-1 };

  while(receiver[0] < receiver[1]) {  // send until two directions converge
    synchronize();

    VRB.Flow(cname, fname, "Spread on S dimension:: %d ==> %d\n", sender[0], receiver[0]);
    VRB.Flow(cname, fname, "Spread on S dimension:: %d ==> %d\n", sender[1], receiver[1]);

    if(qio_arg.Scoor() == sender[0]) {
      socket.Init(data, SCU_SP, SCU_SEND, datablk);
      SCUTrans(&socket);
      SCUTransComplete();
    }
    else if(qio_arg.Scoor() == receiver[0]) {
      socket.Init(data, SCU_SM, SCU_REC, datablk);
      SCUTrans(&socket);
      SCUTransComplete();
    }
    else if(qio_arg.Scoor() == sender[1]) {
      socket.Init(data, SCU_SM, SCU_SEND, datablk);
      SCUTrans(&socket);
      SCUTransComplete();
    }
    else if(qio_arg.Scoor() == receiver[1]) {
      socket.Init(data, SCU_SP, SCU_REC, datablk);
      SCUTrans(&socket);
      SCUTransComplete();
    }

    sender[0] = receiver[0];
    sender[1] = receiver[1];
    
    receiver[0]++;
    receiver[1]--;

  }
}
#else
void SerialIO::sSpread(char * data, const int datablk) const {  

  const char *fname = "sSpread()";

  // clone the lattice from s==0 nodes to all s>0 nodes
  // only used in loading process
  if(qio_arg.Snodes() <= 1)  return;

  // eg. 8 nodes
  // step:  i   ii       iii       iv
  //        0   -->  1   -->   2   -->   3 
  //        |     
  //        \/    
  //        7   -->  6   -->   5   -->   4
  TempBufAlloc snd_buf(datablk);
  TempBufAlloc rcv_buf(datablk);
  int s_nodes=qio_arg.Snodes();
  int s_coor=qio_arg.Scoor();
  if(s_coor == 0) memcpy(snd_buf.FPtr(), data, datablk);
  for(int i =0;i<s_nodes-1;i++){
    getMinusData(rcv_buf.FPtr(),snd_buf.FPtr(),datablk/sizeof(IFloat),4);
    if(i+1 == s_coor) memcpy(data, rcv_buf.FPtr(), datablk);
    memcpy(snd_buf.FPtr(), rcv_buf.FPtr(), datablk);
  }

}
#endif



// Testing functions for *Shift() class

int SerialIO::backForthTest() {
  const char * fname = "backForthTest()";

  int error = 0;

  srand(1234);

  int datablk = 4*18*sizeof(Float) * qio_arg.VolNodeSites();
  int xblk = 4*18*sizeof(Float) * qio_arg.XnodeSites();
  char* data = new char[datablk];
  for(int i=0;i<datablk;i++)
    data[i] = rand() % 256 * uniqueID();

  TempBufAlloc buf(xblk);
  memcpy(buf, data, xblk);

  xShiftNode(data,xblk,-1);
  xShiftNode(data,xblk,1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "xShiftNode() error!!\n";
  }

  xShiftNode(data,xblk,1);
  xShiftNode(data,xblk,-1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "xShiftNode() error!!\n";
  }


  yShift(data,xblk,-1);
  yShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "yShift() error!!\n";
  }

  yShift(data,xblk,1);
  yShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "yShift() error!!" << endl << endl;
  }


  zShift(data,xblk,-1);
  zShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "zShift() error!!" << endl << endl;
  }

  zShift(data,xblk,1);
  zShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "zShift() error!!" << endl << endl;
  }


  tShift(data,xblk,-1);
  tShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "tShift() error!!" << endl << endl;
  }

  tShift(data,xblk,1);
  tShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "tShift() error!!" << endl << endl;
  }

  sShift(data,xblk,-1);
  sShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "sShift() error!!" << endl << endl;
  }

  sShift(data,xblk,1);
  sShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) {
    error = 1;
    cout << "sShift() error!!" << endl << endl;
  }

  delete[] data;

  if(error) return 0;
  return 1;
}


int SerialIO::rotateTest() {
  int error = 0;

  srand(3456);

  int datablk = 4*18*sizeof(Float) * qio_arg.VolNodeSites();
  int xblk = 4*18*sizeof(Float) * qio_arg.XnodeSites();

  char * data = new char[datablk];
  for(int i=0;i<datablk;i++) 
    data[i] = rand() % 256 * uniqueID();

  TempBufAlloc buf(xblk);
  memcpy(buf, data, xblk);

  for(int i=0;i<qio_arg.Xnodes();i++) xShiftNode(data,xblk,-1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "xShiftNode() error!!" << endl << endl;
  }
  for(int i=0;i<qio_arg.Xnodes();i++) xShiftNode(data,xblk,1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "xShiftNode() error!!" << endl << endl;
  }

  for(int i=0;i<qio_arg.Ynodes()*qio_arg.YnodeSites();i++) yShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "yShift() error!!" << endl << endl;
  }
  for(int i=0;i<qio_arg.Ynodes()*qio_arg.YnodeSites();i++) yShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "yShift() error!!" << endl << endl;
  }

  for(int i=0;i<qio_arg.Znodes()*qio_arg.ZnodeSites();i++) zShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "zShift() error!!" << endl << endl;
  }
  for(int i=0;i<qio_arg.Znodes()*qio_arg.ZnodeSites();i++) zShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "zShift() error!!" << endl << endl;
  }

  for(int i=0;i<qio_arg.Tnodes()*qio_arg.TnodeSites();i++) tShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "tShift() error!!" << endl << endl;
  }
  for(int i=0;i<qio_arg.Tnodes()*qio_arg.TnodeSites();i++) tShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "tShift() error!!" << endl << endl;
  }

  for(int i=0;i<qio_arg.Snodes()*qio_arg.SnodeSites();i++) sShift(data,xblk,-1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "sShift() error!!" << endl << endl;
  }
  for(int i=0;i<qio_arg.Snodes()*qio_arg.SnodeSites();i++) sShift(data,xblk,1);
  if(memcmp(buf, data, xblk)) { 
    error = 1;
    cout << "sShift() error!!" << endl << endl;
  }

  delete[] data;

  if(error) return 0;
  return 1;

}


#endif // TARGET == QCDOC



CPS_END_NAMESPACE
