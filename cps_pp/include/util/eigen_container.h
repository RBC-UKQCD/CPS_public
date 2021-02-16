#ifndef _EIGEN_CONTAINER_H_
#define _EIGEN_CONTAINER_H_
#include<config.h>
#ifdef USE_GRID
#include<Grid/Grid.h>
#endif

/*
    Declaration/definition for 
     1.  eigen vectors/values container (EigenContainer class),
     2.  its cache mechanism (EigenCache class),
     3.  with the cache list, which is global instance, just like VRB, ERR, GJP  (EigenCacheList)

   EigenCacheList  uses STL's vector class. So, if there is a problem with the template or STL in the enviornment,
   you have to write (or copy) a limited clone of STL vector.

  FIXME :  Later make  eigen_container.C and move the method functions  AND   the instance of EigenCacheList into that
           Currently EigenCacheList needs to be declared in main.C as I don't want to deal with two files in developping stage.

  FIXME : copmpress / decompress  is still in the middle

  FIXME : do the tests for more than one configuration  and mass, cache may still have a bug

*/

#include <vector>

#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/no_arg.h>
#include <util/time_cps.h>
#include <util/qcdio.h>
#ifdef USE_QIO
#include <util/qio_writeGenericFields.h>
#include <util/qio_readGenericFields.h>
#endif
#include <util/dirac_op.h>
#include <util/eig_io.h>
#include <alg/cg_arg.h>
#include <comms/sysfunc_cps.h>
#include <comms/glb.h>


#include <util/vector.h>

#include <stdlib.h>
#include <cstdlib>


CPS_START_NAMESPACE
#define enum_stringfy( ename ) # ename
//-----------------------------------------------------------
//  Some temporay stuff for compress decompress experiments
//-----------------------------------------------------------
//void lanczos_GramSchm (Float * psi, Float ** vec, int Nvec, size_t f_size, Float * alpha);

//----------------------------------------------------------------------------
//class EigenContainer;		// forward declaration
// Not maintained. Probabaly to be deprecated 

//
// A class for eigenvector cache to reduce I/O
//
enum EigenFormat
{ UNDEFINED, QIO, RBC, RBCcomp };
class EigenCache
{
  friend class EigenContainer;
  friend class AlgLanczos;

  char *cname;
  int neig;
  int neig_b; // blocked vectors, not expanded
  size_t d_size; //data size (float? double?);

//  char fname_root_bc[1024];	// cached root_fname
  //int bc[4]; // boundary conditions in 4 dim, those used in the cached eigen vectors

  size_t f_size; //size of evec in units of Float. Template??

  int alloc_flag;		// if the memory for cache is allocated or not
  int eval_cached;		// if the eval is already cached or not
  int read_interval;		// for testing mostly at the moment

    std::vector < int >index;	//index of eigen vector that is cached.
  // if negative imply it's not cached yet.
  // This is a map between index for eigenv and the index of cache array
  // i.e.   if  I == cached_index[i],  then  evec[I]  holds  i-th eigenvalue
  // Indexing of eigenvalue should be same as the original as it's small.
  // This will be useful for future extion like circular buffer
  // but currently cache size is same as original data ( cacheing all the eigen vectors)

  char cache_name[1024];	// name to identify the cache
  EvecReader evec_reader;

public:
  // Constructer null-ing flags, should be called once in the global scope
    EigenCache ():
    cname ("EigenCache"), read_interval (1000), neig (0), neig_b(0),alloc_flag (0),
    eval_cached (0)
  {
//    *fname_root_bc = 0;		// clear file name
//      index = 0;
  }
  std::vector < Float > evals;
  std::vector < Float * >evecs;

  EigenCache (const char *name)
  : cname ("EigenCache"), read_interval (1000), neig (0), alloc_flag (0),
    eval_cached (0)
  {
    strcpy (cache_name, name);

//    *fname_root_bc = 0;		// clear file name
//      index = 0;
  }
  ~EigenCache() { dealloc();}

  int Neig ()
  {
//    assert(neig <= evals.size());
    VRB.Result(cname,"Neig()","neig=%d neig_b=%d,evals=%d\n",neig,neig_b,evals.size());
    return neig;
  }

  int NeigCoarse ()
  {
    VRB.Result(cname,"NeigCoarse()","neig=%d neig_b=%d,evals=%d\n",neig,neig_b,evals.size());
//    assert(neig_b == evals.size());
    return neig_b;
  }

  char *Name ()
  {
//    return fname_root_bc;
    return cache_name;
  }

  // if the arguments are already cached
  int is_cached (char *a_fname_root_bc, int a_neig)
  {
    int ret= (strcmp (cache_name, a_fname_root_bc) == 0 && ( (neig == a_neig) || (neig_b == a_neig)) );
    VRB.Result (cname, "is_cached", "%p, %s %d %d: %s %di: %d \n", this, cache_name, neig, neig_b, a_fname_root_bc, a_neig, ret);
//    printf ( "Node %d: %s is_cached %p, %s %d: %s %d: %d \n", UniqueID(), cname, this, cache_name, neig, a_fname_root_bc, a_neig, ret);
    return ret;
  }

//  void alloc (char *a_fname_root_bc, int a_neig, size_t a_f_size)
  void alloc (size_t a_neig, size_t a_f_size, size_t a_d_size=sizeof(Float))
  {
    const char *fname = "alloc(C*,I,I)";
    VRB.Func (cname, fname);
    VRB.Flow (cname, fname, "cache_name=%s neig=%d f_size=%ld data_size=%d\n",
		cache_name, a_neig, a_f_size,a_d_size);

    // first deallocate if already allocated
    if(alloc_flag) dealloc ();

    f_size = a_f_size;
    neig = a_neig;
    d_size = a_d_size;

//    strcpy (fname_root_bc, a_fname_root_bc);

    evals.resize (neig);
    evecs.resize (neig);
    for (int i = 0; i < neig; i++) {
      evecs[i] =
	(Float *) smalloc (cname, fname, "evecs[]", f_size * d_size);
      VRB.Debug (cname, fname, "evecs[%d]=%p\n", i, evecs[i]);
    }
    index.resize (neig,-1);
//      index = (int*) smalloc(cname,fname,"index", neig*sizeof(int));

    clear ();
    alloc_flag = 1;
  }

  void clear ()
  {
    VRB.Func (cname, "clear()");
    eval_cached = 0;
    for (int i = 0; i < neig; ++i)
      index[i] = -1;
  }

  void dealloc ()
  {
    const char *fname = "dealloc()";
    VRB.Func (cname, fname);
    if (!alloc_flag)
      return;
//    *fname_root_bc = 0;
#if 1
    resize(0);
#else
    for (int i = 0; i < neig; i++)
      sfree (cname, fname, "evecs[i]", evecs[i]);
    neig = 0;
    evals.resize(0);
    evecs.resize(0);
    index.resize (neig);
#endif
    alloc_flag = 0;
    eval_cached = 0;
  }


  void resize (int new_size)
  {
    const char *fname = "resize(i)";
//    VRB.Func (cname, fname);
    if (!alloc_flag)
      return;
    assert (new_size < neig);
    for(int i=(neig-1); i>=new_size;i--)
    sfree (evecs[i]);
    neig=new_size;
    evecs.resize (neig);
    evals.resize (neig);
    index.resize (neig);
  }

  // save eigenvalues into cache 
  void save (Float * lam)
  {
    VRB.Flow (cname, "save(F*)", "here\n");
    moveFloat (evals.data (), lam, neig);
    eval_cached = 1;
  }

  // load eigenvalues from cache
  // return 0 if it's not in the cache
  int load (Float * lam)
  {
    VRB.Flow (cname, "load(F*)", "%d\n", eval_cached);
    if (!eval_cached)
      return 0;
    moveFloat (lam, evals.data (), neig);
    return 1;
  }

  // save eigenvector into cache
  void savevec (int idx, Vector * v)
  {
    if (index[idx] < 0) {
      index[idx] = idx;		// currently only for a full contents support

      int c_idx = index[idx];	// for future extention, like circular buffer 
      size_t data_size = f_size*d_size/sizeof(Float);
      moveFloat ((Float *) evecs[c_idx], (Float *) v, data_size);
    }
  }

  // load eigenvector from cache
  // return 0 if it's not in the cache
  int loadvec (Vector * v, int idx)
  {
    VRB.Flow (cname, "loadvec(V*,I)", "idx %d index %d\n", idx, index[idx]);
    if (index[idx] < 0)
      return 0;
    int c_idx = index[idx];	// for future extention, like circular buffer 
    size_t data_size = f_size*d_size/sizeof(Float);
    moveFloat ((Float *) v, (Float *) (evecs[c_idx]), data_size);
    return 1;
  }

  // just return the pointer in the cache
  Vector *vec_ptr (int idx)
  {
    VRB.Flow (cname, "vec_ptr(index)", "idx %d index %d\n", idx, index[idx]);
    if (!alloc_flag)
      return 0;
    VRB.Debug (cname, "vec_ptr(index)", "idx %d index %d %p\n", idx,
		index[idx], evecs[idx]);
    assert (idx < neig);
    return (Vector *) ((Float *) evecs[idx]);
  }

  // just return the pointer in the cache
  Vector *set_ptr (int idx, Vector *ptr)
  {
    VRB.Flow (cname, "set_ptr(index)", "idx %d index %d\n", idx, index[idx]);
    if (!alloc_flag)
      return 0;
    assert (idx < neig);
    evecs[idx] = (Float *) ptr;
    VRB.Debug (cname, "set_ptr(index)", "idx %d index %d %p\n", idx,
		index[idx], evecs[idx]);
    return (Vector *) ((Float *) evecs[idx]);
  }

  // just return the pointer in the cache, not copy
  // return 0 if it's not in the cache
  Vector *pointer (int idx)
  {
    VRB.Flow (cname, "pointer(index)", "idx %d index %d\n", idx, index[idx]);
    if (index[idx] < 0)
      return 0;
    int c_idx = index[idx];	// for future extention, like circular buffer 
    return (Vector *) (evecs[c_idx]);
//    Vector* ptr = (Vector*)((Float*)evecs + c_idx*f_size);
    //printf("evec: %d %d %x",idx, c_idx, ptr);
//    return ptr;
  }
  Float *eval_address ()
  {
    VRB.Flow (cname, "eval_address()", "\n");
    return evals.data ();
  }
  void set_neig (int n)
  {
    neig = n;
  }
  void set_index (int n)
  {
    index[n] = n;
  }

  int readCompressed (const char *root_, const char *checksum_dir=NULL, int interval=-1)
  {
    const char *fname="read_compressed()";
    if(interval>0) read_interval = interval;
    VRB.Debug(cname,fname,"dir=%s checksum_dir=%s read_interval=%d\n", root_, checksum_dir, read_interval);
//    exit(-42);
    evec_reader.read_metadata (root_,evals);
    for(int i=0;i< evals.size();i++)
    VRB.Result(cname,fname,"eval[%d]=%e\n",i,evals[i]);
    for (int i = 0; i < evecs.size (); i += read_interval) {
      std::vector < float *>evec_f;
      for (int j = 0; j < read_interval && (i + j) < evecs.size (); j++)
	evec_f.push_back ((float *) evecs[i + j]);
      evec_reader.read_compressed_blocks (root_, checksum_dir);
      VRB.Result(cname,fname,"read_compressed_blocks() done\n");
      evec_reader.build_evecs (evec_f, i, i + evec_f.size ());
      float *temp = (float *) evecs[i];
//      printf("Node  %d: read_compressed: evecs[%d][0]=%g\n",UniqueID(),i,*temp);
    }
    for (int i = 0; i < evecs.size (); i++) {
      float *temp = (float *) evecs[i];
      char *c_tmp = (char *) temp;
      Float sum = 0.;
#pragma omp parallel for reduction (+:sum)
      for (size_t ind = 0; ind < (size_t) f_size ; ind++)
	sum += temp[ind] * temp[ind];
      glb_sum (&sum);
      VRB.Result(cname,fname,"evecs[%d][0]=%g sum[%d]=%g\n", i, *temp, i, sum);
    }
    neig_b=evals.size();
    return evecs.size();
  }

  int readCompressedBlocks (const char *root_, const char *checksum_dir=NULL)
  {
    const char *fname="readCompressedBlocks()";
    VRB.Debug(cname,fname,"dir=%s checksum_dir=%s \n", root_, checksum_dir);
    static Timer timer(cname,fname);
    timer.start();
    evec_reader.read_metadata (root_,evals);
    for(int i=0;i< evals.size();i++)
    VRB.Flow(cname,fname,"eval[%d]=%e\n",i,evals[i]);
      evec_reader.read_compressed_blocks (root_, checksum_dir);
    neig_b = evals.size();
    timer.stop(true);
//    VRB.Result(cname,fname,"%s done\n",fname);
    return 1;
  }

  int blockProj (std::vector < Float * > &_sol, std::vector < Float * >&_rhs, int start, int end, std::vector< Complex > &coef)
  {
    const char *fname="blockProj()";
    static Timer timer(cname,fname);
    timer.start(false);
    evec_reader.proj(_sol,_rhs,start,end,coef);
    timer.stop(true);
    return 1;
  }
  int blockProj2 (std::vector < Float * > &_sol, std::vector < Float * >&_rhs, int start, int end, std::vector< Complex > &coef)
  {
    const char *fname="blockProj2()";
    static Timer timer(cname,fname);
    timer.start(false);
    evec_reader.proj2(_sol,_rhs,start,end,coef);
    timer.stop(true);
    return 1;
  }
};

#ifdef USE_GRID
template < class Field > class EigenCacheGrid:public EigenCache {
private:
  char *cname;
//  const char *fname;
public:

  Grid::GridBase * grid;

  std::vector < Field > evec_grid;
EigenCacheGrid ():grid (NULL), cname ("EigenCacheGrid") {
  }
  EigenCacheGrid (char *name):EigenCache (name)
  {
  }
  ~EigenCacheGrid () { }

};
#endif


//----------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------//
// Dear Chulwoo :
//  Here is the extern declaration  of a global instance, EigenCacheList
//
//  The cache list needes to live after the dirac operator destory, as the same eigen vector
//  would be needed again (at least 12 times for the propagator, and more for other source locations)
//  This is why currently this is a global variable, just like GJP.
//
//  The cleanup (deallocation) of potentially large size memory of cache 
//  is a responsibility of writers of Alg  or main.C, as the appropriate timing for deallocation is unknown. 
//  For example, one would need more than one quark mass or boundary condition at a time in LMA, so it will need
//  more than one instance of a cache  (that's why there is the global vector<Ecache*>, and which timing each of cache
//  becomes useless is hard to know in this low level routine.
//  I don't know how to write a smart automatic garbage collector without perfomance/memory penalties.
//
//  Service routines are provided below also.
//
//------------------------------------------------------------------------------------------------//

extern std::vector < EigenCache * >EigenCacheList;

//Search contents that match to arguments, return 0 if not found
EigenCache *EigenCacheListSearch (char *cache_name, int neig);

// Cleanup list EigenCache, it also destroys contents pointed by the elements.
void EigenCacheListCleanup ();


//----0----------------------------------------------------------------
//
//  I/O etc for eigen value and vectors
//
//  It also contains the cache class defined above
//
//---------------------------------------------------------------------



//extern std::vector<EigenCache*> EigenCacheList;

CPS_END_NAMESPACE
#endif
