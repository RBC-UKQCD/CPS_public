#ifndef _EIGEN_CONTAINER_H_
#define _EIGEN_CONTAINER_H_
#include<config.h>
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
#include <alg/cg_arg.h>

#include <comms/sysfunc_cps.h>
#include <comms/glb.h>


#include <util/vector.h>



CPS_START_NAMESPACE



#define enum_stringfy( ename ) # ename



//-----------------------------------------------------------
//  Some temporay stuff for compress decompress experiments
//-----------------------------------------------------------

void lanczos_GramSchm(Float *psi, Float **vec, int Nvec, int f_size, Float* alpha);

//----------------------------------------------------------------------------
class EigenContainer;  // forward declaration

//
// A class for eigenvector cache to reduce I/O
//
class EigenCache {
  friend class EigenContainer;
  friend class AlgLanczos;

  char* cname;
  int  neig;
 
  char fname_root_bc[1024];  // cached root_fname
  //int bc[4]; // boundary conditions in 4 dim, those used in the cached eigen vectors
  Float* eval; 
  Vector** evec;
  //Float* tmp;

  int f_size;

  int alloc_flag;  // if the memory for cache is allocated or not
  int eval_cached; // if the eval is already cached or not
  int* index; //index of eigen vector that is cached.
              // if negative imply it's not cached yet.
              // This is a map between index for eigenv and the index of cache array
              // i.e.   if  I == cached_index[i],  then  evec[I]  holds  i-th eigenvalue
              // Indexing of eigenvalue should be same as the original as it's small.
              // This will be useful for future extion like circular buffer
              // but currently cache size is same as original data ( cacheing all the eigen vectors)

  char cache_name[1024]; // name to identify the cache

 public:

  // Constructer null-ing flags, should be called once in the global scope
  EigenCache()
    {
      cname="EigenCache";
      *fname_root_bc = 0; // clear file name
      neig = 0;
      alloc_flag = 0;
      eval_cached = 0; // eigen value is not cached yet
      index = 0;
    }

  EigenCache(char* name)
    {
      cname="EigenCache";
      
      strcpy(cache_name,name);

      *fname_root_bc = 0; // clear file name
      neig = 0;
      alloc_flag = 0;
      eval_cached = 0; // eigen value is not cached yet
      index = 0;
    }

  int Neig(){return neig;}
  char* Name(){return fname_root_bc;}
  
  // if the arguments are already cached
  int is_cached( char* a_fname_root_bc, int a_neig )
  {

    return 
      strcmp( fname_root_bc, a_fname_root_bc)==0  &&
      neig == a_neig ;
  }

  void alloc( char* a_fname_root_bc, int a_neig, int a_f_size )
    {
      char* fname = "alloc(C*,I,I)";
      VRB.Func(cname,fname);

      // first deallocate if already allocated
      dealloc();

      f_size = a_f_size;
      neig = a_neig;

      strcpy(fname_root_bc, a_fname_root_bc);
  
      eval = (Float*) smalloc(cname,fname, "eval", neig*sizeof(Float) );
      if(eval==0)ERR.General(cname,fname,"eval could not malloced\n");
      evec = (Vector**) smalloc(cname,fname, "evec[]", neig*sizeof(Vector*));
      if(evec==0)ERR.General(cname,fname,"evec could not malloced\n");
      Float* tmp = (Float*) smalloc(cname,fname, "evec[][]", neig*f_size* sizeof(Float));
      if(tmp==0)ERR.General(cname,fname,"tmp could not malloced\n");
      for(int i=0;i<neig;++i)  evec[i] = (Vector*)( tmp+i*f_size );

      index = (int*) smalloc(cname,fname,"index", neig*sizeof(int));
      
      clear( );
      alloc_flag = 1;
    }

  void clear()
  {
    VRB.Func(cname,"clear()");
    eval_cached = 0;    
    for(int i=0; i< neig;++i) index[i]=-1;
  }

  void dealloc()
  {
    char* fname="dealloc()";
    VRB.Func(cname,fname);
    if(! alloc_flag) return;
    neig = 0;
    *fname_root_bc=0;
    //for(int i=0;i<neig;++i) sfree(cname,fname,"evec[i]", evec[i]);
    //sfree(cname,fname,"tmp", tmp);
    sfree(cname,fname,"evec[0]", evec[0]);
    sfree(cname,fname,"evec", evec);
    sfree(cname,fname,"eval", eval);
    sfree(cname,fname,"index", index);
    alloc_flag = 0;
    eval_cached = 0;
  }

  // save eigenvalues into cache 
  void save( Float* lam )
  {
    VRB.Flow(cname,"save(F*)","here");
    moveFloat( eval, lam, neig);
    eval_cached = 1;
  }

  // load eigenvalues from cache
  // return 0 if it's not in the cache
  int load( Float* lam )
  {
    VRB.Flow(cname,"load(F*)", "%d\n",eval_cached);
    if (! eval_cached ) return 0;
    moveFloat( lam,  eval,  neig );
    return 1;
  }

  // save eigenvector into cache
  void savevec( int idx, Vector* v )
  {
    if ( index[idx] < 0 ) {
      index[idx] = idx;   // currently only for a full contents support

      int c_idx = index[idx]; // for future extention, like circular buffer 
      moveFloat((Float*)(evec[c_idx]), (Float*)v, f_size);
    }
  }

  // load eigenvector from cache
  // return 0 if it's not in the cache
  int loadvec(  Vector* v, int idx )
  {
    VRB.Flow(cname,"loadvec(V*,I)", "idx %d index %d\n",idx, index[idx]);
    if ( index[idx] < 0 )  return 0;
    int c_idx = index[idx]; // for future extention, like circular buffer 
    moveFloat( (Float*)v, (Float*)(evec[c_idx]), f_size);
    return 1;
  }

  // just return the pointer in the cache, not copy
  // return 0 if it's not in the cache
  Vector* pointer(  int idx )
  {
    VRB.Flow(cname,"poiner(I)", "idx %d index %d\n",idx, index[idx]);
    if ( index[idx] < 0 )  return 0;
    int c_idx = index[idx]; // for future extention, like circular buffer 
    return (Vector*)(evec[c_idx]);
  }
  Vector** evec_address()
  {
    VRB.Flow(cname,"evec_address()", "\n");
    return evec;
  }
  Float* eval_address()
  {
    VRB.Flow(cname,"eval_address()", "\n");
    return eval;
  }
  void set_neig(int n)
  {
    neig=n;
  }
};

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

extern std::vector<EigenCache*> EigenCacheList;

//Search contents that match to arguments, return 0 if not found
EigenCache* EigenCacheListSearch( char* fname_root_bc, int neig );

// Cleanup list EigenCache, it also destroies contents pointed by the elements.
void EigenCacheListCleanup( );


//---------------------------------------------------------------------
//
//  I/O etc for eigen value and vectors
//
//  It also contains the cache class defined above
//
//---------------------------------------------------------------------


class EigenContainer {
 private:
  char* cname;
  char* fname;

  int neig;  // number of eigenvalue

  Float* eval; // eigen values (assumes hermitian)

  Vector* evec; // one eigen vector

  int f_size; // size of one vector in the unit of Float
  int f_stride_size; // size of save_stride vectors in the unit of Float

  int f_size_per_site;
  int n_fields;

  // stride is the number of Vectors in one eigenvector
  int stride;

  Lattice* lattice;

  char fname_root_bc[1024]; // the preamble part of the eigen file with the boundary condition appended
  // example
  //    fname_root.bc0001.nev0001  :  1st eigenvector
  //    fname_root.bc0001.eval     :  eigenvalue list
  //  bc0001 is the spacial periodic boundaryconditions with antiperiodic temporal direction

  double mass; // useful for automatic check of eigen value equation

  EigenCache* ecache ;  // pointer to cache class

 public:
  // n_fields : the number of the most outer loop in the file format
  //            it is GJP.Ssites()  for DWF's eigen value.
  // mass : quark mass needed for consistency check for the eigen value equation
  // a_ecache is the cache class provided from outside, if NULL cache menism will be turned off

  // Dear Chulwoo :
  //   You may be puzzled why I needed n_fields.
  //   Since QIO, which is flexible I/O for parallel enviornment, needs to know Float number per site
  //   and the number of these space-time volume chunk in the outer loop, separately
  //   This is why I needed to set up three different integers in the argument. Perhpas there is a better way....
  //

 EigenContainer(Lattice & latt, char* a_fname_root_bc, 
		int neig_, int f_size_per_site_, int n_fields_, 
		EigenCache* a_ecache = 0)
   : cname( "EigenContainer" ),
     fname( "EigenContainer(...)" ),
     lattice( &latt )
    {

      f_size_per_site = f_size_per_site_;
      n_fields = n_fields_;
      f_size = f_size_per_site*n_fields*GJP.VolNodeSites();
      neig = neig_;
      f_stride_size = f_size_per_site*n_fields*GJP.SaveStride()*GJP.VolNodeSites();
      stride = (f_size_per_site / 3) * n_fields * GJP.VolNodeSites() / 2;

      strcpy(fname_root_bc, a_fname_root_bc);

      ecache = 0;
      if( a_ecache ) //if cache system is requested
	{
	  ecache = a_ecache;
	  // if this is not same as already cached file
	  // This should be alloced outside already!
	  if( ! ecache-> alloc_flag ) {
	    ERR.General(cname,fname,"Dummy! requested cache should be alloced!\n");
	    //ecache-> alloc(  fname_root_bc,  neig,  f_size );
	  }
	}
      eval = (Float*) smalloc(cname,fname, "eval", neig*sizeof(Float) );
      if(eval==0)ERR.General(cname,fname,"eval could not malloced\n");
      evec = (Vector*) smalloc(cname,fname, "evec", f_stride_size* sizeof(Float));
      if(evec==0)ERR.General(cname,fname,"evec could not be malloced\n");
    }

  ~EigenContainer()
    {
      fname="~EigenContainer()";
      sfree(cname,fname, "eval",eval);
      sfree(cname,fname, "evec",evec);
    }

  // load eigen values from  fname_root.s-eval
  Float* load_eval( )
  {
    char* fname="load_eval()";

    FILE* fp;
    // try cache fist 
    if( ecache )
      if( ecache-> load( eval ) ) return eval;

    char file[1024];
    // formerly extension was .eval, which contains eigenvalue of MATPCDAG_MATPC, now .evals contains the eigenvalue of MAPC_HERM
    snprintf(file, 1024, "%s.evals",fname_root_bc);
    fp = Fopen(file,"r");
    if(!fp) 
      ERR.General(cname,fname,"failed opening %s, please note .eval (MATPCDAG_MATPC) and .evals (MAT_HERM) have different contents.",file);

    for(int i=0; i< neig; ++i ){
      //printf("NEIG %d %d\n",neig,i);
      int idx;
      Float e=0;      
      if(!UniqueID()) {
	fscanf( fp, "%d %lf", &idx, &e );
	if( idx != i) ERR.General(cname,fname,"In %s file index is wrong %d vs %d\n", file, idx, i);
      }
      glb_sum(&e);
      eval[i]=e;
    }
    Fclose(fp);

    if( ecache )
      ecache-> save( eval );
    return eval;
  }

  // save eigen values from  fname_root.s-eval 
  void save_eval( Float* in_eval )
  {
    char* fname="save_eval()";

      if(!UniqueID()){

	char file[1024];
	// formerly extension was .eval, which contains eigenvalue of MATPCDAG_MATPC, now .evalS contains the eigenvalue of MAPC_HERM
	snprintf(file, 1024, "%s.evals",fname_root_bc); 

	FILE* fp;
	fp = Fopen(file,"w");
	if(!fp) {
	  ERR.General(cname,fname,"failed opening %s, please note .eval (MATPCDAG_MATPC) and .evals (MAT_HERM) have different contents.",file);
	}

	for(int i=0; i< neig; ++i ){
	  Fprintf( fp, "%d %.16e\n", i, in_eval[i] );
	}
	Fclose(fp);

      }
  }

#ifndef USE_QIO
  // load from file with the specified "nev" index
  Vector* nev_load( int index )
  {

    ERR.NotImplemented(cname,"nev_load(I)");
    return NULL;
  }

  // save to file with the specified "nev" index
  void nev_save( int index, Vector* evec_,
		 char* field_type_label,
		 char* ensemble_id="n/a",
		 char* ensemble_label="n/a",
		 int seqNum=777 )
  {
    ERR.NotImplemented(cname,"nev_save(I)");
  }
#else

  // load from file with the specified "nev" index
  Vector* nev_load( int index )
  {

    VRB.Flow(cname,"nev_load(I)","ecache %x \n", ecache);
    if (ecache){ // cached, don't read in again
      VRB.Flow(cname,"nev_load(I)"," index %d ecache->index[index] %d \n",
	       index,ecache->index[index]);
      
      if ( ecache-> index[index] >= 0) {
	if(!UniqueID()) printf("Getting cached eig-vec %d\n",index);
	// copy version
	//ecache-> loadvec ( evec,  index );
	//return evec;
	// no-copy version
	return ecache->pointer( index );
      }
    }

    char file[1024];
    int save_stride = GJP.SaveStride();
    int num = (index/save_stride)*save_stride;
    snprintf(file,1024, "%s.nev%03d", fname_root_bc, num);
    
    qio_readGenericFields readGenField;
    readGenField. read_genericfields( file, save_stride*n_fields, f_size_per_site, evec, QIO_UNKNOWN);
    // convert *to* stag ordering
    if(lattice->Fclass()==F_CLASS_P4){
      for(int i=0;i<save_stride;i++){
	lattice->Fconvert(evec+i*stride,STAG,CANONICAL,1);
      }
    }

    // evec is type Vector, with at least one 3-complex vector at each site and working on checkerboard
    if( ecache ){
      for(int i=0;i<save_stride;i++){
	ecache-> savevec( num+i, evec+i*stride );
	if(!UniqueID()) printf("Cached eig-vec %d (size=%d)\n",num+i,stride);
      }
    }
    printf("returning eig-vec %d\n",index % save_stride);
    return evec+(index % save_stride)*stride;
  }

  // save to file with the specified "nev" index
  void nev_save( int index, Vector* evec_,
		 char* field_type_label,
		 char* ensemble_id="n/a",
		 char* ensemble_label="n/a",
		 int seqNum=777 )
  {
    double time=dclock();
    
    char file[1024];
   
    snprintf(file,1024, "%s.nev%03d", fname_root_bc, index);

    qio_writeGenericFields writeGenField;

    int save_stride = GJP.SaveStride();
    if(lattice->Fclass()==F_CLASS_P4){
      for(int i=0;i<save_stride;i++){
	lattice->Fconvert(evec_+i*stride,CANONICAL,STAG,1);
      }
    }
    writeGenField.setHeader( ensemble_id, ensemble_label, seqNum, field_type_label );
    writeGenField. write_genericfields( file, save_stride*n_fields, f_size_per_site, evec_);

    if(!UniqueID()) printf("nev_save, time to save :%e sec\n",time-dclock());
  }
#endif


  void nev_check( Vector* vtmp, Float mass, Float* residual = 0, Float* eval = 0)
  {
    char* fname="nev_check(V*,F,F*,F*)";

    CgArg cg_arg;
    cg_arg.mass = mass;
    cg_arg.eigen_shift = 0.0;

    Vector* Apsi = (Vector*) smalloc(cname,fname, "Apsi", f_size* sizeof(Float));
    if(Apsi==0)ERR.General(cname,fname,"Apsi could not malloced\n");

    if(lattice->Fclass()==F_CLASS_DWF){
      cg_arg.RitzMatOper = MATPC_HERM; //could be  MATPCDAG_MATPC;
      DiracOpDwf dop( *lattice, 0, 0, &cg_arg, CNV_FRM_NO );
      dop.RitzMat(Apsi, vtmp );
    }else if(lattice->Fclass()==F_CLASS_MOBIUS){
      cg_arg.RitzMatOper = MATPCDAG_MATPC;
      DiracOpMobius dop( *lattice, 0, 0, &cg_arg, CNV_FRM_NO );
      dop.RitzMat(Apsi, vtmp );
    }else if(lattice->Fclass()==F_CLASS_P4){
      cg_arg.RitzMatOper = MATPCDAG_MATPC;
      DiracOpP4 dop( *lattice, 0, 0, &cg_arg, CNV_FRM_NO );
      dop.RitzMat(Apsi, vtmp );
    }else{
      ERR.General(cname,fname,"Error: valid class type is dwf, mobius or p4\n");
    }

    Float alp =  Apsi->ReDotProductGlbSum( vtmp, f_size);

    Apsi->FTimesV1PlusV2(-alp, vtmp, Apsi, f_size);  

#if 0
    Float* ftmp1= (Float*)Apsi;
    //Float phase = atan2( ftmp1[1], ftmp1[0] );
    //Rcomplex rot_phase( cos( phase ), sin(- phase ) );
    
    for(int i=0;i<f_size;i+=2){
      Rcomplex C(ftmp1[i], ftmp1[i+1]);
      //C *= rot_phase;
      if( C.norm() > 1e-10)
	printf("%d  %e %e %e\n", i, C.real(),C.imag(),C.norm());
    }
#endif

    
    Float rnorm = sqrt(Apsi->NormSqGlbSum(f_size ));
    Float norm = sqrt(vtmp->NormSqGlbSum(f_size));
    VRB.Result(cname,fname, "Final True Residual  %e norm %e alpha %.16e\n",rnorm,norm,alp);


    if(residual) *residual = rnorm;
    if(eval) *eval = alp;

    sfree(cname,fname,"Apsi",Apsi);
  }


  void nev_check( int index, Float mass, Float* residual = 0, Float* eval = 0)
  {
    char* fname="nev_check(I,F,F*,F*)";
    Vector* vtmp = nev_load( index );
    nev_check( vtmp, mass, residual, eval );
  }

  
  void compress (Vector* sol, Float mass, int step_eig,  Float* eig1 )
  {
    char* fname = "compress(I)";
    VRB.Func(cname,fname);



    int ncompress = 1+ (neig-1)/step_eig;


    Float* eval = load_eval();

    sol->VecZero(f_size);

    Vector* vtmp = (Vector*) smalloc(cname,fname, "vtmp", f_size* sizeof(Float));
      if(vtmp==0)ERR.General(cname,fname,"vtmp could not malloced\n");
    Vector* Apsi = (Vector*) smalloc(cname,fname, "Apsi", f_size* sizeof(Float));
      if(Apsi==0)ERR.General(cname,fname,"Apsi could not malloced\n");
    // make a linear combination
    for(int iev=0;iev<neig; iev += step_eig) {
    
      Vector* evec = nev_load( iev );

      // test the eigen value of the R gamma_5 MatPc
      { 
	HermicianDWF_ee( vtmp, evec, mass, lattice, Apsi );
	Complex prod = evec -> CompDotProductGlbSum( vtmp, f_size );

	if(!UniqueID())
	printf("%d %e %e eval %e vs %e (%e)\n",
	       iev, prod.real(), prod.imag(),
	       prod.real()*prod.real(),
	       eval[iev],
	       2.0*(prod.real()*prod.real()- eval[iev])/
	       (prod.real()*prod.real()+ eval[iev]));

	eig1[ iev/step_eig ] = prod.real();
      }


    Complex z = 1.0;//evec->CompDotProductGlbSum( src, f_size );

    if(!UniqueID())
      printf("adding %d %d %g %g\n", iev, iev/step_eig, eig1[iev/step_eig], eval[iev]);
    // z /= eval[ iev ];

    sol -> CTimesV1PlusV2(z, evec, sol, f_size );
    }

    sfree(cname,fname,"vtmp",vtmp);
    sfree(cname,fname,"Apsi",Apsi);

  }

  void decompress( Vector** sol, Float mass,  int step_eig,  Float* eig1  )
  {
    char* fname="decompress(F,I,I)";
    VRB.Func(cname,fname);

    Vector* Apsi = (Vector*) smalloc(cname,fname, "Apsi", f_size* sizeof(Float));
      if(Apsi==0)ERR.General(cname,fname,"Apsi could not malloced\n");
    Vector* vtmp = (Vector*) smalloc(cname,fname, "vtmp", f_size* sizeof(Float));
      if(vtmp==0)ERR.General(cname,fname,"vtmp could not malloced\n");

  //-------------------------------------------------------------------
  // Now try to claim the lost informations
  {

    int n_lcon = (neig-1)/step_eig + 1;

    // sol[0] is the linear combination
    Float norm = sqrt(sol[0]->NormSqGlbSum(f_size ));
    sol[0] -> VecTimesEquFloat(1.0/norm, f_size);

    for(int i=0; i< n_lcon-1 ;++i) {
      HermicianDWF_ee( vtmp, sol[i], mass, lattice, Apsi );    
      sol[i+1] -> FTimesV1PlusV2( -1.0/eig1[i], vtmp, sol[i], f_size );
      Float norm = sqrt(sol[i+1]->NormSqGlbSum(f_size ));
      VRB.Result(cname,fname,"extracting %d %g %g\n",i,eig1[i], norm);
      sol[i+1] -> VecTimesEquFloat(1.0/norm, f_size);
    }


    VRB.Result(cname,fname,"claiming %d %g %g\n", 
	       n_lcon-1, eig1[n_lcon-1],1.0);
    nev_check( sol[n_lcon-1],  mass );

    for(int i=n_lcon-2; i>=0; --i){
      lanczos_GramSchm( (Float*)(sol[i]), (Float**)sol+i+1, 
			n_lcon-1-i, f_size,0);
      Float norm = sqrt(sol[i]->NormSqGlbSum(f_size ));
      VRB.Result(cname,fname,"claiming %d %g %g\n",i,eig1[i], norm);
      sol[i] -> VecTimesEquFloat(1.0/norm, f_size);
      nev_check( sol[i],  mass );
    }

  }

  //-------------------------------------------------------------------

    sfree(cname,fname,"Apsi",Apsi);
    sfree(cname,fname,"vtmp",vtmp);
    //    sfree(cname,fname,"sol",sol);
  }


};

CPS_END_NAMESPACE
#endif
