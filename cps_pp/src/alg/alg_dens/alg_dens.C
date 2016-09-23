#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief Methods of the AlgDens class.
  
  $Id: alg_dens.C,v 1.8 2008/02/14 20:45:44 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/02/14 20:45:44 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_dens/alg_dens.C,v 1.8 2008/02/14 20:45:44 chulwoo Exp $
//  $Id: alg_dens.C,v 1.8 2008/02/14 20:45:44 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: alg_dens.C,v $
//  $Revision: 1.8 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_dens/alg_dens.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// alg_dens.C
//
// AlgDens is derived from Alg and is relevant to the 
// stochastic measurement of derivatives of the partition function
// with respect to the chemical potential using the
// Conjugate Gradient algorithm. The type of fermion is
// determined by the argument to the constructor.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <alg/alg_dens.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE

#define POINT
#undef POINT
#define Z2
#undef Z2
//------------------------------------------------------------------
/*!
  \param latt The lattice on which to compute the condensate.
  \param c_arg The common argument structure for all algorithms.
  \param arg The algorithm parameters.
 */
//------------------------------------------------------------------
AlgDens::AlgDens(Lattice& latt, 
	       CommonArg *c_arg,
	       DensArg *arg) : 
	       Alg(latt, c_arg) 
{
  cname = "AlgDens";
  char *fname = "AlgDens(L&,CommonArg*,DensArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_dens_arg = arg;


  // Set the node size of the full (non-checkerboarded) fermion field
  //----------------------------------------------------------------
  f_size = GJP.VolNodeSites() * latt.FsiteSize();


  // Allocate memory for the source.
  //----------------------------------------------------------------
  src = (Vector *) smalloc(f_size * sizeof(Float));
  if(src == 0)
    ERR.Pointer(cname,fname, "src");
  VRB.Smalloc(cname,fname, "src", src, f_size * sizeof(Float));

  srcM = (Vector *) smalloc(f_size * sizeof(Float));
  if(srcM == 0)
    ERR.Pointer(cname,fname, "srcM");
  VRB.Smalloc(cname,fname, "srcM", srcM, f_size * sizeof(Float));


  // Allocate memory for the solutions
  //----------------------------------------------------------------
  save = (Vector *) smalloc(f_size * (alg_dens_arg->max_save) * sizeof(Float));
  if(save == 0)
    ERR.Pointer(cname,fname, "save");
  VRB.Smalloc(cname,fname, "save", save, f_size * (alg_dens_arg->max_save) * sizeof(Float));
  
  sol = (Vector *) smalloc(f_size * sizeof(Float));
  if(sol == 0)
    ERR.Pointer(cname,fname, "sol");
  VRB.Smalloc(cname,fname, "sol", sol, f_size * sizeof(Float));

  solM = (Vector *) smalloc(f_size * sizeof(Float));
  if(solM == 0)
    ERR.Pointer(cname,fname, "solM");
  VRB.Smalloc(cname,fname, "solM", solM, f_size * sizeof(Float));

  //map_table = (int *) smalloc( (alg_dens_arg->max_save) * sizeof(int));
  //if(map_table == 0)
  //  ERR.Pointer(cname,fname, "map_table");
  //VRB.Smalloc(cname,fname, "map_table", map_table, (alg_dens_arg->max_save) * sizeof(int));

  //save_table = (int *) smalloc( (alg_dens_arg->n_obs) * sizeof(int));
  //if(save_table == 0)
  //  ERR.Pointer(cname,fname, "save_table");
  //VRB.Smalloc(cname,fname, "save_table", save_table, (alg_dens_arg->n_obs) * sizeof(int));

  //load_table = (int *) smalloc( (alg_dens_arg->n_obs) * sizeof(int));
  //if(load_table == 0)
  //  ERR.Pointer(cname,fname, "load_table");
  //VRB.Smalloc(cname,fname, "load_table", load_table, (alg_dens_arg->n_obs) * sizeof(int));

  //refresh_table = (int *) smalloc( (alg_dens_arg->n_obs) * sizeof(int));
  //if(refresh_table == 0)
  //  ERR.Pointer(cname,fname, "refresh_table");
  //VRB.Smalloc(cname,fname, "refresh_table", refresh_table, (alg_dens_arg->n_obs) * sizeof(int));

  //a_coor_table = (int *) smalloc( (alg_dens_arg->n_obs) * sizeof(int));
  //if(a_coor_table == 0)
  //  ERR.Pointer(cname,fname, "a_coor_table");
  //VRB.Smalloc(cname,fname, "a_coor_table", a_coor_table, (alg_dens_arg->n_obs) * sizeof(int));

  //b_coor_table = (int *) smalloc( (alg_dens_arg->n_obs) * sizeof(int));
  //if(b_coor_table == 0)
  //  ERR.Pointer(cname,fname, "b_coor_table");
  //VRB.Smalloc(cname,fname, "b_coor_table", b_coor_table, (alg_dens_arg->n_obs) * sizeof(int));

}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgDens::~AlgDens() {
  char *fname = "~AlgDens()";
  VRB.Func(cname,fname);

  // Free memory
  //----------------------------------------------------------------
  VRB.Sfree(cname,fname, "save", save);
  sfree(save);
  VRB.Sfree(cname,fname, "sol", sol);
  sfree(sol);
  VRB.Sfree(cname,fname, "solM", solM);
  sfree(solM);
  VRB.Sfree(cname,fname, "src", src);
  sfree(src);
  VRB.Sfree(cname,fname, "srcM", srcM);
  sfree(srcM);
  //VRB.Sfree(cname,fname, "map_table", map_table);
  //sfree(map_table);
  //VRB.Sfree(cname,fname, "save_table", save_table);
  //sfree(save_table);
  //VRB.Sfree(cname,fname, "load_table", load_table);
  //sfree(load_table);
  //VRB.Sfree(cname,fname, "refresh_table", refresh_table);
  //sfree(refresh_table);
  //VRB.Sfree(cname,fname, "a_coor_table", a_coor_table);
  //sfree(a_coor_table);
  //VRB.Sfree(cname,fname, "b_coor_table", b_coor_table);
  //sfree(b_coor_table);

}


//------------------------------------------------------------------
//! Performs the computation.
/*!
  \post The results are written to the file specified in the common_arg
  structure.
*/
//------------------------------------------------------------------
void AlgDens::run()
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif
  int iter=0;
  int ls;
  int ls_glb;
  int obsID;
  int order;
  Complex TrObs(0.,0.);
  Float norm;
  Float true_res;
  Float true_res2;
  FILE *fp=NULL;
  DensArg *dens_arg;
  CgArg cg_arg_struct;
  CgArg *cg_arg = &cg_arg_struct;
  char *fname = "run()";
  VRB.Func(cname,fname);
  //VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);


  // Set the Lattice pointer dens_arg and cg_arg
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();
  dens_arg = alg_dens_arg;
  cg_arg->max_num_iter = dens_arg->max_num_iter;
  cg_arg->stop_rsd = dens_arg->stop_rsd;
  
  // Make a Float pointer to sol
  //----------------------------------------------------------------
  Float *f_sol = (Float *) sol;

  //----------------------------------------------------------------
  // Initialize source and solution (initial guess)
  // and calculate PsiBarPsi for:
  //----------------------------------------------------------------


  //----------------------------------------------------------------
  // Domain Wall fermions
  //----------------------------------------------------------------
  if(lat.Fclass() == F_CLASS_DWF){

      ERR.General(cname, fname, "not implemented!\n");

  }

  //----------------------------------------------------------------
  // Wilson or Clover fermions
  //----------------------------------------------------------------
  else if (   (lat.Fclass() == F_CLASS_WILSON) 
	   || (lat.Fclass() == F_CLASS_CLOVER)   ) {  

      ERR.General(cname, fname, "not implemented!\n");

  }

  //----------------------------------------------------------------
  // Staggered fermions
  //----------------------------------------------------------------
  else if ( (lat.Fclass() == F_CLASS_STAG)
	   || (lat.Fclass() == F_CLASS_P4)
	   || (lat.Fclass() == F_CLASS_ASQTAD)   ) {  


    /* open file
       --------- */

    if(common_arg->results != 0){
      if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	ERR.FileA(cname,fname, (char *)common_arg->results);
      }
    }

    make_tables();

    /***** loop over sources *****/
    for(int source=0; source < dens_arg->n_src; source++){

      /* set source
	 ---------- */

      #ifdef POINT
//      bzero((char *)src, f_size*sizeof(Float));
      memset((char *)src, 0,f_size*sizeof(Float));
      if(CoorX()==0 && CoorY()==0 && CoorZ()==0 && CoorT()==0)
	{
	  *((IFloat *) &src[0]) = 1.0;
	  printf("Setting point source at origin\n");
	}
      for(int zz = 0; zz < GJP.VolNodeSites(); zz++)
	for(int yy = 0; yy < 6; yy++)
	  if(*(((IFloat *)&src[zz])+yy) > 1e-15)
	    printf("zz = %d, yy = %d, *(src+zz) = %0.16e\n", zz, yy, *(((IFloat *)&src[zz])+yy));
      #else
      lat.RandGaussVector(src, 0.5);
      #endif
      
      #ifdef Z2
      IFloat * tmp1;
      for(int zz = 0; zz < GJP.VolNodeSites(); zz++)
	for(int yy = 0; yy < 6; yy+=2)
	  {
	    tmp1 = (IFloat *)(src + zz);
	    if(*(tmp1+yy) > 0)
	      *(tmp1+yy) = 1.0;
	    else
	      *(tmp1+yy) = -1.0;
	    *(tmp1 +yy+ 1) = 0.0;
	  }
      //printf("Setting Z_2 source.\n");
      //int check_site = 29;
      //printf("Site %d =\n",check_site);
      //tmp1 = (IFloat *)(src + check_site);
      //for(int yy = 0; yy < 6; yy++)
      //  printf("%0.16e  ",*(tmp1+yy));
      //printf("\n");
      #endif

      clear_vector();

      #ifdef POINT
      norm = 1.0;
      #else
      norm = GJP.VolSites() * ( lat.FsiteSize() / 2 );
      #endif

      // Initialize the cg_arg mass, with the first mass we
      // want to compute for:
      switch( dens_arg->pattern_kind ) {
      case ARRAY: 
	cg_arg->mass = dens_arg->mass[0]; 
	break;
      case LIN:   
	cg_arg->mass = dens_arg->mass_start; 
	break;
      case LOG:   
	cg_arg->mass = dens_arg->mass_start; 
	break;
      default: 
	ERR.General(cname, fname,
		    "dens_arg->kind = %d is unrecognized\n", 
		    dens_arg->pattern_kind);
	break;
      }
      
      /***** Loop over masses *****/
      for(int m=0; m<dens_arg->n_masses; m++){
	
	save_vector(src,0);
      
	/***** loop over operators *****/
	for (int iobs=0; iobs<dens_arg->n_obs; iobs++){

	  obsID=dens_arg->obs[iobs];
	  order=(b_coor_table[iobs])%((dens_arg->max_deri)+1);
	  VRB.Flow(cname,fname,"%d: calculate obsID %d (%d,%d) order=%d \n",
		   iobs,obsID,a_coor_table[iobs],b_coor_table[iobs],order);
	  
	  if(refresh_table[iobs]){
	    
	    /* initialize solution 
	       ------------------- */

	    for(int i=0; i< f_size/2; i++){
	      f_sol[2*i] = 1.0;     // real part
	      f_sol[2*i+1] = 0.0;   // imaginary part
	    }
	    
	    /* load vector
	       ----------- */

	    load_vector(srcM,load_table[iobs]);
	    
	    /* do inversion
	       ------------ */

	    iter = lat.FmatInv(sol, srcM, cg_arg, &true_res, CNV_FRM_YES);
	    
	  }

	  /* calculate Tr
	     ------------ */
	  if(order>0){
	    lat.FdMdmu(solM, sol, cg_arg, CNV_FRM_YES, order);
	    TrObs = solM->CompDotProductGlbSum(src, f_size) / norm;
	  }
	  else{
	    TrObs = sol->CompDotProductGlbSum(src, f_size) / norm;
	    //if(obsID==1) TrObs *= 2.0 ;
	    if(b_coor_table[iobs]==0){
	      for(int i=0; i<a_coor_table[iobs]; i++){
		TrObs *= 2.0;
	      }
	    }
	  }
	  
	  /* save vector
	     ----------- */

	  if(save_table[iobs]){
	    if(order>0){
	      save_vector(solM,obsID);
	    }
	    else{
	      save_vector(sol,obsID);
	    }
	  }
	  
	  /* print out results
	     ----------------- */

	  if(refresh_table[iobs]){
	    // Print out obsID, mass, Tr, number of iterations and true residual
	    Fprintf(fp, "%12d %3d %12.6e %23.16e %23.16e %4d %22.16e \n",
		    obsID, source, IFloat(cg_arg->mass), 
		    TrObs.real(), 
		    TrObs.imag(),  
		    iter, IFloat(true_res));
	  }
	  else{
	    // Print out obsID, mass, Tr
	    Fprintf(fp, "%12d %3d %12.6e %23.16e %23.16e \n",
		    obsID, source, IFloat(cg_arg->mass), 
		    TrObs.real(), 
		    TrObs.imag());
	  }
	}
	
	// If there is another mass loop iteration ahead, we should
	// set the cg_arg->mass to it's next desired value
	if( m < dens_arg->n_masses - 1 ) {
	  switch( dens_arg->pattern_kind ) {
	  case ARRAY: 
	    cg_arg->mass = dens_arg->mass[m+1]; 
	    break;
	  case LIN:   
	    cg_arg->mass += dens_arg->mass_step; 
	    break;
	  case LOG:   
	    cg_arg->mass *= dens_arg->mass_step; 
	    break;
	  }
	}
      }
    }
    Fclose(fp);
  }  

  //----------------------------------------------------------------
  // Unknown fermion type
  //----------------------------------------------------------------
  else {
    ERR.General(cname,fname,"Unknown class type %d\n",int(lat.Fclass()));
  }

}


void AlgDens::make_tables()
{
  int obsID;
  int rest, old_rest;
  int a,b,ak,bk,aprev,bprev;
  int brefmin,brefmax,anext;
  int bnextmin,bnextmax;
  int mult;
  DensArg *dens_arg;
  char *fname = "make_tables()";
  VRB.Func(cname,fname);

  dens_arg = alg_dens_arg;

  for(int zz = 0; zz<dens_arg->n_obs; zz++){
    obsID=dens_arg->obs[zz];
    rest=obsID; a=0; b=0; mult=1;
    while(rest>0){
      a++;
      mult*=(dens_arg->max_deri)+1;
      old_rest=rest;
      rest-=mult;
    }
    b=old_rest-1;
    // *(a_coor_table+zz)=a;
    // *(b_coor_table+zz)=b;
    a_coor_table[zz]=a;
    b_coor_table[zz]=b;
    VRB.Flow(cname,fname,"coor_table: %d %d %d\n",zz,a,b);
  }
    
  for(int zz = 0; zz<dens_arg->n_obs; zz++){
    obsID=dens_arg->obs[zz];

    // a = *(a_coor_table+zz);
    // b = *(b_coor_table+zz);
    a=a_coor_table[zz];
    b=b_coor_table[zz];

    // *(load_table+zz)=-1;
    // *(save_table+zz)=0;
    // *(refresh_table+zz)=1;
    load_table[zz]=-1;
    save_table[zz]=0;
    refresh_table[zz]=1;

    if(a==1){
      // *(load_table+zz)=0;
      load_table[zz]=0;
    }

    aprev=a-1;
    mult=(dens_arg->max_deri)+1;
    bprev=b/mult;

    brefmin=b-b%((dens_arg->max_deri)+1);
    brefmax=brefmin+(dens_arg->max_deri);

    anext=a+1;
    bnextmin=b*((dens_arg->max_deri)+1);
    bnextmax=bnextmin+(dens_arg->max_deri);

    for(int kk = 0; kk<dens_arg->n_obs; kk++){
      //ak = *(a_coor_table+kk);
      //bk = *(b_coor_table+kk);
      ak = a_coor_table[kk];
      bk = b_coor_table[kk];
      if((ak==aprev)&&(bk==bprev)){
	// *(load_table+zz)=dens_arg->obs[kk];
	load_table[zz]=dens_arg->obs[kk];
      }
      if(ak==a){
	if((bk>=brefmin)&&(bk<=brefmax)){
	  if(bk<b){
	    // *(refresh_table+zz)=0;
	    refresh_table[zz]=0;
	  }
	}
      }
      if(ak==anext){
	if((bk>=bnextmin)&&(bk<=bnextmax)){
	  // *(save_table+zz)=1;
	  save_table[zz]=1;
	}
      }
    }
    
    VRB.Flow(cname,fname,"                                         \n");
    VRB.Flow(cname,fname,"   load_table: %d %d -> (%d, %d)-(%d, %d)\n",zz,load_table[zz],aprev,bprev,a,b);
    VRB.Flow(cname,fname,"                                         \n");
    VRB.Flow(cname,fname,"                             (%d, %d)    \n",anext,bnextmax);
    VRB.Flow(cname,fname,"                            /            \n");
    VRB.Flow(cname,fname,"   save_table: %d %d -> (%d, %d)         \n",zz,save_table[zz],a,b);
    VRB.Flow(cname,fname,"                            \\           \n");
    VRB.Flow(cname,fname,"                             (%d, %d)    \n",anext,bnextmin);
    VRB.Flow(cname,fname,"                                         \n");
    VRB.Flow(cname,fname,"                             (%d, %d)    \n",a,brefmax);
    VRB.Flow(cname,fname,"                            /            \n");
    VRB.Flow(cname,fname,"refresh_table: %d %d -> (%d, %d)         \n",zz,refresh_table[zz],aprev,bprev);
    VRB.Flow(cname,fname,"                            \\           \n");
    VRB.Flow(cname,fname,"                             (%d, %d)    \n",a,brefmin);
    

    //if(*(load_table+zz)==-1){
    if(load_table[zz]==-1){
      ERR.General(cname,fname,"Operator can not be computed, broken chain.\n");
    }
  }
}

void AlgDens::clear_vector()
{
  int offset;
  DensArg *dens_arg;
  char *fname = "clear_vector()";
  VRB.Func(cname,fname);

  dens_arg = alg_dens_arg;

  for(int zz = 0; zz<dens_arg->max_save; zz++){
    offset=zz*GJP.VolNodeSites();
//    bzero((char *)(save+offset), f_size*sizeof(Float));
    memset((char *)(save+offset), 0, f_size*sizeof(Float));
    // *(map_table+zz)=-1;
    map_table[zz]=-1;
  }
    
}

void AlgDens::save_vector(Vector *v, int obsID)
{
  int entry;
  int offset;
  IFloat *tmp1;
  IFloat *tmp2;
  DensArg *dens_arg;
  char *fname = "save_vector()";
  VRB.Func(cname,fname);

  dens_arg = alg_dens_arg;

  entry=-1;
  for(int zz = 0; zz<dens_arg->max_save; zz++){
    VRB.Flow(cname,fname,"(save) map_table: %d %d \n", zz, map_table[zz]);
    //if(*(map_table+zz)<0){
    if(map_table[zz]<0){
      entry=zz;
      zz=dens_arg->max_save;
    }
  }
  if(entry<0){
    ERR.General(cname,fname,"Vector can not bee saved. Field ran out of memory.\n");
  }
  offset=entry*GJP.VolNodeSites();
  // *(map_table+entry)=obsID;
  map_table[entry]=obsID;
  for(int zz = 0; zz < GJP.VolNodeSites(); zz++){
    tmp1=(IFloat *)(save + zz + offset);
    tmp2=(IFloat *)(v+zz);
    *(tmp1+0)=*(tmp2+0);
    *(tmp1+1)=*(tmp2+1);
    *(tmp1+2)=*(tmp2+2);
    *(tmp1+3)=*(tmp2+3);
    *(tmp1+4)=*(tmp2+4);
    *(tmp1+5)=*(tmp2+5);
  }
      
}

void AlgDens::load_vector(Vector *v, int obsID)
{
  int entry;
  int offset;
  IFloat *tmp1;
  IFloat *tmp2;
  DensArg *dens_arg;
  char *fname = "load_vector()";
  VRB.Func(cname,fname);

  dens_arg = alg_dens_arg;

  entry=-1;
  for(int zz = 0; zz<dens_arg->max_save; zz++){
    VRB.Flow(cname,fname,"(load) map_table: %d %d \n", zz, map_table[zz]);
    //if(*(map_table+zz)==obsID){
    if(map_table[zz]==obsID){
      entry=zz;
      zz=dens_arg->max_save;
    }
  }
  if(entry<0){
    ERR.General(cname,fname,"Vector can not bee loaded. Vector not found.\n");
  }
  offset=entry*GJP.VolNodeSites();
  // *(map_table+entry)=-1;
  map_table[entry]=-1;
  for(int zz = 0; zz < GJP.VolNodeSites(); zz++){
    tmp1=(IFloat *)(save + zz + offset);
    tmp2=(IFloat *)(v+zz);
    *(tmp2+0)=*(tmp1+0);
    *(tmp2+1)=*(tmp1+1);
    *(tmp2+2)=*(tmp1+2);
    *(tmp2+3)=*(tmp1+3);
    *(tmp2+4)=*(tmp1+4);
    *(tmp2+5)=*(tmp1+5);
  }
      
}

CPS_END_NAMESPACE
