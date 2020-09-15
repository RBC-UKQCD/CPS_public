//Compute the low mode part of the W and V vectors. In the Lanczos class you can choose to store the vectors in single precision (despite the overall precision, which is fixed to double here)
//Set 'singleprec_evecs' if this has been done
template< typename mf_Policies>
void A2AvectorW<mf_Policies>::computeVWlow(A2AvectorV<mf_Policies> &V, Lattice &lat, BFM_Krylov::Lanczos_5d<double> &eig, bfm_evo<double> &dwf, bool singleprec_evecs){
  const char *fname = "computeVQlow(....)";
  typedef typename mf_Policies::ComplexType::value_type mf_Float;
  int gparity = GJP.Gparity();
  if(eig.dop.gparity != gparity){ ERR.General("A2AvectorW",fname,"Gparity must be disabled/enabled for *both* CPS and the eigenvectors"); }

  //Double precision temp fields
  CPSfermion4D<ComplexD> afield;  Vector* a = (Vector*)afield.ptr(); //breaks encapsulation, but I can sort this out later.
  CPSfermion5D<ComplexD> bfield;  Vector* b = (Vector*)bfield.ptr();

  int afield_fsize = afield.size()*sizeof(CPSfermion4D<ComplexD>::FieldSiteType)/sizeof(Float); //number of floats in field
  
  const int glb_ls = GJP.SnodeSites() * GJP.Snodes();
    
  Fermion_t tmp[2] = { dwf.allocFermion(), dwf.allocFermion() };
  Fermion_t tmp2[2] = { dwf.allocFermion(), dwf.allocFermion() };
  Fermion_t bq_tmp = dwf.allocCompactFermion(); 
  Fermion_t Mtmp = dwf.allocFermion();
    
  //The general method is described by page 60 of Daiqian's thesis
  for(int i = 0; i < nl; i++) {
    //Step 1) Compute V
    mf_Float* vi = (mf_Float*)V.getVl(i).ptr();

    //Copy bq[i][1] into bq_tmp
    const int len = 24 * eig.dop.node_cbvol * (1 + gparity) * eig.dop.cbLs;
    omp_set_num_threads(bfmarg::threads);		
    if(singleprec_evecs) // eig->bq is in single precision
#pragma omp parallel for  //Bet I could reduce the threading overheads by parallelizing this entire method
      for(int j = 0; j < len; j++) {
	((double*)bq_tmp)[j] = ((float*)(eig.bq[i][1]))[j];
      }
    else // eig.bq is in double precision
#pragma omp parallel 
      { dwf.axpy(bq_tmp, eig.bq[i][1], eig.bq[i][1], 0.);}

    //Do tmp = [ -(Mee)^-1 Meo bq_tmp, bg_tmp ]
#pragma omp parallel
    {		  
      dwf.Meo(bq_tmp,tmp[0],Even,0);	//tmp[0] = Meo bq_tmp 
      dwf.MooeeInv(tmp[0],tmp[1],0);  //tmp[1] = (Mee)^-1 Meo bq_tmp
      dwf.axpy(tmp[0],tmp[1],tmp[1],-2.0); //tmp[0] = -2*tmp[1] + tmp[1] = -tmp[1]
      dwf.axpy(tmp[1],bq_tmp,bq_tmp,0.0);  //tmp[1] = bq_tmp
    }
    //Get 4D part and poke into a
    dwf.cps_impexFermion((Float *)b,tmp,0);
    lat.Ffive2four(a,b,glb_ls-1,0,2); // a[4d] = b[5d walls]
    //Multiply by 1/lambda[i] and copy into v (with precision change if necessary)
    VecTimesEquFloat<mf_Float,Float>(vi, (Float*)a, 1.0 / eig.evals[i], afield_fsize);

    //Step 2) Compute Wl

    //Do tmp = [ -[Mee^-1]^dag [Meo]^dag Doo bq_tmp,  Doo bq_tmp ]    (Note that for the Moe^dag in Daiqian's thesis, the dagger also implies a transpose of the spatial indices, hence the Meo^dag in the code)
#pragma omp parallel
    {
      dwf.Mprec(bq_tmp,tmp[1],Mtmp,DaggerNo);  //tmp[1] = Doo bq_tmp
      dwf.Meo(tmp[1],Mtmp,Even,1); //Mtmp = Meo^dag Doo bq_tmp
      dwf.MooeeInv(Mtmp,tmp[1],1); //tmp[1] = [Mee^-1]^dag Meo^dag Doo bq_tmp
      dwf.axpy(tmp[0],tmp[1],tmp[1],-2.0); //tmp[0] = -tmp[1] = -[Mee^-1]^dag Meo^dag Doo bq_tmp
      dwf.Mprec(bq_tmp,tmp[1],Mtmp,DaggerNo); //tmp[1] = Doo bq_tmp
    }
    //Left-multiply by D-^dag for Mobius
    if(dwf.solver == HmCayleyTanh) {
#pragma omp parallel 
      {
	dwf.G5D_Dminus(tmp,tmp2,1);
	dwf.axpy(tmp, tmp2, tmp2, 0.0);
      }
    }
    //Get 4D part, poke onto a then copy into wl
    dwf.cps_impexFermion((Float *)b,tmp,0);
    lat.Ffive2four(a,b,0,glb_ls-1, 2);
    VecTimesEquFloat<mf_Float,Float>((mf_Float*)wl[i].ptr(), (Float*)a, 1.0, afield_fsize);
  }

  dwf.freeFermion(tmp[0]);
  dwf.freeFermion(tmp[1]);
  dwf.freeFermion(tmp2[0]);
  dwf.freeFermion(tmp2[1]);
  dwf.freeFermion(bq_tmp);
  dwf.freeFermion(Mtmp);
}



//Compute the high mode parts of V and W. 
//singleprec_evecs specifies whether the input eigenvectors are stored in single preciison
//You can optionally pass a single precision bfm instance, which if given will cause the underlying CG to be performed in mixed precision.
//WARNING: if using the mixed precision solve, the eigenvectors *MUST* be in single precision (there is a runtime check)
template< typename mf_Policies>
void A2AvectorW<mf_Policies>::computeVWhigh(A2AvectorV<mf_Policies> &V, BFM_Krylov::Lanczos_5d<double> &eig, bool singleprec_evecs, Lattice &lat, bfm_evo<double> &dwf_d, bfm_evo<float> *dwf_fp){
  const char *fname = "computeVWhigh(....)";
  typedef typename mf_Policies::ComplexType::value_type mf_Float;
  bool mixed_prec_cg = dwf_fp != NULL; 
  if(mixed_prec_cg && !singleprec_evecs){ ERR.General("A2AvectorW",fname,"If using mixed precision CG, input eigenvectors must be stored in single precision"); }

  VRB.Result("A2AvectorW", fname, "Start computing high modes.\n");
    
  //Generate the compact random sources for the high modes
  setWhRandom(args.rand_type);

  //Allocate temp *double precision* storage for fermions
  CPSfermion5D<ComplexD> afield,bfield;
  CPSfermion4D<ComplexD> v4dfield;

  int v4dfield_fsize = v4dfield.size()*sizeof(CPSfermion4D<ComplexD>::FieldSiteType)/sizeof(Float); //number of floats in field
  
  Vector *a = (Vector*)afield.ptr(), *b = (Vector*)bfield.ptr(), *v4d = (Vector*)v4dfield.ptr();

  const int glb_ls = GJP.SnodeSites() * GJP.Snodes();

  Fermion_t src[2] = {dwf_d.allocFermion(),dwf_d.allocFermion()}; 
  Fermion_t V_tmp[2] = {dwf_d.allocFermion(),dwf_d.allocFermion()};
  Fermion_t V_tmp2[2] = {dwf_d.allocFermion(),dwf_d.allocFermion()};
  Fermion_t tmp = dwf_d.allocFermion();

  //Details of this process can be found in Daiqian's thesis, page 60

  //Copy evals into multi1d
  multi1d<float> eval;
  multi1d<double> eval_d;
  if(mixed_prec_cg){
    eval.resize(eig.evals.size());
    for(int i = 0; i < eig.evals.size(); i++) eval[i] = eig.evals[i]; 
  }else{
    eval_d.resize(eig.evals.size());
    for(int i = 0; i < eig.evals.size(); i++) eval_d[i] = eig.evals[i]; 
  }

  for(int i=0; i<nh; i++){
    //Step 1) Get the diluted W vector to invert upon
    getDilutedSource(v4dfield, i);

    //Step 2) Solve V
    lat.Ffour2five(a, v4d, 0, glb_ls-1, 2); // poke the diluted 4D source onto a 5D source (I should write my own code to do this)
    dwf_d.cps_impexFermion((Float *)a,src,1);

    //Multiply src by D_ for Mobius 
    if(dwf_d.solver == HmCayleyTanh) {
#pragma omp parallel 
      {
	dwf_d.G5D_Dminus(src,V_tmp,0); //V_tmp = D_ src
	dwf_d.axpy(src, V_tmp, V_tmp, 0.0); //src = D_ src
	dwf_d.axpby(V_tmp, src, src, 0.0, 0.0); //V_tmp = 0
      }
    }
      
    //We can re-use previously computed solutions to speed up the calculation if rerunning for a second mass by using them as a guess
    //If no previously computed solutions this wastes a few flops, but not enough to care about
    //V vectors default to zero, so this is a zero guess if not reusing existing solutions
    VecTimesEquFloat<Float,mf_Float>((Float*)v4d, (mf_Float*)V.getVh(i).ptr(), 1.0, v4dfield_fsize); // v[i]->v4d to double precision
    lat.Ffour2five(a, v4d, 0, glb_ls-1, 2); // to 5d
    if(dwf_d.solver == HmCayleyTanh) {
      dwf_d.cps_impexFermion((Float*)a, V_tmp, 1); // to bfm
#pragma omp parallel 
      {
	dwf_d.G5D_Dminus(V_tmp,V_tmp2,0);
      }
    }else dwf_d.cps_impexFermion((Float*)a, V_tmp2, 1); // to bfm

    //Do the CG
    if(mixed_prec_cg){
      //Do the mixed precision deflated solve
      dwf_fp->residual = 1e-5; //todo: make this not hardcoded
      int max_cycle = 100;
      omp_set_num_threads(bfmarg::threads);
#pragma omp parallel
      { mixed_cg::threaded_cg_mixed_M(V_tmp2, src, dwf_d, *dwf_fp, max_cycle, CG, &eig.bq, &eval, nl); }
    }else{
      //Do a double precision deflated solve 
#pragma omp parallel
      {
	dwf_d.CGNE_M_high(V_tmp2, src, eig.bq, eval_d, nl, singleprec_evecs);
      }
    }

    //CPSify the solution, including 1/nhit for the hit average
    dwf_d.cps_impexFermion((Float *)b,V_tmp2,0);
    lat.Ffive2four(v4d, b, glb_ls-1, 0, 2);
    VecTimesEquFloat<mf_Float,Float>((mf_Float*)V.getVh(i).ptr(), (Float*)v4d, 1.0 / nhits, v4dfield_fsize);
  }
  dwf_d.freeFermion(src[0]); 
  dwf_d.freeFermion(src[1]);
  dwf_d.freeFermion(V_tmp[0]); 
  dwf_d.freeFermion(V_tmp[1]);
  dwf_d.freeFermion(V_tmp2[0]); 
  dwf_d.freeFermion(V_tmp2[1]);
  dwf_d.freeFermion(tmp); 
}
