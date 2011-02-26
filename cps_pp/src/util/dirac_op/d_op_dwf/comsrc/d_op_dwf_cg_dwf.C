//#if 0
#ifdef USE_CG_DWF
#include<config.h>
#include<util/dirac_op.h>
#include<util/dwf.h>
#include<util/pt.h>
#include<util/time_cps.h>
#include<util/lat_data.h>
#include<dwf-ssed.h>
#include<dwf-ssef.h>

//USING_NAMESPACE_CPS
#define USE_CPS_CG
#undef DIRAC_ONLY

CPS_START_NAMESPACE

static inline void CheckCoor(int ndim, int local[],const int *global){
	for (int i = 0;i<ndim;i++){
		if ( (global[i] / GJP.NodeSites(i) ) != GJP.NodeCoor(i))
		ERR.General("","CheckCoor()","global[%d]=%d,NodeCoor[%d]=%d",i,global[i],i,GJP.NodeCoor(i));
		local[i] = global[i] % GJP.NodeSites(i);	
		if (i == 4) 
		local[i] =  GJP.NodeSites(i) - 1 - local[i];	
	}	
}
static double gauge_reader(
					const void     *OuterGauge,
					void           *env,
					const int       *lattice_addr,
					int		dim,
					int             row,
					int             col,
					int             re_im){

	CPS_NAMESPACE::LatMatrix **mat = (CPS_NAMESPACE::LatMatrix **)OuterGauge;
	CPS_NAMESPACE::Lattice *lat = (CPS_NAMESPACE::Lattice *)env;
	int local[4];
	CheckCoor(4,local,lattice_addr);
#if 0
//	for(int i =0;i<4;i++)
//	printf("gauge_reader %d %d: %d %d %d\n",UniqueID(),i,lattice_addr[i],GJP.NodeCoor(i),local[i]);
	if ( row==col && re_im==0 ) return (double)(1.);
	else return 0.;
#else
	unsigned long index = lat->GsiteOffset(local)/4;
	index = re_im+2*(col+3*(row+3*index));
	double *pointer = mat[dim]->Field();
	return pointer[index];
#endif

}

static MIT_ssed_DWF_gauge_reader gauge_reader_p = &gauge_reader;

double d_fermion_reader(
					const void     *OuterFermion,
					void           *env,
					const int       lattice_addr[5],
					int             color,
					int             dirac,
					int             re_im){

	Lattice *lat = (Lattice *)env;
	int local[5];
	CheckCoor(5,local,lattice_addr);
//	for(int i =0;i<5;i++)

//	printf("fermion_reader %d %d: %d %d\n",UniqueID(),i,lattice_addr[i],local[i]);
	unsigned long index = lat->FsiteOffsetChkb(local);
	index = re_im+2*(color+3*(dirac+4*index));
	double *pointer = (double *)OuterFermion;
#if 0
	if ( fabs( *(pointer+index) ) >1e-8){
	printf("fermion_reader: (%d %d %d %d %d):(%d) (%d %d %d) %g\n",
	lattice_addr[0], lattice_addr[1], lattice_addr[2], lattice_addr[3],
	lattice_addr[4], lat->FsiteOffsetChkb(local), color,dirac,re_im, *(pointer+index));
	}
#endif
	double val =  *(pointer+index);
	if ( dirac >=2 ) val = -val;
	return val ;
}
static MIT_ssed_DWF_fermion_reader fermion_reader_p=&d_fermion_reader;

typedef struct d_fermion_writer_arg{
 Lattice *lat;
 double coeff;
} fermion_writer_arg;

static void d_fermion_writer(
					void     *OuterFermion,
					void           *env,
					const int       lattice_addr[5],
					int             color,
					int             dirac,
					int             re_im,
					double			value){

	fermion_writer_arg *arg = (fermion_writer_arg *)env;
	Lattice *lat = arg->lat;
//	printf("d_fermion_writer(%p %p %d %d %d %g)\n",
//	OuterFermion,env,color,dirac,re_im, value);
	int local[5];
	CheckCoor(5,local,lattice_addr);
	unsigned long index = lat->FsiteOffsetChkb(local);
	index = re_im+2*(color+3*(dirac+4*index));
	double *pointer = (double *)OuterFermion;
#if 0
	if ( fabs(value)  >1e-8){
	printf("fermion_writer: (%d %d %d %d %d):(%d) (%d %d %d) %g\n",
	lattice_addr[0], lattice_addr[1], lattice_addr[2], lattice_addr[3],
	lattice_addr[4], lat->FsiteOffsetChkb(local), color,dirac,re_im, value);
	}
#endif
	if ( dirac >=2 ) value = -value;
	*(pointer+index) = arg->coeff*value;
}
static MIT_ssed_DWF_fermion_writer fermion_writer_p = & d_fermion_writer;

static   LatMatrix *Plus[4];
static   LatMatrix *Minus[4];
static int lat_allocated=0;


//------------------------------------------------------------------
// int MatInv(Vector *out, Vector *in, 
//            Float *true_res, PreserveType prs_in);
// The inverse of the unconditioned Dirac Operator 
// using Conjugate gradient.
// If true_res !=0 the value of the true residual is returned
// in true_res.
// *true_res = |src - MatPcDagMatPc * sol| / |src|
// prs_in is used to specify if the source
// in should be preserved or not. If not the memory usage
// is less by half the size of a fermion vector.
// The function returns the total number of CG iterations.
//------------------------------------------------------------------
int DiracOpDwf::MatInv(Vector *out, 
		       Vector *in, 
		       Float *true_res,
		       PreserveType prs_in) {
  char *fname = "MatInv(V*,V*,F*)";
  VRB.Func(cname,fname);
  VRB.Result(cname,fname,"Using cg-dwf");

  //----------------------------------------------------------------
  // Initialize kappa and ls. This has already been done by the Fdwf
  // and DiracOpDwf constructors but is done again in case the
  // user has modified the GJP.DwfA5Inv(), GJP.DwfHeight() or
  // GJP.SnodeSites() while in the scope of the DiracOpDwf object.
  //----------------------------------------------------------------
  Dwf *dwf_arg = (Dwf *) dwf_lib_arg;
  dwf_arg->ls = GJP.SnodeSites();
  dwf_arg->dwf_kappa = 
    1.0 / ( 2 * (4 + GJP.DwfA5Inv() - GJP.DwfHeight()) );

  //----------------------------------------------------------------
  // Implement routine
  //----------------------------------------------------------------
  Vector *temp2;
  int temp_size = GJP.VolNodeSites() * lat.FsiteSize() / 2;


  Vector *temp = (Vector *) smalloc(cname,fname, "temp",temp_size * sizeof(Float));

//  if(prs_in == PRESERVE_YES){
    temp2 = (Vector *) smalloc(cname,fname,"temp2",2*temp_size * sizeof(Float));
//  }

  Vector *Ax = (Vector *) smalloc(cname,fname, "temp",2*temp_size * sizeof(Float));
  Vector *sol = (Vector *) smalloc(cname,fname, "temp",2*temp_size * sizeof(Float));

    // points to the even part of fermion source 
  Vector *even_in = (Vector *) ( (IFloat *) in + temp_size );

  // points to the even part of fermion solution
  Vector *even_out = (Vector *) ( (IFloat *) out + temp_size );

  Dslash(temp, even_in, CHKB_EVEN, DAG_NO);
  fTimesV1PlusV2((IFloat *)temp, (IFloat) dwf_arg->dwf_kappa, (IFloat *)temp,
    (IFloat *)in, temp_size);


  LatMatrix Unit;
  Unit = 1.;

  {
    ParTransGauge pt_g(lat);
    Matrix *min[4], *mout_p[4], *mout_m[4];
    int dirs[]={0,2,4,6,1,3,5,7};
    for(int i =0;i<4;i++){
	  if (!lat_allocated){
      	Plus[i] = new LatMatrix;
	VRB.Result(cname,fname,"Plus[%d]=%p\n",i,Plus[i]);
      	Minus[i] = new LatMatrix;
	VRB.Result(cname,fname,"Minus[%d]=%p\n",i,Minus[i]);
	  }
      min[i]=min[i+4]=Unit.Mat();
      mout_p[i]=Plus[i]->Mat();
      mout_m[i]=Minus[i]->Mat();
    }
    pt_g.run(4,mout_p,min,dirs);
//    pt_g.shift_link(mout_p,dirs,4);
//    pt_g.shift_field(mout_p,dirs+4,4,1,mout_m);
    pt_g.run(4,mout_m,min,dirs+4);
    for(int i =0;i<4;i++){
	Matrix temp_m;
	Matrix *mat_p = mout_m[i];
	for(int j=0;j<GJP.VolNodeSites();j++){
	    temp_m.Dagger(*mat_p);
	    *mat_p = temp_m;
	    mat_p++;
	}
    }
	lat_allocated=1;
  }


  // save source
//  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)temp2, (IFloat *)in, 
		2*temp_size * sizeof(IFloat) / sizeof(char));
//  }


  int iter;
//    	const int lattice[5]={1,1,1,1,1};
    int lattice[5];
    for(int i=0;i<5;i++) lattice[i]=GJP.Sites(i);

#define DOUBLE(n) MIT_ssed_DWF_##n
#define SINGLE(n) MIT_ssef_DWF_##n

   	int exit_code = SINGLE(init)(lattice,NULL,NULL);
	if (exit_code) exit(-4);

	SINGLE(Gauge) *g_s = SINGLE(load_gauge)(Plus,Minus,&lat,gauge_reader);
	SINGLE(Fermion) *x_s = SINGLE(allocate_fermion)();
	SINGLE(Fermion) *x0_s = NULL;


	double residual;
	double true_res_re,true_res_im;
	double M_0= -2.*(5. -GJP.DwfHeight());
	double m_f=dirac_arg->mass;
	double true_eps_re,true_eps_im;
	true_eps_re=in->NormSqGlbSum(2*temp_size);
	double stp_cnd = 0.0001*dirac_arg->stop_rsd*dirac_arg->stop_rsd *true_eps_re;
	VRB.Result(cname,fname,"M_0=%g m_f=%g,stp_cnd=%e",M_0,m_f,stp_cnd);
	int offset = GJP.VolNodeSites()*lat.FsiteSize()/ (2*6);

  switch (dirac_arg->Inverter) {
  case CG:
//    iter = InvCg(out,in,true_res);

	fermion_writer_arg save_arg;


	if (dirac_arg->max_num_iter>0){

#ifdef DIRAC_ONLY
		VRB.Func(fname,"dirac op");
		DOUBLE(Dirac_Operator)(x,g,M_0,m_f,rhs);
#else
		VRB.Func(fname,"cg_solver");
		int d_iter;
		iter = d_iter=0;
		Float cg_time = 0.;
		const int SINGLE_CG_CUTOFF=3000; //arbitrary number to cut of single CG
		int Nsingle = dirac_arg->max_num_iter/SINGLE_CG_CUTOFF+3;
		cg_time -=dclock();
		for(int i=0;i<Nsingle && true_eps_re >stp_cnd;i++){
			SINGLE(Fermion) *rhs_s = SINGLE(load_fermion)(in,&lat,d_fermion_reader);
			if(i==0) x0_s = SINGLE(load_fermion)(out,&lat,d_fermion_reader);
			else x0_s = SINGLE(allocate_fermion)();
			SINGLE(cg_solver)(x_s, &residual,&d_iter,g_s,M_0,m_f,x0_s,rhs_s,
				stp_cnd,0,SINGLE_CG_CUTOFF);
			iter += d_iter;
	
			save_arg.lat = &lat;
			save_arg.coeff = 1.;
			if (i==0) 
				SINGLE(save_fermion)(out,&save_arg,d_fermion_writer,x_s);
			else{
				SINGLE(save_fermion)(sol,&save_arg,d_fermion_writer,x_s);
				out->FTimesV1PlusV2(1.,sol,out,2*temp_size);
			}
			SINGLE(delete_fermion)(x0_s);
			SINGLE(delete_fermion)(rhs_s);

			Dslash(Ax,out+offset,CHKB_EVEN,DAG_NO);
			Dslash(Ax+offset,out,CHKB_ODD,DAG_NO);
			Ax->FTimesV1PlusV2(M_0,out,Ax,2*temp_size);
			in->FTimesV1PlusV2(-1,Ax,temp2,2*temp_size);
			true_eps_re=in->NormSqGlbSum(2*temp_size);
			VRB.Result(cname,fname,"Norm(in)=%g\n" ,true_eps_re);

			
		} // for(int i=0;i<Nsingle && true_eps_re >stp_cnd;i++)
		cg_time += dclock();

	// lifted from Chroma flop counting by CJ
		{
			int Nc=3,Ns=4;
			unsigned long Ls = lattice[4];
			unsigned long long Ndiag  = (4*Ls+2)*Nc*Ns; 
			unsigned long long NdiagInv = (10*Ls-8)*Nc*Ns;
			unsigned long long Neo    = Ls*(1320+24);
			unsigned long long N_mpsi = 2*Ndiag + 2*Neo + Ls*24;
			unsigned long long Nflops_c = (24*Ls + 2*N_mpsi) + (48*Ls);  /* constant term */
			unsigned long long Nflops_s = (2*N_mpsi + Ls*(2*48+2*24));   /* slope */
			unsigned long long Nflops_per_cbsite = Nflops_c + (iter )*Nflops_s;
			double  Nflops_total = Nflops_per_cbsite*GJP.VolNodeSites()/2;
			print_flops("cg_dwf","SINGLE(cg_solver)",Nflops_total,cg_time);
		}

#endif //DIRAC_ONLY
	}


	VRB.Func(fname,"save_fermion");
#if 0
	save_arg.lat = &lat;
	save_arg.coeff = -2.*(5.0-GJP.DwfHeight());
	SINGLE(save_fermion)(out,&save_arg,d_fermion_writer,x_s);
#endif
//	out->FTimesV1PlusV2((M_0-1.),out,out,2*temp_size);
	out->VecTimesEquFloat(M_0,2*temp_size);
	QMP_barrier();
	VRB.FuncEnd(fname,"save_fermion");
#ifdef USE_CPS_CG
	MatPcDag(in, temp);
	iter += InvCg(out,in,true_res);
#endif

    break;
  default:
    ERR.General(cname,fname,"InverterType %d not implemented\n",
		dirac_arg->Inverter);
  }


	SINGLE(delete_gauge)(g_s);
	SINGLE(delete_fermion)(x_s);
	SINGLE(fini)();
	VRB.FuncEnd(fname,"DWF_fini");

#undef DOUBLE
#undef SINGLE

  // restore source
  if(prs_in == PRESERVE_YES){
    moveMem((IFloat *)in, (IFloat *)temp2, 
		2*temp_size * sizeof(IFloat) / sizeof(char));
  }

#ifdef USE_CPS_CG
  Dslash(temp, out, CHKB_ODD, DAG_NO);

  fTimesV1PlusV2((IFloat *)even_out, (IFloat) dwf_arg->dwf_kappa, (IFloat *)temp,
    (IFloat *) even_in, temp_size);
#endif


  sfree(cname,fname,"temp",temp);
  sfree(cname, fname, "Ax", Ax);
  sfree(cname, fname, "sol", sol);

//  if(prs_in == PRESERVE_YES){
    sfree(cname, fname, "temp2", temp2);
//  }

  return iter;
}
CPS_END_NAMESPACE
#endif
