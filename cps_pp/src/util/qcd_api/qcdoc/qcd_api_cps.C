#include <config.h>
#include <sys/time.h>
#include <qmp.h>
#include <math.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/smalloc.h>
#include <util/error.h>
#include <util/vector.h>
#include <util/time.h>
#include <alg/common_arg.h>
#include <util/dirac_op.h>
#include <util/enum.h>
#include <alg/do_arg.h>
#include <util/qcd_api.h>
#include <util/qcdio.h>
#include <qalloc.h>
CPS_START_NAMESPACE
GlobalJobParameter GJP;
Verbose VRB;
Error	ERR;
LatRanGen LRG;
CPS_END_NAMESPACE

USING_NAMESPACE_CPS

static const int EVEN=0;
static const int ODD=1;
static const int ASQTAD_NUMFLOPS = 1187;
static  Matrix *gauge;
static  Vector *in_vector;
static  Vector *out_vector;
static  Vector *out_odd;
static int vol, size[4];
enum {GAUGE_LEN=18,VECT_LEN=6, VECT_LEN2=8};
static GwilsonFasqtad *lat = NULL;
static DiracOpAsqtad *asqtad = NULL;
static CgArg cg_arg;


static int LexVector(int *x){
	int i = (x[3] + size[3]*(x[0]+size[0]*(x[1]+size[1]*x[2]))) / 2 ;
	if ((x[0]+x[1]+x[2]+x[3])%2==0)
		return i; 
	else
		return  (i+vol/2);
}

static int LexGauge(int *x, int mu){
        int temp =  (x[0] + size[0]*(x[1]+size[1]*(x[2]+size[2]*x[3])));
       return (temp*4 + mu);
}
BndCndType BndCnd(int bnd){
	if (bnd == 0) return BND_CND_PRD;
	else if (bnd == 1) return BND_CND_APRD;
	else {
		printf("QCD_API:BndCnd: not a valid boundary condition (%d)\n",bnd);
		return BND_CND_PRD;
	}
}

extern "C" {

void QOP_init(QcdApiArg *arg)
{
	DoArg do_arg;
	  //----------------------------------------------------------------
	  // Initializes all Global Job Parameters
	  //----------------------------------------------------------------
#if 0
	const QMP_u32_t *size = QMP::getAllocatedDimensions();
	const QMP_u32_t *coors = QMP::getAllocatedCoordinates();
	if(ndims <4 || ndims >5){
		QMP_printf("cps_init()::Number of dimensions should be 4 or 5\n");
	}
#endif
	  int dims[8];
	  int ndims  = arg->ndims;
	  if(arg->x_sites%SizeX()!=0) printf("QOP_init:: x_sites(%d) is not divisible by the number of nodes in X dimension(%d)\n",arg->x_sites,SizeX());
	  if(arg->y_sites%SizeX()!=0) printf("QOP_init:: x_sites(%d) is not divisible by the number of nodes in X dimension(%d)\n",arg->y_sites,SizeY());
	  if(arg->z_sites%SizeX()!=0) printf("QOP_init:: x_sites(%d) is not divisible by the number of nodes in X dimension(%d)\n",arg->z_sites,SizeZ());
	  if(arg->t_sites%SizeX()!=0) printf("QOP_init:: x_sites(%d) is not divisible by the number of nodes in X dimension(%d)\n",arg->t_sites,SizeT());

	  dims[0]= do_arg.x_node_sites = arg->x_sites/SizeX();
	  dims[1]= do_arg.y_node_sites = arg->y_sites/SizeY();
	  dims[2]= do_arg.z_node_sites = arg->z_sites/SizeZ();
	  dims[3]= do_arg.t_node_sites = arg->t_sites/SizeT();
	  do_arg.s_node_sites = 1;
		vol =1;
		for(int i = 0; i<ndims; i++) vol *= dims[i];
		printf("ndims=%d dims={%d %d %d %d}\n",ndims,dims[0],dims[1],dims[2],dims[3]);

		gauge = (Matrix *)qalloc(QNONCACHE,vol*sizeof(IFloat)*GAUGE_LEN*4);
	if(gauge==NULL){
		printf("QOP_init::ran out of qalloc space for gauge field\n");exit(4);
	}


	  do_arg.x_nodes = SizeX();
	  do_arg.y_nodes = SizeY();
	  do_arg.z_nodes = SizeZ();
	  do_arg.t_nodes = SizeT();
	  do_arg.s_nodes = 1;

	  do_arg.x_bc = BndCnd(arg->x_bc);
	  do_arg.y_bc = BndCnd(arg->y_bc);
	  do_arg.z_bc = BndCnd(arg->z_bc);
	  do_arg.t_bc = BndCnd(arg->t_bc);

	  do_arg.start_conf_kind = START_CONF_LOAD;
	  do_arg.start_conf_load_addr       = gauge;
	  do_arg.start_seed_kind = START_SEED_FIXED;
//	  do_arg.colors = 3;
	  do_arg.beta = 6.0;
	  do_arg.dwf_height = 0.9;
//	  do_arg.verbose_level = DEFAULT_VERBOSE_LEVEL;
	  do_arg.asqtad_KS = arg->asqtad_KS;
	  do_arg.asqtad_naik = arg->asqtad_naik;
	  do_arg.asqtad_3staple = arg->asqtad_3staple;
	  do_arg.asqtad_5staple = arg->asqtad_5staple;
	  do_arg.asqtad_7staple = arg->asqtad_7staple;
	  do_arg.asqtad_lepage = arg->asqtad_lepage;
//	  do_arg.verbose_level = -12;

	  GJP.Initialize(do_arg);
//	printf("GJP.VerboseLevel()=%d\n",GJP.VerboseLevel());
//	VRB.Level(GJP.VerboseLevel());

	size[0] = GJP.XnodeSites();
	size[1] = GJP.YnodeSites();
	size[2] = GJP.ZnodeSites();
	size[3] = GJP.TnodeSites();
	  lat = new GwilsonFasqtad ;
}

void QOP_asqtad_dirac_init(
double * (*gauge_pt)( int, int, int, int, int, int, int, int)
){
	double dtime;
	double *gauge_p = (double *)gauge;
	double *double_p;
	int x[4],dir,a,b;

	for(x[3] = 0;x[3]<size[3] ; x[3]++)    // t
	for(x[2] = 0;x[2]<size[2] ; x[2]++)	   // z	
	for(x[1] = 0;x[1]<size[1] ; x[1]++)   // y
	for(x[0] = 0;x[0]<size[0] ; x[0]++){   // x
		for(dir = 0; dir<4;dir++){
			for(a=0;a<3;a++)
			for(b=0;b<3;b++){
				double_p = (double *) gauge_pt(x[0],x[1],x[2],x[3],dir,a,b,0);
				gauge_p[LexGauge(x,dir)*GAUGE_LEN+b*2+a*6] = *double_p;
				double_p = (double *) gauge_pt(x[0],x[1],x[2],x[3],dir,a,b,1);
				gauge_p[LexGauge(x,dir)*GAUGE_LEN+b*2+a*6+1] = (*double_p);
			}
		}
	}
	in_vector = (Vector *)qalloc(QFAST|QCOMMS,vol*sizeof(Vector));
	out_vector = (Vector *)qalloc(QFAST|QCOMMS,vol*sizeof(Vector));
	if(gauge==NULL || in_vector==NULL||out_vector==NULL){
		printf("ran out of qalloc space!\n");exit(4);
	}
	out_odd= &(out_vector[vol/2]);

//	dtime = -dclock();
	if(asqtad != NULL){
		printf("QOP_asqtad_dirac_init: initializing asqtad for the second time \n");
		exit(5);
	}	
	asqtad = new DiracOpAsqtad(*lat, out_vector,in_vector,&cg_arg,CNV_FRM_NO);
//	dtime +=dclock();
//	printf("DiracOpAsqtad:: elapsed time = %e seconds\n",dtime);
}

void QOP_asqtad_dirac_destroy(){
	if(asqtad == NULL){
		printf("QOP_asqtad_dirac_destroy: trying to delete NULL asqtad\n");
		exit(5);
	}	
	delete asqtad;
	asqtad = NULL;
}

int QOP_asqtad_inv_cg(double mass, int niter, double rsqmin, int evenodd,
double *final_rsq,
double * (*in_pt)( int, int, int, int, int, int),
double * (*out_pt)( int, int, int, int, int, int)
){
#if 1
	CommonArg common_arg;
	int x[4], dir;

	size[0] = GJP.XnodeSites();
	size[1] = GJP.YnodeSites();
	size[2] = GJP.ZnodeSites();
	size[3] = GJP.TnodeSites();

	double *in_p = (double *) in_vector;
	double *out_p = (double *) out_vector;
	int index,a,b;
	double *double_p;
	double dtime;
	if (evenodd != EVEN){
		printf("inv_cg(): non-even inversion not implemented");
		fflush(stdout);
		exit(1);
	}

//	dtime = -dclock();
//	FILE *fp = Fopen("/host/chulwoo/qcdocos_qcd_api_cps.out","w");

	for(x[3] = 0;x[3]<size[3] ; x[3]++)    // t
	for(x[2] = 0;x[2]<size[2] ; x[2]++)	   // z	
	for(x[1] = 0;x[1]<size[1] ; x[1]++)   // y
	for(x[0] = 0;x[0]<size[0] ; x[0]++){   // x
	    if( (x[0]+x[1]+x[2]+x[3])%2 == EVEN)
		for(a=0;a<3;a++){

			double_p = (double *) in_pt(x[0],x[1],x[2],x[3],a,0);
			in_p[LexVector(x)*VECT_LEN+a*2] = *double_p;
//			Fprintf(fp,"in[%d][%d]=\t%e+",LexVector(x),a, *double_p);
			double_p = (double *) in_pt(x[0],x[1],x[2],x[3],a,1);
			in_p[LexVector(x)*VECT_LEN+a*2+1] = (*double_p);
//			Fprintf(fp,"\ti%e\n",*double_p);

			double_p = (double *) out_pt(x[0],x[1],x[2],x[3],a,0);
			out_p[LexVector(x)*VECT_LEN+a*2] = *double_p;
//			Fprintf(fp,"out[%d][%d]=\t%e+",LexVector(x),a, *double_p);
			double_p = (double *) out_pt(x[0],x[1],x[2],x[3],a,1);
			out_p[LexVector(x)*VECT_LEN+a*2+1] = (*double_p);
//			Fprintf(fp,"\ti%e\n",*double_p);

		}
	}
//	Fclose(fp);
//	dtime +=dclock();
//	printf("DiracOpAsqtad:: elapsed time = %e seconds\n",dtime);


	cg_arg.mass = mass;
	cg_arg.max_num_iter = niter;
	cg_arg.stop_rsd = sqrt(rsqmin);
	asqtad->DiracArg(&cg_arg);
//common_arg.results = CAST_AWAY_CONST("/host/chulwoo/qcd_api_result");
    Float true_res;
#if 1
//	printf("sys_mhz = %d\n",sys_mhz());
//	printf("GJP.VolNodeSites() = %d\n",GJP.VolNodeSites());
//	printf("inv_cg()\n"); 
//	dtime = -dclock();
	sys_cacheflush(0);
	int iternum = asqtad->InvCg(out_vector,in_vector,0.0,&true_res);
//	dtime +=dclock();
//	long num_flops = iternum * vol* ASQTAD_NUMFLOPS;
//	printf("asqtad.InvCg():: iternum = %d elapsed time = %e seconds\n",iternum, dtime);
//	printf("flops(%ld)/%e seconds = %e MFlops/sec\n",num_flops,dtime,(double)num_flops/(dtime*1000000));
#else
	int iternum = 1;
	sys_cacheflush(0);
	asqtad->Dslash(out_odd,in_vector,(ChkbType)0,DAG_NO);
	asqtad->Dslash(out_vector,out_odd,(ChkbType)1,DAG_NO);
#endif
    *final_rsq = (double)true_res;

//	dtime = -dclock();
//	fp = Fopen("/host/chulwoo/qcdocos_qcd_api_cps.out","a");
	for(x[3] = 0;x[3]<size[3] ; x[3]++)
	for(x[2] = 0;x[2]<size[2] ; x[2]++)
	for(x[1] = 0;x[1]<size[1] ; x[1]++)
	for(x[0] = 0;x[0]<size[0] ; x[0]++){
		for(a=0;a<3;a++){
			double_p = (double *) out_pt(x[0],x[1],x[2],x[3],a,0);
			*double_p = out_p[LexVector(x)*VECT_LEN+a*2];
//			Fprintf(fp,"out[%d][%d]=\t%e+",LexVector(x),a, *double_p);
			double_p = (double *) out_pt(x[0],x[1],x[2],x[3],a,1);
			*double_p = out_p[LexVector(x)*VECT_LEN+a*2+1];
//			Fprintf(fp,"\ti%e\n",*double_p);
		}
	}
//	dtime +=dclock();
//	printf("Convert:: elapsed time = %e seconds\n",dtime);
//	Fclose(fp);

	return iternum;
#else	
	return 0;
#endif

}

void QOP_finalize(){
	if(lat !=NULL)	delete lat;
	if(asqtad !=NULL)	delete asqtad;
	qfree(gauge);
	qfree(in_vector);
	qfree(out_vector);
}

}
