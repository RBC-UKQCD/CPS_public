#include <config.h>
#include <qmp.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/smalloc.h>
#include <util/error.h>
#include <util/vector.h>
#include <alg/common_arg.h>
#include <util/dirac_op.h>
#include <alg/do_arg.h>
#include <stdio.h>

GlobalJobParameter GJP;
Verbose VRB;
Error	ERR;
LatRanGen LRG;

static const int EVEN=0;
static const int ODD=1;

static  Matrix *gauge;
static  Vector *in_vector;
static  Vector *out_vector;
static int vol, size[4];
enum {GAUGE_LEN=18,VECT_LEN=6, VECT_LEN2=8};


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

extern "C" {

void cps_init()
{
	DoArg do_arg;
	  //----------------------------------------------------------------
	  // Initializes all Global Job Parameters
	  //----------------------------------------------------------------
		const QMP_u32_t *dims = QMP::getAllocatedDimensions();
		const QMP_u32_t *coors = QMP::getAllocatedCoordinates();
		QMP_u32_t ndims = QMP::getAllocatedNumberOfDimensions();

	  do_arg.x_node_sites = dims[0];
	  do_arg.y_node_sites = dims[1];
	  do_arg.z_node_sites = dims[2];
	  do_arg.t_node_sites = dims[3];
		if(ndims >4)
	  do_arg.s_node_sites = dims[4];
		vol =1;
		for(int i = 0; i<ndims; i++) vol *= dims[i];
		printf("ndims=%d dims={%d %d %d %d}\n",ndims,dims[0],dims[1],dims[2],dims[3]);

		gauge = (Matrix *)smalloc(vol*sizeof(IFloat)*GAUGE_LEN*4);
		in_vector = (Vector *)smalloc(vol*sizeof(IFloat)*VECT_LEN);
		out_vector = (Vector *)smalloc(vol*sizeof(IFloat)*VECT_LEN);

	  do_arg.x_nodes = dims[0];
	  do_arg.y_nodes = dims[1];
	  do_arg.z_nodes = dims[2];
	  do_arg.t_nodes = dims[3];
	  do_arg.s_nodes = dims[4];

	  printf("start\n");
	  do_arg.x_bc = BND_CND_PRD;
	  do_arg.y_bc = BND_CND_PRD;
	  do_arg.z_bc = BND_CND_PRD;
	  do_arg.t_bc = BND_CND_PRD;
	  do_arg.start_conf_kind = START_CONF_LOAD;
	  do_arg.start_conf_load_addr       = gauge;
	  do_arg.start_seed_kind = START_SEED_FIXED;
	  do_arg.colors = 3;
	  do_arg.beta = 6.0;
	  do_arg.dwf_height = 0.9;
//	  do_arg.verbose_level = DEFAULT_VERBOSE_LEVEL;
	  do_arg.verbose_level = -020304;

	  GJP.Initialize(do_arg);
}

int inv_cg(double mass, int niter, double rsqmin, int evenodd,
double *final_rsq,
double * (*gauge_pt)( int, int, int, int, int, int, int, int),
double * (*in_pt)( int, int, int, int, int, int),
double * (*out_pt)( int, int, int, int, int, int)
){
#if 1
	CommonArg common_arg;
	CgArg cg_arg;
	int x[4], dir;

	printf("inv_cg()\n"); 
	size[0] = GJP.XnodeSites();
	size[1] = GJP.YnodeSites();
	size[2] = GJP.ZnodeSites();
	size[3] = GJP.TnodeSites();

	double *gauge_p = (double *)gauge;
	double *in_p = (double *) in_vector;
	double *out_p = (double *) out_vector;
	int index,a,b;
	double *double_p;
	if (evenodd != EVEN){
		printf("inv_cg(): non-even inversion not implemented");
		fflush(stdout);
		exit(1);
	}

	FILE *fp = fopen("/host/chulwoo/qcdocos_qcd_api_cps.out","w");
	for(x[3] = 0;x[3]<size[3] ; x[3]++)    // t
	for(x[2] = 0;x[2]<size[2] ; x[2]++)	   // z	
	for(x[1] = 0;x[1]<size[1] ; x[1]++)   // y
	for(x[0] = 0;x[0]<size[0] ; x[0]++)   // x
	if( (x[0]+x[1]+x[2]+x[3])%2 == EVEN){
		for(dir = 0; dir<4;dir++){
			for(a=0;a<3;a++)
			for(b=0;b<3;b++){
				double_p = (double *) gauge_pt(x[0],x[1],x[2],x[3],dir,a,b,0);
				gauge_p[LexGauge(x,dir)*GAUGE_LEN+a*2+b*6] = *double_p;
				fprintf(fp,"gauge[%d][%d]=\t%e+",LexGauge(x,dir),a+b*3, *double_p);
				double_p = (double *) gauge_pt(x[0],x[1],x[2],x[3],dir,a,b,1);
				gauge_p[LexGauge(x,dir)*GAUGE_LEN+a*2+b*6+1] = -(*double_p);
				fprintf(fp,"\ti%e\n",*double_p);
			}
		}
		for(a=0;a<3;a++){
			double_p = (double *) in_pt(x[0],x[1],x[2],x[3],a,0);
			in_p[LexVector(x)*VECT_LEN+a*2] = *double_p;
			fprintf(fp,"in[%d][%d]=\t%e+",LexVector(x),a, *double_p);
			double_p = (double *) in_pt(x[0],x[1],x[2],x[3],a,1);
			in_p[LexVector(x)*VECT_LEN+a*2+1] = (*double_p);
			fprintf(fp,"\ti%e\n",*double_p);
		}
	}
	fclose(fp);


	cg_arg.mass = mass;
	cg_arg.max_num_iter = niter;
	cg_arg.stop_rsd = rsqmin;

	GwilsonFasqtad lat;
	DiracOpAsqtad asqtad(lat, out_vector,in_vector,&cg_arg,CNV_FRM_YES);
	printf("DiracOpAsqtad\n");
	
//common_arg.results = CAST_AWAY_CONST("/host/chulwoo/qcd_api_result");
	int iternum = asqtad.InvCg(0.0,0);
	printf("asqtad.InvCg()\n");

	fp = fopen("/host/chulwoo/qcdocos_qcd_api_cps.out","a");
	for(x[3] = 0;x[3]<size[3] ; x[3]++)
	for(x[2] = 0;x[2]<size[2] ; x[2]++)
	for(x[1] = 0;x[1]<size[1] ; x[1]++)
	for(x[0] = 0;x[0]<size[0] ; x[0]++){
		for(a=0;a<3;a++){
			double_p = (double *) out_pt(x[0],x[1],x[2],x[3],a,0);
			*double_p = out_p[LexVector(x)*VECT_LEN+a*2];
			fprintf(fp,"out[%d][%d]=\t%e+",LexVector(x),a, *double_p);
			double_p = (double *) out_pt(x[0],x[1],x[2],x[3],a,1);
			*double_p = out_p[LexVector(x)*VECT_LEN+a*2+1];
			fprintf(fp,"\ti%e\n",*double_p);
		}
	}
	fclose(fp);

	return iternum;
#else	
	return 0;
#endif

}

void cps_end(){
	sfree(gauge);
	sfree(in_vector);
	sfree(out_vector);
}

}
