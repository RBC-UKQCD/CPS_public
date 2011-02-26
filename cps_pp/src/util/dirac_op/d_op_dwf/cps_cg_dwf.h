#include <config.h>
#include <stdio.h>
#ifndef CPS_CG_DWF_H
#define CPS_CG_DWF_H
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/wilson.h>
#include <util/time_cps.h>
#include <util/dwf.h>
#include <mem/p2v.h>
#include <comms/glb.h>

#ifdef USE_CG_DWF
#include<util/pt.h>
#include<util/lat_data.h>
#include<dwf-ssed.h>
#include<dwf-ssef.h>

USING_NAMESPACE_CPS

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
	for(int i =0;i<5;i++)
	printf("fermion_reader %d %d: %d %d\n",UniqueID(),i,lattice_addr[i],local[i]);
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

CPS_START_NAMESPACE

#if 1
class CgDwfWrapper{
	
	private:
	static   LatMatrix *Plus[4];
	static   LatMatrix *Minus[4];
	static const char *cname;
	
	public:
	static   int lat_allocated;
	CgDwfWrapper(){}
	int Init (Lattice &lat){
		
		const char *fname="Init()";
		VRB.Func(cname,fname);
		
		LatMatrix Unit;
		Unit = 1.;
		
		{
			ParTransGauge pt_g(lat);
			Matrix *min[4], *mout_p[4], *mout_m[4];
			int dirs[]={0,2,4,6,1,3,5,7};
			for(int i =0;i<4;i++){
				if (!lat_allocated){
				Plus[i] = new LatMatrix;
//				VRB.Result(cname,fname,"Plus[%d]=%p\n",i,Plus[i]);
				Minus[i] = new LatMatrix;
//				VRB.Result(cname,fname,"Minus[%d]=%p\n",i,Minus[i]);
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
				Matrix temp;
				Matrix *mat_p = mout_m[i];
				for(int j=0;j<GJP.VolNodeSites();j++){
					temp.Dagger(*mat_p);
					*mat_p = temp;
					mat_p++;
				}
			}
			lat_allocated=1;
		}
		VRB.FuncEnd(cname,fname);
	}
	
	int Inv (Lattice *lat,Vector *out, Vector *in, CgArg *dirac_arg,unsigned long long temp_size){

		const char *fname="Inv()";
		VRB.Func(cname,fname);fflush(stdout);
		int iter;
//		unsigned long temp_size = GJP.VolNodeSites() * lat->FsiteSize() / 2;
//		VRB.Result(cname,fname,"temp_size=%d",temp_size);fflush(stdout);
		//    	const int lattice[5]={1,1,1,1,1};
		VRB.Func(cname,fname);fflush(stdout);
		int lattice[5];
		VRB.Func(cname,fname);fflush(stdout);
		for(int i=0;i<5;i++) lattice[i]=GJP.Sites(i);
	
#define DOUBLE(n) MIT_ssed_DWF_##n
#define SINGLE(n) MIT_ssef_DWF_##n
	
		VRB.Func(fname,"SINGLE(init)");
		int exit_code = SINGLE(init)(lattice,NULL,NULL);
			printf("SINGLE(init) %d\n",exit_code);
		if (exit_code) exit(-4);
	
		VRB.Func(fname,"DOUBLE(init)");
		exit_code = DOUBLE(init)(lattice,NULL,NULL);
			printf("DOUBLE(init) %d\n",exit_code);
		if (exit_code) exit(-4);
	
		VRB.Func(fname,"load_gauge");
	
		SINGLE(Gauge) *g_s = SINGLE(load_gauge)(Plus,Minus,lat,gauge_reader);
		VRB.Func(fname,"load_fermion");
		SINGLE(Fermion) *rhs_s = SINGLE(load_fermion)(in,lat,d_fermion_reader);
		VRB.Func(fname,"load_fermion");
		SINGLE(Fermion) *x0_s = SINGLE(load_fermion)(out,lat,d_fermion_reader);
		VRB.Func(fname,"load_fermion");
		SINGLE(Fermion) *x_s = SINGLE(allocate_fermion)();
	
		VRB.Func(fname,"load_gauge");
		DOUBLE(Gauge) *g = DOUBLE(load_gauge)(Plus,Minus,lat,gauge_reader);
		VRB.Func(fname,"load_fermion");
		DOUBLE(Fermion) *rhs = DOUBLE(load_fermion)(in,lat,d_fermion_reader);
		DOUBLE(Fermion) *x0 = NULL;
		VRB.Func(fname,"load_fermion");
		DOUBLE(Fermion) *x = DOUBLE(allocate_fermion)();
		VRB.Func(fname,"load_fermion");
		DOUBLE(Fermion) *Ax = DOUBLE(allocate_fermion)();
		VRB.Func(fname,"load_fermion");
		DOUBLE(Fermion) *Ax_b = DOUBLE(allocate_fermion)();
	
		double residual;
		double true_res_re,true_res_im;
		double M_0=GJP.DwfHeight();
		double m_f=dirac_arg->mass;
		double stp_cnd = dirac_arg->stop_rsd*dirac_arg->stop_rsd *in->NormSqGlbSum(2*temp_size);
		VRB.Result(cname,fname,"M_0=%g m_f=%g,stp_cnd=%e",M_0,m_f,stp_cnd);

#if 0
		VRB.Func(fname,"Dirac_Oper");
		DOUBLE(Dirac_Operator)(x,g,M_0,m_f,rhs);
#endif
	
		fermion_writer_arg save_arg;
	
		VRB.Func(fname,"cg_solver");
		SINGLE(cg_solver)(x_s, &residual,&iter,g_s,M_0,m_f,x0_s,rhs_s,
			dirac_arg->stop_rsd,0,dirac_arg->max_num_iter);
	
		save_arg.lat = lat;
#if 1
		VRB.Func(fname,"save_fermion");
		save_arg.coeff = M_0;
		SINGLE(save_fermion)(out,&save_arg,d_fermion_writer,x_s);
#else
		save_arg.coeff = 1.;
		SINGLE(save_fermion)(out,&save_arg,d_fermion_writer,x_s);

		x0 = DOUBLE(load_fermion)(out,&lat,d_fermion_reader);
		DOUBLE(cg_solver)(x, &residual,&iter,g,M_0,m_f,x0,rhs,
			dirac_arg->stop_rsd,0,dirac_arg->max_num_iter);

	//Compute true residual
		VRB.Func(fname,"Dirac_Oper");
		DOUBLE(Dirac_Operator)(Ax,g,M_0,m_f,x);
		DOUBLE(Add_Fermion)(Ax_b,Ax,-1.,rhs);
		double true_eps_re,true_eps_im;
		DOUBLE(Fermion_Dot_Product)(&true_eps_re,&true_eps_im,Ax_b,Ax_b);
		VRB.Result(cname,fname,"true_eps=%e +i %e\n",true_eps_re,true_eps_im);
	
		VRB.Func(fname,"save_fermion");
		save_arg.coeff = M_0;
		DOUBLE(save_fermion)(out,&save_arg,d_fermion_writer,x);
#endif
		QMP_barrier();
		VRB.FuncEnd(fname,"save_fermion");
	
	}

}; //class
#endif

CPS_END_NAMESPACE

#endif // USE_CG_DWF
#endif // #ifndef CPS_CG_DWF_H
