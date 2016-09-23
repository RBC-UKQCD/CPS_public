#ifndef ALG_WILSONFLOW
#define	ALG_WILSONFLOW
#include <config.h>
#include <alg/alg_base.h>
#include<util/site.h>

//using namespace cps;
CPS_START_NAMESPACE
class AlgWilsonFlow:public Alg
{
private:
	char *cname;
	bool su3_proj;//generally speaking it is not needed
	Float tolerance;
	Float dt;

	Matrix *lat_back;
	Float *Z_Lie;

        int Slab;
	int MatrixSize;
	int GsiteSize;
	int l_node_sites[4];
	int l_dir_offset[4];
	int vol_node_sites;
	int g_node_sites[4];
	int g_dir_offset[4];
	int g_lcl_vol;
	
	void logRun();

        void three_staple( Lattice& latt,  Matrix& link , int *pos, int u);
	void three_staple(Matrix& link, int *pos, int u, Float * gfield, const int * g_dir_offset);
	inline Float * GsiteOffset(Float * p, const int *x, const int *g_dir_offset);
	void PathOrdProdPlus(Matrix & mat, int* x, int* dirs, int n, Float *gfield, const int *g_dir_offset);
	void calculateZ(Lattice &lat, Site &site, int mu, Float Z[8]);
        void calculateZ(int pos[4], int mu, Float Z[8], Float * gfield, const int * g_dir_offset);
        void AssembleGfield(Float* lfield,Float* gfield);
        void DoRK4Step(int rk4_step, int site, Float* lfield, int l_dir_offset[4], Float* gfield, int g_dir_offset[4]);

public:
	AlgWilsonFlow(Lattice& lat, CommonArg *ca, Float dtime=0.01, bool proj=true, Float tol=1e-8);
	virtual ~AlgWilsonFlow();

	void run();
	void smartrun();

	void su3projon(){su3_proj=true;}
	void su3projoff(){su3_proj=false;}
	void set_tol(Float x){tolerance=x;}
	Float get_tol() const {return tolerance;}
	void set_dt(Float dt){this->dt = dt;}
	Float get_dt() const {return dt;}
};

CPS_END_NAMESPACE

#endif

