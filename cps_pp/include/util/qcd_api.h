#ifndef _UTIL_QCD_API_H
#define _UTIL_QCD_API_H
/* Interface between QCDOC CG inverter */

#ifdef __cplusplus
extern "C" {
#endif

struct QcdApiArg{
	int ndims;
	int x_sites;
	int y_sites;
	int z_sites;
	int t_sites;
	int s_sites;
	int x_bc;
	int y_bc;
	int z_bc;
	int t_bc;
	double asqtad_KS;
	double asqtad_naik;
	double asqtad_3staple;
	double asqtad_5staple;
	double asqtad_7staple;
	double asqtad_lepage;
};

void QOP_init(struct QcdApiArg *arg);
void QOP_finalize();

void QOP_asqtad_dirac_init(
double * (*gauge_pt)( int, int, int, int, int, int, int, int)
);
void QOP_asqtad_dirac_destroy();

int QOP_asqtad_inv_cg(double mass, int niter, double rsqmin, int evenodd,
double *final_rsq,
double * (*in_pt)( int, int, int, int, int, int),
double * (*out_pt)( int, int, int, int, int, int)
);

#ifdef __cplusplus
}
#endif

#endif
