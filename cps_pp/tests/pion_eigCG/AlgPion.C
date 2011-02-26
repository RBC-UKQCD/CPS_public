#include"AlgPion.h"
#include<config.h>
#include<util/qcdio.h>
#include<alg/qpropw.h>
#include<alg/wilson_matrix.h>
#include<util/site.h>
#include<alg/common_arg.h>
#include<alg/no_arg.h>
#include<alg/qpropw_arg.h>
#include<util/smalloc.h>
#include<alg/fix_gauge_arg.h>
#include<util/vector.h>
#include<alg/alg_fix_gauge.h>
#include<util/rcomplex.h>
//c++ classes
#include<iostream>
//#include<fstream>
//#include<cassert>
using namespace cps;
using namespace std;

AlgPion::AlgPion(Lattice & latt, CommonArg* comm_arg, QPropWArg *lqprop_arg, EigCGArg *eigcg_arg):Alg(latt,comm_arg)
{
	cname = "AlgPion";
	char *fname = "AlgPion(Lattice &, CommonArg*, QpropWArg *)";
	VRB.Func(cname,fname);

	if(lqprop_arg == NULL)ERR.Pointer(cname,fname,"lqprop_arg");
	this->lqprop_arg = lqprop_arg;
	this->eigcg_arg = eigcg_arg;
}
AlgPion::~AlgPion()
{
	char *fname = "~AlgPion()";
	VRB.Func(cname,fname);
}
void AlgPion::runpion()
{
	char *fname = "runpion()";
	VRB.Func(cname,fname);

	std::cout<<"test: run pion before gauge fixing"<<endl;

	Lattice & lat = AlgLattice();
	//fix gauge
	FixGaugeArg fix_gauge_arg;
	fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
	fix_gauge_arg.hyperplane_start = 0;
	fix_gauge_arg.hyperplane_step = 1;
	fix_gauge_arg.hyperplane_num = GJP.TnodeSites()*GJP.Tnodes();
	fix_gauge_arg.stop_cond = 1e-8;
	fix_gauge_arg.max_iter_num = 2000;
	AlgFixGauge algfixgauge(lat,common_arg,&fix_gauge_arg);
	algfixgauge.run();	

	std::cout<<"test: run pion after gauge fixing"<<endl;
	//create qpropwwallsrc at all the time slices..
	int t_size = GJP.TnodeSites()*GJP.Tnodes();
	QPropWWallSrc *qpropw1, *qpropw2;
	char lcomm_arg_filename[200];
	CommonArg lcomm_arg("lquark","");
	lqprop_arg->t=0;
	sprintf(lcomm_arg_filename,"%s_tsrc_%d_lquark.txt",common_arg->filename,0);
	lcomm_arg.set_filename(lcomm_arg_filename);

	//eigCG para
//	int nev=(lqprop_arg->cg).bicgstab_n;
//	int m=3*nev;
//	Float max_eig=(lqprop_arg->cg).true_rsd;
//	int max_def_len=lqprop_arg->seqNum;
//	int restart_len=3;
//	Float restart[]={1e-3,1e-5,1e-7};
	int nev=eigcg_arg->nev;
	int m=eigcg_arg->m;
	Float max_eig=eigcg_arg->max_eig_cut;
	int max_def_len=eigcg_arg->max_def_len;
	int restart_len=eigcg_arg->restart_len;
	Float *restart=eigcg_arg->restart;
	bool always_restart=eigcg_arg->always_restart;

	int def_len=0;
	int vec_len=GJP.VolNodeSites() * lat.FsiteSize() / (lat.FchkbEvl()+1);
	Vector **V=new Vector*[m];
	for(int i=0;i<m;i++)
		V[i]=(Vector *)smalloc(cname,fname,"V",vec_len*sizeof(Float));
	Float *M=new Float[m];
	float **U=new float*[max_def_len];
	for(int i=0;i<max_def_len;i++)U[i]=(float *)smalloc(cname,fname,"U",vec_len*sizeof(float));
	Rcomplex *H=new Rcomplex[max_def_len*max_def_len];

	qpropw1=new QPropWWallSrc(lat,lqprop_arg,&lcomm_arg);
	qpropw1->eig_Run(V,vec_len,M,max_eig,nev,m,U,H,max_def_len,def_len,restart,restart_len,always_restart,0,1e-8);

	cout<<"deflation length:"<<def_len<<endl;
	for(int i=0;i<m;i++)sfree(V[i]);
	delete [] V;

	//second solver with restarts and stop getting more low modes
	lqprop_arg->t=6; 
	qpropw2=new QPropWWallSrc(lat,lqprop_arg,&lcomm_arg);
	qpropw2->eig_Run(NULL,vec_len,M,max_eig,0,0,U,H,max_def_len,def_len,restart,restart_len,always_restart,0,1e-8);

	algfixgauge.free();	

	Rcomplex *corr=new Rcomplex[t_size];
	WilsonMatrix Mtemp1,Mtemp2;
	char filename[200];
	
  //------------------------------------------
  //***calculate pion correlation function****
  //------------------------------------------
	
	for(int t=0;t<t_size;t++)corr[t]=0.0;
	Site s;
	for(s.Begin(); s.End(); s.nextSite())
	{
	   int t=s.physT();
	   int i=s.Index();
	   Mtemp1=(*qpropw1)[i];
	   Mtemp1.hconj();
	   corr[t]+=Trace( (*qpropw1)[i], Mtemp1 );
	}
	slice_sum((Float*)corr,t_size*2,99);

	sprintf(filename,"%s_pioncorr.txt",common_arg->filename);
	writeCorr(corr,filename);
	delete []corr;
	delete qpropw1;
	delete qpropw2;

	delete []H;
	delete []M;
	for(int i=0;i<max_def_len;i++)sfree(U[i]);
	delete [] U;
}
	
void AlgPion::writeCorr(Rcomplex *corr,char *filename)
{
	char *fname="writeCorr(*,*)";
	int t_size = GJP.TnodeSites()*GJP.Tnodes();
	
	FILE *fp;
	if((fp=Fopen(filename,"a"))==NULL)ERR.FileW(cname,fname,filename);
	for(int t=0;t<t_size;t++)
	{
		Fprintf(fp,"%d\t%e\t%e\n",t,corr[t].real(),corr[t].imag());
	}
	fflush(fp);
	Fclose(fp);
}	
 
 
 
 
 
 
 
 
