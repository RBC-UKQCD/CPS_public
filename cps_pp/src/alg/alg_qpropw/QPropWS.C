//QPropWS.C single precision contructed from QPropW class

#include <config.h> 
#include <alg/qpropw_s.h>
#include <alg/wilson_matrix.h>
#include <omp.h>

CPS_START_NAMESPACE

//constructor from QPropW
QPropWS::QPropWS(QPropW& rhs, bool rand)
{
	char *fname = "QPropWS(const QPropW &)";
	cname = "QPropWS";
	rsrc = NULL;
	int sproplen = GJP.VolNodeSites()*12*12*2;
	prop = (float *)smalloc(sproplen*sizeof(float));
	if(prop==0)ERR.Pointer(cname,fname,"prop");

#pragma omp parallel for
	for(int i=0;i<GJP.VolNodeSites();i++)
	{
		Float *x=(Float *)(&(rhs[i]));
		float *p=prop+i*288;
		for(int j=0;j<288;j++) p[j] = x[j];
	}

	if(rand)	{
		int rsrc_size = 2* GJP.VolNodeSites();
		rsrc = (Float*)smalloc(rsrc_size*sizeof(Float));
		if(rsrc==0)ERR.Pointer(cname,fname,"rsrc");

#pragma omp parallel for
		for(int i=0;i<rsrc_size/2;i++){
			rsrc[i*2]=rhs.rand_src(i).real();
			rsrc[i*2+1]=rhs.rand_src(i).imag();
		}
	}	else{
		rsrc=NULL;
	}
}	
//! copy constructor
//QPropWS(const QPropWS& rhs); 

//! Comunicate Wilson Matrices...
//WilsonMatrix& GetMatrix(const int *, WilsonMatrix&) const;

Complex QPropWS::rand_src(int i) const
{
	if(rsrc!=NULL)return Complex(rsrc[2*i], rsrc[2*i+1]);
	else {ERR.General(cname,"rand_src","rsrc is null");}
}

/*! This is a better name for the WallWallProp */
WilsonMatrix QPropWS::WallSinkProp(int t_sink)
{
	wilson_matrix wmat;
	Float *p=(Float *)(&wmat);
	for(int i=0;i<288;i++)p[i]=0.0;

	for (int i=0; i<GJP.VolNodeSites(); i++) {
		int t = i/(GJP.VolNodeSites()/GJP.TnodeSites());
		t += GJP.TnodeCoor()*GJP.TnodeSites();
		if (t != t_sink) continue;
		for(int j=0;j<288;j++)p[j] += prop[i*288+j];	  
	}

#ifdef PARALLEL
	slice_sum((Float*)&wmat, 288, 99);
	//glb_sum_multi_dir((Float *)&wmat,4,288);
#endif

	return WilsonMatrix(wmat);
}

//QPropW& operator=(const QPropW& rhs);

/*! Returns the prop */
WilsonMatrix QPropWS::operator[](int i) const
{ 
	if(i<0 || i>=GJP.VolNodeSites())ERR.General(cname,"[ ]", "i is out of range!");
	wilson_matrix w;
	Float *p=(Float *)(&w);
	for(int j=0;j<288;j++)p[j]=prop[j+i*288];
	return WilsonMatrix(w);
}

// DESTRUCTORS
QPropWS::~QPropWS()
{
	if(prop!=NULL)sfree(prop);
	prop=NULL;
	if(rsrc!=NULL)sfree(rsrc);
	rsrc=NULL;
}

CPS_END_NAMESPACE
