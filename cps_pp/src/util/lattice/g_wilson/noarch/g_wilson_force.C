#include <config.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/lat_data.h>
#include <util/time_cps.h>
//#include <comms/nga_reg.h>
#include <comms/glb.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE

//const unsigned CBUF_MODE4 = 0xcca52112;

#define PROFILE
//------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float dt):
// It evolves the canonical momentum mom by step_size
// using the pure gauge force.
//------------------------------------------------------------------
ForceArg Gwilson::EvolveMomGforce(Matrix *mom, Float dt){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);
  VRB.Result(cname, fname, "Entering EvolveMomGforce()\n");

  Float L1=0.0;
  Float L2=0.0;
  Float Linf=0.0;

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif
  int if_block = 0;
  if (SigmaBlockSize () > 0) if_block = 1;
  VRB.Result(cname,fname,"if_block=%d\n",if_block);

  LatData Plaqs(1);
  
  int x[4];
if (if_block){

	for (x[0] = 0; x[0] < node_sites[0]; x[0] += sigma_blocks[0]) 
	for (x[1] = 0; x[1] < node_sites[1]; x[1] += sigma_blocks[1]) 
	for (x[2] = 0; x[2] < node_sites[2]; x[2] += sigma_blocks[2]) 
	for (x[3] = 0; x[3] < node_sites[3]; x[3] += sigma_blocks[3]) {
	    int sigma = GetSigma (x, 0, 1);
	  int offset[4],x_tmp[4];
	    Float re_tr_plaq = 0.;
	for (offset[0] = 0; offset[0] < sigma_blocks[0]; offset[0] += 1) 
	for (offset[1] = 0; offset[1] < sigma_blocks[1]; offset[1] += 1) 
	for (offset[2] = 0; offset[2] < sigma_blocks[2]; offset[2] += 1) 
	for (offset[3] = 0; offset[3] < sigma_blocks[3]; offset[3] += 1) {
	    for(int i=0;i<4;i++) x_tmp[i] = x[i]+offset[i];
	    for (int mu = 0; mu < 3; ++mu)
	      for (int nu = mu + 1; nu < 4; ++nu){
		Float re_tr = ReTrPlaqNonlocal (x_tmp, mu, nu);
                re_tr_plaq += re_tr;
if (x_tmp[0]==0)
if (x_tmp[1]==0)
if (x_tmp[2]==0)
if (x_tmp[3]==0)
  VRB.Result(cname,fname,"ReTr(Plaq)[%d][%d][0]=%0.12e\n",mu,nu,re_tr);
}
}
	for (offset[0] = 0; offset[0] < sigma_blocks[0]; offset[0] += 1) 
	for (offset[1] = 0; offset[1] < sigma_blocks[1]; offset[1] += 1) 
	for (offset[2] = 0; offset[2] < sigma_blocks[2]; offset[2] += 1) 
	for (offset[3] = 0; offset[3] < sigma_blocks[3]; offset[3] += 1) {
	    for(int i=0;i<4;i++) x_tmp[i] = x[i]+offset[i];
                Float * tmp_f = (Plaqs.Field(GsiteOffset(x_tmp)/4));
	         *tmp_f =  re_tr_plaq;
	}
  }
}
  
  
  for(x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0]) {
    //printf("x[0] = %d\n", x[0]);
    for(x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1]) {
      //printf("x[1] = %d\n", x[1]);
      for(x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2]) {
	for(x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {
	  
	  int uoff = GsiteOffset(x);
          Matrix mt0;
          Matrix *mp0 = &mt0;
	  
	  for (int mu = 0; mu < 4; ++mu) {
	    if (if_block) GforceSite(*mp0, x, mu,Plaqs.Field());
	    else  GforceSite(*mp0, x, mu);
	    
	    IFloat *ihp = (IFloat *)(mom+uoff+mu);
	    IFloat *dotp = (IFloat *)mp0;
	  if(x[0]<4)
	  if(x[1]==0)
	  if(x[2]==0)
	  if(x[3]==0)
  	VRB.Result(cname,fname,"Gforce[%d][%d]=%0.12e\n",mu,x[0],mp0->norm());
	    fTimesV1PlusV2(ihp, dt, dotp, ihp, 18);
	    Float norm = ((Matrix*)dotp)->norm();
	    Float tmp = sqrt(norm);
	    L1 += tmp;
	    L2 += norm;
	    Linf = (tmp>Linf ? tmp : Linf);
	  }
	}
      }
    }
  }

#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

//  VRB.Result(cname, fname, "Finished EvolveMomGforce()\n");

  return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);

}
CPS_END_NAMESPACE
