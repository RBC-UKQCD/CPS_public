#include <config.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/time.h>
#include <comms/nga_reg.h>
#include <comms/glb.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE

//  CRAM temp buffer
#if TARGET == QCDSP
static Matrix *mp0 = (Matrix *)CRAM_SCRATCH_ADDR;	// ihdot
#else
static Matrix mt0;
static Matrix *mp0 = &mt0;		// ihdot
#endif 

const unsigned CBUF_MODE4 = 0xcca52112;

#define PROFILE
//------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float step_size):
// It evolves the canonical momentum mom by step_size
// using the pure gauge force.
//------------------------------------------------------------------
void Gwilson::EvolveMomGforce(Matrix *mom, Float step_size){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);
#ifdef PROFILE
  Float time = -dclock();
  ForceFlops=0;
#endif

  
  setCbufCntrlReg(4, CBUF_MODE4);

  int x[4];
  
  for(x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0]) {
    for(x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1]) {
      for(x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2]) {
	for(x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {
	  
	  int uoff = GsiteOffset(x);
	  
	  for (int mu = 0; mu < 4; ++mu) {
	    GforceSite(*mp0, x, mu);
	    
	    IFloat *ihp = (IFloat *)(mom+uoff+mu);
	    IFloat *dotp = (IFloat *)mp0;
	    fTimesV1PlusV2(ihp, step_size, dotp, ihp+BANK4_BASE, 18);
	  }
	}
      }
    }
  }
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif
}
CPS_END_NAMESPACE
