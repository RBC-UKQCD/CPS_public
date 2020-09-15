#include <config.h>
//------------------------------------------------------------------
//
// The class functions for FermionVectorTp.
//
//------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>

#include <alg/common_arg.h>
#include <comms/glb.h>
#include <comms/scu.h>
#include <util/site.h>
#include <util/momentum.h>
#include <comms/sysfunc_cps.h>
#include <alg/fermion_vector.h>
#include <util/omp_wrapper.h>


#include <cmath>
#include <cassert>
#include <vector>
//#include <omp.h>
//#include <qmp.h>

CPS_START_NAMESPACE
using namespace std;

const Float& FermionVectorTp::operator[](int i) {
  return fv[i];
}


FermionVectorTp::FermionVectorTp() {
  char *fname = "FermionVectorTp()";
  cname = "FermionVectorTp";
  
  VRB.Func(cname, fname);

  // allocate space for source
  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  if(GJP.Gparity()) fv_size*=2;
  fv = (Float*)smalloc(cname,fname, "fv", fv_size * sizeof(Float));


}

FermionVectorTp::~FermionVectorTp() {

  char *fname = "~FermionVectorTp()";
  
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "fv", fv);
  sfree(fv);
}

// Zero at every color,spin,space-time point
void FermionVectorTp::ZeroSource() {
  char *fname = "ZeroSource()";
  VRB.Func(cname, fname);

  uint64_t fv_size = (uint64_t) GJP.VolNodeSites() * 2 * GJP.Colors() * 4;
  if(GJP.Gparity()) fv_size*=2;

#pragma omp parallel for
  for (uint64_t i=0; i<fv_size; i++) {
	fv[i] = 0.0;
  }
}


static void SetVolSource1fGparity(Float *fv,int color, int spin,  int flavor) {
  const char *cname = "";
  const char *fname = "SetVolSource1fGparity(color,spin,flavor)";

  if(GJP.Gparity1fY() && !GJP.Gparity1fX()) ERR.General(cname,fname,"Designed for X and XY, not Y alone");

  //G-parity 1f model.
  int nsz[4] = {GJP.XnodeSites(),GJP.YnodeSites(),GJP.ZnodeSites(),GJP.TnodeSites()};
  int nodes[4] = {GJP.Xnodes(),GJP.Ynodes(),GJP.Znodes(),GJP.Tnodes()};

  static const Float quadsgn[4] = {1.0,1.0,1.0,-1.0}; //(LL,LR,UL,UR)
  int doquad[4] = {0,0,0,0};
  
  if(flavor == 0){
    doquad[0] = 1; if(GJP.Gparity1fX() && GJP.Gparity1fY()) doquad[3] = 1;
  }else{
    doquad[1] = 1; if(GJP.Gparity1fX() && GJP.Gparity1fY()) doquad[2] = 1; 
  }
  Site s ;
  for(s.Begin();s.End();s.nextSite()){
    int quad = 0;
    if(s.physX() >= nsz[0]*nodes[0]/2){
      if(GJP.Gparity1fX() && GJP.Gparity1fY() && s.physY() >= nsz[1]*nodes[1]/2) quad = 3; //UR
      else quad = 1; //LR
    }else{
      if(GJP.Gparity1fX() && GJP.Gparity1fY() && s.physY() >= nsz[1]*nodes[1]/2) quad = 2; //UL
      else quad = 0; //LL
    }
    if(!doquad[quad]) continue;

    ((Complex *)fv)[color + COLORS*(spin + 4*s.Index())] = Complex(quadsgn[quad],0.0);
  }
}

// Unit color/spin source at every space-time point
void FermionVectorTp::SetVolSource(int color, int spin, int flavor) {

  char *fname = "SetVolSource()";
  
  VRB.Func(cname, fname);
  
  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
		"Color index out of range: color = %d\n", color);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
		"Spin index out of range: spin = %d\n", spin);

  if(GJP.Gparity1fX() || GJP.Gparity1fY()) return SetVolSource1fGparity(fv,color,spin,flavor);

  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;

  int off = 0;
  if(GJP.Gparity() && flavor == 1) off = fv_size;

  for (int i = off; i < fv_size+off; i++) {
    if(i%SPINOR_SIZE == 2*(color + COLORS*spin)) {
      fv[i] = 1.0;
    }
  }
}

// Set the  color/spin source at every space-time point to a given source
void FermionVectorTp::SetVolSource(int color, int spin, Float* src)
{
  char *fname = "SetVolSource()";
  
  VRB.Func(cname, fname);

  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
    "Color index out of range: color = %d\n", color);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
    "Spin index out of range: spin = %d\n", spin);

  //zero source on all nodes
  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  if(GJP.Gparity()) fv_size*=2;

  for (int j = 0; j < fv_size; j++) {
    fv[j] = 0.0;
  }

  for (int i=0; i < GJP.VolNodeSites(); i++) {
	fv[    2*(color + COLORS*(spin + 4*i))] = src[2*i  ];
	fv[1 + 2*(color + COLORS*(spin + 4*i))] = src[2*i+1];
  }
}

static void SetWallSource1fGparity(Float *fv,int color, int spin, int source_time, int flavor){
  //G-parity 1f model.
  int sz[4] = {GJP.XnodeSites(),GJP.YnodeSites(),GJP.ZnodeSites(),GJP.TnodeSites()};
  int toff = + GJP.TnodeCoor()*sz[3];
  if(source_time < toff || source_time >= toff + sz[3]) return;

  int xy_globoff[2] = {GJP.XnodeCoor()*sz[0],GJP.YnodeCoor()*sz[1]};
  int xy_hfglobsz[2] = {GJP.Xnodes()*sz[0]/2,GJP.Ynodes()*sz[1]/2};
  if(!GJP.Gparity1fY()) xy_hfglobsz[1]*=2;

  int t = source_time - toff;

  for(int z=0;z<GJP.ZnodeSites();z++){
    for(int y=0;y<GJP.YnodeSites();y++){
      int yglob = y+xy_globoff[1];
      for(int x=0;x<GJP.XnodeSites();x++){
	int xglob = x+xy_globoff[0];
	if(flavor == 0 && xglob < xy_hfglobsz[0] && yglob < xy_hfglobsz[1]){
	  fv[2 * (color + GJP.Colors() * (spin + 4 * (x + sz[0] * ( y + sz[1] * ( z + sz[2] * t)))))] = 1.0;
	}else if(flavor == 0 && xglob >= xy_hfglobsz[0] && yglob >= xy_hfglobsz[1]){
	  fv[2 * (color + GJP.Colors() * (spin + 4 * (x + sz[0] * ( y + sz[1] * ( z + sz[2] * t)))))] = -1.0;
	}else if(flavor == 1 && xglob >= xy_hfglobsz[0] && yglob < xy_hfglobsz[1]){
	  fv[2 * (color + GJP.Colors() * (spin + 4 * (x + sz[0] * ( y + sz[1] * ( z + sz[2] * t)))))] = 1.0;	    
	}else if(flavor == 1 && xglob < xy_hfglobsz[0] && yglob >= xy_hfglobsz[1]){
	  fv[2 * (color + GJP.Colors() * (spin + 4 * (x + sz[0] * ( y + sz[1] * ( z + sz[2] * t)))))] = 1.0;
	}
      }
    }
  }
  
};


// Unit color/spin source at every space point on time_slice
// Does not zero the rest of the  time slices
void FermionVectorTp::SetWallSource(int color, int spin, int source_time, int flavor) {
  char *fname = "SetWallSource(color,spin,source_time)";
  
  VRB.Func(cname, fname);

  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
    "Color index out of range: color = %d\n", color);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
    "Spin index out of range: spin = %d\n", spin);

  if(GJP.Gparity1fX() || GJP.Gparity1fY()) return SetWallSource1fGparity(fv,color,spin,source_time,flavor);



#ifdef PARALLEL
  int my_node = GJP.TnodeCoor();
  int ts_node = (int)(source_time/GJP.TnodeSites());
#endif
  int node_ts = source_time%GJP.TnodeSites();
  int wall_size = GJP.VolNodeSites()/GJP.TnodeSites()*SPINOR_SIZE;
  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  for (int i = 0; i < fv_size; i++) {
    fv[i] = 0.0;
#ifdef PARALLEL
     if(my_node != ts_node) continue;
#endif
     if(i < wall_size*node_ts || i >= wall_size*(node_ts+1)) continue;
     if(i%SPINOR_SIZE != 2*(color + COLORS*spin)) continue;
     fv[i+flavor*fv_size] = 1.0;
  }
}

void FermionVectorTp::SetBoxSource(int color, 
				   int spin,
				   int start,
				   int end,
				   int source_time,
				   int* src_offset ) {
  
  char *fname = "SetBoxSource(color,spin,start,end,source_time)";
  
//  VRB.Func(cname, fname);
  VRB.Result(cname, fname,"Entered");
  
  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
				"Color index out of range: color = %d\n", color);
  
  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
				"Spin index out of range: spin = %d\n", spin);
  
  //#ifdef PARALLEL
  int tnode = GJP.TnodeCoor();
  int xnode = GJP.XnodeCoor();
  int ynode = GJP.YnodeCoor();
  int znode = GJP.ZnodeCoor();
  int ts_node = (int)(source_time/GJP.TnodeSites());
  //#endif
  int node_ts = source_time%GJP.TnodeSites();
  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  if(GJP.Gparity()) fv_size*=2;

  int xsize = GJP.XnodeSites();
  int ysize = GJP.YnodeSites();
  int zsize = GJP.ZnodeSites();
  //int tsize = GJP.TnodeSites();
  int xsites = GJP.XnodeSites()*GJP.Xnodes();
  int ysites = GJP.YnodeSites()*GJP.Ynodes();
  int zsites = GJP.ZnodeSites()*GJP.Znodes();
  VRB.Result(cname,fname,"sites = %d %d %d\n",xsites,ysites,zsites);
  
  // TIZB offset the start/end locations
  int gsize[3] = { xsize*GJP.Xnodes(), ysize*GJP.Ynodes(), zsize*GJP.Znodes()};
  int src_start[3] = {start,start,start};
  int src_end[3] = {end,end,end};

  if(src_offset){
    for(int i=0;i<3;++i){
      src_start[i] += src_offset[i] + gsize[i]; src_start[i] %= gsize[i];
      src_end[i] += src_offset[i] +gsize[i];  src_end[i] %= gsize[i];
    }
  }
  
  //zero source on all nodes
  for (int i = 0; i < fv_size; i++) {
     fv[i] = 0.0;
  }
  int t = node_ts;
  int src_vol=0;
  for(int x = 0; x < xsize; x++){
    for(int y = 0; y < ysize; y++){
      for(int z = 0; z < zsize; z++){
		for(int j = 0; j < GJP.Colors()*8; j++){
		  int i = j + GJP.Colors()*8*(x + xsize*(y+ysize*(z+zsize*(t))));
		  //#ifdef PARALLEL
		  if(tnode != ts_node)continue;
		  else if( (src_end[0] < xsites) && ( (x + xnode*xsize < src_start[0]) || (x + xnode*xsize > src_end[0]) )) continue;
		  else if( (src_end[1] < ysites) && ( (y + ynode*ysize < src_start[1]) || (y + ynode*ysize > src_end[1]) )) continue;
		  else if( (src_end[2] < zsites) && ( (z + znode*zsize < src_start[2]) || (z + znode*zsize > src_end[2]) )) continue;
		  else if( (src_end[0] >=xsites) && ( (x + xnode*xsize < src_start[0]) && (x + xnode*xsize > src_end[0] - xsites) )) continue;
		  else if( (src_end[1] >=ysites) && ( (y + ynode*ysize < src_start[1]) && (y + ynode*ysize > src_end[1] - ysites) )) continue;
		  else if( (src_end[2] >=zsites) && ( (z + znode*zsize < src_start[2]) && (z + znode*zsize > src_end[2] - zsites) )) continue;
		  //#endif
		  if(i%SPINOR_SIZE != 2*(color + COLORS*spin))continue;
		  //printf("BOXSRC %d %d %d %d %d %d\n", x,y,z,t, spin,color);
		  fv[i] = 1.0;
		  src_vol++;
		}
      }
    }
  }
  Float vol_f = src_vol;
  glb_sum(&vol_f);
  src_vol = vol_f;

  VRB.Result(cname,fname,"src_vol = %d\n",src_vol);
}

inline void compute_coord(int x[4], const int lcl[4], const int low[4], int i)
{
    x[0] = i % lcl[0] + low[0]; i /= lcl[0];
    x[1] = i % lcl[1] + low[1]; i /= lcl[1];
    x[2] = i % lcl[2] + low[2]; i /= lcl[2];
    x[3] = i % lcl[3] + low[3];
}

// Note: The following code sets a 4D box source. If you want to set a
// 3D xyz box, set size[3] == 1 and glb_x[3] to the global time slice
// you want.
void FermionVectorTp::Set4DBoxSource(int color,
                                     int spin,
                                     const int start[4], // global starting location in x, y, z and t directions
                                     const int size[4], // global size in x, y, z and t directions
                                     const Float mom[4]) // momentum
{
    const char *fname = "Set4DBoxSource()";
    VRB.Func(cname,fname);

    if (color < 0 || color >= GJP.Colors())
        ERR.General(cname, fname, "Color index out of range: color = %d\n", color);
  
    if (spin < 0 || spin > 3)
        ERR.General(cname, fname, "Spin index out of range: spin = %d\n", spin);
  
    ZeroSource();

    const int lcl[4] = {
        GJP.XnodeSites(), GJP.YnodeSites(),
        GJP.ZnodeSites(), GJP.TnodeSites(),
    };

    const int shift[4] = {
        GJP.XnodeSites() * GJP.XnodeCoor(), GJP.YnodeSites() * GJP.YnodeCoor(),
        GJP.ZnodeSites() * GJP.ZnodeCoor(), GJP.TnodeSites() * GJP.TnodeCoor(),
    };

    const int glb[4] = {
        GJP.XnodeSites() * GJP.Xnodes(), GJP.YnodeSites() * GJP.Ynodes(),
        GJP.ZnodeSites() * GJP.Znodes(), GJP.TnodeSites() * GJP.Tnodes(),
    };

    // need to check starting point
    for(int mu = 0; mu < 4; ++mu) {
        if(start[mu] >= 0 && start[mu] < glb[mu]) continue;
        ERR.General(cname, fname, "Invalid starting point: start[%d] = %d", mu, start[mu]);
    }

    const int sites = lcl[0] * lcl[1] * lcl[2] * lcl[3];
    // have to use Float since there is no CPS glb_sum() for int.
    Float src_vol = 0;

    VRB.Result(cname, fname, "mom = %.3e %.3e %.3e %.3e\n",
               mom[0], mom[1], mom[2], mom[3]);

#pragma omp parallel for reduction(+:src_vol)
    for(int i = 0; i < sites; ++i) {
        int glb_x[4];
        compute_coord(glb_x, lcl, shift, i);

        bool inbox = true;
        for(int mu = 0; mu < 4; ++mu) {
            if((glb_x[mu] + glb[mu] - start[mu]) % glb[mu] >= size[mu]) {
                inbox = false;
                break;
            }
        }

        if(inbox) {
            src_vol += 1;

            const double PI = 3.1415926535897932384626433832795028842;
            double alpha = 0.;
            for(int mu = 0; mu < 4; ++mu) {
                alpha += mom[mu] * 2.0 * PI * glb_x[mu] / glb[mu];
            }

            fv[i * SPINOR_SIZE + 2 * (color + COLORS * spin)    ] = std::cos(alpha);
            fv[i * SPINOR_SIZE + 2 * (color + COLORS * spin) + 1] = std::sin(alpha);
        }
    }

    glb_sum(&src_vol);
    VRB.Result(cname, fname, "src_vol = %f\n", src_vol);
    VRB.FuncEnd(cname,fname);
}


// Note: The following code sets a Z3 boxed wall source.
void FermionVectorTp::SetZ3BWall(int color, int spin, int t, const int size[3],
                                 const std::vector<Rcomplex> &rand_num)
{
    const char *fname = "SetZ3BWall()";

    if (color < 0 || color >= GJP.Colors())
        ERR.General(cname, fname, "Color index out of range: color = %d\n", color);
  
    if (spin < 0 || spin > 3)
        ERR.General(cname, fname, "Spin index out of range: spin = %d\n", spin);

    for(int mu = 0; mu < 3; ++mu) {
        if(size[mu] > 0) continue;
        ERR.General(cname, fname, "Invalid box size in %d direction: %d\n", mu, size[mu]);
    }
  
    ZeroSource();

    const int lcl[4] = {
        GJP.XnodeSites(), GJP.YnodeSites(),
        GJP.ZnodeSites(), GJP.TnodeSites(),
    };

    const int shift[4] = {
        GJP.XnodeSites() * GJP.XnodeCoor(), GJP.YnodeSites() * GJP.YnodeCoor(),
        GJP.ZnodeSites() * GJP.ZnodeCoor(), GJP.TnodeSites() * GJP.TnodeCoor(),
    };

    const int glb[4] = {
        GJP.XnodeSites() * GJP.Xnodes(), GJP.YnodeSites() * GJP.Ynodes(),
        GJP.ZnodeSites() * GJP.Znodes(), GJP.TnodeSites() * GJP.Tnodes(),
    };
    VRB.Result(cname,fname,"size=%d %d %d t=%d\n",size[0],size[1],size[2],t);
    VRB.Result(cname,fname,"lcl=%d %d %d %d\n",lcl[0],lcl[1],lcl[2],lcl[3]);
    VRB.Result(cname,fname,"shift=%d %d %d %d\n",shift[0],shift[1],shift[2],shift[3]);
    VRB.Result(cname,fname,"glb=%d %d %d %d\n",glb[0],glb[1],glb[2],glb[3]);




// Changed from Hantao's version to avoid putting sources in remainder portion of the 3D slices when size[] does not divide glb[]

    const int rand_grid[3] = {
        (glb[0] ) / size[0], 
        (glb[1] ) / size[1], 
        (glb[2] ) / size[2]
    };

    const int sites = lcl[0] * lcl[1] * lcl[2] * lcl[3];

#pragma omp parallel for
    for(int i = 0; i < sites; ++i) {
        int glb_x[4];
        compute_coord(glb_x, lcl, shift, i);

        if(glb_x[3] != t) continue;
    for(int mu = 0; mu < 3; ++mu) 
        if(glb_x[mu] >= (size[mu]*rand_grid[mu]) ) continue;

        int id = 0;
        for(int k = 0; k < 3; ++k) {
            id = id * rand_grid[k] + glb_x[k] / size[k];
        }

        fv[i * SPINOR_SIZE + 2 * (color + COLORS * spin)    ] = std::real(rand_num[id]);
        fv[i * SPINOR_SIZE + 2 * (color + COLORS * spin) + 1] = std::imag(rand_num[id]);
    }

    // debug code, check the source by printing it.
     for(int i = 0; i < sites; ++i) {
         int glb_x[4];
        compute_coord(glb_x, lcl, shift, i);
         if(glb_x[3] != t) continue;

         printf("Z3B Source: %d %d %d %d = %17.10e %17.10e\n",
                glb_x[0], glb_x[1], glb_x[2], glb_x[3],
                fv[i * SPINOR_SIZE + 2 * (color + COLORS * spin)    ],
                fv[i * SPINOR_SIZE + 2 * (color + COLORS * spin) + 1]);
     }
}

// Set source from previously defined source
// Does not zero the rest of the  time slices
void FermionVectorTp::SetWallSource(int color, int spin, int source_time, 
									Float* src) {
  char *fname = "SetWallSource(color,spin,source_time,src)";
  
  VRB.Func(cname, fname);

  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
				"Color index out of range: color = %d\n", color);
  
  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
				"Spin index out of range: spin = %d\n", spin);
#ifdef PARALLEL
  int my_node = GJP.TnodeCoor();
  int ts_node = (int)(source_time/GJP.TnodeSites());
#endif
  int node_ts = source_time%GJP.TnodeSites();
  int wall_size = GJP.VolNodeSites()/GJP.TnodeSites();

  //zero source on all nodes
  //int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  //for (int i = 0; i < fv_size; i++) {
  //  fv[i] = 0.0;
  //}

#ifdef PARALLEL
  if(my_node == ts_node)
#endif
    {
      for (int i=0; i < GJP.VolNodeSites(); i++) {
		if(i/wall_size != node_ts) continue; // src defined over whole node!
	fv[    2*(color + COLORS*(spin + 4*i))] = src[2*i  ];
	fv[1 + 2*(color + COLORS*(spin + 4*i))] = src[2*i+1];
      }
    }

}

// COULOMB GAUGE ONLY!
void FermionVectorTp::GFWallSource(Lattice &lat, int spin, int dir, int where, int flavor) {
    char *fname = "GFWallSource()";
    VRB.Func(cname, fname);
    VRB.Debug(cname,fname,"lat=%p spin=%d dir=%d  where=%d\n",&lat,spin,dir,where);

//    if(dir != 3) { ERR.NotImplemented(cname, fname);
  if(dir!=3) ERR.General(cname,fname,"Works only for dir=3\n"); //Added by CK, see site loop below to verify

//the local (on processor) length in "dir" direction
    int len= GJP.NodeSites(dir);
// local processor coordinate in d_ direction
  int lproc = GJP.NodeCoor(dir);

  //CK remove nasty switch in favour of improved GJP commands

    // nc: node coordinate
    // lc: local coordinate
    int nc = where / GJP.TnodeSites();
    int lc = where % GJP.TnodeSites();

    if(nc != GJP.TnodeCoor()) return; // nothing to do here.

    Matrix **gm = lat.FixGaugePtr();
    if (!gm) ERR.Pointer(cname,fname,"fix_gauge_ptr");
//    printf("gm(%d)=%p\n",UniqueID(),gm);
#ifdef USE_OMP
    Matrix *pM = gm[lc];
    VRB.Debug(cname,fname,"pM=%p\n",pM);
    int vol_3d = GJP.XnodeSites() * GJP.YnodeSites() * GJP.ZnodeSites();
#pragma omp parallel for
    for(int i = 0; i < vol_3d; ++i) {
        int mid = i;
        int vid = 6 * spin + SPINOR_SIZE * (i + vol_3d * lc);

        Vector *v = (Vector *)(fv + vid);
        Vector vt(*v);

        Matrix mt;
        mt.Dagger(pM[mid]);

        v->DotXEqual(mt, vt);
   }
#else
  Vector temp;
  Matrix tempmat; 
  // find out if this node overlaps with the hyperplane
  // in which the wall source sits
  int has_overlap = 0;
//  if (lproc * len <= where && where < (lproc + 1) * len)
    if (lproc == nc)
    has_overlap = 1;
 
  if (has_overlap) {
    int local = where % len; // on processor coordinate of
    // source hyperplane
    
    Matrix *pM; 
    if(GJP.Gparity()) pM = gm[local + flavor * len]; //CK G-parity, gauge fixing matrices for U* field (flavor==1) are stored after the U hyperplanes
    else pM = gm[local];

    for (int z = 0; z < GJP.ZnodeSites(); z++)
    for (int y = 0; y < GJP.YnodeSites(); y++) 
    for (int x = 0; x < GJP.XnodeSites(); x++)
    {
      // the matrix offset
      int j =  x + GJP.XnodeSites() * ( y + GJP.YnodeSites() * z);
      // the vector offset
      int i = 2 * GJP.Colors() * ( spin + 4 * (
              x + GJP.XnodeSites() * (
              y + GJP.YnodeSites() * (
              z + GJP.ZnodeSites() * local))));
      if(GJP.Gparity() && flavor == 1) i+=GJP.VolNodeSites() * GJP.Colors() * 8; //CK: skip on field

      if(pM == NULL) ERR.General(cname,fname,"Non-zero fermion vector at site %d but no gauge fixing matrix available");

      temp.CopyVec((Vector*)&fv[i], 6);
      tempmat.Dagger((IFloat*)&pM[j]);
      uDotXEqual((IFloat*)&fv[i], (const IFloat*)&tempmat, (const IFloat*)&temp);

      //printf("FT:GFWALL %d %d %d %d : ",x,y,z,spin);
      //for(int col=0;col<6;++col){ printf("%e ", *(col+(Float*)&(fv[i]))); }
      //printf("\n");
    }
  }
#endif
}

//Gauge fix fermion vector, works for any gauge fixing
void FermionVectorTp::GaugeFixVector(Lattice &lat, int spin){
  static const char* fname = "GaugeFixVector(Lattice &lat, int spin)";
  if(lat.FixGaugeKind() == FIX_GAUGE_NONE) return;

  int nflav = 1;
  if(GJP.Gparity()) nflav = 2;

#pragma omp parallel for default(shared)
  for(int site=0;site<GJP.VolNodeSites();site++){
    for(int flavor = 0; flavor < nflav; flavor ++){
      int i = 2 * GJP.Colors() * ( spin + 4 *(site + flavor * GJP.VolNodeSites()) );
      bool iszero(true);
      for(int ii=0;ii<6;ii++) if(fv[i+ii]!=0.0){ iszero=false; break; }
      if(iszero) continue; //skip if fermion vector is zero at this site

      const Matrix* gfmat = lat.FixGaugeMatrix(site,flavor);
      if(gfmat == NULL) ERR.General(cname,fname,"Non-zero fermion vector at site %d but no gauge fixing matrix available");

      Vector temp;
      temp.CopyVec((Vector*)&fv[i], 6);
      Matrix tempmat;
      tempmat.Dagger((const IFloat*)gfmat );
      uDotXEqual((IFloat*)&fv[i], (const IFloat*)&tempmat, (const IFloat*)&temp);
    }
  }
}





// COULOMB GAUGE ONLY!
void FermionVectorTp::GaugeFixSink(Lattice &lat, int dir, int unfix) {
  
  char *fname = "GaugeFixSink()";
  VRB.Func(cname, fname);

    if(dir!=3) {
        ERR.General(cname,fname,"Works only for dir=3\n");
    }
    if(lat.FixGaugeKind()!=FIX_GAUGE_COULOMB_T) {
        ERR.General(cname,fname,"Works only for FIX_GAUGE_COULOMB_T\n");
    }

    Matrix **gm = lat.FixGaugePtr();
  int nflav=1;
  if(GJP.Gparity()) nflav = 2;

  for(int flav = 0; flav < nflav; flav++){
    for (int t=0; t < GJP.TnodeSites(); t++) {

      Matrix *pM = gm[t + flav * GJP.TnodeSites()]; //flavour 1 GF matrices stored after flavour 0 hyperplanes
	Vector temp;

      if(gm[t] != NULL ){
	for (int z = 0; z < GJP.ZnodeSites(); z++)
	  for (int y = 0; y < GJP.YnodeSites(); y++)
	    for (int x = 0; x < GJP.XnodeSites(); x++){
          // the matrix offset
		int j =  x + GJP.XnodeSites() * ( y + GJP.YnodeSites() * z);

	      for (int spin = 0; spin < 4; spin++){
		// the vector offset
		int i= 2 * GJP.Colors() * ( spin + 4 * (
							x + GJP.XnodeSites() * (
										y + GJP.YnodeSites() * (
													z + GJP.ZnodeSites() * t)))) ;
		if(flav) i+= 2 * GJP.Colors() * 4 * GJP.VolNodeSites();

		temp.CopyVec((Vector*)&fv[i], 6);
		if(unfix)
		uDagDotXEqual((IFloat*)&fv[i],(const IFloat*)&pM[j],
				(const IFloat*)&temp);
		else
		uDotXEqual((IFloat*)&fv[i],(const IFloat*)&pM[j],
			     (const IFloat*)&temp);
	      }
        }
      }
    }
  }
}

void FermionVectorTp::LandauGaugeFixSrc(Lattice& lat, int spin)
{
    //Landau gauge fixing the source added by Qi
    char *fname = "LandauGaugeFixSrc()";
    VRB.Func(cname, fname);

    if (lat.FixGaugeKind() != FIX_GAUGE_LANDAU)
	ERR.General(cname, fname, "lattice not in Landau gauge\n");

    Matrix *pM(lat.FixGaugePtr()[0]);

#pragma omp parallel for
    for (int site = 0; site < GJP.VolNodeSites(); site++) {
	// site offset
	Matrix Adj;
	Adj.Dagger(pM[site]);
	const int f_off(2 * GJP.Colors() * (spin + 4 * site));
	Vector *v = (Vector*)(fv + f_off);
	Vector vt(*v);
	uDotXEqual((IFloat*)&fv[f_off], (const IFloat*)&Adj, (const IFloat*)&vt);
    }
}

void FermionVectorTp::LandauGaugeFixSink( Lattice& lat )
{
  char *fname = "LandauGaugeFixSink()";
  VRB.Func(cname, fname);
  
  if(lat.FixGaugeKind()!=FIX_GAUGE_LANDAU)
    ERR.General(cname,fname,"lattice not in Landau gauge\n");

  //For G-parity, Landau gauge fixing matrices are stored at offset Lx*Ly*Lz*Lt
  Matrix *pM(lat.FixGaugePtr()[0]);

  Vector temp;

  int nflav = 1;
  if(GJP.Gparity()) nflav = 2;

  for(int flav = 0; flav < nflav ; flav++){
    Site site;
    while (site.LoopsOverNode()) {
      // site offset 
      const int s_off(site.Index() + flav*GJP.VolNodeSites());
     
      int spin;
      for (spin=0; spin<4; spin++)
	{
	  // the start of the vector offset in floats 
	  // 6 == re/im * colours
	  const int f_off( 2 * GJP.Colors() * ( spin + 4*s_off) );
		
	  temp.CopyVec((Vector*)&fv[f_off], 6 );
		
	  uDotXEqual( (IFloat*)      &fv[f_off],
		      (const IFloat*)&pM [s_off],
		      (const IFloat*)&temp       );
	}
    }
  }

}

void FermionVectorTp::SetLandauGaugeMomentaSource( Lattice& lat,
                                                   int src_colour,
                                                   int src_spin,
                                                   int p[],
						   int flavor) {

  char *fname = "LandauGaugeMomentaSrc()";
  VRB.Func(cname, fname);
  
  if(lat.FixGaugeKind()!=FIX_GAUGE_LANDAU)
    ERR.General(cname,fname,"lattice not in Landau gauge\n");
  
  
  Matrix* pM( lat.FixGaugePtr()[0] );
    
  /*
    to hold the src (colour vector by colour vector to make
    gauge fixing easier
  */
  Vector temp;

  // gauge fixing matrix
  Matrix Adj ;
  
  Site site  ;
  
  const Float PI(3.14159265358979323846264338327950288319716939937510);

  Float pp[4];
  for(int i=0;i<4;i++){
    pp[i] = p[i]*2.0*PI/(GJP.NodeSites(i)*GJP.Nodes(i));
    if(GJP.Bc(i)==BND_CND_GPARITY) pp[i]/=4.0; //CK: units of momentum are pi/2L not 2pi/L in G-parity directions  
  }

  while ( site.LoopsOverNode() )
    {
      // site offset 
      const int s_off( site.Index() + flavor * GJP.VolNodeSites() );
      Adj.Dagger(pM[s_off]);

      // work out the momentum
      const Float pdotx( pp[0]*site.physX() + pp[1]*site.physY() + pp[2]*site.physZ() + pp[3]*site.physT() );      
      const Rcomplex fact( cos(pdotx), sin(pdotx) );
      
      int spin,colour;
      for ( spin=0; spin<4; spin++ )
        {
          const int f_off( 2 * GJP.Colors() * ( spin + 4*s_off) );
          Vector & cvec(*((Vector*)&fv[f_off]));
          cvec.Zero();

          for ( colour=0;colour<3;colour++)
            {
              if ( src_spin == spin && src_colour == colour )
                {
                  cvec[colour] = fact;
                  temp.CopyVec(&cvec,6);
                  uDotXEqual( (IFloat*)      &fv[f_off],
                              (const IFloat*)&Adj       ,
                              (const IFloat*)&temp       );
                }
            } // colour
        } // spin
    }
}

// Gaussian smear previously defined source (or sink)
// * Does not zero the rest of the  time slices
// * Smearing is done on a color vector across a time-slice
//   for a definite value of the spin, so the source is a 12 component
//   object at each site which should have been set previously,
//   that is, we smear fv, the member data of the FermionVector that
//   called this routine.
// * based on Smearing note by K. Orginos, 9/29/2004
// * essentially, hit a vector with the gauge invariant Laplacian:
//                    __
//     S = 1 + w^2/4N \/^2
//     __
//     \/^2 psi = Sum_{i=1,2,3} {U_i(x)psi(x+i) + U^\dagger_i(x-i)psi(x-i)
//                               - 2 psi
// Typical vales: w=4.35, N=40
//
// TB 9/2004
void FermionVectorTp::GaussianSmearVector(Lattice& lat, 
					  int spin, 
					  int Niter, 
					  Float omega, 
					  int source_time)
{
  char *fname = "GaussianSmearedVector(int spin, int source_time)";
  
  VRB.Func(cname, fname);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
		"Spin index out of range: spin = %d\n", spin);
#ifdef PARALLEL
  int my_node=GJP.TnodeCoor();
  int ts_node = (int)(source_time/GJP.TnodeSites());
#endif
  int node_ts = source_time%GJP.TnodeSites();
  
  //Niter = 40;
  //omega = 4.35;
  // the smearing factor 
  Float factor = omega*omega/(4*Niter);

  // the following code is from d_stag_opt_rdm dslash,
  // modified

  const int VECT_LEN=6;

  int nx[4];
  int nb[4];
  int xv[4];
  
  //-----------------------------------------------------------
  //  nx[4]
  //-----------------------------------------------------------
  nx[0] = GJP.XnodeSites();
  nx[1] = GJP.YnodeSites();
  nx[2] = GJP.ZnodeSites();
  nx[3] = GJP.TnodeSites();
  
  //-----------------------------------------------------------
  //  nb[4]
  //-----------------------------------------------------------
  nb[0] = 4;
  nb[1] = nb[0]*nx[0];
  nb[2] = nb[1]*nx[1];
  nb[3] = nb[2]*nx[2];
  
  //-----------------------------------------------------------
  //  xv[4]
  //-----------------------------------------------------------
  xv[0] = 1;
  xv[1] = xv[0]*nx[0];
  xv[2] = xv[1]*nx[1];
  xv[3] = xv[2]*nx[2];

  
  int x[4], offset;

  const Matrix *uoff;
  const Matrix *mp0;

  Vector vtmp0, vtmp1;
  Vector *src;  
  Vector *vp0;
  Vector *vp1;
  // working vector for the smeared source
  Vector* chi;
 
  src = (Vector*) smalloc(cname,fname, "src", sizeof(Vector)*nx[0]*nx[1]*nx[2]);
  chi = (Vector*) smalloc(cname,fname, "chi", sizeof(Vector)*nx[0]*nx[1]*nx[2]);

  // do the smearchig on all nodes since they have to wait anyway.
  // In the end we take only the result for the source time_slice

  // set time-slice
  x[3] = node_ts;

  // do Niter hits on src
  for(int n=0;n<Niter;n++){
    
    //initialize src
    for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
      for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
	for(x[2] = 0; x[2] < nx[2]; ++x[2]) {

	  int j=xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	  int i=j+xv[3]*x[3];
	  *(src+j)=*((Vector *)fv + spin + 4 * i);
	}
      }
    }

    // loop over all sites on time-slice
    for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
      for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
	for(x[2] = 0; x[2] < nx[2]; ++x[2]) {
	  
	  //the offset of the link at lattice site
	  offset = nb[0]*x[0]+nb[1]*x[1]+nb[2]*x[2]+nb[3]*x[3];
	  // the gauge field link, obviously
	  uoff = (Matrix *) lat.GaugeField() + offset;

	  for (int mu = 0; mu < 3; ++mu) {
	    
	    //-------------------------------------------
	    //  calculate U_mu(x) chi(x+mu)
	    //-------------------------------------------
	    
	    // gather chi(x+mu)
	    if(x[mu] == nx[mu]-1) { 	// x+mu off node
	      x[mu] = 0;
	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      getPlusData((IFloat *)&vtmp0, 
			  (IFloat *)src + VECT_LEN*offset, 
			  VECT_LEN, mu);
	      x[mu] = nx[mu]-1;
	      vp0 = &vtmp0;
	      
	    } else { 			// x+mu on node
	      x[mu]++;
	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      vp0 = src + offset;
	      x[mu]--;
	    }
	    
	    mp0 = uoff+mu;
	    
	    // link x chi
	    offset = VECT_LEN*(xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2]);
	    if(mu == 0) uDotXEqual((IFloat *)chi+offset, (const IFloat *)mp0,
				   (const IFloat *)vp0);
	    else uDotXPlus((IFloat *)chi+offset,(const IFloat *)mp0,
			   (const IFloat *)vp0);
	    
	    
	    //-------------------------------------------
	    //  calculate U^dag_mu(x-mu) chi(x-mu)
	    //-------------------------------------------
	    if(x[mu] == 0) { 		// x-mu off node
	      x[mu] = nx[mu]-1;
	      mp0 = uoff+x[mu]*nb[mu]+mu;
	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      vp1 = (Vector *)src + offset;
	      x[mu] = 0;
	      
	      uDagDotXEqual((IFloat *)&vtmp0, (const IFloat *)mp0,
			    (const IFloat *)vp1);
	      
	      getMinusData((IFloat *)&vtmp1, (IFloat *)&vtmp0, VECT_LEN, mu);

	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      *(chi+offset) += vtmp1;
	      
	    } else { 			// x-mu on node
	      
	      x[mu]--;
	      mp0 = uoff-nb[mu]+mu;
	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      vp0 = (Vector *)src + offset;
	      x[mu]++;
	      
	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      uDagDotXPlus((IFloat *)chi+VECT_LEN*offset,(const IFloat *)mp0,
			   (const IFloat *)vp0);
	    }
	  } // mu
	  
	  offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	  *(chi+offset) *= factor;

	  // on site term is just (1-6*w^2/(4N))*src
	  offset *=VECT_LEN;
	  fTimesV1PlusV2((IFloat*)chi+offset, 1.-6.*factor, 
			 (IFloat*)src+offset, 
			 (IFloat*)chi+offset, VECT_LEN);
	  //debug: 
	  /*printf("SRC, CHI (%d,%d,%d)= %e %e %e %e\n", 
		 x[0], x[1], x[2], 
		 *((IFloat*)src+offset), 
		 *((IFloat*)chi+offset),
		 *((IFloat*)chi+offset+2),
		 *((IFloat*)chi+offset+4));*/
	}
      }
    } // end loop over sites on time slice

    //copy chi into fv
    for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
      for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
	for(x[2] = 0; x[2] < nx[2]; ++x[2]) {
	  
	  // only save if we are on the right node
#ifdef PARALLEL
	  if(my_node != ts_node)continue; 
#endif
	  int j=xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	  int i=j+xv[3]*x[3];
	  *((Vector *)fv + spin + 4 * i)=*(chi+j);
	}
      }
    }
  } // iters
  
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "chi", chi);
  sfree(chi);
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "src", src);
  sfree(src); 

  //multiply  fv by 1000. Makes the source a little larger to preven
  //numerical instabilities
#if 0
  for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
    for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
      for(x[2] = 0; x[2] < nx[2]; ++x[2]) {
	// only save if we are on the right node                              
#ifdef PARALLEL
	if(my_node != ts_node)continue;
#endif
	int i=xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2] + xv[3]*x[3];
	*((Vector *)fv + spin + 4 * i)*= 1.0e3 ;
      }
    }
  }
#endif

#ifdef DEBUG_GAUSS_SMEAR
Float nn = 0 ;
for(x[0] = 0; x[0] < nx[0]; ++x[0]) 
   for(x[1] = 0; x[1] < nx[1]; ++x[1]) 
     for(x[2] = 0; x[2] < nx[2]; ++x[2]) 
       for(x[3] = 0; x[3] < nx[3]; ++x[3]) {
	 int i = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2]+xv[3]*x[3];
	 Vector foo = *((Vector *)fv + spin + 4 * i);	
	 nn += foo.NormSqNode(VECT_LEN) ;
	 printf("%i %i %i %i: %g\n",x[0],x[1],x[2],x[3],foo.NormSqNode(VECT_LEN)); 
       }
 printf("SOURCE NORM: %g\n", nn);
#endif
}

// Gaussian smear previously defined source (or sink)
// * Does not zero the rest of the  time slices
// * Smearing is done on a color vector across a time-slice
//   for a definite value of the spin, so the source is a 12 component
//   object at each site which should have been set previously,
//   that is, we smear fv, the member data of the FermionVector that
//   called this routine.
// * based on Smearing note by K. Orginos, 9/29/2004
// * essentially, hit a vector with the gauge invariant Laplacian:
//                    __
//     S = 1 + w^2/4N \/^2
//     __
//     \/^2 psi = Sum_{i=1,2,3} {U_i(x)psi(x+i) + U^\dagger_i(x-i)psi(x-i)
//                               - 2 psi
// Typical vales: w=4.35, N=40
//
// TB 9/2004
void FermionVectorTp::GaussianSmearVector(Lattice& lat, 
					  int spin, 
					  int Niter, 
					  Float omega)
{
  char *fname = "GaussianSmearedVector(int spin, int source_time)";
  
  VRB.Func(cname, fname);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
		"Spin index out of range: spin = %d\n", spin);
  /*
#ifdef PARALLEL
  int my_node=GJP.TnodeCoor();
  int ts_node = (int)(source_time/GJP.TnodeSites());
#endif
  int node_ts = source_time%GJP.TnodeSites();
  */
  //Niter = 40;
  //omega = 4.35;
  // the smearing factor 
  Float factor = omega*omega/(4*Niter);

  // the following code is from d_stag_opt_rdm dslash,
  // modified

  const int VECT_LEN=6;

  int nx[4];
  int nb[4];
  int xv[4];
  
  //-----------------------------------------------------------
  //  nx[4]
  //-----------------------------------------------------------
  nx[0] = GJP.XnodeSites();
  nx[1] = GJP.YnodeSites();
  nx[2] = GJP.ZnodeSites();
  nx[3] = GJP.TnodeSites();
  
  //-----------------------------------------------------------
  //  nb[4]
  //-----------------------------------------------------------
  nb[0] = 4;
  nb[1] = nb[0]*nx[0];
  nb[2] = nb[1]*nx[1];
  nb[3] = nb[2]*nx[2];
  
  //-----------------------------------------------------------
  //  xv[4]
  //-----------------------------------------------------------
  xv[0] = 1;
  xv[1] = xv[0]*nx[0];
  xv[2] = xv[1]*nx[1];
  xv[3] = xv[2]*nx[2];

  
  int x[4], offset;

  const Matrix *uoff;
  const Matrix *mp0;

  Vector vtmp0, vtmp1;
  Vector *src;  
  Vector *vp0;
  Vector *vp1;
  // working vector for the smeared source
  Vector* chi;
 
  src = (Vector*) smalloc(cname,fname, "src", sizeof(Vector)*nx[0]*nx[1]*nx[2]);
  chi = (Vector*) smalloc(cname,fname, "chi", sizeof(Vector)*nx[0]*nx[1]*nx[2]);

  // do the smearchig on all nodes since they have to wait anyway.
  // In the end we take only the result for the source time_slice

  // set time-slice
  //x[3] = node_ts;
  for(x[3]=0;x[3]<GJP.TnodeSites();x[3]++) {
  // do Niter hits on src
  for(int n=0;n<Niter;n++){
    
    //initialize src
    for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
      for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
	for(x[2] = 0; x[2] < nx[2]; ++x[2]) {

	  int j=xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	  int i=j+xv[3]*x[3];
	  *(src+j)=*((Vector *)fv + spin + 4 * i);
	}
      }
    }

    // loop over all sites on time-slice
    for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
      for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
	for(x[2] = 0; x[2] < nx[2]; ++x[2]) {
	  
	  //the offset of the link at lattice site
	  offset = nb[0]*x[0]+nb[1]*x[1]+nb[2]*x[2]+nb[3]*x[3];
	  // the gauge field link, obviously
	  uoff = (Matrix *) lat.GaugeField() + offset;

	  for (int mu = 0; mu < 3; ++mu) {
	    
	    //-------------------------------------------
	    //  calculate U_mu(x) chi(x+mu)
	    //-------------------------------------------
	    
	    // gather chi(x+mu)
	    if(x[mu] == nx[mu]-1) { 	// x+mu off node
	      x[mu] = 0;
	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      getPlusData((IFloat *)&vtmp0, 
			  (IFloat *)src + VECT_LEN*offset, 
			  VECT_LEN, mu);
	      x[mu] = nx[mu]-1;
	      vp0 = &vtmp0;
	      
	    } else { 			// x+mu on node
	      x[mu]++;
	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      vp0 = src + offset;
	      x[mu]--;
	    }
	    
	    mp0 = uoff+mu;
	    
	    // link x chi
	    offset = VECT_LEN*(xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2]);
	    if(mu == 0) uDotXEqual((IFloat *)chi+offset, (const IFloat *)mp0,
				   (const IFloat *)vp0);
	    else uDotXPlus((IFloat *)chi+offset,(const IFloat *)mp0,
			   (const IFloat *)vp0);
	    
	    
	    //-------------------------------------------
	    //  calculate U^dag_mu(x-mu) chi(x-mu)
	    //-------------------------------------------
	    if(x[mu] == 0) { 		// x-mu off node
	      x[mu] = nx[mu]-1;
	      mp0 = uoff+x[mu]*nb[mu]+mu;
	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      vp1 = (Vector *)src + offset;
	      x[mu] = 0;
	      
	      uDagDotXEqual((IFloat *)&vtmp0, (const IFloat *)mp0,
			    (const IFloat *)vp1);
	      
	      getMinusData((IFloat *)&vtmp1, (IFloat *)&vtmp0, VECT_LEN, mu);

	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      *(chi+offset) += vtmp1;
	      
	    } else { 			// x-mu on node
	      
	      x[mu]--;
	      mp0 = uoff-nb[mu]+mu;
	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      vp0 = (Vector *)src + offset;
	      x[mu]++;
	      
	      offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	      uDagDotXPlus((IFloat *)chi+VECT_LEN*offset,(const IFloat *)mp0,
			   (const IFloat *)vp0);
	    }
	  } // mu
	  
	  offset = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	  *(chi+offset) *= factor;

	  // on site term is just (1-6*w^2/(4N))*src
	  offset *=VECT_LEN;
	  fTimesV1PlusV2((IFloat*)chi+offset, 1.-6.*factor, 
			 (IFloat*)src+offset, 
			 (IFloat*)chi+offset, VECT_LEN);
	  //debug: 
	  /*printf("SRC, CHI (%d,%d,%d)= %e %e %e %e\n", 
		 x[0], x[1], x[2], 
		 *((IFloat*)src+offset), 
		 *((IFloat*)chi+offset),
		 *((IFloat*)chi+offset+2),
		 *((IFloat*)chi+offset+4));*/
	}
      }
    } // end loop over sites on time slice

    //copy chi into fv
    for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
      for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
	for(x[2] = 0; x[2] < nx[2]; ++x[2]) {
	  
	  // only save if we are on the right node
	  //#ifdef PARALLEL
	  //if(my_node != ts_node)continue; 
	  //#endif
	  int j=xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2];
	  int i=j+xv[3]*x[3];
	  *((Vector *)fv + spin + 4 * i)=*(chi+j);
	}
      }
    }
  } // iters
  } // time loop
  
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "chi", chi);
  sfree(chi);
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "src", src);
  sfree(src); 

  //multiply  fv by 1000. Makes the source a little larger to preven
  //numerical instabilities
#if 0
  for(x[3]=0;x[3]<GJP.TnodeSites();x[3]++) {
  for(x[0] = 0; x[0] < nx[0]; ++x[0]) {
    for(x[1] = 0; x[1] < nx[1]; ++x[1]) {
      for(x[2] = 0; x[2] < nx[2]; ++x[2]) {
	// only save if we are on the right node                              
	//#ifdef PARALLEL
	//if(my_node != ts_node)continue;
	//#endif
	int i=xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2] + xv[3]*x[3];
	*((Vector *)fv + spin + 4 * i)*= 1.0e3 ;
      }
    }
  }
  }
#endif 

#ifdef DEBUG_GAUSS_SMEAR
Float nn = 0 ;
for(x[0] = 0; x[0] < nx[0]; ++x[0]) 
   for(x[1] = 0; x[1] < nx[1]; ++x[1]) 
     for(x[2] = 0; x[2] < nx[2]; ++x[2]) 
       for(x[3] = 0; x[3] < nx[3]; ++x[3]) {
	 int i = xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2]+xv[3]*x[3];
	 Vector foo = *((Vector *)fv + spin + 4 * i);	
	 nn += foo.NormSqNode(VECT_LEN) ;
	 printf("%i %i %i %i: %g\n",x[0],x[1],x[2],x[3],foo.NormSqNode(VECT_LEN)); 
       }
 printf("SOURCE NORM: %g\n", nn);
#endif
}

Float FermionVectorTp::Norm(){
  Site s;
  Float norm(0);
  for(s.Begin();s.End();s.nextSite()){
    for(int spin(0);spin<4;spin++)
      norm += ((Vector *)fv + spin + 4 * s.Index())->NormSqNode(6) ;
  }
#ifdef PARALLEL
      slice_sum(&norm, 1, 99);
#endif
  return sqrt(norm) ;
}

FermionVectorTp& FermionVectorTp::operator*=(Float f){
  Site s;
  for(s.Begin();s.End();s.nextSite()){
    for(int spin(0);spin<4;spin++)
      *((Vector *)fv + spin + 4 * s.Index()) *= f;
  }
  
  if(GJP.Gparity()){
    Vector* gfv = (Vector*)fv + GJP.VolNodeSites() * 4;
    for(s.Begin();s.End();s.nextSite()){
      for(int spin(0);spin<4;spin++){
	*(gfv + spin + 4 * s.Index()) *= f;
      }
    }
  }
  return *this ;
}

static void SetVolMomSource1fGparity(Float *fv,int color, int spin, ThreeMom& mom, int flavor) {
  const char *cname = "";
  const char *fname = "SetVolMomSource1fGparity(color,spin,x,y,z,t,flavor)";

  if(GJP.Gparity1fY() && !GJP.Gparity1fX()) ERR.General(cname,fname,"Designed for X and XY, not Y alone");

  //G-parity 1f model.
  int nsz[4] = {GJP.XnodeSites(),GJP.YnodeSites(),GJP.ZnodeSites(),GJP.TnodeSites()};
  int nodes[4] = {GJP.Xnodes(),GJP.Ynodes(),GJP.Znodes(),GJP.Tnodes()};

  static const int quadsgn[4] = {1,1,1,-1}; //(LL,LR,UL,UR)
  int doquad[4] = {0,0,0,0};
  
  if(flavor == 0){
    doquad[0] = 1; if(GJP.Gparity1fX() && GJP.Gparity1fY()) doquad[3] = 1;
  }else{
    doquad[1] = 1; if(GJP.Gparity1fX() && GJP.Gparity1fY()) doquad[2] = 1; 
  }
  Site s ;
  for(s.Begin();s.End();s.nextSite()){
    int quad = 0;
    if(s.physX() >= nsz[0]*nodes[0]/2){
      if(GJP.Gparity1fX() && GJP.Gparity1fY() && s.physY() >= nsz[1]*nodes[1]/2) quad = 3; //UR
      else quad = 1; //LR
    }else{
      if(GJP.Gparity1fX() && GJP.Gparity1fY() && s.physY() >= nsz[1]*nodes[1]/2) quad = 2; //UL
      else quad = 0; //LL
    }
    if(!doquad[quad]) continue;

    ((Complex *)fv)[color + COLORS*(spin + 4*s.Index())] = mom.Fact(s) * (double)quadsgn[quad];
  }
}


// Unit color/spin source at every space-time point with momentum phase
void FermionVectorTp::SetVolMomSource(int color, int spin, ThreeMom& mom, int flavor) {

  const char *fname = "SetVolMomSource()";
  
  VRB.Func(cname, fname);
  
  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
    "Color index out of range: color = %d\n", color);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
    "Spin index out of range: spin = %d\n", spin);

  if(GJP.Gparity1fX() || GJP.Gparity1fY()) return SetVolMomSource1fGparity(fv,color,spin,mom,flavor);

  int cmplx_offset = 0;
  if(GJP.Gparity() && flavor == 1) cmplx_offset = GJP.VolNodeSites() * GJP.Colors() * 4;

  Site s ;
  for(s.Begin();s.End();s.nextSite())
    ((Complex *)fv)[color + COLORS*(spin + 4*s.Index()) + cmplx_offset] = mom.Fact(s); //CK: ThreeMom correctly uses units of pi/2L for momenta in G-parity directions
}



// Momentum wall source
// Check this if it is correct...
void FermionVectorTp::SetMomSource(int color, int spin, int source_time,
				   ThreeMom& mom, int flavor) {
  char *fname = "SetMomSource(color,spin,source_time,Mom)";
  
  VRB.Func(cname, fname);
  
  if (color < 0 || color >= GJP.Colors()) {
    ERR.General(cname, fname,
		"Color index out of range: color = %d\n", color);
  }
  if (spin < 0 || spin > 3) {
    ERR.General(cname, fname,
		"Spin index out of range: spin = %d\n", spin);
  }

  int cmplx_offset = 0;
  if(GJP.Gparity() && flavor == 1) cmplx_offset = GJP.VolNodeSites() * GJP.Colors() * 4; //CK: skip on field

  Site s ;
  for(s.Begin();s.End();s.nextSite())
    if( source_time == s.physT()) // continue only physical time is source_time
      ((Complex *)fv)[color + COLORS*(spin + 4*s.Index()) + cmplx_offset] = mom.Fact(s); //CK: ThreeMom correctly uses units of pi/2L for momenta in G-parity directions
}

// Momentum Cosine wall source
// Check this if it is correct...
void FermionVectorTp::SetMomCosSource(int color, int spin, int source_time,
				      ThreeMom& mom, int flavor) {
  char *fname = "SetMomCosSource(color,spin,source_time,Mom)";
  
  VRB.Func(cname, fname);
  
  if (color < 0 || color >= GJP.Colors()) {
    ERR.General(cname, fname,
		"Color index out of range: color = %d\n", color);
  }
  if (spin < 0 || spin > 3) {
    ERR.General(cname, fname,
		"Spin index out of range: spin = %d\n", spin);
  }
  int cmplx_offset = 0;
  if(GJP.Gparity() && flavor == 1) cmplx_offset = GJP.VolNodeSites() * GJP.Colors() * 4; //CK: skip on field

  Site s ;
  for(s.Begin();s.End();s.nextSite())
    if( source_time == s.physT()) // continue only physical time is source_time
      ((Complex *)fv)[color + COLORS*(spin + 4*s.Index()) + cmplx_offset] = mom.FactCos(s); //CK: ThreeMom correctly uses units of pi/2L for momenta in G-parity directions
}

// Momentum Cosine wall source for twisted boundary conditions
// Momenta are given in units of Pi/L rather than 2*Pi/L
// Check this if it is correct...
void FermionVectorTp::SetMomCosTwistSource(int color, int spin, int source_time, ThreeMomTwist& mom) {
  char *fname = "SetMomCosTwistSource(color,spin,source_time,Mom)";
  
  VRB.Func(cname, fname);

  SetMomCosSource(color,spin,source_time,mom);
}

static void SetPointSource1fGparity(Float *fv,int color, int spin, 
				    int x, int y, int z, int t, int flavor) {
  char *cname = "";
  char *fname = "SetPointSource1fGparity(color,spin,x,y,z,t,flavor)";

  //G-parity 1f model.
  int inpos[4] = {x,y,z,t};
  int sz[4] = {GJP.XnodeSites(),GJP.YnodeSites(),GJP.ZnodeSites(),GJP.TnodeSites()};
  int xy_globsz[2] = {GJP.Xnodes()*sz[0],GJP.Ynodes()*sz[1]};  
  int coor[4] = {GJP.XnodeCoor(),GJP.YnodeCoor(),GJP.ZnodeCoor(),GJP.TnodeCoor()};

  if(x >= xy_globsz[0]/2 || (GJP.Gparity1fY() && y >= xy_globsz[1]/2) ) ERR.General(cname,fname,"1f point source location must be in lower-left quadrant of global X-Y plane");

  int npt = 1; //number of point sources needed
  if(GJP.Gparity1fX() && GJP.Gparity1fY()) npt = 2;

  int ptcoord[npt][4];
  for(int i=0;i<4;i++) ptcoord[0][i] = inpos[i];
  ptcoord[0][0] += flavor*xy_globsz[0]/2; //LL or LR

  Float ptsign[npt];
  ptsign[0] = 1.0;

  if(npt == 2){
    if(flavor == 0) ptsign[1] = -1.0;
    else ptsign[1] = 1.0;

    for(int i=0;i<4;i++) ptcoord[1][i] = inpos[i];
    //UR or UL
    ptcoord[1][0] += (1-flavor)*xy_globsz[0]/2 ;  //R / L
    ptcoord[1][1] += xy_globsz[1]/2; //U
  }

  int procCoor[npt][4];
  int localCoor[npt][4];

  for(int p=0;p<npt;p++){
    for(int i=0;i<4;i++){
      procCoor[p][i] = ptcoord[p][i]/sz[i];
      localCoor[p][i] = ptcoord[p][i] % sz[i];
    }
  }
  
  for(int p=0;p<npt;p++){
    bool ptcoormatch = true;
    for(int i=0;i<4;i++) if(coor[i]!=procCoor[p][i]){ ptcoormatch = false; break; }
    if(ptcoormatch) fv[2*(color + GJP.Colors()*(spin + 4*(localCoor[p][0] + sz[0]*(localCoor[p][1] + sz[1]*(localCoor[p][2] + sz[2]*localCoor[p][3])))))] = ptsign[p];
  }
};

// set point source 
void FermionVectorTp::SetPointSource(int color, int spin, 
				     int x, int y, int z, int t, int flavor) {

  char *fname = "SetPointSource(color,spin,x,y,z,t)";
  
  VRB.Func(cname, fname);

  // trap for wrong arguments
  if (x < 0 || x >= GJP.Xnodes() * GJP.XnodeSites() ||
      y < 0 || y >= GJP.Ynodes() * GJP.YnodeSites() ||
      z < 0 || z >= GJP.Znodes() * GJP.ZnodeSites() ||
      t < 0 || t >= GJP.Tnodes() * GJP.TnodeSites())
    ERR.General(cname, fname,
    "Coordonate arguments out of range: x=%d, y=%d, z=%d, t=%d\n",
    x, y, z, t);

  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
    "Color index out of range: color = %d\n", color);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
    "Spin index out of range: spin = %d\n", spin);

  
  if(GJP.Gparity1fX() || GJP.Gparity1fY()) return SetPointSource1fGparity(fv,color,spin,x,y,z,t,flavor);

  // zero the vector
  //int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  //for (int i = 0; i < fv_size; i++)
  //*((Float *)fv + i) = 0;

  // set point source
  int procCoorX = x / GJP.XnodeSites();
  int procCoorY = y / GJP.YnodeSites();
  int procCoorZ = z / GJP.ZnodeSites();
  int procCoorT = t / GJP.TnodeSites();
  int localX = x % GJP.XnodeSites();
  int localY = y % GJP.YnodeSites();
  int localZ = z % GJP.ZnodeSites();
  int localT = t % GJP.TnodeSites();

  int coor_x = 0;
  int coor_y = 0;
  int coor_z = 0;
  int coor_t = 0;
#ifdef PARALLEL
  coor_x = GJP.XnodeCoor();
  coor_y = GJP.YnodeCoor();
  coor_z = GJP.ZnodeCoor();
  coor_t = GJP.TnodeCoor();
#endif
//VRB.Result("","","HH %d %d %d %d\n", coor_x, coor_y, coor_z, coor_t);
  int offset = 0;
  if(GJP.Gparity() && flavor == 1) offset = GJP.VolNodeSites() * GJP.Colors() * 8; //CK: skip on field

  if (coor_x == procCoorX &&
      coor_y == procCoorY &&
      coor_z == procCoorZ &&
      coor_t == procCoorT){
    fv[2 * (color + GJP.Colors() * (spin + 4 * (
    localX + GJP.XnodeSites() * (
    localY + GJP.YnodeSites() * (
    localZ + GJP.ZnodeSites() * localT)))))
    +offset] = 1.0;

    if(GJP.Gparity() && flavor == 2){
      offset = GJP.VolNodeSites() * GJP.Colors() * 8; //CK: temp hack for placing source on both flavors
      
      fv[2 * (color + GJP.Colors() * (spin + 4 * (
      localX + GJP.XnodeSites() * (
      localY + GJP.YnodeSites() * (
      localZ + GJP.ZnodeSites() * localT)))))
      +offset] = 1.0;
    }
  }

}

// Use the gauge fixing matrix as the source
//  Hardwired for Landau Gauge -- 7 Oct 98 -- MBW
void FermionVectorTp::SetGFPointSource(Lattice& lat,int color, int spin,
									   int x, int y, int z, int t) {

  char *fname = "SetGFPointSource()";
  char *cname = "FermionVectorTp";
  VRB.Func(cname, fname);
 
  Matrix **gm = lat.FixGaugePtr();
 
  // trap for wrong arguments
  if (x < 0 || x >= GJP.Xnodes() * GJP.XnodeSites() ||
      y < 0 || y >= GJP.Ynodes() * GJP.YnodeSites() ||
      z < 0 || z >= GJP.Znodes() * GJP.ZnodeSites() ||
      t < 0 || t >= GJP.Tnodes() * GJP.TnodeSites())
    ERR.General(cname, fname,
    "Coordonate arguments out of range: x=%d, y=%d, z=%d, t=%d\n",
    x, y, z, t);
 
  // trap for color index out of range
  //-------------------------------------------------------------
  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
    "Color index out of range: color = %d\n", color);
 
  // trap for spin index out of range
  //-------------------------------------------------------------
  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
    "Spin index out of range: spin = %d\n", spin);
 
  // zero the vector
  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  for (int i = 0; i < fv_size; i++)
    *((Float *)fv + i) = 0;
 
  // set point source
  int procCoorX = x / GJP.XnodeSites();
  int procCoorY = y / GJP.YnodeSites();
  int procCoorZ = z / GJP.ZnodeSites();
  int procCoorT = t / GJP.TnodeSites();
  int localX = x % GJP.XnodeSites();
  int localY = y % GJP.YnodeSites();
  int localZ = z % GJP.ZnodeSites();
  int localT = t % GJP.TnodeSites();
 
  int coor_x = 0;
  int coor_y = 0;
  int coor_z = 0;
  int coor_t = 0;
#ifdef PARALLEL
  coor_x = GJP.XnodeCoor();
  coor_y = GJP.YnodeCoor();
  coor_z = GJP.ZnodeCoor();
  coor_t = GJP.TnodeCoor();
#endif
//VRB.Result("","","HH %d %d %d %d\n", coor_x, coor_y, coor_z, coor_t);
 
  Matrix *pM = gm[0];
 
  if (coor_x == procCoorX &&
      coor_y == procCoorY &&
      coor_z == procCoorZ &&
      coor_t == procCoorT){
   for(int c=0; c < GJP.Colors(); c++){
    // the real part
    *((Float *)fv + 0 + 2 * (c + GJP.Colors() * (spin + 4 * (
    localX + GJP.XnodeSites() * (
    localY + GJP.YnodeSites() * (
    localZ + GJP.ZnodeSites() * localT)))))) =
      *((Float *)pM + 0 + 2 * (c + GJP.Colors() * (color + GJP.Colors() * (
    localX + GJP.XnodeSites() * (
    localY + GJP.YnodeSites() * (
    localZ + GJP.ZnodeSites() * localT))))));
    // the imaginary part
    *((Float *)fv + 1 + 2 * (c + GJP.Colors() * (spin + 4 * (
    localX + GJP.XnodeSites() * (
    localY + GJP.YnodeSites() * (
    localZ + GJP.ZnodeSites() * localT)))))) =
      - *((Float *)pM + 1 + 2 * (c + GJP.Colors() * (color + GJP.Colors() * (
    localX + GJP.XnodeSites() * (
    localY + GJP.YnodeSites() * (
    localZ + GJP.ZnodeSites() * localT))))));
   } //end for
 } // endif
}

void FermionVectorTp::SetGFLfuncSource(Lattice& lat, int color, int spin,
                Float (*func)(int gx,int gy,int gz,int gt)){
	ERR.NotImplemented(cname,"SetGFLfuncSource()");
}
void FermionVectorTp::SetExpSource(int color, int spin, int x, int y, int z, 
int t, Float A, Float B, Float C){
	ERR.NotImplemented(cname,"SetExpSource()");
}


void FermionVectorTp::CopyWilsonVec(int i, WilsonVector& WV) {
  for(int c=0;c<GJP.Colors();c++)
    for(int s=0;s<4;s++) {
	  int ii = 2*(c+GJP.Colors()*(s+4*i));
	  fv[ii  ] = WV.d[s].c[c].real();
	  fv[ii+1] = WV.d[s].c[c].imag();
	}
}

void FermionVectorTp::CopyWilsonMatSink(int i, int spin, int color, 
										WilsonMatrix& W) {
  for(int c=0;c<GJP.Colors();c++)
    for(int s=0;s<4;s++) {
	  int ii = 2*(c+GJP.Colors()*(s+4*i));
        fv[ii  ] = W.wmat().d[s].c[c].d[spin].c[color].real();
        fv[ii+1] = W.wmat().d[s].c[c].d[spin].c[color].imag();
	}
}


/*! Rotates from Chiral to Dirac basis using the WilsonVector :: ChiralToDirac
  routine
*/
void FermionVectorTp::ChiralToDirac() {
  for(int s=0;s<GJP.VolNodeSites();s++)
    ((WilsonVector*)fv + s)->ChiralToDirac() ;
}

/*! Rotates from Dirac to Chiral basis using the WilsonVector :: DiracToChiral
  routine
*/
void FermionVectorTp::DiracToChiral() {
  for(int s=0;s<GJP.VolNodeSites();s++)
    ((WilsonVector*)fv + s)->DiracToChiral() ;
}

CPS_END_NAMESPACE
