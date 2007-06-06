#include <config.h>
//------------------------------------------------------------------
//
// The class functions for FermionVectorTp.
//
//------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <alg/common_arg.h>
#include <comms/glb.h>
#include <comms/scu.h>
#include <util/site.h>
#include <util/momentum.h>
#include <comms/sysfunc.h>
#include <alg/fermion_vector.h>

CPS_START_NAMESPACE

const Float& FermionVectorTp::operator[](int i) {
  return fv[i];
}


FermionVectorTp::FermionVectorTp() {
  char *fname = "FermionVectorTp()";
  cname = "FermionVectorTp";
  
  VRB.Func(cname, fname);

  // allocate space for source
  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  fv = (Float*)smalloc(fv_size * sizeof(Float));
  if(fv == 0) ERR.Pointer(cname, fname, "fv");
  VRB.Smalloc(cname,fname, "fv", fv, fv_size * sizeof(Float));

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

  int fv_size = GJP.VolNodeSites() * 2 * GJP.Colors() * 4;
  for (int i=0; i<fv_size; i++)
	fv[i] = 0.0;
}

// Unit color/spin source at every space-time point
void FermionVectorTp::SetVolSource(int color, int spin) {

  char *fname = "SetVolSource()";
  
  VRB.Func(cname, fname);
  
  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
    "Color index out of range: color = %d\n", color);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
    "Spin index out of range: spin = %d\n", spin);

  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  for (int i = 0; i < fv_size; i++) {
	if(i%SPINOR_SIZE == 2*(color + COLORS*spin)) {
	  fv[i] = 1.0;
	} else {
	  fv[i] = 0.0;
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
  for (int j = 0; j < fv_size; j++) {
    fv[j] = 0.0;
  }

  for (int i=0; i < GJP.VolNodeSites(); i++) {
	fv[    2*(color + COLORS*(spin + 4*i))] = src[2*i  ];
	fv[1 + 2*(color + COLORS*(spin + 4*i))] = src[2*i+1];
  }
}

// Unit color/spin source at every space point on time_slice
// Does not zero the rest of the  time slices
void FermionVectorTp::SetWallSource(int color, int spin, int source_time) {
  char *fname = "SetWallSource(color,spin,source_time)";
  
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
  int wall_size = GJP.VolNodeSites()/GJP.TnodeSites()*SPINOR_SIZE;
  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  for (int i = 0; i < fv_size; i++) {
    fv[i] = 0.0;
#ifdef PARALLEL
     if(my_node != ts_node) continue;
#endif
     if(i < wall_size*node_ts || i >= wall_size*(node_ts+1)) continue;
     if(i%SPINOR_SIZE != 2*(color + COLORS*spin)) continue;
     fv[i] = 1.0;
  }
}

void FermionVectorTp::SetBoxSource(int color, 
								   int spin,
								   int start,
								   int end,
								   int source_time) {

  char *fname = "SetBoxSource(color,spin,start,end,source_time)";
  
  VRB.Func(cname, fname);
  
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
  int xsize = GJP.XnodeSites();
  int ysize = GJP.YnodeSites();
  int zsize = GJP.ZnodeSites();
  //int tsize = GJP.TnodeSites();
  
  //zero source on all nodes
  for (int i = 0; i < fv_size; i++) {
     fv[i] = 0.0;
  }
  int t = node_ts;
  for(int x = 0; x < xsize; x++){
    for(int y = 0; y < ysize; y++){
      for(int z = 0; z < zsize; z++){
		for(int j = 0; j < GJP.Colors()*8; j++){
		  int i = j + GJP.Colors()*8*(x + xsize*(y+ysize*(z+zsize*(t))));
		  //#ifdef PARALLEL
		  if(tnode != ts_node)continue;
		  else if(x + xnode*xsize < start || x + xnode*xsize > end) continue;
		  else if(y + ynode*ysize < start || y + ynode*ysize > end) continue;
		  else if(z + znode*zsize < start || z + znode*zsize > end) continue;
		  //#endif
		  if(i%SPINOR_SIZE != 2*(color + COLORS*spin))continue;
		  fv[i] = 1.0;
		}
      }
    }
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
void FermionVectorTp::GFWallSource(Lattice &lat, int spin, int dir, int where) {

  char *cname = "FermionVectorT";
  char *fname = "GFWallSource()";
  VRB.Func(cname, fname);

  int len;     //the local (on processor) length in "dir" direction
  //int nproc;   // total number of processors in d_ direction
  int lproc;   // local processor coordinate in d_ direction
               // 0 <= lproc <= nproc
 
  switch(dir) {
    case 0:
      len = GJP.XnodeSites();
      //nproc = GJP.Xnodes();
      lproc = GJP.XnodeCoor();
      break;
    case 1:
      len = GJP.YnodeSites();
      //nproc = GJP.Ynodes();
      lproc = GJP.YnodeCoor();
      break;
    case 2:
      len = GJP.ZnodeSites();
      //nproc = GJP.Znodes();
      lproc = GJP.ZnodeCoor();
      break;
    case 3:
      len = GJP.TnodeSites();
      //nproc = GJP.Tnodes();
      lproc = GJP.TnodeCoor();
      break;
    // trap for wrong direction
    //-----------------------------------------------------------
    default:
      {
		len = lproc = 0 ; // Keep GCC happy!
		ERR.General(cname, fname, "bad argument: dir = %d\n", dir);
      }
  }

  Matrix **gm = lat.FixGaugePtr();

  Vector temp;
  Matrix tempmat;

  // find out if this node overlaps with the hyperplane
  // in which the wall source sits
  int has_overlap = 0;
  if (lproc * len <= where && where < (lproc + 1) * len)
    has_overlap = 1;
 
  if (has_overlap) {
    int local = where % len; // on processor coordinate of
                             // source hyperplane
    Matrix *pM = gm[local];

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
      temp.CopyVec((Vector*)&fv[i], 6);
      tempmat.Dagger((IFloat*)&pM[j]);
      uDotXEqual((IFloat*)&fv[i], (const IFloat*)&tempmat, (const IFloat*)&temp);
    }
  }
}

// COULOMB GAUGE ONLY!
void FermionVectorTp::GaugeFixSink(Lattice &lat, int dir) {
  
  char *fname = "GaugeFixSink()";
  VRB.Func(cname, fname);

  Matrix **gm = lat.FixGaugePtr();

  if(dir!=3)
    ERR.General(cname,fname,"Works only for dir=3\n");

  if(lat.FixGaugeKind()!=FIX_GAUGE_COULOMB_T)
    ERR.General(cname,fname,"Works only for FIX_GAUGE_COULOMB_T\n");

  for (int t=0; t < GJP.TnodeSites(); t++) {

    Matrix *pM = gm[t];
    Vector temp;

    if(gm[t] != NULL ){
      for (int z = 0; z < GJP.ZnodeSites(); z++)
	for (int y = 0; y < GJP.YnodeSites(); y++)
	  for (int x = 0; x < GJP.XnodeSites(); x++)
	    {
	      // the matrix offset
	      int j =  x + GJP.XnodeSites() * ( y + GJP.YnodeSites() * z);
	      for (int spin = 0; spin < 4; spin++)
		{
		  // the vector offset
		  int i= 2 * GJP.Colors() * ( spin + 4 * (
                         x + GJP.XnodeSites() * (
                         y + GJP.YnodeSites() * (
                         z + GJP.ZnodeSites() * t)))) ;
		  temp.CopyVec((Vector*)&fv[i], 6);
		  uDotXEqual((IFloat*)&fv[i],(const IFloat*)&pM[j],
			     (const IFloat*)&temp);
		}
	    }
    }
  }
}

void FermionVectorTp::LandauGaugeFixSink( Lattice& lat ) {

  char *fname = "LandauGaugeFixSink()";
  VRB.Func(cname, fname);
  
  if(lat.FixGaugeKind()!=FIX_GAUGE_LANDAU)
    ERR.General(cname,fname,"lattice not in Landau gauge\n");

  Matrix *pM(lat.FixGaugePtr()[0]);

  Vector temp;
  Site site;

  while (site.LoopsOverNode()) {
	// site offset 
	const int s_off(site.Index());
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

void FermionVectorTp::SetLandauGaugeMomentaSource( Lattice& lat,
                                                   int src_colour,
                                                   int src_spin,
                                                   int p[]) {

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
  
  Float   p1( p[0] );
  Float   p2( p[1] );
  Float   p3( p[2] );
  Float   p4( p[3] );
 
  const Float PI(3.14159265358979323846264338327950288319716939937510);

  p1 *= 2.0*PI/(GJP.XnodeSites()*GJP.Xnodes());
  p2 *= 2.0*PI/(GJP.YnodeSites()*GJP.Ynodes());
  p3 *= 2.0*PI/(GJP.ZnodeSites()*GJP.Znodes());
  p4 *= 2.0*PI/(GJP.TnodeSites()*GJP.Tnodes());
  
  while ( site.LoopsOverNode() )
    {
      // site offset 
      const int s_off( site.Index() );
      Adj.Dagger(pM[s_off]);

      // work out the momentum
      
      const Float px( site.physX()*p1 );
      const Float py( site.physY()*p2 );
      const Float pz( site.physZ()*p3 );
      const Float pt( site.physT()*p4 );
      const Float pdotx( px + py + pz + pt );
      
      const Rcomplex fact( cos(pdotx), sin(pdotx) );
      
      int spin,colour;
      for ( spin=0; spin<4; spin++ )
        {
          const int f_off( 2 * GJP.Colors() * ( spin + 4*s_off) );
          Vector & cvec(*((Vector*)&fv[f_off]));
   
 
    /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik
    -------------------- Quarantine starts --------------------------


// Opened by Sam to do F.T on disconnected prop
*/        

          cvec.Zero();
 
/*

    -------------------- Quarantine ends ---------------------------*/

          for ( colour=0;colour<3;colour++)
            {
              if ( src_spin == spin && src_colour == colour )
                {

    /****************************************************************
     The code below is temporarily isolated for purposes of merging
     with CPS main branch,  12/09/04, Oleg Loktik
    -------------------- Quarantine starts --------------------------

// Opened by Sam to do F.T on disconnected prop
 */

                  cvec[colour] = fact;

/*

    -------------------- Quarantine ends ---------------------------*/
                  
                  temp.CopyVec(&cvec,6);
                  uDotXEqual( (IFloat*)      &fv[f_off],
                              (const IFloat*)&Adj       ,
                              (const IFloat*)&temp       );
                }
              else
                {
                  
                }
            } // colour
        } // spin
    }
}

// Momentum wall source
// Check this if it is correct...
void FermionVectorTp::SetMomSource(int color, int spin, int source_time,
				    ThreeMom& mom) {
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

  Site s ;
  for(s.Begin();s.End();s.nextSite())
    if( source_time == s.physT()) // continue only physical time is source_time
      ((Complex *)fv)[color + COLORS*(spin + 4*s.Index())] = mom.Fact(s);
}


// set point source 
void FermionVectorTp::SetPointSource(int color, int spin, 
				     int x, int y, int z, int t) {

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

  if (coor_x == procCoorX &&
      coor_y == procCoorY &&
      coor_z == procCoorZ &&
      coor_t == procCoorT)
    fv[2 * (color + GJP.Colors() * (spin + 4 * (
    localX + GJP.XnodeSites() * (
    localY + GJP.YnodeSites() * (
    localZ + GJP.ZnodeSites() * localT)))))] = 1.0;

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
