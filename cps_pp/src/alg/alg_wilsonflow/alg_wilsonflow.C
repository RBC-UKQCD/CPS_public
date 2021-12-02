#include<alg/alg_wilsonflow.h>

#include<config.h>
#include<math.h>
#include<util/qcdio.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/error.h>
#include<util/site.h>
#include<util/link_buffer.h>
#include<util/smalloc.h>
#include <comms/scu.h>
#include <util/omp_wrapper.h>


CPS_START_NAMESPACE
//using namespace cps;
//#include<iostream>
//using namespace std;

extern void generateSU3(Matrix &mat, Float Q[8]);
extern const Float SU3_lambda[8][18];

AlgWilsonFlow::AlgWilsonFlow(Lattice &lat, CommonArg *ca, Float dtime, bool proj, Float tol):
	Alg(lat,ca),su3_proj(proj),tolerance(tol),dt(dtime)
{
  cname = "AlgWilsonFlow";

  lat_back = new Matrix[GJP.VolNodeSites()*4];
  if(lat_back==NULL){ERR.Pointer(cname,cname, "lat_back");}
  Z_Lie = new Float[GJP.VolNodeSites()*4*8];
  if(Z_Lie==NULL){ERR.Pointer(cname,cname, "Z_Lie");}
  
  //set up some data used by smartrun:
  Slab = 1; //Expansion in each direction
  MatrixSize = 2 * lat.Colors() * lat.Colors();
  GsiteSize = 4 * MatrixSize;
  
  l_node_sites[0] = GJP.XnodeSites();
  l_node_sites[1] = GJP.YnodeSites();
  l_node_sites[2] = GJP.ZnodeSites();
  l_node_sites[3] = GJP.TnodeSites();
  
  l_dir_offset[0] = GsiteSize;
  l_dir_offset[1] = l_dir_offset[0] * l_node_sites[0];
  l_dir_offset[2] = l_dir_offset[1] * l_node_sites[1];
  l_dir_offset[3] = l_dir_offset[2] * l_node_sites[2];

  vol_node_sites = GJP.VolNodeSites();

  g_node_sites[0] = l_node_sites[0] + 2 * Slab;
  g_node_sites[1] = l_node_sites[1] + 2 * Slab;
  g_node_sites[2] = l_node_sites[2] + 2 * Slab;
  g_node_sites[3] = l_node_sites[3] + 2 * Slab;

  g_dir_offset[0] = GsiteSize;
  g_dir_offset[1] = g_dir_offset[0] * g_node_sites[0];
  g_dir_offset[2] = g_dir_offset[1] * g_node_sites[1];
  g_dir_offset[3] = g_dir_offset[2] * g_node_sites[2];

  g_lcl_vol = g_node_sites[0]
    * g_node_sites[1]
    * g_node_sites[2]
    * g_node_sites[3];
}

AlgWilsonFlow::~AlgWilsonFlow()
{
	delete [] lat_back;
	delete [] Z_Lie;
}

void AlgWilsonFlow::logRun()
{
  const char* fname = "logRun";

  if(common_arg->filename != 0) {
    FILE* f = Fopen(common_arg->filename, "a");
    if(!f) ERR.FileA(cname, fname, common_arg->filename);
    Fprintf(f, "AlgWilsonFlow: dt = %e\n", dt);
    Fclose(f);
  } 
}

void AlgWilsonFlow::run()
{
        logRun();

	Lattice &lattice(AlgLattice());
	Site nloop;

	Float QZ[8];
	Matrix expQ;
	//Follow Luscher's paper: arXiv:1006.4518v1 (June 2010)
	
	//W1 = exp(dt*Z(w0)/4) * W0
	for(nloop.Begin();nloop.End();nloop.nextSite())
	{
		for(int mu=0;mu<4;mu++)
		{
		        //QZ[i] = derivative of S with respect to changes in U_mu in the direction of the generator T^a
			calculateZ(lattice, nloop, mu, QZ);
			for(int i=0;i<8;i++)Z_Lie[8*(4*nloop.Index()+mu)+i] = QZ[i]/4.0*dt;
			//expQ = exp(-dt/4 * dS/dU)   [dS/dU is the SU(3) Lie algebra element whose components in the T^a basis are QZ[a]]
			generateSU3(expQ, &(Z_Lie[8*(4*nloop.Index()+mu)]));
			//new U_mu = expQ * U_mu
			mDotMEqual((Float*)(&lat_back[4*nloop.Index()+mu]), (Float *)&expQ, (Float *)(lattice.GetLink(nloop.pos(),mu)));
		}
	}
	lattice.GaugeField(lat_back);
	lattice.ClearAllBufferedLink();
	//W2 = exp{dt * (Z(W1)*8/9 - Z(W0)*17/36) } W1
	for(nloop.Begin();nloop.End();nloop.nextSite())
	{
		for(int mu=0;mu<4;mu++)
		{
			calculateZ(lattice, nloop, mu, QZ);
			for(int i=0;i<8;i++)Z_Lie[8*(4*nloop.Index()+mu)+i] = QZ[i]*8.0/9.0*dt - 17.0/9.0*Z_Lie[8*(4*nloop.Index()+mu)+i];
			generateSU3(expQ, &(Z_Lie[8*(4*nloop.Index()+mu)]));
			mDotMEqual((Float*)(&lat_back[4*nloop.Index()+mu]), (Float *)&expQ, (Float *)(lattice.GetLink(nloop.pos(),mu)));
		}
	}
	lattice.GaugeField(lat_back);
	lattice.ClearAllBufferedLink();
	//lattice = exp{ dt * (3/4*Z(W2) - 8/9*Z(W1) + 17/36 Z(W0))} W2
	for(nloop.Begin();nloop.End();nloop.nextSite())
	{
		for(int mu=0;mu<4;mu++)
		{
			calculateZ(lattice, nloop, mu, QZ);
			for(int i=0;i<8;i++)Z_Lie[8*(4*nloop.Index()+mu)+i] = QZ[i]*3.0/4.0*dt - Z_Lie[8*(4*nloop.Index()+mu)+i];
			generateSU3(expQ, &(Z_Lie[8*(4*nloop.Index()+mu)]));
			mDotMEqual((Float*)(&lat_back[4*nloop.Index()+mu]), (Float *)&expQ, (Float *)(lattice.GetLink(nloop.pos(),mu)));
		}
	}
	lattice.GaugeField(lat_back);
	lattice.ClearAllBufferedLink();

	if(GJP.Gparity()) lattice.CopyConjGaugeField();
}

void AlgWilsonFlow::calculateZ(Lattice &lat, Site &site, int mu, Float Z[8])
{
	Matrix stap;
	three_staple(lat,stap,site.pos(),mu);
	Matrix stapdag;
	stapdag.Dagger(stap);

    const Matrix *link=lat.GetLink(site.pos(),mu);
	Matrix loop;
	mDotMEqual((Float *)(&loop),(Float*)link, (Float*)(&stapdag));

	Float tmp[18];
	for(int i=0;i<8;i++)
	{
		mDotMEqual(tmp, SU3_lambda[i],(Float*)(&loop));
		Z[i]=-(tmp[1]+tmp[9]+tmp[17]); //-ImTr(tmp)
	}
	
#if 0
	cout<<"Z : ";
	for(int i=0;i<8;i++)cout<<Z[i]<<'\t';
	cout<<endl;
#endif
}

void AlgWilsonFlow::three_staple( Lattice& latt,  Matrix& link , int *pos, int u)
{
	Matrix acumulate_mp; acumulate_mp.ZeroMatrix();
	int dir[3];
	//loop over all directions (+ve and -ve)

	for ( int v=0; v<8; v++ )
	{      
		// except the ones the link is aligned with 
		if((v&3)==u)  continue; 

		const int v1((v+4)&7); // direction opposite to v

		dir[0] = v;
		dir[1] = u;
		dir[2] = v1;

		latt.PathOrdProdPlus(acumulate_mp, pos, dir, 3); 
                char message[128];
                sprintf(message, "three_staple, acumulate_mp, pos = (%d, %d, %d, %d), dir = (%d, %d, %d)", pos[0], pos[1], pos[2], pos[3], dir[0], dir[1], dir[2]);
	}
	moveMem((Float *) &link, (Float*)&acumulate_mp, 18*sizeof(Float));
}



inline Float * AlgWilsonFlow::GsiteOffset(Float * p, const int *x, const int *g_dir_offset)
{
  return p +
    x[0] * g_dir_offset[0] +
    x[1] * g_dir_offset[1] +
    x[2] * g_dir_offset[2] +
    x[3] * g_dir_offset[3];
}

void AlgWilsonFlow::PathOrdProdPlus(Matrix & mat, int* x, int* dirs, int n,
    Float *gfield, const int *g_dir_offset)
{
  //CK: For G-parity we do not have to modify this version as it uses a locally buffered set of gauge links for which
  //    we perform the appropriate complex-conjugation of links crossing the boundary when the buffered field is set up

  int abs_dir;
  int dir_sign;
  int link_site[4];
  const int MatrixSize = 18;

  const Matrix * p1;

  int i;
  for(i = 0; i < 4; ++i)
    link_site[i] = x[i];

  //deal with the first link
  //abs_dir is {0,1,2,3}
  //dir_sign is {0 = positive, 1 = negative}
  //------------------------------
  abs_dir = dirs[0] & 3;
  dir_sign = dirs[0] >> 2;

  //if dir_sign == 1, the link is at x-n_v
  link_site[abs_dir] -= dir_sign;

  //p1 = GetBufferedLink(link_site, abs_dir);
  //get the first link
  p1 = (Matrix*) GsiteOffset(gfield, link_site, g_dir_offset) + abs_dir;

  //if dir_sign == 0, march on to the next site,
  //if dir_sign == 1, we have already moved.
  link_site[abs_dir] += 1 - dir_sign;
  //
  //Temporary matrices and ther pointers
  Matrix ma, mb, mc;
  Matrix *pma  = &ma;
  Matrix *pmb  = &mb;
  Matrix *pmc  = &mc;

  //if dir_sign==1 the link is going backward so get its dagger
  if(dir_sign)
    pma -> Dagger((IFloat*)p1);
  else
    memcpy((IFloat*)pma, (IFloat*)p1, MatrixSize * sizeof(IFloat));

  for(i = 1; i < n; ++i)
  {
    abs_dir = dirs[i]&3;
    dir_sign = dirs[i]>>2;

    link_site[abs_dir] -= dir_sign;

    p1 = (Matrix*) GsiteOffset(gfield, link_site, g_dir_offset) + abs_dir;

    link_site[abs_dir] += 1 - dir_sign;

    //put the next link on the path in mb
    //--------------------------------------
    if(dir_sign)
      mb.Dagger((IFloat*)p1);
    else
      memcpy((IFloat*)pmb, (IFloat*)p1, MatrixSize * sizeof(IFloat));

    // if not the last link on the path, just multiply to the earlier result
    // mc = ma * mb
    // if the last link, multiply and add to mat
    if(i != n - 1)
      mDotMEqual((IFloat*)pmc, (IFloat*)pma, (IFloat*)pmb);
    else
      mDotMPlus((IFloat*)&mat, (IFloat*)pma, (IFloat*)pmb);

    //swap pma and pmc;
    Matrix * tmp_p = pma; pma = pmc; pmc = tmp_p;
  }
}

void AlgWilsonFlow::three_staple(Matrix& link, int *pos, int u, Float * gfield, const int * g_dir_offset)
{
  Matrix acumulate_mp;
  acumulate_mp.ZeroMatrix();
  int dir[3];
  //loop over all directions (+ve and -ve)

  for(int v = 0; v < 8; ++v)
  {
    // except the ones the link is aligned with
    if((v & 3) == u)  continue;

    const int v1((v + 4) & 7); // direction opposite to v

    dir[0] = v;
    dir[1] = u;
    dir[2] = v1;

    PathOrdProdPlus(acumulate_mp, pos, dir, 3, gfield, g_dir_offset);
  }

  // 18 is the matrix size
  memcpy((Float *) &link, (Float*) &acumulate_mp, 18 * sizeof(Float));
}


void AlgWilsonFlow::calculateZ(int pos[4], int mu, Float Z[8], Float * gfield, const int * g_dir_offset)
{
	Matrix stap;
	three_staple(stap, pos, mu, gfield, g_dir_offset);
        Matrix stapdag;
	stapdag.Dagger(stap);

	const Matrix *link = (Matrix *) (GsiteOffset(gfield, pos, g_dir_offset) + mu*18);
        Matrix loop;
	mDotMEqual((Float *)(&loop),(Float*)link, (Float*)(&stapdag));

	Float tmp[18];
	for(int i=0;i<8;i++)
	{
		mDotMEqual(tmp, SU3_lambda[i],(Float*)(&loop));
		Z[i]=-(tmp[1]+tmp[9]+tmp[17]); //-ImTr(tmp)
	}
}


void AlgWilsonFlow::AssembleGfield(Float* lfield, Float* gfield) {
  char fname[] = "AssembleGField";  

  Lattice& lat(AlgLattice());

  int x[4], y[4];

  Float *g_offset;  //offset on the expanded hypercube
  Float *l_offset;  //offset on the orignial hypercube

  int xsta[4] = {0, 0, 0, 0};
  int xend[4] = {l_node_sites[0],
                 l_node_sites[1],
                 l_node_sites[2],
                 l_node_sites[3]};
  int ysta[4] = {Slab, Slab, Slab, Slab};

  int yend[4] = {l_node_sites[0] + Slab,
                 l_node_sites[1] + Slab,
                 l_node_sites[2] + Slab,
                 l_node_sites[3] + Slab};

  memcpy(lfield, (Float *)lat.GaugeField(), GsiteSize * vol_node_sites * sizeof(Float));

  // ------------------------------------------------------------------
  // Move the whole local gauge field to a new location
  // Every dimension is shifted by Slab
  // From L: 0    < x[k] < l_node_sites[k]
  // To   G: Slab < y[k] < l_node_sites[k] + Slab
  // ------------------------------------------------------------------
  x[0] = xsta[0];
  y[0] = ysta[0];
  for(x[1] = xsta[1], y[1] = ysta[1]; x[1] < xend[1]; ++x[1], ++y[1])
    for(x[2] = xsta[2], y[2] = ysta[2]; x[2] < xend[2]; ++x[2], ++y[2])
      for(x[3] = xsta[3], y[3] = ysta[3]; x[3] < xend[3]; ++x[3], ++y[3])
      {
        g_offset = GsiteOffset(gfield, y, g_dir_offset);
        l_offset = GsiteOffset(lfield, x, l_dir_offset);
        memcpy(g_offset, l_offset, l_dir_offset[1] * sizeof(Float));
      }

  Float *surf0;
  Float *surf1;
  int SurfSize;
  int s_dir_offset[4];

  // ------------------------------------------------------------------
  // Propagate the surface slab to neighboring nodes
  // Loop over direction
  // ------------------------------------------------------------------
  for(int i = 0; i < 4; ++i)
  {
    // ------------------------------------------------------------------
    // Set offset on each direction
    // 0 < x[i] < Slab
    // 0 < x[k != i] < l_node_sites[k]
    // ------------------------------------------------------------------
    for(int k = 0; k < 4; ++k)
      s_dir_offset[k] = l_dir_offset[k];
    for(int k = i + 1; k < 4; ++k)
      s_dir_offset[k] = s_dir_offset[k] * Slab / l_node_sites[i];
    xend[i] = Slab;

    // ------------------------------------------------------------------
    // Allocate momory for the data
    // ------------------------------------------------------------------
    SurfSize = vol_node_sites * Slab / l_node_sites[i];
    surf0 = (Float*) smalloc(cname, fname, "surf0", GsiteSize * SurfSize * sizeof(Float));
    surf1 = (Float*) smalloc(cname, fname, "surf1", GsiteSize * SurfSize * sizeof(Float));

    // ------------------------------------------------------------------
    // Assemble the data to propagate in the negative "i" direction
    // From L: 0 < x[i] < Slab
    //         0 < x[k != i] < l_node_sites[k]
    // To   S: 0 < x[i] < Slab
    //         0 < x[k != i] < l_node_sites[k]
    // ------------------------------------------------------------------
    x[0] = xsta[0];
    for(x[1] = xsta[1]; x[1] < xend[1]; ++x[1])
      for(x[2] = xsta[2]; x[2] < xend[2]; ++x[2])
        for(x[3] = xsta[3]; x[3] < xend[3]; ++x[3])
        {
          g_offset = GsiteOffset(surf0, x, s_dir_offset);
          l_offset = GsiteOffset(lfield, x, l_dir_offset);
          memcpy(g_offset, l_offset, s_dir_offset[1] * sizeof(Float));	 
        }
    if(GJP.Bc(i) == BND_CND_GPARITY && GJP.NodeCoor(i)==0) for(int f=1;f<GsiteSize*SurfSize;f+=2) surf0[f]*=-1; //complex conjugate links we are sending across -ve G-parity boundary

    getPlusData(surf1, surf0, GsiteSize * SurfSize, i);

    // ------------------------------------------------------------------
    // Dissemble the received data to the allocated memory
    // From S: 0 < x[i] < Slab
    //         0 < x[k != i] < l_node_sites[k]
    // To   G: l_node_sites[i] + Slab < y[i] < l_node_sites[i] + 2 * Slab
    //         Slab < y[k != i] < l_node_sites[k] + Slab
    // ------------------------------------------------------------------
    ysta[i] = l_node_sites[i] + Slab;
    x[0] = xsta[0];
    y[0] = ysta[0];
    for(x[1] = xsta[1], y[1] = ysta[1]; x[1] < xend[1]; ++x[1], ++y[1])
      for(x[2] = xsta[2], y[2] = ysta[2]; x[2] < xend[2]; ++x[2], ++y[2])
        for(x[3] = xsta[3], y[3] = ysta[3]; x[3] < xend[3]; ++x[3], ++y[3])
        {
          l_offset = GsiteOffset(surf1, x, s_dir_offset);
          g_offset = GsiteOffset(gfield, y, g_dir_offset);
          memcpy(g_offset, l_offset, s_dir_offset[1] * sizeof(Float));
        }
    ysta[i] = Slab;

    // ------------------------------------------------------------------
    // If l_node_site[i] == Slab,  we only need to pass the same slab
    // to the opposite direcction
    // Otherwise we need to assemble the appropriate slab again
    // From L: l_node_sites[i] - Slab < x[i] < l_node_sites[i]
    //         0 < x[k != i] < l_node_sites[k]
    // To   S: 0 < x[i] < Slab
    //         0 < x[k != i] < l_node_sites[k]
    // ------------------------------------------------------------------
    if(l_node_sites[i] != Slab)
    {
      ysta[i] = l_node_sites[i];
      x[0] = xsta[0];
      y[0] = ysta[0];
      for(x[1] = xsta[1], y[1] = ysta[1]; x[1] < xend[1]; ++x[1], ++y[1])
        for(x[2] = xsta[2], y[2] = ysta[2]; x[2] < xend[2]; ++x[2], ++y[2])
          for(x[3] = xsta[3], y[3] = ysta[3]; x[3] < xend[3]; ++x[3], ++y[3])
          {
            g_offset = GsiteOffset(surf0, x, s_dir_offset);
            l_offset = GsiteOffset(gfield, y, g_dir_offset);
            memcpy(g_offset, l_offset, s_dir_offset[1] * sizeof(Float));
          }
      ysta[i] = Slab;
    }
    if(GJP.Bc(i) == BND_CND_GPARITY && GJP.NodeCoor(i)==GJP.Nodes(i)-1) for(int f=1;f<GsiteSize*SurfSize;f+=2) surf0[f]*=-1; //complex conjugate links we are sending across +ve G-parity boundary

    getMinusData(surf1, surf0, GsiteSize * SurfSize, i);

    // ------------------------------------------------------------------
    // Dissemble the received data to the allocated memory
    // From S: 0 < x[i] < Slab
    //         0 < x[k != i] < l_node_sites[k]
    // To   G: 0 < y[i] < Slab
    //         Slab < y[k != i] < l_node_sites[k] + Slab
    // ------------------------------------------------------------------
    ysta[i] = 0;
    x[0] = xsta[0];
    y[0] = ysta[0];
    for(x[1] = xsta[1], y[1] = ysta[1]; x[1] < xend[1]; ++x[1], ++y[1])
      for(x[2] = xsta[2], y[2] = ysta[2]; x[2] < xend[2]; ++x[2], ++y[2])
        for(x[3] = xsta[3], y[3] = ysta[3]; x[3] < xend[3]; ++x[3], ++y[3])
        {
          l_offset = GsiteOffset(surf1, x, s_dir_offset);
          g_offset = GsiteOffset(gfield, y, g_dir_offset);
          memcpy(g_offset, l_offset, s_dir_offset[1] * sizeof(Float));
        }
    ysta[i] = Slab;

    sfree(cname, fname, "surf0", surf0);
    sfree(cname, fname, "surf1", surf1);

    // ------------------------------------------------------------------
    // Propagate the cornered chunk to neighboring nodes
    // Loop over direction
    // ------------------------------------------------------------------
    for(int j = i + 1; j < 4; ++j)
    {
      // ------------------------------------------------------------------
      // 0 < x[i, j] < Slab
      // 0 < x[k != i, j] < l_node_sites[k]
      // ------------------------------------------------------------------
      for(int k = i + 1; k < 4; ++k)
        s_dir_offset[k] = l_dir_offset[k] * Slab / l_node_sites[i];
      for(int k = j + 1; k < 4; ++k)
        s_dir_offset[k] = s_dir_offset[k] * Slab / l_node_sites[j];

      xend[j] = Slab;

      Float * surf2;

      // ------------------------------------------------------------------
      // Every chunk contains two corner pieces
      // ------------------------------------------------------------------
      SurfSize = vol_node_sites * Slab  * Slab / l_node_sites[i] / l_node_sites[j];
      surf0 = (Float*) smalloc(cname, fname, "surf0", 2 * GsiteSize * SurfSize * sizeof(Float));
      surf1 = (Float*) smalloc(cname, fname, "surf1", 2 * GsiteSize * SurfSize * sizeof(Float));

      // ------------------------------------------------------------------
      // Assemble the data to propagate in the negative "j" direction
      // From G: 0 < y[i] < Slab
      //         Slab < y[j] <  2 * Slab
      //         Slab < y[k != i, j] < l_node_sites[k] + Slab
      // To   S: 0 < x[i, j] < Slab
      //         0 < x[k != i,j] < l_node_sites[k]
      // Also:
      // From G: l_node_sites[i] + Slab < y[i] < l_node_sites[i] + 2 * Slab
      //         Slab < y[j] <  2 * Slab
      //         Slab < y[k != i, j] < l_node_sites[k] + Slab
      // To   S: 0 < x[i, j] < Slab
      //         0 < x[k != i,j] < l_node_sites[k]
      // ------------------------------------------------------------------
      surf2 = surf0 + GsiteSize * SurfSize;
      ysta[i] = 0;
      x[0] = xsta[0];
      y[0] = ysta[0];
      for(x[1] = xsta[1], y[1] = ysta[1]; x[1] < xend[1]; ++x[1], ++y[1])
        for(x[2] = xsta[2], y[2] = ysta[2]; x[2] < xend[2]; ++x[2], ++y[2])
          for(x[3] = xsta[3], y[3] = ysta[3]; x[3] < xend[3]; ++x[3], ++y[3])
          {
            l_offset = GsiteOffset(surf0, x, s_dir_offset);
            g_offset = GsiteOffset(gfield, y, g_dir_offset);
            memcpy(l_offset, g_offset, s_dir_offset[1] * sizeof(Float));
            y[i] += l_node_sites[i] + Slab;
            l_offset = GsiteOffset(surf2, x, s_dir_offset);
            g_offset = GsiteOffset(gfield, y, g_dir_offset);
            memcpy(l_offset, g_offset, s_dir_offset[1] * sizeof(Float));
            y[i] -= l_node_sites[i] + Slab;
          }
      ysta[i] = Slab;

      if(GJP.Bc(j) == BND_CND_GPARITY && GJP.NodeCoor(j)==0) for(int f=1;f<2*GsiteSize*SurfSize;f+=2) surf0[f]*=-1; //complex conjugate links we are sending across -ve G-parity boundary

      getPlusData(surf1, surf0, 2 * GsiteSize * SurfSize, j);

      // ------------------------------------------------------------------
      // Dissemble the received data to the allocated memory
      // From S: 0 < x[i, j] < Slab
      //         0 < x[k != i,j] < l_node_sites[k]
      // To   G: 0 < y[i] < Slab
      //         l_node_sites[j] + Slab < y[j] < l_node_sites[j] + 2 * Slab
      //         Slab < y[k != i, j] < l_node_sites[k] + Slab
      // Also:
      // From S: 0 < x[i, j] < Slab
      //         0 < x[k != i,j] < l_node_sites[k]
      // To   G: l_node_sites[i] + Slab < y[i] < l_node_sites[i] + 2 * Slab
      //         l_node_sites[j] + Slab < y[j] < l_node_sites[j] + 2 * Slab
      //         Slab < y[k != i, j] < l_node_sites[k] + Slab
      // ------------------------------------------------------------------
      surf2 = surf1 + GsiteSize * SurfSize;
      ysta[i] = 0;
      ysta[j] = l_node_sites[j] + Slab;
      x[0] = xsta[0];
      y[0] = ysta[0];
      for(x[1] = xsta[1], y[1] = ysta[1]; x[1] < xend[1]; ++x[1], ++y[1])
        for(x[2] = xsta[2], y[2] = ysta[2]; x[2] < xend[2]; ++x[2], ++y[2])
          for(x[3] = xsta[3], y[3] = ysta[3]; x[3] < xend[3]; ++x[3], ++y[3])
          {
            l_offset = GsiteOffset(surf1, x, s_dir_offset);
            g_offset = GsiteOffset(gfield, y, g_dir_offset);
            memcpy(g_offset, l_offset, s_dir_offset[1] * sizeof(Float));
            y[i] += l_node_sites[i] + Slab;
            l_offset = GsiteOffset(surf2, x, s_dir_offset);
            g_offset = GsiteOffset(gfield, y, g_dir_offset);
            memcpy(g_offset, l_offset, s_dir_offset[1] * sizeof(Float));
            y[i] -= l_node_sites[i] + Slab;
          }
      ysta[i] = Slab;
      ysta[j] = Slab;

      // ------------------------------------------------------------------
      // If l_node_site[j] == Slab,  We only need to pass the same slab
      // to the opposite direcction
      // Otherwise we need to assemble the appropriate slab
      // Assemble the data to propagate in the positive "j" direction
      // From G: 0 < y[i] < Slab
      //         Slab < y[j] <  2 * Slab
      //         l_node_sites[j] < y[j] < l_node_sites[j] + Slab
      //         Slab < y[k != i, j] < l_node_sites[k] + Slab
      // To   S: 0 < x[i, j] < Slab
      //         0 < x[k != i,j] < l_node_sites[k]
      // Also:
      // From G: l_node_sites[i] + Slab < y[i] < l_node_sites[i] + 2 * Slab
      //         l_node_sites[j] < y[j] < l_node_sites[j] + Slab
      //         Slab < y[k != i, j] < l_node_sites[k] + Slab
      // To   S: 0 < x[i, j] < Slab
      //         0 < x[k != i,j] < l_node_sites[k]
      // ------------------------------------------------------------------
      if(l_node_sites[j] != Slab)
      {
        surf2 = surf0 + GsiteSize * SurfSize;
        ysta[i] = 0;
        ysta[j] = l_node_sites[j];
        x[0] = xsta[0];
        y[0] = ysta[0];
        for(x[1] = xsta[1], y[1] = ysta[1]; x[1] < xend[1]; ++x[1], ++y[1])
          for(x[2] = xsta[2], y[2] = ysta[2]; x[2] < xend[2]; ++x[2], ++y[2])
            for(x[3] = xsta[3], y[3] = ysta[3]; x[3] < xend[3]; ++x[3], ++y[3])
            {
              l_offset = GsiteOffset(surf0, x, s_dir_offset);
              g_offset = GsiteOffset(gfield, y, g_dir_offset);
              memcpy(l_offset, g_offset, s_dir_offset[1] * sizeof(Float));
              y[i] += l_node_sites[i] + Slab;
              l_offset = GsiteOffset(surf2, x, s_dir_offset);
              g_offset = GsiteOffset(gfield, y, g_dir_offset);
              memcpy(l_offset, g_offset, s_dir_offset[1] * sizeof(Float));
              y[i] -= l_node_sites[i] + Slab;
            }
        ysta[i] = Slab;
        ysta[j] = Slab;
      }
      if(GJP.Bc(j) == BND_CND_GPARITY && GJP.NodeCoor(j)==GJP.Nodes(j)-1) for(int f=1;f<2*GsiteSize*SurfSize;f+=2) surf0[f]*=-1; //complex conjugate links we are sending across +ve G-parity boundary

      getMinusData(surf1, surf0, 2 * GsiteSize * SurfSize, j);

      // ------------------------------------------------------------------
      // Dissemble the received data to the allocated memory
      // From S: 0 < x[i, j] < Slab
      //         0 < x[k != i,j] < l_node_sites[k]
      // To   G: 0 < y[i] < Slab
      //         0 < y[j] < Slab
      //         Slab < y[k != i, j] < l_node_sites[k] + Slab
      // Also:
      // From S: 0 < x[i, j] < Slab
      //         0 < x[k != i,j] < l_node_sites[k]
      // To   G: l_node_sites[i] + Slab < y[i] < l_node_sites[i] + 2 * Slab
      //         0 < y[j] < Slab
      //         Slab < y[k != i, j] < l_node_sites[k] + Slab
      // ------------------------------------------------------------------
      surf2 = surf1 + GsiteSize * SurfSize;
      ysta[i] = 0;
      ysta[j] = 0;
      x[0] = xsta[0];
      y[0] = ysta[0];
      for(x[1] = xsta[1], y[1] = ysta[1]; x[1] < xend[1]; ++x[1], ++y[1])
        for(x[2] = xsta[2], y[2] = ysta[2]; x[2] < xend[2]; ++x[2], ++y[2])
          for(x[3] = xsta[3], y[3] = ysta[3]; x[3] < xend[3]; ++x[3], ++y[3])
          {
            l_offset = GsiteOffset(surf1, x, s_dir_offset);
            g_offset = GsiteOffset(gfield, y, g_dir_offset);
            memcpy(g_offset, l_offset, s_dir_offset[1] * sizeof(Float));
            y[i] += l_node_sites[i] + Slab;
            l_offset = GsiteOffset(surf2, x, s_dir_offset);
            g_offset = GsiteOffset(gfield, y, g_dir_offset);
            memcpy(g_offset, l_offset, s_dir_offset[1] * sizeof(Float));
            y[i] -= l_node_sites[i] + Slab;
          }
      ysta[i] = Slab;
      ysta[j] = Slab;

      // ------------------------------------------------------------------
      // Release memory and  set xend[j] = l_node_sites[j]
      // ------------------------------------------------------------------
      sfree(cname, fname, "surf0", surf0);
      sfree(cname, fname, "surf1", surf1);
      xend[j] = l_node_sites[j];
    }
    // ------------------------------------------------------------------
    // Before go on to the next direction, set xend[i] = l_node_sites[i]
    // ------------------------------------------------------------------
    xend[i] = l_node_sites[i];
  }
}


void AlgWilsonFlow::DoRK4Step(int rk4_step, int site, Float* lfield, int l_dir_offset[4], Float* gfield, int g_dir_offset[4]) 
{
  const char* fname = "DoRK4Step";
  //Follow Luscher's paper: arXiv:1006.4518v1 (June 2010)
  
  int Slab = 1;

  int j = site;
  int x[4] = {0};
  x[0] = j % GJP.XnodeSites(); j /= GJP.XnodeSites();
  x[1] = j % GJP.YnodeSites(); j /= GJP.YnodeSites();
  x[2] = j % GJP.ZnodeSites(); j /= GJP.ZnodeSites();
  x[3] = j;
  int y[4] = {x[0] + Slab,
    x[1] + Slab,
    x[2] + Slab,
    x[3] + Slab};
	
  Float QZ[8];
  Matrix expQ;
      
  for(int mu = 0; mu < 4; mu++) 
  {
    //QZ[a] = derivative of S with respect to changes in U_mu in the direction of the generator T^a
    calculateZ(y, mu, QZ, gfield, g_dir_offset);      

    for(int i = 0; i < 8; i++) {
      if(QZ[i] != QZ[i]) ERR.General(cname, fname, "QZ[%d] is null at mu = %d, rk4_step = %d, site = %d, x = (%d, %d, %d, %d), y = (%d, %d, %d, %d)", i, mu, rk4_step, site, x[0], x[1], x[2], x[3], y[0], y[1], y[2], y[3]);
    }
        
    switch(rk4_step) 
    {
      case 0:
        for(int i = 0; i < 8; i++) Z_Lie[8*(4*site+mu)+i] = QZ[i]/4.0 * dt;
        break;

      case 1:
        for(int i = 0; i < 8; i++) Z_Lie[8*(4*site+mu)+i] = QZ[i]*8.0/9.0*dt - 17.0/9.0*Z_Lie[8*(4*site+mu)+i];
        break;

      case 2:
	for(int i = 0; i < 8; i++) Z_Lie[8*(4*site+mu)+i] = QZ[i]*3.0/4.0*dt - Z_Lie[8*(4*site+mu)+i];
        break;
    }	

    generateSU3(expQ, &(Z_Lie[8*(4*site+mu)]));
      
    //new U_mu = expQ * U_mu
    mDotMEqual(GsiteOffset(lfield, x, l_dir_offset) + mu*18, (Float *)&expQ, GsiteOffset(gfield, y, g_dir_offset) + mu*18);
  }
    
}

void AlgWilsonFlow::smartrun()
{
  logRun();

  const char fname[] = "smartrun";

  Lattice& lat(AlgLattice());
  
  Float * lfield = (Float*) smalloc(cname, fname, "lfield", GsiteSize * vol_node_sites * sizeof(Float));

  Float * gfield = (Float*) smalloc(cname, fname, "gfield", GsiteSize * g_lcl_vol * sizeof(Float)); 

  for(int rk4_step = 0; rk4_step < 3; rk4_step++) 
  {
    //Grab all the needed links from neighboring nodes
    AssembleGfield(lfield, gfield); 

    //Do each of the four steps of the RK4 integrator in turn:
    omp_set_num_threads(GJP.Nthreads());

    #pragma omp parallel for
    for(int site = 0; site < GJP.VolNodeSites(); ++site)
    {
      DoRK4Step(rk4_step, site, lfield, l_dir_offset, gfield, g_dir_offset);
    }

    memcpy((Float *)lat.GaugeField(), lfield, GsiteSize * vol_node_sites * sizeof(Float));
  }
  if(GJP.Gparity()) lat.CopyConjGaugeField();

  sfree(cname, fname, "gfield", gfield);
  sfree(cname, fname, "lfield", lfield);
}



CPS_END_NAMESPACE

