#include <chroma.h>
#include <actions/ferm/invert/syssolver_linop_cg_array.h>
#include <bfm_qdp.h>
//#include <qmp.h>

typedef bfm_qdp<float> bfm_t;
//typedef commQMPbagel1<double> bfm_t;

using namespace Chroma;

typedef LatticeFermion T;
typedef multi1d<LatticeFermion> T5;
typedef multi1d<LatticeColorMatrix> U;

//#define Printf if ( QMP_is_primary_node() ) printf
#define Printf printf
#define NITER 1

// //put fermion flavs onto doubled lattice
// void imp_ferm_doublelatt(QDPdouble* to, QDPdouble *f0_from, QDPdouble *f1_from, int lx, int ly, int lz, int lt){
//   //lattice sizes are of those of the doubled lattice
//   int x[4];

//   for(x[0]=0;x[0]<lx/2;x[0]++){
//     for(x[1]=0;x[1]<ly;x[1]++){
//       for(x[2]=0;x[2]<lz;x[2]++){
// 	for(x[3]=0;x[3]<lt;x[3]++){
	  
// 	  for ( int spco=0;spco<12;spco++ ) { 
// 	    for ( int reim=0;reim<2;reim++ ) { 
// 	      int orig_idx;
// 	      {//calculate ptr offset in original, undoubled field
// 		int ccb   = ((x[0]+x[1]+x[2]+x[3])&0x1);
// 		int csite= x[0] + lx/2*(x[1]+ly*(x[2]+lz*x[3]));
// 		csite = csite/2;
// 		int cbvol = (lx/2*ly*lz*lt)/2;
// 		orig_idx = (ccb*cbvol+csite)*12*2 + spco*2 + reim;
// 	      }
// 	      int doubled_idx_f0;
// 	      {//calculate ptr offset doubled field
// 		int ccb   = ((x[0]+x[1]+x[2]+x[3])&0x1);
// 		int csite= x[0] + lx*(x[1]+ly*(x[2]+lz*x[3]));
// 		csite = csite/2;
// 		int cbvol = (lx*ly*lz*lt)/2;
// 		doubled_idx_f0 = (ccb*cbvol+csite)*12*2 + spco*2 + reim;
// 	      }
// 	      int doubled_idx_f1;
// 	      {//calculate ptr offset doubled field
// 		int ccb   = ((x[0]+lx/2+x[1]+x[2]+x[3])&0x1);
// 		int csite= x[0] + lx/2 + lx*(x[1] + ly*(x[2]+lz*x[3]));
// 		csite = csite/2;
// 		int cbvol = (lx*ly*lz*lt)/2;
// 		doubled_idx_f1 = (ccb*cbvol+csite)*12*2 + spco*2 + reim;
// 	      }
// 	      to[doubled_idx_f0] = f0_from[orig_idx];
// 	      to[doubled_idx_f1] = f1_from[orig_idx];
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
  
// }

  //Before doubling each node contains lx*ly*lz*lt sites
  //After doubling each node contains 2lx*ly*lz*lt sites
  //We need to make the second 4-volume accessible to each node for the doubling

    /*
      | 0 1         | 2 3         | 4 5         | 6 7         |     ----> x    original layout
      | {0 1} {2 3} | {4 5} {6 7} | {0 1} {2 3} | {4 5} {6 7} |   curly brackets are the data in the 2 blocks
    */

void create_blocks_gauge(multi1d<multi1d<LatticeColorMatrix> > &u_xblock, multi1d<LatticeColorMatrix> &u_orig){
  //u_xblock is 2 x Nd
  multi1d<int> ncoor = QDP::Layout::nodeCoord();
  multi1d<int> nodes = QDP::Layout::logicalSize();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  int nrow[] = {nodes[0]*lx,nodes[1]*ly,nodes[2]*lz,nodes[3]*lt};

  for(int mu = 0;mu<Nd;mu++){
    QDPIO::cout << "Mu = " << mu << std::endl;
    //loop over global coords
    multi1d<int> xx(4);
    for(xx[0] = 0; xx[0] < nrow[0]; xx[0]++){
      QDPIO::cout << "x = " << xx[0] << std::endl;
      //what xnode is this?
      int xnode = xx[0]/lx;
      //how far is this into the block
      int xinto_block = xx[0] % lx;
      //what x coord do block0 and block1 need?
      int x0 = (xnode*2*lx + xinto_block) % (nodes[0]*lx);
      int x1 = (xnode*2*lx + xinto_block + lx) % (nodes[0]*lx);

      // if(mu==0){
      // 	printf("x=%d in block0->%d block1->%d. Node is %d and distance into node is %d. Pre-wrap %d %d, wrap amount %d\n",xx[0],x0,x1,xnode,xinto_block,xnode*2*lx + xinto_block,xnode*2*lx + xinto_block + lx,nodes[0]*lx);
      // }

      for(xx[1] = 0; xx[1] < nrow[1]; xx[1]++){
	for(xx[2] = 0; xx[2] < nrow[2]; xx[2]++){
	  for(xx[3] = 0; xx[3] < nrow[3]; xx[3]++){
	    QDPIO::cout << "y,z,t = " << xx[1] << "," << xx[2]<< "," << xx[3] << std::endl;

	    multi1d<int> yy(xx); 
	    yy[0] = x0;
	    ColorMatrix cm0 = peekSite(u_orig[mu],yy);
	    yy[0] = x1;
	    ColorMatrix cm1 = peekSite(u_orig[mu],yy);

	    pokeSite(u_xblock[0][mu],cm0,xx);
	    pokeSite(u_xblock[1][mu],cm1,xx);

	  }
	}
      }
    }
  }

}

void create_blocks_ferm(multi1d<multi1d<LatticeFermion> > &f_xblock, multi1d<LatticeFermion> &f0_orig, multi1d<LatticeFermion> &f1_orig, int Ls){
  //f_xblock is 2 x Ls
  multi1d<int> ncoor = QDP::Layout::nodeCoord();
  multi1d<int> nodes = QDP::Layout::logicalSize();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  int nrow[] = {nodes[0]*lx,nodes[1]*ly,nodes[2]*lz,nodes[3]*lt};

  for(int s = 0;s<Ls;s++){
    //loop over global coords
    multi1d<int> xx(4);
    for(xx[0] = 0; xx[0] < nrow[0]; xx[0]++){
      //what xnode is this?
      int xnode = xx[0]/lx;
      //how far is this into the block
      int xinto_block = xx[0] % lx;
      //what x coord do block0 and block1 need?
      int x0 = (xnode*2*lx + xinto_block) % (nodes[0]*lx);
      int x1 = (xnode*2*lx + xinto_block + lx) % (nodes[0]*lx);

      LatticeFermion* block0_from = & f0_orig[s];
      LatticeFermion* block1_from = & f0_orig[s];
      if(xnode*2*lx + xinto_block >= nodes[0]*lx){
	//if(s==0) printf("block0 target x-coord %d, pull from %d, f1\n",xnode*2*lx + xinto_block,x0);
	block0_from = & f1_orig[s];
      }//else if(s==0) printf("block0 target x-coord %d, pull from %d, f0\n",xnode*2*lx + xinto_block,x0);

      if(xnode*2*lx + xinto_block + lx >= nodes[0]*lx){
	//if(s==0) printf("block1 target x-coord %d, pull from %d, f1\n",xnode*2*lx + xinto_block +lx,x1);
	block1_from = & f1_orig[s];
      }//else if(s==0) printf("block1 target x-coord %d, pull from %d, f0\n",xnode*2*lx + xinto_block +lx,x0);

      for(xx[1] = 0; xx[1] < nrow[1]; xx[1]++){
	for(xx[2] = 0; xx[2] < nrow[2]; xx[2]++){
	  for(xx[3] = 0; xx[3] < nrow[3]; xx[3]++){

	    multi1d<int> yy(xx); 
	    yy[0] = x0;
	    Fermion cm0 = peekSite(*block0_from,yy);
	    yy[0] = x1;
	    Fermion cm1 = peekSite(*block1_from,yy);

	    pokeSite(f_xblock[0][s],cm0,xx);
	    pokeSite(f_xblock[1][s],cm1,xx);

	  }
	}
      }
    }
  }

}




//blocks each contain lx space-time slices that are to be imported onto the doubled lattice
void import_blocks_to_double_gauge(QDPdouble* doubled, QDPdouble* block0, QDPdouble* block1){
  int lx = QDP::Layout::subgridLattSize()[0];  //after lattice has been doubled, this lx = 2* original lx
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  int x[4]; 
  for(x[0]=0;x[0]<lx;x[0]++){
    for(x[1]=0;x[1]<ly;x[1]++){
      for(x[2]=0;x[2]<lz;x[2]++){
	for(x[3]=0;x[3]<lt;x[3]++){
	  QDPdouble* block = block0;
	  int x_orig = x[0];
	  if(x_orig >= lx/2){
	    x_orig -= lx/2;
	    block = block1;
	  }

	  for ( int coco=0;coco<9;coco++ ) { 
	    for ( int reim=0;reim<2;reim++ ) { 
	      int orig_idx;
	      {//calculate ptr offset in original, undoubled field
		int ccb   = ((x_orig+x[1]+x[2]+x[3])&0x1);
		int csite= x_orig + lx/2*(x[1]+ly*(x[2]+lz*x[3]));
		csite = csite/2;
		int cbvol = (lx/2*ly*lz*lt)/2;
		orig_idx = (ccb*cbvol+csite)*9*2 + coco*2 + reim;
	      }
	      int doubled_idx;
	      {//calculate ptr offset doubled field
		int ccb   = ((x[0]+x[1]+x[2]+x[3])&0x1);
		int csite= x[0] + lx*(x[1]+ly*(x[2]+lz*x[3]));
		csite = csite/2;
		int cbvol = (lx*ly*lz*lt)/2;
		doubled_idx = (ccb*cbvol+csite)*9*2 + coco*2 + reim;
	      }
	      doubled[doubled_idx] = block[orig_idx];
	    }
	  }
	}
      }
    }
  }
}

//blocks each contain lx space-time slices that are to be imported onto the doubled lattice
void import_blocks_to_double_ferm(QDPdouble* doubled, QDPdouble* block0, QDPdouble* block1){
  int lx = QDP::Layout::subgridLattSize()[0]; //after lattice has been doubled, this lx = 2* original lx
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  int x[4]; 
  for(x[0]=0;x[0]<lx;x[0]++){
    for(x[1]=0;x[1]<ly;x[1]++){
      for(x[2]=0;x[2]<lz;x[2]++){
	for(x[3]=0;x[3]<lt;x[3]++){
	  QDPdouble* block = block0;
	  int x_orig = x[0];
	  if(x_orig >= lx/2){
	    x_orig -= lx/2;
	    block = block1;
	  }

	  for ( int spco=0;spco<12;spco++ ) { 
	    for ( int reim=0;reim<2;reim++ ) { 
	      int orig_idx;
	      {//calculate ptr offset in original, undoubled field
		int ccb   = ((x_orig+x[1]+x[2]+x[3])&0x1);
		int csite= x_orig + lx/2*(x[1]+ly*(x[2]+lz*x[3]));
		csite = csite/2;
		int cbvol = (lx/2*ly*lz*lt)/2;
		orig_idx = (ccb*cbvol+csite)*12*2 + spco*2 + reim;
	      }
	      int doubled_idx;
	      {//calculate ptr offset doubled field
		int ccb   = ((x[0]+x[1]+x[2]+x[3])&0x1);
		int csite= x[0] + lx*(x[1]+ly*(x[2]+lz*x[3]));
		csite = csite/2;
		int cbvol = (lx*ly*lz*lt)/2;
		doubled_idx = (ccb*cbvol+csite)*12*2 + spco*2 + reim;
	      }
	      doubled[doubled_idx] = block[orig_idx];
	    }
	  }
	}
      }
    }
  }
}


int main (int argc,char **argv )
{
  Chroma::initialize(&argc,&argv);


  /********************************************************
   * Command line parsing
   ********************************************************
   */
#define COMMANDLINE
#ifdef COMMANDLINE
  if ( argc != 6 ) { 
    QDPIO::cout << "Usage: " << argv[0] << " lx ly lz lt Ls\n All must be even\n";
    QDPIO::cout << "argc is " << argc << "\n";
   for ( int i=0;i<argc;i++){
     QDPIO::cout << i << " " << argv[i] << "\n";
   }
   exit(-1);

  }
#endif
  /********************************************************
   * Setup QDP
   ********************************************************
   */
  multi1d<int> nrow(Nd);
#ifdef COMMANDLINE
  nrow[0] = atoi(argv[1]);
  nrow[1] = atoi(argv[2]);
  nrow[2] = atoi(argv[3]);
  nrow[3] = atoi(argv[4]);
  int Ls = atoi(argv[5]);
#else
  nrow[0] = 4;
  nrow[1] = 4;
  nrow[2] = 4;
  nrow[3] = 4;
  int Ls  = 4;
#endif

  Layout::setLattSize(nrow);
  Layout::create();

  int lx = QDP::Layout::subgridLattSize()[0];
  int ly = QDP::Layout::subgridLattSize()[1];
  int lz = QDP::Layout::subgridLattSize()[2];
  int lt = QDP::Layout::subgridLattSize()[3];

  /********************************************************
   * Setup DWF operator
   ********************************************************
   */
  bfmarg dwfa;

  dwfa.solver = DWF;
  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;
  dwfa.verbose=1;
  dwfa.reproduce=0;

  int gp[] = {1,0,0};
  dwfa.EnableGparity(gp);

  multi1d<int> procs = QDP::Layout::logicalSize();
  QDPIO::cout << procs.size() << " dim machine\n\t";
  for(int mu=0;mu<4;mu++){
    QDPIO::cout << procs[mu] << " ";
    if ( procs[mu]>1 ) {
      dwfa.local_comm[mu] = 0;
    } else { 
      dwfa.local_comm[mu] = 1;
    }
  }
  QDPIO::cout << "\nLocal comm = ";
  for(int mu=0;mu<4;mu++){
    QDPIO::cout << dwfa.local_comm[mu] << " ";
  }
  QDPIO::cout << "\n";
  
  multi1d<int> ncoor = QDP::Layout::nodeCoord();
  multi1d<int> nodes = QDP::Layout::logicalSize();
  for(int i=0;i<4;i++){ dwfa.ncoor[i] = ncoor[i]; dwfa.nodes[i] = nodes[i]; }
  QDPIO::cout << "\nNode coor: " << ncoor[0] << " " << ncoor[1] << " " << ncoor[2] << " " << ncoor[3];
  QDPIO::cout << "\nGrid size: " << nodes[0] << " " << nodes[1] << " " << nodes[2] << " " << nodes[3] << "\n";

  Real M5(1.8);
  Real mq(0.1);
  
  bfmarg::Threads(64);
  bfmarg::Verbose(1);

  dwfa.precon_5d = 1;
  dwfa.Ls   = Ls;
  dwfa.M5   = toDouble(M5);
  dwfa.mass = toDouble(mq);
  dwfa.Csw  = 0.0;
  dwfa.max_iter = 1000;
  dwfa.residual = 1.e-8;
  QDPIO::cout << "Initialising Gparity bfm operator\n";


  /********************************************************
   * Gaussian gauge field
   ********************************************************
   */


  multi1d<LatticeColorMatrix>  u_nogp(Nd);
  multi1d<LatticeColorMatrix>  ustar_nogp(Nd);
  multi1d<LatticeColorMatrix>  u(2*Nd);
  HotSt(u_nogp);
  for(int m=0;m<Nd;m++) ustar_nogp[m] = conj(u_nogp[m]);
  
  for(int m=0;m<Nd;m++) u[m] = u_nogp[m];
  for(int m=Nd; m< 2*Nd; ++m) u[m] = ustar_nogp[m-Nd];

  //need to put - sign on outwards-facing conjugate links in x-direction
  for(int y=0;y<nrow[1];y++){
    for(int z=0;z<nrow[2];z++){
      for(int t=0;t<nrow[3];t++){
	multi1d<int> yy(4); yy[0] = nrow[0]-1; yy[1] = y; yy[2] = z; yy[3] = t;
	ColorMatrix cj = Real(-1.0) * peekSite(u[Nd],yy);
	pokeSite(u[Nd],cj,yy);
      }
    }
  }
  

  /********************************************************
   * Gaussian source and result vectors
   ********************************************************
   */
  multi1d<LatticeFermion> psi_nogp_f0(Ls);
  multi1d<LatticeFermion> chi_nogp_f0(Ls);

  multi1d<LatticeFermion> psi_nogp_f1(Ls);
  multi1d<LatticeFermion> chi_nogp_f1(Ls);

  for(int s=0;s<Ls;s++) gaussian(psi_nogp_f0[s]);
  for(int s=0;s<Ls;s++) gaussian(chi_nogp_f0[s]);

  for(int s=0;s<Ls;s++) gaussian(psi_nogp_f1[s]);
  for(int s=0;s<Ls;s++) gaussian(chi_nogp_f1[s]);

  multi1d<LatticeFermion> psi(2*Ls);
  multi1d<LatticeFermion> chi(2*Ls);

  //currently expects 4d fields f0 and f1 on each s-slice
  int s=0;
  for(int i=0;i<2*Ls;i+=2){
    psi[i] = psi_nogp_f0[s];
    psi[i+1] = psi_nogp_f1[s];

    chi[i] = chi_nogp_f0[s];
    chi[i+1] = chi_nogp_f1[s];

    s++;
  }

  if(1){    
    bfm_t  dwf;
    QDPIO::cout << "Initialising Gparity bfm operator\n";
    dwf.init(dwfa);
    
    Fermion_t psi_h[2]; //even/odd
    psi_h[0] = dwf.allocFermion();
    psi_h[1] = dwf.allocFermion();
    Fermion_t chi_h[2];
    chi_h[0] = dwf.allocFermion();
    chi_h[1] = dwf.allocFermion();
    
    QDPIO::cout << "Importing gauge and fermion vectors\n";
    dwf.importGauge(u);

    dwf.importFermion(psi,psi_h[0],0);
    dwf.importFermion(psi,psi_h[1],1);
    
    dwf.importFermion(chi,chi_h[0],0);
    dwf.importFermion(chi,chi_h[1],1);
    
    dwf.inv_type=CG_PREC_MDAGM;
    dwf.qdp_chi_h[0]=chi_h[0];
    dwf.qdp_chi_h[1]=chi_h[1];
    dwf.qdp_psi_h[0]=psi_h[0];
    dwf.qdp_psi_h[1]=psi_h[1];

    QDPIO::cout << "Inverting\n";
    bfm_spawn_cg(dwf);


    // for(int i=0;i<NITER;i++) {
    //   QDPIO::cout << "Starting CGNE on iter " << i << "\n";
    //   dwf.CGNE(chi_h,psi_h);
    // }
    
    QDPIO::cout << "Finishing up\n";
    dwf.exportFermion(chi,chi_h[0],0);
    dwf.exportFermion(chi,chi_h[1],1);
    dwf.end();
  }
 
  //now do double lattice
  QDPIO::cout << "Preparing doubled lattice\n";
  int x[4];

  Double n2gp;
  for(int sg=0;sg<2*Ls;sg++) {
    n2gp += norm2(chi[sg]);
  }

  multi1d<multi1d<LatticeColorMatrix> > u_xblock(2);
  multi1d<multi1d<LatticeFermion> > psi_xblock(2);
  multi1d<multi1d<LatticeFermion> > chi_xblock(2);

  for(int i=0;i<2;i++){
    u_xblock[i].resize(Nd);
    psi_xblock[i].resize(Ls);
    chi_xblock[i].resize(Ls);
  }
  QDPIO::cout << "Preparing double gauge\n";
  create_blocks_gauge(u_xblock,u_nogp);
  QDPIO::cout << "Preparing double psi\n";
  create_blocks_ferm(psi_xblock,psi_nogp_f0,psi_nogp_f1,Ls);
  QDPIO::cout << "Preparing double chi\n";
  create_blocks_ferm(chi_xblock,chi_nogp_f0,chi_nogp_f1,Ls);

  QDPdouble* block_ptrs[2][Nd];
  QDPdouble* psi_ptrs[2][Ls];
  QDPdouble* chi_ptrs[2][Ls];
  for(int i=0;i<2;i++){
    for(int mu=0;mu<Nd;mu++){
      block_ptrs[i][mu] = (QDPdouble *)&u_xblock[i][mu].elem(0).elem();
    }
    for(int s=0;s<Ls;s++){
      psi_ptrs[i][s] = (QDPdouble*)&psi_xblock[i][s].elem(0).elem(0).elem(0).real();
      chi_ptrs[i][s] = (QDPdouble*)&chi_xblock[i][s].elem(0).elem(0).elem(0).real();
    }
  }

  //Chroma::finalize();
  //QDP_finalize();
  //QMP_finalize_msg_passing();
  //Do a doubled-lattice G-parity inversion
  nrow[0]*=2;

  //Chroma::initialize(&argc,&argv);
    
  Layout::setLattSize(nrow);
  Layout::create();

  lx = QDP::Layout::subgridLattSize()[0];
  ly = QDP::Layout::subgridLattSize()[1];
  lz = QDP::Layout::subgridLattSize()[2];
  lt = QDP::Layout::subgridLattSize()[3];

  ncoor = QDP::Layout::nodeCoord();
  nodes = QDP::Layout::logicalSize();

  QDPIO::cout << "Doubled lattice " << lx << " " << ly << " " << lz << " " << lt << endl;

  multi1d<LatticeColorMatrix>  u_double(Nd);
  for(int mu=0;mu<4;mu++){
    //copy gauge fields onto doubled lattice from stored blocks. All is local to node as we did comms before when setting up the blocks
    QDPdouble* doubled = (QDPdouble *)&u_double[mu].elem(0).elem();
    import_blocks_to_double_gauge(doubled,block_ptrs[0][mu],block_ptrs[1][mu]);

    multi1d<int> xx(4);
    //conjugate second half (which contains copies of first half)
    for(xx[0]=nrow[0]/2;xx[0]<nrow[0];xx[0]++){
      for(xx[1]=0;xx[1]<nrow[1];xx[1]++){
    	for(xx[2]=0;xx[2]<nrow[2];xx[2]++){
    	  for(xx[3]=0;xx[3]<nrow[3];xx[3]++){

	    ColorMatrix cj = conj(peekSite(u_double[mu],xx));
    	    pokeSite(u_double[mu],cj,xx);
    	  }
    	}
      }
    }

  }
  //put APBC on outwards facing links on G-parity boundary
  for(x[1]=0;x[1]<nrow[1];x[1]++){
    for(x[2]=0;x[2]<nrow[2];x[2]++){
      for(x[3]=0;x[3]<nrow[3];x[3]++){
	  multi1d<int> xm(4);
	  xm[0]=nrow[0]-1; xm[1] = x[1]; xm[2] = x[2]; xm[3] = x[3];
	  ColorMatrix cj = Real(-1.0) * peekSite(u_double[0],xm);
	  pokeSite(u_double[0],cj,xm);
      }
    }
  }



  bfm_t  dwf_nogp;
  dwfa.gparity=0;
  for(int i=0;i<3;i++) dwfa.gparity_dir[i]=0;

  dwfa.node_latt[0]  = lx;
  dwfa.node_latt[1]  = ly;
  dwfa.node_latt[2]  = lz;
  dwfa.node_latt[3]  = lt;

  QDPIO::cout << "Initialising Standard bfm operator\n";
  dwf_nogp.init(dwfa);

  Fermion_t psi_h_nogp[2];
  psi_h_nogp[0] = dwf_nogp.allocFermion();
  psi_h_nogp[1] = dwf_nogp.allocFermion();
  Fermion_t chi_h_nogp[2];
  chi_h_nogp[0] = dwf_nogp.allocFermion();
  chi_h_nogp[1] = dwf_nogp.allocFermion();

  multi1d<LatticeFermion> psi_nogp_both(Ls);
  multi1d<LatticeFermion> chi_nogp_both(Ls);

  for(int s=0;s<Ls;s++){
    import_blocks_to_double_ferm((QDPdouble*)&psi_nogp_both[s].elem(0).elem(0).elem(0).real(),psi_ptrs[0][s],psi_ptrs[1][s]);
    import_blocks_to_double_ferm((QDPdouble*)&chi_nogp_both[s].elem(0).elem(0).elem(0).real(),chi_ptrs[0][s],chi_ptrs[1][s]);
  }

  dwf_nogp.importGauge(u_double);

  /*Import this checkerboard of source field to bagel*/
  dwf_nogp.importFermion(psi_nogp_both,psi_h_nogp[0],0);
  dwf_nogp.importFermion(psi_nogp_both,psi_h_nogp[1],1);
    
  dwf_nogp.importFermion(chi_nogp_both,chi_h_nogp[0],0);
  dwf_nogp.importFermion(chi_nogp_both,chi_h_nogp[1],1);

  dwf_nogp.inv_type=CG_PREC_MDAGM;
  dwf_nogp.qdp_chi_h[0]=chi_h_nogp[0];
  dwf_nogp.qdp_chi_h[1]=chi_h_nogp[1];
  dwf_nogp.qdp_psi_h[0]=psi_h_nogp[0];
  dwf_nogp.qdp_psi_h[1]=psi_h_nogp[1];

  bfm_spawn_cg(dwf_nogp);


  // for(int i=0;i<NITER;i++) {
  //   dwf_nogp.CGNE(chi_h_nogp,psi_h_nogp);
  // }

  dwf_nogp.exportFermion(chi_nogp_both,chi_h_nogp[0],0);
  dwf_nogp.exportFermion(chi_nogp_both,chi_h_nogp[1],1);
  dwf_nogp.end();


  Double n2 = 0.0;
  for(int s=0;s<Ls;s++){
    Double sn2 = norm2(chi_nogp_both[s]);
    n2+=sn2;
  }
  QDPIO::cout << "|| (1flav) - (2flav) || = "<< n2-n2gp << " (" << n2 << "," << n2gp << ")\n";


  Chroma::finalize();

  QDPIO::cout << "Done\n";

}
