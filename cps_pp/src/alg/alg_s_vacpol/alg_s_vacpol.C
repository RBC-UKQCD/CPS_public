//------------------------------------------------------------------
//
// alg_s_vacpol.C
//
// T. Blum (June 2012)
//
// routines to calculate 2 pt correlation function
//------------------------------------------------------------------

#include <stdlib.h>	// exit()
#include <stdio.h>
#include <alg/alg_s_vacpol.h>
#include <alg/alg_s_spect.h>
#include <util/eigen_container.h>
#include <comms/glb.h>
#include <comms/scu.h>
#include <util/site.h>;
#include <util/fft.h>;
#include <util/qcdio.h>;
#include <util/time_cps.h>;

#define PI 3.141592654
#define NXMOM 5
#define NYMOM 5
#define NZMOM 5
#define NTMOM 17

CPS_START_NAMESPACE

QuarkPropSMng AlgStagQuark::sqpm;

//Sets a pointer to the lattice class to access smeared links
//void  set_pt (Fp4 *lat);
static Fp4 *lat_pt;

Rcomplex Zero(0.0,0.0);

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgVacPolStag::AlgVacPolStag(Fp4& latt, CommonArg *c_arg, StagQuarkArg *arg): 
  Alg(latt, c_arg)
{
  cname = "AlgVacPolStag";
  char *fname = "AlgVacPolStag(L&,CommonArg*,MuonArg*)";
  VRB.Func(cname,fname);

  lat_pt = (Fp4*)&latt;

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0) ERR.Pointer(cname,fname, "arg");
  sq_arg = *arg;
  // common_arg set in AlgBase
}
//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
// for twisted, or non-degenate, or ...
AlgVacPolStag::AlgVacPolStag(Fp4& latt, CommonArg *c_arg, StagQuarkArg *arg, StagQuarkArg *arg2):
  Alg(latt, c_arg)
{
  cname = "AlgVacPolStag";
  char *fname = "AlgVacPolStag(L&,CommonArg*,MuonArg*)";
  VRB.Func(cname,fname);

  lat_pt = (Fp4*)&latt;

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0) ERR.Pointer(cname,fname, "arg");
  sq_arg = *arg;
  if(arg2 == 0) ERR.Pointer(cname,fname, "arg2");
  sq_arg_twisted = *arg2;
  // common_arg set in AlgBase
}

//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgVacPolStag::~AlgVacPolStag() {
  char *fname = "~AlgVacPolStag()";
  VRB.Func(cname,fname);
  //sqpm.destroyQuarkPropS(sq_arg->qid);
}

//------------------------------------------------------------------
// Conserved-conserved vector current 2 pt function
//------------------------------------------------------------------
void AlgVacPolStag::VacPolStagConsConsAMA(int Start, int Inc, int tStart, int tInc)
{

  int n[4], p[4], x[4], nu, mu;

  // local sizes
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();

  // shifts to get global coords/mom
  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();
  
  int size = 4*4; // 4 dirs * 4 dirs is # of complex to FFT/site

  // mom
  Float q[4];
  // global sizes
  Float L[4];
  L[0] = n[0]*GJP.Xnodes();
  L[1] = n[1]*GJP.Ynodes();
  L[2] = n[2]*GJP.Znodes();
  L[3] = n[3]*GJP.Tnodes();
  // mom factors
  Float PIL0 = PI/L[0];
  Float PIL1 = PI/L[1];
  Float PIL2 = PI/L[2];
  Float PIL3 = PI/L[3];
  Float TWOPIL0 = 2*PIL0;
  Float TWOPIL1 = 2*PIL1;
  Float TWOPIL2 = 2*PIL2;
  Float TWOPIL3 = 2*PIL3;

  // twist mom. need for the 1/2 mom factors from point-splitting
  // note that this is 1/2*theta
  Float thetaover2[4]; 
  for(int mu=0; mu<4; ++mu){
    thetaover2[mu] = PI*GJP.TwistBc(mu)/GJP.Sites(mu);
  }
  
  Rcomplex *VacPol; // V(q,y)
  VacPol = (Rcomplex*)smalloc(4*4*GJP.VolNodeSites()*sizeof(Rcomplex));
  Rcomplex *VacPolSrc; // V(q,q) (usually won't do all momenta on source end)
  VacPolSrc = (Rcomplex*)smalloc(4*4*GJP.VolNodeSites()*sizeof(Rcomplex));

  
  // zero HVP
  for(p[3] = 0; p[3] < n[3]; p[3]++){
    for(p[2] = 0; p[2] < n[2]; p[2]++){
      for(p[1] = 0; p[1] < n[1]; p[1]++){
	for(p[0] = 0; p[0] < n[0]; p[0]++){
	  
	  int site = p[0] + n[0]*(p[1] +n[1]*(p[2] + n[2]*(p[3])));
	  
	  for(nu=0;nu<4;nu++){
	    for(mu=0;mu<4;mu++){
	      VacPolSrc[nu+4*(mu+4*site)]=Zero;
	    }
	  }
	}
      }
    }
  }
  
  int cnt=0;
  //loop over point sources
  for(x[3] = tStart; x[3] < L[3]; x[3]+=tInc){
    for(x[2] = Start; x[2] < L[2]; x[2]+=Inc){
      for(x[1] = Start; x[1] < L[1]; x[1]+=Inc){
	for(x[0] = Start; x[0] < L[0]; x[0]+=Inc){
	  
	  sq_arg.src.origin[0] = sq_arg_twisted.src.origin[0] = x[0];
	  sq_arg.src.origin[1] = sq_arg_twisted.src.origin[1] = x[1];
	  sq_arg.src.origin[2] = sq_arg_twisted.src.origin[2] = x[2];
	  sq_arg.src.origin[3] = sq_arg_twisted.src.origin[3] = x[3];
	  
	  //VacPolStagConsCons(VacPol);
	  VacPolStagConsConsTwistedBC(VacPol);

	  //add 1/2 point split phases and print
	  // source momenta
	  for(p[3] = 0; p[3] < n[3]; p[3]++){
	    int ptg = p[3]+shiftT;
	    for(p[2] = 0; p[2] < n[2]; p[2]++){
	      int pzg = p[2]+shiftZ;
	      for(p[1] = 0; p[1] < n[1]; p[1]++){
		int pyg = p[1]+shiftY;
		for(p[0] = 0; p[0] < n[0]; p[0]++){
		  int pxg = p[0]+shiftX;
		  
		  if(pxg > NXMOM && pxg < L[0] - NXMOM ) continue;
		  if(pyg > NYMOM && pyg < L[1] - NYMOM ) continue;
		  if(pzg > NZMOM && pzg < L[2] - NZMOM ) continue;
		  if(ptg > NTMOM && ptg < L[3] - NTMOM ) continue;
		  
		  int site = p[0] + n[0]*(p[1] +n[1]*(p[2] + n[2]*(p[3])));
		
		  // phase p.y
		  q[0] = TWOPIL0*pxg*x[0];
		  q[1] = TWOPIL1*pyg*x[1];
		  q[2] = TWOPIL2*pzg*x[2];
		  q[3] = TWOPIL3*ptg*x[3];
		  Float theta_src = (q[0]+q[1]+q[2]+q[3]);

		  for(nu=0;nu<4;nu++){
		    // full point-split source phase
		    q[0] = nu==0 ? PIL0*pxg+thetaover2[0] : 0;
		    q[1] = nu==1 ? PIL1*pyg+thetaover2[1] : 0;
		    q[2] = nu==2 ? PIL2*pzg+thetaover2[2] : 0; 
		    q[3] = nu==3 ? PIL3*ptg+thetaover2[3] : 0;
		    Float theta_src2 = (q[0]+q[1]+q[2]+q[3]);
		    for(mu=0;mu<4;mu++){
		      // extra 1/2 point-split sink phase (FFT already done on sink)
		      q[0] = mu==0 ? PIL0*pxg+thetaover2[0] : 0;
		      q[1] = mu==1 ? PIL1*pyg+thetaover2[1] : 0;
		      q[2] = mu==2 ? PIL2*pzg+thetaover2[2] : 0; 
		      q[3] = mu==3 ? PIL3*ptg+thetaover2[3] : 0;
		      Float theta = theta_src + theta_src2 -(q[0]+q[1]+q[2]+q[3]);
		      Rcomplex phase(cos(theta),sin(theta));
		      
		      VacPolSrc[nu+4*(mu+4*site)] += phase*VacPol[nu+4*(mu+4*site)];
		    }
		  }
		}
	      }
	    }
	  }
	  cnt++;
	}
      }
    }
  }//src points

  //char filename[128];
  //sprintf(filename,"hvp-cons-cons.dat");
  FileIoType ft= ADD_ID;
  FILE *fp=Fopen(ft,(char*)common_arg->results, "w");

  //print
  for(p[3] = 0; p[3] < n[3]; p[3]++){
    int tg = p[3]+shiftT;
    for(p[2] = 0; p[2] < n[2]; p[2]++){
      int zg = p[2]+shiftZ;
      for(p[1] = 0; p[1] < n[1]; p[1]++){
	int yg = p[1]+shiftY;
	for(p[0] = 0; p[0] < n[0]; p[0]++){
	  int xg = p[0]+shiftX;
	  
          if(xg > NXMOM && xg < L[0] - NXMOM ) continue;
          if(yg > NYMOM && yg < L[1] - NYMOM ) continue;
          if(zg > NZMOM && zg < L[2] - NZMOM ) continue;
          if(tg > NTMOM && tg < L[3] - NTMOM ) continue;
	  
	  int site = p[0] + n[0]*(p[1] +n[1]*(p[2] + n[2]*(p[3])));
	  
	  for(mu=0;mu<4;mu++){
	    for(nu=0;nu<4;nu++){
	      Fprintf(ft,fp,"VACPOL0 %d %d  mom %d %d %d %d  %e %e\n",
		      mu, nu,
		      xg, yg, zg, tg,
		      VacPolSrc[nu+4*(mu+4*site)].real()/cnt,VacPolSrc[nu+4*(mu+4*site)].imag()/cnt
		      );
	    }
	  }
	}
      }
    }
  } //momentum
  Fclose(ft,fp);

  sfree(VacPolSrc);
  sfree(VacPol);
 
}





void AlgVacPolStag::VacPolStagConsConsTwistedBC(Rcomplex* VacPol)
{

  char *fname = "VacPolStagConsCons()";
  VRB.Func(cname,fname);
  
  int n[4], x[4], x_n[4];
  
  int i,j,k;
  int source_i;
  int icol,cgn;
  int mu, nu;
  Rcomplex contact[4];
  Rcomplex contact2[4];
  Rcomplex contactJ[4];
  Rcomplex vacpol0;
  Float x1,x2;
  Float phase;
  Float finalrsq;
  Rcomplex cc;
  Matrix temp, temp2;
  Matrix sourcelink[4];
  Matrix *link;
  Matrix propnu_neigh;
  Matrix prop_neigh;
  
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();

  // global sizes
  int L[4];
  L[0] = n[0]*GJP.Xnodes();
  L[1] = n[1]*GJP.Ynodes();
  L[2] = n[2]*GJP.Znodes();
  L[3] = n[3]*GJP.Tnodes();

  int size = 4*4; // 4 dirs * 4 dirs
  
  Rcomplex I(0.0,1.0);
 
  // The quark props
  //------------------------------------

  // set the point-split locations. This is global!
  StagQuarkArg sq_arg0 = sq_arg;
  sq_arg0.src.origin[0]= (sq_arg0.src.origin[0]+1)%L[0];
  sq_arg0.qid=0;
  StagQuarkArg sq_arg1 = sq_arg;
  sq_arg1.src.origin[1]= (sq_arg1.src.origin[1]+1)%L[1];
  sq_arg1.qid=1;
  StagQuarkArg sq_arg2 = sq_arg;
  sq_arg2.src.origin[2]= (sq_arg2.src.origin[2]+1)%L[2];
  sq_arg2.qid=2;
  StagQuarkArg sq_arg3 = sq_arg;
  sq_arg3.src.origin[3]= (sq_arg3.src.origin[3]+1)%L[3];
  sq_arg3.qid=3;
  sq_arg.qid=4;

  // for conserved point split current at source: 1 centered, 4 shifted
  QuarkPropS quark[5] = {
    QuarkPropS (AlgLattice(), sq_arg0),
    QuarkPropS (AlgLattice(), sq_arg1),
    QuarkPropS (AlgLattice(), sq_arg2),
    QuarkPropS (AlgLattice(), sq_arg3),
    QuarkPropS (AlgLattice(), sq_arg) 
  };

  sq_arg0.qid+=5;
  sq_arg0. cg. fname_eigen =   sq_arg_twisted. cg. fname_eigen ;
  sq_arg1.qid+=5;
  sq_arg1. cg. fname_eigen =   sq_arg_twisted. cg. fname_eigen ;
  sq_arg2.qid+=5;
  sq_arg2. cg. fname_eigen =   sq_arg_twisted. cg. fname_eigen ;
  sq_arg3.qid+=5;
  sq_arg3. cg. fname_eigen =   sq_arg_twisted. cg. fname_eigen ;
  sq_arg_twisted.qid=9;
  // for conserved point split current at source: 1 centered, 4 shifted
  // twisted b.c.
  QuarkPropS quarkt[5] = { 
    QuarkPropS (AlgLattice(), sq_arg0),
    QuarkPropS (AlgLattice(), sq_arg1),
    QuarkPropS (AlgLattice(), sq_arg2),
    QuarkPropS (AlgLattice(), sq_arg3),
    QuarkPropS (AlgLattice(), sq_arg_twisted) 
  };
  
  Float etime = time_elapse();

  Float theta[4];
  Rcomplex etheta[4];
  for(int mu=0; mu<4; ++mu){
    theta[mu] = 2.0*PI*GJP.TwistBc(mu)/GJP.Sites(mu);
    Rcomplex cc(cos(theta[mu]),sin(theta[mu]));
    etheta[mu] = cc;
  }

  // get props

  // twisted first so smeared links are untwisted at the end (for cons-curr)
  if(!UniqueID()) printf("Twisting the lattice in %s\n",fname);
  AlgLattice().twist_bc(1);
  for(int nu=0; nu<=4; nu++){
    quarkt[nu].setupQuarkPropS();
    quarkt[nu].getQuarkPropS((char*)common_arg->results);
  }
  //ecache->dealloc();
  // un-twist
  if(!UniqueID()) printf("Un-Twisting the lattice in %s\n",fname);
  AlgLattice().twist_bc(-1);
  etime = time_elapse();
  if(!UniqueID())printf("Time for twisted quark props %g secs\n",etime);

  //snprintf(fname2,1024, "%s.bc%d%d%d%d", sq_arg.cg.fname_eigen, GJP.Bc(0),GJP.Bc(1),GJP.Bc(2),GJP.Bc(3));
  //ecache->alloc( fname2, neig, fsize );

  sq_arg0. cg. fname_eigen =   sq_arg. cg. fname_eigen ;
  sq_arg1. cg. fname_eigen =   sq_arg. cg. fname_eigen ;
  sq_arg2. cg. fname_eigen =   sq_arg. cg. fname_eigen ;
  sq_arg3. cg. fname_eigen =   sq_arg. cg. fname_eigen ;
  for(int nu=0; nu<=4; nu++){
    quark[nu].setupQuarkPropS();
    quark[nu].getQuarkPropS((char*)common_arg->results);
  }
  //ecache->dealloc();
  etime = time_elapse();
  if(!UniqueID())printf("Time for quark props %g secs\n",etime);


  // later, the point-split correlation function 
  // will need the link at source point, and parity
  x[0] = sq_arg.src.origin[0] ;
  x[1] = sq_arg.src.origin[1] ;
  x[2] = sq_arg.src.origin[2] ;
  x[3] = sq_arg.src.origin[3] ;
  // set point source
  int procCoorX = x[0] / GJP.XnodeSites();
  int procCoorY = x[1] / GJP.YnodeSites();
  int procCoorZ = x[2] / GJP.ZnodeSites();
  int procCoorT = x[3] / GJP.TnodeSites();
  int localX = x[0] % GJP.XnodeSites();
  int localY = x[1] % GJP.YnodeSites();
  int localZ = x[2] % GJP.ZnodeSites();
  int localT = x[3] % GJP.TnodeSites();  
  int coor_x = GJP.XnodeCoor();
  int coor_y = GJP.YnodeCoor();
  int coor_z = GJP.ZnodeCoor();
  int coor_t = GJP.TnodeCoor();
  /* phase from conjugating quark -> antiquark */
  Site site(localX,localY,localZ,localT);
  int src_parity = site.physX()+site.physY()+site.physZ()+site.physT();
  src_parity = src_parity % 2;

  int vol =  GJP.VolNodeSites();

  for(int nu=0;nu<4;nu++){
    sourcelink[nu].ZeroMatrix();
    contact[nu] = Zero;
    if (coor_x == procCoorX && coor_y == procCoorY && coor_z == procCoorZ && coor_t == procCoorT){
      Site site(localX,localY,localZ,localT);
      int id = site.Index();
      // first V_x over lattice, then V_y, etc.
      link = lat_pt->Fields(0) + nu*vol + id;
      sourcelink[nu] = *link;
    }
  }

  /* "broadcast" it */
  for(int nu=0;nu<4;nu++){
    slice_sum((Float*)&sourcelink[nu], 18, 99);
  }

  /* source dir */
  for(nu=0;nu<4;nu++){

    Matrix hsrclink;
    hsrclink.Dagger(sourcelink[nu]);

    /* sink dir */
    for(mu=0;mu<4;mu++){
              
      /* calculate the two-point function */
      Site site;
      while ( site.LoopsOverNode() ){

	Matrix temp;
	Matrix temp2;
	Matrix temp3;

	int i=site.Index();
	int neighbor_site=site.plusIndex(mu);
	
	Matrix *link = lat_pt->Fields(0) + i + vol*mu;

	Matrix hlink;
	hlink.Dagger(*link);

	// get prop in +mu dir.
	// need the transpose (conjugate later)
	for(int color=0;color<3;color++){
	  if( /* offsite */ site.Coor(mu)+1 == n[mu] ){
	    getPlusData( (IFloat *)&propnu_neigh+6*color,
			 quark[nu].AccessQuarkPropS(color,neighbor_site), 6, mu);
	  } else {
	    moveFloat((Float*)&propnu_neigh+6*color, 
		      quark[nu].AccessQuarkPropS(color,neighbor_site),6);
	  }
	}
	propnu_neigh.Transpose();
	temp.DotMEqual(*link,propnu_neigh); 
	temp2.Dagger(temp);
	Matrix prop;
        // this one is twisted
	for(int color=0;color<3;color++)
	  moveFloat((Float*)&prop+6*color,quarkt[4].AccessQuarkPropS(color,i),6);
	prop.Transpose();
	temp.DotMEqual(temp2,prop);
	temp2.DotMEqual(temp,sourcelink[nu]);
	Rcomplex cc= temp2.Tr();
        /* phase from conjugating quark -> antiquark */
	int parity = site.physX()+site.physY()+site.physZ()+site.physT();
	parity = parity % 2;
        phase = (src_parity-parity==0 ? 1 : -1);
        cc *= phase;
        VacPol[nu+4*(mu+4*i)] = cc;

        // get twisted prop in +mu dir.
	for(int color=0;color<3;color++){
	  if( /* offsite */ site.Coor(mu)+1 == n[mu] ){
	    getPlusData( (IFloat *)&prop_neigh+6*color,
			 quarkt[4].AccessQuarkPropS(color,neighbor_site), 6, mu);
	  } else {
	    moveFloat((Float*)&prop_neigh+6*color, 
		      quarkt[4].AccessQuarkPropS(color,neighbor_site),6);
	  }
	}
	prop_neigh.Transpose();
	for(int color=0;color<3;color++)
	  moveFloat((Float*)&temp+6*color,quark[nu].AccessQuarkPropS(color,i),6);
	temp.Transpose();
        temp2.Dagger(temp);
        temp.DotMEqual(temp2,*link);
        temp2.DotMEqual(temp,prop_neigh);
        temp.DotMEqual(temp2,sourcelink[nu]);
        cc= temp.Tr();
        phase = (src_parity-parity==0 ? -1 : 1);
        cc *= phase;
        VacPol[nu+4*(mu+4*i)] += etheta[mu]*cc;


        
         // get prop in +mu dir.
	for(int color=0;color<3;color++){
	  if( /* offsite */ site.Coor(mu)+1 == n[mu] ){
	    getPlusData( (IFloat *)&prop_neigh+6*color,
			 quark[4].AccessQuarkPropS(color,neighbor_site), 6, mu);
	  } else {
	    moveFloat((Float*)&prop_neigh+6*color, 
		      quark[4].AccessQuarkPropS(color,neighbor_site),6);
	  }
	}
	prop_neigh.Transpose();
	temp.DotMEqual(*link,prop_neigh);
        temp2.Dagger(temp);
	for(int color=0;color<3;color++)
	  moveFloat((Float*)&temp+6*color,quarkt[nu].AccessQuarkPropS(color,i),6);
	temp.Transpose();
        temp3.DotMEqual(temp2,temp);
        temp2.DotMEqual(temp3,hsrclink);
        cc= temp2.Tr();
        phase = (src_parity-parity==0 ? -1 : 1);
        cc *= phase;
        VacPol[nu+4*(mu+4*i)] += conj(etheta[nu])*cc;

	// get twisted prop in +mu dir.
	// need the transpose (conjugate later)
	for(int color=0;color<3;color++){
	  if( /* offsite */ site.Coor(mu)+1 == n[mu] ){
	    getPlusData( (IFloat *)&propnu_neigh+6*color,
			 quarkt[nu].AccessQuarkPropS(color,neighbor_site), 6, mu);
	  } else {
	    moveFloat((Float*)&propnu_neigh+6*color, 
		      quarkt[nu].AccessQuarkPropS(color,neighbor_site),6);
	  }
	}
	propnu_neigh.Transpose();
	for(int color=0;color<3;color++)
	  moveFloat((Float*)&temp+6*color,quark[4].AccessQuarkPropS(color,i),6);
	temp.Transpose();
        temp2.Dagger(temp);
        temp.DotMEqual(temp2,*link);
        temp2.DotMEqual(temp,propnu_neigh);
        temp.DotMEqual(temp2,hsrclink);
        cc= temp.Tr();
        phase = (src_parity-parity==0 ? 1 : -1);
	cc *= phase;
        VacPol[nu+4*(mu+4*i)] += etheta[mu]*conj(etheta[nu])*cc;

#if 0
	// twisted fields 
	// this phase is necessary for coordinate space vacpol. But in momentum space it is exactly
        // what is needed for FT w/r to p+theta/L momenta, so do nothing for momentum space Pi
	phase = theta[0]*site.physX()+theta[1]*site.physY()+theta[2]*site.physZ()+theta[3]*site.physT();
	phase -= theta[0]*sq_arg.src.origin[0]+theta[1]*sq_arg.src.origin[1]+theta[2]*sq_arg.src.origin[2]+theta[3]*sq_arg.src.origin[3];
	Rcomplex ccc(cos(phase),sin(phase));
        VacPol[nu+4*(mu+4*i)] *= ccc;
#endif
      } /* sites */
    } /* mu */

  
    /* contact terms */
    // delta(x-y+nu)
    x[0] = sq_arg.src.origin[0] ;
    x[1] = sq_arg.src.origin[1] ;
    x[2] = sq_arg.src.origin[2] ;
    x[3] = sq_arg.src.origin[3] ;
    int x0 = (nu==0) ? (1+x[0])%(GJP.XnodeSites()*GJP.Xnodes()) : x[0];
    int y0 = (nu==1) ? (1+x[1])%(GJP.YnodeSites()*GJP.Ynodes()) : x[1];
    int z0 = (nu==2) ? (1+x[2])%(GJP.ZnodeSites()*GJP.Znodes()) : x[2];
    int t0 = (nu==3) ? (1+x[3])%(GJP.TnodeSites()*GJP.Tnodes()) : x[3];		 
    procCoorX = x0 / GJP.XnodeSites();
    procCoorY = y0 / GJP.YnodeSites();
    procCoorZ = z0 / GJP.ZnodeSites();
    procCoorT = t0 / GJP.TnodeSites();
    localX = x0 % GJP.XnodeSites();
    localY = y0 % GJP.YnodeSites();
    localZ = z0 % GJP.ZnodeSites();
    localT = t0 % GJP.TnodeSites();
    coor_x = GJP.XnodeCoor();
    coor_y = GJP.YnodeCoor();
    coor_z = GJP.ZnodeCoor();
    coor_t = GJP.TnodeCoor();

    //Float phi = theta[0] * x0;
    //phi += theta[1] * y0;
    //phi += theta[2] * z0;
    //phi += theta[3] * t0;
    Float phi = theta[nu];
    Rcomplex cc2(cos(phi),sin(phi));
    //Rcomplex cc2(1.0,0.0);

    if (coor_x == procCoorX && coor_y == procCoorY && coor_z == procCoorZ && coor_t == procCoorT){

      Site site(localX,localY,localZ,localT);
      i=site.Index();
      int i_neigh=site.minusIndex(nu);

      for(int color=0;color<3;color++){
	if( /* offsite */ site.Coor(nu)-1 < 0 ){
	  getMinusData( (IFloat *)&temp+6*color,
		       quark[nu].AccessQuarkPropS(color,i_neigh), 6, nu);
	} else {
	  moveFloat((Float*)&temp+6*color, 
		    quark[nu].AccessQuarkPropS(color,i_neigh),6);
	}
      }
      temp.Transpose();
      Matrix hlink;
      hlink.Dagger(sourcelink[nu]);
      temp2.DotMEqual(temp,hlink);
      cc = conj(etheta[nu])*cc2*temp2.Tr();
#if 0
      printf("CONTACT PER  %d  x %d %d %d %d   %e %e\n",
	     nu,
	     site.physX(),site.physY(),site.physZ(),site.physT(),
	     cc.real(),cc.imag()
	     );
#endif
      contact[nu] = 0.5*cc;
      contact2[nu] = -cc; // sub
      contactJ[nu] = -cc; // sub

      for(int color=0;color<3;color++)
	moveFloat((Float*)&temp+6*color,quarkt[4].AccessQuarkPropS(color,i),6);
      temp.Transpose();
      temp2.DotMEqual(temp,sourcelink[nu]);
      cc = -cc2*temp2.Tr();
#if 0
      printf("CONTACT TWISTED  %d  x %d %d %d %d   %e %e\n",
	     nu,
	     site.physX(),site.physY(),site.physZ(),site.physT(),
	     cc.real(),cc.imag()
	     );
#endif
      contact[nu] += 0.5*cc;
      contact2[nu] -= cc; // add (- -)
      // avg per+twisted bc
    }
    /* "broadcast" it */
    slice_sum((Float*)&contact[nu], 2, 99);
    slice_sum((Float*)&contact2[nu], 2, 99);
    slice_sum((Float*)&contactJ[nu], 2, 99);

    // delta(x-y)
    x0 = sq_arg.src.origin[0] ;
    y0 = sq_arg.src.origin[1] ;
    z0 = sq_arg.src.origin[2] ;
    t0 = sq_arg.src.origin[3] ;
    procCoorX = x0 / GJP.XnodeSites();
    procCoorY = y0 / GJP.YnodeSites();
    procCoorZ = z0 / GJP.ZnodeSites();
    procCoorT = t0 / GJP.TnodeSites();
    localX = x0 % GJP.XnodeSites();
    localY = y0 % GJP.YnodeSites();
    localZ = z0 % GJP.ZnodeSites();
    localT = t0 % GJP.TnodeSites();
    coor_x = GJP.XnodeCoor();
    coor_y = GJP.YnodeCoor();
    coor_z = GJP.ZnodeCoor();
    coor_t = GJP.TnodeCoor();

    phi = theta[0] * x0;
    phi += theta[1] * y0;
    phi += theta[2] * z0;
    phi += theta[3] * t0;
    //Rcomplex cc3(cos(phi),-sin(phi));
    Rcomplex cc3(1.0,0.0);

    if (coor_x == procCoorX && coor_y == procCoorY && coor_z == procCoorZ && coor_t == procCoorT){

      Site site(localX,localY,localZ,localT);
      i=site.Index();
      int i_neigh=site.plusIndex(nu);

      for(int color=0;color<3;color++){
	if( /* offsite */ site.Coor(nu)+1 == n[nu] ){
	  getPlusData( (IFloat *)&temp+6*color,
		       quark[4].AccessQuarkPropS(color,i_neigh), 6, nu);
	} else {
	  moveFloat((Float*)&temp+6*color, 
		    quark[4].AccessQuarkPropS(color,i_neigh),6);
	}
      }
      temp.Transpose();
      temp2.DotMEqual(temp,sourcelink[nu]);
      cc = -cc3*temp2.Tr();

#if 0
      printf("CONTACT PER  %d  x %d %d %d %d   %e %e\n",
	     nu,
	     site.physX(),site.physY(),site.physZ(),site.physT(),
	     cc.real(),cc.imag()
	     );
#endif
      
      /* factor of 1/2 for current, 2 for each prop for
         field renorm. So, do nothing.*/
      contact[nu] += 0.5*cc;
      contact2[nu] += cc; // sub (- above)
      contactJ[nu] += cc; // sub (- above)

      for(int color=0;color<3;color++)
	moveFloat((Float*)&temp+6*color,quarkt[nu].AccessQuarkPropS(color,i),6);
      temp.Transpose();
      Matrix hlink;
      hlink.Dagger(sourcelink[nu]);
      temp2.DotMEqual(temp,hlink);
      cc = cc3*conj(etheta[nu])*temp2.Tr();

#if 0
      printf("CONTACT TWISTED %d  x %d %d %d %d   %e %e\n",
	     nu,
	     site.physX(),site.physY(),site.physZ(),site.physT(),
	     cc.real(),cc.imag()
	     );
#endif
      
      /* factor of 1/2 for current, 2 for each prop for
         field renorm. So, do nothing.*/
      contact[nu] += 0.5*cc;

      VacPol[nu+4*(nu+4*i)] += contact[nu];

      contact2[nu] += cc; // (add)
      contact2[nu] *= 0.5;
      contactJ[nu] *= 0.5;
      printf("CONTACT2 %d  x %d %d %d %d   %e %e\n",
	     nu,
	     site.physX(),site.physY(),site.physZ(),site.physT(),
	     contact2[nu].real(),contact2[nu].imag()
	     );
      printf("CONTACTJ %d  x %d %d %d %d   %e %e\n",
	     nu,
	     site.physX(),site.physY(),site.physZ(),site.physT(),
	     contactJ[nu].real(),contactJ[nu].imag()
	     );
    }
  } /* nu */

  etime = time_elapse();
  if(!UniqueID())printf("Time for HVP contraction %g secs\n",etime);

#if 0
  //print out PI before FFT as check
  {
    Site mysite;
    while ( mysite.LoopsOverNode() ){
      
      int i=mysite.Index();
      
      for(nu=0;nu<4;nu++){
	for(mu=0;mu<4;mu++){
	  if(!UniqueID())printf("PI %d %d  x %d %d %d %d   %e %e\n",
				mu,nu,
				mysite.physX(),mysite.physY(),mysite.physZ(),mysite.physT(),
				VacPol[nu+4*(mu+4*i)].real(),VacPol[nu+4*(mu+4*i)].imag()
				);
	}
      }
    }
  }
#endif

#if 0
  //print out props
  {
    Site mysite;
    while ( mysite.LoopsOverNode() ){
      
      int i=mysite.Index();
      
      for(mu=0;mu<5;mu++){
	Vector vec;
	moveFloat((Float*)&vec, quarkt[mu].AccessQuarkPropS(0,i),6);
	Rcomplex cc(*((Float*)&vec), *((Float*)&vec+1));
	phase = theta[0]*mysite.physX()+theta[1]*mysite.physY()+theta[2]*mysite.physZ()+theta[3]*mysite.physT();
	Rcomplex ccc(cos(phase),sin(phase));
	cc *= ccc;
	if(!UniqueID())printf("prop y %d  x %d %d %d %d   %e %e\n",
			      mu,
			      mysite.physX(),mysite.physY(),mysite.physZ(),mysite.physT(),
			      cc.real(), cc.imag());
      }
    }
  }

#endif


#if 0

  //print out divergence
  {
    Site mysite;
    while ( mysite.LoopsOverNode() ){
      
      int i=mysite.Index();
      
      for(nu=0;nu<4;nu++){
	Rcomplex Div(0.,0.);
	for(mu=0;mu<4;mu++){

	  int i_neigh=mysite.minusIndex(mu);

	  Rcomplex cc;
	  if( /* offsite */ mysite.Coor(mu)-1 < 0 ){
	    getMinusData( (IFloat *)&cc, (Float*)&VacPol[nu+4*(mu+4*i_neigh)], 2, mu);
	  } else {
	    cc = VacPol[nu+4*(mu+4*i_neigh)];
	  }
	  if( /* 0 boundary */ mysite.physCoor(mu)== 0 ){
	    phase = theta[mu]*GJP.Sites(mu);
	    Rcomplex ccc(cos(phase),-sin(phase));
	    cc*=ccc;
	  }
	  
	  //Div += VacPol[nu+4*(mu+4*i)]-conj(etheta[mu])*cc;
	  // Twisted difference
	  Div += VacPol[nu+4*(mu+4*i)]-cc;
	  
	}
	
	printf("DIV-VP nu %d  x %d %d %d %d   %e %e\n",
	       nu,
	       mysite.physX(),mysite.physY(),mysite.physZ(),mysite.physT(),
	       Div.real(),Div.imag()
	       );
      }
    }
  }
#endif

  // Fourier transform. "size" is number of complex
  // "1" last arg. is for "forward" fft, i.e. exp -ipx
  // to match our milc code
  FFT4((Float*)VacPol, 1, 2*size, 1, 1);
  etime = time_elapse();
  if(!UniqueID())printf("Time for HVP FFT %g secs\n",etime);

  quarkt[0].destroyQuarkPropS(sq_arg0.qid);
  quarkt[1].destroyQuarkPropS(sq_arg1.qid);
  quarkt[2].destroyQuarkPropS(sq_arg2.qid);
  quarkt[3].destroyQuarkPropS(sq_arg3.qid);
  quarkt[4].destroyQuarkPropS(sq_arg_twisted.qid);

  sq_arg0.qid-=5;
  sq_arg1.qid-=5;
  sq_arg2.qid-=5;
  sq_arg3.qid-=5;

  quark[0].destroyQuarkPropS(sq_arg0.qid);
  quark[1].destroyQuarkPropS(sq_arg1.qid);
  quark[2].destroyQuarkPropS(sq_arg2.qid);
  quark[3].destroyQuarkPropS(sq_arg3.qid);
  quark[4].destroyQuarkPropS(sq_arg.qid);

}



//------------------------------------------------------------------
// Conserved-conserved vector current 2 pt function
//------------------------------------------------------------------
void AlgVacPolStag::VacPolStagConsCons(Rcomplex* VacPol)
{

  char *fname = "VacPolStagConsCons()";
  VRB.Func(cname,fname);
  
  int n[4], x[4], x_n[4];
  
  int i,j,k;
  int source_i;
  int icol,cgn;
  int mu, nu;
  Rcomplex contact[4];
  Rcomplex vacpol0;
  Float x1,x2;
  Float phase;
  Float finalrsq;
  Rcomplex cc;
  Matrix temp, temp2;
  Matrix sourcelink[4];
  Matrix *link;
  Matrix propnu_neigh;
  Matrix prop_neigh;
  
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();

  // global sizes
  int L[4];
  L[0] = n[0]*GJP.Xnodes();
  L[1] = n[1]*GJP.Ynodes();
  L[2] = n[2]*GJP.Znodes();
  L[3] = n[3]*GJP.Tnodes();

  int size = 4*4; // 4 dirs * 4 dirs
  
  Rcomplex I(0.0,1.0);
 
  // The quark props
  //------------------------------------
  // set the point-split locations. This is global!
  StagQuarkArg sq_arg0 = sq_arg;
  sq_arg0.src.origin[0]= (sq_arg0.src.origin[0]+1)%L[0];
  sq_arg0.qid=0;
  StagQuarkArg sq_arg1 = sq_arg;
  sq_arg1.src.origin[1]= (sq_arg1.src.origin[1]+1)%L[1];
  sq_arg1.qid=1;
  StagQuarkArg sq_arg2 = sq_arg;
  sq_arg2.src.origin[2]= (sq_arg2.src.origin[2]+1)%L[2];
  sq_arg2.qid=2;
  StagQuarkArg sq_arg3 = sq_arg;
  sq_arg3.src.origin[3]= (sq_arg3.src.origin[3]+1)%L[3];
  sq_arg3.qid=3;
  sq_arg.qid=4;
  // for conserved point split current at source: 1 centered, 4 shifted
  QuarkPropS quark[5] = { 
    QuarkPropS (AlgLattice(), sq_arg0),
    QuarkPropS (AlgLattice(), sq_arg1),
    QuarkPropS (AlgLattice(), sq_arg2),
    QuarkPropS (AlgLattice(), sq_arg3),
    QuarkPropS (AlgLattice(), sq_arg) 
  };


  // get props
  Float etime = time_elapse();
  for(int nu=0; nu<=4; nu++){
    quark[nu].setupQuarkPropS();
    quark[nu].getQuarkPropS((char*)common_arg->results);
    //for(int c=0;c<3;c++)
      //sol = (Vector*)quark[nu].AccessQuarkPropS(c,0);
    //AlgLattice().Convert(CANONICAL,sol);
  }
  etime = time_elapse();
  if(!UniqueID())printf("Time for quark props %g secs\n",etime);

  // later, the point-split correlation function 
  // will need the link at source point, and parity
  x[0] = sq_arg.src.origin[0] ;
  x[1] = sq_arg.src.origin[1] ;
  x[2] = sq_arg.src.origin[2] ;
  x[3] = sq_arg.src.origin[3] ;
  // set point source
  int procCoorX = x[0] / GJP.XnodeSites();
  int procCoorY = x[1] / GJP.YnodeSites();
  int procCoorZ = x[2] / GJP.ZnodeSites();
  int procCoorT = x[3] / GJP.TnodeSites();
  int localX = x[0] % GJP.XnodeSites();
  int localY = x[1] % GJP.YnodeSites();
  int localZ = x[2] % GJP.ZnodeSites();
  int localT = x[3] % GJP.TnodeSites();  
  int coor_x = GJP.XnodeCoor();
  int coor_y = GJP.YnodeCoor();
  int coor_z = GJP.ZnodeCoor();
  int coor_t = GJP.TnodeCoor();
  /* phase from conjugating quark -> antiquark */
  Site site(localX,localY,localZ,localT);
  int src_parity = site.physX()+site.physY()+site.physZ()+site.physT();
  src_parity = src_parity % 2;

  int vol =  GJP.VolNodeSites();

  for(int nu=0;nu<4;nu++){
    sourcelink[nu].ZeroMatrix();
    contact[nu] = Zero;
    if (coor_x == procCoorX && coor_y == procCoorY && coor_z == procCoorZ && coor_t == procCoorT){
      Site site(localX,localY,localZ,localT);
      int id = site.Index();
      // first V_x over lattice, then V_y, etc.
      link = lat_pt->Fields(0) + nu*vol + id;
      sourcelink[nu] = *link;
#if 0
      for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
	  printf("link %d nu %d  %d %d  %g %g\n",id,nu,i,j,sourcelink[nu](i,j).real(),sourcelink[nu](i,j).imag());
	  //printf("link %d nu %d  %d %d  %g %g\n",id,nu,i,j,link[nu](i,j).real(),link[nu](i,j).imag());
	  //printf("link %d nu %d  %d %d  %g %g\n",id,nu,i,j,*((Float*)link+2*i+6*j),*((Float*)link+2*i+6*j+1));
#endif
    }
  }
#if 0
  {
    Site site;
    while ( site.LoopsOverNode() ){
      
      int id=site.Index();
      Matrix prop;
      for(int nu=0;nu<4;nu++){
	for(int color=0;color<3;color++)
	  moveFloat((Float*)&prop+6*color,quark[nu].AccessQuarkPropS(color,id),6);
	Matrix mat = *(lat_pt->Fields(0) + nu*vol + id);
	for(int i=0; i<3; i++)
	  for(int j=0; j<3; j++){
	    printf("Prop%d %d(%d,%d,%d,%d)  %d %d  %g %g\n",
		   nu,id,site.X(),site.Y(),site.Z(),site.T(),
		   i,j,prop(i,j).real(),prop(i,j).imag());
	    printf("Link %d(%d,%d,%d,%d) nu %d  %d %d  %g %g\n",
	       id,site.X(),site.Y(),site.Z(),site.T(),nu,
	       i,j,mat(i,j).real(),mat(i,j).imag());
	   }
	}
    }
  }
  //exit(0);
#endif  
  /* "broadcast" it */
  for(int nu=0;nu<4;nu++){
    slice_sum((Float*)&sourcelink[nu], 18, 99);
  }

  /* source dir */
  for(nu=0;nu<4;nu++){

    Matrix hsrclink;
    hsrclink.Dagger(sourcelink[nu]);

    /* sink dir */
    for(mu=0;mu<4;mu++){
              
      /* calculate the two-point function */
      Site site;
      while ( site.LoopsOverNode() ){

	int i=site.Index();
	int neighbor_site=site.plusIndex(mu);
	
	Matrix *link = lat_pt->Fields(0) + i + vol*mu;

	Matrix hlink;
	hlink.Dagger(*link);

	// get prop in +mu dir.
	// need the transpose (conjugate later)
	for(int color=0;color<3;color++){
	  if( /* offsite */ site.Coor(mu)+1 == n[mu] ){
	    getPlusData( (IFloat *)&propnu_neigh+6*color,
			 quark[nu].AccessQuarkPropS(color,neighbor_site), 6, mu);
	  } else {
	    moveFloat((Float*)&propnu_neigh+6*color, 
		      quark[nu].AccessQuarkPropS(color,neighbor_site),6);
	  }
	}
	Matrix temp;
	Matrix temp2;
	Matrix temp3;
	propnu_neigh.Transpose();
	temp.DotMEqual(*link,propnu_neigh); 
	temp2.Dagger(temp);
	Matrix prop;
	for(int color=0;color<3;color++)
	  moveFloat((Float*)&prop+6*color,quark[4].AccessQuarkPropS(color,i),6);
	prop.Transpose();
	temp.DotMEqual(temp2,prop);
	temp2.DotMEqual(temp,sourcelink[nu]);
	Rcomplex cc= temp2.Tr();
	
        /* phase from conjugating quark -> antiquark */
	int parity = site.physX()+site.physY()+site.physZ()+site.physT();
	parity = parity % 2;
        phase = (src_parity-parity==0 ? 1 : -1);
        cc *= phase;
        VacPol[nu+4*(mu+4*i)] = cc;

        // get prop in +mu dir.
        // need the transpose (conjugate later)
	for(int color=0;color<3;color++){
	  if( /* offsite */ site.Coor(mu)+1 == n[mu] ){
	    getPlusData( (IFloat *)&prop_neigh+6*color,
			 quark[4].AccessQuarkPropS(color,neighbor_site), 6, mu);
	  } else {
	    moveFloat((Float*)&prop_neigh+6*color, 
		      quark[4].AccessQuarkPropS(color,neighbor_site),6);
	  }
	}
	prop_neigh.Transpose();

	for(int color=0;color<3;color++)
	  moveFloat((Float*)&temp+6*color,quark[nu].AccessQuarkPropS(color,i),6);
	temp.Transpose();

        temp2.Dagger(temp);
        temp.DotMEqual(temp2,*link);
        temp2.DotMEqual(temp,prop_neigh);
        temp.DotMEqual(temp2,sourcelink[nu]);
        cc= temp.Tr();
        phase = (src_parity-parity==0 ? -1 : 1);
        cc *= phase;
        VacPol[nu+4*(mu+4*i)] += cc;
        
        temp.DotMEqual(*link,prop_neigh);
        temp2.Dagger(temp);
	for(int color=0;color<3;color++)
	  moveFloat((Float*)&temp+6*color,quark[nu].AccessQuarkPropS(color,i),6);
	temp.Transpose();
        temp3.DotMEqual(temp2,temp);
        temp2.DotMEqual(temp3,hsrclink);
        cc= temp2.Tr();
        phase = (src_parity-parity==0 ? -1 : 1);
        cc *= phase;
        VacPol[nu+4*(mu+4*i)] += cc;

	for(int color=0;color<3;color++)
	  moveFloat((Float*)&temp+6*color,quark[4].AccessQuarkPropS(color,i),6);
	temp.Transpose();
        temp2.Dagger(temp);
        temp.DotMEqual(temp2,*link);
        temp2.DotMEqual(temp,propnu_neigh);
        temp.DotMEqual(temp2,hsrclink);
        cc= temp.Tr();
        phase = (src_parity-parity==0 ? 1 : -1);
	cc *= phase;
        VacPol[nu+4*(mu+4*i)] += cc;

      } /* sites */
    } /* mu */

  
    /* contact terms */
    x[0] = sq_arg.src.origin[0] ;
    x[1] = sq_arg.src.origin[1] ;
    x[2] = sq_arg.src.origin[2] ;
    x[3] = sq_arg.src.origin[3] ;
    int x0 = (nu==0) ? (1+x[0])%(GJP.XnodeSites()*GJP.Xnodes()) : x[0];
    int y0 = (nu==1) ? (1+x[1])%(GJP.YnodeSites()*GJP.Ynodes()) : x[1];
    int z0 = (nu==2) ? (1+x[2])%(GJP.ZnodeSites()*GJP.Znodes()) : x[2];
    int t0 = (nu==3) ? (1+x[3])%(GJP.TnodeSites()*GJP.Tnodes()) : x[3];
    procCoorX = x0 / GJP.XnodeSites();
    procCoorY = y0 / GJP.YnodeSites();
    procCoorZ = z0 / GJP.ZnodeSites();
    procCoorT = t0 / GJP.TnodeSites();
    localX = x0 % GJP.XnodeSites();
    localY = y0 % GJP.YnodeSites();
    localZ = z0 % GJP.ZnodeSites();
    localT = t0 % GJP.TnodeSites();
    coor_x = GJP.XnodeCoor();
    coor_y = GJP.YnodeCoor();
    coor_z = GJP.ZnodeCoor();
    coor_t = GJP.TnodeCoor();
    if (coor_x == procCoorX && coor_y == procCoorY && coor_z == procCoorZ && coor_t == procCoorT){
      Site site(localX,localY,localZ,localT);
      i=site.Index();
      //mult_su3_nn(&(lattice[i].prop), &(sourcelink[nu]), &temp);
      for(int color=0;color<3;color++)
	moveFloat((Float*)&temp+6*color,quark[4].AccessQuarkPropS(color,i),6);
      temp.Transpose();
      temp2.DotMEqual(temp,sourcelink[nu]);
      cc = temp2.Tr();
       /* factor of 1/2 for current, 2 for each prop for
         field renorm. So, do nothing.*/
      contact[nu] = cc;
    }
    /* "broadcast" it */
    slice_sum((Float*)&contact[nu], 2, 99);

    x0 = sq_arg.src.origin[0] ;
    y0 = sq_arg.src.origin[1] ;
    z0 = sq_arg.src.origin[2] ;
    t0 = sq_arg.src.origin[3] ;
    procCoorX = x0 / GJP.XnodeSites();
    procCoorY = y0 / GJP.YnodeSites();
    procCoorZ = z0 / GJP.ZnodeSites();
    procCoorT = t0 / GJP.TnodeSites();
    localX = x0 % GJP.XnodeSites();
    localY = y0 % GJP.YnodeSites();
    localZ = z0 % GJP.ZnodeSites();
    localT = t0 % GJP.TnodeSites();
    coor_x = GJP.XnodeCoor();
    coor_y = GJP.YnodeCoor();
    coor_z = GJP.ZnodeCoor();
    coor_t = GJP.TnodeCoor();
    if (coor_x == procCoorX && coor_y == procCoorY && coor_z == procCoorZ && coor_t == procCoorT){
      Site site(localX,localY,localZ,localT);
      i=site.Index();
      //mult_su3_na(&(lattice[i].propnu[nu]), &(link[4*i+nu]), &temp);
      for(int color=0;color<3;color++)
	moveFloat((Float*)&temp+6*color,quark[nu].AccessQuarkPropS(color,i),6);
      temp.Transpose();
      Matrix *link = lat_pt->Fields(0) + i + vol*nu;
      Matrix hlink;
      hlink.Dagger(*link);
      temp2.DotMEqual(temp,hlink);

      cc = temp2.Tr();
      /* factor of 1/2 for current, 2 for each prop for
         field renorm. So, do nothing.*/
      contact[nu] -= cc;
      //VacPol[nu+4*(nu+4*i)] -= contact[nu];
    }
  } /* nu */
  etime = time_elapse();
  if(!UniqueID())printf("Time for HVP contraction %g secs\n",etime);

#if 0
  //print out PI
  Site mysite;
  while ( mysite.LoopsOverNode() ){
    
    int i=mysite.Index();
    
    for(nu=0;nu<4;nu++){
      for(mu=0;mu<4;mu++){
	if(!UniqueID())printf("PI %d %d  x %d %d %d %d   %e %e\n",
			      mu,nu,
			      mysite.physX(),mysite.physY(),mysite.physZ(),mysite.physT(),
			      VacPol[nu+4*(mu+4*i)].real(),VacPol[nu+4*(mu+4*i)].imag()
			      );
      }
    }
  }
#endif

#if 1

  //print out divergence
  {
    Site mysite;
    while ( mysite.LoopsOverNode() ){
      
      int i=mysite.Index();
      
      for(nu=0;nu<4;nu++){
	Rcomplex Div(0.,0.);
	for(mu=0;mu<4;mu++){

	  int i_neigh=mysite.minusIndex(mu);

	  Rcomplex cc;
	  if( /* offsite */ mysite.Coor(mu)-1 < 0 ){
	    getMinusData( (IFloat *)&cc, (Float*)&VacPol[nu+4*(mu+4*i_neigh)], 2, mu);
	  } else {
	      cc = VacPol[nu+4*(mu+4*i_neigh)];
	  }

	  Div += VacPol[nu+4*(mu+4*i)]-cc;

	}

	printf("DIV-VP nu %d  x %d %d %d %d   %e %e\n",
	       nu,
	       mysite.physX(),mysite.physY(),mysite.physZ(),mysite.physT(),
	       Div.real(),Div.imag()
	       );
      }
    }
  }
#endif



  // Fourier transform. "size" is number of complex
  // "1" last arg. is for "forward" fft, i.e. exp -ipx
  // to match our milc code
  FFT4((Float*)VacPol, 1, 2*size, 1, 1);
  etime = time_elapse();
  if(!UniqueID())printf("Time for HVP FFT %g secs\n",etime);

  quark[0].destroyQuarkPropS(sq_arg0.qid);
  quark[1].destroyQuarkPropS(sq_arg1.qid);
  quark[2].destroyQuarkPropS(sq_arg2.qid);
  quark[3].destroyQuarkPropS(sq_arg3.qid);
  quark[4].destroyQuarkPropS(sq_arg.qid);

}




//------------------------------------------------------------------
// Conserved-conserved vector current 2 pt function
//------------------------------------------------------------------
void AlgVacPolStag::VacPolStagConsCons()
{

  char *fname = "VacPolStagConsCons()";
  VRB.Func(cname,fname);
  
  int n[4], x[4], x_n[4];
  
  int i,j,k;
  int source_i;
  int icol,cgn;
  int mu, nu;
  Rcomplex contact[4];
  Rcomplex vacpol0;
  Float x1,x2;
  Float phase;
  Float finalrsq;
  Rcomplex cc;
  Matrix temp, temp2;
  Matrix sourcelink[4];
  Matrix *link;
  Matrix propnu_neigh;
  Matrix prop_neigh;
  
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();

  int size = 4*4; // 4 dirs * 4 dirs
  
  Rcomplex I(0.0,1.0);
  Rcomplex Zero(0.0,0.0);
 
  Rcomplex *VacPol;
  VacPol = (Rcomplex*)smalloc(4*4*GJP.VolNodeSites()*sizeof(Rcomplex));

  //loop over source location
  //for(int xs=ptStart;xs<n[0]*GJP.Xnodes();xs+=ptINC){
  // for(int ys=ptStart;ys<n[1]*GJP.Ynodes();ys+=ptINC){
  //  for(int zs=ptStart;zs<n[2]*GJP.Znodes();zs+=ptINC){
  //   for(int ts=ptStart;ts<n[3]*GJP.Tnodes();ts+=ptINC){

  // The quark props
  //------------------------------------
  StagQuarkArg sq_arg0 = sq_arg;
  sq_arg0.src.origin[0]++;
  StagQuarkArg sq_arg1 = sq_arg;
  sq_arg1.src.origin[1]++;
  StagQuarkArg sq_arg2 = sq_arg;
  sq_arg2.src.origin[2]++;
  StagQuarkArg sq_arg3 = sq_arg;
  sq_arg3.src.origin[3]++;
  // for conserved point split current at source: 1 centered, 4 shifted
  QuarkPropS quark[5] = { 
    QuarkPropS (AlgLattice(), sq_arg0),
    QuarkPropS (AlgLattice(), sq_arg1),
    QuarkPropS (AlgLattice(), sq_arg2),
    QuarkPropS (AlgLattice(), sq_arg3),
    QuarkPropS (AlgLattice(), sq_arg) 
  };
  
  // get props
  Vector *sol;
  for(int nu=0; nu<=4; nu++){
    quark[nu].setupQuarkPropS();
    quark[nu].getQuarkPropS((char*)common_arg->results);
    for(int c=0;c<3;c++)
      sol = (Vector*)quark[nu].AccessQuarkPropS(c,0);
    //AlgLattice().Convert(CANONICAL,sol);
  }

  // later, the point-split correlation function 
  // will need the link at source point
  x[0] = sq_arg.src.origin[0] ;
  x[1] = sq_arg.src.origin[1] ;
  x[2] = sq_arg.src.origin[2] ;
  x[3] = sq_arg.src.origin[3] ;
  // set point source
  int procCoorX = x[0] / GJP.XnodeSites();
  int procCoorY = x[1] / GJP.YnodeSites();
  int procCoorZ = x[2] / GJP.ZnodeSites();
  int procCoorT = x[3] / GJP.TnodeSites();
  int localX = x[0] % GJP.XnodeSites();
  int localY = x[1] % GJP.YnodeSites();
  int localZ = x[2] % GJP.ZnodeSites();
  int localT = x[3] % GJP.TnodeSites();
  
  int coor_x = GJP.XnodeCoor();
  int coor_y = GJP.YnodeCoor();
  int coor_z = GJP.ZnodeCoor();
  int coor_t = GJP.TnodeCoor();

  int vol =  GJP.VolNodeSites();

  for(int nu=0;nu<4;nu++){
    sourcelink[nu].ZeroMatrix();
    contact[nu] = Zero;
    if (coor_x == procCoorX && coor_y == procCoorY && coor_z == procCoorZ && coor_t == procCoorT){
      Site site(localX,localY,localZ,localT);
      int id = site.Index();
      // first V_x over lattice, then V_y, etc.
      link = lat_pt->Fields(0) + nu*vol + id;
      sourcelink[nu] = *link;
#if 0
      for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
	  printf("link %d nu %d  %d %d  %g %g\n",id,nu,i,j,sourcelink[nu](i,j).real(),sourcelink[nu](i,j).imag());
	  //printf("link %d nu %d  %d %d  %g %g\n",id,nu,i,j,link[nu](i,j).real(),link[nu](i,j).imag());
	  //printf("link %d nu %d  %d %d  %g %g\n",id,nu,i,j,*((Float*)link+2*i+6*j),*((Float*)link+2*i+6*j+1));
#endif
    }
  }
#if 0
  {
    Site site;
    while ( site.LoopsOverNode() ){
      
      int id=site.Index();
      Matrix prop;
      for(int nu=0;nu<4;nu++){
	for(int color=0;color<3;color++)
	  moveFloat((Float*)&prop+6*color,quark[nu].AccessQuarkPropS(color,id),6);
	Matrix mat = *(lat_pt->Fields(0) + nu*vol + id);
	for(int i=0; i<3; i++)
	  for(int j=0; j<3; j++){
	    printf("Prop%d %d(%d,%d,%d,%d)  %d %d  %g %g\n",
		   nu,id,site.X(),site.Y(),site.Z(),site.T(),
		   i,j,prop(i,j).real(),prop(i,j).imag());
	    printf("Link %d(%d,%d,%d,%d) nu %d  %d %d  %g %g\n",
	       id,site.X(),site.Y(),site.Z(),site.T(),nu,
	       i,j,mat(i,j).real(),mat(i,j).imag());
	   }
	}
    }
  }
  //exit(0);
#endif  
  /* "broadcast" it */
  for(int nu=0;nu<4;nu++){
    slice_sum((Float*)&sourcelink[nu], 18, 99);
  }

  /* source dir */
  for(nu=0;nu<4;nu++){

    Matrix hsrclink;
    hsrclink.Dagger(sourcelink[nu]);

    /* sink dir */
    for(mu=0;mu<4;mu++){
              
      /* calculate the two-point function */
      Site site;
      while ( site.LoopsOverNode() ){

	int i=site.Index();
	int neighbor_site=site.plusIndex(mu);
	
	Matrix *link = lat_pt->Fields(0) + i + vol*mu;

	Matrix hlink;
	hlink.Dagger(*link);

	// get prop in +mu dir.
	// need the transpose (conjugate later)
	for(int color=0;color<3;color++){
	  if( /* offsite */ site.Coor(mu)+1 == n[mu] ){
	    getPlusData( (IFloat *)&propnu_neigh+6*color,
			 quark[nu].AccessQuarkPropS(color,neighbor_site), 6, mu);
	  } else {
	    moveFloat((Float*)&propnu_neigh+6*color, 
		      quark[nu].AccessQuarkPropS(color,neighbor_site),6);
	  }
	}
	Matrix temp;
	Matrix temp2;
	Matrix temp3;
	propnu_neigh.Transpose();
	temp.DotMEqual(*link,propnu_neigh); 
	temp2.Dagger(temp);
	Matrix prop;
	for(int color=0;color<3;color++)
	  moveFloat((Float*)&prop+6*color,quark[4].AccessQuarkPropS(color,i),6);
	prop.Transpose();
	temp.DotMEqual(temp2,prop);
	temp2.DotMEqual(temp,sourcelink[nu]);
	Rcomplex cc= temp2.Tr();
	
        /* phase from conjugating quark -> antiquark */
	int parity = site.physX()+site.physY()+site.physZ()+site.physT();
	parity = parity % 2;
        phase = (parity==0 ? 1 : -1);
        cc *= phase;
        VacPol[nu+4*(mu+4*i)] = cc;

        // get prop in +mu dir.
        // need the transpose (conjugate later)
	for(int color=0;color<3;color++){
	  if( /* offsite */ site.Coor(mu)+1 == n[mu] ){
	    getPlusData( (IFloat *)&prop_neigh+6*color,
			 quark[4].AccessQuarkPropS(color,neighbor_site), 6, mu);
	  } else {
	    moveFloat((Float*)&prop_neigh+6*color, 
		      quark[4].AccessQuarkPropS(color,neighbor_site),6);
	  }
	}
	prop_neigh.Transpose();

	for(int color=0;color<3;color++)
	  moveFloat((Float*)&temp+6*color,quark[nu].AccessQuarkPropS(color,i),6);
	temp.Transpose();

        temp2.Dagger(temp);
        temp.DotMEqual(temp2,*link);
        temp2.DotMEqual(temp,prop_neigh);
        temp.DotMEqual(temp2,sourcelink[nu]);
        cc= temp.Tr();
        phase = (parity==0 ? -1 : 1);
        cc *= phase;
        VacPol[nu+4*(mu+4*i)] += cc;
        
        temp.DotMEqual(*link,prop_neigh);
        temp2.Dagger(temp);
	for(int color=0;color<3;color++)
	  moveFloat((Float*)&temp+6*color,quark[nu].AccessQuarkPropS(color,i),6);
	temp.Transpose();
        temp3.DotMEqual(temp2,temp);
        temp2.DotMEqual(temp3,hsrclink);
        cc= temp2.Tr();
        phase = (parity==0 ? -1 : 1);
        cc *= phase;
        VacPol[nu+4*(mu+4*i)] += cc;

	for(int color=0;color<3;color++)
	  moveFloat((Float*)&temp+6*color,quark[4].AccessQuarkPropS(color,i),6);
	temp.Transpose();
        temp2.Dagger(temp);
        temp.DotMEqual(temp2,*link);
        temp2.DotMEqual(temp,propnu_neigh);
        temp.DotMEqual(temp2,hsrclink);
        cc= temp.Tr();
        phase = (parity==0 ? 1 : -1);
	cc *= phase;
        VacPol[nu+4*(mu+4*i)] += cc;

      } /* sites */
    } /* mu */

    
    /* contact terms */
    x[0] = sq_arg.src.origin[0] ;
    x[1] = sq_arg.src.origin[1] ;
    x[2] = sq_arg.src.origin[2] ;
    x[3] = sq_arg.src.origin[3] ;
    int x0 = (nu==0) ? (1+x[0])%(GJP.XnodeSites()*GJP.Xnodes()) : x[0];
    int y0 = (nu==1) ? (1+x[1])%(GJP.YnodeSites()*GJP.Ynodes()) : x[1];
    int z0 = (nu==2) ? (1+x[2])%(GJP.ZnodeSites()*GJP.Znodes()) : x[2];
    int t0 = (nu==3) ? (1+x[3])%(GJP.TnodeSites()*GJP.Tnodes()) : x[3];
    procCoorX = x0 / GJP.XnodeSites();
    procCoorY = y0 / GJP.YnodeSites();
    procCoorZ = z0 / GJP.ZnodeSites();
    procCoorT = t0 / GJP.TnodeSites();
    localX = x0 % GJP.XnodeSites();
    localY = y0 % GJP.YnodeSites();
    localZ = z0 % GJP.ZnodeSites();
    localT = t0 % GJP.TnodeSites();
    coor_x = GJP.XnodeCoor();
    coor_y = GJP.YnodeCoor();
    coor_z = GJP.ZnodeCoor();
    coor_t = GJP.TnodeCoor();
    if (coor_x == procCoorX && coor_y == procCoorY && coor_z == procCoorZ && coor_t == procCoorT){
      Site site(localX,localY,localZ,localT);
      i=site.Index();
      //mult_su3_nn(&(lattice[i].prop), &(sourcelink[nu]), &temp);
      for(int color=0;color<3;color++)
	moveFloat((Float*)&temp+6*color,quark[4].AccessQuarkPropS(color,i),6);
      temp.Transpose();
      temp2.DotMEqual(temp,sourcelink[nu]);
      cc = temp2.Tr();
      /* factor of 1/2 for current, 2 for each prop for
         field renorm. So, do nothing.*/
      contact[nu] = cc;
    }
    /* "broadcast" it */
    slice_sum((Float*)&contact[nu], 2, 99);

    x0 = sq_arg.src.origin[0] ;
    y0 = sq_arg.src.origin[1] ;
    z0 = sq_arg.src.origin[2] ;
    t0 = sq_arg.src.origin[3] ;
    procCoorX = x0 / GJP.XnodeSites();
    procCoorY = y0 / GJP.YnodeSites();
    procCoorZ = z0 / GJP.ZnodeSites();
    procCoorT = t0 / GJP.TnodeSites();
    localX = x0 % GJP.XnodeSites();
    localY = y0 % GJP.YnodeSites();
    localZ = z0 % GJP.ZnodeSites();
    localT = t0 % GJP.TnodeSites();
    coor_x = GJP.XnodeCoor();
    coor_y = GJP.YnodeCoor();
    coor_z = GJP.ZnodeCoor();
    coor_t = GJP.TnodeCoor();
    if (coor_x == procCoorX && coor_y == procCoorY && coor_z == procCoorZ && coor_t == procCoorT){
      Site site(localX,localY,localZ,localT);
      i=site.Index();
      //mult_su3_na(&(lattice[i].propnu[nu]), &(link[4*i+nu]), &temp);
      for(int color=0;color<3;color++)
	moveFloat((Float*)&temp+6*color,quark[nu].AccessQuarkPropS(color,i),6);
      temp.Transpose();
      Matrix *link = lat_pt->Fields(0) + i + vol*nu;
      Matrix hlink;
      hlink.Dagger(*link);
      temp2.DotMEqual(temp,hlink);

      cc = temp2.Tr();
      /* factor of 1/2 for current, 2 for each prop for
         field renorm. So, do nothing.*/
      contact[nu] -= cc;
      VacPol[nu+4*(nu+4*i)] -= contact[nu];
    }
  } /* nu */

  // Fourier transform. "size" is number of complex
  // "1" last arg. is for "forward" fft, i.e. exp -ipx
  // to match our milc code
  FFT4((Float*)VacPol, 1, 2*size, 1, 1);

 // multiply by extra phases from source, point-split sink
 //-------------------------------------------------------
 
 const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
 const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
 const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
 const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();
 //FileIoType T=ADD_ID;
 char filename[128];
 sprintf(filename,"hvp-cons-cons.dat");
 FILE *fp=fopen(filename, "w");
 Float q[4];
 Float L[4];
 L[0] = n[0]*GJP.Xnodes();
 L[1] = n[1]*GJP.Ynodes();
 L[2] = n[2]*GJP.Znodes();
 L[3] = n[3]*GJP.Tnodes();
 
 for(x[3] = 0; x[3] < n[3]; x[3]++){
   int ptg = x[3]+shiftT;
   for(x[2] = 0; x[2] < n[2]; x[2]++){
     int pzg = x[2]+shiftZ;
     for(x[1] = 0; x[1] < n[1]; x[1]++){
       int pyg = x[1]+shiftY;
       for(x[0] = 0; x[0] < n[0]; x[0]++){
	 int pxg = x[0]+shiftX;
	
if(pxg > 6 || pyg > 6 || pzg > 6 || ptg > 16 ) continue;
 
	 int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
	 
	 for(nu=0;nu<4;nu++){
	   // point-split sink phase
	   q[0] = nu==0 ? PI*pxg/L[0] : 0;
	   q[1] = nu==1 ? PI*pyg/L[1] : 0;
	   q[2] = nu==2 ? PI*pzg/L[2] : 0; 
	   q[3] = nu==3 ? PI*ptg/L[3] : 0;
	   Float theta_src = (q[0]+q[1]+q[2]+q[3]);
	   for(mu=0;mu<4;mu++){
	     // point-split sink phase
	     q[0] = mu==0 ? PI*pxg/L[0] : 0;
	     q[1] = mu==1 ? PI*pyg/L[1] : 0;
	     q[2] = mu==2 ? PI*pzg/L[2] : 0; 
	     q[3] = mu==3 ? PI*ptg/L[3] : 0;
	     Float theta = (q[0]+q[1]+q[2]+q[3]);
	     Rcomplex phase(cos(theta_src-theta),sin(theta_src-theta));
	     VacPol[nu+4*(mu+4*site)] *= phase;
	     printf("VACPOL0 %d %d  mom %d %d %d %d  %e %e\n",
		     mu, nu,
		     pxg, pyg, pzg, ptg,
		     VacPol[nu+4*(mu+4*site)].real(),VacPol[nu+4*(mu+4*site)].imag()
		     );
	   }
	 }
       }
     }
   }
 }
 fclose(fp);

 /*int conf = alg_muon_arg->conf;
 sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.xyz%d%d%d.%d",
	 DIRVP,VPTAG,
	 alg_muon_arg->loop_mass,
	 t_op,
	 qp_arg->x,
	 qp_arg->y,
	 qp_arg->z,
	 alg_muon_arg->conf);
 qio_writeMuon writeQio(argc,argv);
 writeQio.write(filename,AlgLattice(),VacPol,size);
 sfree(sol);*/

 sfree(VacPol);

}



CPS_END_NAMESPACE
