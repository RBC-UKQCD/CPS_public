//------------------------------------------------------------------
//
// alg_muon.C
//
// T. Blum (June 2007)
//
// routines to calculate 2 and 3 pt muon correlation functions
//------------------------------------------------------------------

#include <stdlib.h>	// exit()
#include <stdio.h>
#include <sys/stat.h> 
#include <alg/alg_muon.h>
#include <util/qio_readMuon.h>
#include <util/qio_writeMuon.h>
#include <util/qcdio.h>
#include <util/momentum.h>
#include <comms/glb.h>
#include <comms/scu.h>
#include <alg/alg_eig.h>
#include <alg/eig_arg.h>
#include <util/lat_cont.h>

#define PI 3.141592654

bool FileExists(char * strFilename);

CPS_START_NAMESPACE

// stupid f..cking argc and argv for QMP!
int argc=0;
char *argv[1];

int NPROJ;
Float EPS; 
int XMOM;
int YMOM;
int ZMOM;
int NC;
int NHITS;
char* DIRVP; 
char* DIRML;
char* VPTAG;
int DO_MUONLINE; 
int DO_VACPOL; 
int tINC;
char* EIGTAG;
int ptINC;
int ptStart;
int Nmom;
int MaxMomSq;
ThreeMom *ext_mom;

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgMuon::AlgMuon(Lattice& latt, CommonArg *c_arg, MuonArg *arg, QPropWArg *qarg): 
  Alg(latt, c_arg)
{
  cname = "AlgMuon";
  char *fname = "AlgMuon(L&,CommonArg*,MuonArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0) ERR.Pointer(cname,fname, "arg");
  alg_muon_arg = arg;
  qp_arg = qarg;

  NPROJ=alg_muon_arg->NPROJ;
  EPS=alg_muon_arg->EPS; 
  XMOM=alg_muon_arg->XMOM;
  YMOM=alg_muon_arg->YMOM;
  ZMOM=alg_muon_arg->ZMOM;
  Nmom=alg_muon_arg->Nmom ;
  MaxMomSq=alg_muon_arg->MaxMomSq ;
  NC=alg_muon_arg->NConfs;
  NHITS=alg_muon_arg->NHITS;
  DIRVP=alg_muon_arg->DIRVP; 
  DIRML=alg_muon_arg->DIRML;
  VPTAG=alg_muon_arg->VPTAG;
  DO_MUONLINE=alg_muon_arg->DO_MUONLINE; 
  DO_VACPOL=alg_muon_arg->DO_VACPOL; 
  tINC=alg_muon_arg->tINC;
  EIGTAG=alg_muon_arg->EIGTAG;
  ptINC=alg_muon_arg->ptINC;
  ptStart=alg_muon_arg->ptStart;

  ext_mom = (ThreeMom*)smalloc(Nmom*sizeof(ThreeMom));
  int count(0);
  for(int p1=-XMOM;p1<=XMOM;p1++)
    for(int p2=-YMOM;p2<=YMOM;p2++)
      for(int p3=-ZMOM;p3<=ZMOM;p3++)
        if((p1*p1+p2*p2+p3*p3)<=MaxMomSq)
          //if((p1*p1+p2*p2+p3*p3)!=0)// eliminate the p=0
            if(count<Nmom)
             {
               if(!UniqueID())printf("External Momenta %d  %d %d %d\n",count,p1,p2,p3);
               ext_mom[count] = ThreeMom(p1,p2,p3);
               count++ ;
             }

  prop=NULL ;
  // must initialize
  argv[0]="dummy";
}

AlgMuon::AlgMuon(Lattice& latt, CommonArg *c_arg, MuonArg *arg, QPropWArg *qarg, EigCGArg *ecgarg): 
  Alg(latt, c_arg)
{
  cname = "AlgMuon";
  char *fname = "AlgMuon(L&,CommonArg*,MuonArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0) ERR.Pointer(cname,fname, "arg");
  alg_muon_arg = arg;
  qp_arg = qarg;
  eigcg_arg = ecgarg;

  NPROJ=alg_muon_arg->NPROJ;
  EPS=alg_muon_arg->EPS; 
  XMOM=alg_muon_arg->XMOM;
  YMOM=alg_muon_arg->YMOM;
  ZMOM=alg_muon_arg->ZMOM;
  Nmom=alg_muon_arg->Nmom ;
  MaxMomSq=alg_muon_arg->MaxMomSq ;
  NC=alg_muon_arg->NConfs;
  NHITS=alg_muon_arg->NHITS;
  DIRVP=alg_muon_arg->DIRVP; 
  DIRML=alg_muon_arg->DIRML;
  VPTAG=alg_muon_arg->VPTAG;
  DO_MUONLINE=alg_muon_arg->DO_MUONLINE; 
  DO_VACPOL=alg_muon_arg->DO_VACPOL; 
  tINC=alg_muon_arg->tINC;
  EIGTAG=alg_muon_arg->EIGTAG;
  ptINC=alg_muon_arg->ptINC;
  ptStart=alg_muon_arg->ptStart;

  ext_mom = (ThreeMom*)smalloc(Nmom*sizeof(ThreeMom));
  int count(0);
  for(int p1=-XMOM;p1<=XMOM;p1++)
    for(int p2=-YMOM;p2<=YMOM;p2++)
      for(int p3=-ZMOM;p3<=ZMOM;p3++)
        if((p1*p1+p2*p2+p3*p3)<=MaxMomSq)
          //if((p1*p1+p2*p2+p3*p3)!=0)// eliminate the p=0
            if(count<Nmom)
             {
               if(!UniqueID())printf("External Momenta %d  %d %d %d\n",count,p1,p2,p3);
               ext_mom[count] = ThreeMom(p1,p2,p3);
               count++ ;
             }

  prop=NULL ;
  // must initialize
  argv[0]="dummy";
}

//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgMuon::~AlgMuon() {
  char *fname = "~AlgMuon()";
  VRB.Func(cname,fname);
  //delete prop if it is still alive
  if(prop!=NULL) delete prop ;
  sfree(ext_mom);
}

//------------------------------------------------------------------
// local-conserved
//------------------------------------------------------------------
void AlgMuon::VacPolConsLoc(int t_op, Rcomplex* VacPol)
{

 char *fname = "VacPolConsLoc(int t_op, Rcomplex* VacPol)";
 VRB.Func(cname,fname);

 // assumes NO spreadout Ls
 if(GJP.Snodes() != 1){
   ERR.General(cname,fname,"Dummy: Snodes must be 1!\n");
 }

 int i, j;
 int mu, nu;
 int n[5], x[5], x_n[5];

 n[0] = GJP.XnodeSites();
 n[1] = GJP.YnodeSites();
 n[2] = GJP.ZnodeSites();
 n[3] = GJP.TnodeSites();
 n[4] = GJP.SnodeSites();
 int size = 4*4; // 4 dirs * 4 dirs

 // The lepton prop
 //------------------------------------
 qp_arg->cg.mass=alg_muon_arg->loop_mass;
 qp_arg->t=t_op;
 // these should be set properly in vml file
 //qp_arg->x =  0;
 //qp_arg->y =  0;
 //qp_arg->z =  0;
 // use midpoint prop to store 4d prop for 2nd s-slice
 qp_arg->store_midprop = 1;
 //Float save = GJP.DwfHeight();
 //GJP.SetDwfHeight(1.8);
 // This needs to be changed to an all-to-all

 QPropWPointSrc lepton(AlgLattice(), qp_arg, common_arg);
 //GJP.SetDwfHeight(save);
 
 int fsize = 3*4;
 Rcomplex *sol;
 sol = (Rcomplex*)smalloc(fsize*GJP.VolNodeSites()*sizeof(Rcomplex));
 char filename[100];
 Rcomplex I(0.0,1.0);
 Rcomplex Zero(0.0,0.0);
 int ls = GJP.SnodeSites();

 // initialize to zero since we are accumulating sum over s
 for(mu=0; mu< 4; mu++){
   for(nu=0; nu< 4; nu++){
     for (int site=0; site<GJP.VolNodeSites(); site++) {
       VacPol[nu+4*(mu+4*site)] = Zero;
     }
   }
 }

 for(x[4] = 0; x[4] < n[4]; x[4]++){
   int s=x[4];
   for(mu=0; mu< 4; mu++){//source gamma
     for(nu=0; nu< 4; nu++){//sink gamma (conserved current)
       for(x[3] = 0; x[3] < n[3]; x[3]++){
	 for(x[2] = 0; x[2] < n[2]; x[2]++){
	   for(x[1] = 0; x[1] < n[1]; x[1]++){
	     for(x[0] = 0; x[0] < n[0]; x[0]++){
	       
	       // coordinates for neighbor
	       for (i = 0; i < 4; i++ )
		 x_n[i] = (i == nu) ? (x[i]+1)%n[i] : x[i];               // +hop check
	       
	       // offsets
	       int site = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3])));
               int neighbor_site = x_n[0]+n[0]*(x_n[1]+n[1]*(x_n[2]+n[2]*(x_n[3])));
	       
	       // U_mu(x) where mu = prop_dir
	       Matrix *link = AlgLattice().GaugeField() + 4*site + nu;
               Matrix tempmat;
               tempmat.Dagger(*link);                                     // check dag
	       
               // mult prop on source by gamma_mu * gamma_5
               WilsonMatrix mat = lepton(s,site);                         // s check
               mat.gr(mu).gr(-5);                                         // gamma check
	       
	       WilsonMatrix wmat_neigh;
	       
	       // get prop in +nu dir.
	       // need the transpose (conjugate later)
	       if( /* offsite */ x[nu]+1 == n[nu] ){
		 
		 getPlusData( (IFloat *)&wmat_neigh,
			      (IFloat *)&lepton(ls-1-s,neighbor_site),    // s,x check
			      288, nu);
	       } else {
		 
		 wmat_neigh = lepton(ls-1-s,neighbor_site);               // s,x check
	       }
	       
	       // mult on sink by gamma_5 (1+gamma_nu) U^dagger_mu
	       WilsonMatrix temp = mat;
	       WilsonMatrix temp2 = temp;
	       temp2 += temp.gl(nu);
	       temp2.gl(-5);                                              // gamma check
	       temp2.LeftTimesEqual(tempmat);                             // link check
	       // spin-color trace
	       WilsonMatrix cc = wmat_neigh;
	       cc.hconj();                                                // conj check
	       VacPol[nu+4*(mu+4*site)] += 0.5*Trace(cc,temp2);           // tr. check

	       
               // the other contraction...
	       
               // mult prop on source by gamma_mu * gamma_5
               mat = lepton(ls-1-s,site);
               mat.gr(-5).gr(mu); // check gamma (conj below)
	       
	       // get the fields in plus direction
	       if( /* offsite */ x[nu]+1 == n[nu] ){
		 
		 getPlusData( (IFloat *)&wmat_neigh,
			      (IFloat *)&lepton(s,neighbor_site), 
			      288, nu);
	       } else {
		 
		 wmat_neigh = lepton(s,neighbor_site);               // check s,x
	       }
	       
	       temp = wmat_neigh;
	       temp2 = temp;
	       temp2 -= temp.gl(nu);
	       temp2.gl(-5);
	       temp2.LeftTimesEqual(*link);                     // check gamma, U_mu
		 
	       // spin-color trace
	       cc=mat;
	       cc.hconj();                                      // check s,x, conj
	       VacPol[nu+4*(mu+4*site)] -=  0.5*Trace(cc, temp2);
             }
           }
         }
       }
     }
   }
 }

 // Fourier transform. "size" is number of complex
 // "0" last arg. is for "backward" fft, i.e. exp +ipx
 // which we must do for either the loop or the muon line, so 
 // do it here. i.e. we have to have exp i p(x-y)
 FFT4((Float*)VacPol, 1, 2*size, 1, 0);

 // multiply by extra phases from source, point-split sink
 //-------------------------------------------------------
 
 const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
 const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
 const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
 const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();
 //FileIoType T=ADD_ID;
 //sprintf(filename,"vp-cons-loc.t%d.dat",t_op);
 //FILE *fp=Fopen(T, filename, "w");
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
	 

	 int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
	 
	 for(nu=0;nu<4;nu++){
	   for(mu=0;mu<4;mu++){
	     // point-split sink phase
	     q[0] = nu==0 ? PI*pxg/L[0] : 0;
	     q[1] = nu==1 ? PI*pyg/L[1] : 0;
	     q[2] = nu==2 ? PI*pzg/L[2] : 0; 
	     q[3] = nu==3 ? PI*ptg/L[3] : 0;
	     // source phase (opposite sign to sink (FFT) phase)
	     // could be this is actually needed to put muons on-shell!
	     //q[0] -= 2.*PI*pxg/L[0]*qp_arg->x;
	     //q[1] -= 2.*PI*pyg/L[1]*qp_arg->y;
	     //q[2] -= 2.*PI*pzg/L[2]*qp_arg->z;
	     //q[3] -= 2.*PI*ptg/L[3]*qp_arg->t;
	     Float theta = (q[0]+q[1]+q[2]+q[3]);
	     Rcomplex phase(cos(theta),sin(theta));
	     VacPol[nu+4*(mu+4*site)] *= phase;
	     //debug:
	     /*Rcomplex cc(0.,1.);
	     int fact = (ptg*pxg*pyg*pzg);
	     if(fact==0){fact=1;}
	     VacPol[nu+4*(mu+4*site)] = cc/fact;*/
	     //Fprintf(fp,"t_op= %d mu= %d nu= %d mom= %d %d %d %d %e %e\n", 
	//	     t_op, mu, nu, pxg, pyg, pzg, ptg,
	//	     VacPol[nu+4*(mu+4*site)].real(), 
	//	     VacPol[nu+4*(mu+4*site)].imag());
	   }
	 }
       }
     }
   }
 }
 //Fclose(T,fp);

 int conf = alg_muon_arg->conf;
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
 sfree(sol);

}


//------------------------------------------------------------------
// local-conserved with 1 twisted bc and 1 periodic bc prop
//------------------------------------------------------------------
void AlgMuon::VacPolConsLocTwistedBC(int t_op, Rcomplex* VacPol)
{

 char *fname = "VacPolConsLoc(int t_op, Rcomplex* VacPol)";
 VRB.Func(cname,fname);

 // assumes NO spreadout Ls
 if(GJP.Snodes() != 1){
   ERR.General(cname,fname,"Dummy: Snodes must be 1!\n");
 }

 int i, j;
 int mu, nu;
 int n[5], x[5], x_n[5];

 n[0] = GJP.XnodeSites();
 n[1] = GJP.YnodeSites();
 n[2] = GJP.ZnodeSites();
 n[3] = GJP.TnodeSites();
 n[4] = GJP.SnodeSites();
 int size = 4*4; // 4 dirs * 4 dirs

 // The lepton prop
 //------------------------------------
 qp_arg->cg.mass=alg_muon_arg->loop_mass;
 qp_arg->t=t_op;
 // these should be set properly in vml file
 //qp_arg->x =  0;
 //qp_arg->y =  0;
 //qp_arg->z =  0;
 // use midpoint prop to store 4d prop for 2nd s-slice
 qp_arg->store_midprop = 1;
 //Float save = GJP.DwfHeight();
 //GJP.SetDwfHeight(1.8);
 // This needs to be changed to an all-to-all

 QPropWPointSrc lepton(AlgLattice(), qp_arg, common_arg);
 if(!UniqueID()) printf("Twisting the lattice\n");
 AlgLattice().twist_bc(1);
 QPropWPointSrc Twlepton(AlgLattice(), qp_arg, common_arg);
 if(!UniqueID()) printf("Un-Twisting the lattice\n");
 AlgLattice().twist_bc(-1);

 int fsize = 3*4;
 Rcomplex *sol;
 sol = (Rcomplex*)smalloc(fsize*GJP.VolNodeSites()*sizeof(Rcomplex));
 char filename[100];
 Rcomplex I(0.0,1.0);
 Rcomplex Zero(0.0,0.0);
 int ls = GJP.SnodeSites();

 // initialize to zero since we are accumulating sum over s
 for(mu=0; mu< 4; mu++){
   for(nu=0; nu< 4; nu++){
     for (int site=0; site<GJP.VolNodeSites(); site++) {
       VacPol[nu+4*(mu+4*site)] = Zero;
     }
   }
 }

 for(x[4] = 0; x[4] < n[4]; x[4]++){
   int s=x[4];
   for(mu=0; mu< 4; mu++){//source gamma
     for(nu=0; nu< 4; nu++){//sink gamma (conserved current)
       for(x[3] = 0; x[3] < n[3]; x[3]++){
	 for(x[2] = 0; x[2] < n[2]; x[2]++){
	   for(x[1] = 0; x[1] < n[1]; x[1]++){
	     for(x[0] = 0; x[0] < n[0]; x[0]++){
	       
	       // coordinates for neighbor
	       for (i = 0; i < 4; i++ )
		 x_n[i] = (i == nu) ? (x[i]+1)%n[i] : x[i];               // +hop check
	       
	       // offsets
	       int site = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3])));
               int neighbor_site = x_n[0]+n[0]*(x_n[1]+n[1]*(x_n[2]+n[2]*(x_n[3])));
	       
	       // U_mu(x) where mu = prop_dir
	       Matrix *link = AlgLattice().GaugeField() + 4*site + nu;
               Matrix tempmat;
               tempmat.Dagger(*link);                                     // check dag
	       
               // mult prop on source by gamma_mu * gamma_5
               WilsonMatrix mat = Twlepton(s,site);                         // s check
               mat.gr(mu).gr(-5);                                         // gamma check
	       
	       WilsonMatrix wmat_neigh;
	       
	       // get prop in +nu dir.
	       // need the transpose (conjugate later)
	       if( /* offsite */ x[nu]+1 == n[nu] ){
		 
		 getPlusData( (IFloat *)&wmat_neigh,
			      (IFloat *)&lepton(ls-1-s,neighbor_site),    // s,x check
			      288, nu);
	       } else {
		 
		 wmat_neigh = lepton(ls-1-s,neighbor_site);               // s,x check
	       }
	       
	       // mult on sink by gamma_5 (1+gamma_nu) U^dagger_mu
	       WilsonMatrix temp = mat;
	       WilsonMatrix temp2 = temp;
	       temp2 += temp.gl(nu);
	       temp2.gl(-5);                                              // gamma check
	       temp2.LeftTimesEqual(tempmat);                             // link check
	       // spin-color trace
	       WilsonMatrix cc = wmat_neigh;
	       cc.hconj();                                                // conj check
	       VacPol[nu+4*(mu+4*site)] += 0.5*Trace(cc,temp2);           // tr. check

	       
               // the other contraction...
	       
               // mult prop on source by gamma_mu * gamma_5
               mat = lepton(ls-1-s,site);
               mat.gr(-5).gr(mu); // check gamma (conj below)
	       
	       // get the fields in plus direction
	       if( /* offsite */ x[nu]+1 == n[nu] ){
		 
		 getPlusData( (IFloat *)&wmat_neigh,
			      (IFloat *)&Twlepton(s,neighbor_site), 
			      288, nu);
	       } else {
		 
		 wmat_neigh = Twlepton(s,neighbor_site);               // check s,x
	       }
	       
	       temp = wmat_neigh;
	       temp2 = temp;
	       temp2 -= temp.gl(nu);
	       temp2.gl(-5);
	       temp2.LeftTimesEqual(*link);                     // check gamma, U_mu
		 
	       // spin-color trace
	       cc=mat;
	       cc.hconj();                                      // check s,x, conj
	       VacPol[nu+4*(mu+4*site)] -=  0.5*Trace(cc, temp2);
             }
           }
         }
       }
     }
   }
 }


#if 1
  //print out divergence
  {
    Site mysite;
    while ( mysite.LoopsOverNode() ){
      
      int i=mysite.Index();
      
      for(mu=0;mu<4;mu++){
	Rcomplex Div(0.,0.);
	for(nu=0;nu<4;nu++){
	  int i_neigh=mysite.minusIndex(nu);
	  Div += VacPol[nu+4*(mu+4*i)]-VacPol[nu+4*(mu+4*i_neigh)];
	}
	if(!UniqueID())printf("DIV-VP mu %d  x %d %d %d %d   %e %e\n",
			      mu,
			      mysite.physX(),mysite.physY(),mysite.physZ(),mysite.physT(),
			      Div.real(),Div.imag()
			      );
      }
    }
  }
#endif



 // Fourier transform. "size" is number of complex
 // "0" last arg. is for "backward" fft, i.e. exp +ipx
 // which we must do for either the loop or the muon line, so 
 // do it here. i.e. we have to have exp i p(x-y)
 FFT4((Float*)VacPol, 1, 2*size, 1, 0);

 // multiply by extra phases from source, point-split sink
 //-------------------------------------------------------
 
 const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
 const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
 const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
 const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();
 //FileIoType T=ADD_ID;
 //sprintf(filename,"vp-cons-loc.t%d.dat",t_op);
 //FILE *fp=Fopen(T, filename, "w");
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
	 

	 int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
	 
	 for(nu=0;nu<4;nu++){
	   for(mu=0;mu<4;mu++){
	     // point-split sink phase
	     q[0] = nu==0 ? PI*pxg/L[0] : 0;
	     q[1] = nu==1 ? PI*pyg/L[1] : 0;
	     q[2] = nu==2 ? PI*pzg/L[2] : 0; 
	     q[3] = nu==3 ? PI*ptg/L[3] : 0;
	     // source phase (opposite sign to sink (FFT) phase)
	     // could be this is actually needed to put muons on-shell!
	     //q[0] -= 2.*PI*pxg/L[0]*qp_arg->x;
	     //q[1] -= 2.*PI*pyg/L[1]*qp_arg->y;
	     //q[2] -= 2.*PI*pzg/L[2]*qp_arg->z;
	     //q[3] -= 2.*PI*ptg/L[3]*qp_arg->t;
	     Float theta = (q[0]+q[1]+q[2]+q[3]);
	     Rcomplex phase(cos(theta),sin(theta));
	     VacPol[nu+4*(mu+4*site)] *= phase;
	     //debug:
	     /*Rcomplex cc(0.,1.);
	     int fact = (ptg*pxg*pyg*pzg);
	     if(fact==0){fact=1;}
	     VacPol[nu+4*(mu+4*site)] = cc/fact;*/
	     //Fprintf(fp,"t_op= %d mu= %d nu= %d mom= %d %d %d %d %e %e\n", 
	//	     t_op, mu, nu, pxg, pyg, pzg, ptg,
	//	     VacPol[nu+4*(mu+4*site)].real(), 
	//	     VacPol[nu+4*(mu+4*site)].imag());
	   }
	 }
       }
     }
   }
 }
 //Fclose(T,fp);

 int conf = alg_muon_arg->conf;
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
 sfree(sol);

}

//------------------------------------------------------------------
// local-conserved, loop over point source to test spectral decomp.
// this only works on scalar machine for now.
//------------------------------------------------------------------
void AlgMuon::VacPolConsLocTest(int t_op, Rcomplex* VacPol, int* mom)
{

 char *fname = "VacPolConsLocTest(int t_op, Rcomplex* VacPol)";
 VRB.Func(cname,fname);

 // assumes NO spreadout Ls
 if(GJP.Snodes() != 1){
   ERR.General(cname,fname,"Dummy: Snodes must be 1!\n");
 }

 int i, j;
 int mu, nu, spin, color;
 int n[5], x[5], x_n[5];

 n[0] = GJP.XnodeSites();
 n[1] = GJP.YnodeSites();
 n[2] = GJP.ZnodeSites();
 n[3] = GJP.TnodeSites();
 n[4] = GJP.SnodeSites();
 int size = 4*4; // 4 dirs * 4 dirs

 Float q[4];
 Float L[4];
 L[0] = n[0]*GJP.Xnodes();
 L[1] = n[1]*GJP.Ynodes();
 L[2] = n[2]*GJP.Znodes();
 L[3] = n[3]*GJP.Tnodes();

 const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
 const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
 const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
 const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();

 char filename[100];
 Rcomplex I(0.0,1.0);
 Rcomplex Zero(0.0,0.0);
 int ls = GJP.SnodeSites();


 // initialize to zero since we are accumulating sum over s, 
 // and FT at source
 for(mu=0; mu< 4; mu++){
   for(nu=0; nu< 4; nu++){
     for (int site=0; site<GJP.VolNodeSites(); site++) {
       VacPol[nu+4*(mu+4*site)] = Zero;
     }
   }
 }
 // loop over point sources
 for(int xs=0;xs<n[0]*GJP.Xnodes();xs+=ptINC){
   for(int ys=0;ys<n[1]*GJP.Ynodes();ys+=ptINC){
     for(int zs=0;zs<n[2]*GJP.Znodes();zs+=ptINC){
       // The lepton prop
       //------------------------------------
       qp_arg->cg.mass=alg_muon_arg->loop_mass;
       qp_arg->t=t_op;
       qp_arg->x =  xs;
       qp_arg->y =  ys;
       qp_arg->z =  zs;

       // source phase (opposite sign to sink (FFT) phase)
       q[0] = 2.*PI/L[0]*mom[0]*qp_arg->x;
       q[1] = 2.*PI/L[1]*mom[1]*qp_arg->y;
       q[2] = 2.*PI/L[2]*mom[2]*qp_arg->z;
       Float theta = (q[0]+q[1]+q[2]);
       //theta=0.0; //debug
       Rcomplex phase(cos(theta),-sin(theta));
       //printf("MOM= %d %d %d %e %e %e %e %e\n",mom[0],mom[1],mom[2],q[0],q[1],q[2],phase.real(),phase.imag());
       // This needs to be changed to an all-to-all
       QPropWPointSrc lepton(AlgLattice(), qp_arg, common_arg);       
       
       color=0; // Pure QED for now
       // load in prop, one s-slice at a time
       // "spin" and "color" are source indices.
       // 5th dimension loop is cut in half since
       // sum is same for both halves (mult by 2)
       for(x[4] = 0; x[4] < n[4]; x[4]++){
	 int s=x[4];
	 for(mu=0; mu< 4; mu++){//source gamma
	   for(nu=0; nu< 4; nu++){//sink gamma (conserved current)
	     for(x[3] = 0; x[3] < n[3]; x[3]++){
	       for(x[2] = 0; x[2] < n[2]; x[2]++){
		 for(x[1] = 0; x[1] < n[1]; x[1]++){
		   for(x[0] = 0; x[0] < n[0]; x[0]++){
		     
		     // coordinates for neighbor
		     for (i = 0; i < 4; i++ )
		       x_n[i] = (i == nu) ? (x[i]+1)%n[i] : x[i];               // +hop check
		     
		     // offsets
		     int site = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3])));
		     int neighbor_site = x_n[0]+n[0]*(x_n[1]+n[1]*(x_n[2]+n[2]*(x_n[3])));
		     
		     // U_mu(x) where mu = prop_dir
		     Matrix *link = AlgLattice().GaugeField() + 4*site + nu;
		     Matrix tempmat;
		     tempmat.Dagger(*link);                                     // check dag
 		     
		     // mult prop on source by gamma_mu * gamma_5
		     WilsonMatrix mat = lepton(s,site);                         // s check
 		     mat.gr(mu).gr(-5);                                         // gamma check
		     
		     WilsonMatrix wmat_neigh;
		     
		     // get prop in +nu dir.
		     // need the transpose (conjugate later)
		     if( /* offsite */ x[nu]+1 == n[nu] ){
		       
		       getPlusData( (IFloat *)&wmat_neigh,
				    (IFloat *)&lepton(ls-1-s,neighbor_site),    // s,x check
				    288, nu);
		     } else {
		       
		       wmat_neigh = lepton(ls-1-s,neighbor_site);               // s,x check
		     }
		     
		     // mult on sink by gamma_5 (1+gamma_nu) U^dagger_mu
		     WilsonMatrix temp = mat;
		     WilsonMatrix temp2 = temp;
		     temp2 += temp.gl(nu);
		     temp2.gl(-5);                                              // gamma check
		     temp2.LeftTimesEqual(tempmat);                             // link check
		     // spin-color trace
		     WilsonMatrix cc = wmat_neigh;
		     cc.hconj();                                                // conj check
		     VacPol[nu+4*(mu+4*site)] += 0.5*Trace(cc,temp2)*phase;           // tr. check
		     
		     
		     // the other contraction...
		     
		     // mult prop on source by gamma_mu * gamma_5
		     mat = lepton(ls-1-s,site);
		     mat.gr(-5).gr(mu); // check gamma (conj below)
		     
		     // get the fields in plus direction
		     if( /* offsite */ x[nu]+1 == n[nu] ){
		       
		       getPlusData( (IFloat *)&wmat_neigh,
				    (IFloat *)&lepton(s,neighbor_site), 
				    288, nu);
		     } else {
		       
		       wmat_neigh = lepton(s,neighbor_site);               // check s,x
		     }
		     
		     temp = wmat_neigh;
		     temp2 = temp;
		     temp2 -= temp.gl(nu);
		     temp2.gl(-5);
		     temp2.LeftTimesEqual(*link);                     // check gamma, U_mu
		     
		     // spin-color trace
		     cc=mat;
		     cc.hconj();                                      // check s,x, conj
		     VacPol[nu+4*(mu+4*site)] -=  0.5*Trace(cc, temp2)*phase;
		   }
		 }
	       }
	     }
	   }
	 }
       }// s loop
     }
   }
 } // source loop end
       
 // Fourier transform. "size" is number of complex
 // "0" last arg. is for "backward" fft, i.e. exp +ipx
 // which we must do for either the loop or the muon line, so 
 // do it here. i.e. we have to have exp i p(x-y)
 //FFT4((Float*)VacPol, 1, 2*size, 1, 0);
 
      
 // multiply by extra phases from source, point-split sink
 //-------------------------------------------------------
 
 //FileIoType T=ADD_ID;
 //sprintf(filename,"vp-cons-loc-Test.t%d.dat",t_op);
 //FILE *fp=Fopen(T, filename, "w");
 
 for(x[3] = 0; x[3] < n[3]; x[3]++){
   int ptg = x[3]+shiftT;
   for(x[2] = 0; x[2] < n[2]; x[2]++){
     int pzg = x[2]+shiftZ;
     for(x[1] = 0; x[1] < n[1]; x[1]++){
       int pyg = x[1]+shiftY;
       for(x[0] = 0; x[0] < n[0]; x[0]++){
	 int pxg = x[0]+shiftX;
	 
	 
	 int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
	 
	 for(nu=0;nu<4;nu++){
	   for(mu=0;mu<4;mu++){
	     // point-split sink phase
	     q[0] = nu==0 ? PI*pxg/L[0] : 0;
	     q[1] = nu==1 ? PI*pyg/L[1] : 0;
	     q[2] = nu==2 ? PI*pzg/L[2] : 0; 
	     q[3] = nu==3 ? PI*ptg/L[3] : 0;
	     // source phase (opposite sign to sink (FFT) phase)
	     // could be this is actually needed to put muons on-shell!
	     //q[0] -= 2.*PI*pxg/L[0]*qp_arg->x;
	     //q[1] -= 2.*PI*pyg/L[1]*qp_arg->y;
	     //q[2] -= 2.*PI*pzg/L[2]*qp_arg->z;
	     //q[3] -= 2.*PI*ptg/L[3]*qp_arg->t;
	     Float theta = (q[0]+q[1]+q[2]+q[3]);
	     Rcomplex phase(cos(theta),sin(theta));
	     VacPol[nu+4*(mu+4*site)] *= phase;
	     //if(!UniqueID())printf("VP t_op= %d mu= %d nu= %d mom= %d %d %d %d %.15e %.15e\n", 
		     //t_op, mu, nu, pxg, pyg, pzg, ptg,
		     //VacPol[nu+4*(mu+4*site)].real(), 
		     //VacPol[nu+4*(mu+4*site)].imag());
	   }
	 }
       }
     }
   }
 }
 //Fclose(T,fp);

 int conf = alg_muon_arg->conf;
 sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.mom%d%d%d.dat.%d",
	 DIRVP,VPTAG,
	 alg_muon_arg->loop_mass,
	 t_op,
	 mom[0],mom[1],mom[2],
	 alg_muon_arg->conf);
 qio_writeMuon writeQio(argc,argv);
 writeQio.write(filename,AlgLattice(),VacPol,size);
}



//------------------------------------------------------------------
// local-conserved, loop over point source to test spectral decomp.
// this only works on scalar machine for now.
//------------------------------------------------------------------
void AlgMuon::VacPolConsLocLoopPtSrc(int t_op)
{

 char *fname = "VacPolConsLocLoopPtSrc(int t_op)";
 VRB.Func(cname,fname);

 // assumes NO spreadout Ls
 if(GJP.Snodes() != 1){
   ERR.General(cname,fname,"Dummy: Snodes must be 1!\n");
 }

 int i, j;
 int mu, nu, spin, color;
 int n[5], x[5], x_n[5];

 n[0] = GJP.XnodeSites();
 n[1] = GJP.YnodeSites();
 n[2] = GJP.ZnodeSites();
 n[3] = GJP.TnodeSites();
 n[4] = GJP.SnodeSites();
 int size = 4*4; // 4 dirs * 4 dirs

 Float q[4];
 Float L[4];
 L[0] = n[0]*GJP.Xnodes();
 L[1] = n[1]*GJP.Ynodes();
 L[2] = n[2]*GJP.Znodes();
 L[3] = n[3]*GJP.Tnodes();

 const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
 const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
 const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
 const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();

 char filename[100];
 Rcomplex I(0.0,1.0);
 Rcomplex Zero(0.0,0.0);
 int ls = GJP.SnodeSites();

 // external momenta * mu * nu * internal momenta (4-vol) 
 Rcomplex *ext_phase;
 ext_phase = (Rcomplex*)smalloc(Nmom*sizeof(Rcomplex));
 Rcomplex **VacPol = (Rcomplex **)smalloc(Nmom*sizeof(Rcomplex*));
 for(int p=0;p<Nmom;p++) VacPol[p] = (Rcomplex*)smalloc(16*GJP.VolNodeSites()*sizeof(Rcomplex));

 // initialize to zero since we are accumulating sum over s, 
 // and FT at source
 for(mu=0; mu< 4; mu++){
   for(nu=0; nu< 4; nu++){
     for (int site=0; site<GJP.VolNodeSites(); site++) {
       for(int n=0;n<Nmom;n++)VacPol[n][nu+4*(mu+4*site)] = Zero;
     }
   }
 }
#if 0
 int nev=eigcg_arg->nev;
 int m=eigcg_arg->m;
 Float max_eig=eigcg_arg->max_eig_cut;
 int max_def_len=eigcg_arg->max_def_len;
 int restart_len=eigcg_arg->restart_len;
 Float *restart=eigcg_arg->restart;
 bool always_restart=eigcg_arg->always_restart;
 
 int def_len=0;
 int vec_len=GJP.VolNodeSites() * AlgLattice().FsiteSize() / (AlgLattice().FchkbEvl()+1);
 Vector **V=new Vector*[m];
 for(int i=0;i<m;i++)
   V[i]=(Vector *)smalloc(cname,fname,"V",vec_len*sizeof(Float));
 Float *M=new Float[m];
 float **U=new float*[max_def_len];
 for(int i=0;i<max_def_len;i++)U[i]=(float *)smalloc(cname,fname,"U",vec_len*sizeof(float));
 Rcomplex *H=new Rcomplex[max_def_len*max_def_len];
#endif

 Float b_coeff = 0.5;
 if(AlgLattice().Fclass()==F_CLASS_MOBIUS) b_coeff *= GJP.Mobius_b();

 // loop over point sources
 int nsrc=0;
 for(int xs=ptStart;xs<n[0]*GJP.Xnodes();xs+=ptINC){
   for(int ys=ptStart;ys<n[1]*GJP.Ynodes();ys+=ptINC){
     for(int zs=ptStart;zs<n[2]*GJP.Znodes();zs+=ptINC){
       //if(zs>3*ptINC) continue;
       // The lepton prop
       //------------------------------------
       qp_arg->cg.mass=alg_muon_arg->loop_mass;
       qp_arg->t=t_op;
       qp_arg->x =  xs;
       qp_arg->y =  ys;
       qp_arg->z =  zs;
       QPropWPointSrc lepton(AlgLattice(), qp_arg, common_arg);
       /*if(nsrc==0){
	 lepton.eig_Run(V,vec_len,M,max_eig,nev,m,U,H,max_def_len,def_len,restart,restart_len,always_restart,0,1e-8);
	 cout<<"deflation length:"<<def_len<<endl;
	 for(int i=0;i<m;i++)sfree(V[i]);
	 delete [] V;
       }else{
	 //solver with restarts 
	 lepton.eig_Run(NULL,vec_len,M,max_eig,0,0,U,H,max_def_len,def_len,restart,restart_len,always_restart,0,1e-8);
       }*/
#if 0
       lepton.eig_Run(V,vec_len,M,max_eig,nev,m,U,H,max_def_len,def_len,restart,restart_len,always_restart,0,1e-8);
#endif
       // source phase (opposite sign to sink (FFT) phase)
       for(int im=0;im<Nmom;im++){
         q[0] = 2.*PI/L[0]*ext_mom[im].cmp(0)*qp_arg->x;
         q[1] = 2.*PI/L[1]*ext_mom[im].cmp(1)*qp_arg->y;
         q[2] = 2.*PI/L[2]*ext_mom[im].cmp(2)*qp_arg->z;
         Float theta = (q[0]+q[1]+q[2]);
         Rcomplex cc(cos(theta),-sin(theta));
         ext_phase[im]=cc;
         //printf("MOM= %d %d %d %e %e %e %e %e\n",mom[0],mom[1],mom[2],q[0],q[1],q[2],phase.real(),phase.imag());
       }
       
       color=0; // Pure QED for now
       // load in prop, one s-slice at a time
       // "spin" and "color" are source indices.
       // 5th dimension loop is cut in half since
       // sum is same for both halves (mult by 2)
       for(x[4] = 0; x[4] < n[4]; x[4]++){
	 int s=x[4];
	 Float coeff = 0.5*GJP.Mobius_c();
	 Float coeffp1 = coeff;
	 if(s==0) coeff *= -alg_muon_arg->loop_mass;
	 if(s==ls-1) coeffp1 *= -alg_muon_arg->loop_mass;
	 for(mu=0; mu< 4; mu++){//source gamma
	   for(nu=0; nu< 4; nu++){//sink gamma (conserved current)
	     for(x[3] = 0; x[3] < n[3]; x[3]++){
	       for(x[2] = 0; x[2] < n[2]; x[2]++){
		 for(x[1] = 0; x[1] < n[1]; x[1]++){
		   for(x[0] = 0; x[0] < n[0]; x[0]++){
		     
		     // coordinates for neighbor
		     for (i = 0; i < 4; i++ )
		       x_n[i] = (i == nu) ? (x[i]+1)%n[i] : x[i];               // +hop check
		     
		     // offsets
		     int site = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3])));
		     int neighbor_site = x_n[0]+n[0]*(x_n[1]+n[1]*(x_n[2]+n[2]*(x_n[3])));
		     
		     // U_nu(x) 
		     Matrix *link = AlgLattice().GaugeField() + 4*site + nu;
		     Matrix tempmat;
		     tempmat.Dagger(*link);                                     // check dag
 		     
		     // mult prop on source by gamma_mu * gamma_5
		     WilsonMatrix mat = lepton(s,site);                         // s check
 		     mat.gr(mu).gr(-5);                                         // gamma check
		     
		     WilsonMatrix wmat_neigh;
		     
		     // get prop in +nu dir.
		     // need the transpose (conjugate later)
		     if( /* offsite */ x[nu]+1 == n[nu] ){
		       
		       getPlusData( (IFloat *)&wmat_neigh,
				    (IFloat *)&lepton(ls-1-s,neighbor_site),    // s,x check
				    288, nu);
		     } else {
		       
		       wmat_neigh = lepton(ls-1-s,neighbor_site);               // s,x check
		     }
		     
		     // mult on sink by gamma_5 (1+gamma_nu) U^dagger_mu
		     WilsonMatrix temp = mat;
		     WilsonMatrix temp2 = temp;
		     temp2 += temp.gl(nu);
		     temp2.gl(-5);                                              // gamma check
		     temp2.LeftTimesEqual(tempmat);                             // link check
		     // spin-color trace
		     WilsonMatrix cc = wmat_neigh;
		     cc.hconj();                                                // conj check
		     for(int im=0;im<Nmom;im++)VacPol[im][nu+4*(mu+4*site)] += b_coeff*Trace(cc,temp2)*ext_phase[im];   // tr. check
		     
		     
		     // the other contraction...
		     
		     // mult prop on source by gamma_mu * gamma_5
		     mat = lepton(ls-1-s,site);
		     mat.gr(-5).gr(mu); // check gamma (conj below)
		     
		     // get the fields in plus direction
		     if( /* offsite */ x[nu]+1 == n[nu] ){
		       
		       getPlusData( (IFloat *)&wmat_neigh,
				    (IFloat *)&lepton(s,neighbor_site), 
				    288, nu);
		     } else {
		       
		       wmat_neigh = lepton(s,neighbor_site);               // check s,x
		     }
		     
		     temp = wmat_neigh;
		     temp2 = temp;
		     temp2 -= temp.gl(nu);
		     temp2.gl(-5);
		     temp2.LeftTimesEqual(*link);                     // check gamma, U_mu
		     
		     // spin-color trace
		     cc=mat;
		     cc.hconj();                                      // check s,x, conj
		     for(int im=0;im<Nmom;im++) VacPol[im][nu+4*(mu+4*site)] -=  b_coeff*Trace(cc, temp2)*ext_phase[im];

		     if(AlgLattice().Fclass()==F_CLASS_MOBIUS){

		       int sminus1 = (s-1+ls)%ls; // s-1 check
		       // mult prop on source by gamma_mu * gamma_5
		       mat = lepton(sminus1,site);
		       mat.gr(mu).gr(-5);
		       		       
		       // get prop in +nu dir.
		       // need the transpose (conjugate later)
		       if( /* offsite */ x[nu]+1 == n[nu] ){
			 
			 getPlusData( (IFloat *)&wmat_neigh,
				      (IFloat *)&lepton(ls-1-s,neighbor_site),
				      288, nu);
		       } else {
			 
			 wmat_neigh = lepton(ls-1-s,neighbor_site);
		       }
		       
		       // mult on sink by gamma_5 (1+gamma_nu) U^dagger_nu
		       temp = mat;
		       temp.gl(PR); // (1+g5)/2 check
		       temp2 = temp;
		       temp2 += temp.gl(nu);
		       temp2.gl(-5);
		       temp2.LeftTimesEqual(tempmat);
		       // spin-color trace
		       cc = wmat_neigh;
		       cc.hconj();
		       for(int im=0;im<Nmom;im++)VacPol[im][nu+4*(mu+4*site)] += coeff*Trace(cc,temp2)*ext_phase[im];
		       
		       
		       // the other contraction...
		       
		       // mult prop on source by gamma_mu * gamma_5
		       mat = lepton(ls-1-s,site);
		       mat.gr(-5).gr(mu); // check gamma (conj below)
		       
		       // get the fields in plus direction
		       if( /* offsite */ x[nu]+1 == n[nu] ){
			 
			 getPlusData( (IFloat *)&wmat_neigh,
				      (IFloat *)&lepton(sminus1,neighbor_site), 
				      288, nu);
		       } else {
			 
			 wmat_neigh = lepton(sminus1,neighbor_site);
		       }
		       
		       temp = wmat_neigh;
		       temp.gl(PR);
		       temp2 = temp;
		       temp2 -= temp.gl(nu);
		       temp2.gl(-5);
		       temp2.LeftTimesEqual(*link);                 
		       
		       // spin-color trace
		       cc=mat;
		       cc.hconj();                                  
		       for(int im=0;im<Nmom;im++) VacPol[im][nu+4*(mu+4*site)] -=  coeff*Trace(cc, temp2)*ext_phase[im];
		       
		       // other bit
		       
		       int splus1 = (s+1+ls)%ls; // s+1 check
		       // mult prop on source by gamma_mu * gamma_5
		       mat = lepton(splus1,site);
		       mat.gr(mu).gr(-5);
		       
		       // get prop in +nu dir.
		       // need the transpose (conjugate later)
		       if( /* offsite */ x[nu]+1 == n[nu] ){
			 
			 getPlusData( (IFloat *)&wmat_neigh,
				      (IFloat *)&lepton(ls-1-s,neighbor_site),
				      288, nu);
		       } else {
			 
			 wmat_neigh = lepton(ls-1-s,neighbor_site);
		       }
		       
		       // mult on sink by gamma_5 (1+gamma_nu) U^dagger_mu
		       temp = mat;
		       temp.gl(PL); // (1-g5)/2 check
		       temp2 = temp;
		       temp2 += temp.gl(nu);
		       temp2.gl(-5);
		       temp2.LeftTimesEqual(tempmat);
		       // spin-color trace
		       cc = wmat_neigh;
		       cc.hconj();
		       for(int im=0;im<Nmom;im++)VacPol[im][nu+4*(mu+4*site)] += coeffp1*Trace(cc,temp2)*ext_phase[im];
		       
		       
		       // the other contraction...
		       
		       // mult prop on source by gamma_mu * gamma_5
		       mat = lepton(ls-1-s,site);
		       mat.gr(-5).gr(mu); // check gamma (conj below)
		       
		       // get the fields in plus direction
		       if( /* offsite */ x[nu]+1 == n[nu] ){
			 
			 getPlusData( (IFloat *)&wmat_neigh,
				      (IFloat *)&lepton(splus1,neighbor_site), 
				      288, nu);
		       } else {
			 
			 wmat_neigh = lepton(splus1,neighbor_site);
		       }
		       
		       temp = wmat_neigh;
		       temp.gl(PL);
		       temp2 = temp;
		       temp2 -= temp.gl(nu);
		       temp2.gl(-5);
		       temp2.LeftTimesEqual(*link);                 
		       
		       // spin-color trace
		       cc=mat;
		       cc.hconj();                                  
		       for(int im=0;im<Nmom;im++) VacPol[im][nu+4*(mu+4*site)] -=  coeffp1*Trace(cc, temp2)*ext_phase[im];
		     }
#if 0
for(int im=0;im<Nmom;im++) 
printf("VP(mom %d) %d %d   %d %d %d %d %e %e\n",
        im, mu, nu, x[0],x[1],x[2],x[3],
        VacPol[im][nu+4*(mu+4*site)].real(),
        VacPol[im][nu+4*(mu+4*site)].imag());
#endif
		   }
		 }
	       }
	     }
	   }
	 }
       }// s loop
       nsrc++; //count srcs
     }
   }
 } // source loop end
 // free eigcg vectors
#if 0
 delete [] H;
 for(int iu=0;iu<max_def_len;iu++)sfree(U[iu]);
 delete [] U;
 delete [] M;
 for(int i=0;i<m;i++)sfree(V[i]);
 delete [] V;
#endif
      
 // Fourier transform. "size" is number of complex
 // "0" last arg. is for "backward" fft, i.e. exp +ipx
 // which we must do for either the loop or the muon line, so 
 // do it here. i.e. we have to have exp i p(x-y)
 for(int im=0;im<Nmom;im++)FFT4((Float*)VacPol[im], 1, 2*size, 1, 0);
 
      
 // multiply by extra phases from source, point-split sink
 //-------------------------------------------------------
 
 FileIoType T=ADD_ID;
 sprintf(filename,"vp-cons-loc-Test.t%d.dat",t_op);
 FILE *fp=Fopen(T, filename, "w");
 
 for(x[3] = 0; x[3] < n[3]; x[3]++){
   int ptg = x[3]+shiftT;
   for(x[2] = 0; x[2] < n[2]; x[2]++){
     int pzg = x[2]+shiftZ;
     for(x[1] = 0; x[1] < n[1]; x[1]++){
       int pyg = x[1]+shiftY;
       for(x[0] = 0; x[0] < n[0]; x[0]++){
	 int pxg = x[0]+shiftX;
	 
	 
	 int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
	 
	 for(nu=0;nu<4;nu++){
	   for(mu=0;mu<4;mu++){
	     // point-split sink phase
	     q[0] = nu==0 ? PI*pxg/L[0] : 0;
	     q[1] = nu==1 ? PI*pyg/L[1] : 0;
	     q[2] = nu==2 ? PI*pzg/L[2] : 0; 
	     q[3] = nu==3 ? PI*ptg/L[3] : 0;
	     // source phase (opposite sign to sink (FFT) phase)
	     // could be this is actually needed to put muons on-shell!
	     //q[0] -= 2.*PI*pxg/L[0]*qp_arg->x;
	     //q[1] -= 2.*PI*pyg/L[1]*qp_arg->y;
	     //q[2] -= 2.*PI*pzg/L[2]*qp_arg->z;
	     //q[3] -= 2.*PI*ptg/L[3]*qp_arg->t;
	     Float theta = (q[0]+q[1]+q[2]+q[3]);
	     Rcomplex phase(cos(theta),sin(theta));
	     for(int im=0;im<Nmom;im++) {
		VacPol[im][nu+4*(mu+4*site)] *= phase;
	        if(!UniqueID())printf("VP exmom %d t_op= %d mu= %d nu= %d mom= %d %d %d %d %.15e %.15e\n", 
		     		       im, t_op, mu, nu, pxg, pyg, pzg, ptg,
		                       VacPol[im][nu+4*(mu+4*site)].real(), 
		                       VacPol[im][nu+4*(mu+4*site)].imag());
	     }
	   }
	 }
       }
     }
   }
 }
 Fclose(T,fp);

 int conf = alg_muon_arg->conf;
 for(int im=0;im<Nmom;im++){
   sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.mom%d%d%d.dat.%d",
  	   DIRVP,VPTAG,
  	   alg_muon_arg->loop_mass,
  	   t_op,
  	   ext_mom[im].cmp(0),ext_mom[im].cmp(1),ext_mom[im].cmp(2),
  	   alg_muon_arg->conf);
   qio_writeMuon writeQio(argc,argv);
   writeQio.write(filename,AlgLattice(),VacPol[im],size);
 }
 for(int p=0;p<Nmom;p++) sfree(VacPol[p]);
 sfree(VacPol);
 sfree(ext_phase);
}


//------------------------------------------------------------------
// local-conserved, loop over point source to test spectral decomp.
// this only works on scalar machine for now.
//------------------------------------------------------------------
void AlgMuon::VacPolConsLocPtSrc(int t_op, int *xsrc)
{

 char *fname = "VacPolConsLocPtSrc(int t_op)";
 VRB.Func(cname,fname);

 // assumes NO spreadout Ls
 if(GJP.Snodes() != 1){
   ERR.General(cname,fname,"Dummy: Snodes must be 1!\n");
 }

 int i, j;
 int mu, nu, spin, color;
 int n[5], x[5], x_n[5];

 n[0] = GJP.XnodeSites();
 n[1] = GJP.YnodeSites();
 n[2] = GJP.ZnodeSites();
 n[3] = GJP.TnodeSites();
 n[4] = GJP.SnodeSites();
 int size = 4*4; // 4 dirs * 4 dirs

 Float q[4];
 Float L[4];
 L[0] = n[0]*GJP.Xnodes();
 L[1] = n[1]*GJP.Ynodes();
 L[2] = n[2]*GJP.Znodes();
 L[3] = n[3]*GJP.Tnodes();

 const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
 const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
 const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
 const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();

 char filename[100];
 Rcomplex I(0.0,1.0);
 Rcomplex Zero(0.0,0.0);
 int ls = GJP.SnodeSites();

 // external momenta * mu * nu * internal momenta (4-vol) 
 Rcomplex *ext_phase;
 ext_phase = (Rcomplex*)smalloc(Nmom*sizeof(Rcomplex));
 Rcomplex **VacPol = (Rcomplex **)smalloc(Nmom*sizeof(Rcomplex*));
 for(int p=0;p<Nmom;p++) VacPol[p] = (Rcomplex*)smalloc(16*GJP.VolNodeSites()*sizeof(Rcomplex));

 // initialize to zero since we are accumulating sum over s, 
 // and FT at source
 for(mu=0; mu< 4; mu++){
   for(nu=0; nu< 4; nu++){
     for (int site=0; site<GJP.VolNodeSites(); site++) {
       for(int n=0;n<Nmom;n++)VacPol[n][nu+4*(mu+4*site)] = Zero;
     }
   }
 }
#if 0
 int nev=eigcg_arg->nev;
 int m=eigcg_arg->m;
 Float max_eig=eigcg_arg->max_eig_cut;
 int max_def_len=eigcg_arg->max_def_len;
 int restart_len=eigcg_arg->restart_len;
 Float *restart=eigcg_arg->restart;
 bool always_restart=eigcg_arg->always_restart;
 
 int def_len=0;
 int vec_len=GJP.VolNodeSites() * AlgLattice().FsiteSize() / (AlgLattice().FchkbEvl()+1);
 Vector **V=new Vector*[m];
 for(int i=0;i<m;i++)
   V[i]=(Vector *)smalloc(cname,fname,"V",vec_len*sizeof(Float));
 Float *M=new Float[m];
 float **U=new float*[max_def_len];
 for(int i=0;i<max_def_len;i++)U[i]=(float *)smalloc(cname,fname,"U",vec_len*sizeof(float));
 Rcomplex *H=new Rcomplex[max_def_len*max_def_len];
#endif

       //if(zs>3*ptINC) continue;
       // The lepton prop
       //------------------------------------
       qp_arg->cg.mass=alg_muon_arg->loop_mass;
       qp_arg->t=t_op;
       qp_arg->x =  xsrc[0];
       qp_arg->y =  xsrc[1];
       qp_arg->z =  xsrc[2];
       QPropWPointSrc lepton(AlgLattice(), qp_arg, common_arg);
       /*if(nsrc==0){
	 lepton.eig_Run(V,vec_len,M,max_eig,nev,m,U,H,max_def_len,def_len,restart,restart_len,always_restart,0,1e-8);
	 cout<<"deflation length:"<<def_len<<endl;
	 for(int i=0;i<m;i++)sfree(V[i]);
	 delete [] V;
       }else{
	 //solver with restarts 
	 lepton.eig_Run(NULL,vec_len,M,max_eig,0,0,U,H,max_def_len,def_len,restart,restart_len,always_restart,0,1e-8);
       }*/
#if 0
       lepton.eig_Run(V,vec_len,M,max_eig,nev,m,U,H,max_def_len,def_len,restart,restart_len,always_restart,0,1e-8);
#endif
       // source phase (opposite sign to sink (FFT) phase)
       for(int im=0;im<Nmom;im++){
         q[0] = 2.*PI/L[0]*ext_mom[im].cmp(0)*qp_arg->x;
         q[1] = 2.*PI/L[1]*ext_mom[im].cmp(1)*qp_arg->y;
         q[2] = 2.*PI/L[2]*ext_mom[im].cmp(2)*qp_arg->z;
         Float theta = (q[0]+q[1]+q[2]);
         Rcomplex cc(cos(theta),-sin(theta));
         ext_phase[im]=cc;
         //printf("MOM= %d %d %d %e %e %e %e %e\n",mom[0],mom[1],mom[2],q[0],q[1],q[2],phase.real(),phase.imag());
       }
       
       color=0; // Pure QED for now
       // load in prop, one s-slice at a time
       // "spin" and "color" are source indices.
       // 5th dimension loop is cut in half since
       // sum is same for both halves (mult by 2)
       for(x[4] = 0; x[4] < n[4]; x[4]++){
	 int s=x[4];
	 for(mu=0; mu< 4; mu++){//source gamma
	   for(nu=0; nu< 4; nu++){//sink gamma (conserved current)
	     for(x[3] = 0; x[3] < n[3]; x[3]++){
	       for(x[2] = 0; x[2] < n[2]; x[2]++){
		 for(x[1] = 0; x[1] < n[1]; x[1]++){
		   for(x[0] = 0; x[0] < n[0]; x[0]++){
		     
		     // coordinates for neighbor
		     for (i = 0; i < 4; i++ )
		       x_n[i] = (i == nu) ? (x[i]+1)%n[i] : x[i];               // +hop check
		     
		     // offsets
		     int site = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3])));
		     int neighbor_site = x_n[0]+n[0]*(x_n[1]+n[1]*(x_n[2]+n[2]*(x_n[3])));
		     
		     // U_mu(x) where mu = prop_dir
		     Matrix *link = AlgLattice().GaugeField() + 4*site + nu;
		     Matrix tempmat;
		     tempmat.Dagger(*link);                                     // check dag
 		     
		     // mult prop on source by gamma_mu * gamma_5
		     WilsonMatrix mat = lepton(s,site);                         // s check
 		     mat.gr(mu).gr(-5);                                         // gamma check
		     
		     WilsonMatrix wmat_neigh;
		     
		     // get prop in +nu dir.
		     // need the transpose (conjugate later)
		     if( /* offsite */ x[nu]+1 == n[nu] ){
		       
		       getPlusData( (IFloat *)&wmat_neigh,
				    (IFloat *)&lepton(ls-1-s,neighbor_site),    // s,x check
				    288, nu);
		     } else {
		       
		       wmat_neigh = lepton(ls-1-s,neighbor_site);               // s,x check
		     }
		     
		     // mult on sink by gamma_5 (1+gamma_nu) U^dagger_mu
		     WilsonMatrix temp = mat;
		     WilsonMatrix temp2 = temp;
		     temp2 += temp.gl(nu);
		     temp2.gl(-5);                                              // gamma check
		     temp2.LeftTimesEqual(tempmat);                             // link check
		     // spin-color trace
		     WilsonMatrix cc = wmat_neigh;
		     cc.hconj();                                                // conj check
		     for(int im=0;im<Nmom;im++)VacPol[im][nu+4*(mu+4*site)] += 0.5*Trace(cc,temp2)*ext_phase[im];   // tr. check
		     
		     
		     // the other contraction...
		     
		     // mult prop on source by gamma_mu * gamma_5
		     mat = lepton(ls-1-s,site);
		     mat.gr(-5).gr(mu); // check gamma (conj below)
		     
		     // get the fields in plus direction
		     if( /* offsite */ x[nu]+1 == n[nu] ){
		       
		       getPlusData( (IFloat *)&wmat_neigh,
				    (IFloat *)&lepton(s,neighbor_site), 
				    288, nu);
		     } else {
		       
		       wmat_neigh = lepton(s,neighbor_site);               // check s,x
		     }
		     
		     temp = wmat_neigh;
		     temp2 = temp;
		     temp2 -= temp.gl(nu);
		     temp2.gl(-5);
		     temp2.LeftTimesEqual(*link);                     // check gamma, U_mu
		     
		     // spin-color trace
		     cc=mat;
		     cc.hconj();                                      // check s,x, conj
		     for(int im=0;im<Nmom;im++) VacPol[im][nu+4*(mu+4*site)] -=  0.5*Trace(cc, temp2)*ext_phase[im];
#if 0
for(int im=0;im<Nmom;im++) 
printf("VP(mom %d) %d %d   %d %d %d %d %e %e\n",
        im, mu, nu, x[0],x[1],x[2],x[3],
        VacPol[im][nu+4*(mu+4*site)].real(),
        VacPol[im][nu+4*(mu+4*site)].imag());
#endif
		   }
		 }
	       }
	     }
	   }
	 }
       }// s loop
 // free eigcg vectors
#if 0
 delete [] H;
 for(int iu=0;iu<max_def_len;iu++)sfree(U[iu]);
 delete [] U;
 delete [] M;
 for(int i=0;i<m;i++)sfree(V[i]);
 delete [] V;
#endif
      
 // Fourier transform. "size" is number of complex
 // "0" last arg. is for "backward" fft, i.e. exp +ipx
 // which we must do for either the loop or the muon line, so 
 // do it here. i.e. we have to have exp i p(x-y)
 for(int im=0;im<Nmom;im++)FFT4((Float*)VacPol[im], 1, 2*size, 1, 0);
 
      
 // multiply by extra phases from source, point-split sink
 //-------------------------------------------------------
 
 //FileIoType T=ADD_ID;
 //sprintf(filename,"vp-cons-loc-Test.t%d.dat",t_op);
 //FILE *fp=Fopen(T, filename, "w");
 
 for(x[3] = 0; x[3] < n[3]; x[3]++){
   int ptg = x[3]+shiftT;
   for(x[2] = 0; x[2] < n[2]; x[2]++){
     int pzg = x[2]+shiftZ;
     for(x[1] = 0; x[1] < n[1]; x[1]++){
       int pyg = x[1]+shiftY;
       for(x[0] = 0; x[0] < n[0]; x[0]++){
	 int pxg = x[0]+shiftX;
	 
	 
	 int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
	 
	 for(nu=0;nu<4;nu++){
	   for(mu=0;mu<4;mu++){
	     // point-split sink phase
	     q[0] = nu==0 ? PI*pxg/L[0] : 0;
	     q[1] = nu==1 ? PI*pyg/L[1] : 0;
	     q[2] = nu==2 ? PI*pzg/L[2] : 0; 
	     q[3] = nu==3 ? PI*ptg/L[3] : 0;
	     // source phase (opposite sign to sink (FFT) phase)
	     // could be this is actually needed to put muons on-shell!
	     //q[0] -= 2.*PI*pxg/L[0]*qp_arg->x;
	     //q[1] -= 2.*PI*pyg/L[1]*qp_arg->y;
	     //q[2] -= 2.*PI*pzg/L[2]*qp_arg->z;
	     //q[3] -= 2.*PI*ptg/L[3]*qp_arg->t;
	     Float theta = (q[0]+q[1]+q[2]+q[3]);
	     Rcomplex phase(cos(theta),sin(theta));
	     for(int im=0;im<Nmom;im++) {
		VacPol[im][nu+4*(mu+4*site)] *= phase;
	        //if(!UniqueID())printf("VP exmom %d t_op= %d mu= %d nu= %d mom= %d %d %d %d %.15e %.15e\n", 
		     		       //im, t_op, mu, nu, pxg, pyg, pzg, ptg,
		                       //VacPol[im][nu+4*(mu+4*site)].real(), 
		                       //VacPol[im][nu+4*(mu+4*site)].imag());
	     }
	   }
	 }
       }
     }
   }
 }
 //Fclose(T,fp);

 int conf = alg_muon_arg->conf;
 for(int im=0;im<Nmom;im++){
   sprintf(filename,"%s/vacpol-consloc%s%d%d%d.m%g.t%d.mom%d%d%d.dat.%d",
  	   DIRVP,VPTAG,
	   xsrc[0], xsrc[1], xsrc[2],
  	   alg_muon_arg->loop_mass,
  	   t_op,
  	   ext_mom[im].cmp(0),ext_mom[im].cmp(1),ext_mom[im].cmp(2),
  	   alg_muon_arg->conf);
   qio_writeMuon writeQio(argc,argv);
   writeQio.write(filename,AlgLattice(),VacPol[im],size);
 }
 for(int p=0;p<Nmom;p++) sfree(VacPol[p]);
 sfree(VacPol);
 sfree(ext_phase);
}

//------------------------------------------------------------------
// local-conserved, axial current sink, gamma five source
//------------------------------------------------------------------
void AlgMuon::AxialConsLoc(int t_op)
{

 char *fname = "AxialConsLoc(int t_op)";
 VRB.Func(cname,fname);

 // assumes NO spreadout Ls
 if(GJP.Snodes() != 1){
   ERR.General(cname,fname,"Dummy: Snodes must be 1!\n");
 }

 int i, j;
 int mu, nu, spin, color;
 int n[5], x[5], x_n[5];

 n[0] = GJP.XnodeSites();
 n[1] = GJP.YnodeSites();
 n[2] = GJP.ZnodeSites();
 n[3] = GJP.TnodeSites();
 n[4] = GJP.SnodeSites();
 int size = 4*4; // 4 dirs * 4 dirs

 Float q[4];
 Float L[4];
 L[0] = n[0]*GJP.Xnodes();
 L[1] = n[1]*GJP.Ynodes();
 L[2] = n[2]*GJP.Znodes();
 L[3] = n[3]*GJP.Tnodes();

 const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
 const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
 const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
 const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();

 char filename[100];
 Rcomplex I(0.0,1.0);
 Rcomplex Zero(0.0,0.0);
 int ls = GJP.SnodeSites();

 // external momenta * mu * nu * internal momenta (4-vol) 
 Rcomplex *ext_phase;
 ext_phase = (Rcomplex*)smalloc(Nmom*sizeof(Rcomplex));
 Rcomplex **VacPol = (Rcomplex **)smalloc(Nmom*sizeof(Rcomplex*));
 for(int p=0;p<Nmom;p++) VacPol[p] = (Rcomplex*)smalloc(16*GJP.VolNodeSites()*sizeof(Rcomplex));

 // initialize to zero since we are accumulating sum over s, 
 // and FT at source
 for(mu=0; mu< 4; mu++){
   for(nu=0; nu< 4; nu++){
     for (int site=0; site<GJP.VolNodeSites(); site++) {
       for(int n=0;n<Nmom;n++)VacPol[n][nu+4*(mu+4*site)] = Zero;
     }
   }
 }

 // loop over point sources
 int nsrc=0;
 for(int xs=ptStart;xs<n[0]*GJP.Xnodes();xs+=ptINC){
   for(int ys=ptStart;ys<n[1]*GJP.Ynodes();ys+=ptINC){
     for(int zs=ptStart;zs<n[2]*GJP.Znodes();zs+=ptINC){
       //if(zs>3*ptINC) continue;
       // The lepton prop
       //------------------------------------
       qp_arg->cg.mass=alg_muon_arg->loop_mass;
       qp_arg->t=t_op;
       qp_arg->x =  xs;
       qp_arg->y =  ys;
       qp_arg->z =  zs;
       QPropWPointSrc lepton(AlgLattice(), qp_arg, common_arg);

       // source phase (opposite sign to sink (FFT) phase)
       for(int im=0;im<Nmom;im++){
         q[0] = 2.*PI/L[0]*ext_mom[im].cmp(0)*qp_arg->x;
         q[1] = 2.*PI/L[1]*ext_mom[im].cmp(1)*qp_arg->y;
         q[2] = 2.*PI/L[2]*ext_mom[im].cmp(2)*qp_arg->z;
         Float theta = (q[0]+q[1]+q[2]);
         Rcomplex cc(cos(theta),-sin(theta));
         ext_phase[im]=cc;
         printf("MOM(%d)= %e %e %e PHASE= %e %e\n",im,q[0],q[1],q[2], cc.real(),cc.imag());
       }
       
       // load in prop, one s-slice at a time
       for(x[4] = 0; x[4] < n[4]; x[4]++){
	 int s=x[4];
	 Float sign=1.0;
	 if(s < n[4]/2)sign = -1.0;
	 //for(mu=0; mu< 4; mu++){//source gamma mu
	 for(mu=0; mu< 1; mu++){//source gamma 5
	   for(nu=0; nu< 4; nu++){//sink gamma (conserved current)
	     for(x[3] = 0; x[3] < n[3]; x[3]++){
	       for(x[2] = 0; x[2] < n[2]; x[2]++){
		 for(x[1] = 0; x[1] < n[1]; x[1]++){
		   for(x[0] = 0; x[0] < n[0]; x[0]++){
		     
		     // coordinates for neighbor
		     for (i = 0; i < 4; i++ )
		       x_n[i] = (i == nu) ? (x[i]+1)%n[i] : x[i];               // +hop check
		     
		     // offsets
		     int site = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3])));
		     int neighbor_site = x_n[0]+n[0]*(x_n[1]+n[1]*(x_n[2]+n[2]*(x_n[3])));
		     
		     // U_mu(x) where mu = prop_dir
		     Matrix *link = AlgLattice().GaugeField() + 4*site + nu;
		     Matrix tempmat;
		     tempmat.Dagger(*link);                                     // check dag
 		     
		     // mult prop on source by gamma_5 * gamma_5
		     WilsonMatrix mat = lepton(s,site);                         // s check
 		     //mat.gr(mu).gr(-5);                                         // gamma check
		     
		     WilsonMatrix wmat_neigh;
		     
		     // get prop in +nu dir.
		     // need the transpose (conjugate later)
		     if( /* offsite */ x[nu]+1 == n[nu] ){
		       
		       getPlusData( (IFloat *)&wmat_neigh,
				    (IFloat *)&lepton(ls-1-s,neighbor_site),    // s,x check
				    288, nu);
		     } else {
		       
		       wmat_neigh = lepton(ls-1-s,neighbor_site);               // s,x check
		     }
		     
		     // mult on sink by gamma_5 (1+gamma_nu) U^dagger_mu
		     WilsonMatrix temp = mat;
		     WilsonMatrix temp2 = temp;
		     temp2 += temp.gl(nu);
		     temp2.gl(-5);                                              // gamma check
		     temp2.LeftTimesEqual(tempmat);                             // link check
		     // spin-color trace
		     WilsonMatrix cc = wmat_neigh;
		     cc.hconj();                                                // conj check
		     for(int im=0;im<Nmom;im++)VacPol[im][nu+4*(mu+4*site)] += sign*0.5*Trace(cc,temp2)*ext_phase[im];   // tr. check
		     
		     
		     // the other contraction...
		     
		     // mult prop on source by gamma_5 * gamma_5
		     mat = lepton(ls-1-s,site);
		     //mat.gr(-5).gr(mu); // check gamma (conj below)
		     
		     // get the fields in plus direction
		     if( /* offsite */ x[nu]+1 == n[nu] ){
		       
		       getPlusData( (IFloat *)&wmat_neigh,
				    (IFloat *)&lepton(s,neighbor_site), 
				    288, nu);
		     } else {
		       
		       wmat_neigh = lepton(s,neighbor_site);               // check s,x
		     }
		     
		     temp = wmat_neigh;
		     temp2 = temp;
		     temp2 -= temp.gl(nu);
		     temp2.gl(-5);
		     temp2.LeftTimesEqual(*link);                     // check gamma, U_mu
		     
		     // spin-color trace
		     cc=mat;
		     cc.hconj();                                                // conj check
		     for(int im=0;im<Nmom;im++) VacPol[im][nu+4*(mu+4*site)] -=  sign*0.5*Trace(cc, temp2)*ext_phase[im];
#if 1
if(s==n[4]-1)
for(int im=0;im<Nmom;im++) 
printf("AXIAL(mom %d) %d %d   %d %d %d %d %e %e\n",
        im, mu, nu, x[0],x[1],x[2],x[3],
        VacPol[im][nu+4*(mu+4*site)].real(),
        VacPol[im][nu+4*(mu+4*site)].imag());
#endif
		   }
		 }
	       }
	     }
	   }
	 }
       }// s loop
       nsrc++; //count srcs
     }
   }
 } // source loop end      

 for(int p=0;p<Nmom;p++) sfree(VacPol[p]);
 sfree(VacPol);
 sfree(ext_phase);
}


//------------------------------------------------------------------
// local-conserved, Low-mode averaged, connected diagram
//------------------------------------------------------------------
void AlgMuon::VacPolConsLocLowMode(int t_op, int* mom, Rcomplex* VacPol, 
				   Vector* eval, Vector** evec, int Neig)
{

 char *fname = "VacPolConsLocLowMode(int t_op, Rcomplex* VacPol)";
 VRB.Func(cname,fname);

 // assumes NO spreadout Ls
 if(GJP.Snodes() != 1){
   ERR.General(cname,fname,"Dummy: Snodes must be 1!\n");
 }

 int i, j;
 int mu, nu, spin, color;
 int n[5], x[5], x_n[5];
 Rcomplex I(0.0,1.0);
 Rcomplex Zero(0.0,0.0);
 char filename[128];

 n[0] = GJP.XnodeSites();
 n[1] = GJP.YnodeSites();
 n[2] = GJP.ZnodeSites();
 n[3] = GJP.TnodeSites();
 n[4] = GJP.SnodeSites();

 Float q[4];
 Float L[4];
 L[0] = n[0]*GJP.Xnodes();
 L[1] = n[1]*GJP.Ynodes();
 L[2] = n[2]*GJP.Znodes();
 L[3] = n[3]*GJP.Tnodes();
 const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
 const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
 const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
 const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();

 int size = 4*4; // 4 dirs * 4 dirs

 // The lepton prop
 //------------------------------------
 qp_arg->cg.mass=alg_muon_arg->loop_mass;
 qp_arg->t=t_op;
 // these should be set properly in vml file
 //qp_arg->x =  0;
 //qp_arg->y =  0;
 //qp_arg->z =  0;

 int my_node = GJP.TnodeCoor();
 int ts_node = (int)(qp_arg->t/GJP.TnodeSites()); 
 int node_ts = qp_arg->t % GJP.TnodeSites();
 
 Rcomplex **srcBilinear=(Rcomplex **)smalloc(cname,fname,"srcBilinear",(4)*sizeof(Rcomplex*));
 for(int i=0;i<4;i++){
   srcBilinear[i]=(Rcomplex*)smalloc(cname,fname,"srcBilinear[dir]",(Neig*Neig)*sizeof(Rcomplex));
   if(srcBilinear[i]==NULL)
     {ERR.General(cname,fname,"src bilinear %d not initialized\n",i);}
 }

 // construct 4d src bilinear on 1 time slice, ubar gamma_mu u.
 int f_size_5d= GJP.VolNodeSites()*AlgLattice().FsiteSize();
 int f_size_4d= f_size_5d/n[4];
 Vector *u=(Vector *)smalloc(cname,fname,"u",f_size_4d*sizeof(Float));
 if(u==NULL)ERR.General(cname,fname,"4d vector u was not initialized\n");
 Vector *v=(Vector *)smalloc(cname,fname,"v",f_size_4d*sizeof(Float));
 if(v==NULL)ERR.General(cname,fname,"4d vector v was not initialized\n");

 // init to zero
 for(int i=0;i<Neig;i++){
   for(int j=0;j<Neig;j++){
     for(int dir=0;dir<4;dir++) { 
       srcBilinear[dir][i+Neig*j]=Zero;
     }
   }
 }

 // loop over eigenvalues, compute inner product for each i,j pair
 for(int i=0;i<Neig;i++){

   // 4d eigenvector, source, x R, so same as sink
   u->VecZero(f_size_4d);
   AlgLattice().Ffive2four(u, evec[i], n[4]-1, 0);
   
   for(int j=0;j<Neig;j++){
     
     // 4d eignevector, sink, v
     v->VecZero(f_size_4d);
     AlgLattice().Ffive2four(v, evec[j], n[4]-1, 0);
     
     for(x[3] = 0; x[3] < n[3]; x[3]++){
       for(x[2] = 0; x[2] < n[2]; x[2]++){
	 int zg = x[2]+shiftZ;
	 for(x[1] = 0; x[1] < n[1]; x[1]++){
	   int yg = x[1]+shiftY;
	   for(x[0] = 0; x[0] < n[0]; x[0]++){
	     int xg = x[0]+shiftX;
	     
	     int site = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3])));
	     
	     if(my_node != ts_node) continue;
	     if(x[3] != node_ts) continue;
	     
	     // cp Vector at site for each spin into WilsonVector
	     // so we can multiply by gamma matrix
	     int stride = site*24;
	     WilsonVector ttu;
	     WilsonVector ttv;
	     moveFloat((Float*)&ttu, (Float*)u+stride, 24);
	     moveFloat((Float*)&ttv, (Float*)v+stride, 24);
	     // extra g5 for Herm. Dirac Op. We only need one
	     ttv.gamma(-5);
	     
	     // q.x
	     q[0] = 2.*PI*mom[0]/L[0]*xg;
	     q[1] = 2.*PI*mom[1]/L[1]*yg;
	     q[2] = 2.*PI*mom[2]/L[2]*zg; 
	     Float theta = (q[0]+q[1]+q[2]);
	     //theta=0.0; //debug
	     Rcomplex phase(cos(theta),-sin(theta));
	     //printf("MOM= %d %d %d %e %e %e %e %e\n",mom[0],mom[1],mom[2],q[0],q[1],q[2],phase.real(),phase.imag());fflush(stdout);
	     
	     for(int dir=0;dir<4;dir++) {
	       
	       ttv.gamma(dir);
	       // inner product of two WilsonVectors
	       Rcomplex cc;
	       for(int s=0;s<4;s++){
		 for(int c=0;c<3;c++){
		       cc = conj(ttu.d[s].c[c]);
		       srcBilinear[dir][i+Neig*j] += cc * ttv.d[s].c[c] * phase;
		 }
	       }
	       ttv.gamma(dir);
	     }
	   }
	 }
       }
     }
     //end loops over x
   }
 }

 // Global sum and print
 for(int dir=0;dir<4;dir++) {
   for(int i=0;i<Neig;i++){
     for(int j=0;j<Neig;j++){
       
       Float re = srcBilinear[dir][i+Neig*j].real();
       Float im = srcBilinear[dir][i+Neig*j].imag();
       slice_sum(&re,1,99);
       slice_sum(&im,1,99);
       Rcomplex cc(re,im);
       srcBilinear[dir][i+Neig*j] = cc;
       if(!UniqueID())printf("Source Bilinear dir= %d i= %d j= %d %e %e\n",
			     dir, i,j, 
			     srcBilinear[dir][i+Neig*j].real(), 
			     srcBilinear[dir][i+Neig*j].imag());
       fflush(stdout);
     }
   }
 }

 sfree(v);
 sfree(u);

 // calc vacpol at each site.
 // Construct sink bi-linear and mult by source bi-linear
 for(x[3] = 0; x[3] < n[3]; x[3]++){
   for(x[2] = 0; x[2] < n[2]; x[2]++){
     for(x[1] = 0; x[1] < n[1]; x[1]++){
       for(x[0] = 0; x[0] < n[0]; x[0]++){

	 int site = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3])));
	 
	 for(mu=0; mu< 4; mu++){//source gamma
	   for(nu=0; nu< 4; nu++){//sink gamma (conserved current)

	     VacPol[nu+4*(mu+4*site)] = Zero;
	     Matrix *link = AlgLattice().GaugeField() + 4*site + nu;
	     Matrix cclink;
	     cclink.Dagger(*link);

	     //sum over s, eig value
	     for(x[4] = 0; x[4] < n[4]; x[4]++){

	       // s
	       int site5d = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3]+n[3]*x[4])));
	       // ls-1-s
	       int site5dR = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3]+n[3]*(n[4]-1-x[4]))));

	       // only every 3rd one, others are degenerate for pure QED
	       for(int i=0;i<Neig;i++){

		 Float eval_i = *((Float*)eval + i);

		 for(int j=0;j<Neig;j++){
		   
		   Float eval_j = *((Float*)eval + j);
		   
		   // coordinates for neighbor
		   for (int ii = 0; ii < 4; ii++ )
		     x_n[ii] = (ii == nu) ? (x[ii]+1)%n[ii] : x[ii];
		   int neighbor_site5d = x_n[0]+n[0]*(x_n[1]+n[1]*(x_n[2]+n[2]*(x_n[3]+n[3]*x[4])));
		   int neighbor_site5dR = x_n[0]+n[0]*(x_n[1]+n[1]*(x_n[2]+n[2]*(x_n[3]+n[3]*(n[4]-1-x[4]))));
		   
		   // cp Vector at site for each spin into WilsonVector
		   // so we can multiply by gamma matrix
		   WilsonVector psi; // psi(x,s)
		   int stride = site5d*24;
		   moveFloat((Float*)&psi,(Float*)(evec[j])+stride, 24);
		   WilsonVector psi_p; // psi(x+nu,s)
		   int neighbor_stride = neighbor_site5dR*24;
		   if( /*offsite*/  x[nu]+1 == n[nu] ){
		     getPlusData( (Float *)&psi_p, (Float *)(evec[i])+neighbor_stride, 24, nu);
		   } else {
		     moveFloat((Float*)&psi_p,(Float*)(evec[i])+neighbor_stride, 24);
		   }
		   
		   // compute psibar(x+nu,s) U^dagger(x)_nu g5 (1+gamma_nu) psi(x,s):
		   
		   // extra g5 for Herm Dirac Op.
		   psi.gamma(-5);
		   WilsonVector temp = psi;
		   temp.gamma(nu);
		   psi -= temp;
		   psi.LeftTimesEqual(cclink);
		   Rcomplex cc=Zero;
		   // dot product of two WilsonVectors
		   for(int s=0;s<4;s++){
		     for(int c=0;c<3;c++){
		       cc +=  conj(psi_p.d[s].c[c])*psi.d[s].c[c];
		     }
		   }
		   
		   VacPol[nu+4*(mu+4*site)] += 0.5 * cc / eval_i / eval_j * srcBilinear[mu][j+Neig*i];
		   //VacPol[nu+4*(mu+4*site)] += 0.5*cc / eig_val;
		   
		   // compute -psibar(x,s) U^dagger(x)_nu g5 (1-gamma_nu) psi(x+nu,s):
		   
		   // cp Vector at site for each spin into WilsonVector
		   // so we can multiply by gam\ma matrix
		   // psi(x,s)
		   stride = site5dR*24;
		   moveFloat((Float*)&psi, (Float*)(evec[i])+stride, 24);
		   
		   // psi(x+nu,s)
		   neighbor_stride = neighbor_site5d*24;
		   if(  /*offsite*/  x[nu]+1 == n[nu] ){
		     getPlusData( (IFloat *)&psi_p, (IFloat *)(evec[j])+neighbor_stride, 24, nu);
		   } else {
		     moveFloat((Float*)&psi_p,(Float*)(evec[j])+neighbor_stride, 24);
		   }
		   
		   psi_p.gamma(-5);
		   temp = psi_p;
		   temp.gamma(nu);
		   psi_p += temp;
		   psi_p.LeftTimesEqual(*link);
		   cc=Zero;
		   // dot product of two WilsonVectors
		   for(int s=0;s<4;s++){
		     for(int c=0;c<3;c++){
		       cc +=  conj(psi.d[s].c[c])*psi_p.d[s].c[c];
		     }
		   }
		   
		   VacPol[nu+4*(mu+4*site)] -=  0.5 * cc / eval_i / eval_j * srcBilinear[mu][j+Neig*i];
		 }
	       }
	       //end loop on eigenvalues
	     }
	   }
	 }
       }
     }
   }
 }
 for(int i=0;i<4;i++){
     sfree(srcBilinear[i]);
 }
 sfree(srcBilinear);

 // Fourier transform. "size" is number of complex
 // "0" last arg. is for "backward" fft, i.e. exp +ipx
 // which we must do for either the loop or the muon line, so 
 // do it here. i.e. we have to have exp i p(x-y)
 FFT4((Float*)VacPol, 1, 2*size, 1, 0);

 // multiply by extra phases from source, point-split sink
 //-------------------------------------------------------
 
 FileIoType T=ADD_ID;
 sprintf(filename,"vp-cons-loc-LowMode.mom%d%d%d.t%d.dat",mom[0],mom[1],mom[2],t_op);
 FILE *fp=Fopen(T, filename, "w");
 
 for(x[3] = 0; x[3] < n[3]; x[3]++){
   int ptg = x[3]+shiftT;
   for(x[2] = 0; x[2] < n[2]; x[2]++){
     int pzg = x[2]+shiftZ;
     for(x[1] = 0; x[1] < n[1]; x[1]++){
       int pyg = x[1]+shiftY;
       for(x[0] = 0; x[0] < n[0]; x[0]++){
	 int pxg = x[0]+shiftX;
	 

	 int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
	 
	 for(nu=0;nu<4;nu++){
	   for(mu=0;mu<4;mu++){
	     // point-split sink phase
	     q[0] = nu==0 ? PI*pxg/L[0] : 0;
	     q[1] = nu==1 ? PI*pyg/L[1] : 0;
	     q[2] = nu==2 ? PI*pzg/L[2] : 0; 
	     q[3] = nu==3 ? PI*ptg/L[3] : 0;
	     // source phase (opposite sign to sink (FFT) phase)
	     // could be this is actually needed to put muons on-shell!
	     //q[0] -= 2.*PI*pxg/L[0]*qp_arg->x;
	     //q[1] -= 2.*PI*pyg/L[1]*qp_arg->y;
	     //q[2] -= 2.*PI*pzg/L[2]*qp_arg->z;
	     //q[3] -= 2.*PI*ptg/L[3]*qp_arg->t;
	     Float theta = (q[0]+q[1]+q[2]+q[3]);
	     Rcomplex phase(cos(theta),sin(theta));
	     VacPol[nu+4*(mu+4*site)] *= phase;
	     //Fprintf(fp,"t_op= %d mu= %d nu= %d mom= %d %d %d %d %e %e\n", 
		     //t_op, mu, nu, pxg, pyg, pzg, ptg,
		     //VacPol[nu+4*(mu+4*site)].real(), 
		     //VacPol[nu+4*(mu+4*site)].imag());
             
	   }
	 }
       }
     }
   }
 }
 Fclose(T,fp);

 int conf = alg_muon_arg->conf;
 sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.mom%d%d%d.%s.dat.%d",
	 DIRVP,VPTAG,
	 alg_muon_arg->loop_mass,
	 t_op,
	 mom[0],mom[1],mom[2],
	 EIGTAG,
	 alg_muon_arg->conf);
 qio_writeMuon writeQio(argc,argv);
 writeQio.write(filename,AlgLattice(),VacPol,size);

}


// calc the muon line with an insertion of the
// conserved current
void AlgMuon::MuonLineCons(int *out_mom, QPropWMomSrc &in, Rcomplex* Muon)
{
  
  char *fname = "MuonLineCons(int *out_mom, QPropWMomSrc &min, Rcomplex* Muon))";
  VRB.Func(cname,fname);
  int x[5], x_n[4];
  int n[5];
  int nu, mu, proj;
  int i, spin, color;
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  Rcomplex Ihalf(0.0,0.5);
  Rcomplex Zero(0.0,0.0);
  WilsonMatrix temp;  
  int size = 4*NPROJ; // 4 dirs * 4 projectors
  int fsize = 3*4;
  char filename[100];
  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();

  // outgoing muon with momentum p
  //------------------------------
  qp_arg->cg.mass = alg_muon_arg->line_mass ;
  qp_arg->t = alg_muon_arg->source_time;
  // U1 links only
  LatticeContainer lat_cont;
  // save lattice
  lat_cont.Get(AlgLattice());
  AlgLattice().SetGfieldOrd();
  AlgLattice().mult_su3_links_by_u1(alg_muon_arg->charge);
  QPropWMomSrc out(AlgLattice(), qp_arg, out_mom, common_arg);
  out.Run();
  // restore lattice
  lat_cont.Set(AlgLattice());

  // initialize to zero since we are accumulating sum over s
  for(mu=0; mu< 4; mu++){
    for(proj=0; proj< NPROJ; proj++){
      for (int site=0; site<GJP.VolNodeSites(); site++) {
	Muon[mu+4*(proj+NPROJ*site)] = Zero;
      }
    }
  }

  //FileIoType T=ADD_ID;
  //sprintf(filename,"ml-consloc-cs-sumS.mom%d%d%d.dat",out_mom[0],out_mom[1],out_mom[2]);
  //FILE *fp=Fopen(T, filename, "w");
  
  // 5th dimension loop is cut in half since
  // sum is same for both halves (mult by 2)
  // check first
  int ls = n[4];
  for(x[4] = 0; x[4] < n[4]; x[4]++){ 
    int s=x[4];
    for(nu=0; nu< 4; nu++){//sink gamma (conserved current)
      for(x[3] = 0; x[3] < n[3]; x[3]++){
	int tg = x[3]+shiftT;
	for(x[2] = 0; x[2] < n[2]; x[2]++){
	  int zg = x[2]+shiftZ;
	  for(x[1] = 0; x[1] < n[1]; x[1]++){
	    int yg = x[1]+shiftY;
	    for(x[0] = 0; x[0] < n[0]; x[0]++){
	      int xg = x[0]+shiftX;

	      
	      // coordinates for neighbor
	      for (i = 0; i < 4; i++ )
		x_n[i] = (i == nu) ? (x[i]+1)%n[i] : x[i];               // +hop check
	      
	      // offsets
	      int site = x[0]+n[0]*(x[1]+n[1]*(x[2]+n[2]*(x[3])));
	      int neighbor_site = x_n[0]+n[0]*(x_n[1]+n[1]*(x_n[2]+n[2]*(x_n[3])));
	      
	      // U_nu(x) where nu = prop_dir
	      Matrix *link = AlgLattice().GaugeField() + 4*site + nu;
	      Matrix tempmat;
	      tempmat.Dagger(*link);                                     // check dag
	      
	      // start with incoming
	      WilsonMatrix mat = in(s,site);                           // s check
	      //can't, proj: mat.gr(-5);                              // gamma check
	      
	      WilsonMatrix wmat_neigh;
	      
	      // get prop in +nu dir.
	      // need the transpose (conjugate later)
	      if( /* offsite */ x[nu]+1 == n[nu] ){
		
		getPlusData( (IFloat *)&wmat_neigh,
			     (IFloat *)&out(ls-1-s,neighbor_site),     // s,x check
			     288, nu);
	      } else {
		
		wmat_neigh = out(ls-1-s,neighbor_site);                      // s,x check
	      }
	      
	      // mult on sink by gamma_5 (1+gamma_nu) U^dagger_nu
	      WilsonMatrix temp = mat;
	      WilsonMatrix temp2 = temp;
	      temp2 += temp.gl(nu);
	      temp2.gl(-5);                                              // gamma check
	      temp2.LeftTimesEqual(tempmat);                             // link check
	      // spin-color trace
	      WilsonMatrix cc = wmat_neigh;
	      cc.hconj();                                                // conj check
	      cc.gl(-5);                                                // conj check
	      WilsonMatrix line = cc*temp2;                              // tr. check
	      WilsonMatrix line1 = cc*temp2;                              // debug

	      // the other contraction...
	      
	      // mult sink by gamma5
	      mat = out(ls-1-s,site);
	      mat.gr(-5); // (conj below)
	      
	      // get the fields in plus direction
	      if( /* offsite */ x[nu]+1 == n[nu] ){
		
		getPlusData( (IFloat *)&wmat_neigh,
			     (IFloat *)&in(s,neighbor_site), 
			     288, nu);
	      } else {
		
		wmat_neigh = in(s,neighbor_site);               // check s,x
	      }
	      
	      temp = wmat_neigh;
	      temp2 = temp;
	      temp2 -= temp.gl(nu);
	      temp2.gl(-5);
	      temp2.LeftTimesEqual(*link);                     // check gamma, U_nu
	      
	      // spin-color trace
	      cc=mat;
	      cc.hconj();                                      // check s,x, conj
	      line -=  cc * temp2;
	      WilsonMatrix line2 = cc*temp2; //debug
	      
	      // project with (1+gamma_t)/2 (1, gamma_x gamma_y, gamma_x gamma_z, and perms.)
	      //                 proj=       0, 1              , 2
	      line.PParProjectSink(); //(multiplies on the left)
	      // project some more, and trace
	      Muon[nu+4*(0+NPROJ*site)] +=   0.5 * line.Trace(); 
	      line.gl(0).gl(1); // i gy gx * line
	      Muon[nu+4*(1+NPROJ*site)] += Ihalf * line.Trace(); 
	      line.gl(1).gl(2); // i gz gx ...
	      Muon[nu+4*(2+NPROJ*site)] += Ihalf * line.Trace();
	      line.gl(0).gl(1); // i gz gy ...
	      Muon[nu+4*(3+NPROJ*site)] += Ihalf * line.Trace(); 
	    }
	  }
	}
      }
    } // nu
  } // s
  //Fclose(T,fp);

  // Fourier transform the muon line (vertex coordinate->vertex momentum)
  FFT4((Float*)Muon, 1, 2*size, 1, 1);

  // multiply by extra phases from source, point-split sink
  //-------------------------------------------------------

  int conf = alg_muon_arg->conf;
 
  FileIoType T=ADD_ID;
  /*sprintf(filename,"%s/muonline-consloc.m%g.i000.o%d%d%d.dat.%d",
	  DIRML,
	  alg_muon_arg->line_mass,
	  out_mom[0],out_mom[1],out_mom[2],
	  alg_muon_arg->conf);*/
  //FILE *fp=Fopen(T, filename, "w");
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
	  
	  
	  int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
	  
	  for(nu=0;nu<4;nu++){
	    // point-split sink phase (these should be opp. to vacpol ones)
	    q[0] = nu==0 ? -PI*pxg/L[0] : 0;
	    q[1] = nu==1 ? -PI*pyg/L[1] : 0;
	    q[2] = nu==2 ? -PI*pzg/L[2] : 0;
	    q[3] = nu==3 ? -PI*ptg/L[3] : 0;
	    // this is the source phase from both props. Since they were already
	    // mom src props, only need the t component. And the incoming is 
	    // hard-wired at t=0!
	    //q[3] += 2.*PI*ptg/L[3]*qp_arg->t;
	    Float theta = (q[0]+q[1]+q[2]+q[3]);
	    Rcomplex phase(cos(theta),sin(theta));
	    for(proj=0;proj<NPROJ;proj++){
	      Muon[nu+4*(proj+NPROJ*site)] *= phase;
	      //debug:
	      /*Rcomplex cc(1.,0.);
	      int fact = (ptg*pxg*pyg*pzg);
	      if(fact==0){fact=1;}
	      Muon[nu+4*(proj+NPROJ*site)] = cc*fact;*/
	      //Fprintf(fp,"proj= %d nu= %d mom= %d %d %d %d %e %e\n", 
	      //      proj, nu, pxg, pyg, pzg, ptg,
	      //Fprintf(fp,"%le %le\n", 
	      //      Muon[nu+4*(proj+NPROJ*site)].real(), 
	      //      Muon[nu+4*(proj+NPROJ*site)].imag());
	    }
	  }
	}
      }
    }
  }
  //Fclose(T,fp);

  qio_writeMuon writeQio(argc,argv);
  sprintf(filename,"%s/muonline-consloc.m%g.i000.o%d%d%d.dat.%d",
	  DIRML,
	  alg_muon_arg->line_mass,
	  out_mom[0],out_mom[1],out_mom[2],
	  alg_muon_arg->conf);
  writeQio.write(filename,AlgLattice(),Muon,size); fflush(stdout); 
  
  // fourier transform the sink of out-going prop.
  two_point(out_mom,out);

}




double AlgMuon::photon_prop(int *q,int mu,double momX,double momT)
{
  double qsq;
  qsq = q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
  qsq *= momX*momX;
  qsq += q[3]*q[3]*momT*momT;
  return 1./qsq;
}



void AlgMuon::run()
{

  char *fname = "run()";
  VRB.Func(cname,fname);
  int mu;
  int i,j,k;
  int size = 4*GJP.VolNodeSites();

  // set up for fourier transforms
  //----------------------------------------------------------------
  int in_mom[3];
  in_mom[0]=in_mom[1]=in_mom[2]=0;
  int out_mom[3];
  // fourier transform the external vertex 
  // up to and including this momentum

  if(DO_MUONLINE){
    // holds muon line (4 * 4: gamma_mu times 4 projectors)
    Rcomplex *Muon;
    Muon = (Rcomplex*)smalloc(NPROJ*size*sizeof(Rcomplex));
    {
      // incoming muon prop, only need it once
      qp_arg->cg.mass = alg_muon_arg->line_mass ;
      qp_arg->t =  0;
      int saveSrcColorStart=qp_arg->StartSrcColor;
      int saveSrcColorEnd=qp_arg->EndSrcColor;
      qp_arg->StartSrcColor=0;
      qp_arg->EndSrcColor=1;
      QPropWMomSrc min(AlgLattice(), qp_arg, in_mom, common_arg);
      min.Run();
      // fft and print trace of in-coming prop
      two_point(in_mom, min); // 1=conj the prop
      // Loop over external photon momentum
      for(i=0; i < Nmom; i++){
	// Muonline
	out_mom[0]=ext_mom[i].cmp(0);
	out_mom[1]=ext_mom[i].cmp(1);
	out_mom[2]=ext_mom[i].cmp(2);
	MuonLineCons(out_mom, min, Muon);
      }
      qp_arg->StartSrcColor=saveSrcColorStart;
      qp_arg->EndSrcColor=saveSrcColorEnd;
    }
    sfree(Muon);
  }


  if(DO_VACPOL){

    int f_size_5d= GJP.VolNodeSites()*AlgLattice().FsiteSize();
    int f_size_4d= f_size_5d/GJP.SnodeSites();

#if 0    
    EigArg eig_arg;
    eig_arg.Decode("eig_arg.vml","eig_arg");
    char filename[128];
    // set mass
    eig_arg.mass=alg_muon_arg->loop_mass;
    eig_arg.Mass_init=alg_muon_arg->loop_mass;
    sprintf(filename,"%s/eig_arg.dat",DIRVP);
    eig_arg.Encode(filename,"eig_arg");
    
    int Neig=eig_arg.N_eig;
    Vector *eval=(Vector *)smalloc(cname,fname,"eigval",Neig*sizeof(Float));
    Vector **evec=(Vector **)smalloc(cname,fname,"evec",(Neig)*sizeof(Vector*));
    for(i=0;i<Neig;i++){
      evec[i]=(Vector*)smalloc(cname,fname,"evec[i]",(f_size_5d)*sizeof(Float));
      if(evec[i]==NULL)
	{ERR.General(cname,fname,"vec %d not initialized\n",i);}
      evec[i]->VecZero(f_size_5d);
    }
    
    // Get the Eigenvectors and values!
    AlgEig eig( AlgLattice(), common_arg, &eig_arg );
    if(Neig>0) eig.run(evec,eval);
#endif

    for(int t=alg_muon_arg->oper_time_start; t <= alg_muon_arg->oper_time_end; t++){

      if(alg_muon_arg->vp_kind==CONSERVED_LOCAL){
        Rcomplex *vacpol;
        vacpol = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));
	VacPolConsLoc(t,vacpol);
	sfree(vacpol);
      }else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_TWISTED){
        Rcomplex *vacpol;
        vacpol = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));
	VacPolConsLocTwistedBC(t,vacpol);
	sfree(vacpol);
      }else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOOP_PTSRC){
	VacPolConsLocLoopPtSrc(t);
	//AxialConsLoc(t);
      }else if( alg_muon_arg->vp_kind==CONSERVED_LOCAL_PTSRC_TEST){
        Rcomplex *vacpol;
        vacpol = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));
	// Loop over external photon momentum
	for(i=0; i <= XMOM; i++){
	  for(j=0; j <= YMOM; j++){
	    for(k=0; k <= ZMOM; k++){
	      
	      //skip zero momentum, and others...
	      if(i==0 && j==0 && k==0)continue;
	      int momsq = i*i+j*j+k*k;
	      if(momsq>1)continue;
	      
	      out_mom[0]=i;out_mom[1]=j;out_mom[2]=k;	      
	      VacPolConsLocTest(t,vacpol,out_mom);
	    }
	  }
	}
	sfree(vacpol);
#if 0
      }else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOWMODE && Neig>0 ){

        Rcomplex *vacpol;
        vacpol = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));
	
	for(int n=Neig;n<=Neig;n+=Neig/3){ // use different number of eigvec's in loop
	  
	  // Loop over external photon momentum
	  for(i=0; i <= XMOM; i++){
	    for(j=0; j <= YMOM; j++){
	      for(k=0; k <= ZMOM; k++){
		
		//skip zero momentum, and others...
		if(i==0 && j==0 && k==0)continue;
		int momsq = i*i+j*j+k*k;
		if(momsq>1)continue;
		
		out_mom[0]=i;out_mom[1]=j;out_mom[2]=k;		
		VacPolConsLocLowMode(t,out_mom,vacpol,eval,evec,n);
	      }
	    }
	  }
	} // eigvec's
        sfree(vacpol);
#endif
      }else{
	//if(!UniqueID())ERR.General(cname,fname, "Unrecognized Vac Pol type (%d) or Neig<=0 (%d).\n", alg_muon_arg->vp_kind, Neig);	
	if(!UniqueID())ERR.General(cname,fname, "Unrecognized Vac Pol type (%d)\n", alg_muon_arg->vp_kind);
      }
    }// tOp
#if 0
    sfree(eval);
    for(int i=0;i<Neig;i++){
      sfree(evec[i]);
    }
    sfree(evec);
#endif
  }//vp
  //check LBL qio read/write
  //out_mom[0]=1;out_mom[1]=0;out_mom[2]=0;
  //run_lbl(out_mom, alg_muon_arg->oper_time_start, Muon, vacpol);
}



// check lbl "on the fly", ie don't read from disk
void AlgMuon::run_lbl(int *out_mom, 
		      int t_op, 
		      Rcomplex *Muon, 
		      Rcomplex *vacpol)
{

  char *fname = "run_lbl(i*,i,*Rc,*Rc)";
  VRB.Func(cname,fname);
  int mu,nu,rho;
  int config = alg_muon_arg->conf;
  //----------------------------------------------------------------
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();

  Rcomplex Zero(0.0,0.0);
  int t;
  int i,j,k;
  int x[4];
  // fourier transform the external vertex 
  // up to and including this momentum
  Float momX = 2.*PI/(GJP.Xnodes()*GJP.XnodeSites());
  Float momT = 2.*PI/(GJP.Tnodes()*GJP.TnodeSites());

  // projector
  for(int proj=0;proj<NPROJ;proj++){
    // ext. vertex lorentz index
    for(rho=0;rho<4;rho++){
      //sum over q^2:
      Rcomplex sum=Zero;
      for(x[3] = 0; x[3] < n[3]; x[3]++){
	int tg = x[3]+shiftT;
	Float pt = sin(0.5*tg*momT);
	for(x[2] = 0; x[2] < n[2]; x[2]++){
	  int zg = x[2]+shiftZ;
	  Float pz = sin(0.5*zg*momX);
	  for(x[1] = 0; x[1] < n[1]; x[1]++){
	    int yg = x[1]+shiftY;
	    Float py = sin(0.5*yg*momX);
	    for(x[0] = 0; x[0] < n[0]; x[0]++){
	      int xg = x[0]+shiftX;
	      Float px= sin(0.5*xg*momX);
	      
	      Float qsq = 4.0*(px*px+py*py+pz*pz+pt*pt);
	      if(qsq<EPS)continue;//avoid zero momentum
	      int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
	      qsq = 1.0/qsq;
	      // internal Lorentz indices
	      for(nu=0;nu<4;nu++){
		for(mu=0;mu<4;mu++){
		  //for now, take mu=nu
		  if(mu!=nu)continue;
		  // new convention for cons-loc vacpol computed after 6/26/08!
		  sum += qsq*Muon[mu+4*(proj+NPROJ*site)]*vacpol[nu+4*(rho+4*site)];
		}
	      }
	    }
	  }
	}
      }
      // global sum
      slice_sum((Float*)&sum,2,99);
      if(!UniqueID())printf("CHECK LBL conf= %d t_op= %d proj= %d rho= %d mom= %d %d %d  %le %le \n",
			    config,t_op,proj,rho,
			    out_mom[0],out_mom[1],out_mom[2],
			    sum.real(),sum.imag());
    }
  }
}


bool FileExists(string strFilename) { 
  struct stat stFileInfo; 
  bool blnReturn; 
  int intStat; 

  // Attempt to get the file attributes 
  intStat = stat(strFilename.c_str(),&stFileInfo); 
  if(intStat == 0) { 
    // We were able to get the file attributes 
    // so the file obviously exists. 
    blnReturn = true; 
  } else { 
    // We were not able to get the file attributes. 
    // This may mean that we don't have permission to 
    // access the folder which contains this file. If you 
    // need to do that level of checking, lookup the 
    // return values of stat which will give you 
    // more details on why stat failed. 
    blnReturn = false; 
  } 
   
  return(blnReturn); 
}

void AlgMuon::run_lbl()
{

  char *fname = "run_lbl()";
  VRB.Func(cname,fname);
  int mu,nu,rho;
  // total number of configurations
  int nconfs = NC;
  int start_config = alg_muon_arg->conf;
  int size = 4*GJP.VolNodeSites();
  // holds muon line
  Rcomplex *Muon;
  Muon = (Rcomplex*)smalloc(NPROJ*size*sizeof(Rcomplex));
  // holds Vacuum Polarization
  Rcomplex *vacpol;
  vacpol = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));

   //----------------------------------------------------------------
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();

  Rcomplex Zero(0.0,0.0);
  int t;
  int i,j,k;
  int x[4];
  // fourier transform the external vertex 
  // up to and including this momentum
  Float momX = 2.*PI/(GJP.Xnodes()*GJP.XnodeSites());
  Float momT = 2.*PI/(GJP.Tnodes()*GJP.TnodeSites());

  if(alg_muon_arg->vp_kind!=CONSERVED_LOCAL &&
     alg_muon_arg->vp_kind!=CONSERVED_LOCAL_LOWMODE &&
     alg_muon_arg->vp_kind!=CONSERVED_LOCAL_PTSRC_TEST &&
     alg_muon_arg->vp_kind!=CONSERVED_LOCAL_LOOP_PTSRC){
    ERR.General(cname,fname, "VacPol type %d invalid.\n", alg_muon_arg->vp_kind);
  }

  // configs
  for(int nc=start_config;nc<nconfs;nc++){
    
    for(t=alg_muon_arg->oper_time_start; t <= alg_muon_arg->oper_time_end; t+=tINC){

      // read in data and
      // average over equivalent momenta
      char filename[200];
      // vacpol: use this for point source
      if(alg_muon_arg->vp_kind==CONSERVED_LOCAL){
	qio_readMuon readVac(argc,argv);
        sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.xyz%d%d%d.%d",
         	DIRVP,VPTAG,
	        alg_muon_arg->loop_mass,
		t,
	        qp_arg->x,
         	qp_arg->y,
         	qp_arg->z,
		nc);
	if(FileExists(filename)){
	  readVac.read(filename,AlgLattice(),vacpol,16);
	}else{
          if(!UniqueID())printf("%s does not exist, skipping\n",filename);
	  continue;
	}
      }
      for(i = 0; i < Nmom; i++){
	    
 	    // vacpol (use this code if loop source is fft'd,
	    // i.e., for low modes or random src's.
	    if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOWMODE){
	      qio_readMuon readVac(argc,argv);
	      sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.mom%d%d%d.%s.dat.%d",
		      DIRVP,VPTAG,
		      alg_muon_arg->loop_mass,
		      t,
		      ext_mom[i].cmp(0),
		      ext_mom[i].cmp(1),
		      ext_mom[i].cmp(2),
		      EIGTAG,
		      nc);
	      if(!UniqueID())printf("trying %s\n",filename);
	      if(FileExists(filename)){
		if(!UniqueID())printf("qio trying %s\n",filename);
		readVac.read(filename,AlgLattice(),vacpol,16);
	      }else{
                if(!UniqueID())printf("%s does not exist, skipping\n",filename);
		continue;
	      }
	    } else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_PTSRC_TEST || alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOOP_PTSRC){
	      qio_readMuon readVac(argc,argv);
	      sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.mom%d%d%d.dat.%d",
		      DIRVP,VPTAG,
		      alg_muon_arg->loop_mass,
		      t,
		      ext_mom[i].cmp(0),
		      ext_mom[i].cmp(1),
		      ext_mom[i].cmp(2),
		      nc);
	      if(!UniqueID())printf("trying %s\n",filename);
	      if(FileExists(filename)){
		if(!UniqueID())printf("qio trying %s\n",filename);
		readVac.read(filename,AlgLattice(),vacpol,16);
	      }else{
                if(!UniqueID())printf("%s does not exist, skipping\n",filename);
		continue;
	      }
	    }else{
    		ERR.General(cname,fname, "VacPol type %d invalid.\n", alg_muon_arg->vp_kind);
  	    }

	    
	    // get the Muonline
	    qio_readMuon readMu(argc,argv);
	    sprintf(filename,"%s/muonline-consloc.m%g.i000.o%d%d%d.dat.%d",
		    DIRML,
		    alg_muon_arg->line_mass,
		    ext_mom[i].cmp(0),
		    ext_mom[i].cmp(1),
		    ext_mom[i].cmp(2),
		    nc);
	    if(FileExists(filename)){
	      readMu.read(filename,AlgLattice(),Muon,4*NPROJ);
	    }else{
              if(!UniqueID())printf("%s does not exist, skipping\n",filename);
	      continue;
	    }

	    // projector
	    for(int proj=0;proj<NPROJ;proj++){
	      // ext. vertex lorentz index
	      for(rho=0;rho<4;rho++){
		
		//sum over q^2:
		Rcomplex sum=Zero;
		for(x[3] = 0; x[3] < n[3]; x[3]++){
		  int tg = x[3]+shiftT;
		  //Float pt = cos(tg*momT);
		  Float pt = sin(0.5*tg*momT);
		  for(x[2] = 0; x[2] < n[2]; x[2]++){
		    int zg = x[2]+shiftZ;
		    //Float pz = cos(zg*momX);
		    Float pz = sin(0.5*zg*momX);
		    for(x[1] = 0; x[1] < n[1]; x[1]++){
		      int yg = x[1]+shiftY;
		      //Float py = cos(yg*momX);
		      Float py = sin(0.5*yg*momX);
		      for(x[0] = 0; x[0] < n[0]; x[0]++){
			int xg = x[0]+shiftX;
			//Float px = cos(xg*momX);
			Float px= sin(0.5*xg*momX);

			//Float qsq = 8.0-2.0*(px+py+pz+pt);
			Float qsq = 4.0*(px*px+py*py+pz*pz+pt*pt);
			if(qsq<EPS)continue;//avoid zero momentum
			int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
			qsq = 1.0/qsq;
			// internal Lorentz indices
			for(nu=0;nu<4;nu++){
			  for(mu=0;mu<4;mu++){
			    //for now, take mu=nu
			    if(mu!=nu)continue;
			    // new convention for cons-loc vacpol computed after 6/26/08!
			    sum += qsq*Muon[mu+4*(proj+NPROJ*site)]*vacpol[nu+4*(rho+4*site)];
			  }
			}
		      }
		    }
		  }
		}
		// global sum
		slice_sum((Float*)&sum,2,99);
		if(!UniqueID())printf("LBL conf= %d t_op= %d proj= %d rho= %d mom= %d %d %d  %.15e %.15e \n",
				      nc,t,proj,rho,ext_mom[i].cmp(0),ext_mom[i].cmp(1),ext_mom[i].cmp(2),sum.real(),sum.imag());
	      }
	    }
      }//mom
    }//t
  }//confs
  sfree(vacpol);
  sfree(Muon);
  //close(fp);
}




// calc the subtraction term
void AlgMuon::run_sub()
{

  char *fname = "run_sub()";
  VRB.Func(cname,fname);
  int mu,rho;
  // total number of configurations
  int nconfs = 30; //debug =1, fix this!
  int size = 4*GJP.VolNodeSites();
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  // holds muon line
  Rcomplex *Muon;
  Muon = (Rcomplex*)smalloc(NPROJ*size*sizeof(Rcomplex));
  // holds Vacuum Polarization
  Rcomplex *vacpol;
  vacpol = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));
  // holds data for global sum
  int TSize = n[3]*GJP.Tnodes();
  // holds muon line
  Rcomplex *temp;
  if(4*NPROJ > 16)ERR.General(cname,fname, "size of temp is not enough!");
  temp = (Rcomplex*)smalloc(16*TSize*sizeof(Rcomplex));

  //----------------------------------------------------------------

  // this is my node coordinate
  int myXcoor = GJP.XnodeCoor();
  int myYcoor = GJP.YnodeCoor();
  int myZcoor = GJP.ZnodeCoor();
  int myTcoor = GJP.TnodeCoor();
  
  Rcomplex Zero(0.0,0.0);
  int t;
  int i,j,k,l;
  // fourier transform the external vertex 
  // up to and including this momentum

  char filename[200];

  // configs
  for(int nc=0;nc<nconfs;nc++){
    
    // do the vacuum polarization, loop over ext. vert. insertion time slice
    for(t=alg_muon_arg->oper_time_start; t <= alg_muon_arg->oper_time_end; t++){
      
      // read in data
      qio_readMuon readVac(argc,argv);
      sprintf(filename,"%d/vacpol.t%d.dat",nc,t);
      readVac.read(filename,AlgLattice(),vacpol,16);

      // only need these 3 momenta
      for(i = 0; i <= XMOM; i++){
	for(j = 0; j <= YMOM; j++){
	  for(k = 0; k <= ZMOM; k++){

	    //skip zero momentum, and others...
	    int momsq = i*i+j*j+k*k;
	    if(momsq!=1)continue;

	    // init to 0 for glb sum
	    for(l = 0; l < TSize; l++){
	      for(rho=0;rho<4;rho++){
		for(mu=0;mu<4;mu++){
		  temp[mu+4*(rho+4*l)]=Zero;
		}
	      }
	    }
	    // need stuff for all time components of 4-momentum
	    for(l = 0; l < TSize; l++){

	      // this is node coordinate where momentum lives
	      int procCoorX = i / n[0];
	      int procCoorY = j / n[1];
	      int procCoorZ = k / n[2];
	      int procCoorT = l / n[3];
	      // this is local coordinate where momentum lives
	      int localX = i % n[0];
	      int localY = j % n[1];
	      int localZ = k % n[2];
	      int localT = l % n[3];
	      	      
	      // get stuff if we are on the right node
	      if (myXcoor == procCoorX &&
		  myYcoor == procCoorY &&
		  myZcoor == procCoorZ &&
		  myTcoor == procCoorT){
		
		int site = localX + n[0]*(localY +n[1]*(localZ + n[2]*localT));
		//external, internal lorentz index
		for(rho=0;rho<4;rho++){
		  for(mu=0;mu<4;mu++){
		    temp[mu+4*(rho+4*l)]=vacpol[mu+4*(rho+4*site)];
		  }
		}
	      }
	    }
	    // global sum (send it to all nodes, especially node 0)
	    slice_sum((Float*)temp,32*TSize,99);
	    for(l = 0; l < TSize; l++){
	      for(rho=0;rho<4;rho++){
		for(mu=0;mu<4;mu++){
		  printf("VACP conf= %d t_op= %d rho= %d mu= %d mom= %d %d %d %d  %e %e \n",
			 nc,t,rho,mu,i,j,k,l,
			 temp[mu+4*(rho+4*l)].real(),
			 temp[mu+4*(rho+4*l)].imag());
		}
	      } 
	    }
	  }
	}
      }//mom
    }//t op

    // now the muon line
    for(i = 0; i <= XMOM; i++){
      for(j = 0; j <= YMOM; j++){
	for(k = 0; k <= ZMOM; k++){
	  
	  //skip zero momentum, and others...
	  int momsq = i*i+j*j+k*k;
	  if(momsq!=1)continue;

	  // init to 0 for glb sum
	  for(l = 0; l < TSize; l++){
	    for(int proj=0;proj<NPROJ;proj++){
	      for(mu=0;mu<4;mu++){
		temp[mu+4*(proj+NPROJ*l)]=Zero;
	      }
	    }
	  }
	  
	  // get the Muonline
	  qio_readMuon readMu(argc,argv);
	  sprintf(filename,"%d/muonline.i%d%d%d.o%d%d%d.dat",nc,0,0,0,i,j,k);
	  readMu.read(filename,AlgLattice(),Muon,4*NPROJ);
	  
	  // need stuff for all time components of 4-momentum
	  for(l = 0; l < TSize; l++){
	    
	    // this is node coordinate where momentum lives
	    int procCoorX = i / n[0];
	    int procCoorY = j / n[1];
	    int procCoorZ = k / n[2];
	    int procCoorT = l / n[3];
	    // this is local coordinate where momentum lives
	    int localX = i % n[0];
	    int localY = j % n[1];
	    int localZ = k % n[2];
	    int localT = l % n[3];
	    
	    // print stuff out if we are one the right node
	    if (myXcoor == procCoorX &&
		myYcoor == procCoorY &&
		myZcoor == procCoorZ &&
		myTcoor == procCoorT){
	      
	      int site = localX + n[0]*(localY +n[1]*(localZ + n[2]*localT));
	      //projector, internal lorentz index
	      for(int proj=0;proj<NPROJ;proj++){
		for(mu=0;mu<4;mu++){
		  temp[mu+4*(proj+NPROJ*l)]=Muon[mu+4*(proj+NPROJ*site)];
		}
	      }
	    }
	  }
	  // global sum
	  slice_sum((Float*)temp,2*4*NPROJ*TSize,99);
	  for(l = 0; l < TSize; l++){
	    for(int proj=0;proj<NPROJ;proj++){
	      for(mu=0;mu<4;mu++){
		printf("Muon conf= %d proj= %d mu= %d mom= %d %d %d %d  %e %e \n",
		       nc,proj,mu,i,j,k,l,
		       temp[mu+4*(proj+NPROJ*l)].real(),
		       temp[mu+4*(proj+NPROJ*l)].imag()); 
	      }
	    }
	  }
	}
      }
    }// mom


  }//confs

  // clean up
  sfree(temp);
  sfree(vacpol);
  sfree(Muon);

}


void AlgMuon::run_jk_sub()
{

  char *fname = "run_jk_sub()";
  VRB.Func(cname,fname);
  int mu,nu,rho;
  // total number of configurations
  int nconfs = NC; //debug =1, fix this!
  int jksize = 1;  // omit 1 lattice at a time
  int NT = alg_muon_arg->oper_time_end-alg_muon_arg->oper_time_start+1; 
            // number of time slices for vertex
  int NMOM =3;
  int size = 4*GJP.VolNodeSites();
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  // holds muon line
  Rcomplex *Muon;
  Muon = (Rcomplex*)smalloc(NPROJ*NMOM*size*sizeof(Rcomplex));
  // holds Vacuum Polarization
  Rcomplex *vacpol;
  vacpol = (Rcomplex*)smalloc(4*NT*size*sizeof(Rcomplex));
  // holds data for global sum
  //int TSize = n[3]*GJP.Tnodes();
  // holds muon line
  Rcomplex *temp;
  temp = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));
  Rcomplex *temp2;
  temp2 = (Rcomplex*)smalloc(NPROJ*size*sizeof(Rcomplex));

  //----------------------------------------------------------------

  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();
  Float momX = 2.*PI/(GJP.Xnodes()*GJP.XnodeSites());
  Float momT = 2.*PI/(GJP.Tnodes()*GJP.TnodeSites());

  Rcomplex Zero(0.0,0.0);
  int t;
  int i,j,k;
  // fourier transform the external vertex 
  // up to and including this momentum

  int x[4];
  char filename[200];

  // configs
  Float norm=1.0/(nconfs-jksize);
  for(int nc=0;nc<nconfs;nc++){

    printf("CONFIG: %d\n",nc);

    // avg vacpol first
    int tt=-1;
    for(t=alg_muon_arg->oper_time_start; t <= alg_muon_arg->oper_time_end; t++){

      tt++;

      qio_readMuon readVac(argc,argv);
      sprintf(filename,"%d.b/vacpol.t%d.dat",nc,t);
      readVac.read(filename,AlgLattice(),temp,16);

      for(x[3] = 0; x[3] < n[3]; x[3]++){
	for(x[2] = 0; x[2] < n[2]; x[2]++){
	  for(x[1] = 0; x[1] < n[1]; x[1]++){
	    for(x[0] = 0; x[0] < n[0]; x[0]++){
	      int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*x[3]));
	      // internal Lorentz indices
	      for(nu=0;nu<4;nu++){
		for(mu=0;mu<4;mu++){
		  if(nc==0)
		    vacpol[nu+4*(mu+4*(tt+NT*site))]=temp[nu+4*(mu+4*site)]*norm;
		  else 
		    vacpol[nu+4*(mu+4*(tt+NT*site))]+=temp[nu+4*(mu+4*site)]*norm;
		}
	      }
	    }
	  }
	}
      }
    }

    // now the avg muonline
    int mom=0;
    for(i = 0; i <= XMOM; i++){
      for(j = 0; j <= YMOM; j++){
	for(k = 0; k <= ZMOM; k++){
	  
	  //skip zero momentum, and others...
	  int momsq = i*i+j*j+k*k;
	  if(momsq!=1)continue;

	  if(mom >= NMOM)ERR.General(cname,fname, "too many momenta!");
	  
	  // get the Muonline
	  qio_readMuon readMu(argc,argv);
	  sprintf(filename,"%d.b/muonline.i%d%d%d.o%d%d%d.dat",nc,0,0,0,i,j,k);
	  readMu.read(filename,AlgLattice(),temp2,4*NPROJ);
	  
	  // projector
	  for(int proj=0;proj<NPROJ;proj++){
	    for(x[3] = 0; x[3] < n[3]; x[3]++){
	      for(x[2] = 0; x[2] < n[2]; x[2]++){
		for(x[1] = 0; x[1] < n[1]; x[1]++){
		  for(x[0] = 0; x[0] < n[0]; x[0]++){
		    int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*x[3]));
		    for(mu=0;mu<4;mu++){
		      if(nc==0)
			Muon[mu+4*(proj+NPROJ*(mom+NMOM*site))] = 
			  temp2[mu+4*(proj+NPROJ*site)]*norm; 
		      else
			Muon[mu+4*(proj+NPROJ*(mom+NMOM*site))] += 
			  temp2[mu+4*(proj+NPROJ*site)]*norm; 
		    }
		  }
		}
	      }
	    }
	  }
	  mom++;
	}
      }
    }

  }//confs

  for(int nc=0;nc<nconfs;nc++){
    // read in data
    
    int tt=-1;
    for(t=alg_muon_arg->oper_time_start; t <= alg_muon_arg->oper_time_end; t++){

      tt++;

      qio_readMuon readVac(argc,argv);
      sprintf(filename,"%d.b/vacpol.t%d.dat",nc,t);
      readVac.read(filename,AlgLattice(),temp,16);

      int mom=0;
      for(i = 0; i <= XMOM; i++){
	for(j = 0; j <= YMOM; j++){
	  for(k = 0; k <= ZMOM; k++){

	    //skip zero momentum, and others...
	    int momsq = i*i+j*j+k*k;
	    if(momsq!=1)continue;
	    
	    if(mom == NMOM)ERR.General(cname,fname, "too many momenta!");

	    // get the Muonline
	    qio_readMuon readMu(argc,argv);
	    sprintf(filename,"%d.b/muonline.i%d%d%d.o%d%d%d.dat",nc,0,0,0,i,j,k);
	    readMu.read(filename,AlgLattice(),temp2,4*NPROJ);
	    
	    // projector
	    for(int proj=0;proj<NPROJ;proj++){
	      // ext. vertex lorentz index
	      for(rho=0;rho<4;rho++){
		
		//sum over q^2:
		Rcomplex sum=Zero;
		for(x[3] = 0; x[3] < n[3]; x[3]++){
		  int tg = x[3]+shiftT;
		  //Float pt = cos(tg*momT);
		  Float pt= sin(0.5*tg*momT);
		  for(x[2] = 0; x[2] < n[2]; x[2]++){
		    int zg = x[2]+shiftZ;
		    //Float pz = cos(zg*momX);
		    Float pz= sin(0.5*zg*momX);
		    for(x[1] = 0; x[1] < n[1]; x[1]++){
		      int yg = x[1]+shiftY;
		      //Float py = cos(yg*momX);
		      Float py= sin(0.5*yg*momX);
		      for(x[0] = 0; x[0] < n[0]; x[0]++){
			int xg = x[0]+shiftX;
			//Float px = cos(xg*momX);
			Float px= sin(0.5*xg*momX);

			int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
			//Float qsq = 8.0-2.0*(px+py+pz+pt);
			Float qsq = 4.0*(px*px+py*py+pz*pz+pt*pt);
			
			if(qsq<EPS)continue;//avoid zero momentum
			qsq = 1.0/qsq;
			
			// internal Lorentz indices
			for(nu=0;nu<4;nu++){
			  for(mu=0;mu<4;mu++){
				//for now, take mu=nu
			    if(mu!=nu)continue;
			    // jk block nc
			    Rcomplex ttvp = vacpol[nu+4*(rho+4*(tt+NT*site))]
			      -   temp[nu+4*(rho+4*site)]*norm;
			    Rcomplex ttmu = Muon[mu+4*(proj+NPROJ*(mom+NMOM*site))]
			      -temp2[mu+4*(proj+NPROJ*site)]*norm;
			    sum+= qsq*(ttmu*ttvp);
			  }
			}
		      }
		    }
		  }
		}
		// global sum
		slice_sum((Float*)&sum,2,99);
		printf("LBL jkblk= %d t_op= %d proj= %d rho= %d mom= %d %d %d  %e %e \n",
		       nc,t,proj,rho,i,j,k,sum.real(),sum.imag());		
	      }
	    }
	    mom++;
	  }
	}
      }
    }
  }//confs

  // clean up
  sfree(temp2);
  sfree(temp);
  sfree(vacpol);
  sfree(Muon);

}


// compute and save averages
// Jacknife done in step2()
int AlgMuon::run_jk_sub_step1()
{

  int cnt;
  char *fname = "run_jk_sub_step1()";
  VRB.Func(cname,fname);
  int mu,nu;
  // total number of configurations
  int nconfs = NC;
  int jksize = 1;  // omit 1 lattice at a time
  int NT = alg_muon_arg->oper_time_end-alg_muon_arg->oper_time_start+1; 
            // number of time slices for vertex
  int NMOM =3;
  int size = 4*GJP.VolNodeSites();
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  // holds muon line
  Rcomplex *Muon;
  Muon = (Rcomplex*)smalloc(NPROJ*NMOM*size*sizeof(Rcomplex));
  // holds Vacuum Polarization
  Rcomplex *vacpol;
  vacpol = (Rcomplex*)smalloc(4*NT*size*sizeof(Rcomplex));
  // holds Vacuum Polarization
  Rcomplex *temp;
  temp = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));
  // holds muon line
  Rcomplex *temp2;
  temp2 = (Rcomplex*)smalloc(NPROJ*size*sizeof(Rcomplex));
  //------------------------------------------------------------------

  int t;
  int i,j,k;
  // fourier transform the external vertex 
  // up to and including this momentum

  int x[4];
  char filename[200];


  if(alg_muon_arg->vp_kind!=CONSERVED_LOCAL && alg_muon_arg->vp_kind!=CONSERVED_LOCAL_LOWMODE ){
    ERR.General(cname,fname, "VacPol type %d invalid.\n", alg_muon_arg->vp_kind);
  }

  qio_readMuon readVac(argc,argv);

  if(alg_muon_arg->vp_kind==CONSERVED_LOCAL){
    cnt=0;
    for(int nc=alg_muon_arg->conf;nc<nconfs;nc++){
      
      printf("CONFIG: %d\n",nc);
      
      // avg vacpol first
      int tt=-1;
      for(t=alg_muon_arg->oper_time_start; t <= alg_muon_arg->oper_time_end; t+=tINC){
	tt++;
	//qio_readMuon readVac(argc,argv);
        sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.xyz%d%d%d.%d",
                DIRVP,VPTAG,
                alg_muon_arg->loop_mass,
                t,
                qp_arg->x,
                qp_arg->y,
                qp_arg->z,
                nc);
	readVac.read(filename,AlgLattice(),temp,16);
	
	for(x[3] = 0; x[3] < n[3]; x[3]++){
	  for(x[2] = 0; x[2] < n[2]; x[2]++){
	    for(x[1] = 0; x[1] < n[1]; x[1]++){
	      for(x[0] = 0; x[0] < n[0]; x[0]++){
		int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*x[3]));
		// internal Lorentz indices
		for(nu=0;nu<4;nu++){
		  for(mu=0;mu<4;mu++){
		    if(nc==0)
		      vacpol[nu+4*(mu+4*(tt+NT*site))]=temp[nu+4*(mu+4*site)];
		    else 
		      vacpol[nu+4*(mu+4*(tt+NT*site))]+=temp[nu+4*(mu+4*site)];
		  }
		}
	      }
	    }
	  }
	}
      }
      cnt++;
    }
    // write it out 
    qio_writeMuon writeVac(argc,argv);
    sprintf(filename,"%s/sum-vacpol%s.dat",DIRVP,VPTAG);
    writeVac.write(filename,AlgLattice(),vacpol,NT*16);
  }else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOWMODE){

    int mom=0;
    for(i = 0; i <= XMOM; i++){
      for(j = 0; j <= YMOM; j++){
	for(k = 0; k <= ZMOM; k++){
	  
	  //skip zero momentum, and others...
	  int momsq = i*i+j*j+k*k;
	  if(momsq!=1)continue;
	  
	  if(mom >= NMOM)ERR.General(cname,fname, "too many momenta!");
	  
	  cnt=0;	  
	  for(int nc=0;nc<nconfs;nc++){
	    
	    printf("CONFIG: %d\n",nc);
	    
	    // avg vacpol first
	    int tt=-1;
	    for(t=alg_muon_arg->oper_time_start; t <= alg_muon_arg->oper_time_end; t+=tINC){
	      tt++;
	      //qio_readMuon readVac(argc,argv);
	      sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.mom%d%d%d.dat.%d",
		      DIRVP,VPTAG,
		      alg_muon_arg->loop_mass,
		      t,
		      i,j,k,
		      nc);
	      readVac.read(filename,AlgLattice(),temp,16);
	      
	      for(x[3] = 0; x[3] < n[3]; x[3]++){
		for(x[2] = 0; x[2] < n[2]; x[2]++){
		  for(x[1] = 0; x[1] < n[1]; x[1]++){
		    for(x[0] = 0; x[0] < n[0]; x[0]++){
		      int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*x[3]));
		      // internal Lorentz indices
		      for(nu=0;nu<4;nu++){
			for(mu=0;mu<4;mu++){
			  if(nc==0)
			    vacpol[nu+4*(mu+4*(tt+NT*site))]=temp[nu+4*(mu+4*site)];
			  else 
			    vacpol[nu+4*(mu+4*(tt+NT*site))]+=temp[nu+4*(mu+4*site)];
			}
		      }
		    }
		  }
		}
	      }
	    }
	    cnt++;
	  }
	  // write it out 
	  qio_writeMuon writeVac(argc,argv);
	  sprintf(filename,"%s/sum-vacpol%s.mom%d%d%d.dat",DIRVP,VPTAG,i,j,k);
	  writeVac.write(filename,AlgLattice(),vacpol,NT*16);
	}
      }
    }
  }
  // now the avg muonline
  for(int nc=0;nc<nconfs;nc++){

    printf("(LINE)CONFIG: %d\n",nc);

    int mom=0;
    for(i = 0; i <= XMOM; i++){
      for(j = 0; j <= YMOM; j++){
	for(k = 0; k <= ZMOM; k++){
	  
	  //skip zero momentum, and others...
	  int momsq = i*i+j*j+k*k;
	  if(momsq!=1)continue;

	  if(mom >= NMOM)ERR.General(cname,fname, "too many momenta!");
	  
	  // get the Muonline
	  qio_readMuon readMu(argc,argv);
	  sprintf(filename,"%s/muonline-consloc.m%g.i000.o%d%d%d.dat.%d",
		  DIRML,
		  alg_muon_arg->line_mass,
		  i,j,k,
		  nc);
	  readMu.read(filename,AlgLattice(),temp2,4*NPROJ);
	  
	  // projector
	  for(x[3] = 0; x[3] < n[3]; x[3]++){
	    for(x[2] = 0; x[2] < n[2]; x[2]++){
	      for(x[1] = 0; x[1] < n[1]; x[1]++){
		for(x[0] = 0; x[0] < n[0]; x[0]++){
		  int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*x[3]));
		  for(mu=0;mu<4;mu++){
		    for(int proj=0;proj<NPROJ;proj++){

		      if(nc==0)
			Muon[mu+4*(proj+NPROJ*(mom+NMOM*site))] = 
			  temp2[mu+4*(proj+NPROJ*site)]; 
		      else
			Muon[mu+4*(proj+NPROJ*(mom+NMOM*site))] += 
			  temp2[mu+4*(proj+NPROJ*site)]; 
		    }
		  }
		}
	      }
	    }
	  }
	  mom++;
	}
      }
    }
  }//confs
  qio_writeMuon writeMuon(argc,argv);
  sprintf(filename,"%s/sum-muon.dat",DIRML);
  writeMuon.write(filename,AlgLattice(),Muon,NMOM*4*NPROJ);
  // clean up
  sfree(temp2);
  sfree(temp);
  sfree(vacpol);
  sfree(Muon);
  return(cnt);
}




// compute and save averages
// Jacknife done in step2()
// calcs sums of vacpol, lines, returns # of configs
int AlgMuon::run_jk_sub_step1_mom()
{

  char *fname = "run_jk_sub_step1_mom()";
  VRB.Func(cname,fname);
  int mu,nu;
  // total number of configurations
  int nconfs = NC;
  int jksize = 1;  // omit 1 lattice at a time
  int NT = alg_muon_arg->oper_time_end-alg_muon_arg->oper_time_start+1; 
            // number of time slices for vertex
  int size = 4*GJP.VolNodeSites();
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  // holds muon line
  Rcomplex *Muon;
  Muon = (Rcomplex*)smalloc(NPROJ*Nmom*size*sizeof(Rcomplex));
  // holds Vacuum Polarization
  Rcomplex *vacpol;
  vacpol = (Rcomplex*)smalloc(4*NT*Nmom*size*sizeof(Rcomplex));
  // holds Vacuum Polarization
  Rcomplex *temp;
  temp = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));
  // holds muon line
  Rcomplex *temp2;
  temp2 = (Rcomplex*)smalloc(NPROJ*size*sizeof(Rcomplex));
  //------------------------------------------------------------------

  int t;
  int i,j,k;
  Rcomplex Zero(0.0,0.0);
  // fourier transform the external vertex 
  // up to and including this momentum

  int x[4];
  char filename[200];

  for(i = 0; i < Nmom; i++){
    for(int proj=0;proj<NPROJ;proj++){
      for(x[3] = 0; x[3] < n[3]; x[3]++){
        for(x[2] = 0; x[2] < n[2]; x[2]++){
          for(x[1] = 0; x[1] < n[1]; x[1]++){
	    for(x[0] = 0; x[0] < n[0]; x[0]++){
	      int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*x[3]));
	      for(mu=0;mu<4;mu++){
	        Muon[mu+4*(proj+NPROJ*(i+Nmom*site))] = Zero; 
	      }
	    }
          }
        }
      }
    }
  }
  
  int cnt=0;
  for(int nc=0;nc<nconfs;nc++){

    for(i = 0; i < Nmom; i++){

	  // get the Muonline
	  qio_readMuon readMu(argc,argv);
	  sprintf(filename,"%s/muonline-consloc.m%g.i000.o%d%d%d.dat.%d",
		  DIRML,
		  alg_muon_arg->line_mass,
		  ext_mom[i].cmp(0),
		  ext_mom[i].cmp(1),
		  ext_mom[i].cmp(2),
		  nc);
	  if(FileExists(filename)){
	    readMu.read(filename,AlgLattice(),temp2,4*NPROJ);
	  }else{
	    if(!UniqueID())printf("%s does not exist, skipping\n",filename);
	    continue;
	  }

	  // avg vacpol first
	  int tt=-1;
	  int skip=0;
	  for(t=alg_muon_arg->oper_time_start; t <= alg_muon_arg->oper_time_end; t+=tINC){
	    
	    tt++;
	    
	    if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOWMODE ){
	      qio_readMuon readVac(argc,argv);
	      sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.mom%d%d%d.%s.dat.%d",
		      DIRVP,VPTAG,
		      alg_muon_arg->loop_mass,
		      t,
		      ext_mom[i].cmp(0),
		      ext_mom[i].cmp(1),
		      ext_mom[i].cmp(2),
		      EIGTAG,
		      nc);
	      if(FileExists(filename)){
		readVac.read(filename,AlgLattice(),temp,16);
	      }else{
	        if(!UniqueID())printf("%s does not exist, skipping\n",filename);
		skip=1;
		continue;
	      }
	    }else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_PTSRC_TEST || alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOOP_PTSRC){
	      qio_readMuon readVac(argc,argv);
	      sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.mom%d%d%d.dat.%d",
		      DIRVP,VPTAG,
		      alg_muon_arg->loop_mass,
		      t,
		      ext_mom[i].cmp(0),
                      ext_mom[i].cmp(1),
                      ext_mom[i].cmp(2),
                      nc);
	      if(FileExists(filename)){
		readVac.read(filename,AlgLattice(),temp,16);
	      }else{
	        if(!UniqueID())printf("%s does not exist, skipping\n",filename);
		skip=1;
		continue;
	      }
	    }else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL){
	      qio_readMuon readVac(argc,argv);
	      sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.dat.%d",
		      DIRVP,VPTAG,
		      alg_muon_arg->loop_mass,
		      t,
		      nc);
	      if(FileExists(filename)){
		readVac.read(filename,AlgLattice(),temp,16);
	      }else{
	        if(!UniqueID())printf("%s does not exist, skipping\n",filename);
		skip=1;
		continue;
	      }
	    }

	    for(x[3] = 0; x[3] < n[3]; x[3]++){
	      for(x[2] = 0; x[2] < n[2]; x[2]++){
		for(x[1] = 0; x[1] < n[1]; x[1]++){
		  for(x[0] = 0; x[0] < n[0]; x[0]++){
		    int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*x[3]));
		    // internal Lorentz indices
		    for(nu=0;nu<4;nu++){
		      for(mu=0;mu<4;mu++){
			vacpol[nu+4*(mu+4*(tt+NT*(i+Nmom*site)))]+=temp[nu+4*(mu+4*site)];
		      }
		    }
		  }
		}
	      }
	    }

	  } // t_op

          if(skip==0){	    
	    // projector
	    for(int proj=0;proj<NPROJ;proj++){
	      for(x[3] = 0; x[3] < n[3]; x[3]++){
		for(x[2] = 0; x[2] < n[2]; x[2]++){
		  for(x[1] = 0; x[1] < n[1]; x[1]++){
		    for(x[0] = 0; x[0] < n[0]; x[0]++){
		      int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*x[3]));
		      for(mu=0;mu<4;mu++){
			  Muon[mu+4*(proj+NPROJ*(i+Nmom*site))] += temp2[mu+4*(proj+NPROJ*site)]; 
		      }
		    }
		  }
		}
	      }
	    }
	    cnt++;
	    if(!UniqueID())printf("CONFIG: %d\n",nc);
	  }

    }
  }//confs

  cnt/=Nmom;
  if(cnt==0){
     ERR.General(cname,fname, "no data to jack!\n");
  }

  // write it out 
  if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOWMODE ){
    qio_writeMuon writeVac(argc,argv);
    sprintf(filename,"%s/sum-vacpol%s%s.dat",DIRVP,VPTAG,EIGTAG);
    writeVac.write(filename,AlgLattice(),vacpol,NT*Nmom*16);
  }else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_PTSRC_TEST || alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOOP_PTSRC){
    qio_writeMuon writeVac(argc,argv);
    sprintf(filename,"%s/sum-vacpol%s.dat",DIRVP,VPTAG);
    writeVac.write(filename,AlgLattice(),vacpol,NT*Nmom*16);
  }
  // write it out 
  qio_writeMuon writeMuon(argc,argv);
  sprintf(filename,"%s/sum-muon.dat",DIRML);
  writeMuon.write(filename,AlgLattice(),Muon,Nmom*4*NPROJ);
  // clean up
  sfree(temp2);
  sfree(temp);
  sfree(vacpol);
  sfree(Muon);
  return(cnt);

}

// read averages, do jacknife
void AlgMuon::run_jk_sub_step2(int cnt)
{

  char *fname = "run_jk_sub_step2(int)";
  VRB.Func(cname,fname);
  int mu,nu,rho;
  // total number of configurations
  int nconfs = NC;
  int jksize = 1;  // omit 1 lattice at a time
  int NT = alg_muon_arg->oper_time_end-alg_muon_arg->oper_time_start+1; 
            // number of time slices for vertex
  int NMOM =3;
  int size = 4*GJP.VolNodeSites();
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  // holds muon line
  Rcomplex *Muon;
  Muon = (Rcomplex*)smalloc(NPROJ*NMOM*size*sizeof(Rcomplex));
  // holds Vacuum Polarization
  Rcomplex *vacpol;
  vacpol = (Rcomplex*)smalloc(4*NT*size*sizeof(Rcomplex));
  // holds data for global sum
  //int TSize = n[3]*GJP.Tnodes();
  // holds muon line
  Rcomplex *temp;
  temp = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));
  Rcomplex *temp2;
  temp2 = (Rcomplex*)smalloc(NPROJ*size*sizeof(Rcomplex));

  //----------------------------------------------------------------

  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();
  Float momX = 2.*PI/(GJP.Xnodes()*GJP.XnodeSites());
  Float momT = 2.*PI/(GJP.Tnodes()*GJP.TnodeSites());

  Rcomplex Zero(0.0,0.0);
  int t;
  int i,j,k;
  // fourier transform the external vertex 
  // up to and including this momentum

  int x[4];
  char filename[200];

  // this is my node coordinate
  int myXcoor = GJP.XnodeCoor();
  int myYcoor = GJP.YnodeCoor();
  int myZcoor = GJP.ZnodeCoor();
  int myTcoor = GJP.TnodeCoor();

  if(alg_muon_arg->vp_kind!=CONSERVED_LOCAL){
    ERR.General(cname,fname, "VacPol type %d invalid.\n", alg_muon_arg->vp_kind);
  }

  // read in the averages
  if(alg_muon_arg->vp_kind==CONSERVED_LOCAL){
    qio_readMuon readVac(argc,argv);
    sprintf(filename,"%s/sum-vacpol%s.dat",DIRVP,VPTAG);
    readVac.read(filename,AlgLattice(),vacpol,NT*16);
  }
  qio_readMuon readMuon(argc,argv);
  sprintf(filename,"%s/sum-muon.dat",DIRML);
  readMuon.read(filename,AlgLattice(),Muon,NMOM*4*NPROJ);

  // # of configs
  Float norm=1.0/((Float)(cnt-jksize));

  for(int nc=alg_muon_arg->conf;nc<nconfs;nc++){
    int tt=-1;
    for(t=alg_muon_arg->oper_time_start; t <= alg_muon_arg->oper_time_end; t+=tINC){

      tt++;

      qio_readMuon readVac(argc,argv);
      if(alg_muon_arg->vp_kind==CONSERVED_LOCAL){
        sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.xyz%d%d%d.%d",
                DIRVP,VPTAG,
                alg_muon_arg->loop_mass,
                t,
                qp_arg->x,
                qp_arg->y,
                qp_arg->z,
                nc);
	readVac.read(filename,AlgLattice(),temp,16);
      }

      int mom=0;
      for(i = 0; i <= XMOM; i++){
	for(j = 0; j <= YMOM; j++){
	  for(k = 0; k <= ZMOM; k++){

	    //skip zero momentum, and others...
	    int momsq = i*i+j*j+k*k;
	    if(momsq!=1)continue;
	    
	    if(mom == NMOM)ERR.General(cname,fname, "too many momenta!");
	    // get the Muonline
	    qio_readMuon readMu(argc,argv);
	    sprintf(filename,"%s/muonline-consloc.m%g.i000.o%d%d%d.dat.%d",
		    DIRML,
		    alg_muon_arg->line_mass,
		    i,j,k,
		    nc);
	    readMu.read(filename,AlgLattice(),temp2,4*NPROJ);

	    // projector
	    for(int proj=0;proj<NPROJ;proj++){
	      // ext. vertex lorentz index
	      for(rho=0;rho<4;rho++){
		
		//sum over q^2:
		Rcomplex sum=Zero;
		for(x[3] = 0; x[3] < n[3]; x[3]++){
		  int tg = x[3]+shiftT;
		  //Float pt = cos(tg*momT);
		  Float pt= sin(0.5*tg*momT);
		  for(x[2] = 0; x[2] < n[2]; x[2]++){
		    int zg = x[2]+shiftZ;
		    //Float pz = cos(zg*momX);
		    Float pz= sin(0.5*zg*momX);
		    for(x[1] = 0; x[1] < n[1]; x[1]++){
		      int yg = x[1]+shiftY;
		      //Float py = cos(yg*momX);
		      Float py= sin(0.5*yg*momX);
		      for(x[0] = 0; x[0] < n[0]; x[0]++){
			int xg = x[0]+shiftX;
			//Float px = cos(xg*momX);
			Float px= sin(0.5*xg*momX);

			int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
			//Float qsq = 8.0-2.0*(px+py+pz+pt);
			Float qsq = 4.0*(px*px+py*py+pz*pz+pt*pt);
			
			if(qsq<EPS)continue;//avoid zero momentum
			qsq = 1.0/qsq;
			
			// internal Lorentz indices
			for(nu=0;nu<4;nu++){
			  for(mu=0;mu<4;mu++){
			    //for now, take mu=nu
			    if(mu!=nu)continue;
			    // new convention for cons-loc vacpol computed after 6/26/08!
			    Rcomplex ttvp = (vacpol[nu+4*(rho+4*(tt+NT*site))]
					     - temp[nu+4*(rho+4*site)])*norm;
			    Rcomplex ttmu = (Muon[mu+4*(proj+NPROJ*(mom+NMOM*site))]
					     - temp2[mu+4*(proj+NPROJ*site)])*norm;
			    sum+= qsq*(ttmu*ttvp);
			  }
			}
		      }
		    }
		  }
		}
		// global sum
		slice_sum((Float*)&sum,2,99);
		if(!UniqueID())printf("LBL jkblk= %d t_op= %d proj= %d rho= %d mom= %d %d %d  %.15e %.15e \n",
		       nc,t,proj,rho,i,j,k,sum.real(),sum.imag());		
	      }
	    }
	    mom++;
	  }
	}
      }
    }
  }//confs

  // clean up
  sfree(temp2);
  sfree(temp);
  sfree(vacpol);
  sfree(Muon);

}

void AlgMuon::run_jk_sub_step2(int* out_mom,
			       int t_op,
			       Rcomplex* Muon, 
			       Rcomplex* vacpol,
			       Rcomplex* avgMuon, 
			       Rcomplex* avgVac)
{

  char *fname = "run_jk_sub_step2(out_mom,t,Muon,vacpol)";
  VRB.Func(cname,fname);
  int mu,nu,rho;
  // total number of configurations
  int nconfs = NC;
  int nc  = alg_muon_arg->conf;
  int jksize = 1;  // omit 1 lattice at a time
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();

  //----------------------------------------------------------------

  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();
  Float momX = 2.*PI/(GJP.Xnodes()*GJP.XnodeSites());
  Float momT = 2.*PI/(GJP.Tnodes()*GJP.TnodeSites());

  Rcomplex Zero(0.0,0.0);
  int t;
  int i,j,k;
  // fourier transform the external vertex 
  // up to and including this momentum

  int x[4];

  Float norm=1.0;
  if(nconfs > 1){
    norm=1.0/(nconfs-jksize);
  }

  // projector
  for(int proj=0;proj<NPROJ;proj++){
    // ext. vertex lorentz index
    for(rho=0;rho<4;rho++){
      
      //sum over q^2:
      Rcomplex sum=Zero;
      for(x[3] = 0; x[3] < n[3]; x[3]++){
	int tg = x[3]+shiftT;
	Float pt= sin(0.5*tg*momT);
	for(x[2] = 0; x[2] < n[2]; x[2]++){
	  int zg = x[2]+shiftZ;
	  Float pz= sin(0.5*zg*momX);
	  for(x[1] = 0; x[1] < n[1]; x[1]++){
	    int yg = x[1]+shiftY;
	    Float py= sin(0.5*yg*momX);
	    for(x[0] = 0; x[0] < n[0]; x[0]++){
	      int xg = x[0]+shiftX;
	      Float px= sin(0.5*xg*momX);
	      
	      int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
	      //Float qsq = 8.0-2.0*(px+py+pz+pt);
	      Float qsq = 4.0*(px*px+py*py+pz*pz+pt*pt);
	      
	      if(qsq<EPS)continue;//avoid zero momentum
	      qsq = 1.0/qsq;
	      
	      // internal Lorentz indices
	      for(nu=0;nu<4;nu++){
		for(mu=0;mu<4;mu++){
		  //for now, take mu=nu
		  if(mu!=nu)continue;
		  // jk block nc
		  Rcomplex ttvp = (avgVac[nu+4*(rho+4*site)]
				   - vacpol[nu+4*(rho+4*site)])*norm;
		  Rcomplex ttmu = (avgMuon[mu+4*(proj+NPROJ*site)]
				   -  Muon[mu+4*(proj+NPROJ*site)])*norm;
		  sum+= qsq*(ttmu*ttvp);
		}
	      }
	    }
	  }
	}
      }
      // global sum
      slice_sum((Float*)&sum,2,99);
      if(!UniqueID())printf("LBL jkblk= %d t_op= %d proj= %d rho= %d mom= %d %d %d  %e %e \n",
	     nc,t_op,proj,rho,
	     out_mom[0],out_mom[1],out_mom[2],
	     sum.real(),sum.imag());		
    }
  }
}



void AlgMuon::run_jk_sub_step2_mom(int cnt)
{

  char *fname = "run_jk_sub_step2_mom()";
  VRB.Func(cname,fname);

  if(cnt==0){
     ERR.General(cname,fname, "no data to jack!\n");
  }else if(cnt==1){
     ERR.General(cname,fname, "only 1 block-- no jack\n");
  }

  int mu,nu,rho;
  // total number of configurations
  int nconfs = NC;
  int jksize = 1;  // omit 1 lattice at a time
  int NT = alg_muon_arg->oper_time_end-alg_muon_arg->oper_time_start+1; 
            // number of time slices for vertex
  int size = 4*GJP.VolNodeSites();
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  // holds muon line
  Rcomplex *Muon;
  Muon = (Rcomplex*)smalloc(NPROJ*Nmom*size*sizeof(Rcomplex));
  // holds Vacuum Polarization
  Rcomplex *vacpol;
  vacpol = (Rcomplex*)smalloc(4*NT*Nmom*size*sizeof(Rcomplex));
  // holds data for global sum
  //int TSize = n[3]*GJP.Tnodes();
  // holds muon line
  Rcomplex *temp;
  temp = (Rcomplex*)smalloc(4*size*sizeof(Rcomplex));
  Rcomplex *temp2;
  temp2 = (Rcomplex*)smalloc(NPROJ*size*sizeof(Rcomplex));

  //----------------------------------------------------------------

  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();
  Float momX = 2.*PI/(GJP.Xnodes()*GJP.XnodeSites());
  Float momT = 2.*PI/(GJP.Tnodes()*GJP.TnodeSites());

  Rcomplex Zero(0.0,0.0);
  int t;
  int i,j,k;
  // fourier transform the external vertex 
  // up to and including this momentum

  int x[4];
  char filename[200];

  // this is my node coordinate
  int myXcoor = GJP.XnodeCoor();
  int myYcoor = GJP.YnodeCoor();
  int myZcoor = GJP.ZnodeCoor();
  int myTcoor = GJP.TnodeCoor();

  // read in the averages
  if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOWMODE ){
    qio_readMuon readVac(argc,argv);
    sprintf(filename,"%s/sum-vacpol%s%s.dat",DIRVP,VPTAG,EIGTAG);
    readVac.read(filename,AlgLattice(),vacpol,NT*Nmom*16);
  }else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_PTSRC_TEST||alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOOP_PTSRC){
    qio_readMuon readVac(argc,argv);
    sprintf(filename,"%s/sum-vacpol%s.dat",DIRVP,VPTAG);
    readVac.read(filename,AlgLattice(),vacpol,NT*Nmom*16);
  }else{
    ERR.General(cname,fname, "VacPol type %d invalid.\n", alg_muon_arg->vp_kind);
  }
  // 
  qio_readMuon readMuon(argc,argv);
  sprintf(filename,"%s/sum-muon.dat",DIRML);
  readMuon.read(filename,AlgLattice(),Muon,Nmom*4*NPROJ);

  // configs
  Float norm=1.0/((Float)(cnt-jksize));

  for(int nc=alg_muon_arg->conf;nc<nconfs;nc++){
    // read in data
    int tt=-1;
    for(t=alg_muon_arg->oper_time_start; t <= alg_muon_arg->oper_time_end; t+=tINC){

      tt++;

      for(i = 0; i < Nmom; i++){

	    if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOWMODE){
	      qio_readMuon readVac(argc,argv);
	      sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.mom%d%d%d.%s.dat.%d",
		      DIRVP,VPTAG,
		      alg_muon_arg->loop_mass,
		      t,
		      ext_mom[i].cmp(0),
		      ext_mom[i].cmp(1),
		      ext_mom[i].cmp(2),
		      EIGTAG,
		      nc);
	      if(FileExists(filename)){
		readVac.read(filename,AlgLattice(),temp,16);
	      }else{
                if(!UniqueID())printf("%s does not exist, skipping\n",filename);
		continue;
	      }
	    }else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL_PTSRC_TEST || alg_muon_arg->vp_kind==CONSERVED_LOCAL_LOOP_PTSRC){
	      qio_readMuon readVac(argc,argv);
	      sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.mom%d%d%d.dat.%d",
		      DIRVP,VPTAG,
		      alg_muon_arg->loop_mass,
		      t,
		      ext_mom[i].cmp(0),
		      ext_mom[i].cmp(1),
		      ext_mom[i].cmp(2),
		      nc);
	      if(FileExists(filename)){
		readVac.read(filename,AlgLattice(),temp,16);
	      }else{
                if(!UniqueID())printf("%s does not exist, skipping\n",filename);
		continue;
	      }
	    }else if(alg_muon_arg->vp_kind==CONSERVED_LOCAL){
	      qio_readMuon readVac(argc,argv);
	      sprintf(filename,"%s/vacpol-consloc%s.m%g.t%d.dat.%d",
		      DIRVP,VPTAG,
		      alg_muon_arg->loop_mass,
		      t,
		      nc);
	      if(FileExists(filename)){
		readVac.read(filename,AlgLattice(),temp,16);
	      }else{
                if(!UniqueID())printf("%s does not exist, skipping\n",filename);
		continue;
	      }
	    }

	    // get the Muonline
	    qio_readMuon readMu(argc,argv);
	    sprintf(filename,"%s/muonline-consloc.m%g.i000.o%d%d%d.dat.%d",
		    DIRML,
		    alg_muon_arg->line_mass,
		    ext_mom[i].cmp(0),
		    ext_mom[i].cmp(1),
		    ext_mom[i].cmp(2),
		    nc);
	    if(FileExists(filename)){
	      readMu.read(filename,AlgLattice(),temp2,4*NPROJ);
	    }else{
              if(!UniqueID())printf("%s does not exist, skipping\n",filename);
	      continue;
	    }
	    
	    // projector
	    for(int proj=0;proj<NPROJ;proj++){
	      // ext. vertex lorentz index
	      for(rho=0;rho<4;rho++){
		
		//sum over q^2:
		Rcomplex sum=Zero;
		for(x[3] = 0; x[3] < n[3]; x[3]++){
		  int tg = x[3]+shiftT;
		  //Float pt = cos(tg*momT);
		  Float pt= sin(0.5*tg*momT);
		  for(x[2] = 0; x[2] < n[2]; x[2]++){
		    int zg = x[2]+shiftZ;
		    //Float pz = cos(zg*momX);
		    Float pz= sin(0.5*zg*momX);
		    for(x[1] = 0; x[1] < n[1]; x[1]++){
		      int yg = x[1]+shiftY;
		      //Float py = cos(yg*momX);
		      Float py= sin(0.5*yg*momX);
		      for(x[0] = 0; x[0] < n[0]; x[0]++){
			int xg = x[0]+shiftX;
			//Float px = cos(xg*momX);
			Float px= sin(0.5*xg*momX);

			int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
			//Float qsq = 8.0-2.0*(px+py+pz+pt);
			Float qsq = 4.0*(px*px+py*py+pz*pz+pt*pt);
			
			if(qsq<EPS)continue;//avoid zero momentum
			qsq = 1.0/qsq;
			
			// internal Lorentz indices
			for(nu=0;nu<4;nu++){
			  for(mu=0;mu<4;mu++){
			    //for now, take mu=nu
			    if(mu!=nu)continue;
			    // new convention for cons-loc vacpol computed after 6/26/08!
			    Rcomplex ttvp = (vacpol[nu+4*(rho+4*(tt+NT*(i+Nmom*site)))]
					     - temp[nu+4*(rho+4*site)])*norm;
			    Rcomplex ttmu = (Muon[mu+4*(proj+NPROJ*(i+Nmom*site))]
					     - temp2[mu+4*(proj+NPROJ*site)])*norm;
			    sum+= qsq*(ttmu*ttvp);
			  }
			}
		      }
		    }
		  }
		}
		// global sum
		slice_sum((Float*)&sum,2,99);
		if(!UniqueID())printf("LBL jkblk= %d t_op= %d proj= %d rho= %d mom= %d %d %d  %.15e %.15e \n",
		       nc,t,proj,rho,ext_mom[i].cmp(0),ext_mom[i].cmp(1),ext_mom[i].cmp(2),sum.real(),sum.imag());
	      }
	    }
      }
    }
  }//confs

  // clean up
  sfree(temp2);
  sfree(temp);
  sfree(vacpol);
  sfree(Muon);

}


void AlgMuon::two_point()
{

  char *fname = "two_point()";
  VRB.Func(cname,fname);
  WilsonMatrix temp;

  // set up for fourier transforms
  //----------------------------------------------------------------
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  int nx = n[0]*GJP.Xnodes();
  int ny = n[1]*GJP.Ynodes();
  int nz = n[2]*GJP.Znodes();
  int nt = n[3]*GJP.Tnodes();
  // this is my node coordinate
  int myXcoor = GJP.XnodeCoor();
  int myYcoor = GJP.YnodeCoor();
  int myZcoor = GJP.ZnodeCoor();
  int myTcoor = GJP.TnodeCoor();

  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();

  int t;
  int i,j,k;
  int x[4];
 
  Rcomplex ZERO(0.0,0.0);
  Rcomplex I(0.0,1.0);

  // two-point functions
  //----------------------------------------------------------------
  Rcomplex *two_pt = (Rcomplex*) smalloc(GJP.VolNodeSites()*sizeof(Rcomplex));
  Rcomplex *two_pt_np = (Rcomplex*) smalloc(GJP.VolNodeSites()*sizeof(Rcomplex));

  // get the prop
  //----------------------------------------------------------------
  int *mom = (int*)smalloc(3*sizeof(int));
  mom[0]=mom[1]=mom[2]=0;
  // incoming muon prop, only need it once
  qp_arg->cg.mass = alg_muon_arg->line_mass ;
  qp_arg->t =  0  ;
  //QPropWMomSrc in(AlgLattice(), qp_arg, mom, common_arg);
  // test for normalization:
  qp_arg->x =  0;
  qp_arg->y =  0;
  qp_arg->z =  0;
  QPropWPointSrc in(AlgLattice(), qp_arg, common_arg);

  // do the incoming muon two-point function (trace the prop)
  for(x[3] = 0; x[3] < n[3]; x[3]++){
    for(x[2] = 0; x[2] < n[2]; x[2]++){
      for(x[1] = 0; x[1] < n[1]; x[1]++){
	for(x[0] = 0; x[0] < n[0]; x[0]++){
	  
	  int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
          
	  // two pt functions (muon propagator)
	  //project with (1+-gamma_t/2)
	  temp = in[site];
	  temp.PParProjectSink();
	  two_pt[site] = temp.Trace() ;
	  temp = in[site];
	  temp.NParProjectSink();
	  two_pt_np[site] = temp.Trace() ;
	}
      }
    }
  }
  // Fourier transform spatial comps. (2 is size in Floats per site)
  // Since the source is zero mom, doesn't matter which, 0,1 is chosen
  // for "sign" of fft (last arg)
  FFTXYZ((Float*)two_pt, 1, 2, 1, 1);	
  FFTXYZ((Float*)two_pt_np, 1, 2, 1, 1);	
  // Print out results
  //------------------
  FileIoType T=ADD_ID;
  char filename[100];
  sprintf(filename,"prop-in.exmom000.dat");
  FILE *fp=Fopen(T, filename, "w");
  x[2] = 0;
  x[1] = 0;
  x[0] = 0;
  if(myXcoor==0 && myYcoor==0 && myZcoor==0){

    for(x[3] = 0; x[3] < n[3]; x[3]++){
      int tg = x[3]+shiftT;
      int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
      Fprintf(fp,"muon 2pt FUNCTION (+-parity proj) tsrc= %d tsink= %d mom= %d %d %d %e %e  %e %e\n", 
	      qp_arg->t, tg, 0, 0, 0,
	      two_pt[site].real(), 
	      two_pt[site].imag(),
	      two_pt_np[site].real(), 
	      two_pt_np[site].imag());
    }
  }
  Fclose(T,fp);
  
  // source location of the out-going muon
  qp_arg->t =  alg_muon_arg->source_time ;
  qp_arg->x = 0 ; // fix this!
  qp_arg->y = 0 ; // fix this!
  qp_arg->z = 0 ; // fix this!

  // Loop over momenta
  for(i = 0; i <= XMOM; i++){
    for(j = 0; j <= YMOM; j++){
      for(k = 0; k <= ZMOM; k++){
	
	if(i*i+j*j+k*k != 1){continue;}
        
	// outgoing muon with momentum p
	mom[0]=i; mom[1]=j; mom[2]=k;
	if(!UniqueID())printf("PROP: mom= %d %d %d\n",mom[0],mom[1],mom[2]);
	//QPropWMomSrc out(AlgLattice(), qp_arg, mom, common_arg);
	QPropWPointSrc out(AlgLattice(), qp_arg, common_arg);
	
	// loop over sink position
	for(x[3] = 0; x[3] < n[3]; x[3]++){
	  for(x[2] = 0; x[2] < n[2]; x[2]++){
	    for(x[1] = 0; x[1] < n[1]; x[1]++){
	      for(x[0] = 0; x[0] < n[0]; x[0]++){
		
		int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
                
		// two pt functions (muon propagator)
		//project with (1+-gamma_t/2)
		temp = out[site];
		temp.PParProjectSink();
		two_pt[site] = temp.Trace() ;
		temp = out[site];
		temp.NParProjectSink();
		two_pt_np[site] = temp.Trace() ;
	      }
	    }
	  }
	}
	// Fourier transform spatial comps. (2 is size in Floats per site)
	// since we will take anti-muon -> muon, flip sign in fft, or 0 last arg
	FFTXYZ((Float*)two_pt, 1, 2, 1, 0);	
	FFTXYZ((Float*)two_pt_np, 1, 2, 1, 0);	
	// Print out results
	//------------------
	FileIoType T=ADD_ID;
	char filename[100];
	sprintf(filename,"prop-out.exmom%d%d%d.dat",i,j,k);
	FILE *fp=Fopen(T, filename, "w");
	for(x[3] = 0; x[3] < n[3]; x[3]++){
	  int tg = x[3]+shiftT;
	  for(x[2] = 0; x[2] < n[2]; x[2]++){
	    int pzg = x[2]+shiftZ;
	    if(pzg != k && pzg != nz-1){continue;}
	    for(x[1] = 0; x[1] < n[1]; x[1]++){
	      int pyg = x[1]+shiftY;
	      if(pyg != j && pyg != ny-1){continue;}
	      for(x[0] = 0; x[0] < n[0]; x[0]++){
		int pxg = x[0]+shiftX;
		if(pxg != i && pxg != nx-1){continue;}

		int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
		Fprintf(fp,"muon 2pt FUNCTION (+-parity proj) tsrc= %d tsink= %d mom= %d %d %d %e %e  %e %e\n", 
			qp_arg->t, tg, pxg, pyg, pzg,
			two_pt[site].real(), 
			two_pt[site].imag(),
			two_pt_np[site].real(), 
			two_pt_np[site].imag());
	      }
	    }
	  }
	}
	Fclose(T,fp);
      }
    }
  }// end mom

  sfree(two_pt);
}

void AlgMuon::two_point(int *mom, QPropWMomSrc &prop, int conj)
{

  char *fname = "two_point(QPropWMomSrc &prop)";
  VRB.Func(cname,fname);
  WilsonMatrix temp;

  // set up for fourier transforms
  //----------------------------------------------------------------
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();

  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();

  const int Lt = GJP.Tnodes()*GJP.TnodeSites();
  const int Lx = GJP.Xnodes()*GJP.XnodeSites();
  const int Ly = GJP.Ynodes()*GJP.YnodeSites();
  const int Lz = GJP.Znodes()*GJP.ZnodeSites();

  int x[4];
  
  // two-point functions
  //----------------------------------------------------------------
  Rcomplex *two_pt = (Rcomplex*) smalloc(GJP.VolNodeSites()*sizeof(Rcomplex));
  Rcomplex *two_pt_np = (Rcomplex*) smalloc(GJP.VolNodeSites()*sizeof(Rcomplex));
  Rcomplex *two_pt_p = (Rcomplex*) smalloc(GJP.VolNodeSites()*sizeof(Rcomplex));

  int dir=3;
  if(mom[0]!=0)dir=0;
  if(mom[1]!=0)dir=1;
  if(mom[2]!=0)dir=2;

  // do the incoming muon two-point function (trace the prop)
  for(x[3] = 0; x[3] < n[3]; x[3]++){
    for(x[2] = 0; x[2] < n[2]; x[2]++){
      for(x[1] = 0; x[1] < n[1]; x[1]++){
	for(x[0] = 0; x[0] < n[0]; x[0]++){
	  
	  int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
          
	  // two pt functions (muon propagator)
	  //project with (1+-gamma_t/2)
	  temp = prop[site];
	  if(conj){
	    temp.gl(-5).gr(-5);
	    temp.hconj();
	  }
	  temp.PParProjectSink();
	  two_pt[site] = temp.Trace() ;
	  temp = prop[site];
	  if(conj){
	    temp = prop[site];
	    temp.gl(-5).gr(-5);
	    temp.hconj();
	  }
	  temp.NParProjectSink();
	  two_pt_np[site] = temp.Trace() ;
          temp = prop[site];
          if(conj){
            temp = prop[site];
            temp.gl(-5).gr(-5);
            temp.hconj();
          }
          temp.gl(dir);
          two_pt_p[site] = temp.Trace() ;
	}
      }
    }
  }
  
  // Fourier transform spatial comps. (2 is size in Floats per site)
  // since we will take anti-muon -> muon for two_pt_np, flip sign in fft, 
  // or 0 last arg. ok since for "incoming" prop (two_pt), we only need 0 mom case
  if(conj){
    FFTXYZ((Float*)two_pt, 1, 2, 1, 1);	
    FFTXYZ((Float*)two_pt_np, 1, 2, 1, 1);
    FFTXYZ((Float*)two_pt_p, 1, 2, 1, 1);
  }else{
    FFTXYZ((Float*)two_pt, 1, 2, 1, 0);	
    FFTXYZ((Float*)two_pt_np, 1, 2, 1, 0);
    FFTXYZ((Float*)two_pt_p, 1, 2, 1, 0);
  }
  // debug WI, need 4d FT	
  //FFT4((Float*)two_pt, 1, 2, 1, 0);	
  //FFT4((Float*)two_pt_np, 1, 2, 1, 0);	
  // Print out results
  //------------------
  FileIoType T=ADD_ID;
  char filename[100];
  sprintf(filename,"%s/prop-m%g.mom%d%d%d.%d.dat",
	  DIRML,
	  alg_muon_arg->line_mass,
	  mom[0],mom[1],mom[2],
	  alg_muon_arg->conf);
  for(x[3] = 0; x[3] < n[3]; x[3]++){
    int tg = x[3]+shiftT;
    for(x[2] = 0; x[2] < n[2]; x[2]++){
      int pzg = x[2]+shiftZ;
      if(pzg != mom[2] && pzg != Lz+mom[2]){continue;}
      for(x[1] = 0; x[1] < n[1]; x[1]++){
	int pyg = x[1]+shiftY;
	if(pyg != mom[1] && pyg != Ly+mom[1]){continue;}
	for(x[0] = 0; x[0] < n[0]; x[0]++){
	  int pxg = x[0]+shiftX;
	  if(pxg != mom[0] && pxg != Lx+mom[0]){continue;}
	  
	  int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
          FILE *fp=Fopen(T, filename, "a");
          Fprintf(fp,"2pt FUNC (1+-g4), g(mom) tsink= %d mom= %d %d %d %le %le  %le %le  %le %le\n",
                  tg, pxg, pyg, pzg,
                  two_pt[site].real(),
                  two_pt[site].imag(),
                  two_pt_np[site].real(),
                  two_pt_np[site].imag(),
                  two_pt_p[site].real(),
                  two_pt_p[site].imag());
	  Fclose(T,fp);
	}
      }
    }
  }

  sfree(two_pt_np);
  sfree(two_pt);
}


void AlgMuon::two_pointv2()
{

  char *fname = "two_pointv2";
  VRB.Func(cname,fname);
  WilsonMatrix temp;

  // set up for fourier transforms
  //----------------------------------------------------------------
  int n[5];
  n[0] = GJP.XnodeSites();
  n[1] = GJP.YnodeSites();
  n[2] = GJP.ZnodeSites();
  n[3] = GJP.TnodeSites();
  n[4] = GJP.SnodeSites();
  const int shiftT = GJP.TnodeCoor()*GJP.TnodeSites();
  const int shiftX = GJP.XnodeCoor()*GJP.XnodeSites();
  const int shiftY = GJP.YnodeCoor()*GJP.YnodeSites();
  const int shiftZ = GJP.ZnodeCoor()*GJP.ZnodeSites();
  const Float PX = (2.0*PI)/(GJP.XnodeSites()*GJP.Xnodes());
  const Float PY = (2.0*PI)/(GJP.YnodeSites()*GJP.Ynodes());
  const Float PZ = (2.0*PI)/(GJP.ZnodeSites()*GJP.Znodes());
  const Float PT = (2.0*PI)/(GJP.TnodeSites()*GJP.Tnodes());

  // get the prop
  //----------------------------------------------------------------
  int mom[3]; 
  mom[0]=mom[1]=mom[2]=0;
  // incoming muon prop, only need it once
  qp_arg->cg.mass = alg_muon_arg->line_mass ;
  qp_arg->t =  0  ;
  QPropWMomSrc in(AlgLattice(), qp_arg, mom, common_arg);

  // two-point function
  //----------------------------------------------------------------
  Rcomplex *two_pt_corr_func_pp = (Rcomplex*) smalloc(GJP.Tnodes()*GJP.TnodeSites()*sizeof(Rcomplex));
  Rcomplex *two_pt_corr_func_np = (Rcomplex*) smalloc(GJP.Tnodes()*GJP.TnodeSites()*sizeof(Rcomplex));
  int t;
  int i,j,k;
  int x[4];
  
  Rcomplex ZERO(0.0,0.0);
  Rcomplex I(0.0,1.0);
  
  // source location of the out-going muon
  qp_arg->t = alg_muon_arg->source_time ; // fix this!

  // Loop over momenta
  for(i = 0; i <= XMOM; i++){
    Float p1 = i * PX;
    for(j = 0; j <= YMOM; j++){
      Float p2 = j * PY;
      for(k = 0; k <= ZMOM; k++){
	Float p3 = k * PZ;
	
	// if(i==0 && j==0 && k==0){continue;}
        
	// outgoing muon with momentum p
	mom[0]=i; mom[1]=j; mom[2]=k;
	// exp{-ip.x} source, but will conjugate this below
	QPropWMomSrc out(AlgLattice(), qp_arg, mom, common_arg);
	
	//zero correlation functions
	for(t=0; t< GJP.Tnodes()*GJP.TnodeSites(); t++){
	  two_pt_corr_func_pp[t]= two_pt_corr_func_np[t]=ZERO;
	}          
	
	// loop over sink position
	for(x[3] = 0; x[3] < n[3]; x[3]++){
	  int tglobal = x[3]+shiftT;
	  for(x[2] = 0; x[2] < n[2]; x[2]++){
	    int zglobal = x[2]+shiftZ;
	    for(x[1] = 0; x[1] < n[1]; x[1]++){
	      int yglobal = x[1]+shiftY;
	      for(x[0] = 0; x[0] < n[0]; x[0]++){
		int xglobal = x[0]+shiftX;
		
		int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
                
		Float px = (xglobal)*p1;
		Float py = (yglobal)*p2;
		Float pz = (zglobal)*p3;
		Float pdotx = (px + py + pz);
		// exp{ -i p . x}
		Rcomplex phase(cos(pdotx), -sin(pdotx));
		
		// two pt functions (muon propagator)
		//project with (1+-gamma_t/2)
		temp = out[site];
		temp.PParProjectSink();
		Rcomplex cc_phase = conj(phase);
		two_pt_corr_func_pp[tglobal] += cc_phase * temp.Trace() ;
		//temp.hconj() ;
		//temp.gr(-5).gl(-5);
		temp = out[site];
		temp.NParProjectSink();
		two_pt_corr_func_np[tglobal] += cc_phase * temp.Trace() ;
	      }
	    }
	  }
	}
	// Global sum
	//----------------------------------------------------------------
	for(t=0; t< GJP.Tnodes()*GJP.TnodeSites(); t++){
	  slice_sum((Float*)&two_pt_corr_func_pp[t], 2, 99);
	  slice_sum((Float*)&two_pt_corr_func_np[t], 2, 99);
	}
	// Print out results
	//----------------------------------------------------------------
	if(common_arg->results != 0){
	  FILE *fp;
	  if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
	    ERR.FileA(cname,fname, (char *)common_arg->results);
	  }
	  for(t=0; t< GJP.Tnodes()*GJP.TnodeSites(); t++){
	    Fprintf(fp,"muon 2pt FUNCTION (+/- parity proj) tsrc= %d tsink= %d mom= %d %d %d %e %e  %e %e\n", 
		    qp_arg->t, t, i, j, k,
		    two_pt_corr_func_pp[t].real(), 
		    two_pt_corr_func_pp[t].imag(),
		    two_pt_corr_func_np[t].real(), 
		    two_pt_corr_func_np[t].imag());
	  }
	  Fclose(fp);
	}
      }
    }
  }// end fourier transform
  
  // the zero momentum results
  for(t=0; t< GJP.Tnodes()*GJP.TnodeSites(); t++){
    two_pt_corr_func_np[t]=ZERO;
    two_pt_corr_func_pp[t]=ZERO;
  }          
  for(x[3] = 0; x[3] < n[3]; x[3]++){
    int tglobal = x[3]+shiftT;
    for(x[2] = 0; x[2] < n[2]; x[2]++){
      for(x[1] = 0; x[1] < n[1]; x[1]++){
	for(x[0] = 0; x[0] < n[0]; x[0]++){
	  
	  int site = x[0] + n[0]*(x[1] +n[1]*(x[2] + n[2]*(x[3])));
          
	  //project with (1+gamma_t)/2
	  temp = in[site];
	  temp.PParProjectSink();
	  two_pt_corr_func_pp[tglobal] += temp.Trace() ;
	  //project with (1-gamma_t)/2
	  temp = in[site];
	  temp.NParProjectSink();
	  two_pt_corr_func_np[tglobal] += temp.Trace() ;
	}
      }
    }
  }
  // Global sum
  //----------------------------------------------------------------
  for(t=0; t< GJP.Tnodes()*GJP.TnodeSites(); t++){
    slice_sum((Float*)&two_pt_corr_func_pp[t], 2, 99);
    slice_sum((Float*)&two_pt_corr_func_np[t], 2, 99);
  }
  // Print out results
  //----------------------------------------------------------------
  if(common_arg->results != 0){
    FILE *fp;
    if( (fp = Fopen((char *)common_arg->results, "a")) == NULL ) {
      ERR.FileA(cname,fname, (char *)common_arg->results);
    }
    for(t=0; t< GJP.Tnodes()*GJP.TnodeSites(); t++){
      Fprintf(fp,"muon ZeroMom 2pt FUNCTION (+/- parity proj) tsrc= %d t= %d  %e %e  %e %e\n", 
	      0, t,
	      two_pt_corr_func_pp[t].real(), 
	      two_pt_corr_func_pp[t].imag(),
	      two_pt_corr_func_np[t].real(), 
	      two_pt_corr_func_np[t].imag());
    }
    Fclose(fp);
  }
  sfree(two_pt_corr_func_np);
  sfree(two_pt_corr_func_pp);
}



char * LatticeContainer::cname = "LatticeContainer";
LatticeContainer::LatticeContainer(){
  size_t mat_size = GJP.VolNodeSites()*4;
  gauge_p = new Matrix[mat_size];
}
LatticeContainer::~LatticeContainer(){
  delete[] gauge_p;
}

void LatticeContainer::Get(Lattice &lat){
  str_ord = lat.StrOrd();
  lat.CopyGaugeField(gauge_p);
}
void LatticeContainer::Set(Lattice &lat){
  if (str_ord != lat.StrOrd())
    ERR.General(cname,"Set()","Storage ordering of LatticeContainer(%d) doesn't agree with lattice ordering(%d)"
, str_ord,lat.StrOrd());
  lat.GaugeField(gauge_p);
}

CPS_END_NAMESPACE
