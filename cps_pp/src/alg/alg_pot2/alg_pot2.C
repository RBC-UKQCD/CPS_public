#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-10 14:05:33 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_pot2/alg_pot2.C,v 1.3 2012-08-10 14:05:33 chulwoo Exp $
//  $Id: alg_pot2.C,v 1.3 2012-08-10 14:05:33 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_pot2.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_pot2/alg_pot2.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// last modification: manke, Feb 1, 2001
// 1. changes Vst --> Vts (for consistency)
// 2. explicitly transfer max_XYZT in pot_arg
// 3. print out potentials V(x,y,z,t) separately for x,y,z
// 4. changed V?? (potential) --> W?? (wilson loop)
//------------------------------------------------------------------
//
// alg_pot.C
//
// AlgPot is derived from Alg and it measures the potential
// into different propagation directions
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <alg/alg_pot2.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>           // Tr() and ReTr()
#include <util/verbose.h> 
#include <util/error.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

using namespace std;
//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgPot2::AlgPot2(Lattice& latt, 
	     CommonArg *c_arg,
	     PotArg *arg) : 
	     Alg(latt, c_arg) 
{
  cname = "AlgPot2";
  char *fname = "AlgPot2(L&,CommonArg*,PotArg*)";
  VRB.Func(cname,fname);

  // Initialize the argument pointer
  //----------------------------------------------------------------
  if(arg == 0)
    ERR.Pointer(cname,fname, "arg");
  alg_pot_arg = arg;

  // Calculate normalization factor for each processor
  //----------------------------------------------------------------

  int total_sites    = GJP.VolNodeSites() * GJP.Xnodes() *
                       GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  int int_norm       = GJP.Colors()*total_sites;

  norm_fac           = 1.0 / ( Float(int_norm) );
  xiB2               = GJP.XiBare()*GJP.XiBare();
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
AlgPot2::~AlgPot2() {
  char *fname = "~AlgPot2()";
  VRB.Func(cname,fname);
}





//------------------------------------------------------------------
// run()
//   check is used to set r2max
//------------------------------------------------------------------
void AlgPot2::run()
{
#if TARGET==cpsMPI
  using MPISCU::fprintf;
#endif
  
  char *fname = "run()";
  VRB.Func(cname,fname);
  
  int sweep    = alg_pot_arg->sweep;
  int prop_dir = alg_pot_arg->prop_dir;
  int check    = alg_pot_arg->check;

  // maximal loop extent in (X,Y,Z,T) direction 
  // running over all the biggest loops would often be very expensive

  int max_ext[4];
 
  int max_X = alg_pot_arg->max_X;
  int max_Y = alg_pot_arg->max_Y;
  int max_Z = alg_pot_arg->max_Z;
  int max_T = alg_pot_arg->max_T;

  max_ext[0] = max_X;
  max_ext[1] = max_Y;
  max_ext[2] = max_Z;
  max_ext[3] = max_T;

  int max_prop = max_ext[prop_dir];

  /*
  // to use pot_arg, r2max is defined here 
  int r2max = 1;
  for (int i = 0; i < 4; i++)
    if( i != prop_dir ) {
      if( r2max <= max_ext[i] ) r2max = max_ext[i];
    }
  r2max = r2max*r2max;
  */
  int r2max = check*check;

  // orthogonal directions to prop_dir 
  int dir_orth[3];
  int idir = 0;
  for (int i = 0; i < 4; i++)
    if( i != prop_dir ) {
      dir_orth[idir] = i;
      idir++;
    }

  int max_1 = max_ext[dir_orth[0]];
  int max_2 = max_ext[dir_orth[1]];
  int max_3 = max_ext[dir_orth[2]];

  int num_loop = 0;
  for (int ext1 = 0; ext1 <= max_1; ext1++){
    for (int ext2 = 0; ext2 <= max_2; ext2++){
      for (int ext3 = 0; ext3 <= max_3; ext3++){  

	int cmax = ext3; 
	if (cmax <= ext2) {cmax = ext2;}
	if (cmax <= ext1) {cmax = ext1;}
	int cmin = ext1;
	if (cmin >= ext2) {cmin = ext2;}
	if (cmin >= ext3) {cmin = ext3;}
	  
	int cmid = ext1+ext2+ext3 - cmax - cmin;
	  
	int r2 = ext1*ext1 + ext2*ext2 + ext3*ext3;
	if( r2 <= r2max || cmid == 0 ) num_loop+=max_prop;
      }
    }
  }
			      
  Float aniso_factor;
  Float pot_real;
  Float pot_imag;

  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  // static QQ-pair propagates in temporal direction
  // calculate potential for all possible distances
  // on-axis and off-axis
  
  // allocate space for potential and multiplicity
  // actually the following is bigger than necessary
  // as the loops below are only UP TO < MAX_R2
  // its safe for now though
  //  Complex *Wts = new Complex[max_r2];
  //  int *number  = new int[max_r2];
  Complex *Wloop = new Complex[num_loop];
  //  int (*qq_sep)[6];
  //  qq_sep = new int[num_loop][6];
  int *qq_sep0 = new int[num_loop];
  int *qq_sep1 = new int[num_loop];
  int *qq_sep2 = new int[num_loop];
  int *qq_sep3 = new int[num_loop];
  int *qq_sep4 = new int[num_loop];
  int *qq_sep5 = new int[num_loop];


  for (int i = 0; i < num_loop; i++){
    Wloop[i] = 0.;
    qq_sep0[i] = -1;  // x
    qq_sep1[i] = -1;  // y
    qq_sep2[i] = -1;  // z
    qq_sep3[i] = -1;  // t
    qq_sep4[i] =  0;  // prop 
    qq_sep5[i] =  0;  // r2
  }

  int iloop = 0;
  for (int ext_prop=1; ext_prop <= max_prop ; ext_prop++){

    // aniso_factor = xi0^{-2* (number of temporal links) }
    aniso_factor = power(xiB2,-1.*ext_prop);
    
    // vary the extent of static Q\bar Q source along the three
    // different spatial directions
    // there is a certain restrcition here in that we only
    // consider spatial paths like (x,x,x,x,y,y,z,z,z,z,z)
    // but not the Oh-related ones
    //  e.g. not: (x,y,z,z,x,x,y,z,z,x,z)
    
    for (int ext1 = 0; ext1 <= max_1; ext1++){
      for (int ext2 = 0; ext2 <= max_2; ext2++){
	for (int ext3 = 0; ext3 <= max_3; ext3++){  

	  int cmax = ext3; 
	  int dir_max = 2;
	  if (cmax <= ext2) {cmax = ext2; dir_max=1;}
	  if (cmax <= ext1) {cmax = ext1; dir_max=0;}
	  int cmin = ext1;
	  int dir_min = 0;
	  if (cmin >= ext2) {cmin = ext2; dir_min=1;}
	  if (cmin >= ext3) {cmin = ext3; dir_min=2;}
	  
	  int cmid = ext1+ext2+ext3 - cmax - cmin;
	  int dir_mid = 0+1+2 - dir_max - dir_min;
	  
	  int r2 = ext1*ext1 + ext2*ext2 + ext3*ext3;
	  //	  if( r2 <= r2max ) {
      	  if( r2 <= r2max || cmid == 0 ) {
	    
	    int check1 = 
	      dir_max*dir_max + dir_mid*dir_mid + dir_min*dir_min;
	    if ( check1 != 5 ) {
	      cout << " sum dir*dir should be 5 not " << check1 << endl;
	      exit(-1);
	    }

	    int cmax2 = 2*cmax;
	    int cmid2 = 2*cmid;
	    int cmin2 = 2*cmin;
	    int chi_mid = cmid2 - cmax;
	    int chi_min = cmin2 - cmax;

	    // squared distance between poitn and origin
	    // int r2 = ext1*ext1 + ext2*ext2 + ext3*ext3;
	    // number of links in Wilson Loop
	    int nlinks = 2*(ext1+ext2+ext3+ext_prop);
	    
	    // allocate space for path
	    int *path = new int[nlinks]; 
	    int ipath = 0;

	    // path between Q and Qbar
	    for (int i = 1; i <= cmax; i++) {
	      path[ipath] = dir_orth[dir_max];
	      ipath++;
	      if (chi_mid >= 0) {
		chi_mid -= cmax2;
		path[ipath] = dir_orth[dir_mid];
		ipath++;
	      }
	      if (chi_min >= 0) {
		chi_min -= cmax2;
		path[ipath] = dir_orth[dir_min];
		ipath++;
	      }
	      chi_mid += cmid2;
	      chi_min += cmin2;
	    }
	    
	    // path for propagation
	    for (int i = 0; i < ext_prop; i++) {
	      path[ipath] = prop_dir;
	      ipath++;
	    }
	    
	    // x = 0 , -x = 4 = x + 4
	    // y = 1 , -y = 5 = y + 4
	    // z = 2 , -z = 6 = z + 4
	    // t = 3 , -t = 7 = t + 4
	    // inverse path between Q and Qbar
	    int qq_path = ext1 + ext2 + ext3;
	    int inv_path = qq_path - 1;
	    for (int i = 0; i < qq_path; i++) {
	      path[ipath] = path[inv_path] + 4;  // inverse direction
	      ipath++;
	      inv_path--;
	    }
	    
	    // inverse path for propagation
	    for (int i = 0; i < ext_prop; i++) {
	      path[ipath] = prop_dir + 4;
	      ipath++;
	    }
	    
	    // check for path length
	    if (ipath != nlinks) {
	      cout << "illegal path !!!" << endl;
	      cout << "ipath  = " << ipath << endl;
	      cout << "nlinks = " << nlinks << endl;
	      exit(-1);
	    }
	    
	    // calculated and accumulate the LOOP
	    // over all possible origins on the lattice
	    // printf("Looping over all sources \n");
	    
	    int x[4]={0,0,0,0};
	    Matrix LOOP;
	    //	  LOOP.ZeroMatrix();
	    Complex pot=0.;
	    
	    int isite = 0;
	    for(x[0]=0;x[0]<GJP.XnodeSites();x[0]++)
	      for(x[1]=0;x[1]<GJP.YnodeSites();x[1]++)
		for(x[2]=0;x[2]<GJP.ZnodeSites();x[2]++)
		  for(x[3]=0;x[3]<GJP.TnodeSites();x[3]++){
		    
		    LOOP.ZeroMatrix();
		    lat.PathOrdProdPlus(LOOP, x, path, nlinks);
		    pot += LOOP.Tr();
		    isite++;
		    
		  } // end loop over all local origins
	    // accumulate result for each distance and its multiplicity.
	    // Notice that some distances can not be
	    // represented on the 3D spatial lattice ==
	    // some integers can not be represented as
	    // a sum of 3 squares, in this case number[r2] remains zero
	    Wloop[iloop] = pot;
	    int tmp_sep[4] = {-1,-1,-1,-1};
	    tmp_sep[dir_orth[0]] = ext1;
	    tmp_sep[dir_orth[1]] = ext2;
	    tmp_sep[dir_orth[2]] = ext3;
	    qq_sep0[iloop] = tmp_sep[0];
	    qq_sep1[iloop] = tmp_sep[1];
	    qq_sep2[iloop] = tmp_sep[2];
	    qq_sep3[iloop] = tmp_sep[3];
	    qq_sep4[iloop] = ext_prop;
	    qq_sep5[iloop] = r2;
	    iloop++;
	    
	    // deallocate space for path
	    delete [] path;
	    //	  cout << "after delete path" << endl;
	  } // end for if(r2<=r2max||cmid==0)
	} // end for ext3
      } // end for ext2
    } // end for ext1
    
  } // end for ext_prop=0; ext_prop <= max_T
  
  if (iloop != num_loop ) {
    cout << "something is wrong about loop counting" << endl;
    cout << "iloop = " << iloop << endl;
    cout << "num_loop = " << num_loop << endl;
    exit(-1);
  }
  
  // dump the averaged potential for all r2 at this timeslice
  char wfname[100];
  sprintf(wfname, "Wloop%04dsme%03d.dat", sweep/1000, sweep%1000);
  FILE *fp= Fopen(wfname, "w");
  //  FILE *fp= Fopen("Wloop_new.dat", "w");
  
  Fprintf(fp, "# max_X= %d max_Y= %d max_Z= %d max_T= %d \n",
	  max_X, max_Y, max_Z, max_T);
  Fprintf(fp, "# sweep= %d  prop_dir= %d  num_loop= %d\n",
	  sweep/1000, prop_dir, num_loop);
  Fprintf(fp, "# smear= %d \n", sweep%1000);
  Fprintf(fp, "# x  y  z  t prop   pot_real       pot_imag\n");
  
  for (int i = 0; i < num_loop; i++){
    Complex W =  Wloop[i];
    pot_real = real(W);
    pot_imag = imag(W);
    glb_sum(&pot_real) ;
    glb_sum(&pot_imag) ;
    
    // normalise wrt colour and total space-time
    // take into account the anisotropy factor for each temporal link
    pot_real *= aniso_factor*norm_fac;
    pot_imag *= aniso_factor*norm_fac;

    Fprintf(fp," %d  %d  %d  %d  %d    %e   %e  \n",
	    //	    qq_sep[i][0],qq_sep[i][1],qq_sep[i][2],qq_sep[i][3],
	    //	    qq_sep[i][4], pot_real, pot_imag);
	    qq_sep0[i],qq_sep1[i],qq_sep2[i],qq_sep3[i],
	    qq_sep4[i], pot_real, pot_imag);
    //    Fprintf(fp," %d  %d  %d  %d  %d  %d   %e   %e  \n",
    //	    qq_sep[i][0],qq_sep[i][1],qq_sep[i][2],qq_sep[i][3],
    //	    qq_sep[i][4], qq_sep[i][5], pot_real, pot_imag);
  }
  
  Fclose(fp);
  //  cout << "total_site = " << (1.0 / norm_fac)/3.0 << endl;

  delete[] Wloop;
  //  delete[] qq_sep;
  delete[] qq_sep0;
  delete[] qq_sep1;
  delete[] qq_sep2;
  delete[] qq_sep3;
  delete[] qq_sep4;
  delete[] qq_sep5;

}
CPS_END_NAMESPACE
