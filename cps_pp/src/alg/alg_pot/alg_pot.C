#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:39 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_pot/alg_pot.C,v 1.8 2004-08-18 11:57:39 zs Exp $
//  $Id: alg_pot.C,v 1.8 2004-08-18 11:57:39 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: alg_pot.C,v $
//  $Revision: 1.8 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_pot/alg_pot.C,v $
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
#include <alg/alg_pot.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/smalloc.h>
#include <util/vector.h>           // Tr() and ReTr()
#include <util/verbose.h> 
#include <util/error.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

// #define max(A, B) ((A) > (B) ? (A) : (B))
// #define min(A, B) ((A) < (B) ? (A) : (B))

//------------------------------------------------------------------
// Constructor 
//------------------------------------------------------------------
AlgPot::AlgPot(Lattice& latt, 
	     CommonArg *c_arg,
	     PotArg *arg) : 
	     Alg(latt, c_arg) 
{
  cname = "AlgPot";
  char *fname = "AlgPot(L&,CommonArg*,PotArg*)";
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
AlgPot::~AlgPot() {
  char *fname = "~AlgPot()";
  VRB.Func(cname,fname);
}





//------------------------------------------------------------------
// run()
//
// this calculates 
// 1. the regular potential from the Wilson Loop in temporal direction
//    and with all squared spatial distances r2=x*x + y*y + z*z 
//    --> Wts
// 2. the sideways potential from wilson Loop in z-direction (=2)
//    and with QQ separaton on-axis along
//    a) x  --> Wzx
//    b) t  --> Wzt
//------------------------------------------------------------------
void AlgPot::run()
{
#if TARGET==cpsMPI
    using MPISCU::fprintf;
#endif

  char *fname = "run()";
  VRB.Func(cname,fname);

  int sweep    = alg_pot_arg->sweep;
  int prop_dir = alg_pot_arg->prop_dir;
  // int check    = alg_pot_arg->check;

  // maximal loop extent in (X,Y,Z,T) direction 
  // running over all the biggest loops would often be very expensive
  int max_X    = alg_pot_arg->max_X;
  int max_Y    = alg_pot_arg->max_Y;
  int max_Z    = alg_pot_arg->max_Z;
  int max_T    = alg_pot_arg->max_T;
  // maximal squared distance from origin (0,0,0)
  int max_r2 = (max_X*max_X + max_Y*max_Y + max_Z*max_Z);

  // printf(" max_XYZT = %d %d %d %d max_r2 = %d \n",max_X,max_Y,max_Z,max_T,max_r2);
  
  Float aniso_factor;
  Float pot_real;
  Float pot_imag;



  // Set the Lattice pointer
  //----------------------------------------------------------------
  Lattice& lat = AlgLattice();

  printf("prop_dir = %d sweep = %d \n",prop_dir,sweep);
  // printf("norm_fac = %e \n", norm_fac);
  // printf("XI0 = %e XI02 = %e \n", GJP.XiBare(),xiB2);

  // if propagation direction is NOT t=3 or z=2 --> exit
  if ((prop_dir != 3) && (prop_dir !=2) ){
    printf("Lattice: prop_dir= %d -- not yet implemented\n",prop_dir);
    //    exit(1);
  }
  
  if (prop_dir == 3) {
    // static QQ-pair propagates in temporal direction
    // calculate potential for all possible distances
    // on-axis and off-axis

   // printf(" calculate regular potential in temporal direction\n");

    // allocate space for potential and multiplicity
    // actually the following is bigger than necessary
    // as the loops below are only UP TO < MAX_R2
    // its safe for now though
    Complex *Wts = new Complex[max_r2];
    int *number  = new int[max_r2];

    for (int ext_prop=1; ext_prop < max_T ; ext_prop++){
      //printf("ext_prop = %d max_T = %d \n",ext_prop,max_T);

      // aniso_factor = xi0^{-2* (number of temporal links) }
      aniso_factor = power(xiB2,-1.*ext_prop);
      // printf("XiBare = %e -> temporal links needs factor %e \n",GJP.XiBare(),aniso_factor);

      // initialisation of potential and the number of points
      // with distance r2 (new for every different ext_prop)
      int r2;
      for (r2=0; r2 < max_r2; r2++){
	Wts[r2] = 0.;
	number[r2] = 0;
      }
    
      // vary the extent of static Q\bar Q source along the three
      // different spatial directions
      // there is a certain restrcition here in that we only
      // consider spatial paths like (x,x,x,x,y,y,z,z,z,z,z)
      // but not the Oh-related ones
      //  e.g. not: (x,y,z,z,x,x,y,z,z,x,z)

      for (int ext1 = 0; ext1 < max_X; ext1++){
	for (int ext2 = 0; ext2 < max_Y; ext2++){
	  for (int ext3 = 0; ext3 < max_Z; ext3++){  
	    // printf("ext1-3 = %d %d %d  ",ext1,ext2,ext3);
      
	    int offset=0;
	    // squared distance between poitn and origin
	    r2 = ext1*ext1 + ext2*ext2 + ext3*ext3;
	    // number of links in Wilson Loop
	    int nlinks = 2*(ext1+ext2+ext3+ext_prop);
	    // printf("r2 = %d nlinks = %d ",r2,nlinks);

	    // allocate space for path
	    int *path = new int[nlinks]; 
	    int i;
	    // set up path
	    for (i=0; i < ext1; i++) path[i] = 0;                 // x = 0
	    offset = ext1;
	    for (i=offset; i < offset+ext2; i++) path[i] = 1;     // y = 1
	    offset = offset + ext2;
	    for (i=offset; i < offset+ext3; i++) path[i] = 2;     // z = 2
	    offset = offset + ext3;
	    for (i=offset; i < offset+ext_prop; i++) path[i] = 3; // t = 3
	    offset = offset + ext_prop;

	    // now the whole way back
	    for (i=offset; i < offset+ext3; i++) path[i] = 6;     // -z = 6
	    offset = offset + ext3;
	    for (i=offset; i < offset+ext2; i++) path[i] = 5;     // -y = 5
	    offset = offset + ext2;
	    for (i=offset; i < offset+ext1; i++) path[i] = 4;     // -x = 4
	    offset = offset + ext1;
	    for (i=offset; i < offset+ext_prop; i++) path[i] = 7; // -t=7

	    /*
	    printf("path = ");
	    for (i=0;i<nlinks;i++) printf("%d ",path[i]);
	    printf("\n ");
	    */

	    // calculated and accumulate the LOOP
	    // over all possible origins on the lattice
	    // printf("Looping over all sources \n");

	    int x[4]={0,0,0,0};
	    Matrix LOOP;
	    LOOP.ZeroMatrix();
	    Complex pot=0.;

	    for(x[0]=0;x[0]<GJP.XnodeSites();x[0]++)
	      for(x[1]=0;x[1]<GJP.YnodeSites();x[1]++)
		for(x[2]=0;x[2]<GJP.ZnodeSites();x[2]++)
		  for(x[3]=0;x[3]<GJP.TnodeSites();x[3]++){
		    //  lat.PathOrdProdPlus(LOOP, x, path, nlinks);
		    	   
		    lat.PathOrdProd(LOOP, x, path, nlinks);
		    pot += LOOP.Tr();
		   
		  } // end loop over all local origins

	    // for PathOrdProdPlus calculate trace of accumulated LOOP
	    // pot = LOOP.Tr();
 
	    // printf("Trace = %e %e  \n", real(pot),imag(pot));

	    
	    // accumulate result for each distance and its multiplicity.
	    // Notice that some distances can not be
	    // represented on the 3D spatial lattice ==
	    // some integers can not be represented as
	    // a sum of 3 squares, in this case number[r2] remains zero
	    number[r2] += 1;
	    Wts[r2]    += pot;

	    // printf("r2 = %d, number = %d, pot = %e %e Wts = %e \n",
	    //	   r2, number[r2], real(pot),imag(pot),real(Wts[r2]));
	    

	    // write this potential to file
	    // i.e before summing over all other pot with the same
	    // value of r2 
	    // e.g.  r=5 is either (5,0,0) or (4,3,0)
	    // treat those seperately here

	    char *filename="Wt_x_y_z"; // example name
	    sprintf(filename,"Wt_%d_%d_%d",ext1,ext2,ext3);
	    FILE *fp= Fopen(filename, "a");
	    Float pot_tmp_real = real(pot);
	    Float pot_tmp_imag = imag(pot);
	    glb_sum(&pot_tmp_real) ;
	    glb_sum(&pot_tmp_imag) ;
	    pot_tmp_real *= aniso_factor*norm_fac;
	    pot_tmp_imag *= aniso_factor*norm_fac;
	    Fprintf(fp,"%d %d %d %e %e \n",
		    sweep, ext_prop, r2, pot_tmp_real,pot_tmp_imag);
	    Fclose(fp);

	    // deallocate space for path
	    delete [] path;
	  } // end for ext3
	} // end for ext2
      } // end for ext1
      
      // dump the averaged potential for all r2 at this timeslice
      FILE *fp= Fopen("Wts", "a");
      
      for (r2=0; r2 < max_r2; r2++){
	if (number[r2]!=0) {
	  Complex W =  Wts[r2]/number[r2];
	  pot_real = real(W);
	  pot_imag = imag(W);
	  // printf("before global sum %d %d %d %e %e \n",
	  //	 sweep, ext_prop, r2, pot_real,pot_imag);
	  glb_sum(&pot_real) ;
	  glb_sum(&pot_imag) ;

	  // normalise wrt colour and total space-time
	  // take into account the anisotropy factor for each temporal link
	  pot_real *= aniso_factor*norm_fac;
	  pot_imag *= aniso_factor*norm_fac;


	  //printf("after global sum %d %d %d %e %e \n",
	  //sweep, ext_prop, r2, pot_real,pot_imag);
	  
	  Fprintf(fp,"%d %d %d %e %e \n",
		  sweep, ext_prop, r2, pot_real,pot_imag);

	} // end if (number[r2] !=0)
      } // end for r2=0; r2 < max_r2
    
      Fclose(fp);

    } // end for ext_prop=0; ext_prop < max_T
 
    // deallocate space for potential and multiplicity
    delete[] Wts;
    delete[] number;
  } // endif (prop_dir == 3 := TIME)
  

  // ======================================================

  if (prop_dir == 2) {

    // calculate "side-ways" potential where the 
    // Q-Q pair propagates only into the spatial z-direction
    // in this case we calculate only the TWO ON-AXIS Potential
    // with the quark separation into the X- or T-direction

    //printf(" calculate sideways potential in direction %d \n",prop_dir);

    // loop over all possible "time" extents
    
    for (int ext_prop=1; ext_prop < max_Z ; ext_prop++){
      //printf("ext_prop = %d of max_Z = %d\n",ext_prop,max_Z);

      // loop over the two different choices for the spatial
      // separation, i.e. X and T
      for (int choice=0; choice < 2 ; choice++){
	int sep_dir; // direction of Q-Q separation
	int max_ext; // maximal extent of Q-Q separation
	char *filename; // filename for potential data
		
	//printf("choice = %d\n",choice);
	switch (choice) {
	case 0 :  sep_dir = 0; max_ext = max_X; filename = "Wzx"; break; 
	case 1 :  sep_dir = 3; max_ext = max_T; filename = "Wzt"; break;
	default :  sep_dir = 3; max_ext = max_T; filename = "Wzt"; break;
	}
	//printf(" calculate quark separation in %d -> %s \n",sep_dir,filename);

	// this is now simpler than above as we have only
	// on-axis separation and there is a unique squared distance
	// for each extent "ext" in this particular direction (x or t)

	for (int ext=1; ext < max_ext ; ext++){	

	  int i;
	  int offset;
	  // squared distance (on-axis)
	  // int r2 = ext*ext;
	  // number of links in loop based on this separation
	  int nlinks = 2*(ext+ext_prop);
	  
	  // allocate space for path
	  int *path = new int[nlinks]; 
	  
	  // set up path
	  // QQ-separation (
	  for (i=0; i < ext; i++) path[i] = sep_dir;
	  offset = ext;
	  // QQ-propagation
	  for (i=offset; i < offset+ext_prop; i++) path[i] = prop_dir;
	  offset = offset + ext_prop;

	  // now the whole way back
	  // -sep_dir (-x=4, -t=7 ) --> sep_dir + 4
	  for (i=offset; i < offset+ext; i++) path[i] = sep_dir+4;
	  offset = offset + ext;

	  // now go backward in "time"=-z=6
	  for (i=offset; i < offset+ext_prop; i++) path[i] = prop_dir+4;

	  /*
	  printf("path = ");
	  for (int j=0; j<nlinks; j++) printf("%d",path[j]);
	  printf("\n");
	  */

	  // calculated and accumulate the LOOP
	  // over all possible origins on the lattice
	  int x[4];
	  Matrix LOOP;
	  Complex pot = 0.;

	  for(x[0]=0;x[0]<GJP.XnodeSites();x[0]++)
	    for(x[1]=0;x[1]<GJP.YnodeSites();x[1]++)
	      for(x[2]=0;x[2]<GJP.ZnodeSites();x[2]++)
		for(x[3]=0;x[3]<GJP.TnodeSites();x[3]++){

		  lat.PathOrdProd(LOOP, x, path, nlinks);
		  pot += LOOP.Tr();

		}

	  
	  
	  //printf("write to file = %s \n",filename);
	  //printf("%d %d %d %e %e \n",
	  // sweep, ext_prop, r2, real(pot),imag(pot));

	  pot_real = real(pot);
	  pot_imag = imag(pot);

	  
	  // printf("before global sum %d %d %d %e %e \n",
	  //	 sweep, ext_prop, r2, pot_real,pot_imag);

	  // global sum over all processors
	  glb_sum(&pot_real) ;
	  glb_sum(&pot_imag) ;

	  
	  // normalise wrt colour and total space-time 
	  // for anisotropic lattices this is dependent on 
	  // whether this loop extented into the temporal direction

	  if (sep_dir==3){

	    // each temporal link should be divided by xi0
	    // aniso_factor = xi0^{-2*temp_ext}
	    aniso_factor = power(xiB2,-1.*ext);

	    //printf("ext = %d, xiB2 = %e, aniso_factor = %e pot before = %e\n",
	    //   ext,xiB2,aniso_factor,pot*norm_fac);

	    pot_real *= aniso_factor*norm_fac;
	    pot_imag *= aniso_factor*norm_fac;

	    // if (ext_prop==1) Wzt[ext] = pot_real; // for external use
	    
	    /*
	    // write W(z,t) to file
	    char *fn="Wz99t99";
	    sprintf(fn,"Wz%dt%d",ext_prop,ext);
	    FILE *fp= Fopen(fn, "a");
	    Fprintf(fp,"%d %e %e \n",sweep, pot_real,pot_imag);
	    Fclose(fp);
	    */
	  } else {

	    // if seperation is in spatial direction
	    // no extra care is necessary --> loop is completely spatial

	    pot_real *= norm_fac;
	    pot_imag *= norm_fac;

	    // if (ext_prop==1) Wzx[ext] = pot_real; // for external use

	    /*
	    // write W(z,x) to file
	    char *fn = "Wz99x99";
	    sprintf(fn,"Wz%dx%d",ext_prop,ext);
	    FILE *fp= Fopen(fn, "a");
	    Fprintf(fp,"%d %e %e \n",sweep, pot_real,pot_imag);
	    Fclose(fp);
	    */
	  }

	  // printf("after global sum %d %d %d %e %e \n",
	  //	 sweep, ext_prop, r2, pot_real,pot_imag);
	  
	  FILE *fp= Fopen(filename, "a");
	  Fprintf(fp,"%d %d %d %e %e \n",
		  sweep, ext_prop, ext, pot_real,pot_imag);
	  Fclose(fp);

	  // deallocate space for path
	  delete [] path;

	} // end for ext=0; ext< max_ext;
      } // end for choice=0; choice < 2 ; 
    } // end for ext_prop =0; ext_prop < max_Z ;
  } // endif (prop_dir == 2 := Z-direction)
}

CPS_END_NAMESPACE
