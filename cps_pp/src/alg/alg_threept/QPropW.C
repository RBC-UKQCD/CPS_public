#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:45 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_threept/QPropW.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Id: QPropW.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.8  2002/12/04 17:16:27  zs
//  Merged the new 2^4 RNG into the code.
//  This new RNG is implemented in the LatRanGen class.
//  The following algorithm and utility classes are affected:
//
//  AlgEig                  Fdwf
//  AlgGheatBath            Fstag
//  AlgHmd                  GlobalJobParameter
//  AlgNoise                Lattice
//  AlgPbp                  Matrix
//  AlgThreept              RandomGenerator
//                          Vector
//
//  Revision 1.7  2001/08/16 10:49:41  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.5  2001/07/03 17:00:46  anj
//
//  Multiple minor alterations to change some #include's from referring to
//  files relative to the top-level source directory to referring to files
//  relative to the source-file positions.  This alteration makes the code
//  backwards compatable with the make structure of QCDSP, although this
//  may have to be changed to a more usual form in the future. Anj.
//
//  Revision 1.4  2001/06/28 14:34:10  anj
//
//  The core ANSIfication should now be complete.  There are a few
//  remaining issues, but this version should compile anywhere and be
//  backward compatable with QCDSP (although this requires the top source
//  directory (.../phys/ to be added to the include path).
//
//  The serial GCC version has also been tested, and all test programs
//  appear to behave as they should (not to imply that they all work, but
//  I believe those that should work are ok).  There are minor differences
//  in the results due to rounding, (see example pbp_gccsun.dat files),
//  but that is all.
//
//  Anj.
//
//  Revision 1.3  2001/06/21 15:40:10  anj
//  Updated the _TARTAN ifdefs, using PARALLEL instead (where appropriate).Anj
//
//  Revision 1.2  2001/06/19 18:11:31  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:00  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: QPropW.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_threept/QPropW.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
//
// The class functions for QpropW.
//
// For now this is specific to three colors. The constructor
// will exit if the number of colors is not equal to three.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>     // exit()
#include <stdio.h>
#include<config.h>
#include<alg/qpropw.h>
#include<alg/common_arg.h>
CPS_START_NAMESPACE
#ifdef PARALLEL
CPS_END_NAMESPACE
#include<comms/glb.h>
#include <sysfunc.h>
CPS_START_NAMESPACE
#endif


CPS_END_NAMESPACE
#include<util/dirac_op.h>
#include<util/dwf.h>
#include<util/random.h>
CPS_START_NAMESPACE

QPropW::QPropW(void)
{};
 
QPropWWallSrc::QPropWWallSrc(void)
{};

// "equal" operator for QPropW
//QPropW& QPropW::operator=(const QPropW& rhs)
//{
//  VRB.Sfree("QPropW","operator=", "prop", prop);
//   sfree(prop);
//printf("operator= mallocing prop\n");
//   prop = (WilsonMatrix*) smalloc(GJP.VolNodeSites() * sizeof(WilsonMatrix));
//   if (prop == 0) ERR.Pointer("QPropW","operator=", "prop");
//   VRB.Smalloc("QPropW","operator=", "prop", prop, GJP.VolNodeSites() * sizeof(WilsonMatrix));
//   for(int i=0;i<GJP.VolNodeSites(); i++){
//	prop[i]=rhs.prop[i];
//  }
//   return *this;
//}
// "equal" operator for QPropW
//QPropW& QPropW::operator=(QPropW& rhs)
//{
//        return rhs;
//}
 
 
// Generate the prop over the whole lattice
void QPropW::CG(Lattice &lat, CgArg *arg, FermionVectorTp& source, 
	FermionVectorTp& sol , int& iter, Float& true_res)
{
  char *fname = "CG(L&, CgArg*, source&, sol&, int&, Float&)";
  char *cname = "QPropW";
  VRB.Func(cname, fname);

  // Set the node size of the full (non-checkerboarded) fermion field
  //----------------------------------------------------------------
  int ls = GJP.SnodeSites();
  int ls_glb = GJP.SnodeSites()*GJP.Snodes();
  int f_size = GJP.VolNodeSites() * lat.FsiteSize()/GJP.SnodeSites();
  int f_size_5d = f_size * ls;


  // Do inversion
  //----------------------------------------------------------------
  if(lat.Fclass() == F_CLASS_DWF){
    Vector *src_4d = (Vector *)source.data();
    Vector *sol_4d = (Vector *)sol.data();
    Vector *src_5d = (Vector *)smalloc(f_size_5d * sizeof(IFloat));
    if(src_5d == 0)
      ERR.Pointer(cname,fname, "src_5d");
    VRB.Smalloc(cname,fname, "src_5d", src_5d, f_size_5d * sizeof(IFloat));
    Vector *sol_5d = (Vector *) smalloc(f_size_5d * sizeof(IFloat));
    if(sol_5d == 0)
      ERR.Pointer(cname,fname, "sol_5d");
    VRB.Smalloc(cname,fname, "sol_5d", sol_5d, f_size_5d * sizeof(IFloat));
    lat.Ffour2five(src_5d, src_4d, 0, ls_glb-1);
    lat.Ffour2five(sol_5d, sol_4d, ls_glb-1, 0);

    iter = lat.FmatInv(sol_5d, src_5d, arg, &true_res, CNV_FRM_YES, PRESERVE_NO);

    lat.Ffive2four(sol_4d, sol_5d, ls_glb-1, 0);


    VRB.Sfree(cname,fname, "sol_5d", sol_5d);
    sfree(sol_5d);
    VRB.Sfree(cname,fname, "src_5d", src_5d);
    sfree(src_5d);

  }
  else {
    iter = lat.FmatInv((Vector*)sol.data(),(Vector*)source.data(), arg, CNV_FRM_YES, PRESERVE_NO);
  }

}

// Generate the prop over the whole lattice
// Use previous sol as initial guess
void QPropW::CGDwf(Lattice &lat, CgArg *arg, 
	FermionVectorTp& source, Vector* sol_5d)
{
  char *fname = "CGDwf(L&, CgArg*, source&, sol_5d&)";
  char *cname = "QPropW";
  VRB.Func(cname, fname);

  // Set the node size of the full (non-checkerboarded) fermion field
  //----------------------------------------------------------------
  int ls = GJP.SnodeSites();
  int ls_glb = GJP.SnodeSites()*GJP.Snodes();
  int f_size = GJP.VolNodeSites() * lat.FsiteSize()/GJP.SnodeSites();
  int f_size_5d = f_size * ls;


  // Do inversion
  //----------------------------------------------------------------
    Vector *src_4d = (Vector *)source.data();
    Vector *src_5d = (Vector *)smalloc(f_size_5d * sizeof(IFloat));
    if(src_5d == 0)
      ERR.Pointer(cname,fname, "src_5d");
    VRB.Smalloc(cname,fname, "src_5d", src_5d, f_size_5d * sizeof(IFloat));
    lat.Ffour2five(src_5d, src_4d, 0, ls_glb-1);

    int iter = lat.FmatInv(sol_5d, src_5d, arg, CNV_FRM_YES, PRESERVE_NO);


    VRB.Sfree(cname,fname, "src_5d", src_5d);
    sfree(src_5d);

}

// Generate a prop from a wall source
QPropWWallSrc::QPropWWallSrc(Lattice& lat, CgArg *arg, int source_time,
	CommonArg* common_arg)
{

   char *fname = "QPropWWallSrc(L&, CgArg*, source_time)";
   char *cname = "QPropWWallSrc";
   VRB.Func(cname, fname);

   // Set the node size of the full (non-checkerboarded) fermion field
   //----------------------------------------------------------------
   int f_size = GJP.VolNodeSites() * lat.FsiteSize()/GJP.SnodeSites();
   int iter;
   Float true_res;

   // allocate space for the quark propagator
   //----------------------------------------------------------------
   prop = (WilsonMatrix*) smalloc(GJP.VolNodeSites() * sizeof(WilsonMatrix));
   if (prop == 0) ERR.Pointer(cname, fname, "prop");
   VRB.Smalloc(cname, fname, "prop", prop, 12 * f_size * sizeof(Float));
   FermionVectorTp src;
   FermionVectorTp sol;

   for (int spin = 0; spin < 4; spin++)
     for (int color = 0; color < GJP.Colors(); color++){

        // initial guess
        sol.setVolSource(color, spin);
        // set the source
        //src.setWallSource(color, spin, source_time);
	// 3=t dir
        src.setGFWallSource(lat, color, spin, 3, source_time);
        //src.setLandauWallSource(lat, color, spin, source_time);
        // Get the prop
        CG(lat, arg, src, sol, iter, true_res);
	sol.GaugeFixSink(lat, 3);
	//sol.LandauGaugeFixSink(lat, 3);
        // Collect solutions in propagator.
        int j;
        for(int i=0; i<f_size; i+=SPINOR_SIZE){
           j=i/SPINOR_SIZE; // lattice site
           prop[j].load_vec(spin, color, (wilson_vector &)sol[i]);
        }

    	if(common_arg->results != 0){
      	  FILE *fp;
      	  if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
       	     ERR.FileA(cname,fname, (char *)common_arg->results);
      	  }
      	  fprintf(fp, "Cg iters = %d true residual = %e\n",
			iter, true_res);
      	  fclose(fp);
    	}


   } // End spin-color loop


}

// Generate a prop from all 4d points 
QPropWVolSrc::QPropWVolSrc(Lattice& lat, CgArg *arg)
{
   char *fname = "QPropWVolSrc(L&, CgArg*)";
   char *cname = "QPropWVolSrc";
   VRB.Func(cname, fname);

   // Set the node size of the full (non-checkerboarded) fermion field
   //----------------------------------------------------------------
   int f_size = GJP.VolNodeSites() * lat.FsiteSize()/GJP.SnodeSites();
   int iter;
   Float true_res;


   // allocate space for the quark propagator
   //----------------------------------------------------------------
   prop = (WilsonMatrix*) smalloc(GJP.VolNodeSites() * sizeof(WilsonMatrix));
   if (prop == 0) ERR.Pointer(cname, fname, "prop");
   VRB.Smalloc(cname, fname, "prop", prop, 12 * f_size * sizeof(Float));

   FermionVectorTp src;
   FermionVectorTp sol;

   for (int color = 0; color < 1; color++)
     for (int spin = 0; spin < 1; spin++) {

        // initial guess
        sol.setVolSource();
        // set the source
        src.setVolSource(color, spin);
        // Get the prop
        CG(lat, arg, src, sol, iter, true_res);
        // Collect solutions in propagator.
        int j;
        for(int i=0; i<f_size; i+=SPINOR_SIZE){
           j=i/SPINOR_SIZE; // lattice site
           prop[j].load_vec(spin, color, (wilson_vector &)sol[i]);
        }

   } // End spin-color loop

}

// Generate a prop from two QPropWWallSrc's
QPropWWallSrc::QPropWWallSrc(Lattice& lat, CgArg *arg, 
		QPropWWallSrc& prop1, QPropWWallSrc& prop2 )
{
   char *fname = "QPropWWallSrc(L&, CgArg*, source_time)";
   char *cname = "QPropWWallSrc";
   VRB.Func(cname, fname);

   // Set the node size of the full (non-checkerboarded) fermion field
   //----------------------------------------------------------------
   int f_size = GJP.VolNodeSites() * lat.FsiteSize()/GJP.SnodeSites();


   // allocate space for the quark propagator
   //----------------------------------------------------------------
   prop = (WilsonMatrix*) smalloc(GJP.VolNodeSites() * sizeof(WilsonMatrix));
   if (prop == 0) ERR.Pointer(cname, fname, "prop");
   VRB.Smalloc(cname, fname, "prop", prop, 12 * f_size * sizeof(Float));
 
   for(int i=0; i< GJP.VolNodeSites(); i++){
       prop[i]=((Float)0.5)*(prop1.prop[i]+prop2.prop[i]);
   }


}


// Generate a prop from two QPropWWallSrc's
QPropWWallSrc::QPropWWallSrc(Lattice& lat, CgArg *arg,
                QPropWWallSrc* prop1, QPropWWallSrc* prop2 )
{
   char *fname = "QPropWWallSrc(L&, CgArg*, source_time)";
   char *cname = "QPropWWallSrc";
   VRB.Func(cname, fname);
 
   // Set the node size of the full (non-checkerboarded) fermion field
   //----------------------------------------------------------------
   int f_size = GJP.VolNodeSites() * lat.FsiteSize()/GJP.SnodeSites();
 
   // allocate space for the quark propagator
   //----------------------------------------------------------------
   prop = (WilsonMatrix*) smalloc(GJP.VolNodeSites() * sizeof(WilsonMatrix));
   if (prop == 0) ERR.Pointer(cname, fname, "prop");
   VRB.Smalloc(cname, fname, "prop", prop, 12 * f_size * sizeof(Float));
 
   for(int i=0; i< GJP.VolNodeSites(); i++){
       prop[i]=((Float)0.5)*((*prop1).prop[i]+(*prop2).prop[i]);
   }
 
 
}

// Generate a prop from two QPropWRandWallSrc's
QPropWRandWallSrc::QPropWRandWallSrc(Lattice& lat, CgArg *arg, 
		QPropWRandWallSrc* prop1, QPropWRandWallSrc* prop2 )
{
   char *fname = "QPropWRandWallSrc(L&, CgArg*, source_time)";
   char *cname = "QPropWRandWallSrc";
   VRB.Func(cname, fname);

   // Set the node size of the full (non-checkerboarded) fermion field
   //----------------------------------------------------------------
   int f_size = GJP.VolNodeSites() * lat.FsiteSize()/GJP.SnodeSites();


   // allocate space for the quark propagator
   //----------------------------------------------------------------
   prop = (WilsonMatrix*) smalloc(GJP.VolNodeSites() * sizeof(WilsonMatrix));
   if (prop == 0) ERR.Pointer(cname, fname, "prop");
   VRB.Smalloc(cname, fname, "prop", prop, 12 * f_size * sizeof(Float));
 
   int i;
   for(i=0; i< GJP.VolNodeSites(); i++){
       prop[i]=((Float)0.5)*((*prop1).prop[i]+(*prop2).prop[i]);
   }
   int wall_size = GJP.VolNodeSites()/GJP.TnodeSites();
   rsrc = (Float *) smalloc(2 * wall_size * sizeof(Float) );
   for(i=0; i< 2*wall_size; i++){
      prop[i]=0.5*((*prop1).rsrc[i]+(*prop2).rsrc[i]);
   }

}

QPropWVolSrc::~QPropWVolSrc()
{
  char *fname = "~QPropWVolSrc()";
  char *cname = "QPropWVolSrc";
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "prop", prop);
  sfree(prop);
}


QPropWWallSrc::~QPropWWallSrc()
{
  char *fname = "~QPropWWallSrc()";
  char *cname = "QPropWWallSrc";
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "prop", prop);
  sfree(prop);
}



// Generate a prop from a random wall source

QPropWRandWallSrc::QPropWRandWallSrc(Lattice &lat, CgArg *arg, 
		int source_time, int seed, CommonArg* common_arg)
{
   char *fname = "QPropWRandWallSrc(L&, CgArg*, source_time, int seed)";
   char *cname = "QPropWRandWallSrc";
   VRB.Func(cname, fname);

   // Set the node size of the full (non-checkerboarded) fermion field
   //----------------------------------------------------------------
   int f_size = GJP.VolNodeSites() * lat.FsiteSize()/GJP.SnodeSites();
   int iter;
   Float true_res;


   // allocate space for the quark propagator
   //----------------------------------------------------------------
   prop = (WilsonMatrix*) smalloc(GJP.VolNodeSites() * sizeof(WilsonMatrix));
   if (prop == 0) ERR.Pointer(cname, fname, "prop");
   VRB.Smalloc(cname, fname, "prop", prop, 12 * f_size * sizeof(Float));

   FermionVectorTp source;
   FermionVectorTp sol;

   // make random source
   LRG.SetSigma(0.5);

   int wall_size = GJP.VolNodeSites()/GJP.TnodeSites();
   rsrc = (Float *) smalloc(2 * wall_size * sizeof(Float) );
   int i;
   for(i=0; i< 2*wall_size; i++)rsrc[i] = 0.0;

   int num_src=1, x[4];
   Float norm=1/sqrt(num_src);


   x[3] = source_time % GJP.TnodeSites();
   for(int n=0; n< num_src; n++)        {
     i = 0;
     for(x[2] = 0; x[2] < GJP.ZnodeSites(); x[2]++)
     for(x[1] = 0; x[1] < GJP.YnodeSites(); x[1]++)
     for(x[0] = 0; x[0] < GJP.XnodeSites(); x[0]++) {
       LRG.AssignGenerator(x);
       rsrc[i++] += LRG.Grand()*norm;     // real element
       rsrc[i++] += LRG.Grand()*norm;     // complex element
     }
   }

   for (int color = 0; color < GJP.Colors(); color++)
     for (int spin = 0; spin < 4; spin++){

        // initial guess
        sol.setVolSource();
        // set the source
        source.setWallSource(color, spin, source_time, rsrc);
        // do inversion
        //CG(lat, arg, source, sol);
        CG(lat, arg, source, sol, iter, true_res);

        // Collect solutions in propagator.
        int j;
        for(int i=0; i<f_size; i+=SPINOR_SIZE){
           j=i/SPINOR_SIZE; // lattice site
           prop[j].load_vec(spin, color, (wilson_vector &)sol[i]);
        }

        if(common_arg->results != 0){
          FILE *fp;
          if( (fp = fopen((char *)common_arg->results, "a")) == NULL ) {
             ERR.FileA(cname,fname, (char *)common_arg->results);
          }
          fprintf(fp, "Cg iters = %d true residual = %e\n",
                        iter, true_res);
          fclose(fp);
        }

   } // End spin-color loop

}

QPropWRandWallSrc::~QPropWRandWallSrc()
{
  char *fname = "~QPropWRandWallSrc()";
  char *cname = "QPropWRandWallSrc";
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "prop", prop);
  VRB.Sfree(cname, fname, "rsrc", rsrc);
  sfree(prop);
  sfree(rsrc);
}

void QPropWRandWallSrc::Delete()
{
  char *fname = "Delete()";
  char *cname = "QPropWRandWallSrc";
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "prop", prop);
  VRB.Sfree(cname, fname, "rsrc", rsrc);
  sfree(prop); prop = NULL;
  sfree(rsrc); rsrc = NULL;
}

WilsonMatrix& QPropW::operator[](int i)
{
     return prop[i];
}

const Rcomplex& QPropWRandWallSrc::rand_src(int i)
{
     return (const Rcomplex&)rsrc[2*i];
}

const Float& FermionVectorTp::operator[](int i)
{
     return fv_[i];
}

FermionVectorTp::FermionVectorTp()
{
  char *fname = "FermionVectorTp()";
  char *cname = "FermionVectorTp";
  VRB.Func(cname, fname);

  // allocate space for source
  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  fv_ = (Float *) smalloc(fv_size * sizeof(Float));
  if(fv_ == 0) ERR.Pointer(cname, fname, "fv_");
  VRB.Smalloc(cname,fname, "fv_", fv_, fv_size * sizeof(Float));

}

FermionVectorTp::~FermionVectorTp()
{
  char *fname = "~FermionVectorTp()";
  char *cname = "FermionVectorTp";
  VRB.Func(cname, fname);
  VRB.Sfree(cname, fname, "fv_", fv_);
  sfree(fv_);
}

// here x, y, z, t run over the whole lattice, which means
// 0 <= x <= x_nodes * x_node_sites
// ...
// 0 <= t <= t_nodes * t_node_sites
void FermionVectorTp::setPointSource(int color, int spin,
int x, int y, int z, int t)
{
  char *fname = "setPointSource()";
  char *cname = "FermionVectorTp";
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
    *((rfloat *)fv_ + i) = 0;

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
    *((rfloat *)fv_ + 2 * (color + GJP.Colors() * (spin + 4 * (
    localX + GJP.XnodeSites() * (
    localY + GJP.YnodeSites() * (
    localZ + GJP.ZnodeSites() * localT)))))) = 1.0;
}

// Zero at every color,spin,space-time point
void FermionVectorTp::setVolSource()
{
  char *fname = "setVolSource()";
  char *cname = "FermionVectorTp";
  VRB.Func(cname, fname);

  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  for (int i = 0; i < fv_size; i++) {
     *((rfloat *)fv_ + i)= 0.0;
  }
}

// Unit color/spin source at every space-time point
void FermionVectorTp::setVolSource(int color, int spin)
{
  char *fname = "setVolSource()";
  char *cname = "FermionVectorTp";
  VRB.Func(cname, fname);

  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
    "Color index out of range: color = %d\n", color);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
    "Spin index out of range: spin = %d\n", spin);

  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  for (int i = 0; i < fv_size; i++) {
     *((rfloat *)fv_ + i)= 0.0;
     if(i%SPINOR_SIZE != 2*( color + COLORS*spin ) )continue;
     *((rfloat *)fv_ + i) = 1.0;
  }
}


// Unit color/spin source at every space point on time_slice
void FermionVectorTp::setWallSource(int color, int spin, int source_time)
{
  char *fname = "setWallSource(color,spin,source_time)";
  char *cname = "FermionVectorTp";
  VRB.Func(cname, fname);

  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
    "Color index out of range: color = %d\n", color);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
    "Spin index out of range: spin = %d\n", spin);

#ifdef PARALLEL
  int my_node=CoorT();
#endif
  int ts_node = (int)(source_time/GJP.TnodeSites());
  int node_ts = source_time%GJP.TnodeSites();
  int wall_size = GJP.VolNodeSites()/GJP.TnodeSites()*SPINOR_SIZE;
  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  for (int i = 0; i < fv_size; i++) {
     *((rfloat *)fv_ + i)= 0.0;
#ifdef PARALLEL
     if(my_node != ts_node)continue; 
#endif
     if(i< wall_size*node_ts || i>=wall_size*(node_ts+1) )continue;
     if(i%SPINOR_SIZE != 2*( color + COLORS*spin ) )continue;
     *((rfloat *)fv_ + i) = 1.0;
  }
}

// Unit color/spin source at every space point on time_slice
void FermionVectorTp::setWallSource(int color, int spin, int source_time, Float* src)
{
  char *fname = "setWallSource(color,spin,source_time,Float*)";
  char *cname = "FermionVectorTp";
  VRB.Func(cname, fname);

  if (color < 0 || color >= GJP.Colors())
    ERR.General(cname, fname,
    "Color index out of range: color = %d\n", color);

  if (spin < 0 || spin > 3)
    ERR.General(cname, fname,
    "Spin index out of range: spin = %d\n", spin);
#ifdef PARALLEL
  int my_node=CoorT();
#endif
  int ts_node = (int)(source_time/GJP.TnodeSites());
  int node_ts = source_time%GJP.TnodeSites();
  int wall_size = GJP.VolNodeSites()/GJP.TnodeSites()*SPINOR_SIZE;
  int fv_size = GJP.VolNodeSites() * GJP.Colors() * 8;
  int i, ii;
  for (i=0,ii=0; i < fv_size; i+=2) {
     *((rfloat *)fv_ + i)= 0.0;
     *((rfloat *)fv_ + i+1)= 0.0;
#ifdef PARALLEL
     if(my_node != ts_node)continue;
#endif
     if(i< wall_size*node_ts || i>=wall_size*(node_ts+1) )continue;
     if(i%SPINOR_SIZE != 2*( color + COLORS*spin ) )continue;
     *((rfloat *)fv_ + i) = src[ii];
     *((rfloat *)fv_ + i+1) = src[ii+1];
     ii+=2;
  }
}

void FermionVectorTp::print() const
{
  char *fname = "print()";
  char *cname = "FermionVectorTp";
  VRB.Func(cname, fname);
  for (int i = 0; i <  GJP.VolNodeSites() * GJP.Colors() * 8; i++)
    VRB.Result(cname, fname,
   "fv_[%d] = %g\n", i, IFloat(*((rfloat *)fv_ + i)) );
}




void FermionVectorTp::setGFWallSource(Lattice &lat, int color,
  int spin, int dir, int where)
{

  char *cname = "QPropW";
  char *fname = "setGFWallSource()";
  VRB.Func(cname, fname);

  Matrix ** gm = lat.FixGaugePtr();

  int len;     //the local (on processor) length in "dir" direction
  int nproc;   // total number of processors in d_ direction
  int lproc;   // local processor coordinate in d_ direction
               // 0 <= lproc <= nproc

  switch(dir) {
    case 0: 
      len = GJP.XnodeSites();
      nproc = GJP.Xnodes();
      lproc = GJP.XnodeCoor();
      break;
    case 1: 
      len = GJP.YnodeSites();
      nproc = GJP.Ynodes();
      lproc = GJP.YnodeCoor();
      break;
    case 2: 
      len = GJP.ZnodeSites();
      nproc = GJP.Znodes();
      lproc = GJP.ZnodeCoor();
      break;
    case 3: 
      len = GJP.TnodeSites();
      nproc = GJP.Tnodes();
      lproc = GJP.TnodeCoor();
      break;
    // trap for wrong direction
    //-----------------------------------------------------------
    default:
      ERR.General(cname, fname, "bad argument: dir = %d\n", dir);
  }
    
  // trap for wrong position
  //-------------------------------------------------------------
  if (where < 0 || where >= len * nproc)
    ERR.General(cname, fname, 
    "argument out of range: where = %d\n", where);
  
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
    *((IFloat *)fv_ + i) = 0.0;

  // find out if this node overlaps with the hyperplane
  // in which the wall source sits
  int has_overlap = 0;
  if (lproc * len <= where && where < (lproc + 1) * len)
    has_overlap = 1;

  if (has_overlap) {
    int local = where % len; // on processor coordinate of 
                             // source hyperplane
    int x, y, z, t, c;
    Matrix *pM = gm[local];

    switch (dir) {
      case 0:
      //for (x = 0; x < GJP.XnodeSites(); x++)
        for (y = 0; y < GJP.YnodeSites(); y++)
        for (z = 0; z < GJP.ZnodeSites(); z++)
        for (t = 0; t < GJP.TnodeSites(); t++)
        for (c = 0; c < GJP.Colors(); c++)
        {
          // the real part
          *((IFloat *)fv_ + 
          0 + 2 * (
          c + GJP.Colors() * (
          spin + 4 * (
          local + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * t)))))) = 
          *((IFloat *)pM + 
          0 + 2 * (
          c + GJP.Colors() * (
          color + GJP.Colors() * (
          y + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * t)))));
          // the imaginary part
          *((IFloat *)fv_ + 
          1 + 2 * (
          c + GJP.Colors() * (
          spin + 4 * (
          local + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * t)))))) = 
          - *((IFloat *)pM + 
          1 + 2 * (
          c + GJP.Colors() * (
          color + GJP.Colors() * (
          y + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * t)))));
        }
        break;

      case 1:
        for (x = 0; x < GJP.XnodeSites(); x++)
      //for (y = 0; y < GJP.YnodeSites(); y++)
        for (z = 0; z < GJP.ZnodeSites(); z++)
        for (t = 0; t < GJP.TnodeSites(); t++)
        for (c = 0; c < GJP.Colors(); c++)
        {
          // the real part
          *((IFloat *)fv_ + 
          0 + 2 * (
          c + GJP.Colors() * (
          spin + 4 * (
          x + GJP.XnodeSites() * (
          local + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * t)))))) = 
          *((IFloat *)pM + 
          0 + 2 * (
          c + GJP.Colors() * (
          color + GJP.Colors() * (
          x + GJP.XnodeSites() * (
          z + GJP.ZnodeSites() * t)))));

          // the imaginary part
          *((IFloat *)fv_ + 
          1 + 2 * (
          c + GJP.Colors() * (
          spin + 4 * (
          x + GJP.XnodeSites() * (
          local + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * t)))))) = 
          - *((IFloat *)pM + 
          1 + 2 * (
          c + GJP.Colors() * (
          color + GJP.Colors() * (
          x + GJP.XnodeSites() * (
          z + GJP.ZnodeSites() * t)))));
        }
        break;

      case 2:
        for (x = 0; x < GJP.XnodeSites(); x++)
        for (y = 0; y < GJP.YnodeSites(); y++)
      //for (z = 0; z < GJP.ZnodeSites(); z++)
        for (t = 0; t < GJP.TnodeSites(); t++)
        for (c = 0; c < GJP.Colors(); c++)
        {
          // the real part
          *((IFloat *)fv_ + 
          0 + 2 * (
          c + GJP.Colors() * (
          spin + 4 * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * (
          local + GJP.ZnodeSites() * t)))))) = 
          *((IFloat *)pM + 
          0 + 2 * (
          c + GJP.Colors() * (
          color + GJP.Colors() * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * t)))));

          // the imaginary part
          *((IFloat *)fv_ + 
          1 + 2 * (
          c + GJP.Colors() * (
          spin + 4 * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * (
          local + GJP.ZnodeSites() * t)))))) = 
          - *((IFloat *)pM + 
          1 + 2 * (
          c + GJP.Colors() * (
          color + GJP.Colors() * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * t)))));
        }
        break;

      case 3:
        for (x = 0; x < GJP.XnodeSites(); x++)
        for (y = 0; y < GJP.YnodeSites(); y++)
        for (z = 0; z < GJP.ZnodeSites(); z++)
      //for (t = 0; t < GJP.TnodeSites(); t++)
        for (c = 0; c < GJP.Colors(); c++)
        {
          // the real part
          *((IFloat *)fv_ + 
          0 + 2 * (
          c + GJP.Colors() * (
          spin + 4 * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * local)))))) = 
          *((IFloat *)pM + 
          0 + 2 * (
          c + GJP.Colors() * (
          color + GJP.Colors() * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * z)))));

          // the imaginary part
          *((IFloat *)fv_ + 
          1 + 2 * (
          c + GJP.Colors() * (
          spin + 4 * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * local)))))) = 
          - *((IFloat *)pM + 
          1 + 2 * (
          c + GJP.Colors() * (
          color + GJP.Colors() * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * z)))));
        }
        break;
    }
  }

}

void FermionVectorTp::setLandauWallSource(Lattice &lat, int color,
  int spin, int where)
{

  char *cname = "QPropW";
  char *fname = "setLandauWallSource()";
  VRB.Func(cname, fname);

  Matrix ** gm = lat.FixGaugePtr();

  int len;     //the local (on processor) length in "dir" direction
  int nproc;   // total number of processors in d_ direction
  int lproc;   // local processor coordinate in d_ direction
               // 0 <= lproc <= nproc

  len = GJP.TnodeSites();
  nproc = GJP.Tnodes();
  lproc = GJP.TnodeCoor();
    
  // trap for wrong position
  //-------------------------------------------------------------
  if (where < 0 || where >= len * nproc)
    ERR.General(cname, fname, 
    "argument out of range: where = %d\n", where);
  
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
    *((IFloat *)fv_ + i) = 0.0;

  // find out if this node overlaps with the hyperplane
  // in which the wall source sits
  int has_overlap = 0;
  if (lproc * len <= where && where < (lproc + 1) * len)
    has_overlap = 1;

  if (has_overlap) {
    int local = where % len; // on processor coordinate of 
                             // source hyperplane
    int x, y, z, c;
    Matrix *pM = gm[0];

        for (z = 0; z < GJP.ZnodeSites(); z++)
        for (y = 0; y < GJP.YnodeSites(); y++)
        for (x = 0; x < GJP.XnodeSites(); x++)
        for (c = 0; c < GJP.Colors(); c++)
        {
          // the real part
          *((IFloat *)fv_ + 
          0 + 2 * (
          c + GJP.Colors() * (
          spin + 4 * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * local)))))) = 
          *((IFloat *)pM + 
          0 + 2 * (
          c + GJP.Colors() * (
          color + GJP.Colors() * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * local))))));

          // the imaginary part
          *((IFloat *)fv_ + 
          1 + 2 * (
          c + GJP.Colors() * (
          spin + 4 * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * local)))))) = 
          - *((IFloat *)pM + 
          1 + 2 * (
          c + GJP.Colors() * (
          color + GJP.Colors() * (
          x + GJP.XnodeSites() * (
          y + GJP.YnodeSites() * (
          z + GJP.ZnodeSites() * local))))));
        }
    }
}


void FermionVectorTp::LandauGaugeFixSink(Lattice &lat, int dir)
{

  char *cname = "QPropW";
  char *fname = "GaugeFixSink()";
  VRB.Func(cname, fname);

  Matrix ** gm = lat.FixGaugePtr();

  Vector temp[1];
  Matrix* pM = gm[0];

  for (int t = 0; t < GJP.TnodeSites(); t++)
  for (int z = 0; z < GJP.ZnodeSites(); z++)
  for (int y = 0; y < GJP.YnodeSites(); y++)
  for (int x = 0; x < GJP.XnodeSites(); x++)
  {
      // the matrix offset
      //int j = 2 * GJP.Colors() * GJP.Colors() * (
      int j =  x + GJP.XnodeSites() * ( y + GJP.YnodeSites() * ( z + GJP.ZnodeSites() * t));
      for (int spin = 0; spin < 4; spin++)
      {
         // the vector offset
         int i= 2 * GJP.Colors() * ( spin + 4 * (
         x + GJP.XnodeSites() * (
         y + GJP.YnodeSites() * (
         z + GJP.ZnodeSites() * t)))) ;
         temp->CopyVec((Vector*)&fv_[i], 6);
         uDotXEqual((IFloat*)&fv_[i], (const IFloat*)&pM[j], (const IFloat*)temp);
      }
  }
}


void FermionVectorTp::GaugeFixSink(Lattice &lat, int dir)
{

  char *cname = "QPropW";
  char *fname = "GaugeFixSink()";
  VRB.Func(cname, fname);

  Matrix ** gm = lat.FixGaugePtr();


  for (int t = 0; t < GJP.TnodeSites(); t++)
  {

    Matrix* pM = gm[t];
    Vector temp[1];

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
         temp->CopyVec((Vector*)&fv_[i], 6);
         uDotXEqual((IFloat*)&fv_[i], (const IFloat*)&pM[j], (const IFloat*)temp);
      }
    }
  }
}



CPS_END_NAMESPACE
