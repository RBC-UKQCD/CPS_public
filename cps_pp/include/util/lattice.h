#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/lattice.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: lattice.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.10  2003/02/10 11:15:24  mcneile
//  I have added a place holder for the asqtad code.
//
//  Revision 1.9  2002/12/04 17:16:27  zs
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
//  Revision 1.8  2002/03/11 22:27:10  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.5.2.1  2002/03/08 16:36:49  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.5  2001/11/08 14:39:38  anj
//  A couple of monir portability tweaks.  Anj
//
//  Revision 1.4  2001/08/16 10:50:30  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:17  anj
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
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: lattice.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/lattice.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// lattice.h
//
// Header file for all lattice classes.
//
//------------------------------------------------------------------


#ifndef INCLUDED_LATTICE_H
#define INCLUDED_LATTICE_H

CPS_END_NAMESPACE
#include<util/enum.h>
#include<util/random.h>
#include<util/vector.h>
#include<util/smalloc.h>
#include<util/pmalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/data_types.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include<alg/cg_arg.h>
#include<alg/ghb_arg.h>
#include<alg/eig_arg.h>
CPS_START_NAMESPACE

class LinkBuffer;

//------------------------------------------------------------------
//
// Lattice is the base abstract class.
//
//------------------------------------------------------------------
class Lattice
{

 private:

    char *cname;    // Class name.

    static Matrix* gauge_field;
       // Pointer to the gauge field configuration.
  
    static int is_allocated;	
       // 0 = gauge field has not been allocated
       // 1 = gauge field has been allocated

    static int is_initialized;	
       // 0 = gauge field has not been initialized
       // 1 = gauge field has been initialized

    static StrOrdType str_ord;
       // The gauge field configuration storage order

    static int *g_upd_cnt;
      // Counts number of heat bath/metropolis update sweeps or
      // HMD/HMC (accepted) trajectories over life of job.
      // Allocated with pmalloc().
      // Manipulated using public functions in Lattice base class.

    static Float md_time;
      // Current molecular dynamics time.  Initialized to 0.0 by
      // the constructor.  Should be manipulated by AlgHmcPhi or AlgHmdR
      // using public functions in Lattice base class.  Time is in units
      // of molecular dynamics step size.

    static Matrix **fix_gauge_ptr;
      // Pointer to an array of pointers that point to the
      // various gauge fixed hyperplanes. It is set to 0
      // if the gauge matrices have not been calculated. 
      // The array pointers are set to 0 if the corresponding 
      // hyperplanes do not have any gauge fixing matrices 
      // calculated.

    static FixGaugeType fix_gauge_kind;
      // The kind of gauge fixing. It is set to FIX_GAUGE_NONE
      // if the gauge is not fixed.

 protected:
    
    static int node_sites[5];
    	// node_sites[0] = GJP.XnodeSite();
    	// node_sites[1] = GJP.YnodeSite();
    	// node_sites[2] = GJP.ZnodeSite();
    	// node_sites[3] = GJP.TnodeSite();
    	// node_sites[4] = GJP.SnodeSite();

    static int g_dir_offset[4];
    	// g_dir_offset[0] = 4;
    	// g_dir_offset[1] = g_dir_offset[0]*GJP.XnodeSites();
    	// g_dir_offset[2] = g_dir_offset[1]*GJP.YnodeSites();
    	// g_dir_offset[3] = g_dir_offset[2]*GJP.ZnodeSites();


    void *f_dirac_op_init_ptr;
      // A pointer that is used by the fermion classes to
      // point to data that need only be initialized by the fermion
      // class constructor and are needed by the relevant dirac
      // operator.

    void *aux0_ptr;
      // General purpose auxiliary pointer 0;

    void *aux1_ptr;
      // General purpose auxiliary pointer 1;


    // Added in by Ping for anisotropic lattices
    //------------------------------------------------------------------
    void MltFloatImpl(Float factor, int dir);
    // U_dir(x) *= factor  where dir = [0,1,2,3] as [x,y,z,t]

    LinkBuffer * link_buffer ;
    //For better performance when getting an offsite link.

    // change from phys_v4.0.0
    // the declarations for GetLink and GetLinkOld
    // used to be protected and listed above (before void MltFloatImpl)

 public:

    const Matrix * GetLink(const int *x, int mu) const;
      // GRF:  returns a reference to the link U_\mu(x), where x is
      // defined relative to the local site (0,0,0,0).
      // If the link is on-node, it returns a reference to
      // to the link.  If the link is off-node, it retrieves the link
      // into a static buffer and returns the reference to the buffer.
      // Since the buffer can be used by other routines as well as other
      // calls to this routine, as a general rule, the link should be
      // used immediately or else copied.

    const Matrix * GetLinkOld(Matrix *g_offset, const int *x,
             int dir, int mu) const;
      // get U_mu(x+dir)
      // GRF: renamed to avoid conflict with the more general
      // purpose function

    // end change from phys_v4.0.0 --> phys_v4.1.0

 public:
    friend class LinkBuffer;

    int LinkBufferIsEnabled(){return ((int) link_buffer);}

    int EnableLinkBuffer(int buf_sz);
      //create the LinkBuffer Object only when requested.

    void DisableLinkBuffer();   
      //delete the LinkBuffer Object when it's not in use. 

    const Matrix * GetBufferedLink(const int *x, int mu);
      //get link through the buffer.

    void ClearBufferedLink(const int * x, int mu);
      //delete the changed links in the buffer

    void ClearAllBufferedLink();
      //delete all the buffered links from the buffer
      //this must be called when lat.Unitarize() is used

    int IsOnNode(const int * x);
      //check if the site x is on node.
      //return 1 if on node, 0 if offnode.

    void PathOrdProdPlus(Matrix & mat, const int *x, const int * dirs, int n);
      //given the starting point x, the directions of each step on the path
      //and the number of steps. calculate the path ordered product of 
      //all the links along the path and add the result to mat 
      //each direction could be {0,1,2,3,4,5,6,7} which coresponds to 
      //the directions {n_x, n_y, n_z, n_t, -n_x, -n_y, -n_z, -n_t}
      //the result is returned in mat.

    void PathOrdProd(Matrix & mat, const int *x, const int * dirs, int n);
      //also calculates the path ordered product, but the result returned 
      //is that product of the path that end on the local node

 public: 

// Base class functions
//------------------------------------------------------------------

    Lattice(void);

    virtual ~Lattice(void);

    Matrix *GaugeField(void) const;
    	// Returns the pointer to the gauge field configuration.

    void GaugeField(Matrix *u);
        // Copies the array pointed to by u to the array
        // pointed to by gauge_field.

    int GsiteOffset(int *x) const
        { return x[0]*g_dir_offset[0]+x[1]*g_dir_offset[1]
	        +x[2]*g_dir_offset[2]+x[3]*g_dir_offset[3];  }
        // Sets the offsets for the canonical storage order
        // of the gauge field. x[i] is the ith coordinate
        // where i = {0,1,2,3} = {x,y,z,t}.

    void CopyGaugeField(Matrix* u);
        // Copies the array pointed to by gauge_field to the
        // array pointed to by u.

    StrOrdType StrOrd(void);
        // Returns the storage order.

    int Colors(void);
        // Returns the number of colors.	  

    int GsiteSize(void);
        // Returns the number of gauge field 
        // components (including real/imaginary) on a
        // site of the 4-D lattice.


    void Staple(Matrix& stap, int *x, int mu);
        // It calculates the staple field at x, mu.
        // The staple field is:
        //
        //      V_u(x) = \sum_v(!=u) {
        //		U_v(x+u) U_u(x+v)~ U_v(x)~
        //	     +  U_v(x+u-v)~ U_u(x-v)~ U_v(x-v)  }
        //
        // GRF: consider changing the function name
        // to Lattice::PlaqStaple() for consistency.

    void BufferedStaple(Matrix & stap, const int *x, int mu);
        //Buffered version of staple

    void RectStaple(Matrix& stap, int *x, int mu) ;
        // It calculates the rectangle field at x, mu.
        // The rectangle field is:
        //
        // \sum_{v != u} {
        //     U_u(x+u)    U_v(x+2u)    U_u(x+u+v)~ U_u(x+v)~  U_v(x)~
        //   + U_u(x+u)    U_v(x+2u-v)~ U_u(x+u-v)~ U_u(x-v)~  U_v(x-v)
        //   + U_v(x+u)    U_u(x+v)~    U_u(x-u+v)~ U_v(x-u)~  U_u(x-u)
        //   + U_v(x+u-v)~ U_u(x-v)~    U_u(x-u-v)~ U_v(x-u-v) U_u(x-u)
        //   + U_v(x+u)    U_v(x+u+v)   U_u(x+2v)~  U_v(x+v)~  U_v(x)~
        //   + U_v(x+u-v)~ U_v(x+u-2v)~ U_u(x-2v)~  U_v(x-2v)  U_v(x-v)
        // }
        //
    void BufferedRectStaple(Matrix& stap, const int *x, int mu);
        //buffered version of RectStaple 

    void RectStaple1(Matrix& stap, int *x, int mu);
        // it calculates the flat 6-1 link staple using the PathOrdProdPlus routine
	// \sum_{ v!= +/-u }{
	//    U_v(x+u) U_v(x+u+v) U_{-u}(x+u+v+v) U_{-v}(x+v+v) U_{-v}(x+v)
	// +  U_v(x+u) U_{-u}(x+v+u) U_{-u}(x+v)  U_{-v}(x+v-u) U_u(x-u)
	// +  U_u(x+u) U_v(x+u+u) U_{-u}(x+u+u+v) U_{-u}(x+u+v) U_{-v}(x+v)

    void ChairStaple(Matrix& stap, int *x, int mu);
         //it calculates the chair shaped 6-1 link staple at x, mu
         // \sum_{w != +/-u, w !=+/-v, w != +/-u}{
         //    U_w(x+u) U_v(x+u+w) U_{-w}(x+u+v+w) U_{-u}(x+u+v) U_{-v}(x+v)  
         //  + U_w(x+u) U_v(x+u+w) U_{-u}(x+u+v+w) U_{-v}(x+v+w) U_{-w}(x+w)
         //  + U_w(x+u) U_{-u}(x+u+w) U_v(x+w) U_{-w}(x+w+v) U_{-v}(x+v)
         // }

    void BufferedChairStaple(Matrix &stap, const int *x, int mu);
        //buffered version of ChairStaple
    
    void CubeStaple(Matrix &stap, const int *x, int mu);
        //it calculates the cube shaped 6-1 link staple at x, mu
        // \sum_{w != +/-u, w !=+/-v, w != +/-u}{
        //   U_v(x+u) U_w(x+u+v) U_{-u}(x+u+v+w) U_{-v}(x+v+w) U_{-w}(x+w)
        // }
            
    void BufferedCubeStaple(Matrix &stap, const int *x, int mu);  
        //buffered version of CubeStaple

    virtual void AllStaple(Matrix &stap, const int *x, int mu)=0;
        //pure virtual function must be implemented for the gauge actions 
        //derived from it.
        //given a link calculate all of its staples, depending on the 
        //lattice type.
        //this is used in heatbath.

    void Plaq(Matrix &plaq, int *x, int mu, int nu) const;
        // Added by Ping for the purpose of debugging now, may be more useful
        // later on. It calculates  the plaquette 
        //   U_u(x) U_v(x+u) U_u(x+v)~ U_v(x)~

    Float ReTrPlaq(int *x, int mu, int nu) const;
        // It calculates the real part of the trace of the plaquette 
        // field at site x, mu, nu with mu < nu.
        // The plaquette field is:
        //
        //   U_u(x) U_v(x+u) U_u(x+v)~ U_v(x)~

    Float SumReTrPlaqNode(void) const;
       // It calculates the sum of the real part of the trace of the 
       // plaquette field at each site of the node sublattice.

    Float SumReTrPlaq(void) const;
       // It calculates the sum of the real part of the trace of the 
       // plaquette field at each site of the whole lattice

    Float ReTrRect(int *x, int mu, int nu) const;
       // It calculates the real part of the trace of the rectangle
       // field at site x, in the (mu, nu) plane with the long axis
       // of the rectangle in the mu direction.
       // The rectangle field is:
       //
       //   U_u(x) U_u(x+u) U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
       //

    Float SumReTrRectNode(void) const;
       // It calculates the sum of the real part of the trace of the
       // rectangle field at each site of the node sublattice.

    Float SumReTrRect(void) const;
       // It calculates the sum of the real part of the trace of the
       // rectangle field at each site of the whole lattice.
    
    Float ReTrLoop(const int *x, const int *dir,  int length) ;
    //-------------------------------------------------------------------
    // ReTrLoop(int *x, int *dir,int length):
    //   It calculates the real part of the trace of the loop at site x
    //   specified by the list of directions in dir
    //   length is the length of the loop 
    //
    // Warning!!:
    //   The user is responcible for handing in directions that close a loop!
    //--------------------------------------------------------------------
    
    Float SumReTrCubeNode(void) ;
    //----------------------------------------------------------------------
    // SumReTrLoopCube()
    //   It calculates the sum of the real part of the trace of the
    //   Cube field at each site of the node sublattice.
    //
    //   The Cube loop is: mu nu rho -mu -nu -rho
    //-----------------------------------------------------------------------

    Float SumReTrCube(void) ;
    //-----------------------------------------------------------------------
    // SumReTrCube()
    //   It calculates the sum of the real part of the trace of the
    //   Cube field at each site of the whole lattice
    //-----------------------------------------------------------------------

  // Added in by Ping for anisotropic lattices
  //------------------------------------------------------------------
  Float AveReTrPlaqNodeNoXi(void) const;
  // Normalization:  1 for ordered links
  // Average over plaq's perpendicular to the special anisotropic dir.

  Float AveReTrPlaqNodeXi(void) const;
  // Normalization:  1 for ordered links
  // Average over plaq's parallel to the special anisotropic dir.

  Float AveReTrPlaqNoXi(void) const;
  // Normalization:  1 for ordered links
  // Average over plaq's perpendicular to the special anisotropic dir.

  Float AveReTrPlaqXi(void) const;
  // Normalization:  1 for ordered links
  // Average over plaq's parallel to the special anisotropic dir.

  void MltFloat(Float factor, int dir)      {
    if (factor != 1.0)    MltFloatImpl(factor, dir);
  }    
  // U_dir(x) *= factor  where dir = [0,1,2,3] as [x,y,z,t]
  // Handle all kinds of storage order correctly.

    void Reunitarize(void);
    	// Re-unitarize the gauge field configuration.

    void Reunitarize(Float &dev, Float &max_diff);
    	// Re-unitarize the gauge field configuration
	// and return:
	// dev = sqrt( Sum_i [ (U(i) - V(i))^2 ] / (Vol*4*18) ),
        // max_diff = Max_i[ |U(i) - V(i)| ]
        // where U(i), V(i) is the gauge field before and after 
        // reunitarization. The index i runs over all components of
        // the gauge field.

    int MetropolisAccept(Float delta_h);
        // 0 reject, 1 accept. If delta_h < 0 it accepts
        // unconditionally.

    void EvolveGfield(Matrix *mom, Float step_size);
        // It evolves the gauge field by step_size using
        // the canonical momentum mom

    Float MomHamiltonNode(Matrix *momentum);
        // The conjugate momentum Hamiltonian of the node sublattice.

    void Convert(StrOrdType new_str_ord,
			Vector *f_field_1,
			Vector *f_field_2);
        // If str_ord is not the same as 
        // new_str_ord then it converts the gauge field
        // configuration and the two fermion fields f_field_1, 
        // f_field_2 to new_str_ord.

    void Convert(StrOrdType new_str_ord);
        // If str_ord is not the same as 
        // new_str_ord then it converts the gauge field
        // configuration to new_str_ord.

    void RandGaussAntiHermMatrix(Matrix *mat, Float sigma2);
        // It produces an anti-Hermitian matrix for each
        // site of the lattice with random
        // entries weighted according to 
	// exp(- Tr(mat^2)) / (2 * sigma2)

    void RandGaussVector(Vector *vect, Float sigma);

    void RandGaussVector(Vector *vect, Float sigma,
                        FermionFieldDimension frm_field_dim);
        // This version assumes the fermion field spans the whole lattice

    void RandGaussVector(Vector *vect, Float sigma, int num_chckbds,
                              FermionFieldDimension frm_field_dim = FIVE_D);
        // It produces 3-vectors in canonical storage order
        // for each site of the lattice with random
        // entries weighted according to
        // exp(- Tr(mat^2)) / (2 * sigma2)
        // It defaults to producing fermionic fields over the entire lattice,
        // but can also produce even and odd checkerboards.
        // If keep_snodes == 1, then we divide the site size by SnodeSites()

    void SetGfieldOrd(void);
    	// Sets the gauge field to the identity

    void SetGfieldDisOrd(void);
    	// Sets the gauge field to disordered (random) values

    int GupdCnt(void);
      // Returns the value of *g_upd_cnt.

    int GupdCnt(int set_val);
      // Sets the value of *g_upd_cnt to set_val and returns that value.

    int GupdCntInc(int inc_val = 1);
      // Increments the *g_upd_cnt by inc_val and returns new value.

    Float MdTime(void);
      // Returns the value of md_time.

    Float MdTime(Float set_val);
      // Sets the value of md_time to set_val and returns that value.

    Float MdTimeInc(Float inc_val = 0.5);
      // Increments md_time by inc_val and returns new value.

    void *FdiracOpInitPtr(void)
      { return f_dirac_op_init_ptr; }
      // Returns a pointer that is used by the fermion classes to
      // point to data that need only be initialized by the fermion
      // class constructor and are needed by the relevant dirac
      // operator.

    void FixGaugeAllocate(FixGaugeType GaugeType,int NHplanes=0,int *Hplanes=0);
        // Allocates memory for the gauge fixing matrices
        //
	// FixGaugeType GaugeKind - gauge type. Numerically equal to the number
	// of the direction (ordered in the "canonical" way X=0,
        // Y=1, Z=2, T=3) orthogonal to the three dimensions used
        // in the Coulomb gauge fixing condition in (2). For the
        // Landau gauge a negative number is used.
        //
        // int NHplanes - number of the hyperplanes to fix on the whole
        // machine (not only current node). Not used in the
        // Landau gauge. If set to zero in Coulomb gauge then
        // treated as a request to use all hyperplanes on all nodes.
        //
        // int *Hplanes - list of NHplanes positions of the hyperplanes to fix
        // along the direction orthogonal to them. Not used in the
        // Landau gauge and when NHplanes is set to zero.

    int FixGauge(Float StopCond, int MaxIterNum);
        // Fixes the gauge and returns number of iterations for this node.
        // FixGaugeAllocate must be called first.
        //
        // Float StopCond is the stopping condition;
        //
        // int MaxIterNum - issues a warning if reached

    void FixGaugeFree(void);
        // Free memory for the gauge fixing matrices

    Matrix **FixGaugePtr(void);
      // Returns fix_gauge_ptr (pointer to an array of pointers 
      // that point to the various gauge fixed hyperplanes.

    FixGaugeType FixGaugeKind(void);
      // Returns fix_gauge_kind (the kind of gauge fixing).

    void *Aux0Ptr(void);
      // Returns the general purpose auxiliary pointer 0;

    void *Aux1Ptr(void);
      // Returns the general purpose auxiliary pointer 1;

    void GsoCheck(void);
      // If GJP.Snodes() == 1 it just returns.
      // If GJP.Snodes() != 1 it checks that the "spread-out"
      // gauge field is identical along all s-slices by comparing 
      // the checksum and plaquette value. This situation arises for
      // DWF with the s direction spread out across many processors.
      // If either the checksum or plaquette are not identical it
      // exits with an error.

    void SoCheck(Float num);
      // If GJP.Snodes() == 1 it just returns.
      // If GJP.Snodes() != 1 it checks that the value of num
      // is identical along all s-slices. This situation arises for
      // DWF with the s direction spread out across many processors.
      // If the num along all s-slices is not identical it exits
      // with an error.


// Gauge action related virtual functions.
//------------------------------------------------------------------


// Fermion action related virtual functions.
//------------------------------------------------------------------
    virtual void Gamma5(Vector *v_out, Vector *v_in, int num_sites);
      // v_out = gamma_5 v_in

    virtual void Ffour2five(Vector *five, Vector *four, 
			    int s_r, int s_l);
        // This is "truly" implemented only in the Fdwf derived class

    virtual void Ffive2four(Vector *four, Vector *five, 
			    int s_r, int s_l);
        // This is "truly" implemented only in the Fdwf derived class

    virtual void Freflex (Vector *out, Vector *in) {}
       // This is "truly" implemented only in the Fdwf derived class

    virtual void Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int dir_flag);
    // dir_flag is flag which takes value 0 when all direction contribute to D
    // 1 - when only the special anisotropic direction contributes to D,
    // 2 - when all  except the special anisotropic direction.
    // Currently this function is implemented only in the Fstag class

// Gauge action related pure virtual functions
//------------------------------------------------------------------
    virtual GclassType Gclass(void) = 0;
        // It returns the type of gauge class

    virtual void GactionGradient(Matrix &grad, int *x, int mu) = 0;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).

    virtual void EvolveMomGforce(Matrix *mom, Float step_size) = 0;
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    virtual Float GhamiltonNode(void) = 0;
        // The pure gauge Hamiltonian of the node sublattice



// Fermion action related pure virtual functions 
//------------------------------------------------------------------
    virtual FclassType Fclass(void) = 0;
        // It returns the type of fermion class

    virtual int FsiteOffsetChkb(const int *x) const = 0;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is not the canonical one but it is particular
        // to the fermion type. This function is not
        // relevant to fermion types that do not
        // use even/odd checkerboarding. x[i] is the 
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    virtual int FsiteOffset(const int *x) const = 0;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is the canonical one. X[I] is the
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    virtual int ExactFlavors(void) = 0;
        // Returns the number of exact flavors of the matrix that
        // is inverted during a molecular dynamics evolution.

    virtual int SpinComponents(void) = 0;
        // Returns the number of spin components.

    virtual int FsiteSize(void) = 0;
        // Returns the number of fermion field 
        // components (including real/imaginary) on a
        // site of the 4-D lattice.

    virtual int FchkbEvl(void) = 0;
        // 0 -> If no checkerboard is used for the evolution
        //      or the CG that inverts the evolution matrix.
	// 1 -> If the fermion fields in the evolution
        //      or the CG that inverts the evolution matrix
	//      are defined on a single checkerboard (half the 
	//      lattice).

    virtual int FmatEvlInv(Vector *f_out, Vector *f_in, 
			   CgArg *cg_arg, 
                           Float *true_res,
			   CnvFrmType cnv_frm = CNV_FRM_YES) = 0;
        // It calculates f_out where A * f_out = f_in and
        // A is the preconditioned (if relevant) fermion matrix that
        // appears in the HMC evolution (typically some preconditioned 
        // form of [Dirac^dag Dirac]). The inversion is done
	// with the conjugate gradient. cg_arg is the structure
        // that contains all the control parameters, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
	// f_in and f_out are defined on a checkerboard.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
	// The function returns the total number of CG iterations.

    virtual int FmatEvlInv(Vector *f_out, Vector *f_in, 
			   CgArg *cg_arg, 
			   CnvFrmType cnv_frm = CNV_FRM_YES) = 0;
        // Same as original but with true_res=0;

    virtual int FmatInv(Vector *f_out, Vector *f_in, 
			CgArg *cg_arg, 
                        Float *true_res,
			CnvFrmType cnv_frm = CNV_FRM_YES,
			PreserveType prs_f_in = PRESERVE_YES) = 0;
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix (Dirac operator). The inversion
	// is done with the conjugate gradient. cg_arg is the 
        // structure that contains all the control parameters, f_in 
        // is the fermion field source vector, f_out should be set 
        // to be the initial guess and on return is the solution.
	// f_in and f_out are defined on the whole lattice.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
        // cnv_frm is used to specify if f_in should be converted 
        // from canonical to fermion order and f_out from fermion 
        // to canonical. 
        // prs_f_in is used to specify if the source
        // f_in should be preserved or not. If not the memory usage
        // is less by the size of one fermion vector or by the size 
        // of one checkerboard fermion vector (half a fermion vector).
        // For staggered fermions f_in is preserved regardles of
        // the value of prs_f_in. 
	// The function returns the total number of CG iterations.

    virtual int FmatInv(Vector *f_out, Vector *f_in, 
			CgArg *cg_arg, 
			CnvFrmType cnv_frm = CNV_FRM_YES,
			PreserveType prs_f_in = PRESERVE_YES) = 0;
        // Same as original but with true_res=0;

    virtual int FeigSolv(Vector **f_eigenv, Float lambda[], 
			 Float chirality[], int valid_eig[],
			 Float **hsum,
			 EigArg *eig_arg, 
			 CnvFrmType cnv_frm = CNV_FRM_YES) = 0;
        // It finds the eigenvectors and eigenvalues of A where
        // A is the fermion matrix (Dirac operator). The solution
	// uses Ritz minimization. eig_arg is the 
        // structure that contains all the control parameters, f_eigenv
        // are the fermion field source vectors which should be
        // defined initially, lambda are the eigenvalues returned 
        // on solution. f_eigenv is defined on the whole lattice.
        // hsum are projected eigenvectors.
	// The function returns the total number of Ritz iterations.

    virtual void SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
			Float mass) = 0;
	// It sets the pseudofermion field phi from frm1, frm2.

    virtual void EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size) = 0;
        // It evolves the canonical momentum mom by step_size
        // using the fermion force. 

    virtual Float FhamiltonNode(Vector *phi, Vector *chi) = 0;
        // The fermion Hamiltonian of the node sublattice.
        // chi must be the solution of Cg with source phi.	       

    virtual void Fconvert(Vector *f_field, 
		   	  StrOrdType to,
			  StrOrdType from) = 0;
        // Convert fermion field f_field from -> to


// Bosonic action related pure virtual functions.
//---------------------------------------------------------1---------

    virtual Float BhamiltonNode(Vector *boson, Float mass) = 0;
        // The boson Hamiltonian of the node sublattice.

};


//------------------------------------------------------------------
//
// Gnone is derived from Lattice. Its functions act
// as if there is no gauge action i.e. beta = 0.
//
//------------------------------------------------------------------
class Gnone : public virtual Lattice
{

 private:
    char *cname;    // Class name.

 public:

    Gnone(void);

    virtual ~Gnone(void);

    GclassType Gclass(void);
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
        // It calculates the gauge force at site x and direction mu.

    void EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode(void);
       // The pure gauge Hamiltonian of the node sublattice.

    void AllStaple(Matrix &stap, const int *x, int mu);

};


//------------------------------------------------------------------
//
// Gwilson is derived from Lattice and is relevant to the 
// standard Wilson single plaquette action.
//
//------------------------------------------------------------------
class Gwilson : public virtual Lattice
{

 private:
    char *cname;    // Class name.

 public:

    Gwilson(void);

    virtual ~Gwilson(void);

    GclassType Gclass(void);
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
        // It calculates the gauge force at site x and direction mu.

    void EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode(void);
       // The pure gauge Hamiltonian of the node sublattice.

    void AllStaple(Matrix &stap, const int *x, int mu);

};


//------------------------------------------------------------------
//
// GpowerPlaq is derived from Lattice and is relevant to the 
// power plaquette action. This action is the same as
// the standard Wilson action with the irrelevant power plaquette
// term added to it. The full action is:
//
// Sum_p [ beta * { -Tr[U_p]/3} + ( {1 - Tr[U_p]/3} / c )^k ]
//
// with c = GJP.PowerPlaqCutoff() and k = GJP.PowerPlaqExponent()
//
// This action supresses plaquettes with {1 - ReTr[U_p]/3} > c 
// and threfore reduces lattice dislocations.
//
//------------------------------------------------------------------
class GpowerPlaq : public virtual Lattice
{

 private:
    char *cname;    // Class name.

 public:

    GpowerPlaq(void);

    virtual ~GpowerPlaq(void);

    GclassType Gclass(void);
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
        // It calculates the gauge force at site x and direction mu.

    void EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode(void);
       // The pure gauge Hamiltonian of the node sublattice.

    void PowerStaple(Matrix& pstap, int *x, int mu);
        // It calculates the staple field at x, mu.
        // The staple field is:
        //
        // V_u(x) = \sum_v(!=u) {
        //      ps(x,u,v)   * [ U_v(x+u) U_u(x+v)~ U_v(x)~     ]
        //    + ps(x-v,u,v) * [ U_v(x+u-v)~ U_u(x-v)~ U_v(x-v) ] }
        //
        // where
        //
        // ps(x,u,v) = 
        // beta + {k/c} * { (1 - ReTr[U_p(x,u,v)]/3) / c }^(k-1)
        //
        // with c = GJP.PowerPlaqCutoff() and
        //      k = GJP.PowerPlaqExponent()

    Float PowerPlaq(int *x, int mu, int nu) const;
        // It calculates the power plaquette 
        // field at site x, mu, nu with mu < nu.
        // The power plaquette field is:
        //
        // pp(x,u,v) = { (1 - ReTr[U_p(x,u,v)]/3) / c }^k
        //
        // with c = GJP.PowerPlaqCutoff() and
        //      k = GJP.PowerPlaqExponent()

    Float SumPowerPlaqNode(void) const;
       // It calculates the sum of the power plaquette  
       // field at each site of the node sublattice.

    Float SumPowerPlaq(void) const;
       // It calculates the sum of the power plaquette 
       // field at each site of the whole lattice

    void AllStaple(Matrix &stap, const int *x, int mu);
};

//------------------------------------------------------------------
//
// GimprRect is derived from Lattice and is relevant to the
// action which contains the standard Wilson plaquette operator
// plus the second order rectangle operator.
//
//------------------------------------------------------------------
class GimprRect : public virtual Lattice
{

 private:
    char *cname;    // Class name.

    Float plaq_coeff; // - GJP.Beta() * ( 1.0 - 8.0 * GJP.C1() ) / 3.0

    Float rect_coeff; // - GJP.Beta() * (             GJP.C1() ) / 3.0

 public:

    GimprRect(void);

    virtual ~GimprRect(void);

    GclassType Gclass(void);
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
        // It calculates the gauge force at site x and direction mu.
        // Typical implementation has this func called with
        // Matrix &force = *mp0.  GactionGradient typically uses
        // mp1 thru mp4, so be careful.

    void EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode(void);
       // The pure gauge Hamiltonian of the node sublattice.

    void AllStaple(Matrix &stap, const int *x, int mu);

};


//------------------------------------------------------------------
//
// GpowerRect is derived from Lattice and is relevant to the
// action which contains the standard Wilson plaquette operator
// plus the second order rectangle operator plus a power
// plaquette term plus a power rectangle term.
//
// The full action is:
// 
// (  Sum_p [ c_0*beta*{ -Tr[U_p]/3} + ( {1 - Tr[U_p]/3} / c )^k ]
//  + Sum_r [ c_1*beta*{ -Tr[U_r]/3} + ( {1 - Tr[U_r]/3} / c )^k ] )
// with c = GJP.PowerPlaqCutoff(), k = GJP.PowerPlaqExponent()
// c_0 = 1 - 8 * c_1, c_1 = GJP.C1()
// This action supresses plaquettes with {1 - ReTr[U_p]/3} > c
// and rectangles with {1 - ReTr[U_r]/3} > c
// and threfore reduces lattice dislocations.
//
//------------------------------------------------------------------
class GpowerRect : public virtual Lattice
{

 private:
    char *cname;    // Class name.

    Float plaq_coeff; // - GJP.Beta() * ( 1.0 - 8.0 * GJP.C1() ) / 3.0

    Float rect_coeff; // - GJP.Beta() * (             GJP.C1() ) / 3.0

 public:

    GpowerRect(void);

    virtual ~GpowerRect(void);

    GclassType Gclass(void);
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
        // It calculates the gauge force at site x and direction mu.
        // Typical implementation has this func called with
        // Matrix &force = *mp0.  GactionGradient typically uses
        // mp1 thru mp4, so be careful.

    void EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode(void);
       // The pure gauge Hamiltonian of the node sublattice.

    void PowerStaple(Matrix& pstap, int *x, int mu);
        // It calculates the staple field at x, mu.
        // The staple field is:
        //
        // V_u(x) = \sum_v(!=u) {
        //      ps(x,u,v)   * [ U_v(x+u) U_u(x+v)~ U_v(x)~     ]
        //    + ps(x-v,u,v) * [ U_v(x+u-v)~ U_u(x-v)~ U_v(x-v) ] }
        //
        // where
        //
        // ps(x,u,v) = 
        // beta + {k/c} * { (1 - ReTr[U_p(x,u,v)]/3) / c }^(k-1)
        //
        // with c = GJP.PowerPlaqCutoff() and
        //      k = GJP.PowerPlaqExponent()

    Float PowerPlaq(int *x, int mu, int nu) const;
        // It calculates the power plaquette 
        // field at site x, mu, nu with mu < nu.
        // The power plaquette field is:
        //
        // pp(x,u,v) = { (1 - ReTr[U_p(x,u,v)]/3) / c }^k
        //
        // with c = GJP.PowerPlaqCutoff() and
        //      k = GJP.PowerPlaqExponent()

    Float SumPowerPlaqNode(void) const;
       // It calculates the sum of the power plaquette  
       // field at each site of the node sublattice.

    Float SumPowerPlaq(void) const;
       // It calculates the sum of the power plaquette 
       // field at each site of the whole lattice

    void PowerRectStaple(Matrix& pstap, int *x, int mu);
        // It calculates the rectangle staple field at x, mu.

    Float PowerRect(int *x, int mu, int nu) const;
        // It calculates the power rectangle 
        // field at site x, mu, nu with mu < nu.
        // The power plaquette rectangle is:
        //
        // pp(x,u,v) = { (1 - ReTr[U_r(x,u,v)]/3) / c }^k
        //
        // with c = GJP.PowerPlaqCutoff() and
        //      k = GJP.PowerPlaqExponent()

    Float SumPowerRectNode(void) const;
       // It calculates the sum of the power rectangle  
       // field at each site of the node sublattice.

    Float SumPowerRect(void) const;
       // It calculates the sum of the power rectangle 
       // field at each site of the whole lattice
    void AllStaple(Matrix &stap, const int *x, int mu);
};

//------------------------------------------------------------------
//
// GimprOLSym is derived from Lattice. It implements the
// One Loop Symanzik improved gauge action.
//
//------------------------------------------------------------------
class GimprOLSym : public virtual Lattice
{

 private:
    char *cname;    // Class name.

    Float rect_coeff; // loop 2x1
    Float cube_coeff; // loop 1x1x1
    Float minus_beta_over_3 ;// coefficient needed to calculate the force 

 public:

    GimprOLSym(void);

    virtual ~GimprOLSym(void);

    GclassType Gclass(void);
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
        // It calculates the gauge force at site x and direction mu.
        // Typical implementation has this func called with
        // Matrix &force = *mp0.  GactionGradient typically uses
        // mp1 thru mp4, so be careful.

    void EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode(void);
       // The pure gauge Hamiltonian of the node sublattice.

    void AllStaple(Matrix &stap, const int *x, int mu);

};

//------------------------------------------------------------------
//
//
// Fnone is derived from Lattice. Its functions do nothing
// and return values as if there is no fermion action or
// fermion fields. The number of spin components is zero
// The site size of the fermion array FsiteSize() is 
// set to 1 so that memory allocation would proceed normally.
//
//------------------------------------------------------------------
class Fnone : public virtual Lattice
{
 private:
    char *cname;    // Class name.
    
 public:

    Fnone(void);

    virtual ~Fnone(void);

    FclassType Fclass(void);
        // It returns the type of fermion class

    int FsiteOffsetChkb(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is not the canonical one but it is particular
        // to the fermion type. This function is not
        // relevant to fermion types that do not
        // use even/odd checkerboarding. x[i] is the 
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int FsiteOffset(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is the canonical one. X[I] is the
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int ExactFlavors(void);
        // Returns the number of exact flavors of the matrix that
        // is inverted during a molecular dynamics evolution.

    int SpinComponents(void);
        // Returns the number of spin components.

    int FsiteSize(void);
        // Returns the number of fermion field 
        // components (including real/imaginary) on a
        // site of the 4-D lattice.

    int FchkbEvl(void);
	// returns 1 => The fermion fields in the evolution
        //      or the CG that inverts the evolution matrix
	//      are defined on a single checkerboard (half the 
	//      lattice).

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // It does nothing and returns 0.

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // Same as original but with true_res=0;


    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // It does nothing and returns 0.

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // Same as original but with true_res=0;

    int FeigSolv(Vector **f_eigenv, Float *lambda,
		 Float chirality[], int valid_eig[],
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It does nothing and returns 0.

    void SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
			Float mass);
	// It sets the pseudofermion field phi from frm1, frm2.

    void FforceSite(Matrix& force, Vector *frm, 
                            int *x, int mu);
        // It calculates the fermion force per site x
        // and direction mu. frm is the fermion field that 
        // resulted from the application of the inverter on 
        // the pseudofermion field.

    void EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force. 

    Float FhamiltonNode(Vector *phi, Vector *chi);
        // The fermion Hamiltonian of the node sublattice.
        // chi must be the solution of Cg with source phi.	       

    void Fconvert(Vector *f_field,
			  StrOrdType to,
			  StrOrdType from);
        // Convert fermion field f_field from -> to

    Float BhamiltonNode(Vector *boson, Float mass);
        // The boson Hamiltonian of the node sublattice.
};


//------------------------------------------------------------------
//
// FstagTypes is derived from Lattice and is relevant to
// all fermion classes with Staggered type fermions 
// These classes are derived from FstagTypes
//
//------------------------------------------------------------------
class FstagTypes : public virtual Lattice
{
 private:
    char *cname;    // Class name.
    
 protected:

 public:

    FstagTypes(void);

    virtual ~FstagTypes(void);
};


//------------------------------------------------------------------
//
// Fstag is derived from FstagTypes and is relevant to
// staggered fermions.
//
//------------------------------------------------------------------
class Fstag : public virtual FstagTypes
{
 private:
    char *cname;    // Class name.

    int e_vsize;	// even(odd) vector size
    int xv[3];

    Vector *f_tmp;

    void getUDagX(Vector& v, const Vector *cvp, int *x, int mu) const;
    
 public:

    Fstag(void);

    virtual ~Fstag(void);

    FclassType Fclass(void);
        // It returns the type of fermion class

    int FsiteOffsetChkb(const int *x) const
        { return (x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2] ; }
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is not the canonical one but it is particular
        // to the Staggered fermion type. x[i] is the 
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int FsiteOffset(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is the canonical one. X[I] is the
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int ExactFlavors(void);
        // Returns the number of exact flavors of the matrix that
        // is inverted during a molecular dynamics evolution.

    int SpinComponents(void);
        // Returns the number of spin components.

    int FsiteSize(void);
        // Returns the number of fermion field 
        // components (including real/imaginary) on a
        // site of the 4-D lattice.

    int FchkbEvl(void);
	// returns 1 => The fermion fields in the evolution
        //      or the CG that inverts the evolution matrix
	//      are defined on a single checkerboard (half the 
	//      lattice).

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix that appears in the HMC 
        // evolution ([Dirac^dag Dirac]). The inversion is done
	// with the conjugate gradient. cg_arg is the structure
        // that contains all the control parameters, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
	// f_in and f_out are defined on a checkerboard.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
	// The function returns the total number of CG iterations.

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // Same as original but with true_res=0;

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix (Dirac operator). The inversion
	// is done with the conjugate gradient. cg_arg is the 
        // structure that contains all the control parameters, f_in 
        // is the fermion field source vector, f_out should be set 
        // to be the initial guess and on return is the solution.
	// f_in and f_out are defined on the whole lattice.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|.
        // cnv_frm is used to specify if f_in should be converted 
        // from canonical to fermion order and f_out from fermion 
        // to canonical. 
        // prs_f_in is not used. The source f_in is always preserved.
	// The function returns the total number of CG iterations.

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // Same as original but with true_res=0;

    int FeigSolv(Vector **f_eigenv, Float lambda[], 
		 Float chirality[], int valid_eig[],
		 Float **hsum,
		 EigArg *eig_arg,
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It finds the eigenvectors and eigenvalues of A where
        // A is the fermion matrix (Dirac operator). The solution
	// uses Ritz minimization. eig_arg is the 
        // structure that contains all the control parameters, f_eigenv
        // are the fermion field source vectors which should be
        // defined initially, lambda are the eigenvalues returned 
        // on solution. f_eigenv is defined on the whole lattice.
	// The function returns the total number of Ritz iterations.

    void SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
			Float mass);
	// It sets the pseudofermion field phi from frm1, frm2.

    void FforceSite(Matrix& force, Vector *frm, 
                            int *x, int mu);
        // It calculates the fermion force per site x
        // and direction mu. frm is the fermion field that 
        // resulted from the application of the inverter on 
        // the pseudofermion field.

    void EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force.

    Float FhamiltonNode(Vector *phi, Vector *chi);
        // The fermion Hamiltonian of the node sublattice.
        // chi must be the solution of Cg with source phi.	       

    void Fconvert(Vector *f_field,
			  StrOrdType to,
			  StrOrdType from);
        // Convert fermion field f_field from -> to

    Float BhamiltonNode(Vector *boson, Float mass);
        // The boson Hamiltonian of the node sublattice.

    void Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int dir_flag);

};





//------------------------------------------------------------------
//
// Fstag is derived from FstagTypes and is relevant to
// staggered fermions.
//
//------------------------------------------------------------------
class FstagAsqtad : public virtual FstagTypes
{
 private:
    char *cname;    // Class name.

    int e_vsize;	// even(odd) vector size
    int xv[3];

    Vector *f_tmp;

    void getUDagX(Vector& v, const Vector *cvp, int *x, int mu) const;
    
 public:

    FstagAsqtad(void);

    virtual ~FstagAsqtad(void);

    FclassType Fclass(void);
        // It returns the type of fermion class

    int FsiteOffsetChkb(const int *x) const
        { return (x[3]>>1)+xv[0]*x[0]+xv[1]*x[1]+xv[2]*x[2] ; }
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is not the canonical one but it is particular
        // to the Staggered fermion type. x[i] is the 
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int FsiteOffset(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is the canonical one. X[I] is the
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int ExactFlavors(void);
        // Returns the number of exact flavors of the matrix that
        // is inverted during a molecular dynamics evolution.

    int SpinComponents(void);
        // Returns the number of spin components.

    int FsiteSize(void);
        // Returns the number of fermion field 
        // components (including real/imaginary) on a
        // site of the 4-D lattice.

    // let us not worry about the eigensolver
    int FeigSolv(Vector **f_eigenv, Float lambda[], 
		 Float chirality[], int valid_eig[],
		 Float **hsum,
		 EigArg *eig_arg,
		 CnvFrmType cnv_frm = CNV_FRM_YES) ;


    int FchkbEvl(void);
	// returns 1 => The fermion fields in the evolution
        //      or the CG that inverts the evolution matrix
	//      are defined on a single checkerboard (half the 
	//      lattice).

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix that appears in the HMC 
        // evolution ([Dirac^dag Dirac]). The inversion is done
	// with the conjugate gradient. cg_arg is the structure
        // that contains all the control parameters, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
	// f_in and f_out are defined on a checkerboard.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
	// The function returns the total number of CG iterations.

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // Same as original but with true_res=0;

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix (Dirac operator). The inversion
	// is done with the conjugate gradient. cg_arg is the 
        // structure that contains all the control parameters, f_in 
        // is the fermion field source vector, f_out should be set 
        // to be the initial guess and on return is the solution.
	// f_in and f_out are defined on the whole lattice.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|.
        // cnv_frm is used to specify if f_in should be converted 
        // from canonical to fermion order and f_out from fermion 
        // to canonical. 
        // prs_f_in is not used. The source f_in is always preserved.
	// The function returns the total number of CG iterations.

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // Same as original but with true_res=0;

    void SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
			Float mass);
	// It sets the pseudofermion field phi from frm1, frm2.

    void FforceSite(Matrix& force, Vector *frm, 
                            int *x, int mu);
        // It calculates the fermion force per site x
        // and direction mu. frm is the fermion field that 
        // resulted from the application of the inverter on 
        // the pseudofermion field.

    void EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force.

    Float FhamiltonNode(Vector *phi, Vector *chi);
        // The fermion Hamiltonian of the node sublattice.
        // chi must be the solution of Cg with source phi.	       

    void Fconvert(Vector *f_field,
			  StrOrdType to,
			  StrOrdType from);
        // Convert fermion field f_field from -> to

    Float BhamiltonNode(Vector *boson, Float mass);
        // The boson Hamiltonian of the node sublattice.

    void Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int dir_flag);

};


//------------------------------------------------------------------
//
// FwilsonTypes is derived from Lattice and is relevant to
// all fermion classes with Wilson type fermions 
// (e.g Fwilson, Fclover, Fdwf, ...). These classes are derived
// from FwilsonTypes
//
//------------------------------------------------------------------
class FwilsonTypes : public virtual Lattice
{
 private:
    char *cname;    // Class name.
    
 protected:
    void (*sproj_tr[8])(IFloat *f, 
			IFloat *v, 
			IFloat *w, 
			int num_blk, 
			int v_stride,
			int w_stride) ;
    // Array with entries that point to 8 non-member functions.
    // These functions are called as follows:
    // sproj_tr[SprojType mu]();
    // For the various SprojTypes see enum.h
    // These functions return a color matrix in f constructed from
    // the spinors v, w using: 
    // f_(i,j) = Tr_spin[ (1 +/- gamma_mu) v_i w^dag_j ]
    //
    // num_blk is the number of spinors v, w. The routines 
    // accumulate the sum over spinors in f.
    //
    // v_stride and w_stride are the number of Floats between spinors
    // (a stride = 0 means that the spinors are consecutive in memory)

    void (*Sigmaproj_tr[12])(IFloat *f, 
			IFloat *v, 
			IFloat *w, 
			int num_blk, 
			int v_stride,
			int w_stride) ;
    // Array with entries that point to 12 non-member functions.
    // These functions are called as follows:
    // Sigmaproj_tr[SigmaprojType mu_nu]();
    // For the various SigmaprojTypes see enum.h 
    // These functions return a color matrix in f constructed from
    // the spinors v, w using: 
    // f_(i,j) = 1/2 Tr_spin[ Sigma_{mu,nu} v_i w^dag_j ]
    //
    // num_blk is the number of spinors v, w. The routines 
    // accumulate the sum over spinors in f.
    //
    // v_stride and w_stride are the number of Floats between spinors
    // (a stride = 0 means that the spinors are consecutive in memory)

 public:

    FwilsonTypes(void);

    virtual ~FwilsonTypes(void);

    void Gamma5(Vector *v_out, Vector *v_in, int num_sites);
    // v_out = Gamma5 * v_in. Gamme5 is in the chiral basis
    //
    //          [ 1  0  0  0]
    // Gamma5 = [ 0  1  0  0]
    //          [ 0  0 -1  0]
    //          [ 0  0  0 -1]
    //
    // num_sites is the number of sites. It is assumed
    // that each site has 24 components.

};


//------------------------------------------------------------------
//
// Fwilson is derived from FwilsonTypes and is relevant to
// wilson fermions.
//
//------------------------------------------------------------------
class Fwilson : public virtual FwilsonTypes
{
 private:
    char *cname;    // Class name.
    
 public:

    Fwilson(void);

    virtual ~Fwilson(void);

    FclassType Fclass(void);
        // It returns the type of fermion class

    int FsiteOffsetChkb(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is not the canonical one but it is particular
        // to the Wilson fermion type. x[i] is the 
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int FsiteOffset(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is the canonical one. X[I] is the
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int ExactFlavors(void);
        // Returns the number of exact flavors of the matrix that
        // is inverted during a molecular dynamics evolution.

    int SpinComponents(void);
        // Returns the number of spin components.

    int FsiteSize(void);
        // Returns the number of fermion field 
        // components (including real/imaginary) on a
        // site of the 4-D lattice.

    int FchkbEvl(void);
	// returns 1 => The fermion fields in the evolution
        //      or the CG that inverts the evolution matrix
	//      are defined on a single checkerboard (half the 
	//      lattice).

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the preconditioned fermion matrix that appears
        // in the HMC evolution (even/odd  preconditioning 
        // of [Dirac^dag Dirac]. The inversion is done
	// with the conjugate gradient. cg_arg is the structure
        // that contains all the control parameters, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
	// f_in and f_out are defined on a checkerboard.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
	// The function returns the total number of CG iterations.

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // Same as original but with true_res=0;

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix (Dirac operator). The inversion
	// is done with the conjugate gradient. cg_arg is the 
        // structure that contains all the control parameters, f_in 
        // is the fermion field source vector, f_out should be set 
        // to be the initial guess and on return is the solution.
	// f_in and f_out are defined on the whole lattice.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
        // cnv_frm is used to specify if f_in should be converted 
        // from canonical to fermion order and f_out from fermion 
        // to canonical. 
        // prs_f_in is used to specify if the source
        // f_in should be preserved or not. If not the memory usage
        // is less by half the size of a fermion vector.
	// The function returns the total number of CG iterations.

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // Same as original but with true_res=0;

    int FeigSolv(Vector **f_eigenv, Float lambda[], 
		 Float chirality[], int valid_eig[],
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It finds the eigenvectors and eigenvalues of A where
        // A is the fermion matrix (Dirac operator). The solution
	// uses Ritz minimization. eig_arg is the 
        // structure that contains all the control parameters, f_eigenv
        // are the fermion field source vectors which should be
        // defined initially, lambda are the eigenvalues returned 
        // on solution. f_eigenv is defined on the whole lattice.
	// The function returns the total number of Ritz iterations.

    void SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
			Float mass);
	// It sets the pseudofermion field phi from frm1, frm2.

    void EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force. 

    Float FhamiltonNode(Vector *phi, Vector *chi);
        // The fermion Hamiltonian of the node sublattice.
        // chi must be the solution of Cg with source phi.	       

    void Fconvert(Vector *f_field,
			  StrOrdType to,
			  StrOrdType from);
        // Convert fermion field f_field from -> to

    Float BhamiltonNode(Vector *boson, Float mass);
        // The boson Hamiltonian of the node sublattice.

};


//------------------------------------------------------------------
//
// Fclover is derived from FwilsonTypes and is relevant to
// clover Wilson fermions.
//
//------------------------------------------------------------------
class Fclover : public virtual FwilsonTypes
{
 private:
    char *cname;    // Class name.

    void EvolveMomFforceSupp(Matrix *mom, Vector *v1, Vector *v2,
		 Vector *v3, Vector *v4, Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the clover contribution of fermion force 

 public:

    Fclover(void);
        // Among other things the constructor allocates
        // memory for the even/odd checkerpoard clover
        // matrices. aux0_ptr of the base class is set
        // to the pointer of the even checkerboard matrices
        // and aux1_ptr to the odd.

    virtual ~Fclover(void);

    FclassType Fclass(void);
        // It returns the type of fermion class

    int FsiteOffsetChkb(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is not the canonical one but it is particular
        // to the Clover fermion type. x[i] is the 
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int FsiteOffset(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is the canonical one. X[I] is the
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int ExactFlavors(void);
        // Returns the number of exact flavors of the matrix that
        // is inverted during a molecular dynamics evolution.

    int SpinComponents(void);
        // Returns the number of spin components.

    int FsiteSize(void);
        // Returns the number of fermion field 
        // components (including real/imaginary) on a
        // site of the 4-D lattice.

    int FchkbEvl(void);
        // returns 0 => The fermion fields in the evolution
        // are defined on ODD-EVEN checkerboard (whole
        // lattice).

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the preconditioned fermion matrix that appears
        // in the HMC evolution (even/odd  preconditioning 
        // of [Dirac^dag Dirac]. The inversion of the odd checkerboard
        // piece is done with the conjugate gradient algorithm
        // while the inversion of the even checkerboard is done
        // using standard explicit hermitian matrix inversion of the
        // clover matrix. cg_arg is the structure that contains
        // all the control parameters for the CG, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
        // f_in and f_out are defined on the whole latice.
        // If true_res !=0 the value of the true residual of the CG and
        // is returned in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
        // The function returns the total number of CG iterations.

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // Same as original but with true_res=0;

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix (Dirac operator) with no 
        // preconditioning. The preconditioned matrix is inverted 
        // and from the result the non-preconditioned f_out is 
        // calculated . The inversion of the odd checkerboard
        // piece is done with the conjugate gradient algorithm
        // while the inversion of the even checkerboard is done
        // using standard explicit hermitian matrix inversion of the
        // clover matrix. cg_arg is the structure that contains
        // all the control parameters for the CG, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
        // f_in and f_out are defined on the whole latice.
        // If true_res !=0 the value of the true residual of the CG and
        // is returned in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
        // cnv_frm is used to specify if f_in should be converted 
        // from canonical to fermion order and f_out from fermion 
        // to canonical. 
        // prs_f_in is used to specify if the source
        // f_in should be preserved or not. If not the memory usage
        // is less by the size of one fermion vector.
        // The function returns the total number of CG iterations.

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // Same as original but with true_res=0;

    int FeigSolv(Vector **f_eigenv, Float lambda[],
		 Float chirality[], int valid_eig[],
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It finds the eigenvectors and eigenvalues of A where
        // A is the fermion matrix (Dirac operator). The solution
	// uses Ritz minimization. eig_arg is the 
        // structure that contains all the control parameters, f_eigenv
        // are the fermion field source vectors which should be
        // defined initially, lambda are the eigenvalues returned 
        // on solution. f_eigenv is defined on the whole lattice.
	// The function returns the total number of Ritz iterations.

    void SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
			Float mass);
	// It sets the pseudofermion field phi from frm1, frm2.

    void EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force.

    Float FhamiltonNode(Vector *phi, Vector *frm);
        // The fermion Hamiltonian of the node sublattice.
        // frm must be the solution of FmatEvlInv with source phi.


    void Fconvert(Vector *f_field,
			  StrOrdType to,
			  StrOrdType from);
        // Convert fermion field f_field from -> to

    Float BhamiltonNode(Vector *boson, Float mass);
        // The boson Hamiltonian of the node sublattice
};


//------------------------------------------------------------------
//
// Fdwf is derived from FwilsonTypes and is relevant to
// domain wall fermions.
//
//------------------------------------------------------------------
class Fdwf : public virtual FwilsonTypes
{
 private:
    char *cname;    // Class name.
    
 public:

    Fdwf(void);

    virtual ~Fdwf(void);

    FclassType Fclass(void);
        // It returns the type of fermion class

    int FsiteOffsetChkb(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is not the canonical one but it is particular
        // to the Dwf fermion type. x[i] is the 
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int FsiteOffset(const int *x) const;
        // Sets the offsets for the fermion fields on a 
        // checkerboard. The fermion field storage order
        // is the canonical one. X[I] is the
        // ith coordinate where i = {0,1,2,3} = {x,y,z,t}.

    int ExactFlavors(void);
        // Returns the number of exact flavors of the matrix that
        // is inverted during a molecular dynamics evolution.

    int SpinComponents(void);
        // Returns the number of spin components.

    int FsiteSize(void);
        // Returns the number of fermion field 
        // components (including real/imaginary) on a
        // site of the 4-D lattice.

    int FchkbEvl(void);
        // Returns 0 => If no checkerboard is used for the evolution
        //      or the CG that inverts the evolution matrix.

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the preconditioned fermion matrix that appears
        // in the HMC evolution (even/odd preconditioning 
        // of [Dirac^dag Dirac]). The inversion is done
	// with the conjugate gradient. cg_arg is the structure
        // that contains all the control parameters, f_in is the
        // fermion field source vector, f_out should be set to be
        // the initial guess and on return is the solution.
	// f_in and f_out are defined on a checkerboard.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
	// The function returns the total number of CG iterations.

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // Same as original but with true_res=0;

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // It calculates f_out where A * f_out = f_in and
        // A is the fermion matrix (Dirac operator). The inversion
	// is done with the conjugate gradient. cg_arg is the 
        // structure that contains all the control parameters, f_in 
        // is the fermion field source vector, f_out should be set 
        // to be the initial guess and on return is the solution.
	// f_in and f_out are defined on the whole lattice.
        // If true_res !=0 the value of the true residual is returned
        // in true_res.
        // *true_res = |src - MatPcDagMatPc * sol| / |src|
        // cnv_frm is used to specify if f_in should be converted 
        // from canonical to fermion order and f_out from fermion 
        // to canonical. 
        // prs_f_in is used to specify if the source
        // f_in should be preserved or not. If not the memory usage
        // is less by half the size of a fermion vector.
	// The function returns the total number of CG iterations.

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // Same as original but with true_res=0;
	
    void Ffour2five(Vector *five, Vector *four, int s_u, int s_l);
        // It transforms a 4-dimensional fermion field
        // to a 5-dimensional field. The 5d field is zero
        // except for the upper two components (right chirality)
        // at s = s_u which are equal to the ones of the 4d field
        // and the lower two components (left chirality) 
        // at s_l, which are equal to the ones of the 4d field
        // where s is the coordinate in the 5th direction.
        // For spread-out DWF s_u, s_l refer to the global
        // s coordinate i.e. their range is from 
        // 0 to [GJP.Snodes() * GJP.SnodeSites() - 1]

    void Ffive2four(Vector *four, Vector *five, int s_u, int s_l);
        // It transforms a 5-dimensional fermion field
        // to a 4-dimensional field. The 4d field has
        // the upper two components (right chirality) equal to the
        // ones of the 5d field at s = s_u and the lower two 
        // components (left chirality) equal to the
        // ones of the 5d field at s = s_l, where s is the 
        // coordinate in the 5th direction.
        // For spread-out DWF s_u, s_l refer to the global
        // s coordinate i.e. their range is from 
        // 0 to [GJP.Snodes() * GJP.SnodeSites() - 1]
        // The same 4D field is generarted in all s node slices.

    int FeigSolv(Vector **f_eigenv, Float lambda[],
		 Float chirality[], int valid_eig[],
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It finds the eigenvectors and eigenvalues of A where
        // A is the fermion matrix (Dirac operator). The solution
	// uses Ritz minimization. eig_arg is the 
        // structure that contains all the control parameters, f_eigenv
        // are the fermion field source vectors which should be
        // defined initially, lambda are the eigenvalues returned 
        // on solution. f_eigenv is defined on the whole lattice.
	// The function returns the total number of Ritz iterations.

    void SetPhi(Vector *phi, Vector *frm1, Vector *frm2,	       
		Float mass);
	// It sets the pseudofermion field phi from frm1, frm2.
	
    void EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force.

    Float FhamiltonNode(Vector *phi, Vector *chi);
        // The fermion Hamiltonian of the node sublattice.
        // chi must be the solution of Cg with source phi.	       

    void Fconvert(Vector *f_field,
			  StrOrdType to,
			  StrOrdType from);
        // Convert fermion field f_field from -> to

    Float BhamiltonNode(Vector *boson, Float mass);
        // The boson Hamiltonian of the node sublattice

    void Freflex (Vector *out, Vector *in);
       // Reflexion in s operator, needed for the hermitian version 
       // of the dirac operator in the Ritz solver.
};


//------------------------------------------------------------------
//
// The following classes have double inheritance. The virtual base 
// class is Lattice with two types of derived classes. One type
// is relevant to the gauge part and has a name that begins with 
// "G". The other type is relevant to the fermion part and has
// a name that begins with "F". The classes below inherit from one 
// gauge class and from one fermion class. All combinations are 
// present.
//
//------------------------------------------------------------------


//------------------------------------------------------------------
// Trivial gauge action -- no fermions
//------------------------------------------------------------------
class GnoneFnone 
    : public virtual Lattice, 
    public Gnone, 
    public Fnone
{
 private:
    char *cname;    // Class name.

 public:
    GnoneFnone(void);
    virtual ~GnoneFnone(void);
};


//------------------------------------------------------------------
// Trivial gauge action -- staggered fermion action
//------------------------------------------------------------------
class GnoneFstag 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public Gnone, 
    public Fstag
{
 private:
    char *cname;    // Class name.

 public:
    GnoneFstag(void);
    virtual ~GnoneFstag(void);
};


//------------------------------------------------------------------
// Trivial gauge action -- wilson fermion action
//------------------------------------------------------------------
class GnoneFwilson 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gnone, 
    public Fwilson
{
 private:
    char *cname;    // Class name.

 public:
    GnoneFwilson(void);
    virtual ~GnoneFwilson(void);
};


//------------------------------------------------------------------
// Trivial gauge action -- clover Wilson fermion action
//------------------------------------------------------------------
class GnoneFclover 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gnone, 
    public Fclover
{
 private:
    char *cname;    // Class name.

 public:
    GnoneFclover(void);
    virtual ~GnoneFclover(void);
};


//------------------------------------------------------------------
// Trivial gauge action -- domain wall fermion action
//------------------------------------------------------------------
class GnoneFdwf 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gnone, 
    public Fdwf
{
 private:
    char *cname;    // Class name.

 public:
    GnoneFdwf(void);
    virtual ~GnoneFdwf(void);
};


//------------------------------------------------------------------
// Wilson gauge action -- no fermions
//------------------------------------------------------------------
class GwilsonFnone 
    : public virtual Lattice, 
    public Gwilson, 
    public Fnone
{
 private:
    char *cname;    // Class name.

 public:
    GwilsonFnone(void);
    virtual ~GwilsonFnone(void);
};


//------------------------------------------------------------------
// Wilson gauge action -- staggered fermion action
//------------------------------------------------------------------
class GwilsonFstag 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public Gwilson, 
    public Fstag
{
 private:
    char *cname;    // Class name.

 public:
    GwilsonFstag(void);
    virtual ~GwilsonFstag(void);
};


//------------------------------------------------------------------
// Wilson gauge action -- wilson fermion action
//------------------------------------------------------------------
class GwilsonFwilson 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gwilson, 
    public Fwilson
{
 private:
    char *cname;    // Class name.

 public:
    GwilsonFwilson(void);
    virtual ~GwilsonFwilson(void);
};


//------------------------------------------------------------------
// Wilson gauge action -- clover Wilson fermion action
//------------------------------------------------------------------
class GwilsonFclover 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gwilson, 
    public Fclover
{
 private:
    char *cname;    // Class name.

 public:
    GwilsonFclover(void);
    virtual ~GwilsonFclover(void);
};


//------------------------------------------------------------------
// Wilson gauge action -- domain wall fermion action
//------------------------------------------------------------------
class GwilsonFdwf 
    : public virtual Lattice, 
    public virtual FwilsonTypes, 
    public Gwilson, 
    public Fdwf
{
 private:
    char *cname;    // Class name.

 public:
    GwilsonFdwf(void);
    virtual ~GwilsonFdwf(void);
};


//------------------------------------------------------------------
// PowerPlaq gauge action -- no fermions
//------------------------------------------------------------------
class GpowerPlaqFnone 
    : public virtual Lattice, 
    public GpowerPlaq, 
    public Fnone
{
 private:
    char *cname;    // Class name.

 public:
    GpowerPlaqFnone(void);
    virtual ~GpowerPlaqFnone(void);
};


//------------------------------------------------------------------
// PowerPlaq gauge action -- staggered fermion action
//------------------------------------------------------------------
class GpowerPlaqFstag 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public GpowerPlaq, 
    public Fstag
{
 private:
    char *cname;    // Class name.

 public:
    GpowerPlaqFstag(void);
    virtual ~GpowerPlaqFstag(void);
};


//------------------------------------------------------------------
// PowerPlaq gauge action -- powerPlaq fermion action
//------------------------------------------------------------------
class GpowerPlaqFwilson 
    : public virtual Lattice, 
    public GpowerPlaq, 
    public Fwilson
{
 private:
    char *cname;    // Class name.

 public:
    GpowerPlaqFwilson(void);
    virtual ~GpowerPlaqFwilson(void);
};


//------------------------------------------------------------------
// PowerPlaq gauge action -- clover PowerPlaq fermion action
//------------------------------------------------------------------
class GpowerPlaqFclover 
    : public virtual Lattice, 
    public GpowerPlaq, 
    public Fclover
{
 private:
    char *cname;    // Class name.

 public:
    GpowerPlaqFclover(void);
    virtual ~GpowerPlaqFclover(void);
};


//------------------------------------------------------------------
// PowerPlaq gauge action -- domain wall fermion action
//------------------------------------------------------------------
class GpowerPlaqFdwf 
    : public virtual Lattice, 
    public GpowerPlaq, 
    public Fdwf
{
 private:
    char *cname;    // Class name.

 public:
    GpowerPlaqFdwf(void);
    virtual ~GpowerPlaqFdwf(void);
};


//------------------------------------------------------------------
// Improved rectangle gauge action -- no fermions
//------------------------------------------------------------------
class GimprRectFnone
    : public virtual Lattice,
    public GimprRect,
    public Fnone
{
 private:
    char *cname;    // Class name.

 public:
    GimprRectFnone(void);
    virtual ~GimprRectFnone(void);
};


//------------------------------------------------------------------
// Improved rectangle gauge action -- staggered fermion action
//------------------------------------------------------------------
class GimprRectFstag
    : public virtual Lattice,
    public virtual FstagTypes,
    public GimprRect,
    public Fstag
{
 private:
    char *cname;    // Class name.

 public:
    GimprRectFstag(void);
    virtual ~GimprRectFstag(void);
};


//------------------------------------------------------------------
// Improved rectangle gauge action -- wilson fermion action
//------------------------------------------------------------------
class GimprRectFwilson
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprRect,
    public Fwilson
{
 private:
    char *cname;    // Class name.

 public:
    GimprRectFwilson(void);
    virtual ~GimprRectFwilson(void);
};


//------------------------------------------------------------------
// Improved rectangle gauge action -- clover Wilson fermion action
//------------------------------------------------------------------
class GimprRectFclover
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprRect,
    public Fclover
{
 private:
    char *cname;    // Class name.

 public:
    GimprRectFclover(void);
    virtual ~GimprRectFclover(void);
};


//------------------------------------------------------------------
// Improved rectangle gauge action -- domain wall fermion action
//------------------------------------------------------------------
class GimprRectFdwf
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprRect,
    public Fdwf
{
 private:
    char *cname;    // Class name.

 public:
    GimprRectFdwf(void);
    virtual ~GimprRectFdwf(void);
};


//------------------------------------------------------------------
// PowerRect gauge action -- no fermions
//------------------------------------------------------------------
class GpowerRectFnone 
    : public virtual Lattice, 
    public GpowerRect, 
    public Fnone
{
 private:
    char *cname;    // Class name.

 public:
    GpowerRectFnone(void);
    virtual ~GpowerRectFnone(void);
};


//------------------------------------------------------------------
// PowerRect gauge action -- staggered fermion action
//------------------------------------------------------------------
class GpowerRectFstag 
    : public virtual Lattice, 
    public virtual FstagTypes, 
    public GpowerRect, 
    public Fstag
{
 private:
    char *cname;    // Class name.

 public:
    GpowerRectFstag(void);
    virtual ~GpowerRectFstag(void);
};


//------------------------------------------------------------------
// PowerRect gauge action -- powerRect fermion action
//------------------------------------------------------------------
class GpowerRectFwilson 
    : public virtual Lattice, 
    public GpowerRect, 
    public Fwilson
{
 private:
    char *cname;    // Class name.

 public:
    GpowerRectFwilson(void);
    virtual ~GpowerRectFwilson(void);
};


//------------------------------------------------------------------
// PowerRect gauge action -- clover PowerRect fermion action
//------------------------------------------------------------------
class GpowerRectFclover 
    : public virtual Lattice, 
    public GpowerRect, 
    public Fclover
{
 private:
    char *cname;    // Class name.

 public:
    GpowerRectFclover(void);
    virtual ~GpowerRectFclover(void);
};


//------------------------------------------------------------------
// PowerRect gauge action -- domain wall fermion action
//------------------------------------------------------------------
class GpowerRectFdwf 
    : public virtual Lattice, 
    public GpowerRect, 
    public Fdwf
{
 private:
    char *cname;    // Class name.

 public:
    GpowerRectFdwf(void);
    virtual ~GpowerRectFdwf(void);
};

//------------------------------------------------------------------
// One Loop Symanzik improved gauge action -- no fermions
//------------------------------------------------------------------
class GimprOLSymFnone
    : public virtual Lattice,
    public GimprOLSym,
    public Fnone
{
 private:
    char *cname;    // Class name.

 public:
    GimprOLSymFnone(void);
    virtual ~GimprOLSymFnone(void);
};


//------------------------------------------------------------------
// One Loop Symanzik improved gauge action -- staggered fermion action
//------------------------------------------------------------------
class GimprOLSymFstag
    : public virtual Lattice,
    public virtual FstagTypes,
    public GimprOLSym,
    public Fstag
{
 private:
    char *cname;    // Class name.

 public:
    GimprOLSymFstag(void);
    virtual ~GimprOLSymFstag(void);
};





//------------------------------------------------------------------
// One Loop Symanzik improved gauge action -- Asqtad staggered action
//------------------------------------------------------------------
class GimprOLSymFstagAsqtad
    : public virtual Lattice,
    public virtual FstagTypes,
    public GimprOLSym,
    public FstagAsqtad
{
 private:
    char *cname;    // Class name.

 public:
    GimprOLSymFstagAsqtad(void);
    virtual ~GimprOLSymFstagAsqtad(void);
};


//------------------------------------------------------------------
// One Loop Symanzik improved gauge action -- wilson fermion action
//------------------------------------------------------------------
class GimprOLSymFwilson
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprOLSym,
    public Fwilson
{
 private:
    char *cname;    // Class name.

 public:
    GimprOLSymFwilson(void);
    virtual ~GimprOLSymFwilson(void);
};


//------------------------------------------------------------------
// One Loop Symanzik improved gauge action -- clover Wilson fermion action
//------------------------------------------------------------------
class GimprOLSymFclover
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprOLSym,
    public Fclover
{
 private:
    char *cname;    // Class name.

 public:
    GimprOLSymFclover(void);
    virtual ~GimprOLSymFclover(void);
};


//------------------------------------------------------------------
// One Loop Symanzik improved gauge action -- domain wall fermion action
//------------------------------------------------------------------
class GimprOLSymFdwf
    : public virtual Lattice,
    public virtual FwilsonTypes,
    public GimprOLSym,
    public Fdwf
{
 private:
    char *cname;    // Class name.

 public:
    GimprOLSymFdwf(void);
    virtual ~GimprOLSymFdwf(void);
};



#endif
CPS_END_NAMESPACE
