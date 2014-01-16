#ifndef INCLUDED_LATTICE_H
#define INCLUDED_LATTICE_H

#include<config.h>
#include<math.h>

//----------------------------------------------------------------------
/*!\file
  \brief  Definitions of the Lattice classes.

  $Id: lattice.h,v 1.69 2013-06-25 12:51:12 chulwoo Exp $
*/
/*----------------------------------------------------------------------
  $Author: chulwoo $
  $Date: 2013-06-25 12:51:12 $
  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/lattice.h,v 1.69 2013-06-25 12:51:12 chulwoo Exp $
  $Id: lattice.h,v 1.69 2013-06-25 12:51:12 chulwoo Exp $
  $Name: not supported by cvs2svn $
  $Revision: 1.69 $
  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/lattice.h,v $
  $State: Exp $
*/  
//------------------------------------------------------------------


#include <util/gjp.h>
#include <util/enum.h>
#include <util/random.h>
#include <util/vector.h>
#include <util/lat_data.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/data_types.h>
#include <comms/glb.h>

#include <alg/cg_arg.h>
#include <alg/eig_arg.h>
#include <alg/lanczos_arg.h>
#include <alg/force_arg.h>
#ifdef PARALLEL
#include <comms/sysfunc_cps.h>
#endif

CPS_START_NAMESPACE

class LinkBuffer;

//------------------------------------------------------------------
//
//! The fundamental abstract base class.
/*! This is the basic class of the CPS from which many are derived.
  The derived classes define the lattice action, and many of the methods
  implementing operations with the action are defined or declared here.
  This class holds the gauge configuration and implements operations
  on the gauge configuration.
  all sorts of other things are actually  defined here too, for instance,
  many of the methods used in the HMD algorithms, operations on spinor
  fields, \e etc.  
*/
//------------------------------------------------------------------
class Lattice
{

 private:

    const char *cname;    // Class name.

    static Matrix* gauge_field;
       // Pointer to the gauge field configuration.
    static Float* u1_gauge_field;
       // Pointer to the u1 gauge field configuration.
  
    static int is_allocated;	
       // 0 = gauge field has not been allocated
       // 1 = gauge field has been allocated

    static int is_initialized;	
       // 0 = gauge field has not been initialized
       // 1 = gauge field has been initialized
    static int u1_is_initialized;

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

    static Float fix_gauge_stpCond;
      // The last stopping condition used.

    static int scope_lock;

 protected:

    int smeared;
    
    static int node_sites[5];
    //!< The local lattice dimensions.
    /*!<
      A reimplementation of GlobalJobParameter::XnodeSite(), \e etc:
-   	 node_sites[0] = GJP.XnodeSite();
-    	 node_sites[1] = GJP.YnodeSite();
-    	 node_sites[2] = GJP.ZnodeSite();
-    	 node_sites[3] = GJP.TnodeSite();
-    	 node_sites[4] = GJP.SnodeSite();
    */

    static int g_dir_offset[4];
    //!< Offsets to help find the array index of gauge links.
    /*!<
      Specifically, \a g_dir_offset[i] is the the internal array index of the
      first of the four gauge field links at the lattice site where the \a i th
      coordinate is 1 and all the rest are 0. The canonical order is assumed,
      so \a i = 0, 1, 2 and 3 corresponds to the X , Y, Z and T direections
      respectively.
    */

    void *f_dirac_op_init_ptr;
      // A pointer that is used by the fermion classes to
      // point to data that need only be initialized by the fermion
      // class constructor and are needed by the relevant dirac
      // operator.

    void *aux0_ptr;
      //!< A  pointer!

    void *aux1_ptr;
      //!< Another pointer!

    friend class LinkBuffer;
    
    // Added in by Ping for anisotropic lattices
    //------------------------------------------------------------------
    void MltFloatImpl(Float factor, int dir);
    //!< Multiplies all gauge links in direction \a dir by a real factor.

    LinkBuffer * link_buffer ;
    //!< The array of off-node links, accessed by methods in link_buffer.C



    

 public:

    const Matrix * GetLink(const int *x, int mu) const;
      //!< Gets the gauge link U_mu(x).
      // defined relative to the local site (0,0,0,0).
      // If the link is on-node, it returns a reference to
      // to the link.  If the link is off-node, it retrieves the link
      // into a static buffer and returns the reference to the buffer.
      // Since the buffer can be used by other routines as well as other
      // calls to this routine, as a general rule, the link should be
      // used immediately or else copied.

    const Matrix * GetLinkOld(Matrix *g_offset, const int *x,
             int dir, int mu) const;
    //!< Gets the gauge link U_mu(x+dir)
      // GRF: renamed to avoid conflict with the more general
      // purpose function
        // U1 field
    const Float* GetLinkOld(Float *g_offset, const int *x, int v, int mu) const;

    int LinkBufferIsEnabled() {
      // Line below casts a pointer to an integer to check NULL ness
      // this is naughty (bj)
      // return ((int) link_buffer)

      // This is more explicit 
      if ( link_buffer == 0x0 ) { 
	return 0;
      }
      else {
	return 1;
      }
    }
    //!< Returns true if there is a buffer for the links, false otherwise.

    int EnableLinkBuffer(int buf_sz);
      //!< Creates a link buffer.

    void DisableLinkBuffer();   
      //!< delete the LinkBuffer Object when it's not in use. 

    const Matrix * GetBufferedLink(const int *x, int mu);
      //!< Gets a link from the buffer.

    void ClearBufferedLink(const int * x, int mu);
     //!< Removes links from the buffer.

    void ClearAllBufferedLink();
      //!< Deletes all the links from the buffer.
      //this must be called when lat.Unitarize() is used

    int IsOnNode(const int * x);
      //!< Checks if a lattice site local to this node.

    void PathOrdProdPlus(Matrix & mat, const int *x, const int * dirs, int n);
    //!< Computes the product of links along a path and adds it to a matrix.
      //given the starting point x, the directions of each step on the path
      //and the number of steps. calculate the path ordered product of 
      //all the links along the path and add the result to mat 
      //each direction could be {0,1,2,3,4,5,6,7} which coresponds to 
      //the directions {n_x, n_y, n_z, n_t, -n_x, -n_y, -n_z, -n_t}
      //the result is returned in mat.

    void PathOrdProd(Matrix & mat, const int *x, const int * dirs, int n);
    //!< Computes the product of links along a path.
      //also calculates the path ordered product, but the result returned 
      //is that product of the path that end on the local node

 public: 

//------------------------------------------------------------------
//------------------------------------------------------------------
// BASE CLASS FUNCTIONS
//------------------------------------------------------------------
//------------------------------------------------------------------

    Lattice();

    virtual ~Lattice();

    Matrix *GaugeField() const {
        return gauge_field;
    }
    //!< Returns the pointer to the gauge field configuration.
    Float *U1GaugeField() const;
         //!< Returns the pointer to the gauge field configuration.
    void twist_bc(int sign);
      // implement twisted boundary conditions
    void mult_su3_links_by_u1(const Float Q);
      //mult su3 links by exp(i Q theta_mu) theta_mu = QED potential
    void mult_su3_links_by_mag_field(Float* B0);
      //mult su3 links by constant B field (linear vector potential)
      // and fix boundary links to give same plaquette (flux)


    void GaugeField(Matrix *u);
    //!< Copies an array into the gauge configuration.
    void U1GaugeField(Float *u);
        //!< Copies an array into the gauge configuration.

    int GsiteOffset(const int *x) const
        { return x[0]*g_dir_offset[0]+x[1]*g_dir_offset[1]
	+x[2]*g_dir_offset[2]+x[3]*g_dir_offset[3];  }

	//!< Gets the array index of a gauge link.
	/*!<
	Specifically, the internal array index of the first of the four gauge
	field links at lattice site x for the canonical storage order.
	\param x The lattice site coordinates.
	\return The array index.
	*/
    virtual unsigned long GsiteOffset(const int *x, const int dir) const;


    void CopyGaugeField(Matrix* u);
        //!< Copies the  gauge configuration into an array;
    void CopyU1GaugeField(Float* u);
        //!< Copies the  gauge configuration into an array;

    int CompareGaugeField(Matrix* u);
        //!< Compares the gauge configuration to the lattice instance;

    StrOrdType StrOrd() const ;
        //!< Returns the storage order.

    int Colors() const;
        //!< Returns the number of colors.	  

    int GsiteSize();
        //!< Gets the number of gauge field  components per lattice site.
    int U1GsiteSize();
        //!< Gets the number of gauge field  components per lattice site.

    void Staple(Matrix& stap, int *x, int mu);
        //!< Calculates the gauge field square staple sum around a link
        // GRF: consider changing the function name
        // to Lattice::PlaqStaple() for consistency.

    void BufferedStaple(Matrix & stap, const int *x, int mu);
        //!< Calculates the gauge field square staple sum  around a link
        //Buffered version of staple

    void RectStaple(Matrix& stap, int *x, int mu) ;
        //!< Calculates the rectangle staple sum around a link.
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
        //!< Calculates the rectangle staple sum around a link.
        //buffered version of RectStaple 

	//! Appears not to be implemented.
    void RectStaple1(Matrix& stap, int *x, int mu);
        // it calculates the flat 6-1 link staple using the PathOrdProdPlus routine
	// \sum_{+/-v; v!= u }{
	//    U_v(x+u) U_v(x+u+v) U_{-u}(x+u+v+v) U_{-v}(x+v+v) U_{-v}(x+v)
	// +  U_v(x+u) U_{-u}(x+v+u) U_{-u}(x+v)  U_{-v}(x+v-u) U_u(x-u)
	// +  U_u(x+u) U_v(x+u+u) U_{-u}(x+u+u+v) U_{-u}(x+u+v) U_{-v}(x+v)

    //! Appears not to be implemented.
    void ChairStaple(Matrix& stap, int *x, int mu);
         //it calculates the chair shaped 6-1 link staple at x, mu
         // \sum_{w != +/-u, w !=+/-v, w != +/-u}{
         //    U_w(x+u) U_v(x+u+w) U_{-w}(x+u+v+w) U_{-u}(x+u+v) U_{-v}(x+v)  
         //  + U_w(x+u) U_v(x+u+w) U_{-u}(x+u+v+w) U_{-v}(x+v+w) U_{-w}(x+w)
         //  + U_w(x+u) U_{-u}(x+u+w) U_v(x+w) U_{-w}(x+w+v) U_{-v}(x+v)
         // }

    void BufferedChairStaple(Matrix &stap, const int *x, int mu);
        //!< Calculates the chair shaped staple sum around a link.
    
    //! Appears not to be implemented.
    void CubeStaple(Matrix &stap, const int *x, int mu);
        //it calculates the cube shaped 6-1 link staple at x, mu
        // \sum_{w != +/-u, w !=+/-v, w != +/-u}{
        //   U_v(x+u) U_w(x+u+v) U_{-u}(x+u+v+w) U_{-v}(x+v+w) U_{-w}(x+w)
        // }
            
    void BufferedCubeStaple(Matrix &stap, const int *x, int mu);  
        //!< Calculates the 5-link cube shaped staple sum around a link
        //buffered version of CubeStaple

    virtual void AllStaple(Matrix &stap, const int *x, int mu)=0;
    //!< Computes all of the staple sums around a link.
        //pure virtual function must be implemented for the gauge actions 
        //derived from it.
        //given a link calculate all of its staples, depending on the 
        //lattice type.
        //this is used in heatbath.

    void Plaq(Matrix &plaq, int *x, int mu, int nu) const;
    //!< Computes a plaquette.
        // Added by Ping for the purpose of debugging now, may be more useful
        // later on. It calculates  the plaquette 
        //   U_u(x) U_v(x+u) U_u(x+v)~ U_v(x)~

    Float ReTrPlaq(int *x, int mu, int nu) const;
        //!< Calculates the real part of the trace of a plaquette.
        // field at site x, mu, nu with mu < nu.
        // The plaquette field is:
        //
        //   U_u(x) U_v(x+u) U_u(x+v)~ U_v(x)~
    Float ReU1Plaq(int *x, int mu, int nu) const;

    Float SumReTrPlaqNode() const;
       //!< Calculates the local sum of the real part of the trace of the plaquette 
    Float SumReU1PlaqNode() const;

    Float SumReTrPlaq() const;
       //!< Calculates the global sum of the real part of the trace of the plaquette 
    Float SumReU1Plaq(void) const;

    Float ReTrRect(int *x, int mu, int nu) const;
        //!< Calculates the real part of the trace of a 6-link rectangle.
       // It calculates the real part of the trace of the rectangle
       // field at site x, in the (mu, nu) plane with the long axis
       // of the rectangle in the mu direction.
       // The rectangle field is:
       //
       //   U_u(x) U_u(x+u) U_v(x+2u) U_u(x+u+v)~ U_u(x+v)~ U_v(x)~
       //

    Float SumReTrRectNode() const;
        //!< Calculates the local sum of the real part of the trace of the 6-link rectangle.
       // It calculates the sum of the real part of the trace of the
       // rectangle field at each site of the node sublattice.

    Float SumReTrRect() const;
        //!< Calculates the global sum of the real part of the trace of the 6-link rectangle.
       // It calculates the sum of the real part of the trace of the
       // rectangle field at each site of the whole lattice.
    
    Float ReTrLoop(const int *x, const int *dir,  int length) ;
    //!< Computes the real trace of the product of links along a path.

    Float SumReTrCubeNode() ;
    //!< Calculates the local sum of the real part of the trace of the cube

    Float SumReTrCube() ;
    //!< Calculates the global sum of the real part of the trace of the cube

  // Added in by Ping for anisotropic lattices
  //------------------------------------------------------------------
  Float AveReTrPlaqNodeNoXi() const;
  //!< Calculates the local average of the real part of the trace of the plaquettes perpendicular to the special anisotropic direction.
  // Normalization:  1 for ordered links

  Float AveReTrPlaqNodeXi() const;
  //!< Calculates the local average of the real part of the trace of the plaquettes parallel to the special anisotropic direction.
  // Normalization:  1 for ordered links

  Float AveReTrPlaqNoXi() const;
  //!< Calculates the global average of the real part of the trace of the plaquettes perpendicular to the special anisotropic direction.  
  // Normalization:  1 for ordered links
  // Average over plaq's perpendicular to the special anisotropic dir.

  Float AveReTrPlaqXi() const;
  //!< Calculates the global average of the real part of the trace of the plaquettes parallel to the special anisotropic direction.
  // Normalization:  1 for ordered links
  // Average over plaq's parallel to the special anisotropic dir.

  // Added in by Schmidt for anisotropic lattices / Thermo
  //------------------------------------------------------------------
  Float AveReTrRectNodeNoXi() const;
  //!< Calculates the local average of the real part of the trace of the rectangle perpendicular to the special anisotropic direction.
  // Normalization:  1 for ordered links

  Float AveReTrRectNodeXi1() const;
  //!< Calculates the local average of the real part of the trace of the rectangle, where the short axix is the special anisotropic direction.
  // Normalization:  1 for ordered links

  Float AveReTrRectNodeXi2() const;
  //!< Calculates the local average of the real part of the trace of the rectangle, where the long axis is the special anisotropic direction.
  // Normalization:  1 for ordered links

  Float AveReTrRectNoXi() const;
  //!< Calculates the global average of the real part of the trace of the rectangle perpendicular to the special anisotropic direction.  
  // Normalization:  1 for ordered links
  // Average over rect's perpendicular to the special anisotropic dir.

  Float AveReTrRectXi1() const;
  //!< Calculates the global average of the real part of the trace of the rectangle, where the short axis is the special anisotropic direction.
  // Normalization:  1 for ordered links
  // Average over rect's with short axis in the special anisotropic dir.

  Float AveReTrRectXi2() const;
  //!< Calculates the global average of the real part of the trace of the rectangle, where the long axis is the special anisotropic direction.
  // Normalization:  1 for ordered links
  // Average over rect's with long axis in the special anisotropic dir.

  void MltFloat(Float factor, int dir)      {
    if (factor != 1.0)    MltFloatImpl(factor, dir);
  }    
  //!< Multiplies all gauge links with direction \a dir by a real factor.
  /*!<
    \param factor The real scale factor.
    \param dir The direction index of the links to be scaled; 
    \a dir = 0, 1, 2 or 3 for direction X, Y, Z and T respectively (all kinds
    of storage order are handled correctly).
    \post The  gauge field  links in direction \a dir are scaled.
  */
 
    void Reunitarize();
    //!< Re-unitarize the gauge field configuration.

    void Reunitarize(Float &dev, Float &max_diff);
    	//!< Test the gauge field for unitarity violation and reunitarize it.
	// and return:
	// dev = sqrt( Sum_i [ (U(i) - V(i))^2 ] / (Vol*4*18) ),
        // max_diff = Max_i[ |U(i) - V(i)| ]
        // where U(i), V(i) is the gauge field before and after 
        // reunitarization. The index i runs over all components of
        // the gauge field.

    int MetropolisAccept(Float delta_h, Float *accept);
    //!< Metropolis algorithm decision.
        // 0 reject, 1 accept. If delta_h < 0 it accepts unconditionally.

    int MetropolisAccept(Float delta_h);
    //!< Metropolis algorithm decision.
        // 0 reject, 1 accept. If delta_h < 0 it accepts unconditionally.

    void EvolveGfield(Matrix *mom, Float step_size);
        //!< Molecular dynamics evolution of the gauge field.

    Float MomHamiltonNode(Matrix *momentum);
    //!< The kinetic energy term of the canonical Hamiltonian on the local lattice.

    void Convert(StrOrdType new_str_ord,
		 Vector *f_field_1,
		 Vector *f_field_2);
    //!< Converts the gauge field and two fermion fields to a new data layout.
    void Convert(StrOrdType new_str_ord,
		 Vector *f_field_1);
    //!< Converts a fermion field to a new data layout.
    

    void Convert(StrOrdType new_str_ord);
    //!< Converts the gauge field to a new data layout.    

    //! A random gaussian anti-Hermitian matrix field.
    void RandGaussAntiHermMatrix(Matrix *mat, Float sigma2);
        // It produces an anti-Hermitian matrix for each
        // site of the lattice with random
        // entries weighted according to 
	// exp(- Tr(mat^2)) / (2 * sigma2)

    //! Creates a random gaussian spin-colour field.
    void RandGaussVector(Vector *vect, Float sigma);

    //! Creates a random gaussian spin-colour field.
    void RandGaussVector(Vector *vect, Float sigma,
                        FermionFieldDimension frm_field_dim);
        // This version assumes the fermion field spans the whole lattice

    //! Creates a random gaussian spin-colour field.
    void RandGaussVector(Vector *vect, Float sigma, int num_chckbds,
	FermionFieldDimension frm_field_dim = FIVE_D );
    void RandGaussVector(Vector *vect, Float sigma, int num_chckbds,
	StrOrdType str, FermionFieldDimension frm_field_dim = FIVE_D );
        // It produces 3-vectors in canonical storage order
        // for each site of the lattice with random
        // entries weighted according to
        // exp(- Tr(mat^2)) / (2 * sigma2)
        // It defaults to producing fermionic fields over the entire lattice,
        // but can also produce even and odd checkerboards.
        // If keep_snodes == 1, then we divide the site size by SnodeSites()

    void SetGfieldOrd();
    //!< Creates a unit gauge field.
    void SetU1GfieldOrd();
    //!< Creates a unit U1 gauge field.

    void SetGfieldDisOrd();
    //!< Creates a random (disordered) gauge field.

    int GupdCnt();
    //!< Reads the gauge field updates counter.

    int GupdCnt(int set_val);
    //!< Sets the gauge field updates counter.

    int GupdCntInc(int inc_val = 1);
    //!< Increments the gauge field updates counter.

    Float MdTime();
      // Returns the value of md_time.

    Float MdTime(Float set_val);
      // Sets the value of md_time to set_val and returns that value.

    Float MdTimeInc(Float inc_val = 0.5);
      // Increments md_time by inc_val and returns new value.

    void *FdiracOpInitPtr()
      { return f_dirac_op_init_ptr; }
      // Returns a pointer that is used by the fermion classes to
      // point to data that need only be initialized by the fermion
      // class constructor and are needed by the relevant dirac
      // operator.

    void FixGaugeAllocate(FixGaugeType GaugeType,int NHplanes=0,int *Hplanes=0);
        //!< Allocates memory for the gauge fixing matrices.

    int FixGauge(Float StopCond, int MaxIterNum);
    //!< Fixes the gauge.

        // FixGaugeAllocate must be called first.
        //
        // Float StopCond is the stopping condition;
        //
        // int MaxIterNum - issues a warning if reached

    void FixGaugeFree();
        //!< Free memory used by the gauge fixing matrices.

    Matrix **FixGaugePtr();
      //!< Returns (a pointer to the first element of) an array of pointers to gauge-fixed hyperplanes.

    FixGaugeType FixGaugeKind();
      //!< Returns the kind of gauge fixing.

    Float FixGaugeStopCond();
      //!< Returns the stopping condition used 


    void *Aux0Ptr();
     //!< Returns a general purpose auxiliary pointer.

    void *Aux1Ptr();
      //!< Returns a general purpose auxiliary pointer.

    void GsoCheck();
    //!< Checks that the gauge field is identical on 5th dimension local lattice slices.

      // If GJP.Snodes() == 1 it just returns.
      // If GJP.Snodes() != 1 it checks that the "spread-out"
      // gauge field is identical along all s-slices by comparing 
      // the checksum and plaquette value. This situation arises for
      // DWF with the s direction spread out across many processors.
      // If either the checksum or plaquette are not identical it
      // exits with an error.

    void SoCheck(Float num);
    //!< Checks that a number is identical on 5th dimension local lattice slices.


//------------------------------------------------------------------
//------------------------------------------------------------------
// GAUGE ACTION RELATED VIRTUAL FUNCTIONS.
//------------------------------------------------------------------
//------------------------------------------------------------------


//------------------------------------------------------------------
//------------------------------------------------------------------
// FERMION ACTION RELATED VIRTAUL FUNCTIONS.
//------------------------------------------------------------------
//------------------------------------------------------------------

    //! Not implemented here.
    virtual void Gamma5(Vector *v_out, Vector *v_in, int num_sites);
    
    //! Not implemented here.
    virtual void Ffour2five(Vector *five, Vector *four, 
			    int s_r, int s_l, int Ncb =2);
    // This is "truly" implemented only in the Fdwf derived class

    //! Not implemented here.
    virtual void Ffive2four(Vector *four, Vector *five, 
			    int s_r, int s_l, int Ncb =2);
        // This is "truly" implemented only in the Fdwf derived class

    //! Not implemented here.
    virtual void Freflex (Vector *out, Vector *in) {};
       // This is "truly" implemented only in the Fdwf derived class

    //! Not implemented here.
    virtual void Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int dir_flag);

    // dir_flag is flag which takes value 0 when all direction contribute to D
    // 1 - when only the special anisotropic direction contributes to D,
    // 2 - when all  except the special anisotropic direction.
    // Currently this function is implemented only in the Fstag class
    //! Not implemented here.

    //! Not implemented here.
    virtual void FdMdmu(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int order);
    // order is the order of the derivative with respect to mu
    // Currently this function is implemented only in the Fstag class

   //~~ Added by DHR for twisted-mass fermions
   //~~------------------------------------------------------------------
   //~~ For most functions, twisted mass parameter passed in cg_arg
   //~~ For functions not using cg_arg, added explicit twisted-mass
   //~~ parameter to function call parameters
   //~~ These functions implemented in lattice_base as an error, and
   //~~ "truly" implemented only in the F_wilsonTm derived class

    virtual Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
    		Float mass, Float epsilon, DagType dag);

    virtual ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
			 Float mass, Float epsilon, Float step_size);
    //~~ Molecular dynamics evolution due to the fermion force

    virtual ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
			 Float mass, Float epsilon, Float step_size);
    //~~ Molecular dynamics evolution due to the boson part of quotient integrator

//------------------------------------------------------------------
//------------------------------------------------------------------
// Bosonic action related virtual functions.
//------------------------------------------------------------------
//------------------------------------------------------------------

   //~~ Added by DHR for twisted-mass fermions
   //~~------------------------------------------------------------------
   //~~ For most functions, twisted mass parameter passed in cg_arg
   //~~ For functions not using cg_arg, added explicit twisted-mass
   //~~ parameter to function call parameters
   //~~ These functions implemented in lattice_base as an error, and
   //~~ "truly" implemented only in the F_wilsonTm derived class

    virtual Float BhamiltonNode(Vector *boson, Float mass, Float epsilon);


//------------------------------------------------------------------
//------------------------------------------------------------------
// GAUGE ACTION RELATED PURE VIRTUAL FUNCTIONS
//------------------------------------------------------------------
//------------------------------------------------------------------

    virtual GclassType Gclass() = 0;
        //!< Returns the type of gauge action

    virtual void GactionGradient(Matrix &grad, int *x, int mu) = 0;
    //!< Calculates the partial derivative of the gauge action w.r.t. the link U_mu(x).
    /*!<
      \param grad The computed gradient.
      \param x The coordinates of the lattice site.
      \param mu The direction of the link.
     */

    static uint64_t ForceFlops;
    //!< Counter for flops in the HMD force calculations.

    virtual ForceArg EvolveMomGforce(Matrix *mom, Float step_size) = 0;
    //!< Molecular dynamics evolution of the conjugate momentum
    /*!<
      The momentum is evolved for a single molecular dynamics timestep
      using the force from the pure gauge action.
      \param mom The momentum matrices on all links.
      \param step_size The molecular dynamics timestep used in the numerical
      integration.
      \post \a mom is assigned the value of the momentum after the molecular
      dynamics evolution.
    */

    virtual Float GhamiltonNode() = 0;
    //!< Computes the pure gauge action on the local sublattice.
    /*!<
      \return The value of the pure gauge action on this node.
    */

//------------------------------------------------------------------
//------------------------------------------------------------------
// Fermion action related pure virtual functions 
//------------------------------------------------------------------
//------------------------------------------------------------------

    virtual FclassType Fclass() const = 0;
        //!< Returns the type of fermion action.
    int FstagType(){
      if (Fclass()==F_CLASS_STAG || Fclass() ==F_CLASS_ASQTAD || 
        Fclass() ==F_CLASS_P4) return 1;
      else return 0;
    }

    //~~ added F_CLASS_WILSON_TM for twisted mass fermions
    int FwilsonType(){
      if (Fclass()==F_CLASS_WILSON || Fclass() ==F_CLASS_CLOVER || 
	  Fclass() ==F_CLASS_DWF || Fclass()==F_CLASS_MOBIUS 
        || Fclass() ==F_CLASS_WILSON_TM || Fclass() ==F_CLASS_NAIVE
        || Fclass() ==F_CLASS_BFM ) return 1;
      else return 0;
    }

    //~~ to distinguish 5D types. Currently exclude BFM, as BFM does all the 5D stuff outside CPS.
    int F5D(){
      if ( Fclass() ==F_CLASS_DWF || Fclass()==F_CLASS_MOBIUS 
        || Fclass() ==F_CLASS_MDWF || Fclass() == F_CLASS_BFM ) return 1;
      else return 0;
    }

    virtual int FsiteOffsetChkb(const int *x) const = 0;
    //!< Gets the lattice site index for the odd-even (checkerboard) order.
    /*!<
      When a field is stored in an odd-even (checkerboard) order,
      this method converts a site's
      cartesian coordinates into its lattice site index.
      \param x The cartesian lattice site coordinates.
      \return The lattice site index.
    */

    virtual int FsiteOffset(const int *x) const = 0;
    //!< Gets the lattice site index for the canonical order.
    /*!<
      When the fermion field is stored in canonical order,
      defined by StrOrdType = CANONICAL, this method converts a sites
      cartesian coordinates into its lattice site index.
      \param x The cartesian lattice site coordinates.
      \return The lattice site index.
    */

    virtual int ExactFlavors() const = 0;
    //!<  The number of dynamical flavors.
    /*!<
      \return The number of flavours defined  when this action is used in
      molecular dynamics dynamical fermion algorithms.
    */

    virtual int SpinComponents() const = 0;
        //!< Returns the number of spin components.

    virtual int FsiteSize() const = 0;
    //!< Gets the size of a fermion field per 4-dim. lattice site.
    /*!<
      \return The total number of all fermionic field components, \e i.e.
      the number of floating point numbers, at each site of the 4-dim. lattice.
    */
    
    virtual int FchkbEvl() const = 0;
    //!< Determines whether one or both parities are used in the molecular dynamics evolution.
    /*!<
      Are the fields used in the molecular dynamics algorithms defined on
      the whole lattice or just on sites on one parity?
      \return 0 if both parities are used, 1 if only one parity is used.
    */

    virtual int FmatEvlInv(Vector *f_out, Vector *f_in, 
			   CgArg *cg_arg, 
                           Float *true_res,
			   CnvFrmType cnv_frm = CNV_FRM_YES) = 0;
    //!< The matrix inversion used in the molecular dynamics algorithms.
    /*!<
      Solves \f$ M^\dagger M f_{out} = f_{in} \f$
      for \f$ f_{out} \f$, where \a M is the 
      (possibly odd-even preconditioned) fermionic matrix.

      \param f_out The initial guess of solution vector.
      \param f_in The source vector
      \param cg_arg The solver parameters
      \param true_res Whether or not to report the true residual. This will
      point to the true residual  if it initially points to a non-zero value.
      \param cnv_frm Whether the lattice fields need to be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_out and \a f_in
  are also converted: This assumes they are initially in the same order as
  the gauge field. Fields that are converted are restored to their original
  order upon exit of this method. \e N.B. If the fields are already in the
  suitable order, then specifying ::CNV_FRM_YES here has no effect.
      \return The number of solver iterations.
      \post \a f_out contains the solution vector.
      \post true_res The true residual, if this was non-zero to start with.
      The residual is \f$ |f_{in} - M^\dagger M f_{out}| / |f_{in}| \f$.
    */

    int FmatEvlInv(Vector *f_out, Vector *f_in,  
			   CgArg *cg_arg,  
			   CnvFrmType cnv_frm = CNV_FRM_YES)
	{ return FmatEvlInv(f_out, f_in, cg_arg, 0, cnv_frm); }

    virtual int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
			    int Nshift, int isz, CgArg **cg_arg, 
			    CnvFrmType cnv_frm, MultiShiftSolveType type, 
			    Float *alpha, Vector **f_out_d) = 0;
    //!< The multishift matrix inversion used in the RHMC molecular dynamics algorithms.
    /*!<
      Solves \f$ (M^\dagger M + shift) f_{out} = f_{in} \f$ for \f$ f_{out}\f$
      for a given number of shifts,
      where \a M is the (possibly odd-even preconditioned) fermionic matrix.

      \param f_out The solution vectors.
      \param f_in The source vector
      \param shift The shifts of the fermion matrix.
      \param Nshift The number of shifts
      \param isz The smallest shift (required by MInvCG)
      \param cg_arg The solver parameters  for each shift
      \param cnv_frm Whether the lattice fields need to be converted to
      to a new storage order appropriate for the type of fermion action.
      If this is ::CNV_FRM_NO, then just the gauge field is converted.
      If this is ::CNV_FRM_YES, then the fields \a f_out and \a f_in
      are also converted: This assumes they are initially in the same order as
      the gauge field. Fields that are converted are restored to their original
      order upon exit of this method. \e N.B. If the fields are already in the
      suitable order, then specifying ::CNV_FRM_YES here has no effect.
      \param type The type of multimass inverter.
      If type is ::MULTI, then regular multishift inversion is performed with
      each solution stored separately.
      If type is ::SINGLE, then each solution is multiplied by an
      amount in parameter \a alpha and summed to a single solution vector.
      \param alpha The contribution of each shifted solution to the total
      solution vector if \a type is SINGLE.
      \param f_out_d Not implemented or ignored with Wilsonesque fermions.
      With staggered fermions, if this is initially non-zero
      and ::MULTI is also specified,
      then <em> D f<sub>out</sub></em>, the solution vectors acted on with the
      Dirac D-slash operator, is written here. 
      \return The number of solver iterations.
      \post \a f_out contains the solution vector(s).
    */

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift,
                            int Nshift, int isz, CgArg *cg_arg,
                            CnvFrmType cnv_frm, MultiShiftSolveType type,
                            Float *alpha, Vector **f_out_d);

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, CnvFrmType cnv_frm,
		    MultiShiftSolveType type, Float *alpha)
	{
	    Vector **f_out_d = 0;
	    return FmatEvlMInv(f_out,f_in,shift,Nshift,isz,cg_arg,cnv_frm,
			       type,alpha, f_out_d);
	}
    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg *cg_arg, CnvFrmType cnv_frm,
		    MultiShiftSolveType type, Float *alpha)
	{
	    Vector **f_out_d = 0;
	    return FmatEvlMInv(f_out,f_in,shift,Nshift,isz,cg_arg,cnv_frm,
			       type,alpha, f_out_d);
	}
   
    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, CnvFrmType cnv_frm,
		    Vector **f_out_d)
	{
	    Float *alpha = 0;
	    return FmatEvlMInv(f_out,f_in,shift,Nshift,isz,cg_arg,cnv_frm,
			       MULTI,alpha, f_out_d);
	}
    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg *cg_arg, CnvFrmType cnv_frm,
		    Vector **f_out_d)
	{
	    Float *alpha = 0;
	    return FmatEvlMInv(f_out,f_in,shift,Nshift,isz,cg_arg,cnv_frm,
			       MULTI,alpha, f_out_d);
	}

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, CnvFrmType cnv_frm)
	{
	    Float *alpha=0;
	    Vector **f_out_d=0; 
	    return FmatEvlMInv(f_out,f_in,shift,Nshift,isz,cg_arg,cnv_frm,
			       MULTI,alpha,f_out_d);
	}

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, CnvFrmType cnv_frm,
		    MultiShiftSolveType type, Vector **f_out_d)
	{
	    Float *alpha=0;
	    return FmatEvlMInv(f_out,f_in,shift,Nshift,isz,cg_arg,cnv_frm,
			       MULTI,alpha,f_out_d);
	}

    virtual void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
			    Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm) = 0;

    //! Chronological initial guess for the solver.
    /*!
      This computes the starting guess for the solver for use in the
      HMD force calculations using a  minimal residual chronological method.
      This computes the guess solution as a linear combination of a given
      number of the previous solutions.

      \param sol The chronological guess for the solution vector.
      \param source The source vector in the equation to be solved.
      \param sol_old The previous solutions.
      \param vm The previous solutions.multiplied by \f$ M^\dagger M \f$,
      computed as a necessary by-product.
      \param degree The number of previous solutions.
      \param cg_arg The solver parameters for the fermion matrix.
      \param cnv_frm Whether the lattice fields need to be converted to
      to a new storage order appropriate for the type of fermion action.
      If this is ::CNV_FRM_NO, then just the gauge field is converted.
      If this is ::CNV_FRM_YES, then the vectors \a sol and \a source 
      are also converted: This assumes they are initially in the same order as
      the gauge field. Fields that are converted are restored to their original
      order upon exit of this method. \e N.B. If the fields are already in the
      suitable order, then specifying ::CNV_FRM_YES here has no effect.
      
      \return The residue of this guess.
     */

    virtual int eig_FmatInv(Vector **V, const int vec_len, Float *M, const int nev, const int m, float **U, Rcomplex *invH, const int def_len, const Float *restart, const int restart_len,
			Vector *f_out, Vector *f_in, 
			CgArg *cg_arg, 
                        Float *true_res,
			CnvFrmType cnv_frm = CNV_FRM_YES,
			PreserveType prs_f_in = PRESERVE_YES);

    virtual int FmatInv(Vector *f_out, Vector *f_in, 
			CgArg *cg_arg, 
                        Float *true_res,
			CnvFrmType cnv_frm = CNV_FRM_YES,
			PreserveType prs_f_in = PRESERVE_YES) = 0;
    //!< Fermion matrix inversion.
    /*!<
      Solves <em> A f_out = f_in </em> for \a f_out, where \a A is the
      fermion matrix. The vectors must be defined on the whole lattice,
      not just on sites of a single parity.

      \param f_out The initial guess of solution vector.
      \param f_in The source vector
      \param cg_arg The solver parameters  for each shift
      \param true_res Whether or not to report the true residual. The true
      residual will be  written here if this is non-zero.
      \param cnv_frm Whether the lattice fields need to be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_out and \a f_in
  are also converted: This assumes they are initially in the same order as
  the gauge field. Fields that are converted are restored to their original
  order upon exit of this method. \e N.B. If the fields are already in the
  suitable order, then specifying ::CNV_FRM_YES here has no effect.
      \param prs_f_in Whether or not the source vector is allowed to be
      overwritten, thereby saving memory. For staggered fermions \a f_in is
      preserved regardless of the value of \a prs_f_in. 
      \return The number of solver iterations.
      \post \a f_out contains the solution vector.
      \post \a true_res contains  the true residual, if it was non-zero
      to start with.
      The residual is  <em>  |f_in - A f_out| / |f_in| </em>
    */

    int FmatInv(Vector *f_out, Vector *f_in,
			CgArg *cg_arg, 
			CnvFrmType cnv_frm = CNV_FRM_YES,
			PreserveType prs_f_in = PRESERVE_YES)
	{ return FmatEvlInv(f_out, f_in, cg_arg, 0, cnv_frm); }

	virtual void Fsolfour2five(Vector *sol_5d, Vector *sol_4d, Vector *src_5d, CgArg *cg_arg){
		char *fname = "Fsolfour2five()";
		ERR.NotImplemented(cname,fname);
	}
    // Recover the 5D solution from the 4D solution, without solve the equation again.
    
  //!< It solves the eigenvectors and eigenvalues of the fermion matrix using implicitly restarting Lanczos with polynomial accerelation with shifts.
    /*!<
      \param f_eigenv The computed eigenvalues
      \param lambda The corresponding eigenvalues
      \param eig_arg 
      \param cnv_frm Whether the lattice fields need to be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_out and \a f_in
  are also converted: This assumes they are initially in the same order as
  the gauge field. Fields that are converted are restored to their original
  order upon exit of this method. \e N.B. If the fields are already in the
  suitable order, then specifying ::CNV_FRM_YES here has no effect.
      \return The number of eigensolver iterations.
      \post f_eigenv contains the eigenvectors.
      \post lambda contains the eigenvalues.

      Currently only impolemented for even/odd preconditioned matrix (extention is easy though).
      
    */
    virtual int FeigSolv(Vector **f_eigenv, Float *lambda, 
			 LanczosArg *eig_arg, 
			 CnvFrmType cnv_frm = CNV_FRM_YES){
		char *fname = "FeigSolv(**V,F*,L*,C)";
		ERR.NotImplemented(cname,fname);
		return -1;
	}

  
    //!< It solves the eigenvectors and eigenvalues of the fermion matrix Using Ritz function minimization.
    /*!<
      \param f_eigenv The computed eigenvalues
      \param lambda The corresponding eigenvalues
      \param chirality eigenvector(i)^dagger gamma_5 eigenvector(i)
      \param valid_eig
      \param hsum
      \param eig_arg 
      \param cnv_frm Whether the lattice fields need to be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_out and \a f_in
  are also converted: This assumes they are initially in the same order as
  the gauge field. Fields that are converted are restored to their original
  order upon exit of this method. \e N.B. If the fields are already in the
  suitable order, then specifying ::CNV_FRM_YES here has no effect.
      \return The number of eigensolver iterations.
      \post f_eigenv contains the eigenvectors.
      \post lambda contains the eigenvalues.
     */
    virtual int FeigSolv(Vector **f_eigenv, Float *lambda, 
			 Float *chirality, int *valid_eig,
			 Float **hsum,
			 EigArg *eig_arg, 
			 CnvFrmType cnv_frm = CNV_FRM_YES) = 0;

  
    virtual Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
			 Float mass, DagType dag) = 0;

    //!< Initialises the pseudofermion field.
    /*!<
      The heatbath initialisation of the pseudofermion field
      is done by setting \f$ \phi = M^\dagger \eta \f$ where \f$ \eta \f$
      is a zero mean, unit variance random gaussian field.
      The pseudofermion field may be computed on a single parity only.	     
    
      \pre The random field must already be initialised.
      \param phi The pseudofermion field.
      \param frm1 A random field, possibly on a single parity.
      \param frm2 Another random field, or possibly workspace.
      \param mass The mass parameter of the fermion matrix.       
      \param dag Whether DAG_YES or DAG_NO is applied
      \return The value of the pseudofermionic action on this node (only 
      works for Wilson type fermions).
    */

    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		 Float mass) {
      return SetPhi(phi,frm1,frm2,mass,DAG_YES);
    }

    virtual ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size) = 0;
    //!< Molecular dynamics evolution due to the fermion force
    /*!<
      The momentum is evolved for a single molecular dynamics timestep
      using the force from the fermion action.
      \param mom The momentum matrices on all links of the local lattice.
      \param frm The solution of the fermion matrix inverse with the
      pseudofermion vector source, as computed by FmatEvlInv.
      \param mass The FERMION mass (not used in staggered fermion classes).
      \param step_size The molecular dynamics timestep used in the numerical
      integration.
      \post \a mom is assigned the value of the momentum after the molecular
      dynamics evolution.
    */

    virtual ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
				  Float mass, Float step_size) = 0;
    //!< Molecular dynamics evolution due to the boson part of quotient integrator
    /*!<
      The momentum is evolved for a single molecular dynamics timestep
      using the force from the bosonic part of a quotient action.
      \param mom The momentum matrices on all links of the local lattice.
      \param frm The solution of the fermion matrix inverse with the
      pseudofermion vector source, as computed by FmatEvlInv.
      \param eta The solution of the fermion matrix inverse applied to 
      bosonic operator times the pseudofermion vector source.
      \param mass The BOSON/FERMION mass (not used in staggered fermion classes).
      \param step_size The molecular dynamics timestep used in the numerical
      integration.
      \post \a mom is assigned the value of the momentum after the molecular
      dynamics evolution.
    */

    virtual ForceArg RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
					  int isz, Float *alpha, Float mass, 
					  Float dt, Vector **sol_d, ForceMeasure measure) = 0;

    virtual Float FhamiltonNode(Vector *phi, Vector *chi) = 0;
    //!< Computes the pseudofermionic action on the local sublattice.
    /*!<
      \pre The equation <em> A chi = phi </em>  needs to have been solved for
      chi, where \e A is the inverse of the fermionic matrix in the molecular
      dynamics hamiltonian.
      of the matrix in the molecular dynamics hamiltonian. 
      \param phi The pseudofermion field
      \param chi The solution of <em> A chi = phi </em> 
      where \e A is the inverse of the fermionic matrix in the molecular
      dynamics hamiltonian.
      \return The value of the  pseudofermionic action on this node.
    */

    virtual void Fconvert(Vector *f_field, 
		   	  StrOrdType to,
			  StrOrdType from, int cb=2){
		char *fname = "Fconvert(*V,O,O)";
		ERR.NotImplemented(cname,fname);
	}
    //!< Converts the field layout.
    /*!
      Exactly which data layouts are supported depends on the type of
      fermion action.
      The field exists on all lattice sites (both parities).
      \param f_field The field to be converted.
      \param to The new order.
      \param from The current order.
    */
    void Fconvert(LatVector *f_field, 
		   	  StrOrdType to,
		  StrOrdType from, int cb=2){
      Fconvert(f_field->Vec(),to,from);
    }


//------------------------------------------------------------------
//------------------------------------------------------------------
// Bosonic action related pure virtual functions.
//------------------------------------------------------------------
//------------------------------------------------------------------

    virtual Float BhamiltonNode(Vector *boson, Float mass) = 0;
        // The boson Hamiltonian of the node sublattice.
  
    void ClearSmeared() {smeared = 0;}

    void Shift();

    //!< Method to ensure bosonic force works (only applies to
    //!< staggered fermion formulations).
    virtual void BforceVector(Vector *in, CgArg *cg_arg) = 0;

    // added by Hantao to facilitate doing force statistics.
    void updateForce(ForceArg &f_arg, const Matrix &m)const {
        Float a2 = m.norm();
        Float a = sqrt(a2);

        f_arg.L1 += a;
        f_arg.L2 += a2;
        f_arg.Linf = f_arg.Linf > a ? f_arg.Linf : a;
    }

    //!< Toggle boundary condition
    //
    //!< Note: Agent classes which needs to import gauge field to
    //!external libraries need to overwrite this function.
    virtual void BondCond();
};

//------------------------------------------------------------------
//------------------------------------------------------------------
// DERIVED CLASSES
//------------------------------------------------------------------
//------------------------------------------------------------------

/*! \defgroup gactions Gauge Actions
  \ingroup latactions */

//------------------------------------------------------------------
//! A class implementing a lattice with a zero gauge action.
/*!
  Whatever that means. Most of the methods do nothing.
  \ingroup gactions
*/
//------------------------------------------------------------------
class Gnone : public virtual Lattice
{

 private:
    const char *cname;    // Class name.

 public:

    Gnone();

    virtual ~Gnone();

    GclassType Gclass();
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
    //!< Calculates the gauge force at site x and direction mu.

    ForceArg EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode();
       // The pure gauge Hamiltonian of the node sublattice.

    void AllStaple(Matrix &stap, const int *x, int mu);
    //!< Not implemented.

};


//------------------------------------------------------------------
//! A class implementing a lattice with the standard Wilson plaquette gauge action.
/*! \ingroup gactions */
//------------------------------------------------------------------
class Gwilson : public virtual Lattice
{

 private:
    const char *cname;    // Class name.

 public:

    Gwilson();

    virtual ~Gwilson();

    GclassType Gclass();
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
    //!< Calculates the gauge force at site x and direction mu.

    ForceArg EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode();
       // The pure gauge Hamiltonian of the node sublattice.

    void AllStaple(Matrix &stap, const int *x, int mu);
    //!< Calculates the gauge field square staple sum  around a link
    
};


//------------------------------------------------------------------
//! A class implementing a lattice with the power-plaquette gauge action.
/*!
  This action is the standard Wilson action with an irrelevant power plaquette
 term added to it. The action is:
\f[
 \sum_p  -\frac{1}{3}\beta \mathrm{ReTr}[U_p] + \left(\frac{1}{c}(1-\frac{1}{3}\mathrm{ReTr}[U_p]) \right)^k 
\f]
 where \f$ U_p \f$ is the plaquette and the sum is over all plaquettes.
 \f$ \beta \f$, \a c and \a k are real parameters.
 (see GlobalJobParameter::PowerPlaqCutoff,
      GlobalJobParameter::PowerPlaqExponent and
      GlobalJobParameter::Beta)

   \ingroup gactions
*/
//------------------------------------------------------------------
class GpowerPlaq : public virtual Lattice
{

 private:
    const char *cname;    // Class name.

 public:

    GpowerPlaq();

    virtual ~GpowerPlaq();

    GclassType Gclass();
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
    //!< Calculates the gauge force at site x and direction mu.

    ForceArg EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode();
       // The pure gauge Hamiltonian of the node sublattice.

    void PowerStaple(Matrix& pstap, int *x, int mu);
    //!< Computes the power-plaquette staple sum around a link.

    Float PowerPlaq(int *x, int mu, int nu) const;
        //!< Calculates the power plaquette term.

    Float SumPowerPlaqNode() const;
        //!< Calculates the local sum of the power plaquette term.

    Float SumPowerPlaq() const;
        //!< Calculates the global sum of the power plaquette term.    

    void AllStaple(Matrix &stap, const int *x, int mu);
};

//------------------------------------------------------------------
//! A class implementing a lattice with the plaquette + rectangle gauge action.
/*!
  The action is
  \f[
  -\frac{1}{3}\beta \sum_x \sum_\mu \sum_{\nu\neq\mu}[
  (1-8 c_1) U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x) 
  +c_1 U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu)
U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
]
\f]

  \ingroup gactions
*/
//------------------------------------------------------------------
class GimprRect : public virtual Lattice
{

 private:
    const char *cname;    // Class name.

    Float plaq_coeff; // - GJP.Beta() * ( 1.0 - 8.0 * GJP.C1() ) / 3.0

    Float rect_coeff; // - GJP.Beta() * (             GJP.C1() ) / 3.0

    static unsigned CBUF_MODE4;
	
 public:

    GimprRect();

    virtual ~GimprRect();

    GclassType Gclass();
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
    //!< Calculates the gauge force at site x and direction mu.
        // Typical implementation has this func called with
        // Matrix &force = *mp0.  GactionGradient typically uses
        // mp1 thru mp4, so be careful.

    ForceArg EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode();
       // The pure gauge Hamiltonian of the node sublattice.

    //! Computes sum of the plaquette and rectangle staples around a link.
    void AllStaple(Matrix &stap, const int *x, int mu);

};

//------------------------------------------------------------------
//! A class implementing a tadpole-improved plaquette + rectangle gauge action.
/*!
  The action is
  \f[
  -\frac{1}{3}\beta \sum_x \sum_\mu \sum_{\nu\neq\mu}[
  (1-8 c_1) U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x) 
  + \frac{c_1}{u_0^2} U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu)
U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
]
\f]

  \ingroup gactions
*/
//------------------------------------------------------------------
class GtadpoleRect : public virtual Lattice
{

 private:
    const char *cname;    // Class name.

    Float plaq_coeff; // - GJP.Beta() * ( 1.0 - 8.0 * GJP.C1() ) / 3.0

    Float rect_coeff; // - GJP.Beta() * ( GJP.C1()/(GJP.u0()*GJP.u0()) ) / 3.0

    static unsigned CBUF_MODE4;

 public:

    GtadpoleRect();

    virtual ~GtadpoleRect();

    GclassType Gclass();
        // It returns the type of gauge class

   void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
    //!< Calculates the gauge force at site x and direction mu.
        // Typical implementation has this func called with
        // Matrix &force = *mp0.  GactionGradient typically uses
        // mp1 thru mp4, so be careful.

    ForceArg EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode();
       // The pure gauge Hamiltonian of the node sublattice.

    //! Computes sum of the plaquette and rectangle staples around a link.
    void AllStaple(Matrix &stap, const int *x, int mu);

};


//------------------------------------------------------------------
//! A class implementing a lattice with the power-rectangle gauge action.
/*!
The action is:
 \f[
-\beta \sum_x \sum_p
[ (1-8c_1)\frac{1}{3}\mathrm{ReTr}[U_p] + \{ \frac{1}{c}(1-\frac{1}{3}\mathrm{ReTr}[U_p]) \}^k ]
+ \sum_r [c_1\frac{1}{3}\mathrm{ReTr}[U_r] + \{ \frac{1}{c}(1-\frac{1}{3}\mathrm{ReTr}[U_r]) \}^k ]
\f]
    where the sum is over all plaquettes
\f[
  U_p(x, \mu, \nu) = U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f]
    and all rectangles
\f[
  W_r(x, \mu, \nu) = U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
\f] 

 This action supresses plaquettes with \f$ 1 - \mathrm{ReTr}[U_p]/3 > c \f$
 and rectangles with \f$ 1 - \mathrm{ReTr}[W_r]/3 > c \f$
 and therefore reduces lattice dislocations.

   \ingroup gactions
 */
//------------------------------------------------------------------
class GpowerRect : public virtual Lattice
{

 private:
    const char *cname;    // Class name.

    Float plaq_coeff; // - GJP.Beta() * ( 1.0 - 8.0 * GJP.C1() ) / 3.0

    Float rect_coeff; // - GJP.Beta() * (             GJP.C1() ) / 3.0

 public:

 GpowerRect();

    virtual ~GpowerRect();

    GclassType Gclass();
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
    //!< Calculates the gauge force at site x and direction mu.
        // Typical implementation has this func called with
        // Matrix &force = *mp0.  GactionGradient typically uses
        // mp1 thru mp4, so be careful.

    ForceArg EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode();
       // The pure gauge Hamiltonian of the node sublattice.

    void PowerStaple(Matrix& pstap, int *x, int mu);
        //!< Calculates the sum of the plaquette staples around a link.

    Float PowerPlaq(int *x, int mu, int nu) const;
        //!< Calculates the real part of the trace of the power plaquette.

    Float SumPowerPlaqNode() const;
        //!< Calculates the local sum of the real part of the trace of the power plaquette.

    Float SumPowerPlaq() const;
        //!< Calculates the global sum of the real part of the trace of the power plaquette.

    void PowerRectStaple(Matrix& pstap, int *x, int mu);
        // ! Calculates the rectangle staple sum around a link.

    Float PowerRect(int *x, int mu, int nu) const;
        //!< Calculates the the real part of the trace of the power rectangle.

    Float SumPowerRectNode() const;
        //!< Calculates the local sum of the real part of the trace of the power rectangle.

    Float SumPowerRect() const;
        //!< Calculates the global sum of the real part of the trace of the power rectangle.

    void AllStaple(Matrix &stap, const int *x, int mu);
    //!< Not implemented.
};

//------------------------------------------------------------------
//! A class implementing a lattice with the 1-loop Symanzik improved gauge action.
/*!
  This is a sum of the plaquette, rectangle and cube loops:
  \f[
  \sum_x\sum_\mu\sum_{\nu>\mu} \left[\right.
U_\mu(x) U_\nu(x+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
  \f]\f[
  -\frac{ 1 + 0.4805 \alpha_s}{20 u_0^2} 
  U_\mu(x) U_\mu(x+\mu) U_\nu(x+2\mu) U^\dagger_\mu(x+\mu+\nu) U^\dagger_\mu(x+\nu) U^\dagger_\nu(x)
  \f]\f[
    -\frac{ 0.03325 \alpha_s}{u_0^2} 
     \sum_{\rho>\nu} U_\mu(x) U_\nu(x+\mu) U_\rho(x+\mu+\nu) U^\dagger_\mu(x+\mu+\nu+\rho)
     U^\dagger_\nu(x+\nu+\rho) U^\dagger_\rho(x+\rho)
     \left.\right]
       \f]
       where \f$u_0\f$ is the tadpole coefficient
       and \f$\alpha_s = -4\log(u_0)/3.06839\f$.       

  \ingroup gactions
*/
//------------------------------------------------------------------
class GimprOLSym : public virtual Lattice
{

 private:
    const char *cname;    // Class name.

    Float rect_coeff; // loop 2x1
    Float cube_coeff; // loop 1x1x1
    Float minus_beta_over_3 ;// coefficient needed to calculate the force 

 public:

    GimprOLSym();

    virtual ~GimprOLSym();

    GclassType Gclass();
        // It returns the type of gauge class

    void GactionGradient(Matrix &grad, int *x, int mu) ;
        // Calculates the partial derivative of the gauge action
        // w.r.t. the link U_mu(x).  Typical implementation has this
        // func called with Matrix &grad = *mp0, so avoid using it.

    void GforceSite(Matrix& force, int *x, int mu);
        //!< Calculates the gauge force at site x and direction mu.
        // Typical implementation has this func called with
        // Matrix &force = *mp0.  GactionGradient typically uses
        // mp1 thru mp4, so be careful.

    ForceArg EvolveMomGforce(Matrix *mom, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the pure gauge force.

    Float GhamiltonNode();
       // The pure gauge Hamiltonian of the node sublattice.

    void AllStaple(Matrix &stap, const int *x, int mu);
    //!< Computes the sum of all the staples around a link.

};

//-------------------------------------------------------------------
//! A class for 'improved' actions using smeared links.
/*!
  This class encapsulates the production of smeared gauge links which
  are part of the definition of certain fermion actions.
  The parameter is the number of different smeared fields required. This
  is action-dependent, as is the nature of the smearings and the order in
  which they are stored,  
 */
//------------------------------------------------------------------

class Fsmear : public virtual Lattice
{
private:
    const char *cname;    // Class name.

protected:
    int n_fields;    //!< The number of smeared fields
    Matrix **fields; //!< The smeared links.

public:

    Fsmear(int n_smear);

    ~Fsmear();

    //! Get (a pointer to) the <em>n</em>th smeared field
    Matrix *Fields(int n);

    //! Make the smeared fields.
    virtual void Smear() = 0 ; //pure virtual function

};

//----------------------------------------------------------------------
/*! \defgroup factions Fermion actions
  \ingroup latactions */
//! A class implementing a lattice with a zero fermion action.
/*!
  Most of the methods do nothing.
  \ingroup factions
*/
//------------------------------------------------------------------
// Fnone is derived from Lattice. Its functions do nothing
// and return values as if there is no fermion action or
// fermion fields. The number of spin components is zero
// The site size of the fermion array FsiteSize() is 
// set to 1 so that memory allocation would proceed normally.
//------------------------------------------------------------------
class Fnone : public virtual Lattice
{
 private:
    const char *cname;    // Class name.
    
 public:

    Fnone();

    virtual ~Fnone();

    FclassType Fclass() const;
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

    int ExactFlavors() const;
        // Returns the number of exact flavors of the matrix that
        // is inverted during a molecular dynamics evolution.

    int SpinComponents() const;
        // Returns the number of spin components.

    int FsiteSize() const;

    int FchkbEvl() const;
	// returns 1 => The fermion fields in the evolution
        //      or the CG that inverts the evolution matrix
	//      are defined on a single checkerboard (half the 
	//      lattice).

   int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);
        // It does nothing and returns 0.

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, 
		    CnvFrmType cnv_frm, MultiShiftSolveType type, 
		    Float *alpha, Vector **f_out_d);
        // It does nothing and returns 0.    

    void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
		     Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm);
        // It does nothing and returns 0.    

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);
        // It does nothing and returns 0.

    int FeigSolv(Vector **f_eigenv, Float *lambda, 
		 LanczosArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES){ return 0; };
    
    int FeigSolv(Vector **f_eigenv, Float *lambda,
		 Float *chirality, int *valid_eig,
		 Float **hsum,
		 EigArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
        // It does nothing and returns 0.

    virtual Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
			 Float mass, DagType dag);
	// It sets the pseudofermion field phi from frm1, frm2.

    void FforceSite(Matrix& force, Vector *frm, 
                            int *x, int mu);
    //!< Calculates the pseudofermion force at site x and direction mu.
        // It calculates the fermion force per site x
        // and direction mu. frm is the fermion field that 
        // resulted from the application of the inverter on 
        // the pseudofermion field.

    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
				 Float mass, Float step_size);
        // It evolves the canonical momentum mom by step_size
        // using the fermion force. 

    ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
			  Float mass, Float step_size);
        // It evolve the canonical momentum mom  by step_size
        // using the bosonic quotient force.

    ForceArg RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
			      int isz, Float *alpha, Float mass, Float dt,
			      Vector **sol_d, ForceMeasure measure);

    Float FhamiltonNode(Vector *phi, Vector *chi);
        // The fermion Hamiltonian of the node sublattice.
        // chi must be the solution of Cg with source phi.	       

    void Fconvert(Vector *f_field,
			  StrOrdType to,
		  StrOrdType from, int cb=2);
        // Convert fermion field f_field from -> to

    Float BhamiltonNode(Vector *boson, Float mass);
        // The boson Hamiltonian of the node sublattice.

    //!< Method to ensure bosonic force works.
    void BforceVector(Vector *in, CgArg *cg_arg);

};


//------------------------------------------------------------------
//! A class containing methods relevant to all staggered fermion actions.
//------------------------------------------------------------------
class FstagTypes : public virtual Lattice
{
	 private:
    const char *cname;    // Class name.
    int xv[3];
    
  protected:

    enum
      {
	VECT_LEN=6,          //!< Number of Floats in  a Vector
	MATRIX_SIZE=18,  //!< Number of Floats in  a Matrix
	SITE_LEN=72      //!< Number of Floats in four Matrix's
      };
    
    int bc[4];	        //!< Boundary conditions
    int e_vsize;	//!< Size of a single parity vector field

    static const unsigned CBUF_MODE1 = 0xcb911548;
    static const unsigned CBUF_MODE2 = 0xcca52112;
    static const unsigned CBUF_MODE3 = 0xc98c6106;
    static const unsigned CBUF_MODE4 = 0xcca52112;

    Vector *f_tmp; 
    
 public:

    FstagTypes();

    virtual ~FstagTypes();

    virtual FclassType Fclass() const = 0;

    int FsiteOffsetChkb(const int*) const;

    //!< Gets the lattice site index for the odd-even (checkerboard) order.
    int FsiteOffsetChkb_all(const int*) const;
	
    int FsiteOffset(const int*) const;

    int ExactFlavors() const;

    int SpinComponents() const;

    int FsiteSize() const;

    int FchkbEvl() const;

    void Fconvert(Vector*, StrOrdType, StrOrdType, int cb=2);
    
    Float FhamiltonNode( Vector*,  Vector*) ;
    
    void massRenormalise(Float *mass, Float *trueMass, int degree, 
			 Float *shift, MassRenormaliseDir direction);
    //!< Renormalise the mass (staggered optimisation)

};

//------------------------------------------------------------------
//! A class implementing staggered fermions.
/*!
  \ingroup factions
*/
//------------------------------------------------------------------
class Fstag : public virtual FstagTypes
{
 private:
    const char *cname;    // Class name.
    void getUDagX(Vector& v, const Vector *cvp, int *x, int mu) const;
    
 public:

    Fstag();

    virtual ~Fstag();

    FclassType Fclass() const; 

    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, 
		    CnvFrmType cnv_frm, MultiShiftSolveType type, Float *alpha,
		    Vector **f_out_d);

    void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
		     Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm);

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);

    int FeigSolv(Vector **f_eigenv, Float *lambda, 
		 LanczosArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);    

    int FeigSolv(Vector **f_eigenv, Float *lambda, 
		 Float *chirality, int *valid_eig,
		 Float **hsum,
		 EigArg *eig_arg,
		 CnvFrmType cnv_frm = CNV_FRM_YES);

    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		 Float mass, DagType dag);

    void FforceSite(Matrix& force, Vector *frm, 
                            int *x, int mu);
    //!< Calculates the pseudofermion force at site x and direction mu.

    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
			 Float mass, Float step_size);

    ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
			  Float mass, Float step_size);
        // It evolve the canonical momentum mom  by step_size
        // using the bosonic quotient force.

    ForceArg RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
			       int isz, Float *alpha, Float mass, Float dt,
			       Vector **sol_d, ForceMeasure measure);

    Float BhamiltonNode(Vector *boson, Float mass);

    //!< Method to ensure bosonic force works.
    void BforceVector(Vector *in, CgArg *cg_arg);

    void Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int dir_flag);

    void FdMdmu(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int order);

};


class ParTransAsqtad; //forward declaration

//------------------------------------------------------------------
//! A class implementing improved staggered fermions (the asqtad action).
/*!
  \ingroup factions
*/
//------------------------------------------------------------------

class Fasqtad : public virtual FstagTypes, public virtual Fsmear
{
 private:
    const char *cname;    // Class name.

 public:

    Fasqtad();
    virtual ~Fasqtad();

    FclassType Fclass() const;
    
    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, 
		    CnvFrmType cnv_frm, MultiShiftSolveType type, Float *alpha,
		    Vector **f_out_d);

    void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
		     Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm);

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);

    int FeigSolv(Vector **f_eigenv, Float *lambda, 
		 LanczosArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
    
    int FeigSolv(Vector **f_eigenv, Float *lambda, 
		 Float *chirality, int *valid_eig,
		 Float **hsum,
		 EigArg *eig_arg,
		 CnvFrmType cnv_frm = CNV_FRM_YES);

    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		 Float mass, DagType dag);

    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
			 Float mass, Float step_size);

    ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
			  Float mass, Float step_size);
        // It evolve the canonical momentum mom  by step_size
        // using the bosonic quotient force.

    Float BhamiltonNode(Vector *boson, Float mass);

    //!< Method to ensure bosonic force works.
    void BforceVector(Vector *in, CgArg *cg_arg);

    void Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int dir_flag);

#if 0
    void FdMdmu(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int order);
#endif

    //! Momentum update in the RHMC algorithm.
    ForceArg RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
			       int isz, Float *alpha, Float mass, Float dt,
			       Vector **sol_d, ForceMeasure measure);

    // Various utility routines for the momentum force computation.

    void Smear();
    
  private:

    ForceArg update_momenta(Matrix**, IFloat, Matrix*);
    
    ChkbType parity(const int*);

    void force_product_sum(const Matrix*,  int, IFloat, Matrix*);
    
    void force_product_sum(const Matrix*, const Matrix*, const Matrix*, 
			   IFloat, Matrix*, Matrix*);

    void force_product_sum(const Matrix*, const Matrix*, IFloat, 
			   Matrix*);

    void force_product_d_sum(const Matrix*, const Matrix*, IFloat, 
			     Matrix*);

    void force_product_sum(const Vector*, const Vector*, IFloat, Matrix*);


    
};


class ParTransP4; //forward declaration

//------------------------------------------------------------------
//! A class implementing improved staggered fermions (the P4 action).
/*!
  \ingroup factions
*/
//------------------------------------------------------------------

class Fp4 : public virtual FstagTypes, public virtual Fsmear
{
 private:
    const char *cname;    // Class name.

 public:

    Fp4();
    virtual ~Fp4();

    FclassType Fclass() const;
    
    int FmatEvlInv(Vector *f_out, Vector *f_in, 
		   CgArg *cg_arg, 
		   Float *true_res,
		   CnvFrmType cnv_frm = CNV_FRM_YES);

    int FmatEvlMInv(Vector **f_out, Vector *f_in, Float *shift, 
		    int Nshift, int isz, CgArg **cg_arg, 
		    CnvFrmType cnv_frm, MultiShiftSolveType type, Float *alpha,
		    Vector **f_out_d);

    void FminResExt(Vector *sol, Vector *source, Vector **sol_old, 
		     Vector **vm, int degree, CgArg *cg_arg, CnvFrmType cnv_frm);

    int FmatInv(Vector *f_out, Vector *f_in, 
		CgArg *cg_arg, 
		Float *true_res,
		CnvFrmType cnv_frm = CNV_FRM_YES,
		PreserveType prs_f_in = PRESERVE_YES);

    int FeigSolv(Vector **f_eigenv, Float *lambda, 
		 LanczosArg *eig_arg, 
		 CnvFrmType cnv_frm = CNV_FRM_YES);
    
    int FeigSolv(Vector **f_eigenv, Float *lambda, 
		 Float *chirality, int *valid_eig,
		 Float **hsum,
		 EigArg *eig_arg,
		 CnvFrmType cnv_frm = CNV_FRM_YES);

    Float SetPhi(Vector *phi, Vector *frm1, Vector *frm2,
		 Float mass, DagType dag);

    ForceArg EvolveMomFforce(Matrix *mom, Vector *frm, 
			 Float mass, Float step_size);

    ForceArg EvolveMomFforce(Matrix *mom, Vector *phi, Vector *eta,
			  Float mass, Float step_size);
        // It evolve the canonical momentum mom  by step_size
        // using the bosonic quotient force.

    Float BhamiltonNode(Vector *boson, Float mass);

    //!< Method to ensure bosonic force works.
    void BforceVector(Vector *in, CgArg *cg_arg);

    void Fdslash(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int dir_flag);

    void FdMdmu(Vector *f_out, Vector *f_in, CgArg *cg_arg, 
		 CnvFrmType cnv_frm, int order);


    //! Momentum update in the RHMC algorithm.
    ForceArg RHMC_EvolveMomFforce(Matrix *mom, Vector **sol, int degree,
			      int isz, Float *alpha, Float mass, Float dt,
			      Vector **sol_d, ForceMeasure measure);

    // Various utility routines for the momentum force computation.

    void Smear();
    
  private:

    ForceArg update_momenta(Matrix**, IFloat, Matrix*);
    
    ChkbType parity(const int*);

    void force_product_sum(const Matrix*,  int, IFloat, Matrix*);
    
    void force_product_sum(const Matrix*, const Matrix*, const Matrix*, 
			   IFloat, Matrix*, Matrix*);

    void force_product_sum(const Matrix*, const Matrix*, IFloat, 
			   Matrix*);

    void force_product_d_sum(const Matrix*, const Matrix*, IFloat, 
			     Matrix*);

    void force_product_sum(const Vector*, const Vector*, IFloat, Matrix*);
};

// specific to dwf , in src/util/dirac_op/d_op_mobius/d_op_mobius.C
void ReflectAndMultGamma5( Vector *out, const Vector *in,  int nodevol, int ls);
void HermicianDWF_ee( Vector* vtmp, Vector* evec, Float mass, Lattice* lattice, Vector* Apsi );

CPS_END_NAMESPACE

#include <util/lattice/f_wilson_types.h>
#include <util/lattice/lattice_types.h>



#endif
