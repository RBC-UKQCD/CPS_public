#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_clover/noarch/d_op_clover_supp.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include<config.h>
#include<util/dirac_op.h>
CPS_START_NAMESPACE

//////////////////////////////////////////////////////////////////////////
// d_op_clover_supp.C:
//        implement all the private auxiliary member functions for
//        class DiracOpClover.
//////////////////////////////////////////////////////////////////////////


CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/vector.h>
#include<util/enum.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/wilson.h>
#include<util/clover.h>
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include<comms/scu.h>
#include<util/dense_matrix.h>
CPS_START_NAMESPACE

//--------------------------------------------------------------------------
// enum
//--------------------------------------------------------------------------
enum {
  LINK_SIZE = 18,
  CLOVER_MAT_SIZE = 72,
  HALF_CLOVER_MAT_SIZE = 36,
  HALF_CLOVER_MAT_RANK = 6,
  HALF_LINK_MAT_SIZE = 9  
  //,  SPINOR_SIZE = 24,       defined in wilson.h
  //,  ND = 4                  defined in wilson.h
};

//--------------------------------------------------------------------------
// CRAM scratch buffer: 
//--------------------------------------------------------------------------
//   mp0, mp1, and mp2 for SiteFuv, 
//   additional mp3 for SiteCloverMat
//--------------------------------------------------------------------------
#ifdef _TARTAN
extern const unsigned int clover_cram_scratch_addr;
static Matrix *mp0 = (Matrix *)clover_cram_scratch_addr;      
#else
static IFloat scratch[96];
static Matrix *mp0 = (Matrix *)scratch;
#endif
static Matrix *mp1 = mp0 + 1;
static Matrix *mp2 = mp1 + 1;
static Matrix *mp3 = mp2 + 1;

//--------------------------------------------------------------------------
// DRAM scratch buffer: 
//--------------------------------------------------------------------------
// Used by GetLink to temporary store an off-node link
//--------------------------------------------------------------------------
static Matrix mtmp;

//--------------------------------------------------------------------------
// const Matrix &
// DiracOpClover::GetLink(const int *site, int dir) const
//--------------------------------------------------------------------------
// Purpose: 
//   get a link at specified coordinates and direction, whether on node
//   or off node. The gauge fields are assumed to be in Wilson order.
// Arguments:
//   site:    lattice coordinates[x,y,z,t]
//            which could be out-of-range, i.e., is located off-node.
//   dir:     0,1,2,3 for Ux, Uy, Uz, Ut
//   return:  a pointer to the link. If off-node, it points to
//            a block of static memory. Be careful to use it!
// Wilson order:
//     1. Canonical order is U[t][z][y][x][x,y,z,t]
//     2. Wilson links:   checkerboarded and  EVEN sites first!
//        Wilson spinors: checkerboarded and  ODD  sites first!
//--------------------------------------------------------------------------
const Matrix &
DiracOpClover::GetLink(const int *site, int dir) const
{ 
  //  char *fname = "GetLink()";
  //  VRB.Func(cname,fname);

  // number of local lattice sites
  //------------------------------------------------------------------------
  const int *nsites = clover_lib_arg->nsites;  
  
  // offset out-of-range coordinates site[] into on_node_site[]
  // in order to locate the link
  //------------------------------------------------------------------------
  int on_node_site[4];  
  int on_node = 1;  
  const Matrix *on_node_link;  
  {  
    const Matrix* Ue = lat.GaugeField();         // checkerboarded     
    const Matrix* Uo = Ue + GJP.VolNodeSites()/2 * ND; 
    const int site_offset[4] = {1, nsites[0], nsites[0]*nsites[1],
				nsites[0]*nsites[1]*nsites[2]  };
    int parity = 0;
    int offset = 0;  
    for (int i = 0; i < 4; ++i) {
      on_node_site[i] = (site[i] + nsites[i]) % nsites[i];    
      if (on_node_site[i] != site[i])
	on_node = 0;      
      parity += on_node_site[i];
      offset += site_offset[i] * on_node_site[i];    
    }
    offset = (offset - (offset%2)) / 2;          // checkerboarded
    on_node_link = ((parity%2) == 1 ? Uo : Ue) + offset * ND + dir;  
  }

#ifndef PARALLEL
  return *on_node_link;
#endif

  // send to the destination node if the site is off-node
  //------------------------------------------------------------------------
  if (on_node) {
    return *on_node_link;
  } else {
    Matrix send = *on_node_link;
    Matrix &recv = mtmp;
    for (int i = 0; i < 4; ++i) {
      while (site[i] != on_node_site[i]) {
	if (site[i] < 0) {
	  getMinusData((IFloat *)&recv, (IFloat *)&send, sizeof(recv)/sizeof(IFloat), i);
	  on_node_site[i] -= nsites[i];
	} else {
	  getPlusData((IFloat *)&recv, (IFloat *)&send, sizeof(recv)/sizeof(IFloat), i);
	  on_node_site[i] += nsites[i];
	}
	send = recv;      
      }
    }
    return mtmp;    
  }
}

//--------------------------------------------------------------------------
// void 
// DiracOpClover::SiteFuv(Matrix &Fuv, const int *site, int mu, int nu) const
//--------------------------------------------------------------------------
// Purpose: 
//    calculate F(x, u, v) = 1/4 1/2 (PLAQs - PLAQs^dag)
// where PLAQs are:
//    U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag +
//    U_nu(x) U_mu(x-mu+nu)^dag U_nu(x-mu)^dag U_mu(x-mu) +
//    U_mu(x-mu)^dag U_nu(x-mu-nu)^dag U_mu(x-mu-nu) U_nu(x-nu) +
//    U_nu(x-nu)^dag U_mu(x-nu) U_nu(x+mu-nu) U_mu(x)^dag
// Note:
//  1. We adopt same definition of clover term as in Xiang-Qian Luo's paper,
//     But slightly different notation:
//        Fuv (here) / i = Fuv (Luo)
//  2. MATRIX STORAGE CONVENTION:
//     2.1. All the comments in this functions choose the canonical
//          convention, that is IFloat[row][col][2].
//     2.2. Fuv returned in the canonical convention.
//     2.3. GetLink(..) gives links in the wilson order, that is
//          IFloat[col][row][2].
//     2.4. In order to use the matrix utilities in the system, 
//          we have to change the formulae into:
//            F(x, u, v) = 1/4 1/2 (NewPLAQs - NewPLAQs^dag)^Transpose
//          where NewPLAQs are
//           [ U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag +
//             U_nu(x) U_mu(x-mu+nu)^dag U_nu(x-mu)^dag U_mu(x-mu) +
//             U_mu(x-mu)^dag U_nu(x-mu-nu)^dag U_mu(x-mu-nu) U_nu(x-nu) +
//             U_nu(x-nu)^dag U_mu(x-nu) U_nu(x+mu-nu) U_mu(x)^dag 
//           ] ^Transpose
//  3. site[] has to be a legal site, and should not be off-node which is OK
//     only in calling GetLink.
//--------------------------------------------------------------------------
void 
DiracOpClover::SiteFuv(Matrix &Fuv, const int *site, int mu, int nu) const
{
  //  char *fname = "SiteFuv()";
  //  VRB.Func(cname,fname);

  // neighbour[] is the local coordinates of a neighbour site
  //------------------------------------------------------------------------
  int neighbor[4] = {site[0], site[1], site[2], site[3]};
  
  // mp0 = [ U_mu(x) U_nu(x+mu) ] ^Transpose
  //------------------------------------------------------------------------
  ++(neighbor[mu]);                               // neighbor = x+mu
  mp0->DotMEqual(GetLink(neighbor, nu), GetLink(site, mu));
 
  // mp1 = [ U_mu(x+nu)^dag U_nu(x)^dag ] ^Transpose
  //------------------------------------------------------------------------
  --(neighbor[mu]);  
  ++(neighbor[nu]);                               // neighbor = x+nu
  mp2->DotMEqual(GetLink(neighbor, mu), GetLink(site, nu));
  mp1->Dagger(*mp2);  
 
  // Fuv = [ U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag ] ^Transpose
  //------------------------------------------------------------------------
  Fuv.DotMEqual(*mp1, *mp0);
   
  // mp0 = [ U_nu(x) U_mu(x-mu+nu)^dag ] ^Transpose
  //------------------------------------------------------------------------
  --(neighbor[mu]);                             // neighbor = x-mu+nu
  mp2->Dagger(GetLink(neighbor, mu));
  mp0->DotMEqual(*mp2, GetLink(site, nu));
  
  // mp1 = [ U_nu(x-mu)^dag U_mu(x-mu) ] ^Transpose
  //------------------------------------------------------------------------
  --(neighbor[nu]);                             // neighbor = x-mu
  mp2->Dagger(GetLink(neighbor, nu));
  mp1->DotMEqual(GetLink(neighbor, mu), *mp2);
 
  // Fuv += [ U_nu(x) U_mu(x-mu+nu)^dag U_nu(x-mu)^dag U_mu(x-mu) ] ^Trans
  //------------------------------------------------------------------------
  Fuv.DotMPlus(*mp1, *mp0);

  // mp0 = [ U_mu(x-mu)^dag U_nu(x-mu-nu)^dag ] ^Trans
  //------------------------------------------------------------------------
  mp1->Dagger(GetLink(neighbor, mu));
  --(neighbor[nu]);                             // neighbor = x-mu-nu
  mp2->Dagger(GetLink(neighbor, nu));
  mp0->DotMEqual(*mp2, *mp1);
  
  // mp1 = [ U_mu(x-mu-nu) U_nu(x-nu) ] ^Transpose
  //------------------------------------------------------------------------
  *mp2 = GetLink(neighbor, mu);
  ++(neighbor[mu]);                             // neighbor = x-nu
  mp1->DotMEqual(GetLink(neighbor, nu), *mp2);
  
  // Fuv += [ U_mu(x-mu)^dag U_nu(x-mu-nu)^dag U_mu(x-mu-nu) U_nu(x-nu) ] ^T
  //------------------------------------------------------------------------
  Fuv.DotMPlus(*mp1, *mp0);
 
  // mp0 = [ U_nu(x-nu)^dag U_mu(x-nu) ] ^Transpose
  //------------------------------------------------------------------------
  mp2->Dagger(GetLink(neighbor, nu));
  mp0->DotMEqual(GetLink(neighbor, mu), *mp2);
 
  // mp1 = [ U_nu(x+mu-nu) U_mu(x)^dag ] ^Transpose
  //------------------------------------------------------------------------
  mp2->Dagger(GetLink(site, mu));
  ++(neighbor[mu]);                             // neighbor = x+mu-nu
  mp1->DotMEqual(*mp2, GetLink(neighbor, nu));
 
  // Fuv += [ U_nu(x-nu)^dag U_mu(x-nu) U_nu(x+mu-nu) U_mu(x)^dag ] ^Trans
  //------------------------------------------------------------------------
  Fuv.DotMPlus(*mp1, *mp0);
  
  // Fuv = (PLAQa - PLAQs^dag)/8
  //------------------------------------------------------------------------
  mp0->Dagger(Fuv);
  *mp0 -= Fuv;
  Fuv.Trans(*mp0);  
  Fuv *= -0.125;
}

//--------------------------------------------------------------------------
// local auxiliary function:  matrix *= i
//--------------------------------------------------------------------------
static void iDotM(Matrix &m, int num_complex) 
{
  for (int i = 0; i < num_complex; ++i) {
    Complex & c = m[i];
    IFloat new_im = c.real();
    IFloat new_re = - c.imag();
    c=Complex(new_re, new_im);
  }
}


//--------------------------------------------------------------------------
// void 
// DiracOpClover::SiteCloverMat(const int *site, IFloat *A, int inv) const
//--------------------------------------------------------------------------
// Purpose: 
//   Calculates the clover matrix at the specified site and stored it as
//   two consecutive compressed 6x6 hermitian matrices.
// Arguments:
//   site[0..3]          coordinate[x,y,z,t]
// Convension:
//   Links and fermionic fields are in Wilson order, that is
//   checkerboarded U[t][y][z][x][x,y,z,t].
//--------------------------------------------------------------------------
// Matrix A in Clover Action:
//    Fermion action S = Sum_x_y { phi(x)^bar M_xy phi(y) }
//       where M_xy = I - kappa D,     for Wilson fermions
//                  = A - kappa D,     for Clover improved Wilson fermions
//       where A is local and we adopt the same definition as in 
//       Xiang-Qian Luo's paper.
//   
//    A = I - kappa * Csw / 2 * Sum_u<v { [gamma_u, gamma_v] * Fuv } 
//       where Fuv (here) = i Fuv (Luo),  and we choose
//             v = 0..3 as x,y,z,t
//             gamma_x_y_z = 0  sigma_x_y_z         gamma_t = 0  I 
//                        -sigma_x_y_z  0                     I  0 
//             sigma_x = 0 i    sigma_y = 0 -1     sigma_z = i  0
//                       i 0              1  0               0 -i
//             [sigma_i, sigma_j] = 2 epsilon_i_j_k sigma_k 
//   ===>      [gamma_i, gamma_j] =-2 epsilon_i_j_k (sigma_k    0)
//                                                  (0    sigma_k)
//             [gamma_i, gamma_t] = 2               (sigma_i    0)
//                                                  (0   -sigma_i)
//    A = I + kappa * Csw * (termB - termE    0            )
//                          (0               termB + termE )
//       where termB =  epsilon_ijk sigma_k Fij
//                   =  sigma_z Fxy + sigma_x Fyz - sigma_y Fxz
//                   = (  i Fxy                i Fyz + Fxz )
//                     (  i Fyz - Fxz        - i Fxy       )
//                   = (  i F01                i F12 + F02 )
//                     (  i F12 - F02        - i F01       )
//                   = (  i F01          ___  dagger       )
//                     (  i F12 - F02 __/    - i F01       )
//             termE = sigma_x Fxt + sigma_y Fyt + sigma_z Fzt
//                   = ( i Fzt            i Fxt - Fyt )
//                     ( i Fxt + Fyt      - i Fzt     )
//                   = ( i F23            dagger      )
//                     ( i F03 + F13      - i F23     )
//    A = A0   A1^dag  0    0
//        A1   A2      0    0
//        0    0      A3    A4^dag
//        0    0      A4    A5
//        where  A0 = I + kappa * Csw * i(F01 - F23)
//               A1 =     kappa * Csw * [(iF12-F02) - (iF03+F13)]
//               A2 = I - kappa * Csw * i(F01 - F23)
//               A3 = I + kappa * Csw * i(F01 + F23)
//               A4 =     kappa * Csw * [(iF12-F02) + (iF03+F13)]
//               A5 = I - kappa * Csw * i(F01 + F23)
//--------------------------------------------------------------------------
void 
DiracOpClover::SiteCloverMat(const int *site, IFloat *A) const
{
  char *fname = "SiteCloverMat()";
  VRB.Func(cname,fname);

  // local buffer A0, ... , A5
  //------------------------------------------------------------------------
  Matrix As[2*COLORS]; 
  
  // As[1] = (i F12 - F02) omega
  //------------------------------------------------------------------------
  SiteFuv(As[1], site, 1, 2); 
  iDotM(As[1], HALF_LINK_MAT_SIZE);
  SiteFuv(*mp3, site, 0, 2); 
  As[1] -= *mp3;
  As[1] *= omega; 
 
  // As[4] = (i F03 + F13) omega_xi
  //------------------------------------------------------------------------
  SiteFuv(As[4], site, 0, 3);  
  iDotM(As[4], HALF_LINK_MAT_SIZE);
  SiteFuv(*mp3, site, 1, 3); 
  As[4] += *mp3;        
  As[4] *= omega_xi; 
                        
  // calcuate As[1] and As[4]
  //------------------------------------------------------------------------
  *mp1 = As[4];  
  As[4] += As[1];   
  As[1] -= *mp1;    
 
  //As[1] *= omega; 
  //As[4] *= omega_xi;
 
  // As[0] = i F01 omega 
  //------------------------------------------------------------------------
  SiteFuv(As[0], site, 0, 1); 
  iDotM(As[0], HALF_LINK_MAT_SIZE);  
  As[0] *= omega; 
  
  // As[3] = i F23 omega_xi
  //------------------------------------------------------------------------
  SiteFuv(As[3], site, 2, 3);  
  iDotM(As[3], HALF_LINK_MAT_SIZE);  
  As[3] *= omega_xi; 

 
  // As[0] = i omega (F01 - F23)  and As[3] = i omega (F01 + F23)
  //------------------------------------------------------------------------
  *mp1 = As[3];
  As[3] += As[0];   
  As[0] -= *mp1;   
 
  //As[0] *= omega;
  //As[3] *= omega_xi;
 
  // calcuate As[0], As[2], As[3] and As[5]
  //------------------------------------------------------------------------
  As[2] = 1.0; As[2] -= As[0];                         // As[2] done
  As[0] += 1.0;                                        // As[0] done
  As[5] = 1.0; As[5] -= As[3];                         // As[5] done
  As[3] += 1.0;                                        // As[3] done
 
  // fill in A with As[0], As[1], As[2], As[3], As[4] and As[5]
  // A is made of two 6x6 hermitian matrices, only the lower triangular
  // half matrix elements are stored to save memory.
  //------------------------------------------------------------------------
  {
    IFloat *dst = A;    
    for (int i = 0; i < 6; /* nothing, i=0or3 */) {
      { // fill in As[0] and As[3]
        IFloat *src = (IFloat *)(&(As[i++])); 
        for (int row = 0; row < COLORS; ++row) {
          for (int k = 0; k < 2*row+1; ++k)
            *dst++ = *src++;                 
          src += 2*(COLORS-row)-1;
        }
      }
      { // fill in As[1,2,4,5]
        IFloat *src1 = (IFloat *)(&(As[i++]));
        IFloat *src2 = (IFloat *)(&(As[i++]));
        for (int row = 0; row < COLORS; ++row) {
	  int k;
          for (k = 0; k < 2*COLORS; ++k)
            *dst++ = *src1++;                 // A[1], A[4] filled
          for (k = 0; k < 2*row+1; ++k)
            *dst++ = *src2++;                 // A[2], A[5] filled
          src2 += 2*(COLORS-row)-1;
        }
      }
    }
  }

  // The following implementation causes wrong data!!!!
  //------------------------------------------------------------------------
  /*
  {
    for (int i = 0; i < 2 ; ++i, A += HALF_CLOVER_MAT_SIZE) {                
      mat_hrm_cmpr(A, (IFloat *)(As+3*i), COLORS);       // A[0], A[3] filled
      IFloat *A1 = (IFloat *)(As+i*3+1);                  // A1 = A[1] & A[4]
      IFloat *A2 = (IFloat *)(As+i*3+2);                  // A2 = A[2] & A[5]
      IFloat *dst = A + HALF_LINK_MAT_SIZE;                    
      for (int row = 0; row < COLORS; ++row, A2 += 2*(COLORS-row)-1) {
	for (int k = 0; k < 2*COLORS; ++k)
	  *dst++ = *A1++;                                // A[1], A[4] filled
	for (k = 0; k < 2*row+1; ++k)
	  *dst++ = *A2++;                                // A[2], A[5] filled
      }
    }
  }
  */

  // The following implementation causes wrong data!!!!
  //------------------------------------------------------------------------
  /*
  {
    IFloat *AA = A;    
    for (int i = 0; i < 2 ; ++i, AA += HALF_CLOVER_MAT_SIZE) {                
      mat_hrm_cmpr(A, (IFloat *)(As+3*i), COLORS);      // A[0], A[3] filled
      IFloat *A1 = (IFloat *)(As+i*3+1);
      IFloat *A2 = (IFloat *)(As+i*3+2);
      IFloat *Ap = AA + HALF_LINK_MAT_SIZE;      
      for (int row = 0; row < COLORS; ++row) {
        int k;
        for (k = 0; k < 2*COLORS; ++k)
          *Ap++ = *A1++;                                // A[1], A[4] filled
        for (k = 0; k < 2*row+1; ++k)
          *Ap++ = *A2++;                                // A[2], A[5] filled
        A2 += 2*(COLORS-row)-1;
      }
    }
  }
  */  
}

//--------------------------------------------------------------------------
// void 
// DiracOpClover::CloverMatChkb(ChkbType chkb, int inv) const
//--------------------------------------------------------------------------
// Purpose: 
//   Calculates the clover matrix for all sites of the specified checkerboard
//   and stores them in the memory pointed by lat.Aux0Ptr() for
//   chkb CHKB_EVEN and by lat.Aux1Ptr() for chkb CHKB_ODD.
// Arguments:
//   site[0..3]          coordinate[x,y,z,t]
// Convension:
//   Links and fermionic fields are in Wilson order, that is
//   checkerboarded U[t][y][z][x][x,y,z,t].
//--------------------------------------------------------------------------
void 
DiracOpClover::CloverMatChkb(ChkbType chkb, int inv) const
{
  char *fname = "CloverMatChkb()";
  VRB.Func(cname,fname);

  // get the pointer to the clover matrices.
  //------------------------------------------------------------------------
  IFloat *cl_mat_p = (IFloat *)(chkb == CHKB_EVEN ? 
			      lat.Aux0Ptr() : lat.Aux1Ptr());
  
  // loop in the order consistent with the Wilson checkerboarded storage
  // Note: the sites[0,1,2,3] is in order of [x,y,z,t] 
  //------------------------------------------------------------------------
  {
    IFloat *Ap = cl_mat_p;  
    const int parity = (chkb == CHKB_ODD ? 1 : 0);
    const int *nsites = clover_lib_arg->nsites;    
    int site[4];                                   
    for (site[3] = 0; site[3] < nsites[3]; ++(site[3])) {
      for (site[2] = 0; site[2] < nsites[2]; ++(site[2])) {
	for (site[1] = 0; site[1] < nsites[1]; ++(site[1])) {
	  site[0] = (site[3] + site[2] + site[1] + parity)%2;	
	  for (; site[0] < nsites[0]; site[0] += 2) {
	    SiteCloverMat(site, Ap);	  
	    Ap += CLOVER_MAT_SIZE;     
	  }
	}
      }
    }
  }


  // invert the clover matrix
  //------------------------------------------------------------------------
  if (inv) {
    IFloat *Ap = cl_mat_p; 
    for (int i = 0; i < GJP.VolNodeSites(); ++i) {
      mat_inv(Ap, Ap, 6, MAT_INV_ALG_LDL_CMPR, 0);    
      Ap += HALF_CLOVER_MAT_SIZE;         
    } 
  }
}









// FOR THE PURPOSE OF DEBUGGING
CPS_END_NAMESPACE
#include <stdio.h>
#include<util/smalloc.h>
CPS_START_NAMESPACE

//----------------------------------------------------------------------
//
//----------------------------------------------------------------------
void DiracOpClover::MatDagMatDbg(Vector *out, 
				 Vector *in, 
				 Float *dot_prd, 
				 int direct) 
{
  char *fname = "MatDagMatDbg(V*,V*,F*)";
  VRB.Func(cname,fname);

  // Have to allocate new buffer, double the size of frm_buf0/1
  //--------------------------------------------------------------------
  int vec_size = lat.FsiteSize()*GJP.VolNodeSites();  
  Vector *temp = (Vector *) smalloc(vec_size * sizeof(IFloat));

  // temp = MatPc in
  //--------------------------------------------------------------------
  MatDagOrNotDbg(temp, in, 0, direct);

  // *dot_prd = <temp, temp>
  //--------------------------------------------------------------------
  if (dot_prd) {
    *dot_prd = temp->NormSqNode(vec_size);    
  }

  // out = MatPcDag MatPc in
  //--------------------------------------------------------------------
  MatDagOrNotDbg(out, temp, 1, direct);
  sfree(temp);  
}

//----------------------------------------------------------------------
// Mat(Vector *out, Vector *in) :
//----------------------------------------------------------------------
// Buffer:
//      use clover_lib_arg->frm_buf1,
//      call MatPc in which clover_lib_arg->frm_buf0 will be used.
//----------------------------------------------------------------------
void DiracOpClover::MatDagOrNotDbg(Vector *out, const Vector *in, 
				   int dag, int direct) const
{
  const int half_sites = GJP.VolNodeSites()/2;
  const int vec_size = lat.FsiteSize() * half_sites;    
  IFloat *in_even = (IFloat *)in + vec_size;
  IFloat *out_even = (IFloat *)out + vec_size;  
  IFloat *A_even = (IFloat *)lat.Aux0Ptr();
  IFloat *A_odd = (IFloat *)lat.Aux1Ptr();
  Vector *frm_buf1 = (Vector *)(clover_lib_arg->frm_buf1);  
  Wilson* wilson_p = clover_lib_arg->wilson_p;
  
  // get out_even = Aee in_even - kappa Deo in_odd
  //---------------------------------------------------------------------
  CloverMatChkb(CHKB_EVEN, 0);
  clover_mat_mlt((IFloat *)out_even, A_even, (IFloat *)in_even, half_sites);
  CloverMatChkb(CHKB_EVEN, 1);
  wilson_dslash((IFloat *)frm_buf1, (IFloat *)gauge_field,  
  		(IFloat *)in, 1, dag, wilson_p);         
  fTimesV1PlusV2((IFloat *)out_even, -kappa, (const IFloat *)frm_buf1, 
    		 (const IFloat *)out_even, vec_size);   

 
  // get out_odd = Aoo in_odd - kappa Doe in_even                 directly
  //      OR     = MatPc(in_odd) - kappa Doe Aee^inv out_even   indirectly
  //---------------------------------------------------------------------
  if (direct) {   // no preconditioning
    clover_mat_mlt((IFloat *)out, A_odd, (IFloat *)in, half_sites);
    wilson_dslash((IFloat *)frm_buf1, (IFloat *)gauge_field,  
		  (IFloat *)in_even, 0, dag, wilson_p);         
    fTimesV1PlusV2((IFloat *)out, -kappa, (const IFloat *)frm_buf1, 
		   (const IFloat *)out, vec_size);   
  } else {        // using even-odd preconditioning
    clover_mat_mlt((IFloat *)frm_buf1, A_even, 
		   (IFloat *)out_even, half_sites);
    wilson_dslash((IFloat *)out, (IFloat *)gauge_field,  
		  (IFloat *)frm_buf1, 0, dag, wilson_p);         
    MatPcDagOrNot(frm_buf1, in, dag);
    fTimesV1PlusV2((IFloat *)out, -kappa, (const IFloat *)out, 
		   (const IFloat*)frm_buf1, vec_size);   
  }
}






CPS_END_NAMESPACE
