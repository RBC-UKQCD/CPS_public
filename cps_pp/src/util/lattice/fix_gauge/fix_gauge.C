#include<config.h>
CPS_START_NAMESPACE
/*!\file
  $Id: fix_gauge.C,v 1.9.30.1.6.4 2013/02/19 22:32:04 yinnht Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: yinnht $
//  $Date: 2013/02/19 22:32:04 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/fix_gauge/fix_gauge.C,v 1.9.30.1.6.4 2013/02/19 22:32:04 yinnht Exp $
//  $Id: fix_gauge.C,v 1.9.30.1.6.4 2013/02/19 22:32:04 yinnht Exp $
//  $Name:  $
//  $Locker:  $
//  $Revision: 1.9.30.1.6.4 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/fix_gauge/fix_gauge.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

/*--------------------------------------------------------------------
 * File fix_gauge.C. Version 14.5. Last modified on 98/11/25 at 19:50:37.
 *        Yuriy Zhestkov, Columbia University.
 * E-mail: zhestkov@phys.columbia.edu
 *--------------------------------------------------------------------
 */

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/smalloc.h>
#include <util/rcomplex.h>
#include <util/time_cps.h>
#include <comms/scu.h>
#include <comms/glb.h>
#include <omp.h>
CPS_START_NAMESPACE

#ifdef _TARTAN
CPS_END_NAMESPACE
#include <math.h>
CPS_START_NAMESPACE
#else
CPS_END_NAMESPACE
#include <math.h>
CPS_START_NAMESPACE
#endif

//! The number of colours, again.
#define COLORS 3  // implemented only for this number of colors

//------------------------------------------------------------------------//
// class XXX - gathers functions that hopefully will be replaced with     //
//             system calls                                               //
//------------------------------------------------------------------------//

typedef unsigned uword;  // DSP word is 32 bits

//! A container class for global parameters used in the gauge fixing routines.
class XXX
{
public:

    static const int CheckFreq;
    static const int Dimension;
    static const int SiteSize;

public:

    static int  Coor4d(int dir);

    static int  Node_Size_in_Dir(int dir);
    static int  Num_Nodes_in_Dir(int dir);
};

const int XXX::CheckFreq = 10; // check stop condition every CHECKFREQ's iter
const int XXX::Dimension = 4;
const int XXX::SiteSize = XXX::Dimension;

int XXX::Coor4d(int dir)
{
    int coor[] = {GJP.XnodeCoor(), GJP.YnodeCoor(), GJP.ZnodeCoor(),
                  GJP.TnodeCoor()};
    return coor[dir];
}

int XXX::Node_Size_in_Dir(int dir)
{
    int node_sites[] = {GJP.XnodeSites(), GJP.YnodeSites(), 
                        GJP.ZnodeSites(), GJP.TnodeSites()};
    return node_sites[dir];
}

int XXX::Num_Nodes_in_Dir(int dir)
{
    int machine_size[] = {GJP.Xnodes(),GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};
    return machine_size[dir];
}



//------------------------------------------------------------------------//
//                                                                        //
//  class HyperPlane contains data for a subspace (hyperplane) of the     //
//  lattice the gauge fixing matrices for which has to be found.          //
//  When there are several hyperplanes to be fixed, they can be treated   //
//  separately, provided that the links themselve are not modified.       //
//                                                                        //
//------------------------------------------------------------------------//

class HyperPlane
{
private:

    char *cname;

protected:

    HyperPlane(const Matrix *L0, int Ind2Dir[], 
               int HplDim, Matrix *Gauge, Float SmallFloat);
    virtual ~HyperPlane();

    IFloat HyperPlane_Sum(IFloat sum_1);

protected:

    // L_loc(int ind_link) returns the specified link if it is local 
    // (not in a buffer).
    // L(int ind_link) returns local link or a copy from a buffer.
    // Used always to access links except one place in
    // HyperPlane::copy_nbr_links(..)
    // The link argument - in index order (so, "real link #"=ind2dir[ind_link])

    const Matrix& L_loc(int ind_link) const;
    const Matrix& L(int ind_link) const;

    // G(..), G_loc() return a gauge transformation matrix
    // G(..) might initiate transmission to get those from neighbor nodes

    Matrix& G_loc();
    Matrix& G(int ind_link); // uses g_buf!

private:

    Matrix g_buf;         // used in G() to store received matrix

protected:

    int  hplane_dim;    // dimension of the hyperplane
    int *ind2dir;       // index field - physical direction correspondence
    int *stride;        // address increment for different index fields
    int *node_size;     // number of sites in index direction
    IFloat small_enough; // for the stop condition

    // index coordinates of site (used to avoid argument passing)
    int *index;

    //Add by Jianglei
    Matrix **nbr_gauge_minus_rec;
    Matrix **nbr_gauge_plus_rec;
    Matrix **nbr_gauge_minus_send;
    Matrix **nbr_gauge_plus_send;
    void copy_nbr_gauge();

private:

    Matrix *gauge;          // pointer to a hyperplane of gauge fixing matrices
    const Matrix *l0;       // pointer to the first element of the hyperplane
    Matrix **neighbor_link; // copies of neighbor node links

    // copy the links necessary for gauge fixing calculations from the 
    // neighbor nodes into the local buffers neighbor_link[]
    void copy_nbr_links(int nbr_num, int recurse = 0);
};



//-------------------------------------------------------------------------
//
//  constructor for class HyperPlane
//
//-------------------------------------------------------------------------

HyperPlane::HyperPlane(const Matrix*L0, int Ind2Dir[], 
		       int HplDim, Matrix *Gauge, Float SmallFloat):
    hplane_dim(HplDim),
    gauge(Gauge),
    l0(L0)
{
    cname = "HyperPlane";
    char *fname = "HyperPlane(M*,i*,i,M*,F)";

    VRB.Func(cname,fname);

    //-----------------------------------------------------------------------
    // create a buffer for passing current index in recursive function calls
    //-----------------------------------------------------------------------
    {
        int mem_size = hplane_dim * sizeof(int);
        index = (int*) smalloc(mem_size);
        if(index == NULL)
            ERR.Pointer(cname, fname, "index");
        VRB.Smalloc(cname, fname, "index", index, mem_size);
    }

    //-----------------------------------------------------------------------
    // copy the index_in_hyperplane -- to -- physical_direction information
    //-----------------------------------------------------------------------
    {
        int mem_size = hplane_dim * sizeof(int);
        ind2dir = (int*) smalloc(mem_size);
        if(ind2dir == NULL)
            ERR.Pointer(cname, fname, "ind2dir");
        VRB.Smalloc(cname, fname, "ind2dir", ind2dir, mem_size);
        for(int i=0; i<hplane_dim; i++)
            ind2dir[i] = Ind2Dir[i];
    }
  
    //-----------------------------------------------------------------------
    // calculating distance in address space corresponding to index
    //-----------------------------------------------------------------------
    {
        int mem_size = hplane_dim * sizeof(int);
        stride = (int*) smalloc(mem_size);
        if(stride == NULL)
            ERR.Pointer(cname, fname, "stride");
        VRB.Smalloc(cname, fname, "stride", stride, mem_size);
        for(int i=0; i<hplane_dim; i++)
            {
                stride[i] = XXX::SiteSize;
                for(int j=0; j<ind2dir[i]; j++)
                    stride[i] *= XXX::Node_Size_in_Dir(j);
            }
    }
  
    //-----------------------------------------------------------------------
    // assigning number of sites in the node for index directions
    //-----------------------------------------------------------------------
    {
        int mem_size = hplane_dim * sizeof(int);
        node_size = (int*) smalloc(mem_size);
        if(node_size == NULL)
            ERR.Pointer(cname, fname, "node_size");
        VRB.Smalloc(cname, fname, "node_size", node_size, mem_size);
        for(int i=0; i<hplane_dim; i++)
            node_size[i] = XXX::Node_Size_in_Dir(ind2dir[i]);
    }

    //-----------------------------------------------------------------------
    // find small_enough
    //-----------------------------------------------------------------------
    {
        int num_reals = 2*COLORS*COLORS;
        for(int i=0; i<hplane_dim; i++)
            num_reals *= 
                XXX::Num_Nodes_in_Dir(ind2dir[i]) * node_size[i];
        small_enough = SmallFloat * num_reals;
    }

    //-----------------------------------------------------------------------
    // create and initialize buffers for copies of neighbor node links
    //-----------------------------------------------------------------------
    {
        int mem_size = hplane_dim * sizeof(Matrix*);
        neighbor_link = (Matrix**) smalloc(mem_size);
        if(neighbor_link == NULL)
            ERR.Pointer(cname, fname, "neighbor_link");
        VRB.Smalloc(cname, fname, "neighbor_link", neighbor_link, mem_size);

        for(int i=0; i<hplane_dim; i++)
            {
                mem_size = sizeof(Matrix);
                for(int j=0; j<hplane_dim; j++)
                    mem_size *= (i==j) ? 1 : node_size[j];

                neighbor_link[i] = (Matrix*) smalloc(mem_size);
                if(neighbor_link[i] == NULL)
                    ERR.Pointer(cname, fname, "neighbor_link[i]");
                VRB.Smalloc(cname,fname,"neighbor_link[i]",neighbor_link[i],mem_size);

                copy_nbr_links(i);
            }
    }
    //-------------------------------------------------------------------------------
    //creat and initialize buffers for copies of neighbor node gaueg matrices
    //, add by Jianglei
    //------------------------------------------------------------------------------
    {
        int mem_size = hplane_dim * sizeof(Matrix*);
        nbr_gauge_minus_rec = (Matrix**) smalloc(mem_size);
        nbr_gauge_minus_send = (Matrix**) smalloc(mem_size);
        nbr_gauge_plus_rec = (Matrix**) smalloc(mem_size);
        nbr_gauge_plus_send = (Matrix**) smalloc(mem_size);

        for(int i=0; i<hplane_dim; i++)
            {
                mem_size = sizeof(Matrix);
                for(int j=0; j<hplane_dim; j++)
                    mem_size *= (i==j) ? 1 : node_size[j];

                nbr_gauge_minus_rec[i] = (Matrix*) smalloc(mem_size);
                nbr_gauge_minus_send[i] = (Matrix*) smalloc(mem_size);
                nbr_gauge_plus_rec[i] = (Matrix*) smalloc(mem_size);
                nbr_gauge_plus_send[i] = (Matrix*) smalloc(mem_size);
            }
        copy_nbr_gauge();
    }
}



//-------------------------------------------------------------------------
//
//  destructor for class HyperPlane
//
//-------------------------------------------------------------------------

HyperPlane::~HyperPlane()
{
    char *fname = "~HyperPlane()";

    VRB.Func(cname, fname);

    for(int i=hplane_dim; i--; )
        {
            VRB.Sfree(cname, fname, "neighbor_link[i]", neighbor_link[i]);
            sfree(neighbor_link[i]);
        }

    VRB.Sfree(cname, fname, "neighbor_link", neighbor_link);
    sfree(neighbor_link);

    for(int i=hplane_dim; i--; )
	{
            sfree(nbr_gauge_minus_rec[i]);
            sfree(nbr_gauge_minus_send[i]);
            sfree(nbr_gauge_plus_rec[i]);
            sfree(nbr_gauge_plus_send[i]);
	}
    sfree(nbr_gauge_minus_rec);
    sfree(nbr_gauge_minus_send);
    sfree(nbr_gauge_plus_rec);
    sfree(nbr_gauge_plus_send);

    VRB.Sfree(cname, fname, "node_size", node_size);
    sfree(node_size);

    VRB.Sfree(cname, fname, "stride", stride);
    sfree(stride);

    VRB.Sfree(cname, fname, "ind2dir", ind2dir);
    sfree(ind2dir);

    VRB.Sfree(cname, fname, "index", index);
    sfree(index);
}



//-------------------------------------------------------------------------
//
//  semi-global sum over a hyperplane
//
//-------------------------------------------------------------------------

IFloat HyperPlane::HyperPlane_Sum(IFloat sum_1)
{
    char *fname = "HyperPlane_Sum(f)";

    VRB.Func(cname, fname);

    IFloat sum = sum_1;

    //-----------------------------------------------------------------------
    // propagate and sum over all hyperplane nodes
    //-----------------------------------------------------------------------
    {
        IFloat tmp;
        for(int i=0; i<hplane_dim; i++)
            {
                sum_1 = sum; 
                for(int j=0; j<XXX::Num_Nodes_in_Dir(ind2dir[i])-1; j++)
                    {
                        tmp = sum_1;
                        getPlusData(&sum_1, &tmp, 1 , ind2dir[i]);
                        sum += sum_1;
                    }
            }
    }

    //-----------------------------------------------------------------------
    // propagate and maximize over hyperplane nodes (to ensure same result)
    //-----------------------------------------------------------------------
    {
        IFloat tmp;
        for(int i=0; i<hplane_dim; i++)
            {
                for(int j=0; j<XXX::Num_Nodes_in_Dir(ind2dir[i])-1; j++)
                    {
                        getPlusData(&tmp, &sum_1, 1, ind2dir[i]);
                        sum = (sum>tmp) ? sum : tmp;
                    }
            }
    }
  
    return sum;
}


//-------------------------------------------------------------------------
//                                                        *) sets index[] !
// copy links from neighbor nodes into neighbor_link[]
//
//-------------------------------------------------------------------------

void HyperPlane::copy_nbr_links(int nbr_num, int recurse)
{
    char *fname = "copy_nbr_links(i,i)";

    VRB.Func(cname, fname);

    if(recurse < hplane_dim) // index has not been completed
        if(recurse == nbr_num)
            {
                copy_nbr_links(nbr_num, recurse + 1); // skip this field in index
            }
        else
            {
                for(int i=node_size[recurse]; i--; )
                    {
                        index[recurse] = i;
                        copy_nbr_links(nbr_num, recurse + 1);
                    }
            }
    else  // index is completed
        {
            // send link # nbr_num via wire in positive nbr_num direction

            index[nbr_num] = - 1;
            void * rec = (void*) &L(nbr_num);

            index[nbr_num] = node_size[nbr_num] - 1;
            getMinusData((IFloat*)rec, (IFloat*)&L_loc(nbr_num),
                         sizeof(Matrix)/sizeof(IFloat), ind2dir[nbr_num]);
        }
}

//Coulumb Gauge only !!!
void HyperPlane::copy_nbr_gauge()
{
    int ind1, ind2;

    for(int i = 0; i < node_size[1]; i++)	
        for(int j = 0; j < node_size[2]; j++)	
            {
		ind1 = i + j * node_size[1];
		ind2 = node_size[0] - 1 + node_size[0] * (i + j * node_size[1]);
		nbr_gauge_minus_send[0][ind1] = gauge[ind2];
		ind2 = 0 + node_size[0] * (i + j * node_size[1]);
		nbr_gauge_plus_send[0][ind1] = gauge[ind2];
            }

    for(int i = 0; i < node_size[0]; i++)	
        for(int j = 0; j < node_size[2]; j++)	
            {
		ind1 = i + j * node_size[0];
		ind2 = i + node_size[0] * ( node_size[1] - 1 + node_size[1] * j );
		nbr_gauge_minus_send[1][ind1] = gauge[ind2];
		ind2 = i + node_size[0] * ( 0 + node_size[1] * j );
		nbr_gauge_plus_send[1][ind1] = gauge[ind2];
            }

    for(int i = 0; i < node_size[0]; i++)	
        for(int j = 0; j < node_size[1]; j++)	
            {
		ind1 = i + j * node_size[0];
		ind2 = i + node_size[0] * ( j + node_size[1] * (node_size[2] - 1));
		nbr_gauge_minus_send[2][ind1] = gauge[ind2];
		ind2 = i + node_size[0] * ( j + node_size[1] * (0));
		nbr_gauge_plus_send[2][ind1] = gauge[ind2];
            }
	
    for(int i = 0; i < 3; i++)
	{
            int memsize = sizeof(Matrix) / sizeof(IFloat);
            for(int j = 0; j < 3; j++) if(j != i) memsize *= node_size[j];
            getMinusData((IFloat*)nbr_gauge_minus_rec[i], (IFloat*)nbr_gauge_minus_send[i], memsize, ind2dir[i]);
            getPlusData((IFloat*)nbr_gauge_plus_rec[i], (IFloat*)nbr_gauge_plus_send[i], memsize, ind2dir[i]);
	}
}	


//-------------------------------------------------------------------------
//
// CRITICAL PERFORMANCE CODE !!!
// access to a local link
//
//-------------------------------------------------------------------------

const Matrix& HyperPlane::L_loc(int ind_link) const
{
    int idx = ind2dir[ind_link];
    for(int i=0; i<hplane_dim; i++)
        idx += index[i]*stride[i];
  
    return l0[idx];
}



//-------------------------------------------------------------------------
//
// CRITICAL PERFORMANCE CODE !!!
// access to either a local link or a copy in neighbor_link
//
//-------------------------------------------------------------------------

const Matrix& HyperPlane::L(int ind_link) const
{
    if(index[ind_link] == -1) // link not on the node - redirect to a link buffer
        {
            int idx = 0;
            for(int i=hplane_dim; i--; )
                if(i != ind_link) // skip index field related to the neighbor direction
                    {
                        idx *= node_size[i];
                        idx += index[i];
                    }
            return neighbor_link[ind_link][idx];
        }
    else
        return L_loc(ind_link);
}



//-------------------------------------------------------------------------
//
// CRITICAL PERFORMANCE CODE !!!
// access to a local gauge fixing matrix
//
//-------------------------------------------------------------------------

Matrix& HyperPlane::G_loc()
{
    int idx = index[hplane_dim-1];
    for(int i=hplane_dim-1; i-- > 0 ; )
        {
            idx *= node_size[i];
            idx += index[i];
        }

    return gauge[idx];  
}


//-------------------------------------------------------------------------
//
// CRITICAL PERFORMANCE CODE !!!
// access to either a local gauge fixing matrix or transmitted copy from
// a neighbor node.
//
//-------------------------------------------------------------------------

Matrix& HyperPlane::G(int ind_link)
{
    if(index[ind_link] == -1)
	{
            //index[ind_link] = node_size[ind_link] - 1;
            //getMinusData((IFloat*)&g_buf, (IFloat*)&G_loc(), 
            // sizeof(Matrix)/sizeof(IFloat), ind2dir[ind_link]);
            //index[ind_link] = -1;
            //return g_buf;
            int idx = 0;
            for(int i=hplane_dim; i--; )
                if(i != ind_link) // skip index field related to the neighbor direction
                    {
                        idx *= node_size[i];
                        idx += index[i];
                    }
            return nbr_gauge_minus_rec[ind_link][idx];
	}
    else if(index[ind_link] == node_size[ind_link])
	{
            //index[ind_link] = 0;
            //getPlusData((IFloat*)&g_buf, (IFloat*)&G_loc(), 
            //		sizeof(Matrix)/sizeof(IFloat), ind2dir[ind_link]);
            //index[ind_link] = node_size[ind_link];
            //return g_buf;
            int idx = 0;
            for(int i=hplane_dim; i--; )
                if(i != ind_link) // skip index field related to the neighbor direction
                    {
                        idx *= node_size[i];
                        idx += index[i];
                    }
            return nbr_gauge_plus_rec[ind_link][idx];
	}
    else
	{
            return G_loc();
	}
}



//-------------------------------------------------------------------------//
//                                                                         //
// class FixHPlane - contains functions for finding the gauge fixing       //
// matrices                                                                //
//                                                                         //
//-------------------------------------------------------------------------//

class FixHPlane: public HyperPlane
{
private:

    char *cname;

protected:

    FixHPlane(const Matrix *L0, int Ind2Dir[], 
              int HplDim, Matrix *Gauge, Float SmallFloat);
    ~FixHPlane();

    void iter();
    IFloat delta() { return delta(0); }
    void unitarize(int recurse = 0);

protected:

    void iter(int recurse, int dist);
    IFloat delta(int recurse);

private:

    Matrix A_buf;        // used in findA()
    Matrix& findA();     // stores result in A_buf;
    void fix_g(Matrix&);

protected:

    Matrix tmp_m, tmp_m1; // for temporary local scope usage

private:
  
    int *dist_max;

    friend int Lattice::FixGauge(Float SmallFloat, int MaxIterNum);
};



//-------------------------------------------------------------------------
//
//  constructor for class FixHPlane
//
//-------------------------------------------------------------------------

FixHPlane::FixHPlane(const Matrix *L0, int Ind2Dir[],
		     int HplDim, Matrix *Gauge, Float SmallFloat) :
    HyperPlane(L0, Ind2Dir, HplDim, Gauge, SmallFloat)
{
    cname = "FixHPlane";
    char *fname = "FixHPlane(M*,i*,i,M*,F)";

    VRB.Func(cname, fname);

    //-----------------------------------------------------------------------
    // dist_max[i] - range of order i "distances"
    //-----------------------------------------------------------------------
    {
        int mem_size = hplane_dim * sizeof(int);
        dist_max = (int*) smalloc(mem_size);
        if(dist_max == NULL)
            ERR.Pointer(cname, fname, "dist_max");
        VRB.Smalloc(cname, fname, "dist_max", dist_max, mem_size);

        dist_max[0] = node_size[0] - 1;
        for(int i = 1; i<hplane_dim; i++)
            dist_max[i] = dist_max[i-1] + node_size[i] - 1;
    }
}


//-------------------------------------------------------------------------
//
//  destructor for class FixHPlane
//
//-------------------------------------------------------------------------

FixHPlane::~FixHPlane()
{
    char *fname = "~FixHPlane()";

    VRB.Func(cname, fname);

    VRB.Sfree(cname, fname, "dist_max", dist_max);
    sfree(dist_max);
}



//-------------------------------------------------------------------------
//
//  upper level iteration routine
//
//-------------------------------------------------------------------------

void FixHPlane::iter()
{
    char *fname = "iter()";

    VRB.Func(cname, fname);
  
    int recurse = hplane_dim - 1;
    //----------------------------------------------------------------------
    // I divided distance loop to even and odd part.
    //
    // Even or Odd is defined by
    // mod( index[0]+...+index[hplane_dim-1], 2 ) = 0 or 1
    //                                        Takeshi Yamazaki
    //----------------------------------------------------------------------
    if (GJP.GfixChkb()==0) {
        ERR.General(cname,fname,"Not implemented, you must turn on the checkboard flag!\n");
        exit(-1);
        // VRB.Flow(cname,fname,"using sequential order gauge fixing");
        for(int dist = dist_max[recurse] + 1; dist-- > 0; )
            iter(recurse, dist);
    } else {
        VRB.Flow(cname,fname,"using checkerboared order gauge fixing");
        //Even or Odd part
        for(int dist = dist_max[recurse] + 1; (dist=dist-2) >= 0; )
            iter(recurse, dist);
        copy_nbr_gauge();
        //Odd or Even part
        for(int dist = dist_max[recurse] + 2; (dist=dist-2) >= 0; )
            iter(recurse, dist);
        copy_nbr_gauge();
    }

}


//-------------------------------------------------------------------------
//                                                        *) sets index[] !
//  recursive iteration routine
//
//-------------------------------------------------------------------------

void FixHPlane::iter(int recurse, int dist)
{
    char *fname = "iter(i,i)";

    VRB.Func(cname, fname);

    // "distance" from origin = x[0]+x[1]+...+x[recurse]

    int dm = (recurse>0) ? dist_max[recurse-1] : 0;
    int xmin = (dist > dm) ? dist-dm : 0;
    int s = node_size[recurse] - 1;
    int xmax = (dist < s) ? dist : s;

    for(int xx = xmin; xx <= xmax; xx++)
        {
            index[recurse] = xx;
            if(recurse > 0)
                iter(recurse-1, dist-xx);
            else
                /*	{
                        int go_sign = 0;
                        for(int dim=0; dim < hplane_dim ; dim++) go_sign+=index[dim];
                        go_sign = go_sign % 2;
                        printf("(%d,%d,%d) = %d\n",index[0],index[1],index[2],go_sign);
                */
                fix_g(G_loc());
            //}
        }
}


//-------------------------------------------------------------------------
//
// ATTENTION!
// Here conventions different from those in Dong's thesis are used.
//
// This function returnes matrix A(n) given by
//
//               3   +        
//       A(n) = sum[Lj(n) + Lj(n-j)],
//               j
// where according to new conventions
//                            +
//       Lj(n) = g(n) L0j(n) g(n+j)
//
// It used to be
//               3           +
//       A(n) = sum[Lj(n) + Lj(n-j)], (see Zhihua Dong's thesis, p. 29)
//               j
//                                  +
//       and Lj(n) = g(n+j) L0j(n) g(n)
//
//-------------------------------------------------------------------------

Matrix& FixHPlane::findA()
{
    char *fname = "findA()";

    VRB.Func(cname, fname);

    tmp_m.ZeroMatrix();

    for(int i=0; i<hplane_dim; i++)
        {
            index[i]--;

            // Modified for anisotropic lattices

            if ( ind2dir[i] == GJP.XiDir()) {
                tmp_m1 = L(i);
                tmp_m1 *= GJP.XiGfix();
                tmp_m.DotMPlus(G(i), tmp_m1);
            } else {
                tmp_m.DotMPlus(G(i), L(i));
            }
            // End modification

            index[i] += 2;
            Matrix& G_plus_i = G(i);

            index[i]--;
            tmp_m1.Dagger(L_loc(i)); // local link

            // Added for anisotropic lattices
            if ( ind2dir[i] == GJP.XiDir() ) tmp_m1 *= GJP.XiGfix();
            // End modification

            tmp_m.DotMPlus(G_plus_i, tmp_m1);

        }

    tmp_m1.Dagger(G_loc());  // local G
    A_buf.DotMEqual(tmp_m, tmp_m1);
    return A_buf;
}


//-------------------------------------------------------------------------
//
// find a g+(n) matrix which maximizes the quantity Re Tr A*g+(n)
// see Zhihua Dong's thesis, p. 29 
// Then apply it to g
//
// See also description of findA() for important changes
//
//-------------------------------------------------------------------------

void FixHPlane::fix_g(Matrix& g)
{
    char *fname = "fix_g(M&)";

    VRB.Func(cname, fname);

    //-------------------------------------------------------------------
    // calculate A
    //-------------------------------------------------------------------

    Matrix& A = findA();

    //-------------------------------------------------------------------
    // su2upleft
    //-------------------------------------------------------------------
    {
        Complex a(A[0]),       b(A[1]),
            c(A[COLORS]),  d(A[COLORS+1]);

        Complex v00(real(a) + real(d), imag(d) - imag(a)), 
            v01(real(c) - real(b), -(imag(c) + imag(b)));

#ifdef _TARTAN
        Float sdet = double( sqrt(norm(v00) + norm(v01)) );
#else
        Float sdet = sqrt(norm(v00) + norm(v01));
#endif


        v00 /= sdet; v01 /= sdet;

        Complex v10(-real(v01), imag(v01)), v11(real(v00), -imag(v00));
  
        A[0] = a*v00 + b*v10;        A[1] = a*v01 + b*v11;
        A[COLORS+1] = c*v01 + d*v11;

        Complex e(A[2*COLORS]), f(A[2*COLORS+1]);
        A[2*COLORS] = e*v00 + f*v10; A[2*COLORS+1] = e*v01 + f*v11;

        tmp_m[0] = v11;       tmp_m[1] = -v01;
        tmp_m[COLORS] = -v10; tmp_m[COLORS+1] = v00;
        tmp_m[2] = tmp_m[COLORS+2] 
            = tmp_m[2*COLORS] = tmp_m[2*COLORS+1] = 0.0;
        tmp_m[2*COLORS+2] = 1.0;
    }

    //-------------------------------------------------------------------
    // su2downright
    //-------------------------------------------------------------------
    {
        Complex c(A[2*COLORS+1]), d(A[2*COLORS+2]);

        Complex v00(real(A[COLORS+1]) + real(d), imag(d) - imag(A[COLORS+1])), 
            v01(real(c) - real(A[COLORS+2]), -(imag(c) + imag(A[COLORS+2])));

#ifdef _TARTAN
        Float sdet = double( sqrt(norm(v00) + norm(v01)) );
#else
        Float sdet = sqrt(norm(v00) + norm(v01));
#endif

        v00 /= sdet; v01 /= sdet;

        Complex v10(real(v01), -imag(v01)), v11(real(v00), -imag(v00));
  
        A[2*COLORS+2] = c*v01 + d*v11;
        A[2] = A[1]*v01 + A[2]*v11;

        tmp_m[2*COLORS] = tmp_m[COLORS]*v10; tmp_m[COLORS] *= v11; 
        tmp_m[2*COLORS+1] = tmp_m[COLORS+1]*v10; 
        tmp_m[COLORS+1] *= v11;
        tmp_m[COLORS+2] = -v01; tmp_m[2*COLORS+2] = v00;
    }

    //-------------------------------------------------------------------
    // apply this intermediate transformation
    //-------------------------------------------------------------------

    tmp_m1.DotMEqual(tmp_m, g);

    //-------------------------------------------------------------------
    // su2corners
    //-------------------------------------------------------------------
    {
        tmp_m[0] = Complex(real(A[0]) + real(A[2*COLORS+2]), 
                           imag(A[0]) - imag(A[2*COLORS+2]));
        tmp_m[2*COLORS] = Complex(real(A[2*COLORS]) - real(A[2]), 
                                  imag(A[2*COLORS]) + imag(A[2]));

#ifdef _TARTAN
        Float sdet = double( sqrt(norm(tmp_m[0]) + norm(tmp_m[2*COLORS])) );
#else
        Float sdet = sqrt(norm(tmp_m[0]) + norm(tmp_m[2*COLORS]));
#endif
        tmp_m[0] /= sdet; tmp_m[2*COLORS] /= sdet;
        tmp_m[2] = Complex(-real(tmp_m[2*COLORS]), imag(tmp_m[2*COLORS]));
        tmp_m[2*COLORS+2] = conj(tmp_m[0]);
        tmp_m[1] = tmp_m[COLORS] = tmp_m[COLORS+2] = tmp_m[2*COLORS+1] = 0.0;
        tmp_m[COLORS+1] = 1.0;
    }

    //-------------------------------------------------------------------
    // apply this last transformation
    //-------------------------------------------------------------------

    g.DotMEqual(tmp_m, tmp_m1);
}


//-------------------------------------------------------------------------
//                                                        *) sets index[] !
//  recursively calculates discrepancy from gauge condition
//  
//-------------------------------------------------------------------------

IFloat FixHPlane::delta(int recurse)
{
    char *fname = "delta(f)";

    VRB.Func(cname, fname);

    IFloat dlt = 0.0;

    for(index[recurse]=0; index[recurse]<node_size[recurse]; index[recurse]++)
        if(recurse < hplane_dim-1)
            {
                dlt += delta(recurse + 1);
            }
        else
            {
                Matrix& A = findA();   // reference to A_buf
                tmp_m.Dagger(A);
                A -= tmp_m;
                A -= A.Tr()/COLORS;
                for(int i=0; i<COLORS*COLORS; i++)
                    dlt += norm(A[i]);
            }

    return recurse ? dlt : HyperPlane_Sum(dlt);
}


//-------------------------------------------------------------------------
//
//  recursively corrects unitarity of gauge fixing matreces
//
//-------------------------------------------------------------------------

void FixHPlane::unitarize(int recurse)
{
    char *fname = "unitarize(i)";

    VRB.Func(cname, fname);

    if(recurse < hplane_dim) // index has not been completed
        for(index[recurse] = 0; index[recurse] < node_size[recurse];
            index[recurse]++ )
            unitarize(recurse + 1);
    else
        G_loc().Unitarize();
}

//-------------------------------------------------------------------------
/*!
  \pre Lattice::FixGaugeAllocate should be called prior to this. That is how
  the requested type of gauge fixing is communicated to this method.

  \param Smallfloat The stopping criterion, determining the accuracy of the
  gauge fixing.
  \param MaxIterNum The maximum number of iterations the gauge fixing
  algorithm is allowed before it gives up.
  \return The total number of iterations used. If the algorithm fails to
  converge within the allowed number of iterations, then the -1 times the
  number of iterations reached is returned (this might not be MaxIterNum,
  oddly enough).

  \post The required parts of the gauge field are fixed to the desired gauge.
*/
//-------------------------------------------------------------------------//
//                                                                         //
//  Function FixGauge calculates gauge transformation matrices that        //
//  would put the links on the lattice to a specified gauge.               //
//  For the Landau gauge gauge fixing is applied to the whole lattice.     //
//  The condition is that the matrix                                       //
//                       +     1         +                                 //
//        B(n) = A(n) - A(n) - - Tr(A - A )                        (1)     //
//                             3                                           //
//  vanishes for each site n. (See Zhihua Dong's thesis, p. 29.) Here      //
//  (with important changes described before function findA())             //
//                    +                                                    //
//        A(n) = sum[L'j(n) + L'j(n-j)],                           (2)     //
//                j                                                        //
//                                                                         //
//  and the links L'j are related to the original links Lj via formula     //
//                              +                                          //
//        L'j(n) = g(n) L0j(n) g(n+j)                              (3)     //
//                                                                         //
//-------------------------------------------------------------------------//
//                                                                         //
//  Float SmallFloat - stopping condition;                                 //
//                                                                         //
//  int MaxIterNum - issues a warning if reached.                          //
//                                                                         //
//-------------------------------------------------------------------------//
int Lattice::FixGauge(Float SmallFloat, int MaxIterNum)
{
    char *fname = "FixGauge(F,i)";

    VRB.Func(cname, fname);

    //------------------------------------------------------------------------
    // Initial checking
    //------------------------------------------------------------------------
    {
        if( (fix_gauge_kind == FIX_GAUGE_NONE) || (fix_gauge_ptr == NULL) )
            ERR.General(cname, fname, "Not initialized");

        if(Colors() != COLORS)
            ERR.NotImplemented(cname, fname, "Implemented only for Colors = %d. "
                               "Called with Colors == %d\n", COLORS, Colors());
    }

    // set last Stopping Condition
    fix_gauge_stpCond = SmallFloat;


    //------------------------------------------------------------------------
    // NormDir is numerically identical to fix_gauge_kind
    //------------------------------------------------------------------------
    int NormDir = int(fix_gauge_kind);

    //------------------------------------------------------------------------
    // calculate hyperplane dimensionality
    //------------------------------------------------------------------------
    int HplDim = (fix_gauge_kind == FIX_GAUGE_LANDAU) ? 
        XXX::Dimension : XXX::Dimension-1;

    //------------------------------------------------------------------------
    // allocate and initialize Ind2Dir with the initial order skipping 
    // the normal direction 
    //------------------------------------------------------------------------

    int *Ind2Dir;
  
    {
        int mem_size = XXX::Dimension * sizeof(int);
        Ind2Dir = (int*) smalloc(mem_size);
        if(Ind2Dir == NULL)
            ERR.Pointer(cname, fname, "Ind2Dir");
        VRB.Smalloc(cname, fname, "Ind2Dir", Ind2Dir, mem_size);
    
        int i, ii;
        for(ii=i=0; i<HplDim; i++, ii++)
            {
                if(ii == NormDir)
                    ii++;
                Ind2Dir[i] = ii;
            }
    }
  
  
    //------------------------------------------------------------------------
    // find stride - between hplanes in the L[]
    //------------------------------------------------------------------------
  
    int stride = XXX::SiteSize;
  
    if(fix_gauge_kind != FIX_GAUGE_LANDAU)
        {
            for(int i=0; i<NormDir; i++)
                stride *= XXX::Node_Size_in_Dir(i);
        }
  
    //------------------------------------------------------------------------
    // for each hyperplane - iterate until the stopping condition is reached
    //------------------------------------------------------------------------

    Float not_converged = 0;
    Float tot_iternum = 0;

    {
        int loop_num = 
            (fix_gauge_kind==FIX_GAUGE_LANDAU) ? 1 : XXX::Node_Size_in_Dir(NormDir);

        for(int i=0; i < loop_num; i++)
            {
                //------------------------------------------------------------------
                // if initialized - fix
                //------------------------------------------------------------------

                if( fix_gauge_ptr[i] != NULL )
                    {
	
                        //--------------------------------------------------------------
                        // construct hyperplane
                        //--------------------------------------------------------------
	    
                        FixHPlane hplane(GaugeField() + i*stride, Ind2Dir, HplDim, 
                                         fix_gauge_ptr[i], SmallFloat);
	
                        //--------------------------------------------------------------
                        // find gauge fixing matrices for hyperplane
                        //--------------------------------------------------------------

                        int iternum = 0;

                        while((hplane.delta() > hplane.small_enough)
                              && (iternum < MaxIterNum))
                            {
                                //----------------------------------------------------------
                                // Make iterations
                                //----------------------------------------------------------

                                for(int j=0; j<XXX::CheckFreq; j++)
                                    {
                                        hplane.iter();
                                        iternum++;
                                    }

                                //----------------------------------------------------------
                                // Make sure the matrices are still unitary
                                //----------------------------------------------------------

                                hplane.unitarize();
                            }

                        if(iternum >= MaxIterNum)
                            not_converged += 1;
	    
                        tot_iternum += iternum;
                    }
            }
    }

    //------------------------------------------------------------------------
    // deallocate Ind2Dir
    //------------------------------------------------------------------------
    VRB.Sfree(cname, fname, "Ind2Dir", Ind2Dir);
    sfree(Ind2Dir);
  
    // Added by Hantao, sum in t direction.
    glb_sum(&tot_iternum);
    tot_iternum /= GJP.Xnodes() * GJP.Ynodes() * GJP.Znodes();
    //Add by Jianglei
    VRB.Result(cname, fname, "Iteration numbers = %f\n", tot_iternum);

    //--------------------------------------------------------------
    // Issue a warning through broadcast if MaxIterNum is reached
    //--------------------------------------------------------------
    glb_sum(&not_converged);
    if (not_converged > 0.5) 
        {
            VRB.Warn(cname,fname, 
                     "Some hyperplanes did not reach accuracy in %d iterations\n",
                     MaxIterNum);
            return -(int)tot_iternum;
        }
    else
        {
            return (int)tot_iternum;
        }
}



//-------------------------------------------------------------------------//
/*!
  If Landau gauge fixing is requested, memory is allocated for a gauge fixing
  matrix at each lattice site. If Coulomb gauge fixing is requested, memory is
  allocated for a gauge fixing matrix at each lattice site on the hyperplanes
  requested, if that hyperplane intersects the local lattice on this node.
  The matrices are initialised to unity.

  \param GaugeKind The  type of gauge fixing.
  \param NHplanes The total number of the hyperplanes on which to fix the 
  gauge. If set to zero when Coulomb gauge is requested, then treated as
  a request to fix all hyperplanes in the global lattice.
  This is ignored for Landau gauge fixing.
  \param Hplanes A list of the positions of the hyperplanes on which to fix
  the gauge.along the direction orthogonal to them. This list should be in
  increasing order.  This is ignored for Landau gauge fixing and when \a
  NHplanes is zero.
*/
//-------------------------------------------------------------------------//

void Lattice::FixGaugeAllocate(FixGaugeType GaugeKind,
			       int NHplanes, int *Hplanes)
{
    char *fname = "FixGaugeAllocate(FixGaugeType,i,i*)";

    VRB.Func(cname, fname);


    //------------------------------------------------------------------------
    // Initial checking
    //------------------------------------------------------------------------
    {
        if( (fix_gauge_kind != FIX_GAUGE_NONE) || (fix_gauge_ptr != NULL) )
            ERR.General(cname, fname, "Already initialized");
    
        if(GaugeKind == FIX_GAUGE_NONE)
            ERR.General(cname, fname, "Can't fix to FIX_GAUGE_NONE");

        if( (GaugeKind != FIX_GAUGE_LANDAU) && (NHplanes != 0) )
            if(Hplanes == NULL)
                {
                    ERR.General(cname, fname, "Hplanes is NULL");
                }
            else
                {
                    //--- Hplanes should be in increasing order
                    //--- This condition is assumed when deallocating
                    //--- fix_gauge_ptr[] entries in FixGaugeFree

                    int pn = -1; //--- all entries should be nonnegative

                    for(int i=0; i<NHplanes; i++)
                        if(Hplanes[i] > pn)
                            pn = Hplanes[i];
                        else
                            ERR.General(cname, fname, "Hplane not in increasing order");
                }

        if(Colors() != COLORS)
            ERR.NotImplemented(cname, fname, "Implemented only for Colors == %d. "
                               "Called with Colors == %d\n", COLORS, Colors());
    }

    //------------------------------------------------------------------------
    // set system parameter
    //------------------------------------------------------------------------
    fix_gauge_kind = GaugeKind;
    fix_gauge_stpCond = 0.;

    //------------------------------------------------------------------------
    // set some local parameters
    //------------------------------------------------------------------------
    int GaugePtrLen;
    int MemBlockLen;
    int NormDir;

    {
        //--- NormDir is numerically equal to fix_gauge_kind
        NormDir = int(fix_gauge_kind);

        //--- number of fields in fix_gauge_ptr
        GaugePtrLen = (fix_gauge_kind == FIX_GAUGE_LANDAU) ?
            1 : XXX::Node_Size_in_Dir(NormDir); 

        //--- Number of matrices in each memory block
        MemBlockLen = 1;
        for(int i=0; i<XXX::Dimension; i++)
            if(i != NormDir)
                MemBlockLen *= XXX::Node_Size_in_Dir(i);
    }

    //------------------------------------------------------------------------
    // allocate fix_gauge_ptr and initialize fields to NULL
    //------------------------------------------------------------------------
    {
        int mem_size = GaugePtrLen * sizeof(Matrix*);
    
        fix_gauge_ptr = (Matrix**) 
        smalloc(cname, fname, "fix_gauge_ptr", mem_size);
//        printf("fix_gauge_ptr(%d)=%p\n",UniqueID(),fix_gauge_ptr);
    
        for(int i=0; i<GaugePtrLen; i++)
            fix_gauge_ptr[i] = NULL;
    }

    //------------------------------------------------------------------------
    // allocate buffers for hyperplanes intersecting with the node
    //------------------------------------------------------------------------
    {
        //--- memory required for each buffer
        int mem_size = MemBlockLen * sizeof(Matrix);

        if( (fix_gauge_kind == FIX_GAUGE_LANDAU) || (NHplanes == 0) )
            {
                for(int i=0; i<GaugePtrLen; i++)
                    {
                        fix_gauge_ptr[i] = (Matrix*) smalloc(mem_size);
                        if(fix_gauge_ptr[i] == NULL)
                            ERR.Pointer(cname, fname, "fix_gauge_ptr[]");
                        VRB.Smalloc(cname, fname, "fix_gauge_ptr[]", 
                                    fix_gauge_ptr[i], mem_size);
                    }
            }
        else
            {
                //--- look through the whole list of global hyperplanes
                for(int i=0; i<NHplanes; i++)
                    {
                        // since NormDir is not Landau, safe to use here
                        int crd_min= XXX::Coor4d(NormDir) * XXX::Node_Size_in_Dir(NormDir);
                        int crd_max = crd_min + XXX::Node_Size_in_Dir(NormDir) - 1;
	  
                        //--- if hyperplane does intersect with the node
                        if( (Hplanes[i] >= crd_min) && (Hplanes[i] <= crd_max))
                            {
                                int crd_loc = Hplanes[i] - crd_min;

                                //--- allocate space for a local hyperplane
                                fix_gauge_ptr[crd_loc] = (Matrix*) smalloc(mem_size);
                                if(fix_gauge_ptr[crd_loc] == NULL)
                                    ERR.Pointer(cname, fname, "fix_gauge_ptr[i]");
                                VRB.Smalloc(cname, fname, "fix_gauge_ptr[i]", 
                                            fix_gauge_ptr[crd_loc], mem_size);
                            }
                    }
            }
    }

    //------------------------------------------------------------------------
    // initialize matrices to unity
    //------------------------------------------------------------------------
    {
        for(int i=0; i<GaugePtrLen; i++)
            if(fix_gauge_ptr[i] != NULL)
                for(int j=0; j<MemBlockLen; j++)
                    fix_gauge_ptr[i][j].UnitMatrix();
    }
}


//-------------------------------------------------------------------------//
//                                                                         //
// Deallocates memory                                                      //
//                                                                         //
//-------------------------------------------------------------------------//
/*!
  \post Lattice::Lattice::FixGaugeKind will now return FIX_GAUGE_NONE.
*/
  
void Lattice::FixGaugeFree()
{
    char *fname = "FixGaugeFree()";

    VRB.Func(cname, fname);


    //------------------------------------------------------------------------
    // Initial checking
    //------------------------------------------------------------------------
    {
        if( (fix_gauge_kind == FIX_GAUGE_NONE) || (fix_gauge_ptr == NULL) )
            ERR.General(cname, fname, "Not initialized");
    }


    //------------------------------------------------------------------------
    // Deallocate hyperplane buffers
    //------------------------------------------------------------------------
    {
        //--- Normal direction is numerically equal to fix_gauge_kind
        int NormDir = int(fix_gauge_kind);

        int GaugePtrLen = (fix_gauge_kind == FIX_GAUGE_LANDAU) ?
            1 : XXX::Node_Size_in_Dir(NormDir);

        //--- Notice: to avoid memory fragmentation here, 
        //--- entries in Hplanes[] should be in increasing order
        for(int i=GaugePtrLen; i--; )
            if(fix_gauge_ptr[i] != NULL)
                {
                    VRB.Sfree(cname, fname, "fix_gauge_ptr[i]", fix_gauge_ptr[i]);
                    sfree(fix_gauge_ptr[i]);
                }
    }

    //------------------------------------------------------------------------
    // deallocate the fix_gauge_ptr buffer itself
    //------------------------------------------------------------------------
    {
        VRB.Sfree(cname, fname, "fix_gauge_ptr", fix_gauge_ptr);
        sfree(fix_gauge_ptr);
    }

    //------------------------------------------------------------------------
    // reset system variables
    //------------------------------------------------------------------------
    {
        fix_gauge_ptr = NULL;
        fix_gauge_kind = FIX_GAUGE_NONE;
        fix_gauge_stpCond = 0;
    }
}

CPS_END_NAMESPACE
