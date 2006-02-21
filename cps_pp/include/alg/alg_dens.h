#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definitions of the AlgDens class.

  $Id: alg_dens.h,v 1.2 2006-02-21 02:44:11 chulwoo Exp $
*/
//---------------------------------------------------------------------------

#ifndef INCLUDED_ALG_DENS_H
#define INCLUDED_ALG_DENS_H           //!< Prevent multiple inclusion.

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
#include <alg/dens_arg.h>
CPS_START_NAMESPACE

//---------------------------------------------------------------------------
//! Class for quark condensate calculation.
/*! 
  The derivatives of the partition function with respect ot the chemical
  potential are computed stochastically using the Conjugate Gradient algorithm. 
  They can be computed for a number of different fermion masses.

\ingroup alg 
 */
//---------------------------------------------------------------------------
class AlgDens : public Alg
{
 private:
    char *cname;

    // The argument structure for the pbp algorithm
    DensArg *alg_dens_arg;
 
    // Node checkerboard size of the fermion field
    int f_size;
    int map_table[100];
    int save_table[100];
    int load_table[100];
    int refresh_table[100];
    int a_coor_table[100];
    int b_coor_table[100];

    // The source vector
    Vector *src;
    Vector *srcM;

    // The solution vectors
    Vector *save;
    Vector *sol;
    Vector *solM;
    //int *map_table;
    //int *save_table;
    //int *load_table;
    //int *refresh_table;
    //int *a_coor_table;
    //int *b_coor_table;

    void make_tables(void);
    void clear_vector(void);
    void save_vector(Vector *v, int obsID);
    void load_vector(Vector *v, int obsID);

 public:
    AlgDens(Lattice & latt, CommonArg *c_arg, DensArg *arg);

    virtual ~AlgDens();

    void run(void);

};

#endif





CPS_END_NAMESPACE
