#include <config.h>
/*!\file
  \brief  Implementation of Fsmear class.

  $Id: f_smear.C,v 1.5 2006/02/01 16:46:08 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006/02/01 16:46:08 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/f_smear/f_smear.C,v 1.5 2006/02/01 16:46:08 chulwoo Exp $
//  $Id: f_smear.C,v 1.5 2006/02/01 16:46:08 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/f_smear/f_smear.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// f_smear.C
//
// Fsmear is derived from Lattice and is relevant to the 
// fermion actions with smeared/improved links
//
//------------------------------------------------------------------

#include <util/lattice.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/gjp.h>
CPS_START_NAMESPACE

/*!
  \param n_smear The number of smeared fields.
*/
Fsmear::Fsmear(int n_smear){
    
    cname = "Fsmear";
    const char *fname = "Fsmear";
    VRB.Func(cname,fname);
    
    n_fields = n_smear;
    fields = (Matrix **)smalloc(n_smear * sizeof(Matrix *));
    int size = GJP.VolNodeSites()*4;
    for(int i = 0;i<n_smear;i++){
	fields[i] = (Matrix *)smalloc(size * sizeof(Matrix));
    }
    smeared = 0;
}

Fsmear::~Fsmear()
{
    const char *fname = "~Fsmear()";
    VRB.Func(cname,fname);
    for(int i = 0;i<n_fields;i++) sfree(fields[i]);
    sfree(fields);
}

/*!
  \param n The number of the smeared field to get.
  \return A pointer to the <em>n</em>th smeared field.
*/
Matrix *Fsmear::Fields(int n){

    const char *fname = "Fields";
    if(n >= n_fields)
	ERR.General(cname, fname,"n(%d) is larger than the number of fields (%d)\n",n,n_fields);
    
    return fields[n];
    
}
CPS_END_NAMESPACE

