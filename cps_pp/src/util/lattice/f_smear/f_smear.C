#include <config.h>
/*!\file
  \brief  Implementation of Fsmear class.

  $Id: f_smear.C,v 1.3 2004-08-18 11:58:03 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:03 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_smear/f_smear.C,v 1.3 2004-08-18 11:58:03 zs Exp $
//  $Id: f_smear.C,v 1.3 2004-08-18 11:58:03 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/f_smear/f_smear.C,v $
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

Fsmear::Fsmear(int n_smear){
	cname = "Fsmear";
	char *fname = "Fsmear()";
	VRB.Func(cname,fname);
	n_fields = n_smear;
	fields = (Matrix **)smalloc(n_smear * sizeof(Matrix *));

 	int size = GJP.VolNodeSites()*4;
    for(int i = 0;i<n_smear;i++){
		fields[i] = (Matrix *)smalloc(size * sizeof(Matrix));
    }
	smeared = 0;
}

Fsmear::~Fsmear(){
	char *fname = "~Fsmear()";
	VRB.Func(cname,fname);
    for(int i = 0;i<n_fields;i++){
		sfree(fields[i]);
    }
	sfree(fields);
}

Matrix *Fsmear::Fields(int n){
   if(n >= n_fields)
   ERR.General(cname,"Fields(i)","n(%d) is larger than the number of fields (%d)\n",n,n_fields);
   return fields[n];
}
CPS_END_NAMESPACE

