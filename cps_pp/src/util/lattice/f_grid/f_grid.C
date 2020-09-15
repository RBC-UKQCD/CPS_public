#include<config.h>
#include<util/lattice/fgrid.h>

#ifdef USE_GRID

CPS_START_NAMESPACE

bool FgridBase::grid_initted=false;


//------------------------------------------------------------------
/*!
  \param five The 5-dimensional field.
  \param four The 4-dimensional field.
  \param s_u The global 5th direction (s) coordinate where the
  upper two components (right chirality) of the 5-dim. field
  take the values of those of the 4-dim. field.
  \param s_l The global 5th direction (s) coordinate where the
  lower two components (left chirality) of the 5-dim. field
  take the values of those of the 4-dim. field.
  \post The 5-dim field is zero everywhere except where the global
  5th direction coordinate (s) is \a s_l or \a s_u, where it takes the values
  explained above.
*/
//------------------------------------------------------------------
void FgridBase::Ffour2five(Vector *five, Vector *four, int s_u, int s_l, int Ncb)
{
  int x;
  int i;
  Float *field_4D;
  Float *field_5D;
  char *fname = "Ffour2five(V*,V*,i,i)";
  VRB.Func(cname,fname);


  //------------------------------------------------------------------
  // Initializations
  //------------------------------------------------------------------
  size_t f_size = GJP.VolNodeSites() * FsiteSize()*Ncb/2;
  if(GJP.Gparity()) f_size*=2;

  int ls = GJP.SnodeSites();
  int vol_4d = GJP.VolNodeSites()*Ncb/2;
  int ls_stride = 24 * vol_4d;
  if(GJP.Gparity()) ls_stride*=2;

  int s_u_local = s_u % GJP.SnodeSites();
  int s_l_local = s_l % GJP.SnodeSites();
  int s_u_node = s_u / GJP.SnodeSites();
  int s_l_node = s_l / GJP.SnodeSites();


  //------------------------------------------------------------------
  // Set *five using the 4D field *four. 
  //------------------------------------------------------------------

  // Set all components of the 5D field to zero.
  //---------------------------------------------------------------
  field_5D  = (Float *) five;
  for(i=0; i<f_size; i++){
    field_5D[i]  = 0.0;
  }

  // Do the two upper spin components if s_u is in the node
  //---------------------------------------------------------------
  if( s_u_node == GJP.SnodeCoor() ){
    field_4D  = (Float *) four;
    field_5D  = (Float *) five;
    field_5D  = field_5D  + s_u_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_5D[i]  = field_4D[i];
      }
      field_4D  = field_4D  + 24;
      field_5D  = field_5D  + 24;
    }
    if(GJP.Gparity()){ //CK:08/11 do second stacked field
      for(x=0; x<vol_4d; x++){
	for(i=0; i<12; i++){
	  field_5D[i]  = field_4D[i];
	}
	field_4D  = field_4D  + 24;
	field_5D  = field_5D  + 24;
      }
    }
  }

  // Do the two lower spin components if s_l is in the node
  //----------------------------------------------------------------
  if( s_l_node == GJP.SnodeCoor() ){
    field_4D  = (Float *) four;
    field_5D  = (Float *) five;
    field_4D  = field_4D  + 12;
    field_5D  = field_5D  + 12 + s_l_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_5D[i]  = field_4D[i];
      }
      field_4D  = field_4D  + 24;
      field_5D  = field_5D  + 24;
    }
    if(GJP.Gparity()){
      for(x=0; x<vol_4d; x++){
	for(i=0; i<12; i++){
	  field_5D[i]  = field_4D[i];
	}
	field_4D  = field_4D  + 24;
	field_5D  = field_5D  + 24;
      }
    }
  }

}


//------------------------------------------------------------------
/*!
  \param four The 4-dimensional field.
  \param five The 5-dimensional field.
  \param s_u The global 5th direction (s) coordinate where 
  the values of the upper two components (right chirality) of the 5-dim. field
  are taken by those of the 4-dim. field.
  \param s_l The global 5th direction coordinate (s) where the values of 
  the lower two components (left chirality) of the 5-dim. field
  are taken by  those of the 4-dim. field.
  \post The 5-dim field is zero everywhere except where the global
  5th direction coordinate is \a s_l or \a s_u, where it takes the values
  explained above.
  \post An identical 4-dim. field is reproduced on all nodes in the s
  direction.
*/
//------------------------------------------------------------------
void FgridBase::Ffive2four(Vector *four, Vector *five, int s_u, int s_l, int Ncb)
{
  int x;
  int i;
  Float *field_4D;
  Float *field_5D;
  char *fname = "Ffive2four(V*,V*,i,i)";
  VRB.Func(cname,fname);


  //------------------------------------------------------------------
  // Initializations
  //------------------------------------------------------------------
  int ls = GJP.SnodeSites();
  size_t f_size = GJP.VolNodeSites() * FsiteSize()*Ncb / (ls*2);
  if(GJP.Gparity()) f_size*=2;

  int vol_4d = GJP.VolNodeSites()*Ncb/2;
  int ls_stride = 24 * vol_4d;
  if(GJP.Gparity()) ls_stride*=2;

  int s_u_local = s_u % GJP.SnodeSites();
  int s_l_local = s_l % GJP.SnodeSites();
  int s_u_node = s_u / GJP.SnodeSites();
  int s_l_node = s_l / GJP.SnodeSites();


  //------------------------------------------------------------------
  // Set *four using the 5D field *five. 
  //------------------------------------------------------------------

  // Set all components of the 4D field to zero.
  //---------------------------------------------------------------
  field_4D  = (Float *) four;
  for(i=0; i<f_size; i++){
    field_4D[i]  = 0.0;
  }

  // Do the two upper spin components if s_u is in the node
  //---------------------------------------------------------------
  if( s_u_node == GJP.SnodeCoor() ){
    field_4D = (Float *) four;
    field_5D = (Float *) five;
    field_5D = field_5D + s_u_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_4D[i] = field_5D[i];
      }
      field_4D = field_4D + 24;
      field_5D = field_5D + 24;
    }
    if(GJP.Gparity()){
      for(x=0; x<vol_4d; x++){
	for(i=0; i<12; i++){
	  field_4D[i] = field_5D[i];
	}
	field_4D = field_4D + 24;
	field_5D = field_5D + 24;
      }
    }

  }
  // Do the two lower spin components if s_l is in the node
  //----------------------------------------------------------------
  if( s_l_node == GJP.SnodeCoor() ){
    field_4D = (Float *) four;
    field_5D = (Float *) five;
    field_4D = field_4D + 12;
    field_5D = field_5D + 12 + s_l_local * ls_stride;
    for(x=0; x<vol_4d; x++){
      for(i=0; i<12; i++){
	field_4D[i] = field_5D[i];
      }
      field_4D = field_4D + 24;
      field_5D = field_5D + 24;
    }
    if(GJP.Gparity()){
      for(x=0; x<vol_4d; x++){
	for(i=0; i<12; i++){
	  field_4D[i] = field_5D[i];
	}
	field_4D = field_4D + 24;
	field_5D = field_5D + 24;
      }
    }
  }

  // Sum along s direction to get the same 4D field in all 
  // s node slices.
  //----------------------------------------------------------------
  if( GJP.Snodes() > 1) {
    Float sum;
    field_4D  = (Float *) four;
    for(i=0; i<f_size; i++){
      sum = field_4D[i];
      glb_sum_dir(&sum, 4);
      field_4D[i] = sum;    
    }
  }

}

CPS_END_NAMESPACE
#endif
