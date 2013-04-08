#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of methods converting between different data layouts.

*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/convert/convert.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/lattice.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <util/gjp.h>
//#include <comms/nga_reg.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE

//! A data container in the layout conversion routines.
typedef struct ConvertArgStruct { 
	unsigned	lx ;
	unsigned	ly ;
	unsigned	lz ;
	unsigned	lt ;
	unsigned	nc ;
	unsigned	ns ;
	unsigned	vol ;
	unsigned	site_size ;
	Float		*start_ptr ;
} CAS, *CAP ;

const unsigned CBUF_MODE4 = 0xcca52112;
#ifdef _TARTAN
static Matrix *mp2 = (Matrix *)CRAM_SCRATCH_ADDR+2;
#else
static Matrix tmp_2;
static Matrix *mp2 = &tmp_2;
#endif

extern void MultStagPhases(CAP cap) ;
extern void RunGConverter(CAP cap, unsigned *site_tbl, unsigned *link_tbl) ;
extern void CanonToAnything(CAP cap, StrOrdType new_str_ord) ;

extern void FcanonToWilson(CAP cap, int number_of_checkerboards = 2) ;
extern void FwilsonToCanon(CAP cap, int number_of_checkerboards = 2) ;
extern void FcanonToStag(CAP cap, int number_of_checkerboards = 2) ;
extern void FstagToCanon(CAP cap, int number_of_checkerboards = 2) ;


//-------------------------------------------------------------------------
// Convert(StrOrdType new_str_ord, 
//         Vector *f_field_1,
//         Vector *f_field_2) 
// If str_ord is not the same as 
// new_str_ord then it converts the gauge field
// configuration and the two fermion fields f_field_1, 
// f_field_2 to new_str_ord.
/*!
  The fields exist on all lattice sites (both parities).
  \param new_str_ord The new order in which the fields will be laid out.
  \param f_field_1 A fermion  field.
  \param f_field_2 A fermion  field.
  \post The fermion fields and the gauge field are laid out in the new order.
*/
//-------------------------------------------------------------------------
void Lattice::Convert(StrOrdType new_str_ord, 
		     Vector *f_field_1,
		     Vector *f_field_2) 
{
	char *fname = "Convert(StrOrdType,V*,V*)";
	VRB.Func(cname,fname);

//-------------------------------------------------------------------------
// Check if conversion is needed
//-------------------------------------------------------------------------
	if (new_str_ord == str_ord) {
		VRB.Flow(cname,fname,
			"No conversion necessary from %d to %d\n",
			int(str_ord), int(new_str_ord));
		return ;
	}
	
//-------------------------------------------------------------------------
// Convert the fermion fields 1 and 2.
//-------------------------------------------------------------------------
	if( f_field_1 == 0 )
		ERR.Pointer(cname,fname, "f_field_1");
	if( f_field_2 == 0 )
		ERR.Pointer(cname,fname, "f_field_2");

	Fconvert(f_field_1, new_str_ord, str_ord);
	
	Fconvert(f_field_2, new_str_ord, str_ord);

//-------------------------------------------------------------------------
// Convert the gauge field
//-------------------------------------------------------------------------
	Convert(new_str_ord) ;

	return ;
}


void Lattice::Convert(StrOrdType new_str_ord, 
		     Vector *f_field_1) 
{
	char *fname = "Convert(StrOrdType,V*)";
	VRB.Func(cname,fname);

//-------------------------------------------------------------------------
// Check if conversion is needed
//-------------------------------------------------------------------------
	if (new_str_ord == str_ord) {
		VRB.Flow(cname,fname,
			"No conversion necessary from %d to %d\n",
			int(str_ord), int(new_str_ord));
		return ;
	}
	
//-------------------------------------------------------------------------
// Convert the fermion fields 1
//-------------------------------------------------------------------------
	if( f_field_1 == 0 )
		ERR.Pointer(cname,fname, "f_field_1");
	Fconvert(f_field_1, new_str_ord, str_ord);
	
	return ;
}




//------------------------------------------------------------------
// Convert(StrOrdType new_str_ord):
/*!
  The field exists on all lattice sites (both parities).
  \param new_str_ord The new order in which the field will be laid out.
  \post The gauge field is arranged in the new order.
*/
void Lattice::Convert(StrOrdType new_str_ord)
{
	char *fname = "Convert(StrOrdType)";
	VRB.Func(cname,fname);

//-------------------------------------------------------------------------
// Check if conversion is needed
//-------------------------------------------------------------------------
	if (new_str_ord == str_ord) {
	  VRB.Flow(cname,fname,
		   "No conversion necessary from %d to %d\n",
		   int(str_ord), int(new_str_ord));
	  return ;
	}


	unsigned *site_sort_tbl ;
	unsigned *link_sort_tbl ;
	unsigned x,y,z,t,cb ;
	unsigned r,row,col,mu ;
	unsigned idx = 0 ;

	CAS	cas ;

	cas.lx		= GJP.XnodeSites() ;
	cas.ly		= GJP.YnodeSites() ;
	cas.lz		= GJP.ZnodeSites() ;
	cas.lt		= GJP.TnodeSites() ;
	cas.nc		= Colors() ;
	cas.site_size	= GsiteSize() ;
	cas.start_ptr	= (Float *) GaugeField() ;
	cas.vol		= GJP.VolNodeSites() ;

//-------------------------------------------------------------------------
// first convert to CANONICAL
//-------------------------------------------------------------------------

	switch(str_ord) {
		case STAG :
	                VRB.Flow(cname,fname,
			"Converting gauge field order: STAG -> CANONICAL\n");

			MultStagPhases(&cas) ; // STAG -> CANONICAL

			//---------------------------------------------------
			//  de-Dagger all links
			//---------------------------------------------------
			{
			  Matrix *p = GaugeField();
			  int n_links = 4 * GJP.VolNodeSites();
			  
			  setCbufCntrlReg(4, CBUF_MODE4);
			  
			  for(int i = 0; i < n_links; ++i) {
			    mp2->Dagger((IFloat *)p);
//			    moveMem((IFloat *)(p++), (IFloat *)mp2,
//				    18*sizeof(IFloat));
			    moveFloat((IFloat *)(p++), (IFloat *)mp2, 18);
			  }
			}

			break ;
		case STAG_BLOCK :
		  {
	                VRB.Flow(cname,fname,
			"Converting gauge field order: STAG_BLOCK -> CANONICAL\n");
			//Copy the links into tmp_gauge into STAG ordering
			//Then, copy tmp_gauge back to GaugeField()
			Matrix * tmp_gauge = (Matrix *)
			fmalloc(cname,fname,"tmp_gauge",cas.vol*4*sizeof(Matrix));
			int x[4];
			int mu,current,new_index;
			for(mu = 0; mu < 4; mu++)
			  for(x[2]=0;x[2]<cas.lz;x[2]++)
			    for(x[1]=0;x[1]<cas.ly;x[1]++)
			      for(x[0]=0;x[0]<cas.lx;x[0]++)
				for(x[3]=0;x[3]<cas.lt;x[3]++)
				  {
				    current = cas.vol*mu+(x[3]+cas.lt*(x[0]+cas.lx*(x[1]+cas.ly*x[2])))/2+((x[3]+x[2]+x[1]+x[0])%2)*cas.vol/2;
				    new_index = 4*(x[0]+cas.lx*(x[1]+cas.ly*(x[2]+cas.lz*x[3])))+mu;
//				    moveMem(tmp_gauge+new_index,cas.start_ptr+18*current,sizeof(Matrix));
				    moveFloat((Float *)(tmp_gauge+new_index),(Float*)(cas.start_ptr+18*current),18);
				  }
//			moveMem(cas.start_ptr,tmp_gauge,cas.vol*4*sizeof(Matrix));
			moveFloat((Float*)cas.start_ptr,(Float*)tmp_gauge,cas.vol*4*18);
			ffree(cname,fname, "tmp_gauge", tmp_gauge);

			MultStagPhases(&cas) ; // STAG -> CANONICAL

			//---------------------------------------------------
			//  de-Dagger all links
			//---------------------------------------------------
			{
			  Matrix *p = GaugeField();
			  int n_links = 4 * GJP.VolNodeSites();
			  
			  setCbufCntrlReg(4, CBUF_MODE4);
			  
			  for(int i = 0; i < n_links; ++i) {
			    mp2->Dagger((IFloat *)p);
//			    moveMem((IFloat *)(p++), (IFloat *)mp2,
//				    18*sizeof(IFloat));
			    moveFloat((IFloat *)(p++), (IFloat *)mp2, 18);
			  }
			}
		  }

			break ;
          	case DWF_4D_EOPREC: 
		case DWF_4D_EOPREC_EE: 
		case WILSON :
			VRB.Flow(cname,fname,
			"Converting gauge field order: WILSON -> CANONICAL\n");

			site_sort_tbl = (unsigned *)
			fmalloc(cname,fname,
				    "site_sort_tbl" , 
				    cas.vol*sizeof(unsigned) );
			

			idx = 0 ;

//-------------------------------------------------------------------------
// Loop over current (WILSON) order of sites in lattice
//-------------------------------------------------------------------------

			for (cb=0; cb<2; cb++)
			for (t=0; t<cas.lt; t++)
			for (z=0; z<cas.lz; z++)
			for (y=0; y<cas.ly; y++)
			for (x=0; x<cas.lx; x++) {
				if ((x+y+z+t)%2 == cb) {

//-------------------------------------------------------------------------
// Use desired (CANONICAL) equation for site sequence number
// LSB = 1 indicates site needs converting
//-------------------------------------------------------------------------

					*(site_sort_tbl+idx) = x +
					cas.lx*(y+cas.ly*(z+cas.lz*t))<<1 | 1 ;

					idx++;
				}
			}

//-------------------------------------------------------------------------
// link_sort_tbl should be in CRAM
//-------------------------------------------------------------------------

			
			link_sort_tbl = (unsigned *)
			fmalloc(cname,fname,
				    "link_sort_tbl" , 
				    cas.site_size*sizeof(unsigned)) ;

//-------------------------------------------------------------------------
// Loop over current (WILSON) order of reals in site
//-------------------------------------------------------------------------

			idx = 0 ;

			for (mu=0; mu<4; mu++)		// 4 is spacetime dim
			for (col=0; col<cas.nc; col++)
			for (row=0; row<cas.nc; row++)
			for (r=0; r<2; r++) {		// 2 is # cplx comp

//-------------------------------------------------------------------------
// Use desired (CANONICAL) equation for reals in site
//-------------------------------------------------------------------------

				*(link_sort_tbl+idx) =
					r+2*(col+cas.nc*(row+cas.nc*mu)) ;
				idx++ ;
			}

			RunGConverter(&cas, site_sort_tbl, link_sort_tbl) ;

			ffree(cname,fname, "link_sort_tbl", link_sort_tbl);
			ffree(cname,fname, "site_sort_tbl", site_sort_tbl);

			break ;

		case G_WILSON_HB :

			VRB.Flow(cname,fname,
				 "Converting gauge field order: %s",
				 "G_WILSON_HB -> CANONICAL\n");

			site_sort_tbl = (unsigned *)
			fmalloc(cname,fname,
				    "site_sort_tbl" , 
				    cas.vol*sizeof(unsigned) );



			idx = 0 ;

//-------------------------------------------------------------------------
// Loop over current (G_WILSON_HB) order of sites in lattice
//-------------------------------------------------------------------------

			for (t=0; t<cas.lt; t++)
			for (x=0; x<cas.lx; x++)
			for (y=0; y<cas.ly; y++)
			for (z=0; z<cas.lz; z++) {

//-------------------------------------------------------------------------
// Use desired (CANONICAL) equation for site sequence number
// LSB = 1 indicates site needs converting
//-------------------------------------------------------------------------

				*(site_sort_tbl+idx) = x +
					cas.lx*(y+cas.ly*(z+cas.lz*t))<<1 | 1 ;

				idx++;
			}

//-------------------------------------------------------------------------
// link_sort_tbl should be in CRAM
//-------------------------------------------------------------------------

			link_sort_tbl = (unsigned *)
			fmalloc(cname,fname,
				    "link_sort_tbl" , 
				    cas.site_size*sizeof(unsigned));

//-------------------------------------------------------------------------
// Loop over current (G_WILSON_HB) order of reals in site
//-------------------------------------------------------------------------

			idx = 0 ;

			for (mu=3; mu<7; mu++)		// 4 is spacetime dim
			for (r=0; r<2; r++) 		// 2 is # cplx comp
			for (row=0; row<cas.nc; row++)
			for (col=0; col<cas.nc; col++) {

//-------------------------------------------------------------------------
// Use desired (CANONICAL) equation for reals in site
//-------------------------------------------------------------------------

				*(link_sort_tbl+idx) = r + 2 *
					(col+cas.nc*(row+cas.nc*(mu%4))) ;
				idx++ ;
			}

			RunGConverter(&cas, site_sort_tbl, link_sort_tbl) ;


			ffree(cname,fname, "link_sort_tbl", link_sort_tbl);
			ffree(cname,fname, "site_sort_tbl", site_sort_tbl);

			break ;
		case CANONICAL :
		default :
			break ;
	}

	str_ord = CANONICAL ;

//-------------------------------------------------------------------------
// Next, convert from CANONICAL to new_str_ord
//-------------------------------------------------------------------------

	CanonToAnything(&cas,new_str_ord) ;

	str_ord = new_str_ord ;

//  End  GRF
}

static char *fname_fconvert = "Fconvert(V*,StrOrdType,StrOrdType, int)";
static char *converting_str = "Converting frm field str ord from %d to %d\n";


//------------------------------------------------------------------
//! Does nothing!
//------------------------------------------------------------------
void Fnone::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from, int cb) 
{
	VRB.Func(cname,fname_fconvert);
}

//-------------------------------------------------------------------------
//! Converts the field layout between the canonical and the odd-even order..
/*!
  \param f_field The fermionic field to be converted.
  \param to The new order; either CANONICAL or WILSON.
  \param from The current order; either CANONICAL or WILSON.
*/
//-------------------------------------------------------------------------

void FwilsonTypes::Fconvert(Vector *f_field,StrOrdType to,StrOrdType from, int cb)
{
	CAS	cas ;

	cas.lx		= GJP.XnodeSites() ;
	cas.ly		= GJP.YnodeSites() ;
	cas.lz		= GJP.ZnodeSites() ;
	cas.lt		= GJP.TnodeSites() ;
	cas.nc		= Colors() ;
	cas.ns		= SpinComponents() ;
	cas.site_size	= FsiteSize() ;
	cas.start_ptr	= (Float *) f_field ;
	cas.vol		= GJP.VolNodeSites() ;

	VRB.Func(cname,fname_fconvert);

	if (from == to) {
		VRB.Flow(cname,fname_fconvert,
			"No conversion necessary from %d to %d\n",
			int(from), int(to));
		return ;
	}

	if ((from == CANONICAL) && (to == WILSON)) {

	  VRB.Flow(cname,fname_fconvert, converting_str,
		   int(from), int(to));

	  FcanonToWilson(&cas) ;

	} else if ((from == WILSON) && (to == CANONICAL)) {
	  
	  VRB.Flow(cname,fname_fconvert, converting_str,
		   int(from), int(to));

	  FwilsonToCanon(&cas) ;

	} else {
		ERR.General(cname,fname_fconvert,
			"Unsupported fermion conversion from %d to %d\n",
			int(from), int(to)) ;
	}

	return ;
}


//------------------------------------------------------------------
//! Converts the field layout between the canonical and the odd-even order..
/*!
  \param f_field The fermionic field to be converted.
  \param to The new order; either CANONICAL or STAG.
  \param from The current order; either CANONICAL or STAG.
*/
//------------------------------------------------------------------
void FstagTypes::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from, int num_checkerboard) 
{
  VRB.Func(cname,fname_fconvert);
  
  CAS     cas ;
  
  
  cas.lx          = GJP.XnodeSites() ;
  cas.ly          = GJP.YnodeSites() ;
  cas.lz          = GJP.ZnodeSites() ;
  cas.lt          = GJP.TnodeSites() ;
  cas.nc          = Colors() ;
  cas.ns          = SpinComponents() ;
  cas.site_size   = FsiteSize() ;
  cas.start_ptr   = (Float *) f_field ;
  cas.vol         = GJP.VolNodeSites() ;
  
  if (from == to) {
    VRB.Flow(cname,fname_fconvert,
	     "No conversion necessary from %d to %d\n",
	     int(from), int(to));
    return ;
  }
  
  if ((from == CANONICAL) && (to == STAG)) {
    
    VRB.Flow(cname,fname_fconvert, converting_str,
	     int(from), int(to));
    FcanonToStag(&cas,num_checkerboard) ;
    
  } else if ((from == STAG) && (to == CANONICAL)) {
    
    VRB.Flow(cname,fname_fconvert, converting_str,
	     int(from), int(to));    
    FstagToCanon(&cas,num_checkerboard) ;
    
  } else {
    ERR.General(cname,fname_fconvert,
		"Unsupported fermion conversion from %d to %d\n",
		int(from), int(to)) ;
  }
  
  return ;
}










//------------------------------------------------------------------
//! Converts the field layout between the canonical and the odd-even order..
/*!
  \param f_field The fermionic field to be converted.
  \param to The new order; either CANONICAL or WILSON.
  \param from The current order; either CANONICAL or WILSON.
*/
//------------------------------------------------------------------
void FdwfBase::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from, int cb) 
{
	Float *field_ptr;
	Float *tmp_field_ptr;
	int i;
	int parity;
	CAS	cas ;

  	int ls = GJP.SnodeSites();

	cas.lx		= GJP.XnodeSites() ;
	cas.ly		= GJP.YnodeSites() ;
	cas.lz		= GJP.ZnodeSites() ;
	cas.lt		= GJP.TnodeSites() ;
	cas.nc		= Colors() ;
	cas.ns		= SpinComponents() ;
	cas.site_size	= FsiteSize() / ls ;
	cas.start_ptr	= (Float *) f_field ;
	cas.vol		= GJP.VolNodeSites() ;

        int f_size = cas.vol * FsiteSize();
        int stride = f_size / ls; 
        int half_stride = stride / 2;

	VRB.Func(cname,fname_fconvert);

	if (from == to) {
		VRB.Flow(cname,fname_fconvert,
			"No conversion necessary from %d to %d\n",
			int(from), int(to));
		return ;
	}


	// Allocate memory for a temporary
	//----------------------------------------------------------
	Vector *tmp_f_field = (Vector *) 
	fmalloc(cname,fname_fconvert, 
		    "tmp_f_field", f_size * sizeof(Float));


	if ((from == CANONICAL) && (to == WILSON)) {

	  VRB.Flow(cname,fname_fconvert, converting_str,
		   int(from), int(to));
	  
	  // convert from canonical to intermediate conversion
	  field_ptr = (Float *) f_field;
	  for(i=0; i<ls; i++){
	    cas.start_ptr = field_ptr;
	    FcanonToWilson(&cas) ;
	    field_ptr = field_ptr + stride;
	  }

	  // copy intermediate converted vector to a buffer
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  moveFloat(tmp_field_ptr, field_ptr, f_size );

	  // Set odd part
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  for(i=0; i<ls; i++){
	    parity = (i+1) % 2;
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + half_stride;
	    tmp_field_ptr = tmp_field_ptr + (2 * parity + 1) * half_stride; 
	  }

	  // Set even part
	  tmp_field_ptr = (Float *) tmp_f_field;
	  tmp_field_ptr = tmp_field_ptr + half_stride;
	  for(i=0; i<ls; i++){
	    parity = i % 2;
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + half_stride;
	    tmp_field_ptr = tmp_field_ptr + (2 * parity + 1) * half_stride; 
	  }
	  
	} else if ((from == CANONICAL) && (to == DWF_4D_EOPREC)) {

	  VRB.Flow(cname,fname_fconvert, converting_str,
		   int(from), int(to));
	  
	  // convert from canonical to odd-even on each s-slice
	  field_ptr = (Float *) f_field;
	  for(i=0; i<ls; i++){
	    cas.start_ptr = field_ptr;
	    FcanonToWilson(&cas) ;
	    field_ptr = field_ptr + stride;
	  }

	  // copy intermediate converted vector to a buffer
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  moveFloat(tmp_field_ptr, field_ptr, f_size );

	  // Set odd part: 4d odd-even.
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  for(i=0; i<ls; i++){
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + half_stride;
	    tmp_field_ptr = tmp_field_ptr + stride; 
	  }

	  // Set even part
	  tmp_field_ptr = (Float *) tmp_f_field;
	  tmp_field_ptr = tmp_field_ptr + half_stride;
	  for(i=0; i<ls; i++){
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + half_stride;
	    tmp_field_ptr = tmp_field_ptr + stride; 
	  }


	} else if ((from == CANONICAL) && (to == DWF_4D_EOPREC_EE )) {

	  VRB.Flow(cname,fname_fconvert, converting_str,
		   int(from), int(to));
	  
	  // convert from canonical to odd-even on each s-slice
	  field_ptr = (Float *) f_field;
	  for(i=0; i<ls; i++){
	    cas.start_ptr = field_ptr;
	    FcanonToWilson(&cas) ;
	    field_ptr = field_ptr + stride;
	  }

	  // copy intermediate converted vector to a buffer
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  moveFloat(tmp_field_ptr, field_ptr, f_size );

	  // Set even part:
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  tmp_field_ptr = tmp_field_ptr + half_stride;
	  for(i=0; i<ls; i++){
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + half_stride;
	    tmp_field_ptr = tmp_field_ptr + stride; 
	  }

	  // Set odd part
	  tmp_field_ptr = (Float *) tmp_f_field;
	  for(i=0; i<ls; i++){
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + half_stride;
	    tmp_field_ptr = tmp_field_ptr + stride; 
	  }


	} else if ((from == WILSON) && (to == CANONICAL)) {
	  
	  VRB.Flow(cname,fname_fconvert, converting_str,
		   int(from), int(to));


	  // copy vector to a buffer for intermediate conversion
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
//	  moveMem(tmp_field_ptr, field_ptr, f_size * sizeof(Float));
	  moveFloat(tmp_field_ptr, field_ptr, f_size);


	  // convert odd part to intermediate conversion
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  for(i=0; i<ls; i++){
	    parity = (i+1) % 2;
//	    moveMem(field_ptr, tmp_field_ptr, half_stride * sizeof(Float));
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + (2 * parity + 1) * half_stride; 
	    tmp_field_ptr = tmp_field_ptr + half_stride;
	  }

	  // convert even part to intermediate conversion
	  field_ptr = (Float *) f_field;
	  field_ptr = field_ptr + half_stride;
	  for(i=0; i<ls; i++){
	    parity = i % 2;
//	    moveMem(field_ptr, tmp_field_ptr, half_stride * sizeof(Float));
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + (2 * parity + 1) * half_stride; 
	    tmp_field_ptr = tmp_field_ptr + half_stride;
	  }

	  // convert from intermediate conversion to canonical
	  field_ptr = (Float *) f_field;
	  for(i=0; i<ls; i++){
	    cas.start_ptr = field_ptr;
	    FwilsonToCanon(&cas) ;
	    field_ptr = field_ptr + stride;
	  }

	} else if ((from == DWF_4D_EOPREC) && (to == CANONICAL)) {
	  
	  VRB.Flow(cname,fname_fconvert, converting_str,
		   int(from), int(to));

	  // copy vector to a buffer for intermediate conversion
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  moveFloat(tmp_field_ptr, field_ptr, f_size);

	  // convert odd part to intermediate conversion
	  for(i=0; i<ls; i++){
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + stride; 
	    tmp_field_ptr = tmp_field_ptr + half_stride;
	  }

	  // convert even part to intermediate conversion
	  field_ptr = (Float *) f_field;
	  field_ptr = field_ptr + half_stride;
	  for(i=0; i<ls; i++){
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + stride; 
	    tmp_field_ptr = tmp_field_ptr + half_stride;
	  }

	  // convert from intermediate conversion to canonical
	  field_ptr = (Float *) f_field;
	  for(i=0; i<ls; i++){
	    cas.start_ptr = field_ptr;
	    FwilsonToCanon(&cas) ;
	    field_ptr = field_ptr + stride;
	  }


	} else if ((from == DWF_4D_EOPREC_EE) && (to == CANONICAL)) {
	  
	  VRB.Flow(cname,fname_fconvert, converting_str,
		   int(from), int(to));

	  // copy vector to a buffer for intermediate conversion
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  moveFloat(tmp_field_ptr, field_ptr, f_size);

	  // convert even part to intermediate conversion
	  field_ptr = field_ptr + half_stride;
	  for(i=0; i<ls; i++){
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + stride; 
	    tmp_field_ptr = tmp_field_ptr + half_stride;
	  }

	  // convert odd part to intermediate conversion
	  field_ptr = (Float *) f_field;
	  for(i=0; i<ls; i++){
	    moveFloat(field_ptr, tmp_field_ptr, half_stride );
	    field_ptr = field_ptr + stride; 
	    tmp_field_ptr = tmp_field_ptr + half_stride;
	  }

	  // convert from intermediate conversion to canonical
	  field_ptr = (Float *) f_field;
	  for(i=0; i<ls; i++){
	    cas.start_ptr = field_ptr;
	    FwilsonToCanon(&cas) ;
	    field_ptr = field_ptr + stride;
	  }


	} else {
		ERR.General(cname,fname_fconvert,
			"Unsupported fermion conversion from %d to %d\n",
			int(from), int(to)) ;
	}


	// Free temporary fermion field memory
	//----------------------------------------------------------
	ffree(cname,fname_fconvert, "tmp_f_field", tmp_f_field);

	return ;

}

CPS_END_NAMESPACE
