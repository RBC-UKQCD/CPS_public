#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/convert/convert.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: convert.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.6  2003/02/10 11:23:22  mcneile
//  I have added code for teh asqtad action.
//
//  Revision 1.5  2002/12/04 17:16:27  zs
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
//  Revision 1.4  2001/08/16 10:50:32  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:13:21  anj
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
//  $RCSfile: convert.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/convert/convert.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// convert.C
//
// The interface routine to the conversion programs.  Maintained by
// GRF.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/lattice.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/gjp.h>
#include<comms/nga_reg.h>
#include<comms/cbuf.h>
CPS_START_NAMESPACE

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




//------------------------------------------------------------------
// Convert(StrOrdType new_str_ord):
// If str_ord is not the same as 
// new_str_ord then it converts the gauge field
// configuration to new_str_ord.
//------------------------------------------------------------------
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
			    mp2->Dagger((IFloat *)p+BANK4_BASE);
			    moveMem((IFloat *)(p++), (IFloat *)mp2,
				    18*sizeof(IFloat));
			  }
			}

			break ;
		case WILSON :
			VRB.Flow(cname,fname,
			"Converting gauge field order: WILSON -> CANONICAL\n");

			site_sort_tbl = (unsigned *)
				smalloc(cas.vol*sizeof(unsigned)) ;
			if(site_sort_tbl == 0)
			  ERR.Pointer(cname,fname, "site_sort_tbl"); 
			VRB.Smalloc(cname,fname,
				    "site_sort_tbl" , site_sort_tbl,
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
				smalloc(cas.site_size*sizeof(unsigned)) ;
			if(link_sort_tbl == 0)
			  ERR.Pointer(cname,fname, "link_sort_tbl"); 
			VRB.Smalloc(cname,fname,
				    "link_sort_tbl" , link_sort_tbl,
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

			VRB.Sfree(cname,fname, "link_sort_tbl", link_sort_tbl);
			sfree(link_sort_tbl);
			VRB.Sfree(cname,fname, "site_sort_tbl", site_sort_tbl);
			sfree(site_sort_tbl) ;

			break ;

		case G_WILSON_HB :

			VRB.Flow(cname,fname,
				 "Converting gauge field order: %s",
				 "G_WILSON_HB -> CANONICAL\n");

			site_sort_tbl = (unsigned *)
				smalloc(cas.vol*sizeof(unsigned)) ;
			if(site_sort_tbl == 0)
			  ERR.Pointer(cname,fname, "site_sort_tbl"); 
			VRB.Smalloc(cname,fname,
				    "site_sort_tbl" , site_sort_tbl,
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
				smalloc(cas.site_size*sizeof(unsigned)) ;
			if(link_sort_tbl == 0)
			  ERR.Pointer(cname,fname, "link_sort_tbl"); 
			VRB.Smalloc(cname,fname,
				    "link_sort_tbl" , link_sort_tbl,
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


			VRB.Sfree(cname,fname, "link_sort_tbl", link_sort_tbl);
			sfree(link_sort_tbl);
			VRB.Sfree(cname,fname, "site_sort_tbl", site_sort_tbl);
			sfree(site_sort_tbl) ;

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

char *fname_fconvert = "Fconvert(V*,StrOrdType,StrOrdType, int)";
char *converting_str = "Converting frm field str ord from %d to %d\n";


//------------------------------------------------------------------
// Fnone::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from)
// Convert None fermion field f_field from -> to
//------------------------------------------------------------------
void Fnone::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from) 
{
	VRB.Func(cname,fname_fconvert);
}

//-------------------------------------------------------------------------
// Fwilson::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from) 
//
// Converts the given fermion field between CANONICAL <-> WILSON.  No 
// other conversions are currently supported in this derived class.
//-------------------------------------------------------------------------

void Fwilson::Fconvert(Vector *f_field,StrOrdType to,StrOrdType from) 
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
// Fstag::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from)
// Convert Staggered fermion field f_field from -> to
//------------------------------------------------------------------
void Fstag::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from) 
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
          FcanonToStag(&cas) ;

        } else if ((from == STAG) && (to == CANONICAL)) {

          VRB.Flow(cname,fname_fconvert, converting_str,
                   int(from), int(to));

          FstagToCanon(&cas) ;

        } else {
                ERR.General(cname,fname_fconvert,
                        "Unsupported fermion conversion from %d to %d\n",
                        int(from), int(to)) ;
        }

        return ;
}




//------------------------------------------------------------------
// Fstag::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from)
// Convert Staggered fermion field f_field from -> to
//------------------------------------------------------------------
void FstagAsqtad::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from) 
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
          FcanonToStag(&cas) ;

        } else if ((from == STAG) && (to == CANONICAL)) {

          VRB.Flow(cname,fname_fconvert, converting_str,
                   int(from), int(to));

          FstagToCanon(&cas) ;

        } else {
                ERR.General(cname,fname_fconvert,
                        "Unsupported fermion conversion from %d to %d\n",
                        int(from), int(to)) ;
        }

        return ;
}




//------------------------------------------------------------------
// Fclover::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from)
// Convert Clover fermion field f_field from -> to
// Since the Clover and Wilson storage orders are the same this 
// routine calls the Wilson conversion routine.
//------------------------------------------------------------------
void Fclover::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from) 
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
// Fdwf::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from)
// Convert dwf fermion field f_field from -> to
//------------------------------------------------------------------
void Fdwf::Fconvert(Vector *f_field, StrOrdType to, StrOrdType from) 
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
	Vector *tmp_f_field = (Vector *) smalloc(f_size * sizeof(Float));
	if(tmp_f_field == 0)
	  ERR.Pointer(cname,fname_fconvert, "tmp_f_field");
	VRB.Smalloc(cname,fname_fconvert, 
		    "tmp_f_field", tmp_f_field, f_size * sizeof(Float));


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
	  moveMem(tmp_field_ptr, field_ptr, f_size * sizeof(Float));

	  // Set odd part
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  for(i=0; i<ls; i++){
	    parity = (i+1) % 2;
	    moveMem(field_ptr, tmp_field_ptr, half_stride * sizeof(Float));
	    field_ptr = field_ptr + half_stride;
	    tmp_field_ptr = tmp_field_ptr + (2 * parity + 1) * half_stride; 
	  }

	  // Set even part
	  tmp_field_ptr = (Float *) tmp_f_field;
	  tmp_field_ptr = tmp_field_ptr + half_stride;
	  for(i=0; i<ls; i++){
	    parity = i % 2;
	    moveMem(field_ptr, tmp_field_ptr, half_stride * sizeof(Float));
	    field_ptr = field_ptr + half_stride;
	    tmp_field_ptr = tmp_field_ptr + (2 * parity + 1) * half_stride; 
	  }

	} else if ((from == WILSON) && (to == CANONICAL)) {
	  
	  VRB.Flow(cname,fname_fconvert, converting_str,
		   int(from), int(to));


	  // copy vector to a buffer for intermediate conversion
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  moveMem(tmp_field_ptr, field_ptr, f_size * sizeof(Float));


	  // convert odd part to intermediate conversion
	  field_ptr = (Float *) f_field;
	  tmp_field_ptr = (Float *) tmp_f_field;
	  for(i=0; i<ls; i++){
	    parity = (i+1) % 2;
	    moveMem(field_ptr, tmp_field_ptr, half_stride * sizeof(Float));
	    field_ptr = field_ptr + (2 * parity + 1) * half_stride; 
	    tmp_field_ptr = tmp_field_ptr + half_stride;
	  }

	  // convert even part to intermediate conversion
	  field_ptr = (Float *) f_field;
	  field_ptr = field_ptr + half_stride;
	  for(i=0; i<ls; i++){
	    parity = i % 2;
	    moveMem(field_ptr, tmp_field_ptr, half_stride * sizeof(Float));
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

	} else {
		ERR.General(cname,fname_fconvert,
			"Unsupported fermion conversion from %d to %d\n",
			int(from), int(to)) ;
	}


	// Free temporary fermion field memory
	//----------------------------------------------------------
	VRB.Sfree(cname,fname_fconvert, "tmp_f_field", tmp_f_field);
	sfree(tmp_f_field);

	return ;

}

CPS_END_NAMESPACE
