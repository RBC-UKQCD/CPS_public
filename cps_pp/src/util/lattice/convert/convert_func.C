#include<config.h>
#include<unistd.h>
//#ifndef HAVE_SYNC
//static void inline CPS_NAMESPACE::sync(){}
//#endif
CPS_START_NAMESPACE
/*!\file
  \brief  Functions used by the data layout conversion routines.

  $Id: convert_func.C,v 1.24 2013-04-05 20:05:48 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 20:05:48 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/convert/convert_func.C,v 1.24 2013-04-05 20:05:48 chulwoo Exp $
//  $Id: convert_func.C,v 1.24 2013-04-05 20:05:48 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: convert_func.C,v $
//  $Revision: 1.24 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/lattice/convert/convert_func.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/vector.h>
#include <util/smalloc.h>
#include <util/verbose.h>
#include <util/lattice.h>
#include <comms/sysfunc_cps.h>
////#include <comms/nga_reg.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE

typedef struct ConvertArgStruct {
        unsigned        lx ;
        unsigned        ly ;
        unsigned        lz ;
        unsigned        lt ;
        unsigned        nc ;
	unsigned	ns ;
        unsigned        vol ;
        unsigned        site_size ;
        Float           *start_ptr ;
} CAS, *CAP ;
 
const unsigned CBUF_MODE4 = 0xcca52112;
#ifdef _TARTAN
static Matrix *mp2 = (Matrix *)CRAM_SCRATCH_ADDR+2;
#else
static Matrix tmp_2;
static Matrix *mp2 = &tmp_2;
#endif

void MultStagPhases(CAP cap) ;
void RunGConverter(CAP cap, unsigned *site_tbl, unsigned *link_tbl) ;
void CanonToAnything(CAP cap, StrOrdType new_str_ord) ;

void FcanonToWilson(CAP cap, int number_of_checkerboards) ;
void FwilsonToCanon(CAP cap, int number_of_checkerboards) ;
void FcanonToStag(CAP cap, int number_of_checkerboards) ;
void FstagToCanon(CAP cap, int number_of_checkerboards) ;


char *cname_none = "(none)" ;

#ifdef __cplusplus
extern "C" {
#endif

extern void negate_link(unsigned link_size, IFloat *link) ;
extern void site2cram(void *src, void *dst, unsigned site_size) ;
extern void site2dram(void *src, void *dst, unsigned *link_tbl, unsigned site_size) ;

#ifdef __cplusplus
}
#endif


void MultStagPhases(CAP cap)
{
	unsigned x,y,t,z,link_size ;
	Float *site ;

//-------------------------------------------------------------------------
//  for all that follows, assume euclidean spacetime dimension is 4
//-------------------------------------------------------------------------

	link_size = cap->site_size>>2 ;

	for (t=0; t<cap->lt; t++)
	for (z=0; z<cap->lz; z++)
	for (y=0; y<cap->ly; y++)
	for (x=0; x<cap->lx; x++) {

		site = cap->start_ptr +
			cap->site_size*(x+cap->lx*(y+cap->ly*(z+cap->lz*t)));

#ifdef MILC_COMPATIBILITY  // MILC phase convention - note that this breaks regression tests
		if (t%2) negate_link(link_size, (IFloat *)site);
		if ((t+x)%2) negate_link(link_size, (IFloat *)site+link_size);
		if ((t+x+y)%2) negate_link(link_size, (IFloat *)site+2*link_size) ;
#else	// usual CPS phase definition
		if (x%2) negate_link(link_size, (IFloat *)site+link_size);
		if ((x+y)%2) negate_link(link_size, (IFloat *)site+2*link_size);
		if ((x+y+z)%2) negate_link(link_size, (IFloat *)site+3*link_size) ;
#endif

	}
}


void RunGConverter(CAP cap, unsigned *site_tbl, unsigned *link_tbl)
{
	char *fname = "RunGConverter";

//-------------------------------------------------------------------------
// cram1, cram2 should be in CRAM
//-------------------------------------------------------------------------
//   VRB.Func("",fname);
//  sync();
  const int GSIZE= 72;
//  if(!UniqueID())printf("%s:cap->site_size=%d\n",fname,cap->site_size);
  if (cap->site_size>GSIZE)
  ERR.General("",fname,"cap->site_size(%d)>GSIZE\n",cap->site_size,GSIZE);
  uint32_t vol = cap->vol;

//#ifndef USE_OMP
#if 1
	uint32_t 	low,
			current,
			desired,
			tmp ;
	Float		*cram1,
			*cram2,
			*cram_tmp ;
  Float cram1_stack[GSIZE], cram2_stack[GSIZE];
  cram1 = cram1_stack;
  cram2 = cram2_stack;
	for (low=0; low<vol; low++) {
#else
#pragma omp parallel for default(shared)
	for (low=0; low<vol; low++) {
	uint32_t low,
			current,
			desired,
			tmp ;
	Float		*cram1,
			*cram2,
			*cram_tmp ;
	  Float cram1_stack[GSIZE], cram2_stack[GSIZE];
	  cram1 = cram1_stack;
	  cram2 = cram2_stack;
#endif
//               VRB.Flow("",fname,"low=%d\n",low);
		current = low ;
		desired = *(site_tbl+low) ;

		if (desired & 1) {	// test if site needs conversion

			site2cram(cap->start_ptr+cap->site_size*current,
				cram1, cap->site_size) ;

			while (low != desired>>1) {
				tmp = *(site_tbl + (desired>>1)) ;
				site2cram(cap->start_ptr +
					cap->site_size*(desired>>1),
					cram2, cap->site_size) ;

				site2dram(cram1, cap->start_ptr +
					cap->site_size*(desired>>1),
					link_tbl,
					cap->site_size) ;

				*(site_tbl+(desired>>1)) =
					desired & 0xFFFFFFFE ;

				current = desired >> 1 ;
				desired = tmp ;

				cram_tmp = cram1 ;
				cram1 = cram2 ;
				cram2 = cram_tmp ;
			} ;

			site2dram(cram1,
				cap->start_ptr+cap->site_size*low,
				link_tbl,
				cap->site_size) ;

			*(site_tbl+low) = desired & 0xFFFFFFFE ;
		}
	}

//  sync();
//   VRB.Func("",fname);
//	sfree(cname_none,fname, "cram2", cram2);
//	sfree(cname_none,fname, "cram1", cram1);
}

void CanonToAnything(CAP cap, StrOrdType new_str_ord)
{
	unsigned offset,x,y,z,t,r,row,col,mu ;
	unsigned *site_sort_tbl, *link_sort_tbl ;

	char *fname = "CanonToAnything(CAP,StrOrdType)" ;

	switch (new_str_ord) {
		case STAG :

			VRB.Flow(cname_none,fname,
			"Converting gauge field order: CANONICAL -> STAG\n");

			MultStagPhases(cap) ; // CANONICAL -> STAG

			//---------------------------------------------------
			//  Dagger all links
			//---------------------------------------------------
			{
			  Matrix *p = (Matrix *) cap->start_ptr;
			  int n_links = 4 * cap->vol;
			  
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

			VRB.Flow(cname_none,fname,
			"Converting gauge field order: CANONICAL -> STAG_BLOCK\n");

			MultStagPhases(cap) ; // CANONICAL -> STAG

			//---------------------------------------------------
			//  Dagger all links
			//---------------------------------------------------
			{
			  Matrix *p = (Matrix *) cap->start_ptr;
			  int n_links = 4 * cap->vol;
			  
			  setCbufCntrlReg(4, CBUF_MODE4);
			  
			  for(int i = 0; i < n_links; ++i) {
			    mp2->Dagger((IFloat *)p);
//			    moveMem((IFloat *)(p++), (IFloat *)mp2,
//				    18*sizeof(IFloat));
			    moveFloat((IFloat *)(p++), (IFloat *)mp2, 18);
			  }
			}
			
			//Copy the links into tmp_gauge into STAG_BLOCK ordering
			//Then, copy tmp_gauge back to GaugeField()

			Matrix * tmp_gauge = (Matrix *)
			fmalloc(cname_none,fname,"tmp_gauge",cap->vol*4*sizeof(Matrix));

			int x[4];
			int mu,current,new_index;
			for(mu = 0; mu < 4; mu++)
			  for(x[2]=0;x[2]<cap->lz;x[2]++)
			    for(x[1]=0;x[1]<cap->ly;x[1]++)
			      for(x[0]=0;x[0]<cap->lx;x[0]++)
				for(x[3]=0;x[3]<cap->lt;x[3]++)
				  {
 				    new_index = cap->vol*mu+(x[3]+cap->lt*(x[0]+cap->lx*(x[1]+cap->ly*x[2])))/2 + ((x[3]+x[2]+x[1]+x[0])%2)*cap->vol/2;
				    current = 4*(x[0]+cap->lx*(x[1]+cap->ly*(x[2]+cap->lz*x[3])))+mu;
//				    moveMem(tmp_gauge+new_index,cap->start_ptr+18*current,sizeof(Matrix));
				    moveFloat((Float*)(tmp_gauge+new_index),(Float*)(cap->start_ptr+18*current),18);
				  }
//			moveMem(cap->start_ptr,tmp_gauge,cap->vol*4*sizeof(Matrix));
			moveFloat((Float*)cap->start_ptr,(Float*)tmp_gauge,cap->vol*4*18);
			ffree(cname_none,fname, "tmp_gauge", tmp_gauge);
		  }

			break ;

		case WILSON :
		case G_WILSON_HB :

			if (new_str_ord == WILSON)
				VRB.Flow(cname_none,fname,
				"Converting gauge field order: %s",
				"CANONICAL -> WILSON\n");
			else
				VRB.Flow(cname_none,fname,
				"Converting gauge field order: %s",
				"CANONICAL -> G_WILSON_HB\n");

			site_sort_tbl = (unsigned *)
			fmalloc(cname_none,fname,
				    "site_sort_tbl" , 
				    cap->vol*sizeof(unsigned)) ;



			offset = 0 ;

//-------------------------------------------------------------------------
// Loop over current (CANONICAL) equation for site sequence number
//-------------------------------------------------------------------------

			for (t=0; t<cap->lt; t++)
			for (z=0; z<cap->lz; z++)
			for (y=0; y<cap->ly; y++)
			for (x=0; x<cap->lx; x++) {

//-------------------------------------------------------------------------
// Use desired (WILSON or G_WILSON_HB) equation for site sequence number
// actual WILSON formula is (x-x%2+lx*(y+ly*(z+lz*t))+vol*((x+y+z+t)%2))/2
// LSB = 1 indicates site needs converting
//-------------------------------------------------------------------------

				if (new_str_ord == WILSON)
					*(site_sort_tbl+offset) = (x - x%2
					+ cap->lx*(y+cap->ly*(z+cap->lz*t))
					+ cap->vol*((x+y+z+t)%2)) | 1 ;
				else		// G_WILSON_HB
					*(site_sort_tbl+offset) = 
					z+cap->lz*(y+cap->ly*(x+cap->lx*t))
					<< 1 | 1 ;

				offset++ ;
			}

			link_sort_tbl = (unsigned *)
			fmalloc(cname_none,fname,
				    "link_sort_tbl" , 
				    cap->site_size*sizeof(unsigned)) ;

//-------------------------------------------------------------------------
// Loop over current (CANONICAL) order of reals in site
//-------------------------------------------------------------------------

			offset = 0 ;

			for (mu=0; mu<4; mu++)          // 4 is spacetime dim
			for (row=0; row<cap->nc; row++)
			for (col=0; col<cap->nc; col++)
			for (r=0; r<2; r++) {           // 2 is # cplx comp

//-------------------------------------------------------------------------
// Use desired (WILSON or G_WILSON_HB) equation for reals in site
//-------------------------------------------------------------------------

				if (new_str_ord == WILSON)
					*(link_sort_tbl+offset) =
					r+2*(row+cap->nc*(col+cap->nc*mu)) ;
				else
					*(link_sort_tbl+offset) = col+cap->nc*
					(row+cap->nc*(r+2*((mu+1)%4))) ;
				
				offset++ ;
			}

			RunGConverter(cap, site_sort_tbl, link_sort_tbl) ;

			sfree(cname_none,fname, "link_sort_tbl", link_sort_tbl);
			sfree(cname_none,fname, "site_sort_tbl", site_sort_tbl);

			break ;
		case CANONICAL :
		default :
			break ;
	}
   VRB.FuncEnd("",fname);
}


void FcanonToWilson(CAP cap, int num_chkbds)
{
	char *fname = "FcanonToWilson(CAP, int)";
	VRB.Func(cname_none,fname);

	unsigned idx,x,y,z,t,r,color,spin ;

	unsigned *site_sort_tbl ;
	unsigned *component_sort_tbl ;

	site_sort_tbl = (unsigned *) 
	fmalloc(cname_none,fname,
		    "site_sort_tbl" , 
		    cap->vol * sizeof(unsigned)) ;

//-------------------------------------------------------------------------
// component_sort_tbl should be in CRAM
//-------------------------------------------------------------------------

	component_sort_tbl = 
		(unsigned *) 
	fmalloc(cname_none,fname,
		    "component_sort_tbl" , 
		    cap->site_size * sizeof(unsigned));

//-------------------------------------------------------------------------
// Loop over current (CANONICAL) order of sites in local volume
// Use desired (WILSON) equation for site sequence number
// actual WILSON formula is (x-x%2+lx*(y+ly*(z+lz*t))+vol*((x+y+z+t+1)%2))/2
// LSB = 1 indicates site needs converting
//-------------------------------------------------------------------------

	idx = 0 ;
        if(num_chkbds == 2)  {
	for (t=0; t<cap->lt; t++)
	for (z=0; z<cap->lz; z++)
	for (y=0; y<cap->ly; y++)
	for (x=0; x<cap->lx; x++) {

		*(site_sort_tbl+idx) = (x - x%2 
			+ cap->lx*(y+cap->ly*(z+cap->lz*t)) 
			+ cap->vol*((x+y+z+t+1)%2)) | 1 ;

		idx++ ;
	  }
	}
	else if(num_chkbds == 1) {
	  return; 	// No conversion is necessary for a half lattice
	}
//-------------------------------------------------------------------------
// Loop over current (CANONICAL) order of reals in site
// Use desired (WILSON) equation for reals in site
//-------------------------------------------------------------------------
	idx = 0 ;

	for (spin=0; spin<cap->ns; spin++)
	for (color=0; color<cap->nc; color++)
	for (r=0; r<2; r++) {			// 2 is # cplx comp
	
		*(component_sort_tbl+idx) = r + 2*(color + cap->nc*spin) ;

		idx++ ;
	}

	RunGConverter(cap, site_sort_tbl, component_sort_tbl);

	sfree(cname_none,fname, "component_sort_tbl", component_sort_tbl);
//	sfree(component_sort_tbl);
	sfree(cname_none,fname, "site_sort_tbl", site_sort_tbl);
//	sfree(site_sort_tbl);
}

void FwilsonToCanon(CAP cap, int num_chkbds)
{
	char *fname = "FwilsonToCanon(CAP, int)";
	VRB.Func(cname_none,fname);

	unsigned idx,x,y,z,t,r,color,spin,cb ;

	unsigned *site_sort_tbl ;
	unsigned *component_sort_tbl ;

	site_sort_tbl = (unsigned *) 
	fmalloc(cname_none,fname,
		    "site_sort_tbl" , 
		    cap->vol * sizeof(unsigned)) ;

//-------------------------------------------------------------------------
// component_sort_tbl should be in CRAM
//-------------------------------------------------------------------------

	component_sort_tbl = (unsigned *)
	fmalloc(cname_none,fname, "component_sort_tbl" , 
		    cap->site_size * sizeof(unsigned));

//-------------------------------------------------------------------------
// Loop over current (WILSON) order of sites in local volume
// Use desired (CANONICAL) equation for site sequence number
// LSB = 1 indicates site needs converting
//-------------------------------------------------------------------------

	idx = 0 ;
	if(num_chkbds == 2) {
	for (cb=0; cb<2; cb++)
	for (t=0; t<cap->lt; t++)
	for (z=0; z<cap->lz; z++)
	for (y=0; y<cap->ly; y++)
	for (x=0; x<cap->lx; x++) {
		if ((x+y+z+t+1)%2 == cb) {
			*(site_sort_tbl+idx) = 
			x + cap->lx*(y+cap->ly*(z+cap->lz*t))<<1 | 1 ;
			idx++ ;
		}
	  }
	}
	else if(num_chkbds == 1) {
	  return;	// No conversion is necessary for a half lattice
	}			// in the WILSON storage ordering
//-------------------------------------------------------------------------
// Loop over current (WILSON) order of reals in site
// Use desired (CANONICAL) equation for reals in site
//-------------------------------------------------------------------------

	idx = 0 ;

	for (spin=0; spin<cap->ns; spin++)
	for (color=0; color<cap->nc; color++)
	for (r=0; r<2; r++) {			// 2 is # cplx comp
	
		*(component_sort_tbl+idx) = r + 2*(color + cap->nc*spin) ;

		idx++ ;
	}

	RunGConverter(cap, site_sort_tbl, component_sort_tbl);

	sfree(cname_none,fname, "component_sort_tbl", component_sort_tbl);
//	sfree(component_sort_tbl);
	sfree(cname_none,fname, "site_sort_tbl", site_sort_tbl);
//	sfree(site_sort_tbl);
}

void FcanonToStag(CAP cap, int num_chkbds)
{
        char *fname = "FcanonToStag(CAP)";
        VRB.Func(cname_none ,fname);

        unsigned idx, x, y, z, t, r, color, spin;

        unsigned *site_sort_tbl ;
        unsigned *component_sort_tbl ;

        site_sort_tbl = (unsigned *) 
        fmalloc(cname_none,fname,
                    "site_sort_tbl" , 
                    cap->vol * sizeof(unsigned)) ;

//-------------------------------------------------------------------------
// component_sort_tbl should be in CRAM
//-------------------------------------------------------------------------

        component_sort_tbl =
                (unsigned *) 
        fmalloc(cname_none,fname,
                    "component_sort_tbl" , 
                    cap->site_size * sizeof(unsigned));

//-------------------------------------------------------------------------
// Loop over current (CANONICAL) order of sites in local volume
// Use desired (STAG) equation for site sequence number
// actual STAG formula is
// {(t + lt*(x + lx*(y + ly * z))) + vol * ((x+y+z+t)%2) } / 2
// LSB = 1 indicates site needs converting
//-------------------------------------------------------------------------

        idx = 0 ;

        if(num_chkbds == 2) {
          for (t=0; t < cap->lt; t++)
          for (z=0; z < cap->lz; z++)
          for (y=0; y < cap->ly; y++)
          for (x=0; x < cap->lx; x++) {
            *(site_sort_tbl + idx) = (t + cap->lt * (x + cap->lx *
                (y + cap->ly * z)) + cap->vol*((x+y+z+t)%2)) | 1 ;
            idx++ ;
          }
        }
        else if(num_chkbds == 1) {
          for (t=0; t < cap->lt; t++)
          for (z=0; z < cap->lz; z++)
          for (y=0; y < cap->ly; y++)
          for (x=0; x < cap->lx; x++)
            if( (x+y+z+t)%2 == 0)     {
              *(site_sort_tbl + idx) =
                (t + cap->lt * (x + cap->lx * (y + cap->ly * z))) | 1 ;
              idx++ ;
            }
        }


//-------------------------------------------------------------------------
// Loop over current (CANONICAL) order of reals in site
// Use desired (STAG) equation for reals in site
//-------------------------------------------------------------------------
        idx = 0 ;

        for (spin=0; spin<cap->ns; spin++)
        for (color=0; color<cap->nc; color++)
        for (r=0; r<2; r++) {                   // 2 is # cplx comp

                *(component_sort_tbl+idx) = r + 2*(color + cap->nc*spin) ;

                idx++ ;
        }

        RunGConverter(cap, site_sort_tbl, component_sort_tbl);

        sfree(cname_none,fname, "component_sort_tbl", component_sort_tbl);
//        sfree(component_sort_tbl);

        sfree(cname_none,fname, "site_sort_tbl", site_sort_tbl);
//        sfree(site_sort_tbl);
}

void FstagToCanon(CAP cap, int num_chkbds)
{
        char * fname = "FstagToCanon(CAP)";
        VRB.Func(cname_none, fname);

        unsigned idx, x, y, z, t, cb, r, color, spin;

        unsigned *site_sort_tbl ;
        unsigned *component_sort_tbl ;

        site_sort_tbl = (unsigned *) 
        fmalloc(cname_none,fname,
                    "site_sort_tbl" , 
                    cap->vol * sizeof(unsigned)) ;

//-------------------------------------------------------------------------
// component_sort_tbl should be in CRAM
//-------------------------------------------------------------------------

        component_sort_tbl =
                (unsigned *) 
        fmalloc(cname_none,fname,
                    "component_sort_tbl" , 
                    cap->site_size * sizeof(unsigned));

//-------------------------------------------------------------------------
// Loop over current (WILSON) order of sites in local volume
// Use desired (CANONICAL) equation for site sequence number
// LSB = 1 indicates site needs converting
//-------------------------------------------------------------------------

        idx = 0 ;

        if(num_chkbds == 2)       {
          for (cb = 0; cb < 2; cb++)
          for (z=0; z < cap->lz; z++)
          for (y=0; y < cap->ly; y++)
          for (x=0; x < cap->lx; x++)
          for (t=0; t < cap->lt; t++)
          if( (x+y+z+t)%2 == cb)      {
             *(site_sort_tbl+idx) =
                  (x + cap->lx*(y+cap->ly*(z+cap->lz*t))) << 1 | 1 ;
             idx++ ;
          }
        }
        else if(num_chkbds == 1)       {
          for (z=0; z < cap->lz; z++)
          for (y=0; y < cap->ly; y++)
          for (x=0; x < cap->lx; x++)
          for (t=0; t < cap->lt; t++)
          if( (x+y+z+t)%2 == 0)      {
             *(site_sort_tbl+idx) =
                  (x + cap->lx*(y+cap->ly*(z+cap->lz*t))) | 1 ;
             idx++ ;
          }
        }
//-------------------------------------------------------------------------
// Loop over current (STAG) order of reals in site
// Use desired (CANONICAL) equation for reals in site
//-------------------------------------------------------------------------

        idx = 0 ;

        for (spin=0; spin<cap->ns; spin++)
        for (color=0; color<cap->nc; color++)
        for (r=0; r<2; r++) {                   // 2 is # cplx comp

                *(component_sort_tbl+idx) = r + 2*(color + cap->nc*spin) ;

                idx++ ;
        }

        RunGConverter(cap, site_sort_tbl, component_sort_tbl);

        sfree(component_sort_tbl, cname_none,fname, "component_sort_tbl");
//        sfree(component_sort_tbl);
        sfree(site_sort_tbl,cname_none,fname, "site_sort_tbl");
//        sfree(site_sort_tbl);
}


CPS_END_NAMESPACE
