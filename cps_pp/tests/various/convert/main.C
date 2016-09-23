#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:58:11 $
//  $Header: /space/cvs/cps/cps++/tests/various/convert/main.C,v 1.5 2004/08/18 11:58:11 zs Exp $
//  $Id: main.C,v 1.5 2004/08/18 11:58:11 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: main.C,v $
//  $Revision: 1.5 $
//  $Source: /space/cvs/cps/cps++/tests/various/convert/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/qcdio.h>
#include <stdlib.h>	// exit()
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
#include<alg/do_arg.h>
CPS_START_NAMESPACE

GlobalJobParameter GJP;
LatRanGen LRG;
Verbose VRB;
Error ERR;

main(void)
{
	char *cname = "(none)" ;
	char *fname = "main()";

//-------------------------------------------------------------------------
// Initializes all Global Job Parameters
//-------------------------------------------------------------------------

	DoArg do_arg;

	do_arg.x_node_sites =  2 ;
	do_arg.y_node_sites =  6 ;
	do_arg.z_node_sites = 10 ;
	do_arg.t_node_sites = 14 ;
	do_arg.s_node_sites =  1 ;
	do_arg.x_nodes = 2;
	do_arg.y_nodes = 2;
	do_arg.z_nodes = 2;
	do_arg.t_nodes = 2;
	do_arg.x_bc = BND_CND_PRD;
	do_arg.y_bc = BND_CND_PRD;
	do_arg.z_bc = BND_CND_PRD;
	do_arg.t_bc = BND_CND_APRD;
	do_arg.start_conf_kind = START_CONF_DISORD;
	do_arg.start_seed_kind = START_SEED_FIXED;
	do_arg.beta = 6.0;


	GJP.Initialize(do_arg);


//-------------------------------------------------------------------------
// Set verbose level
//-------------------------------------------------------------------------


	VRB.DeactivateLevel(VERBOSE_RNGSEED_LEVEL);
	VRB.ActivateLevel(VERBOSE_FUNC_LEVEL);
	    

	{
		GwilsonFnone lat ;

		IFloat *base, comp ;

		unsigned x,y,z,t,mu,row,col,r ;
		unsigned offset ;

		base = (IFloat *)lat.GaugeField() ;

		comp = 0.0 ;	// component label
		offset = 0 ;

//-------------------------------------------------------------------------
// Loop over lattice in Canonical order labelling components serially
//-------------------------------------------------------------------------

		for (t=0; t<GJP.TnodeSites(); t++)
		for (z=0; z<GJP.ZnodeSites(); z++)
		for (y=0; y<GJP.YnodeSites(); y++)
		for (x=0; x<GJP.XnodeSites(); x++)
		for (mu=0; mu<4; mu++)
		for (row=0; row<GJP.Colors(); row++)
		for (col=0; col<GJP.Colors(); col++)
		for (r=0; r<2; r++) {
			*(base + offset) = comp ;
			offset++ ;
			comp += 1.0 ;
		}

//-------------------------------------------------------------------------
// Check that Convert traps on CANONICAL -> CANONICAL
//-------------------------------------------------------------------------

		lat.Convert(CANONICAL) ;

//-------------------------------------------------------------------------
// CANONICAL -> STAG
//
// Loop over components in CANONICAL order but generate addresses using
// STAG formula (actually the same) and check to make sure phases are
// correct.
//-------------------------------------------------------------------------

		lat.Convert(STAG) ;

		comp = 0.0 ;

		for (t=0; t<GJP.TnodeSites(); t++)
		for (z=0; z<GJP.ZnodeSites(); z++)
		for (y=0; y<GJP.YnodeSites(); y++)
		for (x=0; x<GJP.XnodeSites(); x++)
		for (mu=0; mu<4; mu++) {

//-------------------------------------------------------------------------
// must compute the STAG phase (eta) here since depends on mu
//-------------------------------------------------------------------------

			IFloat eta ;

			switch (mu) {
				case 1 :
					eta = x%2 ? -1.0 : 1.0 ;
					break ;
				case 2 :
					eta = (x+y)%2 ? -1.0 : 1.0 ;
					break ;
				case 3 :
					eta = (x+y+z)%2 ? -1.0 : 1.0 ;
					break ;
				case 0 :
				default :
					eta = 1.0 ;
					break ;
			}

			for (row=0; row<GJP.Colors(); row++)
			for (col=0; col<GJP.Colors(); col++)
			for (r=0; r<2; r++) {
				offset = z + GJP.ZnodeSites()*t ;
				offset = y + GJP.YnodeSites()*offset ;
				offset = x + GJP.XnodeSites()*offset ;
				offset = row + GJP.Colors()*(mu + 4*offset) ;
				offset = r + 2*(col + GJP.Colors()*offset) ;
				
				if (*(base+offset) != eta*comp)
					VRB.Warn(cname,fname,
					"(%d,%d,%d,%d)(%d)\t%f != %f*%f\n",
					x,y,z,t,mu,*(base+offset),eta,comp);

				comp += 1.0 ;
			}
		}

//-------------------------------------------------------------------------
// check that Convert() traps on STAG -> STAG
//-------------------------------------------------------------------------

		lat.Convert(STAG) ;

//-------------------------------------------------------------------------
// STAG -> CANONICAL
//
// Loop over components in CANONICAL order and check labelling
//-------------------------------------------------------------------------

		lat.Convert(CANONICAL) ;

		offset = 0 ;
		comp = 0.0 ;

		for (t=0; t<GJP.TnodeSites(); t++)
		for (z=0; z<GJP.ZnodeSites(); z++)
		for (y=0; y<GJP.YnodeSites(); y++)
		for (x=0; x<GJP.XnodeSites(); x++)
		for (mu=0; mu<4; mu++)
		for (row=0; row<GJP.Colors(); row++)
		for (col=0; col<GJP.Colors(); col++)
		for (r=0; r<2; r++) {
			if (*(base + offset) != comp)
				VRB.Warn(cname,fname,
					"(%d,%d,%d,%d)(%d)\t%f != %f\n",
					x,y,z,t,mu,*(base+offset),comp);
			offset++ ;
			comp += 1.0 ;
		}

//-------------------------------------------------------------------------
// CANONICAL -> WILSON
//
// Loop over components in CANONICAL order but generate addresses using
// WILSON formula and check to make sure components are correct.
//-------------------------------------------------------------------------

		lat.Convert(WILSON) ;

		offset = 0 ;
		comp = 0.0 ;

		for (t=0; t<GJP.TnodeSites(); t++)
		for (z=0; z<GJP.ZnodeSites(); z++)
		for (y=0; y<GJP.YnodeSites(); y++)
		for (x=0; x<GJP.XnodeSites(); x++)
		for (mu=0; mu<4; mu++)
		for (row=0; row<GJP.Colors(); row++)
		for (col=0; col<GJP.Colors(); col++)
		for (r=0; r<2; r++) {
			offset = z + GJP.ZnodeSites()*t ;
			offset = y + GJP.YnodeSites()*offset ;
			offset = x - x%2 + GJP.XnodeSites()*offset ;
			offset += GJP.VolNodeSites()*((x+y+t+z)%2) ;

			// next line is tricky, actually 4*(offset/2)

			offset = col + GJP.Colors()*(mu + 2*offset) ;
			offset = r + 2*(row + GJP.Colors()*offset) ;

			if (*(base+offset) != comp)
				VRB.Warn(cname,fname,
				"(%d,%d,%d,%d)(%d)\t%u\t%f != %f\n",
				x,y,z,t,mu,offset,*(base+offset),comp);

			comp += 1.0 ;
		}

//-------------------------------------------------------------------------
// check that Convert() traps on WILSON -> WILSON
//-------------------------------------------------------------------------

		lat.Convert(WILSON) ;

//-------------------------------------------------------------------------
// WILSON -> CANONICAL
//
// Loop over components in CANONICAL order and check labelling
//-------------------------------------------------------------------------

		lat.Convert(CANONICAL) ;

		offset = 0 ;
		comp = 0.0 ;

		for (t=0; t<GJP.TnodeSites(); t++)
		for (z=0; z<GJP.ZnodeSites(); z++)
		for (y=0; y<GJP.YnodeSites(); y++)
		for (x=0; x<GJP.XnodeSites(); x++)
		for (mu=0; mu<4; mu++)
		for (row=0; row<GJP.Colors(); row++)
		for (col=0; col<GJP.Colors(); col++)
		for (r=0; r<2; r++) {
			if (*(base + offset) != comp)
				VRB.Warn(cname,fname,
					"(%d,%d,%d,%d)(%d)\t%f != %f\n",
					x,y,z,t,mu,*(base+offset),comp);
			offset++ ;
			comp += 1.0 ;
		}

//-------------------------------------------------------------------------
// CANONICAL -> G_WILSON_HB
//
// Loop over components in CANONICAL order but generate addresses using
// G_WILSON_HB formula and check to make sure components are correct.
//-------------------------------------------------------------------------

		lat.Convert(G_WILSON_HB) ;

		offset = 0 ;
		comp = 0.0 ;

		for (t=0; t<GJP.TnodeSites(); t++)
		for (z=0; z<GJP.ZnodeSites(); z++)
		for (y=0; y<GJP.YnodeSites(); y++)
		for (x=0; x<GJP.XnodeSites(); x++)
		for (mu=0; mu<4; mu++)
		for (row=0; row<GJP.Colors(); row++)
		for (col=0; col<GJP.Colors(); col++)
		for (r=0; r<2; r++) {
			offset = x + GJP.XnodeSites()*t ;
			offset = y + GJP.YnodeSites()*offset ;
			offset = z + GJP.ZnodeSites()*offset ;

			offset = r + 2*((mu+1)%4 + 4*offset) ;

			offset = row + GJP.Colors()*offset ;
			offset = col + GJP.Colors()*offset ;


			if (*(base+offset) != comp)
				VRB.Warn(cname,fname,
				"(%d,%d,%d,%d)(%d)\t%u\t%f != %f\n",
				x,y,z,t,mu,offset,*(base+offset),comp);

			comp += 1.0 ;
		}

//-------------------------------------------------------------------------
// check that Convert() traps on G_WILSON_HB -> G_WILSON_HB
//-------------------------------------------------------------------------

		lat.Convert(G_WILSON_HB) ;

//-------------------------------------------------------------------------
// G_WILSON_HB -> CANONICAL
//
// Loop over components in CANONICAL order and check labelling
//-------------------------------------------------------------------------

		lat.Convert(CANONICAL) ;

		offset = 0 ;
		comp = 0.0 ;

		for (t=0; t<GJP.TnodeSites(); t++)
		for (z=0; z<GJP.ZnodeSites(); z++)
		for (y=0; y<GJP.YnodeSites(); y++)
		for (x=0; x<GJP.XnodeSites(); x++)
		for (mu=0; mu<4; mu++)
		for (row=0; row<GJP.Colors(); row++)
		for (col=0; col<GJP.Colors(); col++)
		for (r=0; r<2; r++) {
			if (*(base + offset) != comp)
				VRB.Warn(cname,fname,
					"(%d,%d,%d,%d)(%d)\t%f != %f\n",
					x,y,z,t,mu,*(base+offset),comp);
			offset++ ;
			comp += 1.0 ;
		}

	} // end scope for GwilsonFnone

	{
		GwilsonFwilson lat ;

		IFloat *f1, *f2, comp ;

		unsigned idx,r,color,spin,x,y,z,t,vol ;

		vol = 2 * lat.Colors() * lat.SpinComponents()
			* GJP.VolNodeSites() ;

		f1 = (IFloat *)smalloc((size_t)vol*sizeof(IFloat)) ;
		if (f1 == 0)
			ERR.Pointer(cname,fname, "f1") ;
		VRB.Smalloc(cname,fname,"f1",f1,vol) ;

		f2 = (IFloat *)smalloc((size_t)vol*sizeof(IFloat)) ;
		if (f2 == 0)
			ERR.Pointer(cname,fname, "f2") ;
		VRB.Smalloc(cname,fname,"f2",f2,vol) ;

//-------------------------------------------------------------------------
// Loop over fermion fields in Canonical order, labelling serially
//-------------------------------------------------------------------------

		comp = 0.0 ;
		idx = 0 ;

		for (t=0; t<GJP.TnodeSites(); t++)
		for (z=0; z<GJP.ZnodeSites(); z++)
		for (y=0; y<GJP.YnodeSites(); y++)
		for (x=0; x<GJP.XnodeSites(); x++)
		for (spin=0; spin<lat.SpinComponents(); spin++)
		for (color=0; color<lat.Colors(); color++)
		for (r=0; r<2; r++) {
			*(f1+idx) = comp ;
			*(f2+idx) = comp ;
			idx++ ;
			comp += 1.0 ;
		}
		
//-------------------------------------------------------------------------
// Check that Fconvert traps on CANONICAL -> CANONICAL
//-------------------------------------------------------------------------

		lat.Convert(CANONICAL, (Vector *) f1, (Vector *) f2) ;

//-------------------------------------------------------------------------
// CANONICAL -> WILSON
//
// Loop over components in CANONICAL order but generate addresses using
// WILSON formula and check to make sure components are correct.
//-------------------------------------------------------------------------

		lat.Convert(WILSON, (Vector *) f1, (Vector *) f2) ;

		idx = 0 ;
		comp = 0.0 ;

		for (t=0; t<GJP.TnodeSites(); t++)
		for (z=0; z<GJP.ZnodeSites(); z++)
		for (y=0; y<GJP.YnodeSites(); y++)
		for (x=0; x<GJP.XnodeSites(); x++)
		for (spin=0; spin<lat.SpinComponents(); spin++)
		for (color=0; color<lat.Colors(); color++)
		for (r=0; r<2; r++) {
			idx = z + GJP.ZnodeSites()*t ;
			idx = y + GJP.YnodeSites()*idx ;
			idx = x - (x%2) + GJP.XnodeSites()*idx ;
			idx += GJP.VolNodeSites()*((x+y+t+z+1)%2) ;
			idx >>= 1 ;
			idx = spin + lat.SpinComponents()*idx ;
			idx = color + lat.Colors()*idx ;
			idx = r + 2*idx ;

			if (*(f1+idx) != comp)
				VRB.Warn(cname,fname,
				"f1(%d,%d,%d,%d)\t%u\t%f != %f\n",
				x,y,z,t,idx,*(f1+idx),comp) ;

			if (*(f2+idx) != comp)
				VRB.Warn(cname,fname,
				"f2(%d,%d,%d,%d)\t%u\t%f != %f\n",
				x,y,z,t,idx,*(f2+idx),comp) ;

			comp += 1.0 ;

		}

//-------------------------------------------------------------------------
// check that Convert() traps on WILSON -> WILSON
//-------------------------------------------------------------------------

		lat.Convert(WILSON, (Vector *) f1, (Vector *) f2) ;

//-------------------------------------------------------------------------
// WILSON -> CANONICAL
//
// Loop over components in CANONICAL order and check labelling
//-------------------------------------------------------------------------

		lat.Convert(CANONICAL, (Vector *) f1, (Vector *) f2) ;

		idx = 0 ;
		comp = 0.0 ;

		for (t=0; t<GJP.TnodeSites(); t++)
		for (z=0; z<GJP.ZnodeSites(); z++)
		for (y=0; y<GJP.YnodeSites(); y++)
		for (x=0; x<GJP.XnodeSites(); x++)
		for (spin=0; spin<lat.SpinComponents(); spin++)
		for (color=0; color<lat.Colors(); color++)
		for (r=0; r<2; r++) {

			if (*(f1+idx) != comp)
				VRB.Warn(cname,fname,
				"f1(%d,%d,%d,%d)\t%u\t%f != %f\n",
				x,y,z,t,idx,*(f1+idx),comp) ;

			if (*(f2+idx) != comp)
				VRB.Warn(cname,fname,
				"f2(%d,%d,%d,%d)\t%u\t%f != %f\n",
				x,y,z,t,idx,*(f2+idx),comp) ;

			idx++ ;
			comp += 1.0 ;
		}

	} // end scope for GwilsonFwilson
}
CPS_END_NAMESPACE
