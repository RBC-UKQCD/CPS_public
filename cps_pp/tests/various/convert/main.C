#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/convert/main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: main.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.9  2002/12/04 17:16:27  zs
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
//  Revision 1.8  2002/03/11 22:26:57  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.5.2.1  2002/03/08 16:36:28  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.5  2001/08/17 20:03:37  anj
//  Multiple (extra) changes to make the test suite smaller (16CPUs
//  required, not 64) and faster.  Anj
//
//  Revision 1.4  2001/08/16 12:54:19  anj
//  Some fixes follosin the float-> float change, mostly of the (variable
//  anme) float_p -> float_p type.  A few fixes to ensure the test
//  programs use the same level of verbosity throughout, and an update of
//  the regression.pl script to make it more useful. Anj
//
//  Revision 1.3  2001/08/16 10:50:06  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "float".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and float).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:29  anj
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
//  Revision 1.2  2001/05/25 06:16:04  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: main.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/convert/main.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <stdio.h>
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
	do_arg.colors = 3;
	do_arg.beta = 6.0;
	do_arg.verbose_level = -10050402; // = 100 ;

	GJP.Initialize(do_arg);


//-------------------------------------------------------------------------
// Set verbose level
//-------------------------------------------------------------------------

	VRB.Level(GJP.VerboseLevel());

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
