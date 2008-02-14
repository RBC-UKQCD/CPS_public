//--------------------------------------------------------------------
/*!\file
  \brief Implementation of smearing class methods.

  AlgOlegSmear:

  Ape smearing with Thomas' SU(3) projection.

  This provides Ape smearing with Thomas' SU(3) projection, used in
  the static-light matrix elements study by Oleg, Taku, Thomas.
  Needs to set the parameters properly to do the same thing they did:
    asa.coef      = 1;
    asa.orthog    =-1;
    asa.tolerance = ?;
  This class is derived from AlgApeSmear, which uses the build-in SU(3)
  projection that all the other AlgSmear classes use.

  This class is primarily for the test of the new static-light code
  by Yasumichi. For future projects, HYP or stout smearing might be
  better to use.
*/
//--------------------------------------------------------------------
#include <config.h>
#include <math.h>
#include <util/qcdio.h>
#include <alg/alg_smear.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/site.h>
#include <util/link_buffer.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE

AlgOlegSmear::AlgOlegSmear(Lattice&     lat,
			   CommonArg*   ca,
			   ApeSmearArg* asa):
  // does not perform builtin SU(3) projecton, see run().
  AlgApeSmear(lat,ca,asa,0),
  cname("AlgOlegSmear"),
  coef(asa->coef)
{
  cname = "AlgOlegSmear";
  lat_back2 = new Matrix[GJP.VolNodeSites()*4];
  if ( lat_back2 == 0x0 ) { ERR.Pointer(cname, cname,"lat_back2"); }
}
  
AlgOlegSmear::~AlgOlegSmear() {
  delete[] lat_back2;
}

/*!
  If an output file is specified in the CommonArg argument, then
  the smearing coefficients are written to the file.

*/
void AlgOlegSmear::run()
{
  if(common_arg->filename != 0){
    FILE* f = Fopen(common_arg->filename, "a");
    if(!f) ERR.FileA(cname, "run", common_arg->filename);
    Fprintf(f,"AlgOlegSmear: coef = %e \n",coef);
    Fclose(f);
  }

  Lattice& lattice(AlgLattice());

  // backup the original,
  // which needs to be restored when SU(3) projection fails
  lattice.CopyGaugeField(lat_back2);

  AlgSmear::run();

  // do Thomas' SU(3) projection
  Matrix* u;
  u = (Matrix*) lattice.GaugeField();
  int isSingular;
  Site s;
  while( s.LoopsOverNode() ) {
    int offset = s.Index()*4;
    for(int mu=0;mu<4;mu++) {
      if( u->ProjSU3() ) { // returns 1 if the matrix is singular
	printf("Smeared link SU(3) matrix appears to be singular.Substituted original link\n");
	*u = *(lat_back2+offset+mu);
      }
      u++;
    }
  }
}

CPS_END_NAMESPACE
