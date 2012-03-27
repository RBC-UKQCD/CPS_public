// Initialize the Nuc3pt argument structure
#include <config.h>
#include <alg/nuc3pt_arg.h>
#include <util/gjp.h>
#include <util/smalloc.h>
CPS_START_NAMESPACE
Nuc3ptArg::Nuc3ptArg(): 
  num_masses(0),
  t_source(0),
  source_inc(0),
  num_src(1),
  t_sink(GJP.Tnodes()*GJP.TnodeSites()/2),
  BoxStart(0),BoxEnd(1),
  theta(0), 
  load_u1_lat(0), 
  u1_gauge_ptr(NULL),
  src_type(POINT),
  DoUnPolarized(1),
  DoUnPolarizedMom(1),
  DoPolarized(1),
  DoPolarizedMom(1),
  //    DoMagMom(0),
  DoHalfFermion(1),
  DoPerPlusAper(0),
  MaxMom2(0),
  DoSS2ptF(1),
  DoGa1Proj(1),
  DoConserved(1),
  num_mult(1)
{
  cname = "Nuc3ptArg" ;
  cg.mass=.1;
  cg.max_num_iter=1000;
  cg.stop_rsd=1e-8;
  cg.RitzMatOper = NONE;
  x[0]=x[1]=x[2]=0;
}

int Nuc3ptArg::NumMasses(){ return num_masses ; }

Float Nuc3ptArg::Mass(int m){ return mass[m] ; }

void Nuc3ptArg::check_args()
{
  char *fname = "check_args" ;
  char *msg = "src location error\n";
  
  if((t_source < 0) || (t_source >= GJP.Tnodes()*GJP.TnodeSites()) ||
     (t_sink  <  0) || (t_sink   >= GJP.Tnodes()*GJP.TnodeSites()))
    ERR.General(cname,fname, "source/sink location error\n");
  
  switch (src_type)
    {
    case BOX: // Check The box source settings
      {
	for(int i(0);i<4;i++)
	  if((BoxStart>=GJP.Nodes(i)*GJP.NodeSites(i)) || (BoxStart<0)||
	     (BoxEnd  >=GJP.Nodes(i)*GJP.NodeSites(i)) || (BoxEnd  <0))
	    ERR.General(cname,fname, "Box source size error\n");
      }
      break ;
    case GAUSS_GAUGE_INV:
    case POINT: // Check The point source settings
      {
	for(int i(0);i<3;i++)
	  if((x[i]>=GJP.Nodes(i)*GJP.NodeSites(i)) || (x[i]<0))
	    ERR.General(cname,fname, msg);
      }
      break ;
    default:
      ERR.General(cname,fname, "Unsupported source");
    }
}
CPS_END_NAMESPACE
