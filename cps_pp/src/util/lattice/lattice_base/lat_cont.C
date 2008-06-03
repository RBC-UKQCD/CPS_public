#include <config.h>
#include <util/lat_cont.h>
#include <util/error.h>

CPS_START_NAMESPACE
char * LatticeContainer::cname = "LatticeContainer";
LatticeContainer::LatticeContainer(){
	size_t mat_size = GJP.VolNodeSites()*4;
	gauge_p = new Matrix[mat_size];
}
LatticeContainer::~LatticeContainer(){
	delete[] gauge_p;
}

void LatticeContainer::Get(Lattice &lat){
	str_ord = lat.StrOrd();
	lat.CopyGaugeField(gauge_p);
}
void LatticeContainer::Set(Lattice &lat){
	if (str_ord != lat.StrOrd())
	ERR.General(cname,"Set()","Storage ordering of LatticeContainer(%d) doesn't agree with lattice ordering(%d)", str_ord,lat.StrOrd());
	lat.GaugeField(gauge_p);
}
CPS_END_NAMESPACE
