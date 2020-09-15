#include <config.h>
#include <util/lattice.h>

#ifndef INCLUDED_LAT_CONT_H
#define INCLUDED_LAT_CONT_H

CPS_START_NAMESPACE
class LatticeContainer{
	private:
		static char *cname;
		Matrix * gauge_p;
		StrOrdType str_ord;
    public:
		LatticeContainer();
		~LatticeContainer();
		
		void Get(Lattice &lat);
		void Set(Lattice &lat);
		Matrix* GaugeField() const{ return gauge_p; } //Added by CK
		
};
CPS_END_NAMESPACE

#endif
