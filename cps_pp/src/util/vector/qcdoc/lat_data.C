#include <config.h>
#include <util/lat_data.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <qalloc.h>

CPS_START_NAMESPACE
int LatData::DEFAULT_FLAG = QCOMMS;
void LatData::Init(int flags, int len, int volume){
	if (status == INITTED)
        ERR.General(cname,"Init()","not allowed to initialize twice\n");
	size = len;
	vol = volume;
	data = (IFloat *)qalloc(flags, sizeof(IFloat)*size*vol);
	if (data == NULL)
	data = (IFloat *)qalloc(DEFAULT_FLAG, sizeof(IFloat)*size*vol);
	if (data == NULL)
	ERR.General("LatData","Init()","out of memory");
        VRB.Flow(cname,"Init()","flags=%x vol=%d data=%p\n",flags,vol,data);
}

#if 0
IFloat *LatData::Field(int pos, int n){
	IFloat *pointer = data+pos*size+n;
	return pointer;
}
#endif

LatData::~LatData(){
 	if (data!= NULL) qfree(data);
}

CPS_END_NAMESPACE
