#include <config.h>
#include <util/lat_data.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
//#include <qalloc.h>

CPS_START_NAMESPACE
//int LatData::DEFAULT_FLAG = QCOMMS;
void LatData::Init(LatDataAlloc flags, int len, int volume){
	if (status == INITTED)
        ERR.General(cname,"Init()","not allowed to initialize twice\n");
	size = len;
	vol = volume;
	switch (flags){
		case DEFAULT:
			data = (IFloat *)smalloc(sizeof(IFloat)*size*vol);
			break;
		case FAST:
			data = (IFloat *)fmalloc(sizeof(IFloat)*size*vol);
			break;
		default:
			ERR.General("LatData","Init()","invalid allocation flag");
	}
	if (data == NULL)
	ERR.General("LatData","Init()","out of memory");
        VRB.Flow(cname,"Init()","this=%p flags=%x vol=%d data=%p\n",this,flags,vol,data);

    status = INITTED;
}

#if 0
IFloat *LatData::Field(int pos, int n){
	IFloat *pointer = data+pos*size+n;
	return pointer;
}
#endif

LatData::~LatData(){
 	if (data!= NULL) ffree(data);
}

CPS_END_NAMESPACE
