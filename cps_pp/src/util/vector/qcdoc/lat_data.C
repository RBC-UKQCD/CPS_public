#include <config.h>
#include <util/lat_data.h>
#include <util/vector.h>
#include <qalloc.h>

CPS_START_NAMESPACE
int LatData::DEFAULT_FLAG = QCOMMS;
int LatData::Init(int flags, int len, int volume){
	size = len;
	vol = volume;
	data = (IFloat *)qalloc(flags, sizeof(IFloat)*size*vol);
	if (data == NULL)
	data = (IFloat *)qalloc(DEFAULT_FLAG, sizeof(IFloat)*size*vol);
}

IFloat *LatData::Field(int pos, int n){
	IFloat *pointer = data+pos*size+n;
	return pointer;
}

LatData::~LatData(){
 	if (data!= NULL) qfree(data);
}

CPS_END_NAMESPACE
