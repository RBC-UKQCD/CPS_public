#include <config.h>
#include <util/error.h>
#include <util/lat_data.h>
#include <util/vector.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE
int LatData::DEFAULT_FLAG = 0;

int LatData::Init(int flags, int len, int volume){
	if (status == INITTED)
	ERR.General(cname,"Init()","not allowed to initialize twice\n");
	size = len;
	vol = volume;
	data = (IFloat *)smalloc(sizeof(IFloat)*size*vol);
	status = INITTED;
	return status;
}

IFloat *LatData::Field(int pos, int n){
	IFloat *pointer = data+pos*size+n;
	return pointer;
}

LatData::~LatData(){
 	if (data!= NULL) sfree(data);
}

CPS_END_NAMESPACE
