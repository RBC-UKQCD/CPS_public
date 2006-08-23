#include <config.h>
#include <stdio.h>
#include <util/error.h>
#include <util/lat_data.h>
#include <util/vector.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE
//int LatData::DEFAULT_FLAG = 0;

void LatData::Init(LatDataAlloc flags, int len, int volume){
//	printf("LatData::Init(%p)\n",this);
	if (status != NEW)
	ERR.General(cname,"Init()","not allowed to initialize twice\n");
	size = len;
	vol = volume;
	data = (IFloat *)smalloc(sizeof(IFloat)*size*vol);
	if (data == NULL)
        ERR.General("LatData","Init()","out of memory");
//	printf("size=%d vol=%d data(%p)=%p\n",size,vol,&data,data);
	status = INITTED;
//	return status;
}

#if 0
const IFloat *LatData::Field(int pos, int n){
	IFloat *pointer = data+pos*size+n;
	return pointer;
}
#endif

LatData::~LatData(){
//	printf("LatData::~LatData(%p)\n",this);
 	if (data!= NULL) sfree(data);
}

CPS_END_NAMESPACE
