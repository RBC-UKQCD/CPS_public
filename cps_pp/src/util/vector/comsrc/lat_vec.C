//------------------------------------------------------------------
/*!\file
  \brief  Implementation of LatData class methods.

  $Id: lat_vec.C,v 1.4 2004-09-02 16:58:11 zs Exp $
*/
//------------------------------------------------------------------

#include <config.h>
#include <util/lat_data.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE

LatData::LatData(const LatData &lat){
    ERR.General("LatData","LatData(&LatData)","Copy constructor not allowed");

}

LatData &LatData::operator=(const LatData &lat){
    ERR.General("LatData","LatData(&LatData)","Copy constructor not allowed");
}

void LatVector::Init(int flag, int n_vec, int vol)
{
//	printf("LatVector::Init(%p)\n",this);
	vec_size = 2*GJP.Colors();
	if (vol == 0) vol = GJP.VolNodeSites();
	LatData::Init(flag,n_vec*vec_size,vol);
//	printf("data=%p\n",data);
}


LatVector::~LatVector(){
//	printf("LatVector::~LatVector(%p)\n",this);
}

Vector *LatVector::Vec(int pos, int vec_row){
//	printf("data=%p\n",data);
	Vector *pointer = (Vector *)(data+pos*size+vec_row*vec_size);
	return pointer;
}

Float *LatVector::Field(int pos, int vec_row,int n){
	Float *pointer = (Float *)(data+pos*size+vec_row*vec_size+n);
	return pointer;
}

LatMatrix::LatMatrix(int flag, int n_vec, int vol)
: LatData(){
//	printf("LatMatrix::LatMatrix()\n");
	mat_size = 2*GJP.Colors()*GJP.Colors();
	if (vol == 0) vol = GJP.VolNodeSites();
	LatData::Init(flag, n_vec*mat_size,vol);
}
LatMatrix::~LatMatrix(){
//	printf("LatMatrix::~LatMatrix()\n");
}

Matrix *LatMatrix::Mat(int pos, int vec_row){
	Matrix *pointer = (Matrix *)(data+pos*size+vec_row*mat_size);
	return pointer;
}

Float *LatMatrix::Field(int pos, int mat_row, int n){
	Float *pointer = (Float *)(data+pos*size+mat_row*mat_size+n);
	return pointer;
}

CPS_END_NAMESPACE
