//------------------------------------------------------------------
/*!\file
  \brief  Implementation of LatData class methods.

  $Id: lat_vec.C,v 1.7 2006-07-13 20:14:44 chulwoo Exp $
*/
//------------------------------------------------------------------

#include <config.h>
#include <util/lat_data.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/smalloc.h>

CPS_START_NAMESPACE

char *LatData::cname ="LatData";
char *LatMatrix::cname ="LatMatrix";

LatData::LatData(){
	VRB.Flow(cname,"LatData()","this=%p\n",this);
	status = NEW;
	data = NULL;
}

LatData::LatData(const LatData &lat){
    ERR.General("LatData","LatData(&LatData)","Copy constructor not allowed");
}

LatData &LatData::operator=(const LatData &lat){
    ERR.General("LatData","LatData(&LatData)","Copy constructor not allowed");
}

LatMatrix::LatMatrix(const LatMatrix &lat) :LatData(){
    ERR.General("LatMatrix","LatMatrix(&LatMatrix)","Copy constructor not allowed");
}

LatMatrix &LatMatrix::operator=(const LatMatrix &lat){
    ERR.General("LatMatrix","LatMatrix(&LatMatrix)","Copy constructor not allowed");
}

int LatData::Check(const LatData &lat_data){
    if (size != lat_data.size)
      ERR.General(cname,"LatData::Check()","size does not aggree\n");
    if (vol != lat_data.vol)
      ERR.General(cname,"LatData::Check()","vol does not aggree\n");
    return 0;
}

void LatVector::Init(LatDataAlloc flag, int n_vec, int vol)
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
LatMatrix::LatMatrix(){
	const char *fname = "LatMatrix()";
	VRB.Flow(cname,fname,"this(%p) entered\n",this);
	Init(DEFAULT,1,0);
}

void LatMatrix::Init(LatDataAlloc flag, int n_vec, int vol)
{
	const char *fname = "Init(f,i,i)";
	VRB.Flow(cname,fname,"this(%p) entered\n",this);
	mat_size = 2*GJP.Colors()*GJP.Colors();
        n_mat = n_vec;
	if (vol == 0) vol = GJP.VolNodeSites();
	LatData::Init(flag, n_vec*mat_size,vol);
        VRB.Flow(cname,fname,"%p::data(%p)=%p\n",this,&data,data);
	VRB.FuncEnd(cname,fname);
}
LatMatrix::~LatMatrix(){
//	printf("LatMatrix::~LatMatrix()\n");
}

Matrix *LatMatrix::Mat(int pos, int vec_row){
//        printf("%p::data(%p)=%p\n",this,&data,data);
	Matrix *pointer = (Matrix *)(data+pos*size+vec_row*mat_size);
	return pointer;
}

Float *LatMatrix::Field(int pos, int mat_row, int n){
	Float *pointer = (Float *)(data+pos*size+mat_row*mat_size+n);
	return pointer;
}

int LatMatrix::Check(const LatMatrix &lat_mat){
    return LatData::Check(lat_mat);
}

void LatMatrix::operator= (IFloat c){
	const char *fname = "operator=(F)";
//	VRB.Func(cname,fname);
	Matrix *pointer = (Matrix *)data;
//        VRB.Flow(cname,fname,"%p::data(%p)=%p\n",this,&data,data);
	for(int i =0;i<vol*n_mat;i++){
		*(pointer+i) = c;
	}
//	VRB.FuncEnd(cname,fname);
}

void LatMatrix::MulDag (LatMatrix &c){
	LatMatrix::Check(c);
	Matrix *pointer = (Matrix *)data;
	const Matrix *c_mat_p = c.Mat();
        Matrix tmp1,tmp2;
	for(int i =0;i<vol*n_mat;i++){
		tmp1 = (*(pointer+i));
		tmp2.Dagger(*(c_mat_p+i));
		(pointer+i)->DotMEqual( tmp2,tmp1);
	}
}

void LatMatrix::TrLessAntiHermMatrix(){
	Matrix *pointer = (Matrix *)data;
	for(int i =0;i<vol*n_mat;i++){
		(pointer+i)->TrLessAntiHermMatrix();
	}
}

void LatMatrix::operator+= (LatMatrix & lat_mat){
	Check(lat_mat);
	Matrix *pointer = (Matrix *)data;
	const Matrix *lat_mat_p = lat_mat.Mat();
#if TARGET== QCDOC
        Float scale=1.;
        vaxpy3_m(pointer,&scale,pointer,lat_mat_p,n_mat*vol*3);
#else
	for(int i =0;i<vol*n_mat;i++){
		*(pointer+i) += *(lat_mat_p+i);
	}
#endif
//	return *this;
}


CPS_END_NAMESPACE
