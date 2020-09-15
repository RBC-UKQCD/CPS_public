#ifndef _GSL_WRAPPER
#define _GSL_WRAPPER

//#include<gsl/gsl_blas.h>
#include<gsl/gsl_blas.h>

CPS_START_NAMESPACE

template<typename mf_Float>
struct gsl_wrapper{
};

template<>
struct gsl_wrapper<float>{
  typedef gsl_complex_float complex;
  typedef gsl_matrix_complex_float matrix_complex;
  typedef gsl_vector_complex_float vector_complex;
  typedef gsl_block_complex_float_struct block_complex_struct;
  typedef gsl_matrix_complex_float_const_view matrix_complex_const_view;

  inline static matrix_complex* matrix_complex_alloc(const size_t i,const size_t j){ return gsl_matrix_complex_float_alloc(i,j); }
  inline static void matrix_complex_set(matrix_complex * m, const size_t i, const size_t j, const complex x){ gsl_matrix_complex_float_set(m,i,j,x); }
  inline static void matrix_complex_set_zero(matrix_complex *m){ gsl_matrix_complex_float_set_zero(m); }
  inline static int matrix_complex_set_row(matrix_complex * m, const size_t i, const vector_complex * v){ return gsl_matrix_complex_float_set_row(m,i,v); }

  inline static complex matrix_complex_get(const matrix_complex * m, const size_t i, const size_t j){ return gsl_matrix_complex_float_get(m,i,j); }
  inline static complex* matrix_complex_ptr(matrix_complex *m, const size_t i, const size_t j){ return gsl_matrix_complex_float_ptr(m,i,j); }

  inline static void matrix_complex_free (matrix_complex * m){ gsl_matrix_complex_float_free(m); }

  inline static vector_complex* vector_complex_alloc(const size_t i){ return gsl_vector_complex_float_alloc(i); }
  inline static complex vector_complex_get (const vector_complex * v, const size_t i){ return gsl_vector_complex_float_get(v,i); }
  inline static void vector_complex_free (vector_complex * m){ gsl_vector_complex_float_free(m); }
  inline static complex* vector_complex_ptr(vector_complex * v, const size_t i){   gsl_vector_complex_float_ptr(v,i); }


  inline static int blas_gemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, const complex alpha, 
			       const matrix_complex * A, const matrix_complex * B, const complex beta, matrix_complex * C){
    return gsl_blas_cgemm(TransA,TransB,alpha,A,B,beta,C);
  }
  inline static int blas_gemv (CBLAS_TRANSPOSE_t TransA, const complex alpha, const matrix_complex * A, const vector_complex * x, const complex beta, vector_complex * y){
    return gsl_blas_cgemv(TransA,alpha,A,x,beta,y);
  }
  inline static matrix_complex_const_view matrix_complex_const_submatrix(const matrix_complex * m, 
									 const size_t i, const size_t j, 
									 const size_t n1, const size_t n2){
    return gsl_matrix_complex_float_const_submatrix(m,i,j,n1,n2);
  }
  //Compute the scalar product x^T y for the vectors x and y, returning the result in result. 
  inline static int blas_dotu(const vector_complex * x, const vector_complex * y, complex * dotc){
    return gsl_blas_cdotu(x,y,dotc);
  }
  //Compute the complex conjugate scalar product x^\dagger y for the vectors x and y, returning the result in result. 
  inline static int blas_dotc(const vector_complex * x, const vector_complex * y, complex * dotc){
    return gsl_blas_cdotc(x,y,dotc);
  }

};

template<>
struct gsl_wrapper<double>{
  typedef gsl_complex complex;
  typedef gsl_matrix_complex matrix_complex;
  typedef gsl_vector_complex vector_complex;
  typedef gsl_block_complex_struct block_complex_struct;
  typedef gsl_matrix_complex_const_view matrix_complex_const_view;

  inline static matrix_complex* matrix_complex_alloc(const size_t i,const size_t j){ return gsl_matrix_complex_alloc(i,j); }
  inline static void matrix_complex_set(matrix_complex * m, const size_t i, const size_t j, const complex x){ gsl_matrix_complex_set(m,i,j,x); }
  inline static void matrix_complex_set_zero(matrix_complex *m){ gsl_matrix_complex_set_zero(m); }
  inline static int matrix_complex_set_row(matrix_complex * m, const size_t i, const vector_complex * v){ return gsl_matrix_complex_set_row(m,i,v); }

  inline static complex matrix_complex_get(const matrix_complex * m, const size_t i, const size_t j){ return gsl_matrix_complex_get(m,i,j); }
  inline static complex* matrix_complex_ptr(matrix_complex *m, const size_t i, const size_t j){ return gsl_matrix_complex_ptr(m,i,j); }
  inline static void matrix_complex_free (matrix_complex * m){ gsl_matrix_complex_free(m); }

  inline static vector_complex* vector_complex_alloc(const size_t i){ return gsl_vector_complex_alloc(i); }
  inline static complex vector_complex_get (const vector_complex * v, const size_t i){ return gsl_vector_complex_get(v,i); }
  inline static void vector_complex_free (vector_complex * m){ gsl_vector_complex_free(m); }
  inline static complex* vector_complex_ptr(vector_complex * v, const size_t i){   gsl_vector_complex_ptr(v,i); }

  inline static int blas_gemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, const complex alpha, 
			       const matrix_complex * A, const matrix_complex * B, const complex beta, matrix_complex * C){
    return gsl_blas_zgemm(TransA,TransB,alpha,A,B,beta,C);
  }
  inline static int blas_gemv (CBLAS_TRANSPOSE_t TransA, const complex alpha, const matrix_complex * A, const vector_complex * x, const complex beta, vector_complex * y){
    return gsl_blas_zgemv(TransA,alpha,A,x,beta,y);
  }
  inline static matrix_complex_const_view matrix_complex_const_submatrix(const matrix_complex * m, 
									 const size_t i, const size_t j, 
									 const size_t n1, const size_t n2){
    return gsl_matrix_complex_const_submatrix(m,i,j,n1,n2);
  }
  //Compute the scalar product x^T y for the vectors x and y, returning the result in result.
  inline static int blas_dotu(const vector_complex * x, const vector_complex * y, complex * dotc){
    return gsl_blas_zdotu(x,y,dotc);
  }
  //Compute the complex conjugate scalar product x^\dagger y for the vectors x and y, returning the result in result. 
  inline static int blas_dotc(const vector_complex * x, const vector_complex * y, complex * dotc){
    return gsl_blas_zdotc(x,y,dotc);
  }
};

CPS_END_NAMESPACE


#endif
