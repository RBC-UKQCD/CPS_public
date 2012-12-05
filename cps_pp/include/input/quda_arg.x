enum CudaPrecision {
  CUDA_HALF_PRECISION,
  CUDA_SINGLE_PRECISION,
  CUDA_DOUBLE_PRECISION
};

enum CudaReconstructType {
  CUDA_RECONSTRUCT_NO, // store all 18 real numbers epxlicitly
  CUDA_RECONSTRUCT_8, // reconstruct from 8 real numbers
  CUDA_RECONSTRUCT_12 // reconstruct from 12 real numbers
};


/*! A structure to hold the QUDA solver parameters.*/
class QudaArg {
  CudaPrecision gauge_prec; /*!< Precision of gauge field. */
  CudaPrecision spinor_prec; /*!< Precision of spinors. */
  CudaReconstructType reconstruct; /*!< How the gauge field is stored on CUDA device. */
  CudaPrecision gauge_prec_sloppy; /*!< Precision of gauge field. */
  CudaPrecision spinor_prec_sloppy; /*!< Precision of spinors. */
  CudaReconstructType reconstruct_sloppy; /*!< How the gauge field is stored on CUDA device. */
  Float reliable_delta; /*!< delta parameter for reliableupdates. */
  int max_restart; /*< Maximum number of restarts. */
  int device; /*!< Which CUDA device to use. */
  memfun QudaArg();
};
