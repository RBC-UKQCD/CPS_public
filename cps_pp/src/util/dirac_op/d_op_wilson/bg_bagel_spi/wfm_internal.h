#ifndef _WFM_PAB_INTERNAL_H_
#define _WFM_PAB_INTERNAL_H_

/*--------------------------------------------------------------------------*/
/* Wilson assembler routines                                                */
/*--------------------------------------------------------------------------*/

/*
 * FP KERNELS
 */

extern "C" {

void dec_hsu3    (void *psi,void *gauge,void *len,void *tab);
void dec_hsu3_dag(void *psi,void *gauge,void *len,void *tab);

void rec_su3    (void *psi,void *gauge,void *chiin,void *len);
void rec_su3_dag(void *psi,void *gauge,void *chiin,void *len);


  /*
   * PAB Jan 2006, bgl SIMD2 specific variants... these change the 2-spinor ordering a bit
   * and require a modified wfm::interleave_site.
   */
void bgl_dec_hsu3    (void *psi,void *gauge,void *len,void *tab,void *complex_i);
void bgl_dec_hsu3_dag(void *psi,void *gauge,void *len,void *tab,void *complex_i);
void bgl_rec_su3    (void *psi,void *gauge,void *chiin,void *len,void *complex_i);
void bgl_rec_su3_dag(void *psi,void *gauge,void *chiin,void *len,void *complex_i);


/*Feb 2006 Sloppy precision variants - truncate 2spinors to SP for bandwidth*/
void s_dec_hsu3    (void *psi,void *gauge,void *len,void *tab);
void s_dec_hsu3_dag(void *psi,void *gauge,void *len,void *tab);
void s_rec_su3    (void *psi,void *gauge,void *chiin,void *len);
void s_rec_su3_dag(void *psi,void *gauge,void *chiin,void *len);
void s_bgl_dec_hsu3    (void *psi,void *gauge,void *len,void *tab,void *complex_i);
void s_bgl_dec_hsu3_dag(void *psi,void *gauge,void *len,void *tab,void *complex_i);
void s_bgl_rec_su3    (void *psi,void *gauge,void *chiin,void *len,void *complex_i);
void s_bgl_rec_su3_dag(void *psi,void *gauge,void *chiin,void *len,void *complex_i);


void xaxpy_norm (Float *scalep,
		 Float *InOutScale, 
		 Float *Add, 
		 int len,
		 Float *res);

void xaxpy (Float *scalep,Float *InOutScale, Float *Add, int len);


  /*Used to scatter a recv buf onto the appropriate face*/
void qcdoc_face_scatter(Float *TwoSpinor,
		  Float *RcvBuf, 
		  unsigned long*FaceTable,
		  unsigned long);
void qcdoc_s_face_scatter(Float *TwoSpinor,
		  Float *RcvBuf, 
		  unsigned long*FaceTable,
		  unsigned long);
void bgl_face_scatter(Float *TwoSpinor,
		  Float *RcvBuf, 
		  unsigned long*FaceTable,
		  unsigned long);
void bgl_s_face_scatter(Float *TwoSpinor,
		  Float *RcvBuf, 
		  unsigned long*FaceTable,
		  unsigned long);
}


#endif
