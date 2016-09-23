#ifndef _QLA_COMPLEX_H
#define _QLA_COMPLEX_H

typedef int QLA_Int;

typedef float QLA_F_Real;
typedef double QLA_D_Real;
typedef long double QLA_Q_Real;

typedef struct {
   QLA_F_Real real;	
   QLA_F_Real imag;
} QLA_F_Complex;

typedef struct {
   QLA_D_Real real;	
   QLA_D_Real imag;
} QLA_D_Complex;

typedef struct {
   QLA_Q_Real real;	
   QLA_Q_Real imag;
} QLA_Q_Complex;

/* Translation of generic to specific types */
/* Works only if QLA_Precision is defined */
/* (Long double precision is of very limited use and is not accessible
   through a generic choice.) */

#if ( QLA_Precision == 'F' )

typedef QLA_F_Real QLA_Real;
typedef QLA_F_Complex QLA_Complex;

#elif ( QLA_Precision == 'D' )

typedef QLA_D_Real QLA_Real;
typedef QLA_D_Complex QLA_Complex;

#elif ( QLA_Precision == 'Q' )

typedef QLA_Q_Real QLA_Real;
typedef QLA_Q_Complex QLA_Complex;

#endif

/* Accessors for real and imaginary parts */

#define QLA_real(a) (a).real
#define QLA_imag(a) (a).imag

/* Macros for complex arithmetic */

#define QLA_norm2_c(a) (QLA_real(a)*QLA_real(a) + QLA_imag(a)*QLA_imag(a))
#define QLA_norm_c(a)  sqrt(QLA_norm2_c(a))
#define QLA_arg_c(a)   atan2((double)QLA_imag(a), (double)QLA_real(a))
			
/* Unary operations */


#define QLA_c_eq_r(c,a)    {QLA_real(c) = a;\
                            QLA_imag(c) = 0.;}

#define QLA_c_eq_c(c,a)    {QLA_real(c) =  QLA_real(a);\
                            QLA_imag(c) =  QLA_imag(a);}

#define QLA_c_eqm_r(c,a)   {QLA_real(c) = -(a);\
                            QLA_imag(c) = 0.;}

#define QLA_c_eqm_c(c,a)   {QLA_real(c) = -QLA_real(a);\
                            QLA_imag(c) = -QLA_imag(a);}

#define QLA_c_peq_r(c,a)    QLA_real(c) += (a);
			    

#define QLA_c_peq_c(c,a)   {QLA_real(c) += QLA_real(a);\
                            QLA_imag(c) += QLA_imag(a);}

#define QLA_c_meq_r(c,a)    QLA_real(c) -= (a);
			    

#define QLA_c_meq_c(c,a)   {QLA_real(c) -= QLA_real(a);\
                            QLA_imag(c) -= QLA_imag(a);}

#define QLA_c_eq_ca(c,a)   {QLA_real(c)  =  QLA_real(a);\
                            QLA_imag(c)  = -QLA_imag(a);}

#define QLA_c_peq_ca(c,a)  {QLA_real(c) +=  QLA_real(a);\
                            QLA_imag(c) -=  QLA_imag(a);}

#define QLA_c_meq_ca(c,a)  {QLA_real(c) -=  QLA_real(a);\
                            QLA_imag(c) +=  QLA_imag(a);}

#define QLA_c_eqm_ca(c,a)  {QLA_real(c) =  -QLA_real(a);\
                            QLA_imag(c) =   QLA_imag(a);}

#define QLA_r_eq_Re_c(c,a)  c = QLA_real(a);
			    

#define QLA_r_eq_Im_c(c,a)  c = QLA_imag(a);


#define QLA_r_peq_Re_c(c,a) c += QLA_real(a);


#define QLA_r_peq_Im_c(c,a) c += QLA_imag(a);


#define QLA_r_meq_Re_c(c,a) c -= QLA_real(a);


#define QLA_r_meq_Im_c(c,a) c -= QLA_imag(a);


#define QLA_r_eqm_Re_c(c,a) c = -QLA_real(a);


#define QLA_r_eqm_Im_c(c,a) c = -QLA_imag(a);


#define QLA_c_eq_ic(c,a)   {QLA_real(c) = -QLA_imag(a);\
                            QLA_imag(c) =  QLA_real(a);}

#define QLA_c_peq_ic(c,a)  {QLA_real(c) -= QLA_imag(a);\
                            QLA_imag(c) += QLA_real(a);}

#define QLA_c_meq_ic(c,a)  {QLA_real(c) += QLA_imag(a);\
                            QLA_imag(c) -= QLA_real(a);}

#define QLA_c_eqm_ic(c,a)  {QLA_real(c) =  QLA_imag(a);\
                            QLA_imag(c) = -QLA_real(a);}


/* Binary operations */

#define QLA_c_eq_c_plus_c(c,a,b)  {QLA_real(c) = QLA_real(a) + QLA_real(b);\
                                   QLA_imag(c) = QLA_imag(a) + QLA_imag(b);}

#define QLA_c_eqm_c_plus_c(c,a,b)  {QLA_real(c) = -QLA_real(a) - QLA_real(b);\
                                    QLA_imag(c) = -QLA_imag(a) - QLA_imag(b);}

#define QLA_c_peq_c_plus_c(c,a,b)  {QLA_real(c) += QLA_real(a) + QLA_real(b);\
                                    QLA_imag(c) += QLA_imag(a) + QLA_imag(b);}

#define QLA_c_meq_c_plus_c(c,a,b)  {QLA_real(c) -= QLA_real(a) + QLA_real(b);\
                                    QLA_imag(c) -= QLA_imag(a) + QLA_imag(b);}

#define QLA_c_eq_c_plus_ic(c,a,b) {QLA_real(c) = QLA_real(a) - QLA_imag(b);\
                                   QLA_imag(c) = QLA_imag(a) + QLA_real(b);}

#define QLA_c_eqm_c_plus_ic(c,a,b) {QLA_real(c) = -QLA_real(a) + QLA_imag(b);\
                                    QLA_imag(c) = -QLA_imag(a) - QLA_real(b);}

#define QLA_c_peq_c_plus_ic(c,a,b) {QLA_real(c) += QLA_real(a) - QLA_imag(b);\
                                    QLA_imag(c) += QLA_imag(a) + QLA_real(b);}

#define QLA_c_meq_c_plus_ic(c,a,b) {QLA_real(c) -= QLA_real(a) - QLA_imag(b);\
                                    QLA_imag(c) -= QLA_imag(a) + QLA_real(b);}

#define QLA_c_eq_r_plus_ir(c,a,b) {QLA_real(c) = (a); QLA_imag(c) = (b);}


#define QLA_c_peq_r_plus_ir(c,a,b) {QLA_real(c) += (a); QLA_imag(c) += (b);}


#define QLA_c_eqm_r_plus_ir(c,a,b) {QLA_real(c) = -(a); QLA_imag(c) = -(b);}


#define QLA_c_meq_r_plus_ir(c,a,b) {QLA_real(c) -= (a); QLA_imag(c) -= (b);}


#define QLA_c_eq_c_minus_c(c,a,b) {QLA_real(c) = QLA_real(a) - QLA_real(b);\
                                   QLA_imag(c) = QLA_imag(a) - QLA_imag(b);}

#define QLA_c_eqm_c_minus_c(c,a,b) {QLA_real(c) = -QLA_real(a) + QLA_real(b);\
                                    QLA_imag(c) = -QLA_imag(a) + QLA_imag(b);}

#define QLA_c_peq_c_minus_c(c,a,b) {QLA_real(c) += QLA_real(a) - QLA_real(b);\
                                    QLA_imag(c) += QLA_imag(a) - QLA_imag(b);}

#define QLA_c_meq_c_minus_c(c,a,b) {QLA_real(c) -= QLA_real(a) - QLA_real(b);\
                                    QLA_imag(c) -= QLA_imag(a) - QLA_imag(b);}

#define QLA_c_eq_c_minus_ca(c,a,b) {QLA_real(c) = QLA_real(a) - QLA_real(b);\
                                   QLA_imag(c) = QLA_imag(a) + QLA_imag(b);}

#define QLA_c_eq_c_minus_ic(c,a,b){QLA_real(c) = QLA_real(a) + QLA_imag(b);\
                                   QLA_imag(c) = QLA_imag(a) - QLA_real(b);}

#define QLA_c_eqm_c_minus_ic(c,a,b){QLA_real(c) = -QLA_real(a) - QLA_imag(b);\
                                    QLA_imag(c) = -QLA_imag(a) + QLA_real(b);}

#define QLA_c_peq_c_minus_ic(c,a,b){QLA_real(c) += QLA_real(a) + QLA_imag(b);\
                                    QLA_imag(c) += QLA_imag(a) - QLA_real(b);}

#define QLA_c_meq_c_minus_ic(c,a,b){QLA_real(c) -= QLA_real(a) + QLA_imag(b);\
                                    QLA_imag(c) -= QLA_imag(a) - QLA_real(b);}

#define QLA_c_eq_c_times_c(c,a,b)  {QLA_real(c) = QLA_real(a)*QLA_real(b)\
                                                 - QLA_imag(a)*QLA_imag(b); \
		                    QLA_imag(c) = QLA_real(a)*QLA_imag(b)\
                                                 + QLA_imag(a)*QLA_real(b);}

#define QLA_c_peq_c_times_c(c,a,b) {QLA_real(c) += (QLA_real(a)*QLA_real(b)\
                                                 - QLA_imag(a)*QLA_imag(b)); \
		                    QLA_imag(c) += (QLA_real(a)*QLA_imag(b)\
                                                 + QLA_imag(a)*QLA_real(b));}

#define QLA_c_eqm_c_times_c(c,a,b) {QLA_real(c) = -QLA_real(a)*QLA_real(b)\
                                                 + QLA_imag(a)*QLA_imag(b); \
		                    QLA_imag(c) = -QLA_real(a)*QLA_imag(b)\
                                                 - QLA_imag(a)*QLA_real(b);}

#define QLA_c_meq_c_times_c(c,a,b) {QLA_real(c) -= (QLA_real(a)*QLA_real(b)\
                                                 - QLA_imag(a)*QLA_imag(b)); \
		                    QLA_imag(c) -= (QLA_real(a)*QLA_imag(b)\
                                                 + QLA_imag(a)*QLA_real(b));}

#define QLA_r_eq_Re_c_times_c(c,a,b)  { c = QLA_real(a)*QLA_real(b)\
                                          - QLA_imag(a)*QLA_imag(b);}

#define QLA_r_peq_Re_c_times_c(c,a,b) { c += (QLA_real(a)*QLA_real(b)\
                                          - QLA_imag(a)*QLA_imag(b));}

#define QLA_r_eqm_Re_c_times_c(c,a,b) { c = -QLA_real(a)*QLA_real(b)\
                                          + QLA_imag(a)*QLA_imag(b);}

#define QLA_r_meq_Re_c_times_c(c,a,b) { c -= (QLA_real(a)*QLA_real(b)\
                                          - QLA_imag(a)*QLA_imag(b));}

#define QLA_r_eq_Im_c_times_c(c,a,b)  { c = QLA_real(a)*QLA_imag(b)\
                                          + QLA_imag(a)*QLA_real(b);}

#define QLA_r_peq_Im_c_times_c(c,a,b) { c += (QLA_real(a)*QLA_imag(b)\
                                          + QLA_imag(a)*QLA_real(b));}

#define QLA_r_eqm_Im_c_times_c(c,a,b) { c = -QLA_real(a)*QLA_imag(b)\
                                          - QLA_imag(a)*QLA_real(b);}

#define QLA_r_meq_Im_c_times_c(c,a,b) { c -= (QLA_real(a)*QLA_imag(b)\
                                          + QLA_imag(a)*QLA_real(b));}

#define QLA_c_eq_c_div_c(c,a,b) { QLA_D_Real t = QLA_real(b)*QLA_real(b)\
                                          + QLA_imag(b)*QLA_imag(b); \
		                  QLA_real(c) = (QLA_real(a)*QLA_real(b)\
                                          + QLA_imag(a)*QLA_imag(b))/t; \
		                  QLA_imag(c) = (QLA_imag(a)*QLA_real(b)\
                                          - QLA_real(a)*QLA_imag(b))/t; }

#define QLA_c_eq_c_times_ca(c,a,b) {QLA_real(c) = QLA_real(a)*QLA_real(b)\
                                                  + QLA_imag(a)*QLA_imag(b); \
	  	                    QLA_imag(c) = QLA_imag(a)*QLA_real(b)\
                                                 - QLA_real(a)*QLA_imag(b);}

#define QLA_c_peq_c_times_ca(c,a,b) {QLA_real(c) += QLA_real(a)*QLA_real(b)\
                                                  + QLA_imag(a)*QLA_imag(b); \
	  	                    QLA_imag(c) += QLA_imag(a)*QLA_real(b)\
                                                 - QLA_real(a)*QLA_imag(b);}

#define QLA_c_eqm_c_times_ca(c,a,b) {QLA_real(c) = -QLA_real(a)*QLA_real(b)\
                                                  - QLA_imag(a)*QLA_imag(b); \
	  	                    QLA_imag(c) = -QLA_imag(a)*QLA_real(b)\
                                                 + QLA_real(a)*QLA_imag(b);}

#define QLA_c_meq_c_times_ca(c,a,b) {QLA_real(c) -= QLA_real(a)*QLA_real(b)\
                                                  + QLA_imag(a)*QLA_imag(b); \
	  	                    QLA_imag(c) -= QLA_imag(a)*QLA_real(b)\
                                                 - QLA_real(a)*QLA_imag(b);}

#define QLA_r_eq_Re_c_times_ca(c,a,b) {c = QLA_real(a)*QLA_real(b)\
                                         + QLA_imag(a)*QLA_imag(b);}

#define QLA_r_peq_Re_c_times_ca(c,a,b) {c += QLA_real(a)*QLA_real(b)\
                                           + QLA_imag(a)*QLA_imag(b);}

#define QLA_r_eqm_Re_c_times_ca(c,a,b) {c = -QLA_real(a)*QLA_real(b)\
                                          - QLA_imag(a)*QLA_imag(b);}

#define QLA_r_meq_Re_c_times_ca(c,a,b) {c -= QLA_real(a)*QLA_real(b)\
                                           + QLA_imag(a)*QLA_imag(b);}

#define QLA_r_eq_Im_c_times_ca(c,a,b) {c = QLA_imag(a)*QLA_real(b)\
                                         - QLA_real(a)*QLA_imag(b);}

#define QLA_r_peq_Im_c_times_ca(c,a,b) {c += QLA_imag(a)*QLA_real(b)\
                                           - QLA_real(a)*QLA_imag(b);}

#define QLA_r_eqm_Im_c_times_ca(c,a,b) {c = -QLA_imag(a)*QLA_real(b)\
                                          + QLA_real(a)*QLA_imag(b);}

#define QLA_r_meq_Im_c_times_ca(c,a,b) {c -= QLA_imag(a)*QLA_real(b)\
                                           - QLA_real(a)*QLA_imag(b);}

#define QLA_c_eq_ca_times_c(c,a,b) {QLA_real(c) = QLA_real(a)*QLA_real(b)\
                                                 + QLA_imag(a)*QLA_imag(b); \
		                    QLA_imag(c) = QLA_real(a)*QLA_imag(b)\
                                                  - QLA_imag(a)*QLA_real(b);}

#define QLA_c_peq_ca_times_c(c,a,b) {QLA_real(c) += QLA_real(a)*QLA_real(b)\
                                                 + QLA_imag(a)*QLA_imag(b); \
		                    QLA_imag(c) += QLA_real(a)*QLA_imag(b)\
                                                  - QLA_imag(a)*QLA_real(b);}

#define QLA_c_eqm_ca_times_c(c,a,b) {QLA_real(c) = -QLA_real(a)*QLA_real(b)\
                                                 - QLA_imag(a)*QLA_imag(b); \
		                    QLA_imag(c) = -QLA_real(a)*QLA_imag(b)\
                                                  + QLA_imag(a)*QLA_real(b);}

#define QLA_c_meq_ca_times_c(c,a,b) {QLA_real(c) -= QLA_real(a)*QLA_real(b)\
                                                 + QLA_imag(a)*QLA_imag(b); \
		                    QLA_imag(c) -= QLA_real(a)*QLA_imag(b)\
                                                  - QLA_imag(a)*QLA_real(b);}

#define QLA_r_eq_Re_ca_times_c(c,a,b) { c = QLA_real(a)*QLA_real(b)\
                                          + QLA_imag(a)*QLA_imag(b);}

#define QLA_r_peq_Re_ca_times_c(c,a,b) { c += QLA_real(a)*QLA_real(b)\
                                            + QLA_imag(a)*QLA_imag(b);}

#define QLA_r_eqm_Re_ca_times_c(c,a,b) {c = -QLA_real(a)*QLA_real(b)\
                                           - QLA_imag(a)*QLA_imag(b);}

#define QLA_r_meq_Re_ca_times_c(c,a,b) {c -= QLA_real(a)*QLA_real(b)\
                                           + QLA_imag(a)*QLA_imag(b);}

#define QLA_r_eq_Im_ca_times_c(c,a,b) {c = QLA_real(a)*QLA_imag(b)\
                                         - QLA_imag(a)*QLA_real(b);}

#define QLA_r_peq_Im_ca_times_c(c,a,b) {c += QLA_real(a)*QLA_imag(b)\
                                          - QLA_imag(a)*QLA_real(b);}

#define QLA_r_eqm_Im_ca_times_c(c,a,b) {c = -QLA_real(a)*QLA_imag(b)\
                                          + QLA_imag(a)*QLA_real(b);}

#define QLA_r_meq_Im_ca_times_c(c,a,b) {c -= QLA_real(a)*QLA_imag(b)\
                                           - QLA_imag(a)*QLA_real(b);}


#define QLA_c_eq_ca_times_ca(c,a,b) {QLA_real(c) = QLA_real(a)*QLA_real(b)\
                                                  - QLA_imag(a)*QLA_imag(b); \
		                     QLA_imag(c) = -QLA_real(a)*QLA_imag(b)\
                                                   - QLA_imag(a)*QLA_real(b);}

#define QLA_c_peq_ca_times_ca(c,a,b) {QLA_real(c) += QLA_real(a)*QLA_real(b)\
                                                  - QLA_imag(a)*QLA_imag(b); \
		                     QLA_imag(c) += -QLA_real(a)*QLA_imag(b)\
                                                   - QLA_imag(a)*QLA_real(b);}

#define QLA_c_eqm_ca_times_ca(c,a,b) {QLA_real(c) = -QLA_real(a)*QLA_real(b)\
                                                  + QLA_imag(a)*QLA_imag(b); \
		                     QLA_imag(c) = QLA_real(a)*QLA_imag(b)\
                                                   + QLA_imag(a)*QLA_real(b);}

#define QLA_c_meq_ca_times_ca(c,a,b) {QLA_real(c) -= QLA_real(a)*QLA_real(b)\
                                                  - QLA_imag(a)*QLA_imag(b); \
		                     QLA_imag(c) -= -QLA_real(a)*QLA_imag(b)\
                                                   - QLA_imag(a)*QLA_real(b);}

#define QLA_r_eq_Re_ca_times_ca(c,a,b) {c = QLA_real(a)*QLA_real(b)\
                                          - QLA_imag(a)*QLA_imag(b);}

#define QLA_r_peq_Re_ca_times_ca(c,a,b) {c += QLA_real(a)*QLA_real(b)\
                                            - QLA_imag(a)*QLA_imag(b);}

#define QLA_r_eqm_Re_ca_times_ca(c,a,b) {c = -QLA_real(a)*QLA_real(b)\
                                           + QLA_imag(a)*QLA_imag(b);}

#define QLA_r_meq_Re_ca_times_ca(c,a,b) {c -= QLA_real(a)*QLA_real(b)\
                                            - QLA_imag(a)*QLA_imag(b);}

#define QLA_r_eq_Im_ca_times_ca(c,a,b) {c = -QLA_real(a)*QLA_imag(b)\
                                          - QLA_imag(a)*QLA_real(b);}

#define QLA_r_peq_Im_ca_times_ca(c,a,b) {c += -QLA_real(a)*QLA_imag(b)\
                                            - QLA_imag(a)*QLA_real(b);}

#define QLA_r_eqm_Im_ca_times_ca(c,a,b) {c = QLA_real(a)*QLA_imag(b)\
                                           + QLA_imag(a)*QLA_real(b);}

#define QLA_r_meq_Im_ca_times_ca(c,a,b) {c -= -QLA_real(a)*QLA_imag(b)\
                                            - QLA_imag(a)*QLA_real(b);}

#define QLA_c_eq_c_times_r(c,a,b)   {QLA_real(c) = (b) * QLA_real(a);\
                                     QLA_imag(c) = (b) * QLA_imag(a); }

#define QLA_c_peq_c_times_r(c,a,b)   {QLA_real(c) += (b) * QLA_real(a);\
                                     QLA_imag(c) += (b) * QLA_imag(a); }

#define QLA_c_eqm_c_times_r(c,a,b)   {QLA_real(c) = -((b) * QLA_real(a));\
                                     QLA_imag(c) = -((b) * QLA_imag(a)); }

#define QLA_c_meq_c_times_r(c,a,b)   {QLA_real(c) -= (b) * QLA_real(a);\
                                     QLA_imag(c) -= (b) * QLA_imag(a); }

#define QLA_c_eq_r_times_c(c,a,b)    QLA_c_eq_c_times_r(c,b,a)
#define QLA_c_peq_r_times_c(c,a,b)   QLA_c_peq_c_times_r(c,b,a)
#define QLA_c_eqm_r_times_c(c,a,b)   QLA_c_eqm_c_times_r(c,b,a)
#define QLA_c_meq_r_times_c(c,a,b)   QLA_c_meq_c_times_r(c,b,a)

#define QLA_c_peq_c_times_r(c,a,b)  {QLA_real(c) += (b) * QLA_real(a);\
                                     QLA_imag(c) += (b) * QLA_imag(a); }

#define QLA_c_peq_r_times_c(c,a,b)   QLA_c_peq_c_times_r(c,b,a)

#define QLA_c_eq_c_div_r(c,a,b)     {QLA_real(c) = QLA_real(a)/(b);\
                                     QLA_imag(c) = QLA_imag(a)/(b);}


/* Ternary operations */
#define QLA_c_eq_c_times_c_plus_c(c,a,x,b) \
                      {QLA_real(c) = QLA_real(a)*QLA_real(x)\
                                   - QLA_imag(a)*QLA_imag(x) + QLA_real(b); \
		       QLA_imag(c) = QLA_real(a)*QLA_imag(x)\
                                   + QLA_imag(a)*QLA_real(x) + QLA_imag(b);}
                                 

#define QLA_c_eq_c_times_c_minus_c(c,a,x,b) \
                      {QLA_real(c) = QLA_real(a)*QLA_real(x)\
                                   - QLA_imag(a)*QLA_imag(x) - QLA_real(b); \
		       QLA_imag(c) = QLA_real(a)*QLA_imag(x)\
                                   + QLA_imag(a)*QLA_real(x) - QLA_imag(b);}
                                 

#define QLA_c_eq_c_times_r_plus_r(c,a,x,b) \
                      {QLA_real(c) = QLA_real(a)*(x) + (b); \
		       QLA_imag(c) = QLA_imag(a)*(x);}
                                 

#define QLA_c_eq_c_times_r_minus_r(c,a,x,b) \
                      {QLA_real(c) = QLA_real(a)*(x) - (b); \
		       QLA_imag(c) = QLA_imag(a)*(x);}

                                 

#define QLA_c_eq_r_times_c_plus_c(c,a,x,b) \
                      {QLA_real(c) = (a)*QLA_real(x) + QLA_real(b); \
		       QLA_imag(c) = (a)*QLA_imag(x) + QLA_imag(b);}
                                 

#define QLA_c_eq_r_times_c_minus_c(c,a,x,b) \
                      {QLA_real(c) = (a)*QLA_real(x) - QLA_real(b); \
		       QLA_imag(c) = (a)*QLA_imag(x) - QLA_imag(b);}
                                 

#endif	/* _QLA_COMPLEX_H */
