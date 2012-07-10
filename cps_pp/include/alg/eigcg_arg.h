/*
 * Please do not edit this file.
 * It was generated using PAB's VML system.
 */

#ifndef _EIGCG_ARG_H_RPCGEN
#define _EIGCG_ARG_H_RPCGEN

#include <config.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
CPS_START_NAMESPACE

class VML;
class EigCGArg {
public:
	 bool Encode(char *filename,char *instance);
	 bool Decode(char *filename,char *instance);
	 bool Vml(VML *vmls,char *instance);
	int nev;
	int m;
	int max_def_len;
	Float max_eig_cut;
	bool_t always_restart;
	int restart_len;
	Float restart[10];
	   EigCGArg (  ) ;
	   ~EigCGArg (  ) ;
};

/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_EigCGArg (VML *, char *instance, EigCGArg*);

#else /* K&R C */
extern  bool_t vml_EigCGArg (VML *, char *instance, EigCGArg*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !_EIGCG_ARG_H_RPCGEN */