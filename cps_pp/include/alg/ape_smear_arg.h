/*
 * Please do not edit this file.
 * It was generated by CPC make system.
 */

#ifndef _APE_SMEAR_ARG_H_RPCGEN
#define _APE_SMEAR_ARG_H_RPCGEN

#include <config.h>
#include <util/vml/types.h>
#include <util/vml/vml.h>
#include <util/enum.h>
#include <util/defines.h>
CPS_START_NAMESPACE

class ApeSmearArg {
public:
	 ApeSmearArg(char *filename);
	 void Encode(char *filename,char *instance);
	 void Decode(char *filename,char *instance);
	 void Vml(VML *vmls,char *instance);
	Float tolerance;
	int orthog;
	Float coef;
	   ApeSmearArg (  ) ;
};

/* the xdr functions */

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)
extern  bool_t vml_ApeSmearArg (VML *, char *instance, ApeSmearArg*);

#else /* K&R C */
extern  bool_t vml_ApeSmearArg (VML *, char *instance, ApeSmearArg*);

#endif /* K&R C */

#ifdef __cplusplus
}
#endif
CPS_END_NAMESPACE

#endif /* !_APE_SMEAR_ARG_H_RPCGEN */