#include<config.h>
CPS_START_NAMESPACE
/*  inst_arg.h */

/*  The structure type InstArg holds the parameters relevant
    to the generation of an instanton configuration. */

#ifndef INCLUDED_INST_ARG_H
#define INCLUDED_INST_ARG_H

CPS_END_NAMESPACE
#include <util/data_types.h>
CPS_START_NAMESPACE



enum InstType { SINGULAR			= 0,
		REGULAR			        = 1,
		REGULAR_TRANSFORMED		= 2,
		REGULAR_SQUASHED		= 3,
		REGULAR_SQUASHED_TRANSFORMED	= 4,
		CONSTANT_FIELD		        = 5 };

enum InstMethod { ADD				= 0,
		DESTROY				= 1 };

struct InstArg {

  InstType inst_kind;  /* The kind of instanton configuration.  */

  InstMethod inst_method;
		       /* Transcription paradigm:
			ADD = multiply link by link,
			DESTROY = replace exsiting links	*/

  Float charge;        /* The instanton charge.                 */

  Float rho;           /* The instanton radius.                 */

  Float rho_cutoff;    /* The cutoff radius.                    */

  int n1;	       /* Abelian charge index 1 for InstType=5 */

  int n2;	       /* Abelian charge index 2 for InstType=5 */
};

#endif /* !INCLUDED_INST_ARG_H */

CPS_END_NAMESPACE
