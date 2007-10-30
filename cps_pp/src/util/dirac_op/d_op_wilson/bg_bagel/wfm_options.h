#ifndef WFM_OPTIONS_H
#define WFM_OPTIONS_H

/* Undef the unwanted from the environment -- eg the compiler command line */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#include "wfm_options_internal.h"

/* Prefix everything with WFM_ */
static const char* const WFM_PACKAGE(PACKAGE);
static const char* const WFM_PACKAGE_BUGREPORT(PACKAGE_BUGREPORT);
static const char* const WFM_PACKAGE_NAME(PACKAGE_NAME);
static const char* const WFM_PACKAGE_STRING(PACKAGE_STRING);
static const char* const WFM_PACKAGE_TARNAME(PACKAGE_TARNAME);
static const char* const WFM_PACKAGE_VERSION(PACKAGE_VERSION);
static const char* const WFM_VERSION(VERSION);
                                                                                
/* Undef the unwanted */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#endif

