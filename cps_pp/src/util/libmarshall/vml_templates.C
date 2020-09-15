#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

CPS_START_NAMESPACE
#include <util/vml/vml_templates.h>

#ifdef _USE_STDLIB

void rpc_print<char *>::doit(char const * const &what, const int &size, const std::string &prefix){
  std::cout << prefix << "'" << what << "' (" << static_cast<void const * const &>(what) << ")\n";
}

#endif

CPS_END_NAMESPACE
