/* Hacked by Peter Boyle for VML 2004 *//*
 * rpc_pktout.c, UDP packet class outputter for the RPC protocol compiler
 * PAB 2002.
 * This auto-derives classes from the base class RPCPkt which I have written.
 * For use in the Qdaemon with RPCTransporter class.
 */
#include <stdio.h>
#include <ctype.h>
#include "rpc_parse.h"
#include "rpc_util.h"
#include "proto.h"

void print_include(void)
{

}

void print_pkt_derivedclass( definition *def )
{
  version_list *vers;
  proc_list *proc;
  int i;
  const char *ext;
  int vno ;
  decl_list *dl;

  if ( def->def_kind == DEF_PROGRAM ) {

    for ( vers = def->def.pr.versions; vers != NULL; vers = vers->next){

      for (proc = vers->procs; proc != NULL;proc = proc->next){

/* Argument packet */
      dl = proc->args.decls;

      fprintf(stderr,"QCDOC-RPC derived class: ");
      fprintf(stderr,"RPCPkt_%s\n",proc->proc_name);


      f_print (fout, "\n");
      f_print (fout, "class RPCPkt_%s : public RPCPkt {",proc->proc_name);

      /*Input and output arg/result*/
      f_print (fout, "\nprivate:\n");
      f_print (fout, "\t");
      ptype(dl->decl.prefix,dl->decl.type, 1);
      f_print (fout, "in;\n");
      f_print (fout, "\t");
      ptype(proc->res_prefix,proc->res_type, 1);
      f_print (fout,"out; \n");
      f_print (fout,"\tvirtual caddr_t Request(void) { return (caddr_t)&in; } ;\n");
      f_print (fout,"\tvirtual caddr_t Reply  (void) { return (caddr_t)&out; } ;\n\n");

      f_print (fout, "\npublic:\n");

      /*Constructor*/
      f_print (fout, "   RPCPkt_%s( %s *in) :\n",proc->proc_name,dl->decl.type);
      f_print (fout,"\n");
      f_print (fout, "\tRPCPkt( (xdrproc_t)xdr_");
      ptype(dl->decl.prefix,dl->decl.type, 1);
      f_print (fout,",\n\t\t(xdrproc_t)xdr_");
      ptype(proc->res_prefix,proc->res_type, 1);
      f_print (fout,",\n\t\t%s,%s,%s",proc->proc_name,vers->vers_name,def->def_name);
      f_print (fout,") {\n");
      f_print (fout,"       xdrencode((caddr_t)in);\n");
      f_print (fout,"    };\n");

      f_print (fout,"};\n");

      }

    }

   }
}






