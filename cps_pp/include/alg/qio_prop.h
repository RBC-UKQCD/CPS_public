//------------------------------------------------------------------
//
// qio_prop.h
//
//------------------------------------------------------------------


#ifndef INCLUDED_QIO_PROP_H
#define INCLUDED_QIO_PROP_H

#include <config.h>

#include <stdlib.h>	// exit()
#include <stdio.h>

#include <util/lattice.h>
#include <alg/alg_base.h>

//#include <alg/qpropw.h>
#include "qpropw.h"

#if TARGET == QCDOC
#include <comms/sysfunc_cps.h>
#endif

CPS_START_NAMESPACE


class QIO_Prop
{
 private:
    char* cname;
    
 public:
    QIO_Prop();

    ~QIO_Prop();

    void QIO_SaveProp(QPropW& prop, char*, char*, int) ;
    void QIO_ReadProp(QPropW& prop, char*) ;

};

CPS_END_NAMESPACE

#endif
