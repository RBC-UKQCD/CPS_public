#include <config.h>
#ifdef USE_QIO
#define VOLFMT QIO_SINGLEFILE

#include <alg/qio_prop.h>

#if TARGET == QCDOC
#include <util/qio_readPropagator.h>
#include <util/qio_writePropagator.h>
#endif

CPS_START_NAMESPACE

QIO_Prop::QIO_Prop()
{
  cname = "QIO_Prop";
  char *fname = "QIO_Prop()";
    
  VRB.Func(cname,fname);
}

QIO_Prop::~QIO_Prop()
{
}

void QIO_Prop::QIO_SaveProp(QPropW& prop, char* qio_out_prop, char* ensemble_label, int ensemble_num )
{
#if TARGET == QCDOC
  const FP_FORMAT outformat (   FP_IEEE64BIG);

                             /* one of
                                FP_UNKNOWN
                                FP_AUTOMATIC
                                FP_TIDSP32
                                FP_IEEE32
                                FP_IEEE32BIG
                                FP_IEEE32LITTLE
                                FP_IEEE64
                                FP_IEEE64BIG
                                FP_IEEE64LITTLE
                             */
#endif
  printf("   ---- starting QIO-part ----- \n");
	
  {
    char ensemble_id[200];
    char* dummy;
    int idummy;
    sprintf(ensemble_id, "%d", ensemble_num);

#if TARGET == QCDOC
    qio_writePropagator writePropQio(idummy, &dummy);

    writePropQio.setHeader(ensemble_id, ensemble_label, ensemble_num);
#endif
    printf("  writing: %s (QIO-format, double)\n",qio_out_prop);
#if TARGET == QCDOC
//    writePropQio.write(qio_out_prop, &prop[0], 12, VOLFMT);
#endif
  }

}

void QIO_Prop::QIO_ReadProp(QPropW &prop, char* qio_in_prop)
{
#if TARGET == QCDOC
  const FP_FORMAT outformat (   FP_IEEE64BIG);

                             /* one of
                                FP_UNKNOWN
                                FP_AUTOMATIC
                                FP_TIDSP32
                                FP_IEEE32
                                FP_IEEE32BIG
                                FP_IEEE32LITTLE
                                FP_IEEE64
                                FP_IEEE64BIG
                                FP_IEEE64LITTLE
                             */
#endif
  printf("   ---- starting QIO-part ----- \n");
	
  {
    char* dummy;
    int idummy;

#if TARGET == QCDOC
    qio_readPropagator readPropQio(idummy, &dummy);
#endif
    printf("   reading: %s (QIO-format, double)\n",qio_in_prop);
#if TARGET == QCDOC
 //   readPropQio.read(qio_in_prop, &prop[0], 12, VOLFMT);
#endif
  }

}

CPS_END_NAMESPACE
#endif
