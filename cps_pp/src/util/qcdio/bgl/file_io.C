#include <config.h>
#include <stdio.h>
#include <stdarg.h>
#include <util/qcdio.h>
CPS_START_NAMESPACE
static FILE FAKE;
const int MAX_FILENAME =200;
FILE *Fopen( FileIoType type, const char *filename, const char *mode){
  FILE *fp = NULL;
  if ( type == ZERO_ONLY && UniqueID() ) return &FAKE;
  if(type == ADD_ID){
    char fname[MAX_FILENAME];
    if(strlen(filename)+6 >MAX_FILENAME){
	  fprintf(stderr,"Fopen: filename(%s) is too long\n",filename);
      return NULL;
    }
    sprintf(fname,"%s.%d",filename,UniqueID());
    fp = fopen(fname,mode);
  } else {
    fp =  fopen(filename,mode);
  }
  if (fp==NULL)
    ERR.General("","Fopen","cannot open %s",filename);
  return fp;
}

int Fclose( FileIoType type, FILE *stream){
  if ( stream == &FAKE )  return 1;
  return fclose(stream);
}

int Fprintf( FileIoType type, FILE *stream, const char *format,...){
  if ( stream == &FAKE )  return 1;
  va_list args;
  va_start(args,format);
  int nb = vfprintf(stream, format, args);
  va_end(args);
  return nb;
}

int Fprintf( FILE *stream, const char *format,...){
  if ( stream == &FAKE )  return 1;
  va_list args;
  va_start(args,format);
  int nb = vfprintf(stream, format, args);
  va_end(args);
  return nb;
}

int Vfprintf( FileIoType type, FILE *stream, const char *format, va_list ap){
  if ( stream == &FAKE )  return 1;
  return vfprintf(stream, format, ap);
}
CPS_END_NAMESPACE
