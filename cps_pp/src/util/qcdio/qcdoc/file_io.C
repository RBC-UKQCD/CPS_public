#include <config.h>
#include <stdio.h>
#include <stdarg.h>
#include <util/qcdio.h>
#include <qalloc.h>
CPS_START_NAMESPACE
static FILE FAKE;
const int MAX_FILENAME =200;
FILE *Fopen( FileIoType type, const char *filename, const char *mode){
  if ( type == ZERO_ONLY && UniqueID() ) return &FAKE;
  if(type == ADD_ID){
    char fname[MAX_FILENAME];
    if(strlen(filename)+6 >MAX_FILENAME){
	  fprintf(stderr,"Fopen: filename(%s) is too long\n",filename);
      return NULL;
    }
    sprintf(fname,"%s.%d",filename,UniqueID());
    return fopen(fname,mode);
  } else {
    return fopen(filename,mode);
  }
}

int Fclose( FileIoType type, FILE *stream){
//  if ( type == ZERO_ONLY && UniqueID() ) return 1;
  if ( stream == &FAKE )  return 1;
  return fclose(stream);
}

int Fprintf( FileIoType type, FILE *stream, const char *format,...){
//  if ( type == ZERO_ONLY && UniqueID() ) return 1;
  if ( stream == &FAKE )  return 1;
#if 0
  {
     ERR.General("","Fprintf()","Trying to write to a fake file pointer, (did you open without ADD_ID option?)\n");
  }
#endif
  va_list args;
  va_start(args,format);
  return vfprintf(stream,format,args);
}

int Fprintf( FILE *stream, const char *format,...){
//  if ( UniqueID() ) return 1;
  if ( stream == &FAKE )  return 1;
  va_list args;
  va_start(args,format);
  return vfprintf(stream,format,args);
}

int Vfprintf( FileIoType type, FILE *stream, const char *format, va_list ap){
//  if ( type == ZERO_ONLY && UniqueID() ) return 1;
  if ( stream == &FAKE )  return 1;
  return vfprintf(stream, format, ap);
}
CPS_END_NAMESPACE
