#include <config.h>
#include <stdio.h>
#include <stdarg.h>
#include <util/qcdio.h>
CPS_START_NAMESPACE
FILE *Fopen( FileIoType type, const char *filename, const char *mode){
return fopen(filename,mode);
}
int Fclose( FileIoType type, FILE *stream){
return fclose(stream);
}

int Fprintf(FILE *stream, const char *format,...){
  va_list args;
  va_start(args,format);
  vfprintf(stream,format,args);
}

int Fprintf( FileIoType type, FILE *stream, const char *format,...){
  va_list args;
  va_start(args,format);
  vfprintf(stream,format,args);
}

int Vfprintf( FileIoType type, FILE *stream, const char *format, va_list ap){
  return vfprintf(stream, format, ap);
}
CPS_END_NAMESPACE
