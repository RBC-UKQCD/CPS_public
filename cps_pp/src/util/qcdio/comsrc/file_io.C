#include <config.h>
#include <stdio.h>
#include <stdarg.h>
#include <util/qcdio.h>
//#include <qalloc.h>
CPS_START_NAMESPACE
static FILE FAKE;
static FILE *FAKE_P = &FAKE;
const int MAX_FILENAME =200;
FILE *Fopen( FileIoType type, const char *filename, const char *mode){
  FILE *result = NULL;
#ifdef UNIFORM_SEED_NO_COMMS
  return FAKE_P;
#endif
  if ( type == ZERO_ONLY && UniqueID() ) return FAKE_P;
  if(type == ADD_ID){
    char fname[MAX_FILENAME];
    if(strlen(filename)+6 >MAX_FILENAME){
	  fprintf(stderr,"Fopen: filename(%s) is too long\n",filename);
      result = NULL;
    }
    else {
      sprintf(fname,"%s.%d",filename,UniqueID());
      result = fopen(fname,mode);
    }
  } else {
    result = fopen(filename,mode);
  }
  if (result == NULL){
	fprintf(stderr,"Fopen: opening %s failed\n",filename);
  }
  return result;
}

size_t Fwrite( const void *ptr, size_t size, size_t n, FILE *stream){
  if ( stream == FAKE_P )  return n;
  return fwrite( ptr, size, n, stream);
}

size_t Fread(void *ptr, size_t size, size_t n, FILE *stream){
  if ( stream == FAKE_P )
    ERR.General("","Fread()","Trying to read from invalid file stream, possibly from using Fopen with wrong flags");
  return fread( ptr, size, n, stream);
}

int Fflush(FILE *stream)
{
    if(stream == FAKE_P) return 0;
    return fflush(stream);
}

int Fclose( FileIoType type, FILE *stream){
  if ( stream == FAKE_P )  return 1;
  return fclose(stream);
}

int Fprintf( FileIoType type, FILE *stream, const char *format,...){
  if ( stream == FAKE_P )  return 1;
  va_list args;
  va_start(args,format);
  int nb = vfprintf(stream, format, args);
  va_end(args);
  return nb;
}

int Fprintf( FILE *stream, const char *format,...){
  if ( stream == FAKE_P )  return 1;
  va_list args;
  va_start(args,format);
  int nb = vfprintf(stream, format, args);
  va_end(args);
  return nb;
}

int Vfprintf( FileIoType type, FILE *stream, const char *format, va_list ap){
  if ( stream == FAKE_P )  return 1;
  return vfprintf(stream, format, ap);
}
CPS_END_NAMESPACE
