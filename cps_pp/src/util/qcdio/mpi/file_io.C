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

int Fprintf( FileIoType type, FILE *stream, const char *format,...){

    va_list args;
    va_start(args,format);
    int nbytes = vfprintf(stream,format,args);
    va_end(args);
    return nbytes;
}

int Fprintf( FILE *stream, const char *format,...){

    va_list args;
    va_start(args,format);
    int nbytes = vfprintf(stream,format,args);
    va_end(args);
    return nbytes;
}

int Vfprintf( FileIoType type, FILE *stream, const char *format, va_list ap){
    return vfprintf(stream, format, ap);
}
CPS_END_NAMESPACE
