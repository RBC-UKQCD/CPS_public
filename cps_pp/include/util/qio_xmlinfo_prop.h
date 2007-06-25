#ifndef __QIOXMLINFOPROP__
#define __QIOXMLINFOPROP__

#include <util/qio_general.h>
#include <util/fpconv.h>


CPS_START_NAMESPACE
using namespace std;


// these are missing in qio_xml.h
//
//char *QIO_encode_as_string(char *buf, QIO_TagCharValue *tag_value,
//			   int *remainder);
//
//char *QIO_next_tag(char *parse_pt, char *tag, char **left_angle);
//
//char *QIO_get_tag_value(char *parse_pt, char *tag, char *value_string);
//
//void QIO_decode_as_string(char *tag, char *value_string,
//			  QIO_TagCharValue *tag_value);
//
//int QIO_check_string_occur(QIO_TagCharValue *tag_value);
//
// end missing


typedef struct {
  QIO_TagCharValue version ;
  QIO_TagCharValue type;
  QIO_TagCharValue info;
} CPS_QIO_PROP_FileRecordInfo;


typedef struct {
  QIO_TagCharValue version ;
  QIO_TagIntValue spin;
  QIO_TagIntValue color;
  QIO_TagCharValue info;
} CPS_QIO_PROP_UserRecordInfo;

typedef struct {
  QIO_TagCharValue version ;
  QIO_TagCharValue info;
} CPS_QIO_PROP_PAIRS_UserRecordInfo;


typedef struct {
  QIO_TagCharValue version;
  QIO_TagCharValue info;
} CPS_QIO_SOURCE_UserRecordInfo;

typedef struct {
  QIO_TagCharValue version;
  QIO_TagIntValue spin;
  QIO_TagIntValue color;
  QIO_TagCharValue info;
} CPS_QIO_SOURCE_PAIRS_UserRecordInfo;


#define CPS_QIO_PROP_FILERECORDFORMATVERSION "1.0"

#define CPS_QIO_PROP_USERRECORDFORMATVERSION "1.0"
#define CPS_QIO_PROP_PAIRS_USERRECORDFORMATVERSION "1.0"

#define CPS_QIO_SOURCE_USERRECORDFORMATVERSION "1.0"
#define CPS_QIO_SOURCE_PAIRS_USERRECORDFORMATVERSION "1.0"

#define CPS_QIO_PROP_FILE_RECORD_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"type"   , "", "", 0}, \
   {"info"   , "", "", 0}  \
}

#define CPS_QIO_PROP_USER_RECORD_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"spin"   , "", 0, 0}, \
   {"color" , "", 0, 0}, \
   {"info"   , "", "", 0}  \
}

#define CPS_QIO_PROP_PAIRS_USER_RECORD_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"info"   , "", "", 0}  \
}

#define CPS_QIO_SOURCE_USER_RECORD_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"info"   , "", "", 0}  \
}

#define CPS_QIO_SOURCE_PAIRS_USER_RECORD_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"spin"   , "", 0, 0}, \
   {"color" , "", 0, 0}, \
   {"info"   , "", "", 0}  \
}


#define CPS_QIO_PROP_FILE_RECORD_INFO_WRAPPER {\
  {"usqcdPropFile", "", "" , 0}       \
}

#define CPS_QIO_PROP_USER_RECORD_INFO_WRAPPER {\
  {"usqcdPropInfo", "", "" , 0}       \
}

#define CPS_QIO_PROP_PAIRS_USER_RECORD_INFO_WRAPPER {\
  {"usqcdPropInfo", "", "" , 0}       \
}

#define CPS_QIO_SOURCE_USER_RECORD_INFO_WRAPPER {\
  {"usqcdSourceInfo", "", "" , 0}       \
}

#define CPS_QIO_SOURCE_PAIRS_USER_RECORD_INFO_WRAPPER {\
  {"usqcdSourceInfo", "", "" , 0}       \
}





typedef struct {
  QIO_TagCharValue userrecordinfo_tags;
} CPS_QIO_PROP_FileRecordInfoWrapper;

typedef struct {
  QIO_TagCharValue userrecordinfo_tags;
} CPS_QIO_PROP_UserRecordInfoWrapper;

typedef struct {
  QIO_TagCharValue userrecordinfo_tags;
} CPS_QIO_PROP_PAIRS_UserRecordInfoWrapper;

typedef struct {
  QIO_TagCharValue userrecordinfo_tags;
} CPS_QIO_SOURCE_UserRecordInfoWrapper;

typedef struct {
  QIO_TagCharValue userrecordinfo_tags;
} CPS_QIO_SOURCE_PAIRS_UserRecordInfoWrapper;



CPS_QIO_PROP_FileRecordInfo *CPS_QIO_PROP_create_file_record_info(char *type, char *info);
void CPS_QIO_PROP_destroy_file_record_info(CPS_QIO_PROP_FileRecordInfo *record_info);

CPS_QIO_PROP_UserRecordInfo *CPS_QIO_PROP_create_user_record_info(int spin, int color, char *info);
void CPS_QIO_PROP_destroy_user_record_info(CPS_QIO_PROP_UserRecordInfo *record_info);

CPS_QIO_PROP_PAIRS_UserRecordInfo *CPS_QIO_PROP_PAIRS_create_user_record_info(char *info);
void CPS_QIO_PROP_PAIRS_destroy_user_record_info(CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info);

CPS_QIO_SOURCE_UserRecordInfo *CPS_QIO_SOURCE_create_user_record_info(char *info);
void CPS_QIO_SOURCE_destroy_user_record_info(CPS_QIO_SOURCE_UserRecordInfo *record_info);

CPS_QIO_SOURCE_PAIRS_UserRecordInfo *CPS_QIO_SOURCE_PAIRS_create_user_record_info(int spin, int color, char *info);
void CPS_QIO_SOURCE_PAIRS_destroy_user_record_info(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info);




void CPS_QIO_PROP_encode_file_record_info(QIO_String *record_string, 
				     CPS_QIO_PROP_FileRecordInfo *record_info);


void CPS_QIO_PROP_encode_user_record_info(QIO_String *record_string, 
				     CPS_QIO_PROP_UserRecordInfo *record_info);

void CPS_QIO_PROP_PAIRS_encode_user_record_info(QIO_String *record_string, 
				     CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info);

void CPS_QIO_SOURCE_encode_user_record_info(QIO_String *record_string, 
				     CPS_QIO_SOURCE_UserRecordInfo *record_info);

void CPS_QIO_SOURCE_PAIRS_encode_user_record_info(QIO_String *record_string, 
				     CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info);


int CPS_QIO_PROP_decode_file_record_info(CPS_QIO_PROP_FileRecordInfo *record_info,
				    QIO_String *record_string);

int CPS_QIO_PROP_decode_user_record_info(CPS_QIO_PROP_UserRecordInfo *record_info,
				    QIO_String *record_string);

int CPS_QIO_PROP_PAIRS_decode_user_record_info(CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info,
				    QIO_String *record_string);

int CPS_QIO_SOURCE_decode_user_record_info(CPS_QIO_SOURCE_UserRecordInfo *record_info,
				    QIO_String *record_string);

int CPS_QIO_SOURCE_PAIRS_decode_user_record_info(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info,
				    QIO_String *record_string);



char *CPS_QIO_PROP_file_get_type(CPS_QIO_PROP_FileRecordInfo *record_info);
char *CPS_QIO_PROP_file_get_info(CPS_QIO_PROP_FileRecordInfo *record_info);


int CPS_QIO_PROP_user_get_spin(CPS_QIO_PROP_UserRecordInfo *record_info);
int CPS_QIO_PROP_user_get_color(CPS_QIO_PROP_UserRecordInfo *record_info);
char *CPS_QIO_PROP_user_get_info(CPS_QIO_PROP_UserRecordInfo *record_info);

char *CPS_QIO_PROP_PAIRS_user_get_info(CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info);

char *CPS_QIO_SOURCE_user_get_info(CPS_QIO_SOURCE_UserRecordInfo *record_info);

int CPS_QIO_SOURCE_PAIRS_user_get_spin(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info);
int CPS_QIO_SOURCE_PAIRS_user_get_color(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info);
char *CPS_QIO_SOURCE_PAIRS_user_get_info(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info);



int CPS_QIO_PROP_file_defined_type(CPS_QIO_PROP_FileRecordInfo *record_info);
int CPS_QIO_PROP_file_defined_info(CPS_QIO_PROP_FileRecordInfo *record_info);

int CPS_QIO_PROP_user_defined_spin(CPS_QIO_PROP_UserRecordInfo *record_info);
int CPS_QIO_PROP_user_defined_color(CPS_QIO_PROP_UserRecordInfo *record_info);
int CPS_QIO_PROP_user_defined_info(CPS_QIO_PROP_UserRecordInfo *record_info);

int CPS_QIO_PROP_PAIRS_user_defined_info(CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info);

int CPS_QIO_SOURCE_user_defined_info(CPS_QIO_SOURCE_UserRecordInfo *record_info);

int CPS_QIO_SOURCE_PAIRS_user_defined_spin(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info);
int CPS_QIO_SOURCE_PAIRS_user_defined_color(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info);
int CPS_QIO_SOURCE_PAIRS_user_defined_info(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info);




int CPS_QIO_PROP_insert_filerecord_version(CPS_QIO_PROP_FileRecordInfo *record_info, char *version);

int CPS_QIO_PROP_insert_userrecord_version(CPS_QIO_PROP_UserRecordInfo *record_info, char *version);
int CPS_QIO_PROP_PAIRS_insert_userrecord_version(CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info, char *version);
int CPS_QIO_SOURCE_insert_userrecord_version(CPS_QIO_SOURCE_UserRecordInfo *record_info, char *version);
int CPS_QIO_SOURCE_PAIRS_insert_userrecord_version(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info, char *version);


int CPS_QIO_PROP_insert_userrecordinfo_info( CPS_QIO_PROP_UserRecordInfo *record_info, char *info);
int CPS_QIO_PROP_PAIRS_insert_userrecordinfo_info( CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info, char *info);
int CPS_QIO_SOURCE_insert_userrecordinfo_info( CPS_QIO_SOURCE_UserRecordInfo *record_info, char *info);
int CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_info( CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info, char *info);



int CPS_QIO_PROP_insert_userrecordinfo_spin( CPS_QIO_PROP_UserRecordInfo *record_info, int spin);
int CPS_QIO_PROP_insert_userrecordinfo_color( CPS_QIO_PROP_UserRecordInfo *record_info, int color);

int CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info, int spin);
int CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info, int color);




int CPS_QIO_PROP_insert_filerecordinfo_type( CPS_QIO_PROP_FileRecordInfo *record_info, char *type);
int CPS_QIO_PROP_insert_filerecordinfo_info( CPS_QIO_PROP_FileRecordInfo *record_info, char *info);


int CPS_QIO_PROP_insert_filerecord_tag_string(CPS_QIO_PROP_FileRecordInfoWrapper *wrapper,
					 char *recordinfo_tags);

int CPS_QIO_PROP_insert_userrecord_tag_string(CPS_QIO_PROP_UserRecordInfoWrapper *wrapper,
					 char *recordinfo_tags);
int CPS_QIO_PROP_PAIRS_insert_userrecord_tag_string(CPS_QIO_PROP_PAIRS_UserRecordInfoWrapper *wrapper,
					 char *recordinfo_tags);
int CPS_QIO_SOURCE_insert_userrecord_tag_string(CPS_QIO_SOURCE_UserRecordInfoWrapper *wrapper,
					 char *recordinfo_tags);
int CPS_QIO_SOURCE_PAIRS_insert_userrecord_tag_string(CPS_QIO_SOURCE_PAIRS_UserRecordInfoWrapper *wrapper,
					 char *recordinfo_tags);


char *CPS_QIO_PROP_get_file_record_info_tag_string(CPS_QIO_PROP_FileRecordInfoWrapper *wrapper);

char *CPS_QIO_PROP_get_user_record_info_tag_string(CPS_QIO_PROP_UserRecordInfoWrapper *wrapper);
char *CPS_QIO_PROP_PAIRS_get_user_record_info_tag_string(CPS_QIO_PROP_PAIRS_UserRecordInfoWrapper *wrapper);
char *CPS_QIO_SOURCE_get_user_record_info_tag_string(CPS_QIO_SOURCE_UserRecordInfoWrapper *wrapper);
char *CPS_QIO_SOURCE_PAIRS_get_user_record_info_tag_string(CPS_QIO_SOURCE_PAIRS_UserRecordInfoWrapper *wrapper);



CPS_END_NAMESPACE
#endif // __QIOXMLINFOPROP__
