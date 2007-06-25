#ifndef __QIOXMLINFO__
#define __QIOXMLINFO__

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
  QIO_TagCharValue plaq;
  QIO_TagCharValue linktr;
  QIO_TagCharValue info;
} CPS_QIO_UserRecordInfo;


#define CPS_QIO_USERRECORDFORMATVERSION "1.0"

#define CPS_QIO_USER_RECORD_INFO_TEMPLATE {\
   {"version", "", "", 0}, \
   {"plaq"   , "", "", 0}, \
   {"linktr" , "", "", 0}, \
   {"info"   , "", "", 0}  \
}


#define CPS_QIO_USER_RECORD_INFO_WRAPPER {\
  {"usqcdInfo", "", "" , 0}       \
}


typedef struct {
  QIO_TagCharValue userrecordinfo_tags;
} CPS_QIO_UserRecordInfoWrapper;

CPS_QIO_UserRecordInfo *CPS_QIO_create_user_record_info(char *plaq, char *linktr, char *info);
void CPS_QIO_destroy_user_record_info(CPS_QIO_UserRecordInfo *record_info);


void CPS_QIO_encode_user_record_info(QIO_String *record_string, 
				     CPS_QIO_UserRecordInfo *record_info);

int CPS_QIO_decode_user_record_info(CPS_QIO_UserRecordInfo *record_info,
				    QIO_String *record_string);

char *CPS_QIO_get_plaq(CPS_QIO_UserRecordInfo *record_info);
char *CPS_QIO_get_linktr(CPS_QIO_UserRecordInfo *record_info);
char *CPS_QIO_get_info(CPS_QIO_UserRecordInfo *record_info);

int CPS_QIO_defined_plaq(CPS_QIO_UserRecordInfo *record_info);
int CPS_QIO_defined_linktr(CPS_QIO_UserRecordInfo *record_info);
int CPS_QIO_defined_info(CPS_QIO_UserRecordInfo *record_info);


int CPS_QIO_insert_userrecord_version(CPS_QIO_UserRecordInfo *record_info, char *version);

int CPS_QIO_insert_userrecordinfo_plaq( CPS_QIO_UserRecordInfo *record_info, char *plaq);
int CPS_QIO_insert_userrecordinfo_linktr( CPS_QIO_UserRecordInfo *record_info, char *linktr);
int CPS_QIO_insert_userrecordinfo_info( CPS_QIO_UserRecordInfo *record_info, char *info);

int CPS_QIO_insert_userrecord_tag_string(CPS_QIO_UserRecordInfoWrapper *wrapper,
					 char *recordinfo_tags);

char *CPS_QIO_get_user_record_info_tag_string(CPS_QIO_UserRecordInfoWrapper *wrapper);

CPS_END_NAMESPACE
#endif // __QIOXMLINFO__
