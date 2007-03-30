#ifdef USE_QIO
#include <config.h>
#include <util/qio_xmlinfo.h>

CPS_START_NAMESPACE
using namespace std;

CPS_QIO_UserRecordInfo *CPS_QIO_create_user_record_info(char *plaq, char *linktr, char *info)
{
  
  // taken from QIO_create_record_info and modified

  CPS_QIO_UserRecordInfo templ = CPS_QIO_USER_RECORD_INFO_TEMPLATE;
  CPS_QIO_UserRecordInfo *record_info;
  time_t cu_time;

  record_info = (CPS_QIO_UserRecordInfo *)malloc(sizeof(CPS_QIO_UserRecordInfo));
  if(!record_info)return NULL;
  time(&cu_time);

  *record_info = templ;
  CPS_QIO_insert_userrecord_version(record_info,CPS_QIO_USERRECORDFORMATVERSION);  
  
  CPS_QIO_insert_userrecordinfo_plaq( record_info, plaq);
  CPS_QIO_insert_userrecordinfo_linktr( record_info, linktr);
  CPS_QIO_insert_userrecordinfo_info( record_info, info);

  return record_info;

}

void CPS_QIO_destroy_user_record_info(CPS_QIO_UserRecordInfo *record_info){
  free(record_info);
}




void CPS_QIO_encode_user_record_info(QIO_String *record_string, 
				     CPS_QIO_UserRecordInfo *record_info)
{
  // taken from QIO_encode_record_info

  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  CPS_QIO_UserRecordInfoWrapper wrapper = CPS_QIO_USER_RECORD_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->plaq, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->linktr, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  CPS_QIO_insert_userrecord_tag_string(&wrapper, recordinfo_tags);

  /* Now build final XML string */
  QIO_string_realloc(record_string, QIO_STRINGALLOC);
  buf  = QIO_string_ptr(record_string);
  remainder = QIO_string_length(record_string);

  /* Begin with xml info stuff */
  strncpy(buf,QIO_XMLINFO,remainder);
  buf[remainder-1] = '\0';
  n = strlen(buf);
  remainder -= n;
  buf += n;
  if(remainder < 0){
    printf("CPS_QIO_encode_user_record_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.userrecordinfo_tags, &remainder);
  }
}

int CPS_QIO_decode_user_record_info(CPS_QIO_UserRecordInfo *record_info,
				    QIO_String *record_string)
{

  // taken from QIO_decode_record_info

  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  CPS_QIO_UserRecordInfoWrapper wrapper = CPS_QIO_USER_RECORD_INFO_WRAPPER;
  CPS_QIO_UserRecordInfo templ = CPS_QIO_USER_RECORD_INFO_TEMPLATE;
  char *left_angle;

  /* Initialize record info structure from a template */
  *record_info = templ;
  
  /* Start parsing record_string */
  /* Check leading tag, which is probably the info phrase "<?xml ...?>" */
  /* We ignore it if it is there */
  tmp_pt = QIO_next_tag(parse_pt, tag, &left_angle);
  if(strcmp(tag,QIO_QUESTXML)==0){
    /* Found ?xml, so resume parsing after the closing ">", ignoring
       the field. Otherwise, leave the parse_pt at its initial value */
    parse_pt = tmp_pt;
  }

  /* Open top-level tag (wrapper) and extract string containing tags */
  parse_pt = QIO_get_tag_value(parse_pt, tag, tags_string);
  QIO_decode_as_string (tag, tags_string, &wrapper.userrecordinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.userrecordinfo_tags))
    return QIO_BAD_XML;
  /* Otherwise start parsing the string of tags */
  parse_pt = CPS_QIO_get_user_record_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_string(tag,value_string,&record_info->plaq);
    QIO_decode_as_string(tag,value_string,&record_info->linktr);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_string_occur(&record_info->plaq);
  errors += QIO_check_string_occur(&record_info->linktr);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;


}

char *CPS_QIO_get_plaq(CPS_QIO_UserRecordInfo *record_info)
{
  return record_info->plaq.value;
}

char *CPS_QIO_get_linktr(CPS_QIO_UserRecordInfo *record_info)
{
  return record_info->linktr.value;
}

char *CPS_QIO_get_info(CPS_QIO_UserRecordInfo *record_info)
{
  return record_info->info.value;
}

int CPS_QIO_defined_plaq(CPS_QIO_UserRecordInfo *record_info)
{
  return record_info->plaq.occur;
}

int CPS_QIO_defined_linktr(CPS_QIO_UserRecordInfo *record_info)
{
  return record_info->linktr.occur;
}

int CPS_QIO_defined_info(CPS_QIO_UserRecordInfo *record_info)
{
  return record_info->info.occur;
}


int CPS_QIO_insert_userrecord_version(CPS_QIO_UserRecordInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;


}


int CPS_QIO_insert_userrecordinfo_plaq( CPS_QIO_UserRecordInfo *record_info, char *plaq)
{
  record_info->plaq.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->plaq.value, plaq, QIO_MAXVALUESTRING-1);
  record_info->plaq.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->plaq.occur = 1;
  if(strlen(plaq) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int CPS_QIO_insert_userrecordinfo_linktr( CPS_QIO_UserRecordInfo *record_info, char *linktr)
{
  record_info->linktr.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->linktr.value, linktr, QIO_MAXVALUESTRING-1);
  record_info->linktr.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->linktr.occur = 1;
  if(strlen(linktr) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int CPS_QIO_insert_userrecordinfo_info( CPS_QIO_UserRecordInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int CPS_QIO_insert_userrecord_tag_string(CPS_QIO_UserRecordInfoWrapper *wrapper,
                                 char *recordinfo_tags){
  wrapper->userrecordinfo_tags.occur = 0;
  if(!recordinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->userrecordinfo_tags.value, recordinfo_tags,
          QIO_MAXVALUESTRING-1);
  wrapper->userrecordinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->userrecordinfo_tags.occur = 1;
  if(strlen(recordinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

char *CPS_QIO_get_user_record_info_tag_string(CPS_QIO_UserRecordInfoWrapper *wrapper){
  return wrapper->userrecordinfo_tags.value;
}



CPS_END_NAMESPACE
#endif

