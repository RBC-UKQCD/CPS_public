#ifdef USE_QIO
#include <config.h>
#include <util/qio_xmlinfo_prop.h>

CPS_START_NAMESPACE
using namespace std;


// first the user-record-stuff...


// QIO_PROP  (has spin, color)

CPS_QIO_PROP_UserRecordInfo *CPS_QIO_PROP_create_user_record_info(int spin, int color, char *info)
{
  
  // taken from QIO_create_record_info and modified

  CPS_QIO_PROP_UserRecordInfo templ = CPS_QIO_PROP_USER_RECORD_INFO_TEMPLATE;
  CPS_QIO_PROP_UserRecordInfo *record_info;
  time_t cu_time;

  record_info = (CPS_QIO_PROP_UserRecordInfo *)malloc(sizeof(CPS_QIO_PROP_UserRecordInfo));
  if(!record_info)return NULL;
  time(&cu_time);

  *record_info = templ;
  CPS_QIO_PROP_insert_userrecord_version(record_info,CPS_QIO_PROP_USERRECORDFORMATVERSION);  
  
  CPS_QIO_PROP_insert_userrecordinfo_spin( record_info, spin);
  CPS_QIO_PROP_insert_userrecordinfo_color( record_info, color);
  CPS_QIO_PROP_insert_userrecordinfo_info( record_info, info);

  return record_info;

}

void CPS_QIO_PROP_destroy_user_record_info(CPS_QIO_PROP_UserRecordInfo *record_info){
  free(record_info);
}




void CPS_QIO_PROP_encode_user_record_info(QIO_String *record_string, 
				     CPS_QIO_PROP_UserRecordInfo *record_info)
{
  // taken from QIO_encode_record_info

  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  CPS_QIO_PROP_UserRecordInfoWrapper wrapper = CPS_QIO_PROP_USER_RECORD_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_int(buf,&record_info->spin, &remainder);
  buf = QIO_encode_as_int(buf,&record_info->color, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  CPS_QIO_PROP_insert_userrecord_tag_string(&wrapper, recordinfo_tags);

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
    printf("CPS_QIO_PROP_encode_user_record_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.userrecordinfo_tags, &remainder);
  }
}

int CPS_QIO_PROP_decode_user_record_info(CPS_QIO_PROP_UserRecordInfo *record_info,
				    QIO_String *record_string)
{

  // taken from QIO_decode_record_info

  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  CPS_QIO_PROP_UserRecordInfoWrapper wrapper = CPS_QIO_PROP_USER_RECORD_INFO_WRAPPER;
  CPS_QIO_PROP_UserRecordInfo templ = CPS_QIO_PROP_USER_RECORD_INFO_TEMPLATE;
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
  parse_pt = CPS_QIO_PROP_get_user_record_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_int(tag,value_string,&record_info->spin);
    QIO_decode_as_int(tag,value_string,&record_info->color);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_int_occur(&record_info->spin);
  errors += QIO_check_int_occur(&record_info->color);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;


}

int CPS_QIO_PROP_user_get_spin(CPS_QIO_PROP_UserRecordInfo *record_info)
{
  return record_info->spin.value;
}

int CPS_QIO_PROP_user_get_color(CPS_QIO_PROP_UserRecordInfo *record_info)
{
  return record_info->color.value;
}

char *CPS_QIO_PROP_user_get_info(CPS_QIO_PROP_UserRecordInfo *record_info)
{
  return record_info->info.value;
}

int CPS_QIO_PROP_user_defined_spin(CPS_QIO_PROP_UserRecordInfo *record_info)
{
  return record_info->spin.occur;
}

int CPS_QIO_PROP_user_defined_color(CPS_QIO_PROP_UserRecordInfo *record_info)
{
  return record_info->color.occur;
}

int CPS_QIO_PROP_user_defined_info(CPS_QIO_PROP_UserRecordInfo *record_info)
{
  return record_info->info.occur;
}


int CPS_QIO_PROP_insert_userrecord_version(CPS_QIO_PROP_UserRecordInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;


}


int CPS_QIO_PROP_insert_userrecordinfo_spin( CPS_QIO_PROP_UserRecordInfo *record_info, int spin)
{
  record_info->spin.occur = 0;
  if(!record_info)return QIO_BAD_ARG;

  record_info->spin.value = spin;
  record_info->spin.occur = 1;
  return QIO_SUCCESS;

}

int CPS_QIO_PROP_insert_userrecordinfo_color( CPS_QIO_PROP_UserRecordInfo *record_info, int color)
{
  record_info->color.occur = 0;
  if(!record_info)return QIO_BAD_ARG;

  record_info->color.value = color;
  record_info->color.occur = 1;
  return QIO_SUCCESS;

}

int CPS_QIO_PROP_insert_userrecordinfo_info( CPS_QIO_PROP_UserRecordInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int CPS_QIO_PROP_insert_userrecord_tag_string(CPS_QIO_PROP_UserRecordInfoWrapper *wrapper,
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

char *CPS_QIO_PROP_get_user_record_info_tag_string(CPS_QIO_PROP_UserRecordInfoWrapper *wrapper){
  return wrapper->userrecordinfo_tags.value;
}


// copy and pace START

// copy 1 -> PROP_PAIRS

// QIO_PROP_PAIRS ( _no_ color, spin)

CPS_QIO_PROP_PAIRS_UserRecordInfo *CPS_QIO_PROP_PAIRS_create_user_record_info(char *info)
{
  
  // taken from QIO_create_record_info and modified

  CPS_QIO_PROP_PAIRS_UserRecordInfo templ = CPS_QIO_PROP_PAIRS_USER_RECORD_INFO_TEMPLATE;
  CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info;
  time_t cu_time;

  record_info = (CPS_QIO_PROP_PAIRS_UserRecordInfo *)malloc(sizeof(CPS_QIO_PROP_PAIRS_UserRecordInfo));
  if(!record_info)return NULL;
  time(&cu_time);

  *record_info = templ;
  CPS_QIO_PROP_PAIRS_insert_userrecord_version(record_info,CPS_QIO_PROP_PAIRS_USERRECORDFORMATVERSION);  
  
  CPS_QIO_PROP_PAIRS_insert_userrecordinfo_info( record_info, info);

  return record_info;

}

void CPS_QIO_PROP_PAIRS_destroy_user_record_info(CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info){
  free(record_info);
}




void CPS_QIO_PROP_PAIRS_encode_user_record_info(QIO_String *record_string, 
				     CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info)
{
  // taken from QIO_encode_record_info

  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  CPS_QIO_PROP_PAIRS_UserRecordInfoWrapper wrapper = CPS_QIO_PROP_PAIRS_USER_RECORD_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  CPS_QIO_PROP_PAIRS_insert_userrecord_tag_string(&wrapper, recordinfo_tags);

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
    printf("CPS_QIO_PROP_PAIRS_encode_user_record_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.userrecordinfo_tags, &remainder);
  }
}

int CPS_QIO_PROP_PAIRS_decode_user_record_info(CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info,
				    QIO_String *record_string)
{

  // taken from QIO_decode_record_info

  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  CPS_QIO_PROP_PAIRS_UserRecordInfoWrapper wrapper = CPS_QIO_PROP_PAIRS_USER_RECORD_INFO_WRAPPER;
  CPS_QIO_PROP_PAIRS_UserRecordInfo templ = CPS_QIO_PROP_PAIRS_USER_RECORD_INFO_TEMPLATE;
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
  parse_pt = CPS_QIO_PROP_PAIRS_get_user_record_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;


}



char *CPS_QIO_PROP_PAIRS_user_get_info(CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info)
{
  return record_info->info.value;
}


int CPS_QIO_PROP_PAIRS_user_defined_info(CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info)
{
  return record_info->info.occur;
}


int CPS_QIO_PROP_PAIRS_insert_userrecord_version(CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;


}




int CPS_QIO_PROP_PAIRS_insert_userrecordinfo_info( CPS_QIO_PROP_PAIRS_UserRecordInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int CPS_QIO_PROP_PAIRS_insert_userrecord_tag_string(CPS_QIO_PROP_PAIRS_UserRecordInfoWrapper *wrapper,
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

char *CPS_QIO_PROP_PAIRS_get_user_record_info_tag_string(CPS_QIO_PROP_PAIRS_UserRecordInfoWrapper *wrapper){
  return wrapper->userrecordinfo_tags.value;
}



// copy 2 -> SOURCE

// QIO_SOURCE ( _no_ color, spin )

CPS_QIO_SOURCE_UserRecordInfo *CPS_QIO_SOURCE_create_user_record_info(char *info)
{
  
  // taken from QIO_create_record_info and modified

  CPS_QIO_SOURCE_UserRecordInfo templ = CPS_QIO_SOURCE_USER_RECORD_INFO_TEMPLATE;
  CPS_QIO_SOURCE_UserRecordInfo *record_info;
  time_t cu_time;

  record_info = (CPS_QIO_SOURCE_UserRecordInfo *)malloc(sizeof(CPS_QIO_SOURCE_UserRecordInfo));
  if(!record_info)return NULL;
  time(&cu_time);

  *record_info = templ;
  CPS_QIO_SOURCE_insert_userrecord_version(record_info,CPS_QIO_SOURCE_USERRECORDFORMATVERSION);  
  
  CPS_QIO_SOURCE_insert_userrecordinfo_info( record_info, info);

  return record_info;

}

void CPS_QIO_SOURCE_destroy_user_record_info(CPS_QIO_SOURCE_UserRecordInfo *record_info){
  free(record_info);
}




void CPS_QIO_SOURCE_encode_user_record_info(QIO_String *record_string, 
				     CPS_QIO_SOURCE_UserRecordInfo *record_info)
{
  // taken from QIO_encode_record_info

  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  CPS_QIO_SOURCE_UserRecordInfoWrapper wrapper = CPS_QIO_SOURCE_USER_RECORD_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  CPS_QIO_SOURCE_insert_userrecord_tag_string(&wrapper, recordinfo_tags);

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
    printf("CPS_QIO_SOURCE_encode_user_record_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.userrecordinfo_tags, &remainder);
  }
}

int CPS_QIO_SOURCE_decode_user_record_info(CPS_QIO_SOURCE_UserRecordInfo *record_info,
				    QIO_String *record_string)
{

  // taken from QIO_decode_record_info

  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  CPS_QIO_SOURCE_UserRecordInfoWrapper wrapper = CPS_QIO_SOURCE_USER_RECORD_INFO_WRAPPER;
  CPS_QIO_SOURCE_UserRecordInfo templ = CPS_QIO_SOURCE_USER_RECORD_INFO_TEMPLATE;
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
  parse_pt = CPS_QIO_SOURCE_get_user_record_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;


}


char *CPS_QIO_SOURCE_user_get_info(CPS_QIO_SOURCE_UserRecordInfo *record_info)
{
  return record_info->info.value;
}


int CPS_QIO_SOURCE_user_defined_info(CPS_QIO_SOURCE_UserRecordInfo *record_info)
{
  return record_info->info.occur;
}


int CPS_QIO_SOURCE_insert_userrecord_version(CPS_QIO_SOURCE_UserRecordInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;


}






int CPS_QIO_SOURCE_insert_userrecordinfo_info( CPS_QIO_SOURCE_UserRecordInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int CPS_QIO_SOURCE_insert_userrecord_tag_string(CPS_QIO_SOURCE_UserRecordInfoWrapper *wrapper,
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

char *CPS_QIO_SOURCE_get_user_record_info_tag_string(CPS_QIO_SOURCE_UserRecordInfoWrapper *wrapper){
  return wrapper->userrecordinfo_tags.value;
}



// copy 3 -> SOURCE_PAIRS

// QIO_SOURCE_PAIRS (has color, spin)

CPS_QIO_SOURCE_PAIRS_UserRecordInfo *CPS_QIO_SOURCE_PAIRS_create_user_record_info(int spin, int color, char *info)
{
  
  // taken from QIO_create_record_info and modified

  CPS_QIO_SOURCE_PAIRS_UserRecordInfo templ = CPS_QIO_SOURCE_PAIRS_USER_RECORD_INFO_TEMPLATE;
  CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info;
  time_t cu_time;

  record_info = (CPS_QIO_SOURCE_PAIRS_UserRecordInfo *)malloc(sizeof(CPS_QIO_SOURCE_PAIRS_UserRecordInfo));
  if(!record_info)return NULL;
  time(&cu_time);

  *record_info = templ;
  CPS_QIO_SOURCE_PAIRS_insert_userrecord_version(record_info,CPS_QIO_SOURCE_PAIRS_USERRECORDFORMATVERSION);  
  
  CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( record_info, spin);
  CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( record_info, color);
  CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_info( record_info, info);

  return record_info;

}

void CPS_QIO_SOURCE_PAIRS_destroy_user_record_info(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info){
  free(record_info);
}




void CPS_QIO_SOURCE_PAIRS_encode_user_record_info(QIO_String *record_string, 
				     CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info)
{
  // taken from QIO_encode_record_info

  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  CPS_QIO_SOURCE_PAIRS_UserRecordInfoWrapper wrapper = CPS_QIO_SOURCE_PAIRS_USER_RECORD_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_int(buf,&record_info->spin, &remainder);
  buf = QIO_encode_as_int(buf,&record_info->color, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  CPS_QIO_SOURCE_PAIRS_insert_userrecord_tag_string(&wrapper, recordinfo_tags);

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
    printf("CPS_QIO_SOURCE_PAIRS_encode_user_record_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.userrecordinfo_tags, &remainder);
  }
}

int CPS_QIO_SOURCE_PAIRS_decode_user_record_info(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info,
				    QIO_String *record_string)
{

  // taken from QIO_decode_record_info

  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  CPS_QIO_SOURCE_PAIRS_UserRecordInfoWrapper wrapper = CPS_QIO_SOURCE_PAIRS_USER_RECORD_INFO_WRAPPER;
  CPS_QIO_SOURCE_PAIRS_UserRecordInfo templ = CPS_QIO_SOURCE_PAIRS_USER_RECORD_INFO_TEMPLATE;
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
  parse_pt = CPS_QIO_SOURCE_PAIRS_get_user_record_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_int(tag,value_string,&record_info->spin);
    QIO_decode_as_int(tag,value_string,&record_info->color);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_int_occur(&record_info->spin);
  errors += QIO_check_int_occur(&record_info->color);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;


}

int CPS_QIO_SOURCE_PAIRS_user_get_spin(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info)
{
  return record_info->spin.value;
}

int CPS_QIO_SOURCE_PAIRS_user_get_color(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info)
{
  return record_info->color.value;
}

char *CPS_QIO_SOURCE_PAIRS_user_get_info(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info)
{
  return record_info->info.value;
}

int CPS_QIO_SOURCE_PAIRS_user_defined_spin(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info)
{
  return record_info->spin.occur;
}

int CPS_QIO_SOURCE_PAIRS_user_defined_color(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info)
{
  return record_info->color.occur;
}

int CPS_QIO_SOURCE_PAIRS_user_defined_info(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info)
{
  return record_info->info.occur;
}


int CPS_QIO_SOURCE_PAIRS_insert_userrecord_version(CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;


}


int CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_spin( CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info, int spin)
{
  record_info->spin.occur = 0;
  if(!record_info)return QIO_BAD_ARG;

  record_info->spin.value = spin;
  record_info->spin.occur = 1;
  return QIO_SUCCESS;

}

int CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_color( CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info, int color)
{
  record_info->color.occur = 0;
  if(!record_info)return QIO_BAD_ARG;

  record_info->color.value = color;
  record_info->color.occur = 1;
  return QIO_SUCCESS;

}

int CPS_QIO_SOURCE_PAIRS_insert_userrecordinfo_info( CPS_QIO_SOURCE_PAIRS_UserRecordInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int CPS_QIO_SOURCE_PAIRS_insert_userrecord_tag_string(CPS_QIO_SOURCE_PAIRS_UserRecordInfoWrapper *wrapper,
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

char *CPS_QIO_SOURCE_PAIRS_get_user_record_info_tag_string(CPS_QIO_SOURCE_PAIRS_UserRecordInfoWrapper *wrapper){
  return wrapper->userrecordinfo_tags.value;
}



// copy and pace END



// now for the file-record


CPS_QIO_PROP_FileRecordInfo *CPS_QIO_PROP_create_file_record_info(char *type, char *info)
{
  
  // taken from QIO_create_record_info and modified

  CPS_QIO_PROP_FileRecordInfo templ = CPS_QIO_PROP_FILE_RECORD_INFO_TEMPLATE;
  CPS_QIO_PROP_FileRecordInfo *record_info;
  time_t cu_time;

  record_info = (CPS_QIO_PROP_FileRecordInfo *)malloc(sizeof(CPS_QIO_PROP_FileRecordInfo));
  if(!record_info)return NULL;
  time(&cu_time);

  *record_info = templ;
  CPS_QIO_PROP_insert_filerecord_version(record_info,CPS_QIO_PROP_FILERECORDFORMATVERSION);  
  
  CPS_QIO_PROP_insert_filerecordinfo_type( record_info, type);
  CPS_QIO_PROP_insert_filerecordinfo_info( record_info, info);

  return record_info;

}

void CPS_QIO_PROP_destroy_file_record_info(CPS_QIO_PROP_FileRecordInfo *record_info){
  free(record_info);
}




void CPS_QIO_PROP_encode_file_record_info(QIO_String *record_string, 
				     CPS_QIO_PROP_FileRecordInfo *record_info)
{
  // taken from QIO_encode_record_info

  char *buf;
  int remainder,n;
  char recordinfo_tags[QIO_MAXVALUESTRING];
  CPS_QIO_PROP_FileRecordInfoWrapper wrapper = CPS_QIO_PROP_FILE_RECORD_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = recordinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&record_info->version, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->type, &remainder);
  buf = QIO_encode_as_string(buf,&record_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  CPS_QIO_PROP_insert_filerecord_tag_string(&wrapper, recordinfo_tags);

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
    printf("CPS_QIO_PROP_encode_file_record_info: record_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.userrecordinfo_tags, &remainder);
  }
}

int CPS_QIO_PROP_decode_file_record_info(CPS_QIO_PROP_FileRecordInfo *record_info,
				    QIO_String *record_string)
{

  // taken from QIO_decode_record_info

  char *parse_pt = QIO_string_ptr(record_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  CPS_QIO_PROP_FileRecordInfoWrapper wrapper = CPS_QIO_PROP_FILE_RECORD_INFO_WRAPPER;
  CPS_QIO_PROP_FileRecordInfo templ = CPS_QIO_PROP_FILE_RECORD_INFO_TEMPLATE;
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
  parse_pt = CPS_QIO_PROP_get_file_record_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&record_info->version);
    QIO_decode_as_string(tag,value_string,&record_info->type);
    QIO_decode_as_string(tag,value_string,&record_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&record_info->version);
  errors += QIO_check_string_occur(&record_info->type);
  errors += QIO_check_string_occur(&record_info->info);

  return errors;


}

char *CPS_QIO_PROP_file_get_type(CPS_QIO_PROP_FileRecordInfo *record_info)
{
  return record_info->type.value;
}


char *CPS_QIO_PROP_file_get_info(CPS_QIO_PROP_FileRecordInfo *record_info)
{
  return record_info->info.value;
}

int CPS_QIO_PROP_file_defined_type(CPS_QIO_PROP_FileRecordInfo *record_info)
{
  return record_info->type.occur;
}


int CPS_QIO_PROP_file_defined_info(CPS_QIO_PROP_FileRecordInfo *record_info)
{
  return record_info->info.occur;
}


int CPS_QIO_PROP_insert_filerecord_version(CPS_QIO_PROP_FileRecordInfo *record_info, char *version)
{
  record_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(record_info->version.value, version, QIO_MAXVALUESTRING-1);
  record_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;


}


int CPS_QIO_PROP_insert_filerecordinfo_type( CPS_QIO_PROP_FileRecordInfo *record_info, char *type)
{
  record_info->type.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->type.value, type, QIO_MAXVALUESTRING-1);
  record_info->type.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->type.occur = 1;
  if(strlen(type) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}


int CPS_QIO_PROP_insert_filerecordinfo_info( CPS_QIO_PROP_FileRecordInfo *record_info, char *info)
{
  record_info->info.occur = 0;
  if(!record_info)return QIO_BAD_ARG;
  strncpy(record_info->info.value, info, QIO_MAXVALUESTRING-1);
  record_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  record_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

int CPS_QIO_PROP_insert_filerecord_tag_string(CPS_QIO_PROP_FileRecordInfoWrapper *wrapper,
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

char *CPS_QIO_PROP_get_file_record_info_tag_string(CPS_QIO_PROP_FileRecordInfoWrapper *wrapper){
  return wrapper->userrecordinfo_tags.value;
}










CPS_END_NAMESPACE
#endif

