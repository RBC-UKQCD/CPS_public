// these are missing in qio_xml.h

char *QIO_encode_as_string(char *buf, QIO_TagCharValue *tag_value,
			   int *remainder);

char *QIO_encode_as_int(char *buf, QIO_TagIntValue *tag_value,
			int *remainder);

char *QIO_next_tag(char *parse_pt, char *tag, char **left_angle);

char *QIO_get_tag_value(char *parse_pt, char *tag, char *value_string);

void QIO_decode_as_string(char *tag, char *value_string,
			  QIO_TagCharValue *tag_value);

void QIO_decode_as_int(char *tag, char *value_string,
		       QIO_TagIntValue *tag_value);

int QIO_check_string_occur(QIO_TagCharValue *tag_value);

int QIO_check_int_occur(QIO_TagIntValue *tag_value);


// end missing
