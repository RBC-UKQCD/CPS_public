/* Hacked by Peter Boyle for VML 2004 */
char *
stpcpy (dest, src)
     char *dest;
     const char *src;
{
  register char *d = dest;
  register const char *s = src;

  do
    *d++ = *s;
  while (*s++ != '\0');

  return d - 1;
}

