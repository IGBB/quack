#include "svg.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INDENT "  "
int _svg_indent_level = 0;


char* svg_attr(const char* name, const char* fmt, ...){
  char *attr;
  int retval, length, cur_length;
  va_list vl;

  /* Use vsnprintf to get the length of value string before allocating memory*/
  va_start(vl, fmt);
  retval = vsnprintf(NULL, 0, fmt, vl);
  va_end(vl);
  if(retval < 0) return NULL;

  /* Allocate enough memory for entire attribute string: 
       ' $name="$value"'
       length($name) + length($value) + 4 literals [ =""] + null char
   */
  length = retval + strlen(name) + 5;
  attr = malloc(length);
  if(attr == NULL) return NULL;


  /* Print beginning of attribute: ' $name="' ; freeing memory if failed*/
  retval = snprintf(attr, length-1, " %s=\"", name);
  if(retval < 0){
    free(attr);
    return NULL;
  }
  cur_length = retval;

  /* Print attribute value; freeing memory if failed*/
  va_start(vl, fmt);
  retval = vsnprintf(attr + cur_length, length-cur_length-1, fmt, vl);
  va_end(vl);
  if(retval < 0){
    free(attr);
    return NULL;
  }
  cur_length += retval;

  /* Close quote and terminate string */
  attr[cur_length] = '"';
  attr[cur_length+1] = '\0';
  
  return attr;
}



void _svg_tag(const int simple, const char* type, int num, ...){
  int i = 0;
  va_list vl;
  char* attr;

  /* Indent current tag */
  for( i = 0; i < _svg_indent_level; i++)
    printf(INDENT);

  /* Increase indent level if started tag with internal elements */
  if(!simple)
    _svg_indent_level++;

  /* Open tag */
  printf("<%s", type);

  /* Print each attribute, free memory as used */
  va_start(vl, num);
  for(i = 0; i < num; i++){
    attr = va_arg(vl, char*);
    printf("%s", attr);
    free(attr);
  }
  va_end(vl);

  /* print nested closing tag if no internal elements expected */
  if(simple) printf("/");

  /* close tag */
  printf(">\n");
}


void svg_end_tag(const char* type){
  int i = 0;

  if(_svg_indent_level > 0)
    _svg_indent_level--;

  for( i = 0; i < _svg_indent_level; i++)
    printf(INDENT);

  
  printf("</%s>\n", type);
  
}
