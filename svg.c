#include "svg.h"

#include <stdarg.h>
#include <stdio.h>


#define INDENT "  "
int _svg_indent_level = 0;

void svg_vprintf(int num, va_list vl){
  int i;
  
  for(i = 0; i < num; i++){
    vprintf(va_arg(vl, char*), vl);
  }

  
  return;
}


void svg_start_tag(const char* type, int num, ...){
  int i = 0;

  for( i = 0; i < _svg_indent_level; i++)
    printf(INDENT);

  _svg_indent_level++;


  va_list vl;
  va_start(vl, num);

  printf("<%s", type);
  svg_vprintf(num, vl);
  printf(">\n");

  va_end(vl);

}
void svg_simple_tag(const char* type, int num, ...){
  int i = 0;

  for( i = 0; i < _svg_indent_level; i++)
    printf(INDENT);

  va_list vl;
  va_start(vl, num);

  printf("<%s", type);
  svg_vprintf(num, vl);
  printf("/>\n");

  va_end(vl);
}


void svg_end_tag(const char* type){
  int i = 0;

  if(_svg_indent_level > 0)
    _svg_indent_level--;

  for( i = 0; i < _svg_indent_level; i++)
    printf(INDENT);

  
  printf("</%s>\n", type);
  
}
