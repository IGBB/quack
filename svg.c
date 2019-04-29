#include "svg.h"

#include <stdarg.h>
#include <stdio.h>


void svg_vprintf(int num, va_list vl){
  int i;
  
  for(i = 0; i < num; i++){
    vprintf(va_arg(vl, char*), vl);
  }

  
  return;
}


void svg_start_tag(const char* type, int num, ...){
  va_list vl;
  va_start(vl, num);

  printf("<%s", type);
  svg_vprintf(num, vl);
  printf(">\n");

  va_end(vl);

}
void svg_simple_tag(const char* type, int num, ...){
  va_list vl;
  va_start(vl, num);

  printf("<%s", type);
  svg_vprintf(num, vl);
  printf("/>\n");

  va_end(vl);
}


void svg_end_tag(const char* type){
  printf("</%s>", type);
}
