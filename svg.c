#include "svg.h"

#include <stdarg.h>
#include <stdio.h>

void svg_start(const char* type, int num, ...){
  char format[1000];
  int length = 0;
  printf("<%s", type);
  
  va_list vl;
  va_start(vl, num);
  for( int i = 0; i < num; i++){
    vprintf(va_arg(vl, char*), vl);
  }
  va_end(vl);

  printf(">\n");
  
  return;
}

void svg_end(const char* type){
  printf("</%s>", type);
}
