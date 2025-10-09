#ifndef __SVG_H
#define __SVG_H

#include <stdio.h>

/* Add an attribute to the svg tag. 
     - `name` = name of attribute. 
     - `fmt`  = printf format
     - `...`  = are variables passed to printf

  Allocates memory for return string. MUST be free'd after use 
 */
char* svg_attr(const char* name, const char* fmt, ...);

/* Print an svg tag (i.e <type ...>). `svg_start_tag` doesn't add a final '/',
     while `svg_simple_tag` does.
     - `type` = type of element (e.g. 'svg', 'line', 'rect', etc)
     - `num`  = number of attributes 
     - `...`  = a list of svg_attr
 */
void _svg_tag(FILE* svg, const int simple, const char* type, int num, ...);
#define svg_start_tag(output, ...) _svg_tag(output, 0, __VA_ARGS__)
#define svg_simple_tag(output, ...)  _svg_tag(output, 1, __VA_ARGS__)


/* Prints closing svg tag (i.e </type>).
     - `type` = type of element (e.g. 'svg', 'line', 'rect', etc)
 */
void svg_end_tag(FILE* svg, const char* type);


#endif
