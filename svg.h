#ifndef __SVG_H
#define __SVG_H

/* Add an attribute to the svg tag. 
     - `name` = name of attribute. **must be a literal and cannot be quoted**
     - `type` = printf format literal, unquoted string 
     - `...`  = are variables passed to printf
 */
#define svg_attr(name, type, ...) " " #name "=" #type , __VA_ARGS__

/* Print an svg tag (i.e <type ...>). `svg_start_tag` doesn't add a final '/',
     while `svg_simple_tag` does.
     - `type` = type of element (e.g. 'svg', 'line', 'rect', etc)
     - `num`  = number of attributes 
     - `...`  = a list of svg_attr
 */
void svg_start_tag(const char* type, int num, ...);
void svg_simple_tag(const char* type, int num, ...);


/* Prints closing svg tag (i.e </type>).
     - `type` = type of element (e.g. 'svg', 'line', 'rect', etc)
 */
void svg_end_tag(const char* type);


#endif
