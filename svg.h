#ifndef __SVG_H
#define __SVG_H

#define svg_attr(name, type, ...) " " #name "=" #type , __VA_ARGS__

void svg_start_tag(const char* type, int num, ...);
void svg_end_tag(const char* type);
void svg_simple_tag(const char* type, int num, ...);

#endif
