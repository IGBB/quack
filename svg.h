#ifndef __SVG_H
#define __SVG_H

#define svg_attr(name, type, ...) " " #name "=" #type , __VA_ARGS__

void svg_start(const char* type, int num, ...);
void svg_end(const char* type);

#endif
