#include <stdio.h>

#include "svg.h"

int main (int argc, char **argv)
{

  int width=800, height=1024;

  printf("<svg width='%d' height='%d' viewBox='20 0 %d %d' " 
           "xmlns='http://www.w3.org/2000/svg' " 
           "xmlns:xlink='http://www.w3.org/1999/xlink'>\n", 
           width, height, width, height);


  svg_start("svg", 5,
            svg_attr(width,   "%d", width),
            svg_attr(height,  "%d", height),
            svg_attr(viewBox, "%d %d %d %d", 20, 0, width, height),
            svg_attr(xmlns,       "%s", "http://www.w3.org/2000/svg"),
            svg_attr(xmlns:xlink, "%s", "http://www.w3.org/1999/xlink")
            );



  
}
