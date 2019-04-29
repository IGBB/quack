#include <stdio.h>

#include "svg.h"

int main (int argc, char **argv)
{

  int width=800, height=1024;


  svg_start_tag("svg", 5,
                svg_attr(width,   "%d", width),
                svg_attr(height,  "%d", height),
                svg_attr(viewBox, "%d %d %d %d", 20, 0, width, height),
                svg_attr(xmlns,       "%s", "http://www.w3.org/2000/svg"),
                svg_attr(xmlns:xlink, "%s", "http://www.w3.org/1999/xlink")
                );

  svg_simple_tag("line", 5,
                 svg_attr(x1,   "%d", 20),
                 svg_attr(x2,  "%d", 100),
                 svg_attr(y1,   "%d", 20),
                 svg_attr(y2,  "%d", 100),
                 svg_attr(stroke, "%s", "black")
                 );


  svg_end_tag("svg");
  
}
