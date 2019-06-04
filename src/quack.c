#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>

#include "kseq.h"
#include "svg.h"
#include "arg.h"
#include "seq.h"

extern void draw(sequence_data* data, int position, int adapters_used);




int main (int argc, char **argv)
{
    struct arguments arguments;
    arguments = parse_options(argc, argv);

    int paired, unpaired, adapters;
    int width, height;
    int *kmers = NULL;

    paired = (arguments.forward != NULL && arguments.reverse != NULL);
    unpaired = (arguments.unpaired != NULL);
    adapters = (arguments.adapters != NULL);
    
    /* If paired and unparied data are both set or unset, then throw error */
    if(paired == unpaired){
      printf("%s\n", "Usage: quack [OPTION...]\nTry `quack --help' or `quack --usage' for more information.");
      exit(1);
    }

    if(adapters) kmers = read_adapters(arguments.adapters);

    width  = (paired)?1195:615;
    height = (adapters)?610:510;

    if(arguments.name != NULL)
      height += 30;
    
    svg_start_tag("svg", 5,
                  svg_attr("width",   "%d", width),
                  svg_attr("height",  "%d", height),
                  svg_attr("viewBox", "%d %d %d %d", 0, 0, width, height),
                  svg_attr("xmlns",       "%s", "http://www.w3.org/2000/svg"),
                  svg_attr("xmlns:xlink", "%s", "http://www.w3.org/1999/xlink")
                  );

    // If name is given, add to middle of viewBox (half of width + min-x of viewbox)
    if(arguments.name != NULL){
      svg_start_tag("text", 6,
                    svg_attr("x", "%d", (width/2)),
                    svg_attr("y", "%d", 30),
                    svg_attr("font-family", "%s", "sans-serif"),
                    svg_attr("text-anchor", "%s", "middle"),
                    svg_attr("font-size",   "%s", "30px"),
                    svg_attr("fill",        "%s", "black"));
      printf("%s", arguments.name);
      svg_end_tag("text");

      svg_start_tag("g", 1, 
                svg_attr("transform", "translate(%d %d)", 0, 30)
                );

    }
  
    sequence_data *data = read_fastq(((paired)?arguments.forward:arguments.unpaired), kmers);
    sequence_data *transformed_data = transform(data);
    draw(transformed_data, 0, adapters);
    free(data);
    
    if(paired){
      data = read_fastq(arguments.reverse, kmers);
      transformed_data = transform(data);
      draw(transformed_data, 1, adapters);
      free(data);
    }

    if(arguments.name != NULL) svg_end_tag("g");

    svg_end_tag("svg");
    
    exit (0);
}
