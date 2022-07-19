#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>

#include "kseq.h"
#include "svg.h"
#include "arg.h"
#include "seq.h"

extern void draw_txt(
  FILE * output,
  sequence_data* forward, sequence_data* reverse,
  char* name, int adapters_used);
extern void draw_svg(
  FILE * output,
  sequence_data* forward, sequence_data* reverse,
  char* name, int adapters_used, int saturation);


int main (int argc, char **argv)
{
    struct arguments arguments;
    arguments = parse_options(argc, argv);

    int paired, unpaired, adapters;
    int *kmers = NULL;

    paired = (arguments.forward != NULL && arguments.reverse != NULL);
    unpaired = (arguments.unpaired != NULL);
    adapters = (arguments.adapters != NULL);
    
    /* If paired and unparied data are both set or unset, then throw error */
    if(paired == unpaired){
      fprintf(stderr, "Usage: quack [OPTION...]\nTry `quack --help' or `quack --usage' for more information.");
      exit(1);
    }

    if(adapters) kmers = read_adapters(arguments.adapters);

    sequence_data *forward, *reverse;
    if(paired){
      forward = read_fastq(arguments.forward, kmers, arguments.encoding, arguments.saturation);
      reverse = read_fastq(arguments.reverse, kmers, arguments.encoding, arguments.saturation);
    } else {
      forward = read_fastq(arguments.unpaired, kmers, arguments.encoding, arguments.saturation);
      reverse = NULL;
    }

    if(arguments.txt)
      draw_txt(arguments.txt, forward, reverse, arguments.name, adapters);

    if(arguments.svg)
      draw_svg(arguments.svg, forward, reverse, arguments.name, adapters, arguments.saturation);
  

    exit (0);
}
