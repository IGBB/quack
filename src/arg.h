#ifndef __ARG_H
#define __ARG_H

#include <stdio.h>
#include "seq.h"

extern const char* const program_version;

struct arguments {
   char *name,
       *forward,
       *reverse,
       *unpaired,
       *adapters;
    FILE * svg, * txt;
    encoding_t encoding;
    int saturation;
};

struct arguments parse_options(int argc, char **argv);
  
#endif
