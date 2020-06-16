#ifndef __ARG_H
#define __ARG_H

#include "seq.h"

extern const char* const program_version;

typedef enum {
    TXT, SVG
} output_format;

struct arguments {
   char *name,
       *forward,
       *reverse,
       *unpaired,
       *adapters;
   output_format outfmt;
   encoding_t encoding;
};

struct arguments parse_options(int argc, char **argv);
  
#endif
