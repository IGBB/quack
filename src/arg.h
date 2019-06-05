#ifndef __ARG_H
#define __ARG_H

extern const char* const program_version;

struct arguments {
    char *name,
         *forward,
         *reverse,
         *unpaired,
         *adapters;
};

struct arguments parse_options(int argc, char **argv);
  
#endif
