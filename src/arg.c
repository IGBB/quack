#include "arg.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "ketopt.h"

const char* const program_version = "quack 2.0";

const char* const help_message = ""
  "Usage: quack [OPTION...]\n"
  "quack -- A FASTQ quality assessment tool\n\n"
  "  -1, --forward file.1.fq.gz      Forward strand\n"
  "  -2, --reverse file.2.fq.gz      Reverse strand\n"
  "  -a, --adapters adapters.fa.gz   (Optional) Adapters file\n"
  "  -n, --name NAME                 (Optional) Display in output\n"
  "  -u, --unpaired unpaired.fq.gz   Data (only use with -u)\n"
  "  -?, --help                      Give this help list\n"
  "      --usage                     (use alone)\n"
  "  -V, --version                   Print program version (use alone)\n"
  "Report bugs to <thrash@igbb.msstate.edu>.\n";

static ko_longopt_t longopts[] = {
    { "forward",  ko_required_argument, '1' },
    { "reverse",  ko_required_argument, '2' },
    { "unpaired", ko_required_argument, 'u' },
    { "adapters", ko_required_argument, 'a' },

    { "name",     ko_required_argument, 'n' },

    { "help",    ko_no_argument, 'h' },
    { "usage",   ko_no_argument, 'h' },
    { "version", ko_no_argument, 'V' },

    { NULL, 0, 0 }
  };

struct arguments parse_options(int argc, char **argv) {
  struct arguments arguments = {
                                .unpaired = NULL,
                                .forward = NULL,
                                .reverse = NULL,
                                .name = NULL,
                                .adapters = NULL
  };


  ketopt_t opt = KETOPT_INIT;

  int i, c;
  while ((c = ketopt(&opt, argc, argv, 1, "1:2:a:u:n:?V", longopts)) >= 0) {
    switch(c){
    case '1': arguments.forward  = opt.arg; break;
    case '2': arguments.reverse  = opt.arg; break;
    case 'u': arguments.unpaired = opt.arg; break;
    case 'a': arguments.adapters = opt.arg; break;

    case 'n': arguments.name = opt.arg;     break;

    case 'V':
      printf("%s\n", program_version);
      exit(EXIT_SUCCESS);

    case 'h':
    case '?':
    default:
      printf("%s\n", help_message);
      exit(EXIT_SUCCESS);

    };
  }

    return arguments;
}
