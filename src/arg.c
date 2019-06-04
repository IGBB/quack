#include "arg.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

const char* const program_version = "quack 2.0";

struct arguments parse_options(int argc, char **argv) {
  struct arguments arguments = {
                                .unpaired = NULL,
                                .forward = NULL,
                                .reverse = NULL,
                                .name = NULL,
                                .adapters = NULL
  };

    if (argc== 1 || argc == 2)  {
        if (argc == 1 || (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "--usage") == 0 || strcmp(argv[1], "-?") == 0)) {
            printf("Usage: quack [OPTION...]\n"
            "quack -- A FASTQ quality assessment tool\n\n"
            "  -1, --forward file.1.fq.gz      Forward strand\n"
            "  -2, --reverse file.2.fq.gz      Reverse strand\n"
            "  -a, --adapters adapters.fa.gz   (Optional) Adapters file\n"
            "  -n, --name NAME                 (Optional) Display in output\n"
            "  -u, --unpaired unpaired.fq.gz   Data (only use with -u)\n"
            "  -?, --help                      Give this help list\n"
            "      --usage                     (use alone)\n"
            "  -V, --version                   Print program version (use alone)\n"
            "Report bugs to <thrash@igbb.msstate.edu>.\n");
        }

        if (argc == 2 && ((strcmp(argv[1], "-V") == 0 || strcmp(argv[1], "--version") == 0))) {
            printf("%s\n", program_version);
            exit (0);
        }
    }

    if(argc>2 && argc % 2 != 0)
    {
        int counter=1;
        while (counter<argc) {
            if (strcmp(argv[counter], "--forward") == 0 || strcmp(argv[counter], "-1") == 0) {
                arguments.forward = argv[counter+1];
            }

            else if (strcmp(argv[counter], "--reverse") == 0 || strcmp(argv[counter], "-2") == 0) {
                arguments.reverse = argv[counter+1];
            }

            else if (strcmp(argv[counter], "--adapters") == 0 || strcmp(argv[counter], "-a") == 0) {
                arguments.adapters = argv[counter+1];
            }

            else if (strcmp(argv[counter], "--unpaired") == 0 || strcmp(argv[counter], "-u") == 0) {
                arguments.unpaired = argv[counter+1];
            }

            else if (strcmp(argv[counter], "--name") == 0 || strcmp(argv[counter], "-n") == 0) {
                arguments.name = argv[counter+1];
            }
            else {
                printf("Usage: quack [OPTION...]\n"
                "quack -- A FASTQ quality assessment tool\n\n"
                "  -1, --forward file.1.fq.gz      Forward strand\n"
                "  -2, --reverse file.2.fq.gz      Reverse strand\n"
                "  -a, --adapters adapters.fa.gz    Adapters file\n"
                "  -n, --name NAME            Display in output\n"
                "  -u, --unpaired unpaired.fq.gz        Data (only use with -u)\n"
                "  -?, --help                 Give this help list\n"
                "      --usage                (use alone)\n"
                "  -V, --version              Print program version (use alone)\n"
                "Report bugs to <thrash@igbb.msstate.edu>.\n");
            }

            counter = counter+2;
        }
    }
    return arguments;
}
