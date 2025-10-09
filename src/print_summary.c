#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "seq.h"

void print_entry (FILE*output, sequence_data* seq,
                  char* name, char* direction ){

    if(name)
        fprintf(output, "%s\t", name);
    if(direction)
        fprintf(output, "%s\t", direction);

    fprintf(output, "%lu\t", seq->number_of_sequences);

    size_t i;
    uint64_t size = 0;
    for (i = 0; i < seq->max_length; i++)
        size += seq->bases[i].length_count * (i + 1);
    fprintf(output, "%lu\t", size);

    uint64_t tmp = 0;
    for (i = 0; i < seq->max_length && 2 * tmp < size; i++)
        tmp += seq->bases[i].length_count * (i + 1);
    fprintf(output, "%lu\n", i + 1);
}

void print_summary(
    FILE * output,
    sequence_data* forward,
    sequence_data* reverse,
    char* name,
    int adapters_used){

    if(name)
        fprintf(output, "Name\t");
    if(reverse)
        fprintf(output, "Direction\t");
    fprintf(output, "NumReads\t");
    fprintf(output, "NumBases\t");
    fprintf(output, "N50\n");

    if(reverse){
        print_entry(output, forward, name, "Forward");
        print_entry(output, reverse, name, "Reverse");
    }else{
        print_entry(output, forward, name, NULL);
    }

}
