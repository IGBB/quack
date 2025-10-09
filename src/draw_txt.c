#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "seq.h"


void draw_txt_content(FILE* output, sequence_data* data){

    int i, j;

    fprintf(output, "### Per Base Content \n\n");
    fprintf(output, "|   Base |          A          |          T          |          C          |          G          |\n");
    fprintf(output, "| ------ | ------------------- | ------------------- | ------------------- | ------------------- |\n");

    for(i = 0; i < data->max_length; i++){
        fprintf(output, "| % 6d ", i);
        for(j = 0; j < 4; j++){
            fprintf(output, "| % 10d (%5.2f%%) ", data->bases[i].content[j],
                   (data->bases[i].content[j] * 100.0)/data->bases[i].length_count);
        }
        fprintf(output, "|\n");
    }
    fprintf(output, "\n");
}


void draw_txt_quality(FILE* output, sequence_data* data){
    int i, j;
    uint64_t total[91] = {0}, tt = 0;


    fprintf(output, "### Per Base Sequence Quality \n\n");
    fprintf(output, "|   Base ");
    for(i = data->min_score; i <= data->max_score; i++)
        fprintf(output, "| %17d ", i);
    fprintf(output, "|\n");

    fprintf(output, "| ------ ");
    for(i = data->min_score; i <= data->max_score; i++)
        fprintf(output, "| ----------------- ");
    fprintf(output, "|\n");

    for(i = 0; i < data->max_length; i++){
        fprintf(output, "| % 6d ", i);
        for(j = data->min_score; j <= data->max_score; j++){
            tt += data->bases[i].scores[j];
            total[j] += data->bases[i].scores[j];
            fprintf(output, "| % 8d (%5.2f%%) ", data->bases[i].scores[j],
                   (data->bases[i].scores[j] * 100.0)/data->bases[i].length_count);
        }
        fprintf(output, "|\n");
    }

     fprintf(output, "| ------ ");
    for(i = data->min_score; i <= data->max_score; i++)
        fprintf(output, "| ----------------- ");
    fprintf(output, "|\n");

    fprintf(output, "| % 6s ", "Total");
    for(j = data->min_score; j <= data->max_score; j++){
        fprintf(output, "| % 8d (%5.2f%%) ", total[j],
               (total[j] * 100.0)/tt);
    }
    fprintf(output, "|\n");

    fprintf(output, "\n");


}


void draw_txt_length(FILE* output, sequence_data * data){
    int i;

    fprintf(output, "### Length Distribution \n\n");
    fprintf(output, "|   Base |     Terminating     |      Cumulative      |\n");
    fprintf(output, "| ------ | ------------------- | -------------------- |\n");

    for(i = 0; i < data->max_length-1; i++){
        uint64_t term = data->bases[i].length_count-data->bases[i+1].length_count;
        fprintf(output, "| % 6d ", i);
        fprintf(output, "| % 10d (%6.2f%%) ", term, term*100.0/data->number_of_sequences);
        fprintf(output, "| % 10d (%6.2f%%) ", data->bases[i].length_count,
               data->bases[i].length_count*100.0/data->number_of_sequences);
        fprintf(output, "|\n");
    }

    fprintf(output, "| % 6d ", i);
    fprintf(output, "| % 10d (%6.2f%%) ", data->bases[i].length_count,
           data->bases[i].length_count*100.0/data->number_of_sequences);
    fprintf(output, "| % 10d (%6.2f%%) ", data->bases[i].length_count,
           data->bases[i].length_count*100.0/data->number_of_sequences);
    fprintf(output, "|\n");

}

void draw_txt_adapters(FILE* output, sequence_data * data){
    int i;

    fprintf(output, "### Adapter Distribution \n\n");
    fprintf(output, "|   Base |        Count        |\n");
    fprintf(output, "| ------ | ------------------- |\n");

    for(i = 0; i < data->max_length; i++){
        fprintf(output, "| % 6d ", i);
        fprintf(output, "| % 10d (%6.2f%%) ", data->bases[i].kmer_count,
               data->bases[i].kmer_count*100.0/data->bases[i].length_count);
        fprintf(output, "|\n");
    }

}


void draw_txt(
    FILE * output,
    sequence_data* forward,
    sequence_data* reverse,
    char* name,
    int adapters_used){

    if(name)
        fprintf(output, "# %s\n\n", name);

    if(reverse)
        fprintf(output, "## Forward\n\n");

    draw_txt_content(output, forward);
    draw_txt_quality(output, forward);
    draw_txt_length(output, forward);

    if(adapters_used)
        draw_txt_adapters(output, forward);

    if(reverse){
        fprintf(output, "## Reverse\n\n");

        draw_txt_content(output, reverse);
        draw_txt_quality(output, reverse);
        draw_txt_length(output, reverse);

        if(adapters_used)
            draw_txt_adapters(output, reverse);
    }
}
