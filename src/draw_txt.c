#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "seq.h"


void draw_txt_content(sequence_data* data){

    int i, j;

    printf("### Per Base Content \n\n");
    printf("|   Base |          A          |          T          |          C          |          G          |\n");
    printf("| ------ | ------------------- | ------------------- | ------------------- | ------------------- |\n");

    for(i = 0; i < data->max_length; i++){
        printf("| % 6d ", i);
        for(j = 0; j < 4; j++){
            printf("| % 10d (%5.2f%%) ", data->bases[i].content[j],
                   (data->bases[i].content[j] * 100.0)/data->bases[i].length_count);
        }
        printf("|\n");
    }
    printf("\n");
}


void draw_txt_quality(sequence_data* data){
    int i, j;
    uint64_t total[91] = {0}, tt = 0;


    printf("### Per Base Sequence Quality \n\n");
    printf("|   Base ");
    for(i = data->min_score; i <= data->max_score; i++)
        printf("| %17d ", i);
    printf("|\n");

    printf("| ------ ");
    for(i = data->min_score; i <= data->max_score; i++)
        printf("| ----------------- ");
    printf("|\n");

    for(i = 0; i < data->max_length; i++){
        printf("| % 6d ", i);
        for(j = data->min_score; j <= data->max_score; j++){
            tt += data->bases[i].scores[j];
            total[j] += data->bases[i].scores[j];
            printf("| % 8d (%5.2f%%) ", data->bases[i].scores[j],
                   (data->bases[i].scores[j] * 100.0)/data->bases[i].length_count);
        }
        printf("|\n");
    }

     printf("| ------ ");
    for(i = data->min_score; i <= data->max_score; i++)
        printf("| ----------------- ");
    printf("|\n");

    printf("| % 6s ", "Total");
    for(j = data->min_score; j <= data->max_score; j++){
        printf("| % 8d (%5.2f%%) ", total[j],
               (total[j] * 100.0)/tt);
    }
    printf("|\n");

    printf("\n");


}


void draw_txt_length(sequence_data * data){
    int i;

    printf("### Length Distribution \n\n");
    printf("|   Base |     Terminating     |      Cumulative      |\n");
    printf("| ------ | ------------------- | -------------------- |\n");

    for(i = 0; i < data->max_length-1; i++){
        uint64_t term = data->bases[i].length_count-data->bases[i+1].length_count;
        printf("| % 6d ", i);
        printf("| % 10d (%6.2f%%) ", term, term*100.0/data->number_of_sequences);
        printf("| % 10d (%6.2f%%) ", data->bases[i].length_count,
               data->bases[i].length_count*100.0/data->number_of_sequences);
        printf("|\n");
    }

    printf("| % 6d ", i);
    printf("| % 10d (%6.2f%%) ", data->bases[i].length_count,
           data->bases[i].length_count*100.0/data->number_of_sequences);
    printf("| % 10d (%6.2f%%) ", data->bases[i].length_count,
           data->bases[i].length_count*100.0/data->number_of_sequences);
    printf("|\n");

}

void draw_txt_adapters(sequence_data * data){
    int i;

    printf("### Adapter Distribution \n\n");
    printf("|   Base |        Count        |\n");
    printf("| ------ | ------------------- |\n");

    for(i = 0; i < data->max_length; i++){
        printf("| % 6d ", i);
        printf("| % 10d (%6.2f%%) ", data->bases[i].kmer_count,
               data->bases[i].kmer_count*100.0/data->bases[i].length_count);
        printf("|\n");
    }

}


void draw_txt(sequence_data* forward,
              sequence_data* reverse,
              char* name,
              int adapters_used){

    if(name)
        printf("# %s\n\n", name);

    if(reverse)
        printf("## Forward\n\n");

    draw_txt_content(forward);
    draw_txt_quality(forward);
    draw_txt_length(forward);

    if(adapters_used)
        draw_txt_adapters(forward);

    if(reverse){
        printf("## Reverse\n\n");

        draw_txt_content(reverse);
        draw_txt_quality(reverse);
        draw_txt_length(reverse);

        if(adapters_used)
            draw_txt_adapters(reverse);
    }
}
