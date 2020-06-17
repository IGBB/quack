#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "svg.h"
#include "seq.h"

#define GRAPH_WIDTH 400.0

void draw_svg_content(sequence_data* data){

    int i;
    uint64_t x, y, len;

    float *perc[4];
    for(i = 0; i < 4; i++)
        perc[i] = malloc(data->max_length * sizeof(float));


    /* Calculate cumulative percentage of base content, in decending order so
     * they stack */
    for(x = 0; x < data->max_length; x++){
        y = 0;
        len = data->bases[x].length_count;
        for (i = 3; i >= 0; i--) {
            y += data->bases[x].content[i];
            perc[i][x] = y * 100.0 / len;
        }
    }

    /* Allocate 16 characters per base (A,C,T,G) per point.
       1 = ' '
       4 = max length is capped in kb range, will start compressing if larger
       2 = '.5' added to point
       1 = ','
       6 = read percentage length
       1 = ' '
       1 = padding for miscalculation
    */
    size_t ratio_points_length = 16*data->max_length;
    char * ratio_points_start[4];
    char * ratio_points[4];
    for(i = 0; i < 4; i++){
        ratio_points[i] = malloc(ratio_points_length);
        ratio_points_start[i] = ratio_points[i];
    }

    /* Add points to string, adjusting ration_points[i] to the end of the string
     * snprintf restricts the length of the string to 16. should stop
     * overflows*/
    for (x = 0; x < data->max_length; x++) {
        for(i = 0; i < 4; i++){
            len = snprintf(ratio_points[i], 16, " %4d.5,%6.2f ", x, perc[i][x]);
            ratio_points[i] += len;
        }
    }


    svg_start_tag("g", 1,
                  svg_attr("transform", "scale(%f %d)", GRAPH_WIDTH/data->max_length, 1)
                  );

    /* Set background color */
    svg_simple_tag("rect", 3,
                   svg_attr("width",  "%d", data->max_length),
                   svg_attr("height", "%d", 100),
                   svg_attr("fill", "%s", "#CCC")
                   );

    /* Draw each distribution */
    char *ratio_labels[4] = {"%A", "%T", "%C", "%G"};
    char *ratio_colors[4] = {"#648964", "#89bc89", "#84accf", "#5d7992"};
    for(i = 0; i < 4; i++){
        svg_simple_tag("polyline", 3,
     /* Since coordinates for lines and rectangles don't work the same; set the
       first point of each line to start off graph. Then, add 0.5 to the x of each
       point. Finally, end the line off graph. */
                      svg_attr("points", "0,0 0,%6.2f %s %d,%6.2f %d,0",
                                perc[i][0],
                                ratio_points_start[i],
                                data->max_length,
                                perc[i][data->max_length-1],
                                data->max_length),
                       svg_attr("fill", "%s", ratio_colors[i]),
                       svg_attr("stroke", "%s", "none")
                       );

        free(ratio_points_start[i]);
        free(perc[i]);
    }

    svg_end_tag("g");

    /* Add graph label */
    svg_start_tag("text", 5,
                  svg_attr("y",           "%d", 95),
                  svg_attr("fill",        "%s", "#CCC"),
                  svg_attr("x",           "%d", 5),
                  svg_attr("font-family", "%s", "sans-serif"),
                  svg_attr("font-size",   "%s", "15px")
                  );
    printf("%s\n", "Base Content Percentage");
    svg_end_tag("text");


}


void draw_svg_quality(sequence_data* data){
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


void draw_svg_length(sequence_data * data){
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

void draw_svg_adapters(sequence_data * data){
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


void draw_svg(sequence_data* forward,
              sequence_data* reverse,
              char* name,
              int adapters_used){

    /* Set width */
    int width = 615;
    if(reverse)
        width = 1195;

    /* Set height */
    int height = 510;
    if(adapters_used)
        height = 610;

    // Start svg
    svg_start_tag("svg", 5,
                  svg_attr("width",   "%d", width),
                  svg_attr("height",  "%d", height),
                  svg_attr("viewBox", "%d %d %d %d", 0, 0, width, height),
                  svg_attr("xmlns",       "%s", "http://www.w3.org/2000/svg"),
                  svg_attr("xmlns:xlink", "%s", "http://www.w3.org/1999/xlink")
                  );

    /* If name is given, add to middle of viewBox (half of width + min-x of viewbox) */
    if(name != NULL){
        height += 30;
        svg_start_tag("text", 6,
                      svg_attr("x", "%d", (width/2)),
                      svg_attr("y", "%d", 30),
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("font-size",   "%s", "30px"),
                      svg_attr("fill",        "%s", "black"));
        printf("%s\n", name);
        svg_end_tag("text");

        svg_start_tag("g", 1,
                      svg_attr("transform", "translate(%d %d)", 0, 30)
        );

    }

    svg_start_tag("g", 1,
                      svg_attr("transform", "translate(%d %d)", 50, 10)
    );
    draw_svg_content(forward);
    svg_end_tag("g");

    svg_start_tag("g", 1,
                      svg_attr("transform", "translate(%d %d)", 50, 10)
    );
    draw_svg_content(forward);
    svg_end_tag("g");


    /* draw(data, 0, adapters); */
    /* free(data); */

    /* if(paired){ */
    /*   data = read_fastq(arguments.reverse, kmers); */
    /*   transformed_data = transform(data); */
    /*   draw(data, 1, adapters); */
    /*   free(data); */
    /* } */

    if(name != NULL) svg_end_tag("g");

    svg_end_tag("svg");

}
