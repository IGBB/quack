#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#include "svg.h"
#include "seq.h"

#define GRAPH_WIDTH  450.0
#define GRAPH_HEIGHT 250.0
#define PERF_SIZE 100.0


char * point_string(int length, float* x, float* y){
    char * points, *cur_pos;
    size_t max = 0;
    int i = 0;

    /* Find lenght of each point string */
    for(i = 0; i < length; i++)
        max += snprintf(NULL, 0, " %.2f,%.2f ", x[i], y[i]);
    max++;

    /* Alloc enough memory and create string */
    points = malloc(sizeof(char[max]));
    cur_pos = points;
    for(i = 0; i < length; i++){
        max = sprintf(cur_pos, " %.2f,%.2f ", x[i], y[i]);
        cur_pos += max;
    }

    return points;
}

void draw_svg_content(sequence_data* data){

    int i, j;
    uint64_t len, y;
    char * points;

    float *perc[4], *x;
    for(i = 0; i < 4; i++)
        perc[i] = malloc(data->max_length * sizeof(float));
    x = malloc(data->max_length * sizeof(float));

    /* Calculate cumulative percentage of base content, in decending order so
     * they stack */
    for(i = 0; i < data->max_length; i++){
        x[i] = i + 0.5;

        y = 0;
        len = data->bases[i].length_count;
        for (j = 3; j >= 0; j--) {
            y += data->bases[i].content[j];
            perc[j][i] = y * 100.0 / len;
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

        points = point_string(data->max_length, x, perc[i]);

        svg_simple_tag("polyline", 3,
     /* Since coordinates for lines and rectangles don't work the same; set the
       first point of each line to start off graph. Then, add 0.5 to the x of each
       point. Finally, end the line off graph. */
                      svg_attr("points", "0,0 0,%6.2f %s %d,%6.2f %d,0",
                                perc[i][0],
                                points,
                                data->max_length,
                                perc[i][data->max_length-1],
                                data->max_length),
                       svg_attr("fill", "%s", ratio_colors[i]),
                       svg_attr("stroke", "%s", "none")
                       );

        free(points);
    }

    free(x);

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
    float perc, *x;

    x = malloc(sizeof(float[data->max_length]));

    svg_start_tag("g", 1,
                  svg_attr("transform", "translate(0, %f) scale(%f %f)",
                           GRAPH_HEIGHT,
                           GRAPH_WIDTH/(data->max_length+1),
                           -1.0 * GRAPH_HEIGHT/data->max_score)
                  );

   /* Define background for heatmap. Must be in descending order or highest score
     background will overwrite all other backgrounds */
    #define score_back(score, color)                                \
    svg_simple_tag("rect", 6,                                       \
                   svg_attr("x",      "%d", 0),                     \
                   svg_attr("y",      "%d", 0),                     \
                   svg_attr("width",  "%d", data->max_length),      \
                   svg_attr("height", "%d", score),                 \
                   svg_attr("stroke", "%s", "none"),                \
                   svg_attr("fill",   "%s", color)                  \
                   )
    score_back(data->max_score+1, "#ccebc5"); // green
    score_back(28,              "#ffffcc"); // yellow
    score_back(20,              "#fbb4ae"); // red


    /* Draw each heatmap */
    for(i = 0; i < data->max_length; i++){
        x[i] = i + 0.5;
        for(j = data->min_score; j <= data->max_score; j++){
            if(data->bases[i].scores[j] == 0)
                continue;
            perc =  ((float) data->bases[i].scores[j] )/data->bases[i].length_count;
            svg_simple_tag("rect", 8,
                       svg_attr("x",      "%d", i),
                       svg_attr("y",      "%d", j),
                       svg_attr("fill-opacity", "%f", perc),
                       svg_attr("width",  "%d", 1),
                       svg_attr("height", "%d", 1),
                       svg_attr("stroke", "%s", "none"),
                       svg_attr("stroke-width", "%d", 0),
                       svg_attr("fill",   "%s", "black")
                       );
        }
    }


    char* points = point_string(data->max_length, x, data->avg_score);

    svg_simple_tag("polyline", 5,
                   /* Since coordinates for lines and rectangles don't work the same; set the
                      first point of each line to start off graph. Then, add 0.5 to the x of each
                      point. Finally, end the line off graph. */
                   svg_attr("points", "0,%.2f %s %d,%.2f ",
                            data->avg_score[0],
                            points,
                            data->max_length,
                            data->avg_score[data->max_length-1]),
                   svg_attr("fill", "%s", "none"),
                   svg_attr("stroke", "%s", "black"),
                   /* The stroke width needs to be the inverse of the height
                    * scale to keep it at a constant thickness for any length of
                    * sequence*/
                   svg_attr("stroke-width", "%f", 2.0 * data->max_score/GRAPH_HEIGHT ),
                   svg_attr("stroke-opacity", "%f", 0.6 )
    );

    free(points);
    free(x);


    svg_end_tag("g");

    /* Add graph label */
    svg_start_tag("text", 5,
                  svg_attr("y",           "%f", GRAPH_HEIGHT-5),
                  svg_attr("fill",        "%s", "#AAA"),
                  svg_attr("x",           "%d", 5),
                  svg_attr("font-family", "%s", "sans-serif"),
                  svg_attr("font-size",   "%s", "15px")
                  );
    printf("%s\n", "Per Base Sequence Quality");
    svg_end_tag("text");


}


void draw_svg_length(sequence_data * data){
    int i;
    char * points;
    float *y, *x;
    x = malloc(data->max_length * sizeof(float));
    y = malloc(data->max_length * sizeof(float));

    for(i = 0; i < data->max_length; i++){
        x[i] = i + 0.5;
        y[i] = data->bases[i].length_count*100.0/data->number_of_sequences;
    }

    points = point_string(data->max_length, x, y);

    svg_start_tag("g", 1,
                  svg_attr("transform", "scale(%f %d)", GRAPH_WIDTH/data->max_length, 1)
                  );

    /* Set background color */
    svg_simple_tag("rect", 3,
                   svg_attr("width",  "%d", data->max_length),
                   svg_attr("height", "%d", 100),
                   svg_attr("fill", "%s", "#CCC")
                   );

    /* Draw length distribution */

    svg_simple_tag("polyline", 3,
                   /* Since coordinates for lines and rectangles don't work the same; set the
                      first point of each line to start off graph. Then, add 0.5 to the x of each
                      point. Finally, end the line off graph. */
                   svg_attr("points", "0,0 0,%6.2f %s %d,%6.2f %d,0",
                            y[0],
                            points,
                            data->max_length,
                            y[data->max_length-1],
                            data->max_length),
                   svg_attr("fill", "%s", "steelblue"),
                   svg_attr("stroke", "%s", "none")
    );

    free(points);
    free(x);
    free(y);

    svg_end_tag("g");

    /* Add graph label */
    svg_start_tag("text", 5,
                  svg_attr("y",           "%d", 95),
                  svg_attr("fill",        "%s", "#AAA"),
                  svg_attr("x",           "%d", 5),
                  svg_attr("font-family", "%s", "sans-serif"),
                  svg_attr("font-size",   "%s", "15px")
                  );
    printf("%s\n", "Length Distibution");
    svg_end_tag("text");


}

void draw_svg_adapter(sequence_data * data){
    int i;
    char * points;
    float *y, *x;
    x = malloc(data->max_length * sizeof(float));
    y = malloc(data->max_length * sizeof(float));

    for(i = 0; i < data->max_length; i++){
        x[i] = i + 0.5;
        y[i] = data->bases[i].kmer_count*100.0/data->number_of_sequences;
    }

    points = point_string(data->max_length, x, y);

    svg_start_tag("g", 1,
                  svg_attr("transform", "scale(%f %d)", GRAPH_WIDTH/data->max_length, 1)
                  );

    /* Set background color */
    svg_simple_tag("rect", 3,
                   svg_attr("width",  "%d", data->max_length),
                   svg_attr("height", "%d", 100),
                   svg_attr("fill", "%s", "#CCC")
                   );

    /* Draw length distribution */

    svg_simple_tag("polyline", 3,
                   /* Since coordinates for lines and rectangles don't work the same; set the
                      first point of each line to start off graph. Then, add 0.5 to the x of each
                      point. Finally, end the line off graph. */
                   svg_attr("points", "0,0 0,%6.2f %s %d,%6.2f %d,0",
                            y[0],
                            points,
                            data->max_length,
                            y[data->max_length-1],
                            data->max_length),
                   svg_attr("fill", "%s", "steelblue"),
                   svg_attr("stroke", "%s", "none")
    );

    free(points);
    free(x);
    free(y);

    svg_end_tag("g");

    /* Add graph label */
    svg_start_tag("text", 5,
                  svg_attr("y",           "%d", 95),
                  svg_attr("fill",        "%s", "#AAA"),
                  svg_attr("x",           "%d", 5),
                  svg_attr("font-family", "%s", "sans-serif"),
                  svg_attr("font-size",   "%s", "15px")
                  );
    printf("%s\n", "Adapter Distibution");
    svg_end_tag("text");


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
                      svg_attr("transform", "translate(%d %d)", 110, 10)
    );
    draw_svg_content(forward);
    svg_end_tag("g");

    svg_start_tag("g", 1,
                      svg_attr("transform", "translate(%d %d)", 110, 130)
    );
    draw_svg_quality(forward);
    svg_end_tag("g");

    svg_start_tag("g", 1,
                      svg_attr("transform", "translate(%d %d)", 110, 400)
    );
    draw_svg_length(forward);
    svg_end_tag("g");

     svg_start_tag("g", 1,
                      svg_attr("transform", "translate(%d %d)", 110, 520)
    );
    draw_svg_adapter(forward);
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
