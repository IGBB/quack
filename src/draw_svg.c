#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <string.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "svg.h"
#include "seq.h"

#define GRAPH_WIDTH  450.0
#define GRAPH_HEIGHT 250.0
#define PERF_SIZE 100.0
#define GRAPH_PAD 5

char *ratio_labels[4] = {"%A", "%T", "%C", "%G"};
char *ratio_colors[4] = {"#648964", "#89bc89", "#84accf", "#5d7992"};

#define PI 3.14159265

FILE * output;

typedef struct {
    float x,y;
} fpair_t;

/******************************************************************************/
/*            Generic find max function                                       */
/******************************************************************************/
float get_max(sequence_data* f, sequence_data* r,
              int (*stop)(sequence_data*),
              float (*accessor)(sequence_data*, int)){
    float max = 0.0;
    int i,n;

    for( n = 0; n < 2; n++ ){
        sequence_data* data = (n == 0)?f:r;
        if(data == NULL) continue;

        for(i = 0; i < (*stop)(data); i++){
            float current = (*accessor)(data, i);
            if(current > max) max = current;
        }
    }

    // Adjust max to nearest ten
    max = ceil(max / 10.0) * 10.0;

    return max;
}

/******************************************************************************/
/*            Iterator Stops function                                        */
/******************************************************************************/
inline int length_stop(sequence_data* data){
    return data->max_length;
}
inline int score_stop(sequence_data* data ){
    return data->max_score+1;
}
inline int sequence_stop(sequence_data* data ){
    return data->windows+1;
}

/*  Max length functions */
inline float length_accessor(sequence_data* data, int i){
    return data->bases[i].length_count*100.0/data->number_of_sequences;
}
float get_max_length(sequence_data* f, sequence_data* r){
    return get_max(f,r,length_stop,length_accessor);
}

/*  Max adapter functions */
inline float adapter_accessor(sequence_data* data, int i){
    float ret = data->bases[i].kmer_count*100.0/data->number_of_sequences;
    return ret;
}
float get_max_adapter(sequence_data* f, sequence_data* r){
    float max = get_max(f,r,length_stop,adapter_accessor);
    /*  Adjust max to sane level */
    if(max < 20.0) max = 20.0;
    return max;
}

/*  Max score functions */
inline float score_accessor(sequence_data* data, int i){
    unsigned int j;
    float sum = 0;
    for(j = 0; j < data->max_length; j++){
        sum += data->bases[j].scores[i];
    }
    return sum*100.0/((float)data->number_of_sequences * data->max_length);
}
float get_max_score(sequence_data* f, sequence_data* r){
    return get_max(f,r,score_stop,score_accessor);
}

inline float content_accessor(sequence_data* data, int i){
    int j;
    uint64_t max = 0;
    for (j = 3; j >= 0; j--) {
        if(max < data->bases[i].content[j])
            max =  data->bases[i].content[j];
    }
    return max * 100.0 / data->bases[i].length_count;
}
float get_max_content(sequence_data* f, sequence_data* r){
    return get_max(f,r,length_stop,content_accessor);
}

/*  Max saturation functions */
inline float saturation_accessor(sequence_data* data, int i){
    float ret = data->new_kmers[i];
    return ret;
}
float get_max_saturation(sequence_data* f, sequence_data* r){
    float max = get_max(f,r,sequence_stop,saturation_accessor);
    /*  Adjust max to sane level */
    if(max < 20.0) max = 20.0;
    return max;
}


char * point_string(int length, fpair_t* p){
    char * points, *cur_pos;
    size_t max = 0;
    int i = 0;

    /* Find length of each point string */
    for(i = 0; i < length; i++)
        max += snprintf(NULL, 0, " %.2f,%.2f ", p[i].x, p[i].y);
    max++;

    /* Alloc enough memory and create string */
    points = malloc(sizeof(char[max]));
    cur_pos = points;
    for(i = 0; i < length; i++){
        max = sprintf(cur_pos, " %.2f,%.2f ", p[i].x, p[i].y);
        cur_pos += max;
    }

    return points;
}


#define svg_start_label(posx, posy)                                     \
    svg_start_tag(output, "text", 5,                                            \
                  svg_attr("x",           "%f", (float)(posx)),         \
                  svg_attr("y",           "%f", (float)(posy)),         \
                  svg_attr("fill",        "%s", "#AAA"),                \
                  svg_attr("font-family", "%s", "sans-serif"),          \
                  svg_attr("font-size",   "%s", "15px")                 \
    )

#define svg_end_label() svg_end_tag(output, "text")


int calc_step(float min, float max, int breaks){
    float length = max - min;
    /* Find a tick step that is
     * - at least 5
     * - divisible by 5,
     * - creates 5 ticks */
    int step = length/breaks;
    step = (step /5 )*5;
    if (step < 5) step = 5;

    return step;
}

void xaxis_length_ticks( fpair_t size, fpair_t scale){
    int step = calc_step(0,size.x,10);
    int i;

    svg_simple_tag(output, "line", 7,
                   svg_attr("y1", "%d", 0),
                   svg_attr("y2", "%d", 0),
                   svg_attr("x1", "%d", 0),
                   svg_attr("x2", "%f", size.x),
                   svg_attr("stroke", "%s", "black"),
                   svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                   svg_attr("stroke-width", "%f", 0.5));


    for (i = step; i <= size.x - step; i += step){
        svg_simple_tag(output, "line", 7,
                       svg_attr("x1", "%f", i-0.5),
                       svg_attr("y1", "%f", -2.0 / scale.y),
                       svg_attr("x2", "%f", i-0.5),
                       svg_attr("y2", "%f", 2.0 / scale.y),
                       svg_attr("stroke", "%s", "black"),
                       svg_attr("stroke-opacity", "%f", 0.5),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke")
        );
    }

    svg_simple_tag(output, "line", 7,
                       svg_attr("x1", "%f", size.x-0.5),
                       svg_attr("y1", "%f", -2.0 / scale.y),
                       svg_attr("x2", "%f", size.x-0.5),
                       svg_attr("y2", "%f", 2.0 / scale.y),
                       svg_attr("stroke", "%s", "black"),
                       svg_attr("stroke-opacity", "%f", 0.5),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke")
        );
}

void yaxis_score_labels( fpair_t size, fpair_t scale){

    int i, step = calc_step(0, size.y, 10);

    /*  loop through each step. excluded last step if too close to max score */
    for (i = step; i < size.y - (step/2); i += step){

        svg_start_tag(output, "text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %f) scale(%f %f)",
                               -4.0, i+0.5,
                               1.0/scale.x, 1.0/scale.y));

        fprintf(output, "%d", i);
        svg_end_tag(output, "text");


    }
    /* add max score label */
        svg_start_tag(output, "text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %f) scale(%f %f)",
                               -4.0, size.y-0.5,
                               1.0/scale.x, 1.0/scale.y));

        fprintf(output, "%d", (int)size.y-1);
        svg_end_tag(output, "text");

    svg_simple_tag(output, "line", 7,
                   svg_attr("x1", "%d", 0),
                   svg_attr("x2", "%d", 0),
                   svg_attr("y1", "%d", 0),
                   svg_attr("y2", "%f", size.y),
                   svg_attr("stroke", "%s", "black"),
                   svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                   svg_attr("stroke-width", "%f", 0.5));

}


void xaxis_percent_labels( fpair_t size, fpair_t scale){
    int i, step = calc_step(0, size.x, 4);

    for(i = step; i < size.x; i+=step){
        svg_simple_tag(output, "line", 7,
                       svg_attr("x1", "%d", i),
                       svg_attr("x2", "%d", i),
                       svg_attr("y1", "%d", 0),
                       svg_attr("y2", "%f", size.y),
                       svg_attr("stroke", "%s", "white"),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                       svg_attr("stroke-width", "%d", 1)

        );
        svg_start_tag(output, "text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "start"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %f) scale(%f %f) rotate(-45)",
                               (float)i,  size.y + 1,
                               1/scale.x, 1/scale.y ));
        fprintf(output, "%d%%", i);
        svg_end_tag(output, "text");
    }

}


void yaxis_percent_labels( fpair_t size, fpair_t scale){
    int i, step = calc_step(0,size.y,4);

    for(i = step; i < size.y; i+=step){
        svg_simple_tag(output, "line", 7,
                       svg_attr("y1", "%d", i),
                       svg_attr("y2", "%d", i),
                       svg_attr("x1", "%d", 0),
                       svg_attr("x2", "%f", size.x),
                       svg_attr("stroke", "%s", "white"),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                       svg_attr("stroke-width", "%f", 1.0)

                       );
        svg_start_tag(output, "text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "start"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %f) scale(%f %f)",
                               size.x + 1, (float)i,
                               1/scale.x, 1/scale.y ));
        fprintf(output,"%d%%", i);
        svg_end_tag(output, "text");
    }

}

void yaxis_saturation_labels( fpair_t size, fpair_t scale){
    int i, step = calc_step(0,size.y,4);

    for(i = step; i < size.y; i+=step){
        svg_simple_tag(output, "line", 7,
                       svg_attr("y1", "%d", i),
                       svg_attr("y2", "%d", i),
                       svg_attr("x1", "%d", 0),
                       svg_attr("x2", "%f", size.x),
                       svg_attr("stroke", "%s", "white"),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                       svg_attr("stroke-width", "%f", 1.0)

                       );
        svg_start_tag(output, "text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "start"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %f) scale(%f %f)",
                               size.x + 1, (float)i,
                               1/scale.x, 1/scale.y ));
        fprintf(output,"%d", i);
        svg_end_tag(output, "text");
    }

}

void draw_svg_graph(sequence_data* data,
                     fpair_t original, fpair_t final,
                     fpair_t flip, fpair_t translate,
                     void (*graph_points)(sequence_data*),
                     void (*xaxis)(fpair_t, fpair_t),
                     void (*yaxis)(fpair_t, fpair_t),
                     char* label){

    fpair_t scale = (fpair_t){final.x/original.x, final.y/original.y};
    fpair_t scale_translate = (fpair_t){0,0};

    if(flip.x){
        scale.x *= -1;
        scale_translate.x = final.x;
    }
    if(flip.y){
        scale.y *= -1;
        scale_translate.y = final.y;
    }

    svg_start_tag(output, "g", 1, svg_attr("transform", "translate(%f %f)", translate.x, translate.y));
    svg_start_tag(output, "g", 1, svg_attr("transform", "translate(%f %f) scale(%f %f)",
                                   scale_translate.x, scale_translate.y,
                                   scale.x, scale.y
                  )
    );



    (*graph_points)(data);

    if(xaxis != NULL) (*xaxis)(original, scale);
    if(yaxis != NULL) (*yaxis)(original, scale);

    svg_end_tag(output, "g");

    /* Add graph label */
    svg_start_label(GRAPH_PAD, final.y - GRAPH_PAD);
    fprintf(output,"%s\n", label);
    svg_end_label();

    svg_end_tag(output, "g");

}

void draw_histo_length(sequence_data* data){
    unsigned int i;
    float y;
    for(i = 0; i < data->max_length-1; i++){
        y = (data->bases[i].length_count - data->bases[i+1].length_count);
        y = y*100.0/data->number_of_sequences;

        svg_simple_tag(output, "rect", 5,
                       svg_attr("x", "%d", i),
                       svg_attr("height", "%f", y),
                       svg_attr("width", "%f", 1.1),
                       svg_attr("y", "%d", 0),
                       svg_attr("fill", "%s", "steelblue")
                       );


    }

    /*  Draw last bar */
    y = (data->bases[i].length_count)*100.0/data->number_of_sequences;

    svg_simple_tag(output, "rect", 5,
                   svg_attr("x", "%d", i),
                   svg_attr("height", "%f", y),
                   svg_attr("width", "%f", 1.1),
                   svg_attr("y", "%d", 0),
                   svg_attr("fill", "%s", "steelblue")
    );

}

void draw_svg_length(sequence_data * data, fpair_t translate, float max){
    unsigned int i;
    fpair_t original, final, flip;




    original = (fpair_t){data->max_length, max};
    final    = (fpair_t){GRAPH_WIDTH, PERF_SIZE};
    flip     = (fpair_t){0,1};

    draw_svg_graph(data, original, final, flip,
                    translate, draw_histo_length,
                    xaxis_length_ticks, yaxis_percent_labels, "Length Distribution");

    for(i = 0; i < data->max_length - 1; i++){
        /* Stop if current point is different than next */
        if(data->bases[i].length_count != data->bases[i+1].length_count)
            break;
    }

    if(i >= data->max_length -1){
        svg_simple_tag(output, "rect", 6,
                       svg_attr("x", "%f", translate.x),
                       svg_attr("y", "%f", translate.y),
                       svg_attr("width", "%f", final.x),
                       svg_attr("height", "%f", final.y),
                       svg_attr("fill", "%s", "#CCC"),
                       svg_attr("opacity", "%f", 0.5)
                       );

        svg_start_tag(output, "text", 7,
                      svg_attr("x", "%f", translate.x + final.x/2),
                      svg_attr("y", "%f", translate.y + final.y/2),
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("font-variant", "%s", "small-caps"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("font-size",   "%s", "12px"),
                      svg_attr("fill",        "%s", "black"));
        fprintf(output, "All reads are %" PRIu64 "-bp", data->max_length);
        svg_end_tag(output, "text");
    }

}

void draw_line_content(sequence_data* data){
    uint64_t i;
    int j;
    char * points_string;

    fpair_t *points = calloc(data->max_length, sizeof(fpair_t));

    /* Calculate cumulative percentage of base content */
    for(j = 0; j < 4; j++){
        for(i = 0; i < data->max_length; i++){
            points[i].x = i +0.5;
            points[i].y = data->bases[i].content[j] * 100.0 / data->bases[i].length_count;
        }
       
        points_string = point_string(data->max_length, points);

        svg_simple_tag(output, "polyline", 5,
                      svg_attr("points", "%s", points_string),
                       svg_attr("fill", "%s", "none"),
                       svg_attr("stroke", "%s", ratio_colors[j]),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                       svg_attr("stroke-width", "%f", 2.0)
        );

        free(points_string);

    }

    free(points);
}

void xaxis_length_dotted(fpair_t size, fpair_t scale){
    int i, step = calc_step(0, size.x, 10);
    for (i = step; i <= size.x - step; i += step){
        svg_simple_tag(output, "line", 8,
                       svg_attr("x1", "%f", i-0.5),
                       svg_attr("y1", "%d", -1),
                       svg_attr("x2", "%f", i-0.5),
                       svg_attr("y2", "%f", size.y),
                       svg_attr("stroke", "%s", "black"),
                       svg_attr("stroke-opacity", "%f", 0.25),
                       svg_attr("stroke-dasharray", "%d %d", 1, 9),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke")
                       );
    }
}
void xaxis_length_dotted_labels(fpair_t size, fpair_t scale){
    int i, step = calc_step(0, size.x, 10);

    xaxis_length_dotted(size, scale);

    for (i = step; i <= size.x - step; i += step){
        svg_start_tag(output, "text", 7,
                      svg_attr("font-family", "sans-serif"),
                      svg_attr("text-anchor", "middle"),
                      svg_attr("dominant-baseline", "middle"),
                      svg_attr("font-size", "8px"),
                      svg_attr("fill", "black"),
                      svg_attr("fill-opacity", "0.5"),
                      svg_attr("transform", "translate(%f %d) scale(%f %f)",
                               i -0.5, -1,
                               1/scale.x, 1/scale.y));

        fprintf(output, "%d", i);
        svg_end_tag(output, "text");
    }
        svg_start_tag(output, "text", 7,
                      svg_attr("font-family", "sans-serif"),
                      svg_attr("text-anchor", "middle"),
                      svg_attr("dominant-baseline", "middle"),
                      svg_attr("font-size", "8px"),
                      svg_attr("fill", "black"),
                      svg_attr("fill-opacity", "0.5"),
                      svg_attr("transform", "translate(%f %d) scale(%f %f)",
                               size.x -0.5, -1,
                               1/scale.x, 1/scale.y));

        fprintf(output, "%d", (int) size.x);
        svg_end_tag(output, "text");
}

void yaxis_content(fpair_t size, fpair_t scale){
    int i;
    for (i = 10; i < size.y; i += 10){
                svg_simple_tag(output, "line", 7,
                       svg_attr("y1", "%d", i),
                       svg_attr("y2", "%d", i),
                       svg_attr("x1", "%d", 0),
                       svg_attr("x2", "%f", size.x),
                       svg_attr("stroke", "%s", "#AAA"),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                       svg_attr("stroke-width", "%f", 0.5));

                svg_start_tag(output, "text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "start"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %d) scale(%f %f)",
                               size.x, i,
                           1/scale.x, 1/scale.y));
        fprintf(output, "%d%%", i);
        svg_end_tag(output, "text");
    }

    for(i = 0; i<4; i++){
        svg_start_tag(output, "g", 1,
                      svg_attr("transform", "translate(%f %f) scale(%f %f)",
                               -20.0 / scale.x , size.y * (i + 1) / 6.0,
                               1/scale.x, 1/scale.y));
        svg_simple_tag(output, "rect", 5,
                       svg_attr("y", "%d", -15),
                       svg_attr("x", "%d", 15),
                       svg_attr("width", "%d", 4),
                       svg_attr("height", "%d", 16),
                       svg_attr("fill",        "%s", ratio_colors[i]));

        svg_start_tag(output, "text", 6,
                      svg_attr("x", "%d", 8),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("font-size",   "%s", "16px")

        );
        fprintf(output, "%c", rev_lookup[i] );
        svg_end_tag(output, "text");
        svg_end_tag(output, "g");
    }


}
void draw_svg_content(sequence_data* data, fpair_t translate, float max){
    fpair_t original, final, flip;


    original = (fpair_t){data->max_length, max};
    final    = (fpair_t){GRAPH_WIDTH, PERF_SIZE};
    flip     = (fpair_t){0, 1};

    draw_svg_graph(data, original, final, flip,
                    translate, draw_line_content,
                   xaxis_length_dotted, yaxis_content,
                   "Base Content Percentage");

}

void draw_heatmap_quality(sequence_data* data){
    uint64_t i;
    int j;
    /* Draw heatmap background */
    int score_delineation[4] = {0,       20,       28, data->max_score+1};
    char* score_color[4]      = {"", "#d73027","#ffffbf",         "#1a9850"};

    for(i = 1; i <= 3; i++){
        int height = score_delineation[i] - score_delineation[i-1];
        svg_simple_tag(output, "rect", 5,
                       svg_attr("y", "%d", score_delineation[i-1]),
                       svg_attr("width", "%d", data->max_length),
                       svg_attr("height", "%d", height),
                       svg_attr("fill", "%s", score_color[i]),
                       svg_attr("fill-opacity", "%f", 0.1)
                       );
    }

    /* Draw each heatmap */

    fpair_t * p = calloc(data->max_length, sizeof(fpair_t));
    for(i = 0; i < data->max_length; i++){
        p[i].x = i + 0.5;
        p[i].y = data->avg_score[i];

        for(j = data->min_score; j <= data->max_score; j++){
            if(data->bases[i].scores[j] == 0)
                continue;
            float perc =  ((float) data->bases[i].scores[j] )/data->bases[i].length_count;
            svg_simple_tag(output, "rect", 8,
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


    char* points = point_string(data->max_length, p);

    svg_simple_tag(output, "polyline", 6,
                   /* Since coordinates for lines and rectangles don't work the same; set the
                      first point of each line to start off graph. Then, add 0.5 to the x of each
                      point. Finally, end the line off graph. */
                   svg_attr("points", "%s", points),
                   svg_attr("fill", "%s", "none"),
                   svg_attr("stroke", "%s", "black"),
                   svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                   svg_attr("stroke-width", "%f", 1.0),
                   svg_attr("stroke-opacity", "%f", 0.6 )
    );

    free(points);
    free(p);
}

void draw_svg_quality(sequence_data* data, fpair_t translate){
    fpair_t original, final, flip;

    original = (fpair_t){data->max_length, data->max_score+1};
    final    = (fpair_t){GRAPH_WIDTH, GRAPH_HEIGHT};
    flip     = (fpair_t){0, 1};

    draw_svg_graph(data, original, final, flip,
                    translate, draw_heatmap_quality,
                   xaxis_length_dotted_labels, NULL,
                  "Per Base Sequence Quality");

}

void draw_histo_score(sequence_data * data){
    int i;

    for(i = 0; i <= data->max_score; i++){
        float y = score_accessor(data, i);
        svg_simple_tag(output, "rect",5,
                       svg_attr("y",      "%d", i),
                       svg_attr("width",  "%f", y),
                       svg_attr("x",      "%f", 0.0),
                       svg_attr("height", "%f", 1.0),
                       svg_attr("fill",   "%s", "steelblue")
                       );
    }

}

void draw_svg_score(sequence_data * data, fpair_t translate, int flipx, float max){
    fpair_t original, final, flip;


    original = (fpair_t){max, data->max_score+1};
    final    = (fpair_t){PERF_SIZE, GRAPH_HEIGHT};
    flip     = (fpair_t){!flipx, 1};

    draw_svg_graph(data, original, final, flip,
                    translate, draw_histo_score,
                   xaxis_percent_labels, yaxis_score_labels,
                   "<tspan dy=\"-15\">Score</tspan>"
                   "<tspan dy=\"15\" x=\"5\">Distribution</tspan>");

}

void draw_histo_adapter(sequence_data * data){
    unsigned int i;
    float y;

    for(i = 0; i < data->max_length; i++){
        y = data->bases[i].kmer_count*100.0/data->number_of_sequences;

        svg_simple_tag(output, "rect", 5,
                       svg_attr("x", "%d", i),
                       svg_attr("height", "%f", y),
                       svg_attr("width", "%f", 1.1),
                       svg_attr("y", "%d", 0),
                       svg_attr("fill", "%s", "steelblue")
                       );


    }

}

void draw_svg_adapter(sequence_data * data, fpair_t translate, float max){
    fpair_t original_size, final_size, flip;


    original_size = (fpair_t){data->max_length, max};
    final_size    = (fpair_t){GRAPH_WIDTH, PERF_SIZE};
    flip          = (fpair_t){0,1};

    draw_svg_graph(data,
                   original_size, final_size, flip, translate,
                   draw_histo_adapter, xaxis_length_ticks, yaxis_percent_labels,
                   "Adapter Distribution");


}

void draw_histo_saturation(sequence_data * data){
    unsigned int i;
    float y;

    for(i = 0; i < data->windows; i++){
        y = saturation_accessor(data, i);

        svg_simple_tag(output, "rect", 5,
                       svg_attr("x", "%d", i),
                       svg_attr("height", "%f", y),
                       svg_attr("width", "%f", 1.1),
                       svg_attr("y", "%d", 0),
                       svg_attr("fill", "%s", "steelblue")
                       );


    }
    printf("Done\n");

}

void draw_svg_saturation(sequence_data * data, fpair_t translate, float max){
    fpair_t original_size, final_size, flip;

    original_size = (fpair_t){data->windows, max};
    final_size    = (fpair_t){GRAPH_WIDTH, PERF_SIZE};
    flip          = (fpair_t){0,1};

    draw_svg_graph(data,
                   original_size, final_size, flip, translate,
                   draw_histo_saturation, xaxis_length_ticks, yaxis_saturation_labels,
                   "Novel K-mer Saturation");

}

void draw_svg(FILE * out,
              sequence_data* forward,
              sequence_data* reverse,
              char* name,
              int adapters_used,
              int saturation_curve){

    output = out;

    int i;
    sequence_data* list[2] = {forward, reverse};

    float max_length     = get_max_length(forward, reverse);
    float max_score      = get_max_score(forward, reverse);
    float max_adapter    = get_max_adapter(forward, reverse);
    float max_content    = get_max_content(forward, reverse);
    float max_saturation = get_max_saturation(forward, reverse);

    /* Set size */
    fpair_t size = {615, 510};
    if(reverse)
        size.x = 1195;

    if(adapters_used)
        size.y = 610;
    if(name != NULL)
        size.y += 30;

    if(saturation_curve)
        size.y += 100;

    // Start svg
    svg_start_tag(output, "svg", 5,
                  svg_attr("width",   "%f", size.x),
                  svg_attr("height",  "%f", size.y),
                  svg_attr("viewBox", "%f %f %f %f", 0.0, 0.0, size.x, size.y),
                  svg_attr("xmlns",       "%s", "http://www.w3.org/2000/svg"),
                  svg_attr("xmlns:xlink", "%s", "http://www.w3.org/1999/xlink")
                  );

    /* If name is given, add to middle of viewBox (half of width + min-x of viewbox) */
    if(name != NULL){
        svg_start_tag(output, "text", 6,
                      svg_attr("x", "%f", (size.x/2)),
                      svg_attr("y", "%d", 30),
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("font-size",   "%s", "30px"),
                      svg_attr("fill",        "%s", "black"));
        fprintf(output, "%s\n", name);
        svg_end_tag(output, "text");

        svg_start_tag(output, "g", 1,
                      svg_attr("transform", "translate(%d %d)", 0, 50)
        );

    }


    fpair_t translate = {GRAPH_PAD*2 + PERF_SIZE + 10, GRAPH_PAD};
    float score_locations[2] = {GRAPH_PAD, GRAPH_PAD*4 + GRAPH_WIDTH*2 + PERF_SIZE + 60};
    for(i = 0; i < 2; i++){
        sequence_data * data = list[i];
        if(data == NULL) continue;
       
        translate.y = GRAPH_PAD;

        svg_start_tag(output, "text", 6,
                      svg_attr("x", "%f", (GRAPH_WIDTH/2) + translate.x),
                      svg_attr("y", "%d", 0),
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("font-size",   "%s", "15px"),
                      svg_attr("fill",        "%s", "black"));
        fprintf(output, "%" PRId64 " reads encoded in phred%d\n", data->number_of_sequences, data->encoding);
        svg_end_tag(output, "text");


        draw_svg_length(data, translate, max_length);

        translate.y += PERF_SIZE + GRAPH_PAD;

        draw_svg_score(data, (fpair_t){score_locations[i], translate.y}, i, max_score);


        draw_svg_quality(data, translate);

        translate.y += GRAPH_HEIGHT + GRAPH_PAD;
        draw_svg_content(data, translate, max_content);

        if(adapters_used) {
          translate.y += PERF_SIZE + GRAPH_PAD;
          draw_svg_adapter(data, translate, max_adapter);
        }

        if(saturation_curve) {
            translate.y += PERF_SIZE + GRAPH_PAD;
            draw_svg_saturation(data, translate, max_saturation);    
        }

        translate.x += GRAPH_PAD + GRAPH_WIDTH + 40;
    }

    if(name != NULL) svg_end_tag(output, "g");

    svg_end_tag(output, "svg");

}
