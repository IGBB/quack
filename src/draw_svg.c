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
#define GRAPH_PAD 5

char *ratio_labels[4] = {"%A", "%T", "%C", "%G"};
char *ratio_colors[4] = {"#648964", "#89bc89", "#84accf", "#5d7992"};

#define PI 3.14159265

typedef struct {
    float x,y;
} fpair_t;

char * point_string(int length, fpair_t* p){
    char * points, *cur_pos;
    size_t max = 0;
    int i = 0;

    /* Find lenght of each point string */
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

#define draw_svg_background(wdth, hght)                          \
    svg_simple_tag("rect", 3,                                    \
                   svg_attr("width",  "%f", (float)(wdth)),      \
                   svg_attr("height", "%f", (float)(hght)),      \
                   svg_attr("fill", "%s", "none")                \
                   )

#define svg_start_label(posx, posy)                                     \
    svg_start_tag("text", 5,                                            \
                  svg_attr("x",           "%f", (float)(posx)),         \
                  svg_attr("y",           "%f", (float)(posy)),         \
                  svg_attr("fill",        "%s", "#AAA"),                \
                  svg_attr("font-family", "%s", "sans-serif"),          \
                  svg_attr("font-size",   "%s", "15px")                 \
    )

#define svg_end_label() svg_end_tag("text")

void draw_svg_axis_label(fpair_t start, fpair_t end) {

    fpair_t middle = (fpair_t){ (start.x + end.x)/2.0 ,
                                (start.y + end.y)/2.0};

    char * baselines[2] = {"hanging", "baseline"};
    char * anchor[2] = {"end", "end"};
    if( start.y > end.y ){
        baselines[0] = "baseline";
        baselines[1] = "hanging";
    }else if (start.y == end.y){
        baselines[0] = "baseline";
        anchor[1] = "start";
    }

    double rotate = atan2((end.y - start.y), (end.x -start.x)) * 180 / PI;

    /* Convert to int for comparison to avoid odd floating point encoding artifats
     * For example, atan2 returns -90.000000000000001 for a vertical line.
     * */
    if((int)rotate < -90 || (int)rotate >= 90)
        rotate -= 180;

    svg_start_tag("text", 8,
                  svg_attr("x",           "%f", start.x),
                  svg_attr("y",           "%f", start.y),
                  svg_attr("fill",        "%s", "#AAA"),
                  svg_attr("font-family", "%s", "sans-serif"),
                  svg_attr("font-size",   "%s", "10px"),
                  svg_attr("alignment-baseline",   "%s", baselines[0]),
                  svg_attr("dominant-baseline",   "%s", baselines[0]),
                  svg_attr("text-anchor",   "%s", anchor[0])
                  );
    printf("0");
    svg_end_tag("text");

    svg_start_tag("text", 8,
                  svg_attr("x",           "%f", end.x),
                  svg_attr("y",           "%f", end.y),
                  svg_attr("fill",        "%s", "#AAA"),
                  svg_attr("font-family", "%s", "sans-serif"),
                  svg_attr("font-size",   "%s", "10px"),
                  svg_attr("alignment-baseline",   "%s", baselines[1]),
                  svg_attr("dominant-baseline",   "%s", baselines[1]),
                  svg_attr("text-anchor",   "%s", anchor[1])
                  );
    printf("100");
    svg_end_tag("text");

    svg_start_tag("text", 8,
                  svg_attr("x",           "%f", middle.x),
                  svg_attr("y",           "%f", middle.y),
                  svg_attr("fill",        "%s", "#AAA"),
                  svg_attr("font-family", "%s", "sans-serif"),
                  svg_attr("font-size",   "%s", "12px"),
                  svg_attr("dominant-baseline",   "%s", "baseline"),
                  svg_attr("text-anchor",   "%s", "middle"),
                  svg_attr("transform", "rotate(%f %f,%f)", rotate, middle.x, middle.y)
                  );
    printf("Percent");
    svg_end_tag("text");


}



void draw_svg_distro(int length, fpair_t* p,
                     fpair_t original, fpair_t final,
                     fpair_t flip, fpair_t translate,
                     float rotate, char* label){

    fpair_t scale = (fpair_t){final.x/original.x, final.y/original.y};
    fpair_t scale_translate = (fpair_t){0,0};

    int i;

    if(flip.x){
        scale.x *= -1;
        scale_translate.x = final.x;
    }
    if(flip.y){
        scale.y *= -1;
        scale_translate.y = final.y;
    }

    svg_start_tag("g", 1, svg_attr("transform", "translate(%f %f)", translate.x, translate.y));
    svg_start_tag("g", 1, svg_attr("transform", "translate(%f %f) scale(%f %f) rotate(%f)",
                                   scale_translate.x, scale_translate.y,
                                   scale.x, scale.y, rotate
                  )
    );

    draw_svg_background(original.x, original.y);

    /* /\* Draw distribution *\/ */
    /* char* points = point_string(length, p); */

    /* svg_simple_tag("polyline", 5, */
    /*                svg_attr("points", "%s", points), */
    /*                svg_attr("stroke", "%s", "steelblue"), */
    /*                svg_attr("fill", "%s", "none"), */
    /*                svg_attr("stroke-width", "%fpx", 0.5), */
    /*                svg_attr("vector-effect", "%s", "non-scaling-stroke") */
    /* ); */

    /* free(points); */

    for( i = 0; i < length; i++  )
        svg_simple_tag("rect",5,
                       svg_attr("x",       "%f", p[i].x),
                       svg_attr("height",  "%f", p[i].y),
                       svg_attr("y",       "%f", 0.0),
                       svg_attr("width",   "%f", 1.1),
                       svg_attr("fill",   "%s", "steelblue")
                       );

    for(i = 20; i < 100; i+=20){
        svg_simple_tag("line", 7,
                       svg_attr("y1", "%d", i),
                       svg_attr("y2", "%d", i),
                       svg_attr("x1", "%d", 0),
                       svg_attr("x2", "%f", original.x),
                       svg_attr("stroke", "%s", "white"),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                       svg_attr("stroke-width", "%f", 1.0)

                       );
        svg_start_tag("text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "start"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %f) scale(%f %f)",
                               original.x + 1, (float)i,
                               1/scale.x, 1/scale.y ));
        printf("%d%%", i);
        svg_end_tag("text");
    }
    svg_simple_tag("line", 7,
                   svg_attr("y1", "%d", 0),
                   svg_attr("y2", "%d", 0),
                   svg_attr("x1", "%d", 0),
                   svg_attr("x2", "%f", original.x),
                   svg_attr("stroke", "%s", "black"),
                   svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                   svg_attr("stroke-width", "%f", 0.5));



    svg_end_tag("g");

    /* Add graph label */
    svg_start_label(GRAPH_PAD, final.y - GRAPH_PAD);
    printf("%s\n", label);
    svg_end_label();

    svg_end_tag("g");

}



void draw_svg_content(sequence_data* data, int translate_x, int translate_y){

    int i, j;
    uint64_t len, y;
    char * points;
    float max;

    fpair_t *perc[4];
    for(i = 0; i < 4; i++){
        perc[i] = malloc(sizeof(fpair_t[data->max_length]));
    }

    /* Calculate cumulative percentage of base content, in decending order so
     * they stack */
    for(i = 0; i < data->max_length; i++){

        y = 0;
        len = data->bases[i].length_count;
        for (j = 3; j >= 0; j--) {
            perc[j][i].x = i +0.5;
            perc[j][i].y = data->bases[i].content[j] * 100.0 / len;

            if(max < perc[j][i].y) max = perc[j][i].y;
        }
    }

    max = ceil(max / 10.0) * 10.0;

    svg_start_tag("g", 1, svg_attr("transform", "translate(%d %d)", translate_x, translate_y));
    svg_start_tag("g", 1,
                  svg_attr("transform", "translate(%f %f) scale(%f %f)",
                           0.0, PERF_SIZE,
                           GRAPH_WIDTH/data->max_length, -1 * PERF_SIZE/max)

                  );


    draw_svg_background(data->max_length, max);

    for (i = 10; i < max; i += 10){
                svg_simple_tag("line", 7,
                       svg_attr("y1", "%d", i),
                       svg_attr("y2", "%d", i),
                       svg_attr("x1", "%d", 0),
                       svg_attr("x2", "%d", data->max_length),
                       svg_attr("stroke", "%s", "#AAA"),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                       svg_attr("stroke-width", "%f", 0.5));
    }

    /* Draw each distribution */
    for(i = 0; i < 4; i++){

        points = point_string(data->max_length, perc[i]);

        svg_simple_tag("polyline", 5,
                      svg_attr("points", "%s", points),
                       svg_attr("fill", "%s", "none"),
                       svg_attr("stroke", "%s", ratio_colors[i]),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                       svg_attr("stroke-width", "%f", 2.0)
        );

        free(points);

        free(perc[i]);
    }


    for (i = 10; i < max; i += 10){
        svg_start_tag("text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "start"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%d %d) scale(%f %f)",
                               data->max_length + 1, i,
                           1.0/(GRAPH_WIDTH/data->max_length), -1.0/(PERF_SIZE/max)));
        printf("%d%%", i);
        svg_end_tag("text");


    }

        /* Find a tick step that is
     * - at least 5
     * - divisible by 5,
     * - creates 10 ticks */
    int step = data->max_length/10;
    step = (step /5 )*5;
    if (step < 5) step = 5;



    for (i = step; i <= data->max_length - step; i += step){
        svg_simple_tag("line", 8,
                       svg_attr("x1", "%f", i-0.5),
                       svg_attr("y1", "%d", -1),
                       svg_attr("x2", "%f", i-0.5),
                       svg_attr("y2", "%f", max),
                       svg_attr("stroke", "%s", "black"),
                       svg_attr("stroke-opacity", "%f", 0.25),
                       svg_attr("stroke-dasharray", "%d %d", 1, 9),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke")
                       );
    }


    svg_end_tag("g");




    /* Add graph label */
    svg_start_label(GRAPH_PAD, PERF_SIZE - GRAPH_PAD);
    printf("%s\n", "Base Content Percentage");
    svg_end_label();


    for(i = 0; i<4; i++){
        svg_start_tag("text", 5,
                      svg_attr("x",           "%f", -15.0),
                      svg_attr("y",           "%f", PERF_SIZE * (i+1)/ 6.0),
                      svg_attr("fill",        "%s", ratio_colors[i]),
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("font-size",   "%s", "15px")
                      );
        printf("%c", rev_lookup[i] );
        svg_end_tag("text");
    }

    svg_end_tag("g");

}


void draw_svg_quality(sequence_data* data, int translate_x, int translate_y){
    int i, j;
    float perc;
    fpair_t *p;

    p = malloc(sizeof(fpair_t[data->max_length]));

    svg_start_tag("g", 1, svg_attr("transform", "translate(%d %d)", translate_x, translate_y));
    svg_start_tag("g", 1,
                  svg_attr("transform", "translate(0, %f) scale(%f %f)",
                           GRAPH_HEIGHT,
                           GRAPH_WIDTH/data->max_length,
                           -1.0 * GRAPH_HEIGHT/(data->max_score+1))
                  );

    int score_delineation[4] = {0,       20,       28, data->max_score+1};
    char* score_color[4]      = {"", "#d73027","#ffffbf",         "#1a9850"};

    for(i = 1; i <= 3; i++){
        float line_width = 2.5;
        int height = score_delineation[i] - score_delineation[i-1];
        svg_simple_tag("rect", 5,
                       svg_attr("y", "%d", score_delineation[i-1]),
                       svg_attr("width", "%d", data->max_length),
                       svg_attr("height", "%d", height),
                       svg_attr("fill", "%s", score_color[i]),
                       svg_attr("fill-opacity", "%f", 0.1)
                       );
        /* svg_simple_tag("line", 7, */
        /*                svg_attr("x1", "%d", 0), */
        /*                svg_attr("x2", "%d", 0), */
        /*                svg_attr("y1", "%d", score_delineation[i] ), */
        /*                svg_attr("y2", "%d", score_delineation[i-1]), */
        /*                svg_attr("stroke", "%s", score_color[i]), */
        /*                svg_attr("stroke-width", "%f", line_width), */
        /*                svg_attr("vector-effect", "%s", "non-scaling-stroke") */
        /*                ); */
        /* svg_simple_tag("line", 7, */
        /*                svg_attr("x1", "%d", data->max_length), */
        /*                svg_attr("x2", "%d", data->max_length), */
        /*                svg_attr("y1", "%d", score_delineation[i] ), */
        /*                svg_attr("y2", "%d", score_delineation[i-1]), */
        /*                svg_attr("stroke", "%s", score_color[i]), */
        /*                svg_attr("stroke-width", "%f", line_width), */
        /*                svg_attr("vector-effect", "%s", "non-scaling-stroke") */
        /*                ); */
    }
    /* Draw each heatmap */
    for(i = 0; i < data->max_length; i++){
        p[i].x = i + 0.5;
        p[i].y = data->avg_score[i];

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


    char* points = point_string(data->max_length, p);

    svg_simple_tag("polyline", 6,
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


    /* Find a tick step that is
     * - at least 5
     * - divisible by 5,
     * - creates 10 ticks */
    int step = data->max_length/10;
    step = (step /5 )*5;
    if (step < 5) step = 5;

    /*  invert scale */
    fpair_t scale = { 1.0 / (GRAPH_WIDTH/data->max_length),
                     -1.0 / (GRAPH_HEIGHT/(data->max_score+1)) };


    for (i = step; i <= data->max_length - step; i += step){
        svg_simple_tag("line", 8,
                       svg_attr("x1", "%f", i-0.5),
                       svg_attr("y1", "%d", -1),
                       svg_attr("x2", "%f", i-0.5),
                       svg_attr("y2", "%d", data->max_score+2),
                       svg_attr("stroke", "%s", "black"),
                       svg_attr("stroke-opacity", "%f", 0.25),
                       svg_attr("stroke-dasharray", "%d %d", 1, 9),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke")
                       );


        svg_start_tag("text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %f) scale(%f %f)",
                               i-0.5, -2.0,
                               scale.x, scale.y));

        printf("%d", i);
        svg_end_tag("text");


    }
        svg_start_tag("text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%d %f) scale(%f %f)",
                               data->max_length, -2.0,
                               scale.x, scale.y));

        printf("%d", data->max_length);
        svg_end_tag("text");



    svg_end_tag("g");

    /* Add graph label */
    svg_start_label(GRAPH_PAD, GRAPH_HEIGHT - GRAPH_PAD);
    printf("%s\n", "Per Base Sequence Quality");
    svg_end_label();

    svg_end_tag("g");

}


void draw_svg_length(sequence_data * data, int translate_x, int translate_y){
    int i, different;
    fpair_t *points = malloc(sizeof(fpair_t[data->max_length]));
    fpair_t original_size, final_size, translate, flip;
    float max;

    for(i = 0; i < data->max_length; i++){
        points[i].x = (float) i;
        points[i].y = data->bases[i].length_count*100.0/data->number_of_sequences;
    }

    different = 0;
    max = 0.0;
    for(i = 0; i < data->max_length - 1; i++){
        points[i].y -= points[i+1].y;

        /* Set different flag if current point is greater than zero */
        different |= (points[i].y > 0.0);

        if(max < points[i].y) max = points[i].y;
    }

    // check last point
    if(max < points[i].y) max = points[i].y;

    // Adjust max to nearest ten
    max = ceil(max / 10.0) * 10.0;

    original_size = (fpair_t){data->max_length, 100};
    final_size    = (fpair_t){GRAPH_WIDTH, PERF_SIZE};
    translate     = (fpair_t){translate_x, translate_y};
    flip          = (fpair_t){0,1};

    draw_svg_distro(data->max_length, points,
                    original_size, final_size, flip,
                    translate, 0.0, "Length Distribution");

    free(points);

    if(!different){
        svg_simple_tag("rect", 6,
                       svg_attr("x", "%f", translate.x),
                       svg_attr("y", "%f", translate.y),
                       svg_attr("width", "%f", final_size.x),
                       svg_attr("height", "%f", final_size.y),
                       svg_attr("fill", "%s", "#CCC"),
                       svg_attr("opacity", "%f", 0.5)
                       );

        svg_start_tag("text", 7,
                      svg_attr("x", "%f", translate.x + final_size.x/2),
                      svg_attr("y", "%f", translate.y + final_size.y/2),
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("font-variant", "%s", "small-caps"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("font-size",   "%s", "12px"),
                      svg_attr("fill",        "%s", "black"));
        printf("All reads are %d-bp", data->max_length);
        svg_end_tag("text");
}

}

void draw_svg_score(sequence_data * data, int translate_x, int translate_y, int flipx){
    int i,j;
    float total, max;

    fpair_t* points = calloc(data->max_score+1, sizeof(fpair_t));
    fpair_t original, final, translate, scale;

    for(i = 0; i <= data->max_score; i++){
        points[i].x = i;
        points[i].y = 0; 
    }

    for(i = 0; i < data->max_length; i++){
        for(j =0; j <= data->max_score; j++){
            total += data->bases[i].scores[j];
            points[j].y += data->bases[i].scores[j];
        }
    }

    max = 0.0;
    for(i = 0; i <= data->max_score ; i++){
        points[i].y = points[i].y * 100.0 / total;
        if(max < points[i].y) max = points[i].y;
    }

    // Get next ten for max
    max = ceil(max/10.0)*10.0;

    original = (fpair_t){100, data->max_score+1};
    final    = (fpair_t){PERF_SIZE, GRAPH_HEIGHT};
    translate     = (fpair_t){translate_x, translate_y};

    scale = (fpair_t){final.x/original.x, final.y/original.y};
    fpair_t scale_translate = (fpair_t){0,0};


    if(!flipx){
        scale.x *= -1;
        scale_translate.x = final.x;
    }

    scale.y *= -1;
    scale_translate.y = final.y;


    svg_start_tag("g", 1, svg_attr("transform", "translate(%f %f)", translate.x, translate.y));
    svg_start_tag("g", 1, svg_attr("transform", "translate(%f %f) scale(%f %f)",
                                   scale_translate.x, scale_translate.y,
                                   scale.x, scale.y
                  )
    );

    draw_svg_background(original.x, original.y);


    for( i = 0; i <= data->max_score; i++  )
        svg_simple_tag("rect",5,
                       svg_attr("y",      "%f", points[i].x),
                       svg_attr("width",  "%f", points[i].y),
                       svg_attr("x",      "%f", 0.0),
                       svg_attr("height", "%f", 1.0),
                       svg_attr("fill",   "%s", "steelblue")
                       );


    for(i = 20; i < 100; i+=20){
        svg_simple_tag("line", 7,
                       svg_attr("x1", "%d", i),
                       svg_attr("x2", "%d", i),
                       svg_attr("y1", "%d", 0),
                       svg_attr("y2", "%d", data->max_score+1),
                       svg_attr("stroke", "%s", "white"),
                       svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                       svg_attr("stroke-width", "%d", 1)

        );
        svg_start_tag("text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "start"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %f) scale(%f %f) rotate(-45)",
                               (float)i,  original.y + 1,
                               1/scale.x, 1/scale.y ));
        printf("%d%%", i);
        svg_end_tag("text");
    }

    svg_simple_tag("line", 7,
                   svg_attr("x1", "%d", 0),
                   svg_attr("x2", "%d", 0),
                   svg_attr("y1", "%d", 0),
                   svg_attr("y2", "%f", original.y),
                   svg_attr("stroke", "%s", "black"),
                   svg_attr("vector-effect", "%s", "non-scaling-stroke"),
                   svg_attr("stroke-width", "%f", 0.5));


    /* Find a tick step that is
     * - at least 5
     * - divisible by 5,
     * - creates 10 ticks */
    int step = data->max_score/10;
    step = (step /5 )*5;
    if (step < 5) step = 5;

    /*  loop through each step. excluded last step if too close to max score */
    for (i = step; i < data->max_score - (step/2); i += step){
        /* svg_simple_tag("line", 8, */
        /*                svg_attr("x1", "%f", i-0.5), */
        /*                svg_attr("y1", "%d", -1), */
        /*                svg_attr("x2", "%f", i-0.5), */
        /*                svg_attr("y2", "%d", data->max_score+2), */
        /*                svg_attr("stroke", "%s", "black"), */
        /*                svg_attr("stroke-opacity", "%f", 0.25), */
        /*                svg_attr("stroke-dasharray", "%d %d", 1, 9), */
        /*                svg_attr("vector-effect", "%s", "non-scaling-stroke") */
        /*                ); */


        svg_start_tag("text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %f) scale(%f %f)",
                               100.0, i+0.5,
                               1.0/scale.x, 1.0/scale.y));

        printf("%d", i);
        svg_end_tag("text");


    }
    /* add max score label */
        svg_start_tag("text", 8,
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("dominant-baseline", "%s", "middle"),
                      svg_attr("font-size",   "%s", "8px"),
                      svg_attr("vector-effect", "%s", "non-scaling-size"),
                      svg_attr("fill",        "%s", "black"),
                      svg_attr("fill-opacity",        "%f", 0.5),
                      svg_attr("transform", "translate(%f %f) scale(%f %f)",
                               100.0, data->max_score+0.5,
                               1.0/scale.x, 1.0/scale.y));

        printf("%d", data->max_score);
        svg_end_tag("text");




    svg_end_tag("g");

    /* Add graph label */
    svg_start_label(GRAPH_PAD, final.y - GRAPH_PAD);
    printf("<tspan dy=\"-15\">Score</tspan>"
           "<tspan dy=\"15\" x=\"5\">Distribution</tspan>");
    svg_end_label();

    svg_end_tag("g");


    free(points);

}

void draw_svg_adapter(sequence_data * data, int translate_x, int translate_y){
    int i;
    fpair_t *points = malloc(sizeof(fpair_t[data->max_length + 4]));
    fpair_t original_size, final_size, translate, flip;

    for(i = 0; i < data->max_length; i++){
        points[i+2].x = i + 0.5;
        points[i+2].y = data->bases[i].kmer_count*100.0/data->number_of_sequences;
    }

    points[0] = (fpair_t){0,0};
    points[1] = (fpair_t){0, points[2].y};

    points[data->max_length+2] = (fpair_t){data->max_length, points[data->max_length+1].y};
    points[data->max_length+3] = (fpair_t){data->max_length, 0};

    original_size = (fpair_t){data->max_length, 100};
    final_size    = (fpair_t){GRAPH_WIDTH, PERF_SIZE};
    translate     = (fpair_t){translate_x, translate_y};
    flip          = (fpair_t){0,0};

    draw_svg_distro(data->max_length+4, points,
                    original_size, final_size, flip,
                    translate, 0.0, "Adapter Distribution");

    free(points);


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
    if(name != NULL)
        height += 30;

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
                      svg_attr("transform", "translate(%d %d)", 0, 35)
        );

    }

    int current_y = GRAPH_PAD;
    int current_x = GRAPH_PAD*2 + PERF_SIZE;


    int i;



    svg_start_tag("text", 6,
                      svg_attr("x", "%f", (GRAPH_WIDTH/2) + PERF_SIZE + GRAPH_PAD),
                      svg_attr("y", "%d", 0),
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("font-size",   "%s", "15px"),
                      svg_attr("fill",        "%s", "black"));
        printf("%d reads encoded in phred%d\n", forward->number_of_sequences, forward->encoding);
        svg_end_tag("text");


    draw_svg_length(forward, current_x, current_y);

    current_y += PERF_SIZE + GRAPH_PAD;




    draw_svg_score(forward, GRAPH_PAD, current_y, 0);



    draw_svg_quality(forward, current_x, current_y);

    current_y += GRAPH_HEIGHT + GRAPH_PAD;

    draw_svg_content(forward, current_x, current_y);



    current_y += PERF_SIZE + GRAPH_PAD;
    if(adapters_used){
        draw_svg_adapter(forward, current_x, current_y);

    }

    if(reverse){
        current_y = GRAPH_PAD;
        current_x += GRAPH_PAD + GRAPH_WIDTH + 40;


    svg_start_tag("text", 6,
                      svg_attr("x", "%f", (GRAPH_WIDTH/2) + current_x),
                      svg_attr("y", "%d", 0),
                      svg_attr("font-family", "%s", "sans-serif"),
                      svg_attr("text-anchor", "%s", "middle"),
                      svg_attr("font-size",   "%s", "15px"),
                      svg_attr("fill",        "%s", "black"));
        printf("%d reads encoded in phred%d\n", reverse->number_of_sequences, reverse->encoding);
        svg_end_tag("text");


        draw_svg_length(reverse, current_x, current_y);

        current_y += PERF_SIZE + GRAPH_PAD;



        draw_svg_score(reverse, current_x + GRAPH_WIDTH + GRAPH_PAD, current_y, 1);
        draw_svg_quality(reverse, current_x, current_y);

        current_y += GRAPH_HEIGHT + GRAPH_PAD;

        draw_svg_content(reverse, current_x, current_y);

        current_y += PERF_SIZE + GRAPH_PAD;
        if(adapters_used)
            draw_svg_adapter(reverse, current_x, current_y);

    }


    if(name != NULL) svg_end_tag("g");

    svg_end_tag("svg");

}
