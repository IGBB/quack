#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <string.h>


#include "svg.h"
#include "seq.h"


#define svg_axis_label( posx, posy, rot, label)                 \
  svg_start_tag("text", 7,                                      \
                svg_attr("x",           "%d", posx),              \
                svg_attr("fill",        "%s", "#AAA"),            \
                svg_attr("y",           "%d", posy),              \
                svg_attr("font-family", "%s", "sans-serif"),      \
                svg_attr("font-size",   "%s", "15px"),            \
                svg_attr("text-anchor", "%s", "middle"),          \
                svg_attr("transform", "rotate(%d)", rot)          \
                );                                              \
  printf("%s\n", label);                                        \
  svg_end_tag("text");
#define svg_axis_number( posx, posy, a, number)                 \
  svg_start_tag("text", 6,                                      \
                svg_attr("x",           "%d", posx),              \
                svg_attr("fill",        "%s", "#AAA"),            \
                svg_attr("y",           "%d", posy),              \
                svg_attr("font-family", "%s", "sans-serif"),      \
                svg_attr("font-size",   "%s", "10px"),            \
                svg_attr("text-anchor", "%s", a)                  \
                );                                              \
  printf("%d\n", number);                                       \
  svg_end_tag("text");

#define svg_center_label( posx, posy, fillv, label_format, label)       \
  svg_start_tag("text", 7,                                              \
                svg_attr("x",           "%d", posx),                      \
                svg_attr("y",           "%d", posy),                      \
                svg_attr("fill",        "%s", fillv),                     \
                svg_attr("font-family", "%s", "sans-serif"),              \
                svg_attr("font-size",   "%s", "15px"),                    \
                svg_attr("font-weight", "%s", "bold"),                    \
                svg_attr("text-anchor", "%s", "middle")                   \
                );                                                      \
  printf(label_format, label);                                          \
  svg_end_tag("text");


void draw(sequence_data* data, int position, int adapters_used) {
  int i, j, x, y;
  int offset = 0;
  int sum = 0;
  int counter = 0;
  char *encoding = "";
  uint64_t min_score = UINT32_MAX;
  int max_score = 0;
  uint64_t number_of_bases = 0;
  uint64_t total_counts[91] = {0};
  float averages[500];

  // get encoding
  i = 0;
  while (i < data->max_length && strcmp(encoding, "") == 0) {
    if (min_score < 31) {
      encoding = "phred33";
    }
    for (j = 0; j < 91; j++) {
      if (j < min_score) {
        if (data->bases[i].scores[j] != 0) {
          min_score = j;
        }
      }
    }
    i++;
  }

  // get max score and score distribution and average scores at each position
  for (i = 0; i < data->max_length; i++) {
    sum = 0;
    for (j = 0; j < 91; j++) {
      if (data->bases[i].scores[j] > 0 && j > max_score) {
        max_score = j;
      }
      total_counts[j] = total_counts[j] + data->bases[i].scores[j];
      sum = sum + j*data->bases[i].scores[j];
      number_of_bases++;
    }
    averages[i] = sum/100.0;
  }

  if (max_score < 40) {
    max_score = 40;
  }
  else {
    max_score++;
  }

  if (min_score < 31) {
    encoding = "phred33";
  }
  else {
    encoding = "phred64";
    offset = 31;
    max_score = max_score - offset;
  }


  /********** File Stats ***************/
   svg_start_tag("text", 6,
                 svg_attr("x", "%d", (position==0)?355:835),
                 svg_attr("y", "%d", 20),
                 svg_attr("text-anchor", "%s", "middle"),
                 svg_attr("font-family", "%s", "sans-serif"),
                 svg_attr("font-size", "%s", "15px"),
                 svg_attr("fill", "%s", "#555")
                 );
   svg_start_tag("tspan", 0);
   printf("%d", data->number_of_sequences);
   svg_end_tag("tspan");
   svg_start_tag("tspan", 1, svg_attr("fill", "%s", "#888"));
   printf("&#160;reads with endcoding&#160;");
   svg_end_tag("tspan");
   svg_start_tag("tspan", 0);
   printf("%s", encoding);
   svg_end_tag("tspan");
   svg_end_tag("text");
  
  /* Group for rug plot */
  svg_start_tag("g", 1, 
                svg_attr("transform", "translate(%d %d)", 5, 25)
                );

  /* Horizontal Tick marks */
  x = 100;
  if (position==1)
    x = 1000;
  for (i = 10; i < 100; i+=10){
    y = 105 + i * 250 / 100;
    
    svg_simple_tag("line",6,
                 svg_attr("x1", "%d", x),
                 svg_attr("x2", "%d", x+100),
                 svg_attr("y1", "%d", y),
                 svg_attr("y2", "%d", y),
                 svg_attr("stroke", "%s", "black"),
                 svg_attr("stroke-width", "%f", (i%20 == 10)?1:0.5)
                 );
  }


  /* Group for the vertical section of rug plot */
  svg_start_tag("g", 1, 
                svg_attr("transform", "translate(%d 0)", (position == 1)?610:130)
                );
 

  /*************** Vertical Tick Marks ***************/
  y = (adapters_used == 0)?400:500;
  for (i = 10; i < 100; i+=10){
    x = i * 450 / 100;
  svg_simple_tag("line",6,
                 svg_attr("x1", "%d", x),
                 svg_attr("x2", "%d", x),
                 svg_attr("y1", "%d", 10),
                 svg_attr("y2", "%d", y),
                 svg_attr("stroke", "%s", "black"),
                 svg_attr("stroke-width", "%f", (i%20 == 10)?1:0.5)
                 );
  }


  /*************** Base Ratio ***************/

  /* Flip svg to make svg coordinate system match cartesian coordinate. Must use
     a group since imagemagick uses SVG v1.1*/
  svg_start_tag("g", 1, svg_attr("transform", "translate(%d,%d) scale(%d, %d)", 0, 100, 1,-1));
                
  svg_start_tag("svg", 4,
                svg_attr("width",  "%d", 450),
                svg_attr("height", "%d", 100),
                svg_attr("preserveAspectRatio", "%s", "none"),
                svg_attr("viewBox", "0 0 %d %d", data->max_length, data->number_of_sequences)
                );

  /* Set background color */
  svg_simple_tag("rect", 3,
                 svg_attr("width",  "%s", "100%"),
                 svg_attr("height", "%s", "100%"),
                 svg_attr("fill", "%s", "#CCC")
                 );
                 
  
  /* Allocate 25 characters per base (A,C,T,G) per point. 
     4 = max length is capped in kb range, will start compressing if larger
     2 = '.5' added to point
     1 = ','
     15 = number of reads (10 quadrillion reads will break it)
     1 = ' '
     2 = padding for miscalculation
     
  */
  size_t ratio_points_length = 25*data->max_length;
  char * ratio_points[4];
  char tmp[25];
  for(i = 0; i < 4; i++){
    ratio_points[i] = malloc(ratio_points_length);
  }

  /* Since coordinates for lines and rectangles don't work the same; set the
     first point of each line to start off graph. Then, add 0.5 to the x of each
     point. Finally, end the line off graph. */
  y = 0;
  for(i = 0; i < 4; i++){
    y += data->bases[0].content[i];
    snprintf(ratio_points[i], ratio_points_length, "0,%d ", y);
  }

  /* Calculate cumlative sum for each x position and add it to point string */
  for (x = 0; x < data->max_length; x++) {
    y = 0;
    for(i = 0; i < 4; i++ ){
      y += data->bases[x].content[i];
      
      snprintf(tmp, 20, "%d.5,%d ", x, y);
      strncat(ratio_points[i], tmp, ratio_points_length);
    }
  }


  /* Make lines end off graph */
  y = 0;
  for(i = 0; i < 4; i++){
    y += data->bases[data->max_length-1].content[i];

    snprintf(tmp, 20, "%d,%d ", data->max_length, y);
    strncat(ratio_points[i], tmp, ratio_points_length);
  }

  
  /* Draw each distribution, in decending order so they stack */
  char *ratio_labels[4] = {"%A", "%T", "%C", "%G"};
  char *ratio_colors[4] = {"#648964", "#89bc89", "#84accf", "#5d7992"};
  for(i = 3; i >= 0; i--){
    svg_simple_tag("polyline", 3,
                   svg_attr("points",      "0,0 %s %d,0", ratio_points[i], data->max_length),
                   svg_attr("fill", "%s", ratio_colors[i]),
                   svg_attr("stroke", "%s", "none")
                   );
  }

  for(i = 0; i < 4; i++)
    free(ratio_points[i]);

  svg_end_tag("svg"); // Base Ratio
  svg_end_tag("g"); // Base Ratio

  /* Lables */

  svg_start_tag("text", 5,
                svg_attr("y",           "%d", 95),
                svg_attr("fill",        "%s", "#CCC"),
                svg_attr("x",           "%d", 5),
                svg_attr("font-family", "%s", "sans-serif"),
                svg_attr("font-size",   "%s", "15px")
                );
  printf("%s\n", "Base Content Percentage");
  svg_end_tag("text");


  /* Print Ratio labels in center, only needed the first time
     465 = width of vertical section (450) + width of margin (30) halved 
   */
  if(position==0){
    for( i = 0; i < 4; i++){
      svg_center_label(465, 20*(4-i), ratio_colors[i], "%s", ratio_labels[i]);
    }
  }

  if(position == 0){
    svg_axis_label(-50, -5, -90, "Percent");
    svg_axis_number(-5, 100, "end", 0);
    svg_axis_number(-5, 5,   "end", 100);
  }else{
    svg_axis_label(50, -455, 90, "Percent");
    svg_axis_number(455, 100, "start", 0);
    svg_axis_number(455, 5,   "start", 100);
  }
  
  /*************** Heatmap ***************/
  /* Flip svg to make svg coordinate system match cartesian coordinate. The y
     is negative to compensate for the horizontal flip */
    svg_start_tag("g", 1, svg_attr("transform", "translate(%d,%d) scale(%d, %d)", 0, 355, 1,-1));

    svg_start_tag("svg", 4,
                svg_attr("width",  "%d", 450),
                svg_attr("height", "%d", 250),
                svg_attr("preserveAspectRatio", "%s", "none"),
                svg_attr("viewBox", "0 0 %d %d", data->max_length, max_score)
                );
                  

  /* Define background for heatmap. Must be in descending order or highest score
     background will overwrite all other backgrounds */
  #define score_back(score, color)                      \
    svg_simple_tag("rect", 6,                           \
                   svg_attr("x",      "%d", 0),           \
                   svg_attr("y",      "%d", 0),           \
                   svg_attr("width",  "%s", "100%"),      \
                   svg_attr("height", "%d", score),       \
                   svg_attr("stroke", "%s", "none"),      \
                   svg_attr("fill",   "%s", color)        \
                   )
  score_back(max_score, "#ccebc5"); // green
  score_back(28,        "#ffffcc"); // yellow
  score_back(20,        "#fbb4ae"); // red

  /* Allocate 15 characters per point. 
     4 = max length is capped in kb range, will start compressing if larger
     2 = '.5' added to point
     1 = ','
     5 = '##.##' for float average
     1 = ' '
     2 = padding for miscalculation
   */
  size_t mean_line_points_length = 15*data->max_length;
  char * mean_line_points = malloc(mean_line_points_length);

  #define running_average_length 5
  float running_average_data[running_average_length], running_average;

  for (i = 0; i < running_average_length; i++){
    running_average_data[i] = averages[i];
  }
  
  /* Make line start off graph */
  snprintf(mean_line_points, mean_line_points_length, "0,%0.2f ", averages[0]);

  for (x = 0; x < data->max_length; x++) {
    for (y = offset; y < max_score+offset; y++) {
      if(data->bases[x].scores[y] > 0)
        svg_simple_tag("rect", 8,
                       svg_attr("x",      "%d", x),
                       svg_attr("y",      "%d", y),
                       svg_attr("fill-opacity", "%f", (float)(data->bases[x].scores[y])/100.0),
                       svg_attr("width",  "%d", 1),
                       svg_attr("height", "%d", 1),
                       svg_attr("stroke", "%s", "none"),
                       svg_attr("stroke-width", "%d", 0),
                       svg_attr("fill",   "%s", "black")
                       );
    }

    /* Replace the oldest data point with a new one. This will make the first
       $running_average_length average the same; however, after that they should
       start moving. */
    running_average_data[x%running_average_length] = averages[x];
    running_average = 0;
    for (i = 0; i < running_average_length; i++)
      running_average += running_average_data[i];
    running_average /= running_average_length;

    
    snprintf(tmp, 20, "%d.5,%0.2f ", x, running_average);
    strncat(mean_line_points, tmp, mean_line_points_length);
  }

  /* Make line end off graph */
  snprintf(tmp, 20, "%d,%0.2f", data->max_length, running_average);
  strncat(mean_line_points, tmp, mean_line_points_length);

  
  /* Print mean line */
  svg_simple_tag("polyline", 6,
                 svg_attr("points",      "%s", mean_line_points),
                 svg_attr("stroke", "%s", "black"),
                 svg_attr("stroke-width", "%f", 0.5),
                 svg_attr("stroke-opacity", "%f", 0.5),
                 svg_attr("fill",   "%s", "none"),
                 svg_attr("stroke-linejoin", "%s", "round")
                 );
   
  free(mean_line_points);
    
  svg_end_tag("svg"); // Heatmap
  svg_end_tag("g"); // Heatmap


  /* Lables */

  svg_start_tag("text", 5,
                svg_attr("y",           "%d", 350),
                svg_attr("fill",        "%s", "#888"),
                svg_attr("x",           "%d", 5),
                svg_attr("font-family", "%s", "sans-serif"),
                svg_attr("font-size",   "%s", "15px")
                );
  printf("%s\n", "Per Base Sequence Quality");
  svg_end_tag("text");


  /* Print base quality labels in center, only needed the first time
     465 = width of vertical section (450) + width of margin (30) halved 
     112 = start of graph (105) + half size of text (7)
   */
  if(position==0){
    svg_center_label(465, 112, "#888", "%d", max_score);
    svg_center_label(465, 112+(int)((max_score-28)*250/max_score), "#888", "%d", 28);
    svg_center_label(465, 112+(int)((max_score-20)*250/max_score), "#888", "%d", 20);
  }

  
  /*************** Length Distro ***************/

  /* Length Distro graph grows away from heatmap. No need to flip or have
     negative y*/
  svg_start_tag("svg", 6,
                svg_attr("x",      "%d", 0),
                svg_attr("y",      "%d", 360),
                svg_attr("width",  "%d", 450),
                svg_attr("height", "%d", 100),
                svg_attr("preserveAspectRatio", "%s", "none"),
                svg_attr("viewBox", "0 0 %d %d", data->max_length, data->number_of_sequences)
                );

  /* Set background color */
  svg_simple_tag("rect", 3,
                 svg_attr("width",  "%s", "100%"),
                 svg_attr("height", "%s", "100%"),
                 svg_attr("fill", "%s", "#EEE")
                 );
                 
  for (x = 0; x < data->max_length; x++) {
    if( data->bases[x].length_count > 0)
      svg_simple_tag("rect", 6,
                     svg_attr("x",      "%d", x),
                     svg_attr("y",      "%d", 0),
                     svg_attr("width",  "%d", 1),
                     svg_attr("height", "%d", data->bases[x].length_count),
                     svg_attr("stroke", "%s", "none"),
                     svg_attr("fill",   "%s", "steelblue")
                     );
  }

  
  
  svg_end_tag("svg"); // Length Distro


  /* Lables */

  svg_start_tag("text", 5,
                svg_attr("y",           "%d", 455),
                svg_attr("fill",        "%s", "#888"),
                svg_attr("x",           "%d", 5),
                svg_attr("font-family", "%s", "sans-serif"),
                svg_attr("font-size",   "%s", "15px")
                );
  printf("%s\n", "Length Distribution");
  svg_end_tag("text");

  if(position == 0){
    svg_axis_label(-410, -5, -90, "Percent");
    svg_axis_number(-5, 370, "end", 0);
    svg_axis_number(-5, 460, "end", 100);
  }else{
    svg_axis_label(410, -455, 90, "Percent");
    svg_axis_number(455, 370, "start", 0);
    svg_axis_number(455, 460, "start", 100);
  }

  
  /*************** Adapter Distro ***************/

  if(adapters_used==1){
    /* Adapter Distro graph grows away from heatmap. No need to flip or have
       negative y*/
    svg_start_tag("svg", 6,
                  svg_attr("x",      "%d", 0),
                  svg_attr("y",      "%d", 465),
                  svg_attr("width",  "%d", 450),
                  svg_attr("height", "%d", 100),
                  svg_attr("preserveAspectRatio", "%s", "none"),
                  svg_attr("viewBox", "0 0 %d %d", data->max_length, data->number_of_sequences)
                  );

    /* Set background color */
    svg_simple_tag("rect", 3,
                   svg_attr("width",  "%s", "100%"),
                   svg_attr("height", "%s", "100%"),
                   svg_attr("fill", "%s", "#EEE")
                   );
                 
    for (x = 0; x < data->max_length; x++) {
      if( data->bases[x].kmer_count > 0)
        svg_simple_tag("rect", 6,
                       svg_attr("x",      "%d", x),
                       svg_attr("y",      "%d", 0),
                       svg_attr("width",  "%d", 1),
                       svg_attr("height", "%d", data->bases[x].kmer_count),
                       svg_attr("stroke", "%s", "none"),
                       svg_attr("fill",   "%s", "steelblue")
                       );
    }

    svg_end_tag("svg"); // Adapter Distro

    /* Lables */

    svg_start_tag("text", 5,
                  svg_attr("y",           "%d", 560),
                  svg_attr("fill",        "%s", "#888"),
                  svg_attr("x",           "%d", 5),
                  svg_attr("font-family", "%s", "sans-serif"),
                  svg_attr("font-size",   "%s", "15px")
                  );
    printf("%s\n", "Adapter Distribution");
    svg_end_tag("text");

    if(position == 0){
      svg_axis_label(-515, -5, -90, "Percent");
      svg_axis_number(-5, 475, "end", 0);
      svg_axis_number(-5, 565, "end", 100);
    }else{
      svg_axis_label(515, -455, 90, "Percent");
      svg_axis_number(455, 475, "start", 0);
      svg_axis_number(455, 565, "start", 100);
    }
  }

  /*************** Bottom Label ***************/
  y = 470;
  if(adapters_used==1) y+=105;
  
  svg_axis_label(225,  y+5, 0, "Base Pairs");
  svg_axis_number(0,   y, "middle", 0);
  svg_axis_number(450, y, "middle", 100);

  
  svg_end_tag("g"); // rug plot vertical section

  /*************** Score Distro ***************/

  /* Score Distro graph is shown on either side of the heatmap, so it needs to
     flip vertically if position == 0. Flipped horizontally as well. Use
     negative positions in directions svg is flipped*/

  svg_start_tag("g", 1,
                svg_attr("transform", "translate(%d,%d) scale(%d, %d)",
                         (position == 0)?125:1065, 355, (position == 0)?-1:1,-1));

  svg_start_tag("svg", 4,
                  svg_attr("width",  "%d", 100),
                  svg_attr("height", "%d", 250),
                  svg_attr("preserveAspectRatio", "%s", "none"),
                  svg_attr("viewBox", "0 0 100 %d", max_score)
                  );

    /* Set background color */
    svg_simple_tag("rect", 3,
                   svg_attr("width",  "%s", "100%"),
                   svg_attr("height", "%s", "100%"),
                   svg_attr("fill", "%s", "#EEE")
                   );
                 
    for (y = 0; y < max_score; y++) {
      if(total_counts[y] > 0)
        svg_simple_tag("rect", 6,
                       svg_attr("x",      "%d", 0),
                       svg_attr("y",      "%d", y),
                       svg_attr("width",  "%d", total_counts[y]*100/number_of_bases),
                       svg_attr("height", "%d", 1),
                       svg_attr("stroke", "%s", "none"),
                       svg_attr("fill",   "%s", "steelblue")
                       );
    }

    svg_end_tag("svg"); // Score Distro
    svg_end_tag("g"); // Score Distro

    /* Lables */

    if(position == 0){
      svg_start_tag("text", 5,
                    svg_attr("y",           "%d", 335),
                    svg_attr("fill",        "%s", "#888"),
                    svg_attr("x",           "%d", 30),
                    svg_attr("font-family", "%s", "sans-serif"),
                    svg_attr("font-size",   "%s", "15px")
                    );
      svg_start_tag("tspan", 0);
      printf("%s\n", "Score");
      svg_end_tag("tspan");
    
      svg_start_tag("tspan", 2, svg_attr("dy", "%d", 15), svg_attr("x", "%d", 30));
      printf("%s\n", "Distribution");
      svg_end_tag("tspan");

      svg_end_tag("text");

      svg_axis_label(72, 100, 0, "Percent");
      /* svg_axis_number(125, 100, "middle", 0); */
      svg_axis_number(25,  100, "middle", 100);

      svg_axis_label(-230, 20, -90, "Score");
      svg_axis_number(20, 110, "end", max_score);
      svg_axis_number(20, 355, "end", 1);

    }else{
      svg_start_tag("text", 5,
                    svg_attr("y",           "%d", 335),
                    svg_attr("fill",        "%s", "#888"),
                    svg_attr("x",           "%d", 1070),
                    svg_attr("font-family", "%s", "sans-serif"),
                    svg_attr("font-size",   "%s", "15px")
                    );
      svg_start_tag("tspan", 0);
      printf("%s\n", "Score");
      svg_end_tag("tspan");
    
      svg_start_tag("tspan", 2, svg_attr("dy", "%d", 15), svg_attr("x", "%d", 1070));
      printf("%s\n", "Distribution");
      svg_end_tag("tspan");

      svg_end_tag("text");

      svg_axis_label(1115,  100, 0, "Percent");
      /* svg_axis_number(1070, 100, "middle", 0); */
      svg_axis_number(1165, 100, "middle", 100);

      svg_axis_label(230, -1170, 90, "Score");
      svg_axis_number(1170, 110, "start", max_score);
      svg_axis_number(1170, 355, "start", 1);
    }

    svg_end_tag("g"); // rug plot 

}
