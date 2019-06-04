#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>

#include "kseq.h"
#include "svg.h"

#define unlikely(x) __builtin_expect ((x), 0)
#define likely(x)       __builtin_expect((x),1)


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



const char *program_version = "quack 1.1.1";
struct arguments {
    char *name, *forward, *reverse, *unpaired, *adapters;
};

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

typedef struct {
    uint64_t scores[91];
    uint64_t content[4];
    uint64_t length_count;
    uint64_t kmer_count;
} base_information;

typedef struct {
    base_information *bases;
    uint64_t max_length;
    uint64_t original_max_length;
    uint64_t number_of_sequences;
} sequence_data;

/* Convert ASCII to Integer for A T C and G */
/*                A     C           G                                      T */
int lookup[20] = {0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

KSEQ_INIT(gzFile, gzread)

int* read_adapters(char *adapters_file) {
    int kmer_size = 10;
    int array_size = pow(4, kmer_size);
    gzFile fp;
    kseq_t *seq;
    int i, l, index;
    fp = gzopen(adapters_file, "r");
    seq = kseq_init(fp);
    int *kmers = malloc(array_size*sizeof(int));
    memset(kmers, 0, array_size*sizeof(int));
    while ((l = kseq_read(seq)) >= 0) {
        index = 0;
            for (i = 0; i < kmer_size; i++) {
                index = ((index << 2) + (lookup[seq->seq.s[i]-65 & ~32])) & (array_size-1);
            }
        for (; i < seq->seq.l; i++) {
            index = ((index << 2) + (lookup[seq->seq.s[i]-65 & ~32])) & (array_size-1);
            kmers[index] = 1;
        }

    }
    kseq_destroy(seq);
    gzclose(fp);
    return kmers;
}

sequence_data* read_fastq(char *fastq_file, int *kmers) {
    gzFile fp;
    kseq_t *seq;
    int i, l, index;
    int kmer_size = 10;
    int array_size = pow(4, kmer_size);
    int max_length = 0;
    fp = gzopen(fastq_file, "r");
    seq = kseq_init(fp);
    base_information *bases = NULL;
    sequence_data *to_return = malloc(sizeof(sequence_data));
    int number_of_sequences = 0;

    while ((l = kseq_read(seq)) >= 0) {
        if (unlikely(seq->seq.l > max_length)) {
            bases = realloc(bases, seq->seq.l*sizeof(base_information));
            memset(bases+max_length, 0, (seq->seq.l - max_length)*sizeof(base_information));
            max_length = seq->seq.l;
        }
        for (i = 0; i < seq->seq.l; i++) {
            int base = seq->seq.s[i];
            int offset = lookup[base-65 & ~32];
            bases[i].content[offset]++;
            int quality = seq->qual.s[i]-33;
            bases[i].scores[quality]++;
        }
        index = 0;
        for (i = 0; i < kmer_size; i++) {
            index = ((index << 2) + (lookup[seq->seq.s[i]-65 & ~32])) & (array_size-1);
        }
        if (kmers) {
            for (; kmers[index] == 0 && i < seq->seq.l; i++) {
                index = ((index << 2) + (lookup[seq->seq.s[i]-65 & ~32])) & (array_size-1);
            }
        }
        if (i < seq->seq.l) {
            bases[i].kmer_count++;
        }

        bases[seq->seq.l-1].length_count++;
        number_of_sequences++;
    }
    kseq_destroy(seq);
    gzclose(fp);
    to_return->bases = bases;
    to_return->max_length = max_length;
    to_return->number_of_sequences = number_of_sequences;
    return to_return;
}

sequence_data* transform(sequence_data* data) {
    int i, j;
    data->original_max_length = data->max_length;
    // binning
    if (data->max_length > 3000) {
        fprintf(stderr, "Binning...\n");
        int bin_size = 100;
        int unbinned;
        int binned = 0;

        for (unbinned = 1; unbinned < data->max_length; unbinned++) {
             if (unbinned%bin_size == 0){
                binned++;
                for (i = 0; i < 4; i++) {
                    data->bases[binned].content[i] = 0;
                }
                for (i = 0; i < 91; i++) {
                    data->bases[binned].scores[i] = 0;
                }
                data->bases[binned].length_count = 0;
             }
             for (i = 0; i < 4; i++) {
                data->bases[binned].content[i] = data->bases[binned].content[i] + data->bases[unbinned].content[i];
             }
             for (i = 0; i < 91; i++) {
                data->bases[binned].scores[i] = data->bases[binned].scores[i] + data->bases[unbinned].scores[i];
             }
            // fprintf(stderr, "%d\n", data->bases[binned].length_count);
            data->bases[binned].length_count = data->bases[binned].length_count + data->bases[unbinned].length_count;
            data->bases[binned].kmer_count = data->bases[binned].kmer_count + data->bases[unbinned].kmer_count;
        }
        data->max_length = binned;
    }

    for (i = 1; i < data->max_length; i++) {
        data->bases[i].kmer_count = data->bases[i-1].kmer_count + data->bases[i].kmer_count;
    }

    // transforming data for drawing (percentages rather than counts)
    for (i = 0; i < data->max_length; i++) {
        int content_sum = 0;
        int score_sum = 0;
        /* for (j = 0; j < 4; j++) { */
        /*     content_sum = content_sum + data->bases[i].content[j]; */
        /* } */
        /* if (content_sum != 0) { */
        /*     for (j = 0; j < 4; j++) { */
        /*         data->bases[i].content[j] = 100*data->bases[i].content[j]/content_sum; */
        /*     } */
        /* } */
        for (j = 0; j < 91; j++) {
            score_sum = score_sum + data->bases[i].scores[j];
        }
        if (score_sum != 0) {
            for (j = 0; j < 91; j++) {
                data->bases[i].scores[j] = 100*data->bases[i].scores[j]/score_sum;
            }
        }
        data->bases[i].length_count = ceil(100*(float)data->bases[i].length_count/data->number_of_sequences);
        data->bases[i].kmer_count = ceil(100*(float)data->bases[i].kmer_count/(float)data->number_of_sequences);

    }
    return data;
}

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
   printf(" reads with endcoding ");
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
  svg_start_tag("g", 1, svg_attr("transform", "translate(%d,%d),scale(%d, %d)", 0, 100, 1,-1));
                
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
    svg_start_tag("g", 1, svg_attr("transform", "translate(%d,%d),scale(%d, %d)", 0, 355, 1,-1));

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

    snprintf(tmp, 20, "%d.5,%0.2f ", x, averages[x]);
    strncat(mean_line_points, tmp, mean_line_points_length);
  }

  /* Make line end off graph */
  snprintf(tmp, 20, "%d,%0.2f", data->max_length, averages[data->max_length-1]);
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
                svg_attr("viewBox", "0 0 %d 100", data->max_length)
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
                  svg_attr("viewBox", "0 0 %d 100", data->max_length)
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
                svg_attr("transform", "translate(%d,%d),scale(%d, %d)",
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

int main (int argc, char **argv)
{
    struct arguments arguments;
    arguments = parse_options(argc, argv);

    int paired, unpaired, adapters;
    int width, height;
    int *kmers = NULL;

    paired = (arguments.forward != NULL && arguments.reverse != NULL);
    unpaired = (arguments.unpaired != NULL);
    adapters = (arguments.adapters != NULL);
    
    /* If paired and unparied data are both set or unset, then throw error */
    if(paired == unpaired){
      printf("%s\n", "Usage: quack [OPTION...]\nTry `quack --help' or `quack --usage' for more information.");
      exit(1);
    }

    if(adapters) kmers = read_adapters(arguments.adapters);

    width  = (paired)?1195:615;
    height = (adapters)?610:510;

    if(arguments.name != NULL)
      height += 30;
    
    svg_start_tag("svg", 5,
                  svg_attr("width",   "%d", width),
                  svg_attr("height",  "%d", height),
                  svg_attr("viewBox", "%d %d %d %d", 0, 0, width, height),
                  svg_attr("xmlns",       "%s", "http://www.w3.org/2000/svg"),
                  svg_attr("xmlns:xlink", "%s", "http://www.w3.org/1999/xlink")
                  );

    // If name is given, add to middle of viewBox (half of width + min-x of viewbox)
    if(arguments.name != NULL){
      svg_start_tag("text", 6,
                    svg_attr("x", "%d", (width/2)),
                    svg_attr("y", "%d", 30),
                    svg_attr("font-family", "%s", "sans-serif"),
                    svg_attr("text-anchor", "%s", "middle"),
                    svg_attr("font-size",   "%s", "30px"),
                    svg_attr("fill",        "%s", "black"));
      printf("%s", arguments.name);
      svg_end_tag("text");

      svg_start_tag("g", 1, 
                svg_attr("transform", "translate(%d %d)", 0, 30)
                );

    }
  
    sequence_data *data = read_fastq(((paired)?arguments.forward:arguments.unpaired), kmers);
    sequence_data *transformed_data = transform(data);
    draw(transformed_data, 0, adapters);
    free(data);
    
    if(paired){
      data = read_fastq(arguments.reverse, kmers);
      transformed_data = transform(data);
      draw(transformed_data, 1, adapters);
      free(data);
    }

    if(arguments.name != NULL) svg_end_tag("g");

    svg_end_tag("svg");
    
    exit (0);
}
