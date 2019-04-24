#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>

#include "kseq.h"

#define unlikely(x) __builtin_expect ((x), 0)
#define likely(x)       __builtin_expect((x),1)

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

int lookup[20] = {0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3};

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
        for (j = 0; j < 4; j++) {
            content_sum = content_sum + data->bases[i].content[j];
        }
        if (content_sum != 0) {
            for (j = 0; j < 4; j++) {
                data->bases[i].content[j] = 100*data->bases[i].content[j]/content_sum;
            }
        }
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
    int i, j;
    int offset = 0;
    int sum = 0;
    int counter = 0;
    char *encoding = "";
    uint64_t min_score = UINT32_MAX;
    int max_score = 0;
    uint64_t number_of_bases = 0;
    uint64_t total_counts[91] = {0};
    uint64_t averages[500] = {0};

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
        averages[i] = sum/100;
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


    // adjust the x position of the SVG based on whether this is paired end data
    int adjust = (position == 1)*120;
    printf("<g width=\"700\" height=\"700\" xmlns=\"http://www.w3.org/2000/svg\" transform=\"translate(%d 50)\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n", (position == 1)*610);
    if (position == 0) {
        printf("<text x=\"395\" y=\"0\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\">%lu <tspan opacity=\"0.75\"> reads with encoding</tspan> %s <tspan opacity=\"0.75\">in</tspan> forward <tspan opacity=\"0.75\">read</tspan></text>\n", data->number_of_sequences, encoding);
    }
    else if (position == 1) {
        printf("<text x=\"275\" y=\"0\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\">%lu <tspan opacity=\"0.75\"> reads with encoding</tspan> %s <tspan opacity=\"0.75\">in</tspan> reverse <tspan opacity=\"0.75\">read</tspan></text>\n", data->number_of_sequences, encoding);
    }

    // tick marks

    #define tickmarks(x1,x2,y1,y2, stroke_width) printf("<line x1=\"%d\" x2=\"%d\" y1=\"%d\" y2=\"%d\" style=\"stroke:black;stroke-width:%d;stroke-opacity:0.75\" />\n", x1,x2,y1,y2,stroke_width);


        // vertical
    for (i = 1; i < 8; i++) {
        tickmarks(i*450/8+50+120*(position==0), i*450/8+50+120*(position==0), 118, 122, (i+1)%2+1)
        tickmarks(i*450/8+50+120*(position==0), i*450/8+50+120*(position==0), 383, 387, (i+1)%2+1)
        tickmarks(i*450/8+50+120*(position==0), i*450/8+50+120*(position==0), 493, 497, (i+1)%2+1)
    }

        // horizontal
    for (i = 1; i < 8; i++) {
        tickmarks(163+340*(position==1), 167+340*(position==1), i*254/8+125, i*254/8+125, i%2+1)
    }

    // base ratios
    printf("<rect x=\"%d\" y=\"15\" width=\"450\" height=\"100\" style=\"fill:#89bc89;fill-opacity:1.0;stroke-opacity:1.0\" />\n", 170-adjust);
    printf("<svg x=\"%d\" y=\"15\" width=\"450\" height=\"100\" preserveAspectRatio=\"none\" viewBox=\"0 0 %lu 100\">", 170-adjust, data->max_length-1);
    printf("<polyline points=\"0,100 ");
    for (i=0; i<data->max_length; i++) {
        int content = data->bases[i].content[0];
        printf("%d,%d ", i, content);
    }
    printf("%lu, 100 0,100\" style=\"fill:#648964;stroke:none\"/>\n", data->max_length-1);
    printf("<polyline points=\"0,100 ");
    for (i=0; i<data->max_length; i++) {
        int content = data->bases[i].content[0] + data->bases[i].content[3];
        printf("%d,%d ", i, content);
    }
    printf("%lu, 100 0,100\" style=\"fill:#84accf;stroke:none\"/>\n", data->max_length-1);
    printf("<polyline points=\"0,100 ");
    for (i=0; i<data->max_length; i++) {
        int content = data->bases[i].content[0] + data->bases[i].content[3] + data->bases[i].content[1];
        printf("%d,%d ", i, content);
    }
    printf("%lu, 100 0,100\" style=\"fill:#5d7992;stroke:none\"/>\n", data->max_length-1);

    printf("</svg>\n");

    // score distribution
    printf("<rect x=\"%d\" y=\"125\" width=\"100\" height=\"254\" style=\"fill:rgb(245,245,245);fill-opacity:1.0;stroke-opacity:1.0\" />\n", 60+(position == 1)*450);
    printf("<svg x=\"%d\" y=\"125\" width=\"100\" height=\"254\" preserveAspectRatio=\"none\" viewBox=\"0 0 100 %d\">\n", 60+(position == 1)*450, max_score);
    for (i = 0; i < max_score+offset; i++){
        printf("<rect x=\"%lu\" y=\"%d\" width=\"%lu\" height=\"1\" style=\"fill:steelblue;fill-opacity:1\" />", (uint64_t)abs(100*(position == 0)-(100*total_counts[i]/number_of_bases)*(position == 0)), max_score - i -1+offset, (uint64_t)100*total_counts[i]/number_of_bases);
    }
    printf("</svg>\n");

    //length distribution
    printf("<rect x=\"%d\" y=\"390\" width=\"450\" height=\"100\" style=\"fill:rgb(245,245,245);fill-opacity:1.0;stroke-opacity:1.0\" />\n", 170-adjust);
    printf("<svg x=\"%d\" y=\"390\" width=\"450\" height=\"100\" preserveAspectRatio=\"none\" viewBox=\"0 0 %lu 100\">", 170-adjust, data->max_length);
    for (i = 0; i < data->max_length; i++){
        printf("<rect x=\"%d\" y=\"0\" width=\"1\" height=\"%lu\" style=\"stroke:none;fill:steelblue;fill-opacity:1\" />", i, data->bases[i].length_count);
    }
    printf("</svg>\n");

    //kmer distribution
    if (adapters_used == 1) {
        printf("<rect x=\"%d\" y=\"500\" width=\"450\" height=\"100\" style=\"fill:rgb(245,245,245);fill-opacity:1.0;stroke-opacity:1.0\" />\n", 170-adjust);
        printf("<svg x=\"%d\" y=\"500\" width=\"450\" height=\"100\" preserveAspectRatio=\"none\" viewBox=\"0 0 %lu 100\">", 170-adjust, data->max_length);
        for (i = 0; i < data->max_length; i++){
            printf("<rect x=\"%d\" y=\"0\" width=\"1\" height=\"%lu\" style=\"stroke:none;fill:steelblue;fill-opacity:1\" />", i, data->bases[i].kmer_count);
        }
        printf("</svg>\n");
    }

    // heatmap
    printf("<rect x=\"%d\" y=\"125\" width=\"450\" height=\"254\" style=\"fill:rgb(245,245,245);fill-opacity:1.0;stroke-opacity:1.0\" />", 170-adjust);
    printf("<svg x=\"%d\" y=\"125\" width=\"450\" height=\"254\" preserveAspectRatio=\"none\" viewBox=\"0 0 %lu %d\">", 170-adjust, data->max_length, max_score);
    printf("<rect x=\"0\" y=\"0\" width=\"%lu\" height=\"%d\" style=\"stroke:none;stroke-width:1;fill-opacity:0.2;fill:mediumseagreen;stroke-opacity:1.0\" />", data->max_length, max_score-40+12);
    printf("<rect x=\"0\" y=\"%d\" width=\"%lu\" height=\"%d\" style=\"stroke:none;stroke-width:1;fill-opacity:0.2;fill:gold;stroke-opacity:1.0\" />", max_score-40+12, data->max_length, 8);
    printf("<rect x=\"0\" y=\"%d\" width=\"%lu\" height=\"%d\" style=\"stroke:none;stroke-width:1;fill-opacity:0.2;fill:tomato;stroke-opacity:1.0\" />", max_score-20, data->max_length, 20);
    for (i = 0; i < data->max_length; i++) {
        for (j = offset; j < max_score+offset; j++) {
            printf("<rect x=\"%d\" y=\"%d\" width=\"1\" height=\"1\"  style=\"stroke:none;fill:%s;fill-opacity:%f\" />\n", i, max_score-j-1+offset, "#000", (float)data->bases[i].scores[j]/100);
        }
    }

    // mean line
    printf("<polyline points=\" ");
    for (i=0; i< data->max_length; i++) {
        printf("%d,%lu ", i, max_score-averages[i]+offset);
    }
    printf("\" style=\"stroke-width: 0.5; opacity: 0.5; fill:none;stroke:#000\"/>\n");

    printf("</svg>\n");

    // heatmap scale line
    printf("<svg x=\"%d\" y=\"125\" width=\"450\" height=\"254\" preserveAspectRatio=\"none\" viewBox=\"0 0 254 %d\">\n", 41+584*(position == 0), max_score);
    if (position != 1) {
        printf("<text x=\"%d\" y=\"12\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" transform=\"scale(0.5 %f)\">%d</text>\n", 7, (float)(max_score + 1)/254, max_score);
        printf("<text x=\"%d\" y=\"20\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" transform=\"translate(0 %d) scale(0.5 %f)\">28</text>\n", 7, 10+max_score-40, (float)(max_score + 1)/254);
        printf("<text x=\"%d\" y=\"28\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" transform=\"translate(0 %d) scale(0.5 %f)\">20</text>\n", 7, 16+max_score-40, (float)(max_score + 1)/254);
    }
    printf("</svg>\n");

     // labels
    if (position != 1) {
        // score distribution
        printf("<text x=\"-275\" y=\"50\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" transform=\"rotate(-90)\" fill-opacity=\"0.5\">Score</text>");
        printf("<text x=\"60\" y=\"360\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" fill-opacity=\"0.75\">Score</text>");
        printf("<text x=\"60\" y=\"375\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" fill-opacity=\"0.75\">Distribution</text>");

        // base content distribution
        printf("<text x=\"628\" y=\"40\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"#89bc89\" font-weight = \"bold\">%%A</text>");
        printf("<text x=\"628\" y=\"60\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"#648964\" font-weight = \"bold\">%%T</text>");
        printf("<text x=\"628\" y=\"80\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"#84accf\" font-weight = \"bold\">%%G</text>");
        printf("<text x=\"628\" y=\"100\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"#5d7992\" font-weight = \"bold\">%%C</text>");
    }
    else if (position == 1) {
        printf("<text x=\"-275\" y=\"630\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" transform=\"rotate(-90)\" fill-opacity=\"0.5\">Score</text>");
        printf("<text x=\"510\" y=\"360\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" fill-opacity=\"0.75\">Score</text>");
        printf("<text x=\"510\" y=\"375\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" fill-opacity=\"0.75\">Distribution</text>");
    }

    // score distribution
    printf("<text x=\"%d\" y=\"136\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">%d</text>", 30+585*(position==1), max_score);
    printf("<text x=\"%d\" y=\"378\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">1</text>", 30+585*(position==1));

    // base content
    printf("<text x=\"%d\" y=\"110\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"white\" fill-opacity=\"0.75\">Base Content Percentage</text>", 172-adjust);

    // length distribution
    printf("<text x=\"%d\" y=\"486\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" fill-opacity=\"0.75\">Length Distribution</text>", 172-adjust);

    // kmer distribution

    if (adapters_used == 1) {
        printf("<text x=\"%d\" y=\"596\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" fill-opacity=\"0.75\">Adapter Content</text>", 172-adjust);
    }

    // common labels
    printf("<text x=\"%d\" y=\"115\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">100</text>", 60+(position == 1)*530);
    printf("<text x=\"%d\" y=\"24\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">100</text>", 137+(position == 1)*366);
    printf("<text x=\"%d\" y=\"490\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">100</text>", 137+(position == 1)*366);
    printf("<text x=\"%d\" y=\"115\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">0</text>",  154+(position == 1)*356);
    printf("<text x=\"%d\" y=\"400\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">0</text>",  147+(position == 1)*360);
    printf("<text x=\"-465\" y=\"%d\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" transform=\"rotate(-90)\" fill-opacity=\"0.5\">Percent</text>", 157+(position == 1)*360);
    printf("<text x=\"-90\" y=\"%d\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" transform=\"rotate(-90)\" fill-opacity=\"0.5\">Percent</text>", 157+(position == 1)*360);
    printf("<text x=\"%d\" y=\"115\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" fill-opacity=\"0.5\">Percent</text>", 90+(position == 1)*437);


    if (adapters_used == 1) {
        printf("<text x=\"%d\" y=\"617\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">1</text>", 170-adjust);
        printf("<text x=\"%d\" y=\"617\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">%lu</text>", 600-adjust, data->original_max_length);
        printf("<text x=\"%d\" y=\"620\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" fill-opacity=\"0.5\">Base Pairs</text>", 350-adjust);
        printf("<text x=\"%d\" y=\"510\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">0</text>",  147+(position == 1)*360);
        printf("<text x=\"-575\" y=\"%d\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" transform=\"rotate(-90)\" fill-opacity=\"0.5\">Percent</text>", 157+(position == 1)*360);
        printf("<text x=\"%d\" y=\"600\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">100</text>", 137+(position == 1)*366);
    }
    else {
        printf("<text x=\"%d\" y=\"517\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">1</text>", 170-adjust);
        printf("<text x=\"%d\" y=\"517\" font-family=\"sans-serif\" font-size=\"12px\" fill=\"black\" fill-opacity=\"0.5\">%lu</text>", 600-adjust, data->original_max_length);
        printf("<text x=\"%d\" y=\"520\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" fill-opacity=\"0.5\">Base Pairs</text>", 350-adjust);
    }

    // encoding
    printf("<text x=\"%d\" y=\"375\" font-family=\"sans-serif\" font-size=\"15px\" fill=\"black\" fill-opacity=\"0.75\">Per base sequence quality</text>", 172-adjust);

    printf("</g>\n");
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
    if((paired ^ unpaired) == 0){
      printf("%s\n", "Usage: quack [OPTION...]\nTry `quack --help' or `quack --usage' for more information.");
      exit(1);
    }

    if(adapters) kmers = read_adapters(arguments.adapters);

    width  = (paired)?1230:700;
    height = (adapters)?700:600;
    
    printf("<svg width='%d' height='%d' viewBox='20 0 %d %d' " 
           "xmlns='http://www.w3.org/2000/svg' " 
           "xmlns:xlink='http://www.w3.org/1999/xlink'>\n", 
           width, height, width, height);

    // If name is given, add to middle of viewBox (half of width + min-x of viewbox)
    if(arguments.name != NULL)
      printf("<text x='%d' y='25' font-family='sans-serif' " 
             "text-anchor='middle' font-size='30px' " 
             "fill='black'>%s</text>\n", width/2 + 20, arguments.name);
    
    sequence_data *data = read_fastq(((paired)?arguments.forward:arguments.unpaired), kmers);
    sequence_data *transformed_data = transform(data);
    draw(transformed_data, 0, 1);
    free(data);
    
    if(paired){
      data = read_fastq(arguments.reverse, kmers);
      transformed_data = transform(data);
      draw(transformed_data, 1, 1);
      free(data);
    }

    printf("</svg>\n");
    
    exit (0);
}
