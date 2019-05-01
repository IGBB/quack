# quack

A FASTQ quality assessment tool

## Citation
A. Thrash, M. Arick, and D. G. Peterson, “Quack: A quality assurance tool for high throughput sequence data,” Analytical Biochemistry, vol. 548, pp. 38–43, 2018. https://doi.org/10.1016/j.ab.2018.01.028

## Latest Release

The latest release of quack and its binaries can always be found [here](https://github.com/IGBB/quack/releases/latest).

## Dependencies

* zlib
* klib (pulled by the submodule update below)

## Installation from Source

git clone https://github.com/IGBB/quack.git

cd quack/

make && make test

## Binaries

Binaries are available in the bin/ folder. Current testing of these binaries has been limited. If a binary doesn't work, try compiling from source on your system.

## Running Quack

Quack has the following options.

```
  -1, --forward     forward strand data in gzipped FASTQ format, must be used with -2 or --reverse
  -2, --reverse     reverse strand data in gzipped FASTQ format, must be used with -1 or --forward
  -a, --adapters    adapters in gzipped FASTA format (optional)
  -n, --name    a descriptive name to be printed with the output image (optional)
  -u, --unpaired    unpaired data in gzipped FASTQ format
  -?, --help, --usage   prints the help or usage information
  -V, --version prints the program version
```

Quack takes gzipped FASTQ-formatted files as input for data and gzipped As output, quack prints an SVG formatted image to standard output.


### Examples

#### Paired-end with name and adapters
`quack -1 reads.1.fastq.gz -2 reads.2.fastq.gz -n sample_name -a adapters_files.fasta.gz > sample_name.svg`

#### Unpaired with name and adapters
`quack -u reads.fastq.gz -n sample_name -a adapters.fa.gz > sample_name.svg`

#### Paired-end without name and adapters
`quack -1 reads.1.fastq.gz -2 reads.2.fastq.gz > sample_name.svg`

#### Unpaired without name and adapters
`quack -u reads.fastq.gz > sample_name.svg`

### Output

Quack is capable of producing output for single-ended data and paired-end data. Only the singled-ended data is labeled, since the paried-end data has all the same parts.

#### Single-ended Data
<img src="https://ars.els-cdn.com/content/image/1-s2.0-S0003269718300630-gr4.jpg">


A. The base content distribution showing the percentage of each nucleotide in each column of an array.  
B. A heatmap showing the distribution of sequence quality for each column and a line representing mean quality scores across the array  
C. A score distribution graph showing the percentage of bases matching certain scores, with 100% on the left of the graph and 0% on the right. The highest scoring data appears at the top of the graph.  
D. Length distribution graph showing the percentage of reads of a given length  
E. Adapter content distribution graph showing how adapter content is distributed throughout an array  

#### Paired-end Data
<img src="https://ars.els-cdn.com/content/image/1-s2.0-S0003269718300630-gr5.jpg">
