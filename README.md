# quack

A FASTQ quality assessment tool

## Dependencies

* zlib
* klib (pulled by the submodule update below)

## Installation from Source

git clone https://github.com/IGBB/quack.git

cd quack/

git submodule update --init --recursive

gcc quack.c -o quack -lz -lm -I klib -O3

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
