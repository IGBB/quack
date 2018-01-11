# quack

A FASTQ quality assessment tool

## Dependencies

* zlib
* klib (pulled by the submodule update below)

## Installation

git clone https://github.com/IGBB/quack.git

cd quack/

git submodule update --init --recursive

gcc quack.c -o quack -lz -lm -I klib -O3
