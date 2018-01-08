# quack

A FASTQ quality assessment tool

## Installation


### Get quack

git clone https://github.com/IGBB/quack.git

cd quack/

git submodule update --init --recursive


### Linux (most varieties)

gcc quack.c -o quack -lz -lm -I klib -O3

### MacOS (tested on 10.11.6)

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew install argp-standalone

brew install gcc

/usr/local/Cellar/gcc/7.2.0/bin/gcc-7 quack.c -o quack -lz -lm -I klib -I /usr/local/include/ -O3 -L /usr/local/lib/ -largp
