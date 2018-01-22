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

These binaries are also available in the [latest release](https://github.com/IGBB/quack/releases/latest).
