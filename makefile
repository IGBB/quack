src = $(wildcard *.c)
obj = $(src:.c=.o)

LDFLAGS = -lz -lm
CFLAGS = -Iklib -O3

quack: $(obj)
	$(CC) -o $@ $^ $(LDFLAGS)

klib/kseq.h:
	git submodule update --init --recursive

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

%.png: %.svg
	convert $< $@

ERR1438847_%.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/007/ERR1438847/$@

test.sing.adapter.svg : quack ERR1438847_1.fastq.gz
	./quack -u ERR1438847_1.fastq.gz -a all.fa.gz -n ERR1438847 > $@

test.sing.noadapt.svg : quack ERR1438847_1.fastq.gz
	./quack -u ERR1438847_1.fastq.gz -n ERR1438847 > $@

test.pair.adapter.svg : quack ERR1438847_1.fastq.gz ERR1438847_2.fastq.gz
	./quack -1 ERR1438847_1.fastq.gz -2 ERR1438847_2.fastq.gz -a all.fa.gz -n ERR1438847 > $@

test.pair.noadapt.svg : quack ERR1438847_1.fastq.gz ERR1438847_2.fastq.gz
	./quack -1 ERR1438847_1.fastq.gz -2 ERR1438847_2.fastq.gz -n ERR1438847 > $@

.PHONY: clean test
clean:
	rm -f $(obj) quack ERR1438847_1.fastq.gz ERR1438847_2.fastq.gz

test: test.pair.adapter.png \
      test.pair.noadapt.png \
      test.sing.adapter.png \
      test.sing.noadapt.png

