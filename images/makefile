
SRR1168757.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/007/SRR1168757/$@

ERR1438847_%.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/007/ERR1438847/$@

long.fastq.gz:
	wget -O $@.tmp ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR605/007/SRR6059707/SRR6059707_1.fastq.gz
	mv $@.tmp $@

long.adapter.svg long.adapter.txt: long.adapter
long.adapter: long.fastq.gz ../quack ../all.fa.gz
long.noadapt.svg long.noadapt.txt: long.noadapt
long.noadapt: long.fastq.gz ../quack

long.adapter long.noadapt: long.%:
	../quack -u $< $(if $(word 3,$^), -a $(word 3,$^)) -n SRR6059707 --svg $@.svg --txt $@.txt
	touch $@


single.adapter.svg single.adapter.txt: single.adapter
single.adapter: SRR1168757.fastq.gz ../quack
	../quack -u $< -a ../all.fa.gz -n SRR1168757 --svg $@.svg --txt $@.txt
	touch $@

single.noadapt.svg : ../quack SRR1168757.fastq.gz
	../quack -u SRR1168757.fastq.gz -n SRR1168757 > $@

paired.adapter.svg : ../quack ERR1438847_1.fastq.gz ERR1438847_2.fastq.gz
	../quack -1 ERR1438847_1.fastq.gz -2 ERR1438847_2.fastq.gz -a ../all.fa.gz -n ERR1438847 > $@

paired.noadapt.svg : ../quack ERR1438847_1.fastq.gz ERR1438847_2.fastq.gz
	../quack -1 ERR1438847_1.fastq.gz -2 ERR1438847_2.fastq.gz -n ERR1438847 > $@

explanation.overlay.png : explanation.overlay.svg
	convert -background none $< $@

explanation.png :  explanation.overlay.png single.adapter.png
	composite -gravity center $^ $@

%.png: %.svg
	convert $< $@

.PHONY: clean all
clean:
	rm -f $(obj) quack ERR1438847_1.fastq.gz ERR1438847_2.fastq.gz

all: paired.adapter.png \
     paired.noadapt.png \
     single.adapter.png \
     single.noadapt.png \
	 long.adapter.png \
	 long.noadapt.png \
     explanation.png
