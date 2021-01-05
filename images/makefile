
SRR1168757.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR116/007/SRR1168757/$@

ERR1438847_%.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/007/ERR1438847/$@

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
     explanation.png
