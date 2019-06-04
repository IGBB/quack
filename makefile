src = $(wildcard src/*.c)
obj = $(src:.c=.o)

override LDFLAGS := -lz -lm $(LDFLAGS)
override CFLAGS := -Iklib -O3 $(CFLAGS)

all : klib/kseq.h quack

quack: $(obj)
	$(CC) -o $@ $^ $(LDFLAGS)

klib/kseq.h:
	git submodule update --init --recursive

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY: all clean images test
clean:
	rm -f $(obj) quack

images: quack
	$(MAKE) -C images all

test: images
