all: release

clean:
	rm -f jpeg_encoder

release: clean jpeg_encoder.c
	gcc -ggdb -O0 -lm -Wall jpeg_encoder.c -o jpeg_encoder

info: clean jpeg_encoder.c
	gcc -ggdb -DINFO -O0 -lm -Wall jpeg_encoder.c -o jpeg_encoder

debug: clean jpeg_encoder.c
	gcc -ggdb -DDEBUG -O0 -lm -Wall jpeg_encoder.c -o jpeg_encoder
