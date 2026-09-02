# The bibliography is inline, so a single latexmk run suffices.
.PHONY: all clean

all: main.pdf

main.pdf: main.tex
	latexmk -pdf main.tex

clean:
	latexmk -C
