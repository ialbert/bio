SOURCES = $(wildcard *.sh)
TARGETS = $(SOURCES:.sh=.txt)

%.txt: %.sh
	python build.py $< > $@

all: $(TARGETS)

push:
	git commit -am "docs updated by `whoami`"
	git push

# Generate the docs.
docs:
	(cd .. && Rscript -e "bookdown::render_book(input='index.txt', output_dir='.book', output_format='bookdown::gitbook')")
