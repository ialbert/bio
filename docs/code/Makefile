SOURCES = $(wildcard *.sh)
TARGETS = $(SOURCES:.sh=.txt)

%.txt: %.sh
	bash $< > $@

all: $(TARGETS)

push:
	git commit -am "docs updated by `whoami`"
	git push

