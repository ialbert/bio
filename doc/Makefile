all: serve

publish: build sync

init:
	rm -rf _book

serve: init
	Rscript -e "bookdown::serve_book(dir='.', preview=TRUE, output_dir='_book', port=8000)"

push:
	-git commit -am `whoami`
	-git push
