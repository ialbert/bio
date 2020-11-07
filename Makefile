.PHONY: dist build test docs

all: serve

publish: build sync

init:
	rm -rf _book

test:
	pytest

build_data:
	(cd test/data && bash ../test_bio_data.sh)

build_test:
	(cd test && python generate.py)
	pytest

serve: init
	rm -rf doc
	Rscript -e "bookdown::serve_book(dir='doc', preview=TRUE, output_dir='doc/_book', port=8000)"

clean:
	rm -rf dist build bio.egg-info

push:
	git commit -am "code update by `whoami`"
	git push

build:
	python setup.py sdist bdist_wheel

upload: build
	rm -rf dist
	python setup.py sdist bdist_wheel
	#python -m twine upload --repository testpypi dist/*
	python -m twine upload --repository pypi dist/*




