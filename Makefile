.PHONY: dist build test docs

# Location of the documentation branch.
REMOTE=www@bioinfo.help:~/sites/bioinfo.help/

all: serve

# Run the tests.
test:
	pytest

# Generate the test data for the tests.
build_data:
	(cd test/data && bash ../bio-examples.sh)

# Generate test from the example script.
build_test:
	(cd test && python generate.py)
	pytest --ff -x

# Generate the docs.
docs:
	(cd docs && Rscript -e "bookdown::render_book(input='index.txt', output_dir='.book', output_format='bookdown::gitbook')")

# Push out the docs to remote docs.
sync:
	rsync -avz docs/.book/* ${REMOTE}

# Serve the documentation as a webpage.
serve:
	(cd docs && rm -rf .book)
	Rscript -e "bookdown::serve_book(dir='docs', preview=TRUE, output_dir='.book', port=8000)"

# Clean the files.
clean:
	rm -rf dist build bio.egg-info

# Quick push for small changes.
push:
	git commit -am "docs updated by `whoami`"
	git push

# Build Python package.
build:
	python setup.py sdist bdist_wheel

# Upload new version to PyPI.
publish: test build
	rm -rf dist
	python setup.py sdist bdist_wheel
	#python -m twine upload --repository testpypi dist/*
	python -m twine upload --repository pypi dist/*


REMOTE=www@biostarhandbook.com:/home/www/book/data_www/bio

# Upload prebuilt data to distribution site.
upload:
	rsync -avz --progress ~/.bio/taxdb.json ${REMOTE}
	rsync -avz --progress ~/.bio/taxdb.sqlite ${REMOTE}





