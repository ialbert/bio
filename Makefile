.PHONY: dist build test docs

# Location of the documentation branch.
REMOTE = www@bioinfo.help:~/sites/bioinfo.help/

all: serve

# Run the tests.
test:
	pytest -x

# Update the usage data.
testdata:
	(cd test/data && bash ../usage.sh)


# Generate test from the example script.
generate:
	(cd test && python generate.py)
	pytest -x

# Update the data for the test script.
fulltest: testdata generate

# Generate the docs.
docs:
	(cd docs && Rscript -e "bookdown::render_book(input='index.Rmd', output_dir='.book', output_format='bookdown::gitbook')")

# Push out the docs to remote docs.
sync:
	rsync -avz  docs/.book/ ${REMOTE}
	rsync -avz docs/bio-explain.md ~/book/biostar-handbook-2/books/main/doc/ontology/
	rsync -avz docs/bio-search.md ~/book/biostar-handbook-2/books/main/doc/data/

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
pypi: test build
	rm -rf dist
	python setup.py sdist bdist_wheel
	#python -m twine upload --repository testpypi dist/*
	python -m twine upload --repository pypi dist/*


#REMOTE=www@biostarhandbook.com:/home/www/book/data_www/bio

# Upload prebuilt data to distribution site.
upload:
	(cd ~/.bio && GZIP=-9 && tar czvf biodata.tar.gz taxdump.tar.gz *.json *.sqlite)
	rsync -avz --progress ~/.bio/biodata.tar.gz ${REMOTE}data


