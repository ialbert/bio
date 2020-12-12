.PHONY: dist build test docs

# Location of the documentation branch.
DOCBRANCH=../bio-docs

all: serve

# Run the tests.
test:
	pytest

# Generate the test data for the tests.
build_data:
	(cd test/data && bash ../bio_examples.sh)

# Generate test from the example script.
build_test:
	(cd test && python generate.py)
	pytest

# Generate the docs.
docs:
	(cd docs && Rscript -e "bookdown::render_book(input='index.txt', output_dir='.html', output_format='bookdown::gitbook')")

# Synchronize the docs.
sync:
	# Get the curent docs.
	(cd ${DOCBRANCH} && git pull origin gh-pages)
	# Synchronize to documentation branch.
	rsync -avz docs/.html/* ${DOCBRANCH}
	# Commit and push out changes.
	(cd ${DOCBRANCH} && git commit -am 'updated the documentation' && git push origin gh-pages)

# Serve the documentation as a webpage.
serve:
	Rscript -e "bookdown::serve_book(dir='docs', preview=TRUE, output_dir='.html', port=8000)"

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
upload: test build
	rm -rf dist
	python setup.py sdist bdist_wheel
	#python -m twine upload --repository testpypi dist/*
	python -m twine upload --repository pypi dist/*

# Uploads prebuilt data to Google Cloud
upload_prebuilt:
	bash docs/upload_prebuilt_data.sh





