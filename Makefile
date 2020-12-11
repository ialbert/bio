.PHONY: dist build test docs

DOCBRANCH=../bio-docs

all: serve

publish: build sync

test:
	pytest

build_data:
	(cd test/data && bash ../bio_examples.sh)

build_test:
	(cd test && python generate.py)
	pytest

docs:
	(cd docs && Rscript -e "bookdown::render_book(input='index.txt', output_dir='.html', output_format='bookdown::gitbook')")
	# Get the curent docs.
	(cd ${DOCBRANCH} && git pull origin gh-pages)
	# Synchronize to documentation branch.
	rsync -avz docs/.html/* ${DOCBRANCH}
	# Commit and push out changes.
	(cd ${DOCBRANCH} && git commit -am 'updated the documentation' && git push origin gh-pages)


serve:
	Rscript -e "bookdown::serve_book(dir='docs', preview=TRUE, output_dir='.html', port=8000)"

clean:
	rm -rf dist build bio.egg-info

push:
	git commit -am "code update by `whoami`"
	git push

build:
	python setup.py sdist bdist_wheel

upload: test build
	rm -rf dist
	python setup.py sdist bdist_wheel
	#python -m twine upload --repository testpypi dist/*
	python -m twine upload --repository pypi dist/*


# Uploads prebuilt data to Google Cloud
upload_prebuilt:
	bash docs/upload_prebuilt_data.sh





