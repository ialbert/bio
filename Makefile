.PHONY: dist build test docs

# Location of the documentation branch.
REMOTE = www@bioinfo.help:~/sites/bioinfo.help/

all: serve

# Run the tests.
test:
	@python src/biorun/test.py

# Update the usage data.
testdata:
	(cd test/data && bash ../usage.sh)


# Update the data for the test script.
fulltest: testdata test

# Generate the output for docs
code:
	(cd docs/code && bash code.sh)

# Generate the docs.
docs:
	(cd docs && Rscript -e "bookdown::render_book(input='index.Rmd', output_dir='.book', output_format='bookdown::gitbook')")

# Push out the docs to remote docs.
sync: docs
	rsync -avz  docs/.book/ ${REMOTE}
	rsync -avz docs/*.md ~/book/biostar-handbook-2/books/bio/
	rsync -avz docs/images/igv-index.png ~/book/biostar-handbook-2/books/main/images

# Serve the documentation as a webpage.
serve:
	(cd docs && rm -rf .book)
	Rscript -e "bookdown::serve_book(dir='docs', preview=TRUE, output_dir='.book', port=8000)"

# Clean the files.
clean:
	rm -rf dist build bio.egg-info

# Quick push for small changes.
push:
	git commit -am "commit by `whoami`"
	git push

# Build Python package.
build:
	rm -rf build dist
	hatch build
	ls -lh dist

# Publish the package
publish: test build
	hatch publish
	
#REMOTE=www@biostarhandbook.com:/home/www/book/data_www/bio

# Upload prebuilt data to distribution site.
upload:
	(cd ~/.bio && GZIP=-9 && tar czvf biodata.tar.gz taxdump.tar.gz *.json *.sqlite assembly_summary_genbank.txt)
	rsync -avz --progress ~/.bio/biodata.tar.gz ${REMOTE}data


