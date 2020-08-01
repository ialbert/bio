.PHONY: dist build

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




