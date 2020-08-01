import setuptools
import biorun

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bio",
    version=biorun.VERSION,
    author="Istvan Albert",
    author_email="istvan.albert@gmail.com",
    description="bio",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/ialbert/bio",
    packages=["biorun"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'plac',
        'attrs',
        'biopython',
        'intervaltree',
    ],

    python_requires='>=3.6',

)
