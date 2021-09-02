import biorun
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="bio",
    version=biorun.VERSION,
    author="Istvan Albert",
    author_email="istvan.albert@gmail.com",
    description="bio",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ialbert/bio",
    packages=find_packages(include=["biorun", "biorun.*"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],

    install_requires=[
        'biopython >= 1.79',
        'requests',
        'tqdm>=4.62',
        'plac',
        'mygene',
    ],

    entry_points={
        'console_scripts': [
            'bio=biorun.__main__:run',
            'comm.py=biorun.scripts.comm:run',
            'uniq.py=biorun.scripts.uniq:run',
            'fasta_filter.py=biorun.scripts.fasta_filter:run',
        ],
    },

    python_requires='>=3.6',

)
