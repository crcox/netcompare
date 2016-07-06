import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "netcompare",
    version = "0.0.1",
    author = "Chris Cox",
    author_email = "cox.crc@gmail.com",
    description = ("Test if two groups of brain networks are significantly non-overlapping."),
    license = "MIT",
    keywords = "fMRI network RSA",
    scripts = ['netcompare'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
    ],
)
