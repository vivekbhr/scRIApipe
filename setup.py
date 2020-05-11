#!/usr/bin/env python
from setuptools import setup, find_packages
from glob import glob
import os.path

# Set __version__ to keep __init__.py non-empty
exec(open('scriapipe/__init__.py').read())
requires = open("requirements.txt").read().strip().split("\n")

setup(
    name='scRIApipe',
    version=__version__,
    description='A workflow for identifying differantial transcript usage',
    url='https://github.com/vivekbhr/scRIApipe',
    author='Vivek Bhardwaj and Thijs Makaske',
    license='MIT',
    keywords='splicing isoform kallisto',
    packages=find_packages(),
    scripts=['bin/scRIA'],
    install_requires=requires,
    include_package_data=True,
    data_files=[("", ["LICENSE"])]
)
