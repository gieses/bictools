#!/usr/bin/env python

import sys
import os.path

from setuptools import setup, find_packages
#This is a list of files to install, and where
#(relative to the 'root' dir, where setup.py is)
#You could be more specific.



setup(name = "bictools",
    version = "0.1dev",
    description = "bictools - a python package to perform variant detection on \
proteomic MS data using BICEPS (Renard et al. 2012)",
    author = "Sven H. Giese",
    author_email = "sven.giese@tu-berlin.de",
    url = "",
    #Name the folder where your packages live:
    #(If you have other packages (dirs) or modules (py files) then
    #put them into the package directory - they will be found
    #recursively.)
    requires = [ 'HTSeq', 'numpy', 'pyopenms', 'matplotlib', 'pandas', 'brewer2mpl'],
    packages=find_packages(),
    #'package' package must contain files (see list above)
    #I called the package 'package' thus cleverly confusing the whole issue...
    #This dict maps the package name =to=> directories
    #It says, package *needs* these files.
    #'runner' is in the root.
    long_description = open('README').read(),
    py_modules = [ 
	  'bictools.scripts.bic_CSV',
	  'bictools.scripts.bic_DB',
	  'bictools.scripts.bic_evaluate'
	],
    scripts = [
          'scripts/bictools-bic_CSV',
          'scripts/bictools-bic_DB',
          'scripts/bictools-bic_evaluate']
    #
    #This next part it for the Cheese Shop, look a little down the page.
    #classifiers = []
)






