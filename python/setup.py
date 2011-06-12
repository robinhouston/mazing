#!/usr/bin/env python

from distutils.core import setup, Extension

setup(
  name = 'mazing',
  version = '1',
  author = 'Robin Houston',
  author_email = 'robin@boskent.com',
  url = 'http://github.com/robinhouston/mazing',
  description = 'Random access to mazes',
  ext_modules = [Extension(
    "mazing",
    sources = ["mazingmodule.c"],
    libraries = ["gmp", "mazing"],
  )],
  classifiers = [
    "Development Status :: 4 - Beta",
    "License :: OSI Approved :: BSD License",
    "Programming Language :: Python",
    "Programming Language :: C",
  ],
)
