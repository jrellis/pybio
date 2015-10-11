from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize
import os

def setup_package():

    metadata = dict(
            name='pybio',
            version='0.116',
            description='A bioinformatics sequence library',
            url='https://github.com/jrellis/pybio',
            author='Rob Ellis',
            author_email='jrellis@fas.harvard.edu',
            packages=['pybio'],
            ext_modules = cythonize([Extension("*", ["pybio/alignment/*.pyx","pybio/alignment/lib/ssw.c"])]),
            install_requires=['cython', 'numpy'],
            )

    setup(**metadata)

if __name__=='__main__':
    setup_package()
