from setuptools import setup, Extension, find_packages, Command
from setuptools.command.develop import develop
import numpy
import os

cython_globs = ['pybio/alignment/*.pyx']

class cythonize(Command):
    description = 'Builds cython files.'
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        from Cython.Build import cythonize
        cythonize(cython_globs)

class cython_develop(develop):
    description = 'Builds cython files and runs develop.'
    def run(self):
        from Cython.Build import cythonize
        cythonize(cython_globs)
        develop.run(self)

def setup_package():

    metadata = dict(
            name='pybio',
            version='0.128',
            description='A bioinformatics sequence library',
            url='https://github.com/jrellis/pybio',
            author='Rob Ellis',
            author_email='jrellis@fas.harvard.edu',
            packages=['pybio', 'pybio.alignment'],
            ext_modules = [
                Extension('pybio.alignment.smith_waterman',sources=['pybio/alignment/smith_waterman.c', 'pybio/alignment/lib/ssw.c'],depends=['pybio/alignment/lib/ssw.h'], include_dirs=[numpy.get_include()])
                ],
            cmdclass = {'cythonize':cythonize, 'cython_develop':cython_develop},
            install_requires=['numpy', 'requests', 'setuptools'],
            setup_requires=['pytest-runner'],
            tests_require=['pytest'],
            )

    setup(**metadata)

if __name__=='__main__':
    setup_package()
