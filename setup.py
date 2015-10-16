from setuptools import setup, Extension, find_packages, Command
from setuptools.command.develop import develop
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
            version='0.121',
            description='A bioinformatics sequence library',
            url='https://github.com/jrellis/pybio',
            author='Rob Ellis',
            author_email='jrellis@fas.harvard.edu',
            packages=['pybio', 'pybio.alignment'],
            ext_modules = [
                Extension('pybio.alignment.smith_waterman',sources=['pybio/alignment/smith_waterman.c', 'pybio/alignment/lib/ssw.c'],depends=['pybio/alignment/lib/ssw.h'])
                ],
            cmdclass = {'cythonize':cythonize, 'cython_develop':cython_develop},
            install_requires=['numpy', 'setuptools'],
            )

    setup(**metadata)

if __name__=='__main__':
    setup_package()
