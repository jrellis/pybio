from setuptools import setup, Extension, find_packages
import os

def setup_package():

    metadata = dict(
            name='pybio',
            version='0.11',
            description='A bioinformatics sequence library',
            url='https://github.com/jrellis/pybio',
            author='Rob Ellis',
            author_email='jrellis@fas.harvard.edu',
            packages=['pybio'],
            ext_modules = [Extension('hello',['pybio/hello.c'])],
            )

    cwd = os.path.abspath(os.path.dirname(__file__))
    print os.path.join(cwd, 'PKG-INFO')
    if not os.path.exists(os.path.join(cwd, 'PKG-INFO')):
        from Cython.Build import cythonize
        print('Cythonizing...')
        metadata['ext_modules'] = cythonize([Extension("*", ["pybio/alignment/*.pyx","pybio/alignment/lib/ssw.c"])])

    setup(**metadata)

if __name__=='__main__':
    setup_package()
