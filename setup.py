from setuptools import setup

setup(
    name='set_perm',
    version='0.0.1',
    description='genomic set tests',
    py_modules=['set_perm.py'],
    scripts=['bin/set_perm'],
    package_dir={'': 'bin'},
    url='https://github.com/joshuamschmidt/set_perm',
    license='MIT',
    author='joshua schmidt',
    author_email='joshmschmidt1@gmail.com'
)
