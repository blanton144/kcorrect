# encoding: utf-8
#
# setup.py
#


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from setuptools import setup, find_packages

import os
import argparse
import sys


NAME = 'kcorrect'
VERSION = 'dev'
RELEASE = 'dev' in VERSION


def run(packages, install_requires):

    setup(name=NAME,
          version=VERSION,
          license='BSD3',
          description='Description of your project.',
          long_description=open('README.rst').read(),
          author='Michael Blanton',
          author_email='michael.blanton@nyu.edu',
          keywords='astronomy software',
          url='https://github.com/blanton144/kcorrect',
          include_package_data=True,
          packages=packages,
          install_requires=install_requires,
          package_dir={'': 'python'},
          classifiers=[
              'Development Status :: 4 - Beta',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: BSD License',
              'Natural Language :: English',
              'Operating System :: OS Independent',
              'Programming Language :: Python',
              'Programming Language :: Python :: 3.9',
              'Topic :: Documentation :: Sphinx',
              'Topic :: Software Development :: Libraries :: Python Modules',
          ],
          )


def get_requirements(opts):
    ''' Get the proper requirements file based on the optional argument '''

    if opts.dev:
        name = 'requirements.txt'
    elif opts.doc:
        name = 'requirements_doc.txt'
    else:
        name = 'requirements.txt'

    requirements_file = os.path.join(os.path.dirname(__file__), name)
    install_requires = [line.strip().replace('==', '>=') for line in open(requirements_file)
                        if not line.strip().startswith('#') and line.strip() != '']
    return install_requires


def remove_args(parser):
    ''' Remove custom arguments from the parser '''

    arguments = []
    for action in list(parser._get_optional_actions()):
        if '--help' not in action.option_strings:
            arguments += action.option_strings

    for arg in arguments:
        if arg in sys.argv:
            sys.argv.remove(arg)


if __name__ == '__main__':

    # Custom parser to decide whether which requirements to install
    parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]))
    parser.add_argument('-d', '--dev', dest='dev', default=False, action='store_true',
                        help='Install all packages for development')
    parser.add_argument('-o', '--doc', dest='doc', default=False, action='store_true',
                        help='Install only core + documentation packages')

    # We use parse_known_args because we want to leave the remaining args for distutils
    args = parser.parse_known_args()[0]

    # Get the proper requirements file
    install_requires = get_requirements(args)

    # Now we remove all our custom arguments to make sure they don't interfere with distutils
    remove_args(parser)

    # Have distutils find the packages
    packages = find_packages(where='python')

    # Runs distutils
    run(packages, install_requires)
