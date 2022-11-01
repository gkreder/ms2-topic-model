from gettext import install
from setuptools import setup
import os
import sys

lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = os.path.join(lib_folder, "requirements.txt")
install_requires = []
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = [x for x in f.read().splitlines() if not x.startswith('#') and not x.startswith("_")]


setup(name='ms2-topic-model',
      version='1.0',
      description='Code base for LLDA based MS2 substructure prediction',
      url='http://github.com/gkreder/ms2-topic-model',
      author='Gabriel Reder',
      author_email='gkreder@gmail.com',
      license='MIT',
      packages=['ms2-topic-model'],
      install_requires = install_requires)