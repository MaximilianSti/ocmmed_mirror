# -*- coding: utf-8 -*-


from __future__ import absolute_import, print_function

from setuptools import setup, find_packages
import sys

requirements = [
        'dexom-python']

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='ocmmed',
    version='1.1.4',
    packages=find_packages('.'),
    install_requires=requirements,
    include_package_data=True,
    author='Maximilian Stingl',
    author_email='contact-metexplore@inrae.fr',
    description='Obtaining cell-specific metabolic models through enumeration with DEXOM',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://forgemia.inra.fr/metexplore/cbm/ocmmed',
    python_requires=">=3.7,<3.10",
)
sys.path.append('ocmmed')
