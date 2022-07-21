# -*- coding: utf-8 -*-


from __future__ import absolute_import, print_function

from setuptools import setup, find_packages
import sys

requirements = [
        'dexom-python',
        'miom[full]']

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='ocmmed',
    version='0.1',
    packages=find_packages('.'),
    install_requires=requirements,
    include_package_data=True,
    author='Maximilian Stingl',
    author_email='maximilian.h.a.stingl@gmail.com',
    description='Obtaining cell-specific metabolic models through enumeration with DEXOM',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://forgemia.inra.fr/metexplore/cbm/ocmmed',
    python_requires=">=3.7",
)
sys.path.append('ocmmed')
