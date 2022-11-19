#!/usr/bin/env python
"""
# Author: Wenze Huang
# Created Time : Mon 31 Oct 2022 09:42:31 PM CST
# File Name: setup.py
# Description:
"""

from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='TISnet',
      version='0.0.1',
      description='TISnet',
      packages=find_packages(),

      author='Wenze Huang',
      author_email='hwz16@tsinghua.org.cn',
      url='https://github.com/huangwenze/TISnet',
      install_requires=requirements,
      python_requires='>3.6.0',

      classifiers=[
          'Development Status :: 1 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.6',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Operating System :: CentOS :: Linux',
     ],
     )
