#!/usr/bin/env python
from setuptools import setup, find_packages
from PTAS.__init__ import __version__

setup(name='PTAS',
      version=__version__,
      description='Pipeline for Targeted Amplicon Sequencing.',
      author='Koki Chigira',
      author_email='s211905s@st.go.tuat.ac.jp',
      url='https://github.com/KChigira/PTAS/',
      license='MIT',
      packages=find_packages(),
      install_requires=[
        'pandas',
        'matplotlib',
      ],
      entry_points={'console_scripts': [
            'mkvcf = PTAS.mkvcf:main',
            'mkprimer = PTAS.mkprimer:main',
            'mkselect = PTAS.mkselect:main',
            'mkbind = PTAS.mkbind:main',
            ]
      }
    )
