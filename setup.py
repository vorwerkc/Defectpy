from setuptools import setup
from setuptools import find_packages

setup(name='defectpy',
      version='0.1',
      packages=find_packages(),
      description='Utilities for Defect Calculations',
      url='https://github.com/vorwerkc/defectpy.git',
      author='Christian Vorwerk',
      license='GPLv3',
      install_requires=[
          'numpy',
          'scipy'
      ],
      python_requires='>=2.7, >=3.0, !=3.0.*, !=3.1.*, !=3.2.*, <4',
      zip_safe=True)
