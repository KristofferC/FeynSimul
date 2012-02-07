from setuptools import setup, find_packages

setup(name='FeynSimul',
      version = '0.4',
      packages = find_packages(),
      package_data = {'FeynSimul':["GPUsrc/kernel.c"]},
      include_package_data = True,
      )
