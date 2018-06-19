from setuptools import setup, find_packages

setup(name='simtools',
      version='0.1',
      description='Tools for simulating associations with genetic data',
      url='http://github.com/rmporsch/simtools',
      author='Robert Porsch',
      author_email='rmporsch@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'pandas',
          'numpy'
          ],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
