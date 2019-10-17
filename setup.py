from setuptools import setup

setup(name='gpnet',
      version='0.1',
      description='Signomial programming compatible models '
                  'for network design.',
      url='https://github.com/1ozturkbe/EGPAnet',
      author='Berk Ozturk',
      author_email='bozturk@mit.edu',
      license='MIT',
      packages=['gpnet'],
      install_requires=[
          'gpkit',
          'robust'
      ])
