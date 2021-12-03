from setuptools import find_packages, setup

from simplepso import __version__

with open('README.md') as f:
    readme = f.read()


setup(
    name='simplepso',
    packages=find_packages(),
    version=__version__,
    description='Simple usage particle swarm optimization',
    author='James C. Pino',
    author_email='james.ch.pino@gmail.com',
    url='https://github.com/LoLab-VU/simplePSO',
    keywords=['optimization',
              'systems biology'],
    classifiers=['License :: OSI Approved :: BSD License',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 2',
                 'Programming Language :: Python :: 3'],
    include_package_data=True,
    install_requires=['matplotlib >= 1.5.0',
                      'numpy >= 1.11.0',
                      'scipy >= 0.17.1',
                      'pysb >= 1.1.1'],
    )
