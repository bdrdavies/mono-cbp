"""Setup script for mono-cbp."""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Read requirements
with open(os.path.join(this_directory, 'requirements.txt'), encoding='utf-8') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

# Read version
version = {}
with open(os.path.join(this_directory, 'mono_cbp', '__init__.py'), encoding='utf-8') as f:
    exec(f.read(), version)

setup(
    name='mono-cbp',
    version=version['__version__'],
    author='Benjamin Davies',
    author_email='bdrdavies@gmail.com',
    description='Pipeline for detecting circumbinary planets in TESS light curves',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/bdrdavies/mono-cbp',
    packages=find_packages(exclude=['tests', 'docs', 'examples']),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    python_requires='>=3.8',
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'mono-cbp=mono_cbp.cli:main',
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
