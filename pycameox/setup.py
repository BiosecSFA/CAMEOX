"""pycameox setup.py"""
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pycameox',
    version='0.0.3',
    author='Jose Manuel Martí',
    author_email='martimartine1@llnl.gov',
    description='Python library for CAMEOX (CAMEOs eXtended)',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/BiosecSFA/CAMEOX',
    project_urls={
        'Bug Tracker': 'https://github.com/BiosecSFA/CAMEOX/issues',
        },
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Framework :: IPython',
        'Framework :: Jupyter',
        'Intended Audience :: Science/Research',
        'License :: Other/Proprietary License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Julia',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.8",
    install_requires=[
        'matplotlib',
        'numpy',
        'pandas',
        'plotly>4.14.0,<5',
        'scipy',
        ],
)
