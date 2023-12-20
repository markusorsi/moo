from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Molecule Overlap Optimizer (MOO)'
LONG_DESCRIPTION = 'RDKit wrapper to generate 3D structures with maxmium overlap.'

# Setting up
setup(
        name="moo-chem", 
        version=VERSION,
        author="Markus Orsi",
        author_email="<markus.orsi@unibe.ch>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], 
        keywords=['cheminformatics', 'moo'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)