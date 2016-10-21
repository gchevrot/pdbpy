from setuptools import setup, find_packages

with open('README.md') as file:
    long_description = file.read()

version = {}
exec(open('pdbpy/version.py').read(), version)

setup(name='pdbpy',
            version=version['version'],
            description='Structural properties of a pdb file',
            long_description=long_description,
            author='Guillaume Chevrot',
            author_email='guillaume.chevrot@univ-orleans.from',
            url='http://github.com/gchevrot/pdbpy',
            license='new BSD',
            packages=find_packages(),
            platforms='any',
            install_requires=["numpy", "urllib3"],
            provides=["ActivePapers.Py"],
            keywords=['PDB', 'computational biology', 'bioinformatics', 'PDB chemical / physical properties', 'structural biophysics'],
            classifiers=[
                          "Development Status :: 3 - Alpha",
                          "Intended Audience :: Science/Research",
                          "License :: OSI Approved :: BSD License",
                          "Operating System :: OS Independent",
                          "Programming Language :: Python :: 3.3",
                          "Programming Language :: Python :: 3.4",
                          "Programming Language :: Python :: 3.5",
                          "Topic :: Scientific/Engineering",
                        ]
     )

