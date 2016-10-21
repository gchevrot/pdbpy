From a PDB file (local file or through the PDB on internet), calculate the properties of the 1st chain.

Available data and properties
-----------------------------

Number of residues, coordinates of the atoms, center of gravity, radius of
gyration.

For an example, see this [notebook](https://github.com/gchevrot/pdbpy/blob/master/examples/example.ipynb) 


Requirements
------------

pdppy has been tested only with Python 3.5. It requires the following
libraries:

 - urllib3 (tested with version 1.18)
 - numpy 1.6 or later [http://numpy.scipy.org/](http://numpy.scipy.org/)


Installation
------------

 git clone https://github.com/gchevrot/pdbpy.git
 
 python setup.py install
