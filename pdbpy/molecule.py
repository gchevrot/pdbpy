import sys
import numpy as np
from pdbpy.extract import extract_coordinates, extract_calpha
from pdbpy.residues import extract_residues
from pdbpy.data import aa_sidechain_chemical_properties as aa_hydrophobicity

class Molecule:
    def __init__(self, pdb_name, download_from_pdb=True):
        """
        Parameters
        ----------
        pdb_name:
            Name of the pdb file. (ex: 1dpx or 1dpx.pdb) 

        download_from_pdb:
            default is True. Use the download_pdb function (need an internet connection)
            If False, use a local pdb file.

        Return
        ------
        A Molecule object containing the coordinates in nanometer
        """
        # Extracting coordinates from a pdb file - result is a file (pdb_name_coordinates.pdb)
        self.pdb_name = pdb_name
        self.download_from_pdb = download_from_pdb
        extract_coordinates(pdb_name, download_from_pdb=download_from_pdb)
        if pdb_name[-4:] == '.pdb':
            pdb_file = pdb_name[:-4] + '_coordinates.pdb'
        else:
            pdb_file = pdb_name + '_coordinates.pdb'
        # Number of columns
        with open(pdb_file, 'r') as f:
            line = f.readline()
            n_columns = len(line.split())
        # coordinates - divide by 10 for the conversion angstrom ==> nanometer
        if n_columns == 12:
            self.xyz = np.loadtxt(pdb_file, usecols=(-6,-5,-4)) / 10
            self.n_res = len(set(np.loadtxt(pdb_file, usecols=(-7,))))
        elif n_columns == 11: # case whne the before last column is missing
            self.xyz = np.loadtxt(pdb_file, usecols=(-5,-4,-3)) / 10
            self.n_res = len(set(np.loadtxt(pdb_file, usecols=(-6,))))
        else:
            print('''Number of columns in the PDB file is "anormal" and cannot 
                  be treated with the current code.''')
            sys.exit()

    def center_of_gravity(self):
        """
        Return the center of gravity from atomic coordinates
        """
        return self.xyz.mean(axis=0)

    def radius_of_gyration(self):
        """ 
        Return the radius of gyration (in nm)
        """
        dist = (self.xyz - self.center_of_gravity())**2
        # My old wrong definition:
        #dist = (dist.sum(axis=1))**0.5
        #return dist.mean()
        dist = (dist.sum(axis=1)).mean()
        return dist**0.5

    def radius_of_gyration_normalized(self):
        """
        Return the radius of gyration (in nm) normalized with the number of residues
        """
        return self.radius_of_gyration() / self.n_res

    def hydrophobicity(self):
        """
        Return the percentage of hydrophobic residue
        """
        res_sequence = extract_residues(self.pdb_name, download_from_pdb=self.download_from_pdb)
        hydrophilic = 0
        hydrophobic = 0
        for res in res_sequence:
            if aa_hydrophobicity[res] == 'hydrophilic':
                hydrophilic += 1
            else:
                hydrophobic += 1
        return hydrophobic/(hydrophilic+hydrophobic)*100

    def msd(self):



