import sys
import numpy as np
from pdbpy.extract import extract_coordinates, extract_calpha_coordinates
from pdbpy.residues import extract_residues
from pdbpy.data import aa_sidechain_chemical_properties as aa_hydrophobicity
from pdbpy.msd import msd, msd_fft

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
        """
        self.pdb_name = pdb_name
        self.download_from_pdb = download_from_pdb
        # Extracting coordinates from a pdb file - result is a file (pdb_name_coordinates.pdb)
        extract_coordinates(pdb_name, download_from_pdb=download_from_pdb)
        # Extracting coordinates from a pdb file - result is a file (pdb_name_calpha.pdb)
        extract_calpha_coordinates(self.pdb_name, download_from_pdb=self.download_from_pdb)

    def coordinates(self):
        """
        Return
        ------
        The coordinates in nanometer
        """
        if self.pdb_name[-4:] == '.pdb':
            pdb_file = self.pdb_name[:-4] + '_coordinates.pdb'
        else:
            pdb_file = self.pdb_name + '_coordinates.pdb'
        # Number of columns
        with open(pdb_file, 'r') as f:
            line = f.readline()
            n_columns = len(line.split())
        # coordinates - divide by 10 for the conversion angstrom ==> nanometer
        if n_columns == 12:
            coordinates = np.loadtxt(pdb_file, usecols=(-6,-5,-4)) / 10
        elif n_columns == 11: # case whne the before last column is missing
            coordinates = np.loadtxt(pdb_file, usecols=(-5,-4,-3)) / 10
        else:
            print('''{}: Number of columns in the PDB file is "anormal" and cannot 
                  be treated with the current code.'''.format(self.pdb_name))
            sys.exit()
        return coordinates

    def calpha_coordinates(self):
        """
        Return
        ------
        The coordinates in nanometer of the carbon alpha
        """
        if self.pdb_name[-4:] == '.pdb':
            pdb_file = self.pdb_name[:-4] + '_calpha.pdb'
        else:
            pdb_file = self.pdb_name + '_calpha.pdb'
        # Number of columns
        with open(pdb_file, 'r') as f:
            line = f.readline()
            n_columns = len(line.split())
        # coordinates - divide by 10 for the conversion angstrom ==> nanometer
        if n_columns == 12:
            coordinates = np.loadtxt(pdb_file, usecols=(-6,-5,-4)) / 10
        elif n_columns == 11: # case whne the before last column is missing
            coordinates = np.loadtxt(pdb_file, usecols=(-5,-4,-3)) / 10
        else:
            print('''{}: Number of columns in the PDB file is "anormal" and cannot 
                  be treated with the current code.'''.format(self.pdb_name))
            sys.exit()
        return coordinates

    def number_of_residues(self):
        """
        Return
        ------
        The number of residues
        """
        # Extracting coordinates from a pdb file - result is a file (pdb_name_coordinates.pdb)
        if self.pdb_name[-4:] == '.pdb':
            pdb_file = self.pdb_name[:-4] + '_coordinates.pdb'
        else:
            pdb_file = self.pdb_name + '_coordinates.pdb'
        # Number of columns
        with open(pdb_file, 'r') as f:
            line = f.readline()
            n_columns = len(line.split())
        # coordinates - divide by 10 for the conversion angstrom ==> nanometer
        if n_columns == 12:
            n_res = len(set(np.loadtxt(pdb_file, usecols=(-7,))))
        elif n_columns == 11: # case whne the before last column is missing
            n_res = len(set(np.loadtxt(pdb_file, usecols=(-6,))))
        else:
            print('''{}: Number of columns in the PDB file is "anormal" and cannot 
                  be treated with the current code.'''.format(self.pdb_name))
            sys.exit()
        return n_res 

    def center_of_gravity(self):
        """
        Return the center of gravity from atomic coordinates
        """
        return self.coordinates().mean(axis=0)

    def radius_of_gyration(self):
        """ 
        Return the radius of gyration (in nm)
        """
        dist = (self.coordinates() - self.center_of_gravity())**2
        # My old wrong definition:
        #dist = (dist.sum(axis=1))**0.5
        #return dist.mean()
        dist = (dist.sum(axis=1)).mean()
        return dist**0.5

    def radius_of_gyration_normalized(self):
        """
        Return the radius of gyration (in nm) normalized with the number of residues
        """
        return self.radius_of_gyration() / self.number_of_residues()

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

    def msl(self):
        """
        Return the mean square length of the protein calculted with the C-alpha atoms. 
        It is calculated like a MSD.
        """
        #return msd(self.coordinates())
        return msd(self.calpha_coordinates())

    def msl_fft(self):
        """
        Return the mean square length of the protein calculted with the C-alpha atoms. 
        It is calculated like a MSD (using FFT in this case to calculate the MSD).
        """
        #return msd_fft(self.coordinates())
        return msd_fft(self.calpha_coordinates())

