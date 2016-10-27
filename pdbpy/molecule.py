import sys
import numpy as np
from pdbpy.extract import extract_coordinates, extract_calpha_coordinates
from pdbpy.residues import extract_residues
from pdbpy.data import aa_sidechain_chemical_properties as aa_hydrophobicity
from pdbpy.msd import msd, msd_fft
from pdbpy.download import download_pdb
from pdbpy.inspection import is_dna_or_rna


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
        # Verifying that it is not a RNA or DNA molecule
        if is_dna_or_rna(self.pdb_name, self.download_from_pdb):
            #print("{} corresponds to a DNA or RNA molecule. This code cannot analyze DNA or RNA.".format(self.pdb_name))
            sys.exit()
        #if self.download_from_pdb:
        #    download_pdb(self.pdb_name)
        # Extract coordinates within __init__ rather than with dedicated function for performance reasons
        self.coordinates = extract_coordinates(self.pdb_name, download_from_pdb=self.download_from_pdb) 
        self.calpha_coordinates = extract_calpha_coordinates(self.pdb_name, download_from_pdb=self.download_from_pdb) 
    
#    def coordinates(self):
#        """
#        Return
#        ------
#        coordinates: numpy array, dimension: (3, n)
#            The coordinates in nanometer
#        """
#        coordinates = extract_coordinates(self.pdb_name, download_from_pdb=self.download_from_pdb) 
#        return coordinates
#
#    def calpha_coordinates(self):
#        """
#        Return
#        ------
#        coordinates: numpy array, dimension: (3, n)
#            The coordinates in nanometer of the carbon alpha
#        """
#        calpha_coordinates = extract_calpha_coordinates(self.pdb_name, download_from_pdb=self.download_from_pdb) 
#        return calpha_coordinates

    def number_of_residues(self):
        """
        Return
        ------
        The number of residues
        """
        # Extracting coordinates from a pdb file - result is a file (pdb_name_coordinates.pdb)
        if self.pdb_name[-4:] == '.pdb':
            pdb_file = self.pdb_name
        else:
            pdb_file = self.pdb_name + '.pdb'
        
        res_number = []
        with open(pdb_file, 'r') as f:
            for line in f:
                if line[:3] == 'TER':
                    break
                if line[:4] == 'ATOM':
                    res_number.append(int(line[22:26]))
                n_res = len(set(res_number))
        return n_res

    def center_of_gravity(self):
        """
        Return the center of gravity from atomic coordinates
        """
        return self.coordinates.mean(axis=0)

    def radius_of_gyration(self):
        """ 
        Return the radius of gyration (in nm)
        """
        dist = (self.coordinates - self.center_of_gravity())**2
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
            # sometimes residue name is not known - so hydrophobicity cannot be calculated
            if res == 'UNK':
                return 'NaN'
            if aa_hydrophobicity[res] == 'hydrophilic':
                hydrophilic += 1
            else:
                hydrophobic += 1
        return hydrophobic/(hydrophilic+hydrophobic)*100

    def msl(self, calpha=True, all_atoms=False):
        """
        Return the mean square length of the protein calculted with the C-alpha 
        atoms if calpha is True or for all atoms if all_atoms is True
        It is calculated like a MSD.
        """
        # if both are True, there is a problem, you should choose only one 
        if calpha and all_atoms:
            print("Choose calpha or all_atoms! You cannot choose both!")
        if calpha:
            return msd(self.calpha_coordinates)
        if all_atoms:
            return msd(self.coordinates)

    def msl_fft(self, calpha=True, all_atoms=False):
        """
        Return the mean square length of the protein calculted with the C-alpha 
        atoms if calpha is True or for all atoms if all_atoms is True
        It is calculated like a MSD (using FFT in this case to calculate the MSD).
        """
        # if both are True, there is a problem, you should choose only one 
        if calpha and all_atoms:
            print("Choose calpha or all_atoms! You cannot choose both!")
        if calpha:
            return msd_fft(self.calpha_coordinates)
        if all_atoms:
            return msd_fft(self.coordinates)

