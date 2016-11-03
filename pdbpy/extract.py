import numpy as np
from pdbpy.download import download_pdb


def extract_coordinates(pdb_name, download_from_pdb=True):
    """
    Extracting the lines containing the coordinates of the 1st chain

    Parameters
    ----------
    pdb_name:
        Name of the pdb file. (ex: 1dpx or 1dpx.pdb) 

    download_from_pdb:
        default is True. Use the download_pdb function (need an internet connection)
        If False, use a local pdb file.

    Return
    ------
    coordinates: numpy array, dimension: (n, 3)
                coordinates in nanometers 
    """
    if download_from_pdb:
        download_pdb(pdb_name)
    if pdb_name[-4:] == '.pdb':
        pdb_file = pdb_name
    else:
        pdb_file = pdb_name + '.pdb'

    coordinates = []
    with open(pdb_file, 'r') as input:
        for line in input:
            # Save only the 1st chain
            if line[:3] == 'TER':
                break
            if line[:4] == 'ATOM':
                # Sometimes in X-ray cristallography, one sees superposition of 
                # differents positions for each atom. We will extract only the first position
                # The different positions are denoted with a letter preceding the residue name 
                # in column 17. Note the column 55-60 (field "occupancy") correspond to a reduced 
                # electronic density for each atom
                if line[16] == ' ' or line[16] == 'A': 
                    coordinates.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    # Divide by 10, so coordinates are in nanometers
    return np.array(coordinates) / 10


def extract_calpha_coordinates(pdb_name, download_from_pdb=True):
    """
    Extracting the lines containing the carbon alpha coordinates of the 1st chain

    Parameters
    ----------
    pdb_name:
        Name of the pdb file. (ex: 1dpx or 1dpx.pdb) 

    download_from_pdb:
        default is True. Use the download_pdb function (need an internet connection)
        If False, use a local pdb file.

    Return
    ------
    coordinates: numpy array, dimension: (n, 3)
                coordinates in nanometers 
    """
    if download_from_pdb:
        download_pdb(pdb_name)
    if pdb_name[-4:] == '.pdb':
        pdb_file = pdb_name
    else:
        pdb_file = pdb_name + '.pdb'

    coordinates = []
    with open(pdb_file, 'r') as input:
        for line in input:
            # Save only the 1st chain
            if line[:3] == 'TER':
                break
            if line[:4] == 'ATOM':
                if line[13:15] == 'CA':
                    # Sometimes in X-ray cristallography, one sees superposition of 
                    # differents positions for each atom. We will extract only the first position
                    # The different positions are denoted with a letter preceding the residue name 
                    # in column 17. Note the column 55-60 (field "occupancy") correspond to a reduced 
                    # electronic density for each atom
                    if line[16] == ' ' or line[16] == 'A': 
                        coordinates.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    # Divide by 10, so coordinates are in nanometers
    return np.array(coordinates) / 10


def extract_coordinates_to_file(pdb_name, download_from_pdb=True):
    """
    Extracting the lines containing the coordinates of the 1st chain

    Parameters
    ----------
    pdb_name:
        Name of the pdb file. (ex: 1dpx or 1dpx.pdb) 

    download_from_pdb:
        default is True. Use the download_pdb function (need an internet connection)
        If False, use a local pdb file.

    Return
    ------
    A text file (note that the coordinates are in Angstrom)
    """
    if download_from_pdb:
        download_pdb(pdb_name)
    if pdb_name[-4:] == '.pdb':
        pdb_file = pdb_name
        output_name = pdb_name[:-4]
    else:
        pdb_file = pdb_name + '.pdb'
        output_name = pdb_name

    with open(pdb_file, 'r') as input, open (output_name+'_coordinates.pdb', 'w') as output:
        for line in input:
            # Save only the 1st chain
            if line[:3] == 'TER':
                break
            if line[:4] == 'ATOM':
                # Sometimes in X-ray cristallography, one sees superposition of 
                # differents positions for each atom. We will extract only the first position
                # The different positions are denoted with a letter preceding the residue name 
                # in column 17. Note the column 55-60 (field "occupancy") correspond to a reduced 
                # electronic density for each atom
                if line[16] == ' ' or line[16] == 'A': 
                    output.write(line)

def extract_calpha_coordinates_to_file(pdb_name, download_from_pdb=True):
    """
    Extracting the lines containing the carbon alpha coordinates of the 1st chain

    Parameters
    ----------
    pdb_name:
        Name of the pdb file. (ex: 1dpx or 1dpx.pdb) 

    download_from_pdb:
        default is True. Use the download_pdb function (need an internet connection)
        If False, use a local pdb file.

    Return
    ------
    A text file (note that the coordinates are in Angstrom)
    """
    if download_from_pdb:
        download_pdb(pdb_name)
    if pdb_name[-4:] == '.pdb':
        pdb_file = pdb_name
        output_name = pdb_name[:-4]
    else:
        pdb_file = pdb_name + '.pdb'
        output_name = pdb_name

    with open(pdb_file, 'r') as input, open (output_name+'_calpha.pdb', 'w') as output:
        for line in input:
            # Save only the 1st chain
            if line[:3] == 'TER':
                break
            if line[:4] == 'ATOM':
                if line[13:15] == 'CA':
                    # Sometimes in X-ray cristallography, one sees superposition of 
                    # differents positions for each atom. We will extract only the first position
                    # The different positions are denoted with a letter preceding the residue name 
                    # in column 17. Note the column 55-60 (field "occupancy") correspond to a reduced 
                    # electronic density for each atom
                    if line[16] == ' ' or line[16] == 'A': 
                        output.write(line)

