import re
from pdbpy.download import download_pdb


def is_dna_or_rna(pdb_name, download_from_pdb=True):
    """
    Test if the the pdb file corresponds to a DNA or RNA structure 

    Parameters
    ----------
    pdb_name:
        Name of the pdb file. (ex: 1dpx or 1dpx.pdb) 

    download_from_pdb:
        default is True. Use the download_pdb function (need an internet connection)
        If False, use a local pdb file.

    Return
    ------
    result: bool
        True if the it is a DNA or a RNA molecule, False otherwise
    """
    if download_from_pdb:
        download_pdb(pdb_name)
    if pdb_name[-4:] == '.pdb':
        pdb_file = pdb_name
    else:
        pdb_file = pdb_name + '.pdb'

    result = False
    with open(pdb_file, 'r') as input:
        for line in input:
            # Save only the 1st chain
            if line[:10] == 'COMPND   2':
                if re.search(r'DNA', line) or re.search(r'RNA', line):
                    result = True

    return result

