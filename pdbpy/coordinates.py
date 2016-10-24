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
                output.write(line)

