from download import download_pdb


def extract_residues(pdb_name, download_from_pdb=True):
    """
    Extracting the residue sequence 

    Parameters
    ----------
    pdb_name:
        Name of the pdb file. (ex: 1dpx or 1dpx.pdb) 

    download_from_pdb:
        default is True. Use the download_pdb function (need an internet connection)
        If False, use a local pdb file.

    Return
    ------
    res_seq: set
        The sequence of residue
    """
    if download_from_pdb:
        download_pdb(pdb_name)
    if pdb_name[-4:] == '.pdb':
        pdb_file = pdb_name
    else:
        pdb_file = pdb_name + '.pdb'

    res = []
    res_number = []

    with open(pdb_file, 'r') as input:
        for line in input:
            # Save only the 1st chain
            if line[:3] == 'TER':
                break
            if line[:4] == 'ATOM':
                res.append(line[17:20])     # columns 17 to 19 correspond to the residue name
                if len(line.split()) == 11:
                    res_number.append(int(line.split()[4]))
                if len(line.split()) == 12:
                    res_number.append(int(line.split()[5]))

    # extract the residue sequence
    res_seq = []
    temp = -1000
    for i, residue_number in enumerate(res_number):
        if residue_number > temp:
            temp = residue_number
            res_seq.append(res[i]) 
    return res_seq
