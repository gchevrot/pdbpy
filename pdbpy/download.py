import urllib3


def download_pdb(pdb_name):
    """
    Download a pdb file from the PDB. Need an internet connection.

    Parameters
    ----------
    Name of the pdb file. (ex: 1dpx or 1dpx.pdb) 

    Return
    ------
    The pdb file
    """
    # URL of the PDB
    url = 'http://files.rcsb.org/download/'
    if pdb_name[-4:] == '.pdb':
        file_name = pdb_name
    else:
        file_name = pdb_name + '.pdb'
    http = urllib3.PoolManager()
    url = url + file_name
    req = http.request('GET', url, preload_content = False)
    with open(file_name, 'wb') as out:
        while True:
            data = req.read()
            if not data:
                break
            out.write(data)
    req.release_conn()
