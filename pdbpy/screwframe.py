## This module implements the module method
## The code is adaptated from the ActivePaper ScrewFrame (https://zenodo.org/record/21690#.WBtGQZPhDfA)

import numpy as np
from pdbpy.geometry import norm, normalize, screw_motion


def frenet_basis(point, before, after):
    """
    Computes a discrete approximation to the Frenet basis at the
    position of a C-alpha atom in a peptide chain.

    Parameters
    ----------
    point: numpy.ndarray
           Position of the C-alpha atom
    before: numpy.ndarray
            Position of the preceding C-alpha atom
    after: numpy.ndarray
            Position of the following C-alpha atom

    Returns
    -------
    numpy.ndarray
    Tangent and normal vector
    """
    # Tangent vector t
    # See definition in Acta Cryst. (2015), D71, 1411: 
    # t = r'/norm(r')  with r'(i) = [r(i+1) - r(i-1)] / 2\delta\lambda    (r: position, \delta\lambda: dist between C-alpha
    t = normalize(after-before)
    # Normal vector - code from konrad's ActivePaper
    da = after-point
    n = normal(da - np.dot(da, t)*t)
    return np.array([t, n])

def frenet_bases(points):
    """
    tangent and normal vectors. Same algorithm as frenet_basis, the difference is that
    tangent and normal vectors from the all the positions of the C-alpha
    atoms are computed.
    It is equivalent(?) to the frenet_chain in Konrad's ActivePaper.

    Parameters
    ----------
    points: numpy.ndarray
            Positions of the C-alpha atoms

    Returns
    -------
    numpy.ndarray
    tangent and normal vectors
    dimensions: (a, b, c) with 
        a = 0 for the tangent vectors
        a = 1 for the normal vectors
        b: number of C-alpha atoms / frenet basis
        c = 0: x-coordinate
        c = 1: y-coordinate
    It has two elements less than the list of C-alpha positions
    because no useful Frenet bases can be defined for the
    first and last position.
    """
    # tangent vectors ("point_after - point_before")
    t = normalize(points[2:] - points[:-2])
    # normal vectors 
    # "point_after - point_current"
    da = points[2:] - points[1:-1]
    n = [v_da - np.dot(v_da, v_t) * v_t for v_da, v_t in zip(da, t)]
    n = normalize(np.array(n))
    return np.array([t, n])

def frame_rotation_and_distance(fbs):
    """
    Computes the orientation change from one Frenet basis
    to another, and the angular distance between 2 bases.

    Parameters
    ----------
    fbs: np.ndarray
      shape: see frenet_bases
      all the frenet bases obtained from the C-alpha atoms

    Returns
    -------
    list of tuples
        (q, delta): 
            q: quaternion describing the rotation that transforms
               a frenet basis to the next
            delta: angular distance
    """
    results = []   # will contain the tuples q and delta
    tangent = 0    # fbs[0] <==> tangent vectors 
    normal = 1     # fbs[1] <==> normal vectors
    # tangent vectors: sum and difference between the successive frenet bases
    sum_tangents = fbs[tangent, 0:-1] + fbs[tangent, 1:]
    diff_tangents = fbs[tangent, 0:-1] - fbs[tangent, 1:]
    # normal vectors: sum and difference between the successive frenet bases
    sum_normals = fbs[normal, 0:-1] + fbs[normal, 1:]
    diff_normals = fbs[normal, 0:-1] - fbs[normal, 1:]

    for i in range(len(sum_tangents)):
        m = np.zeros((4, 4), np.float)

        s = sum_tangents[i]
        d = diff_tangents[i]
        k = np.array([[0,      d[0],  d[1],  d[2]],
                      [-d[0],  0,     s[2], -s[1]],
                      [-d[1], -s[2],  0,     s[0]],
                      [-d[2],  s[1], -s[0],  0]])
        m += np.dot(k.T, k)

        s = sum_normals[i]
        d = diff_normals[i]
        k = np.array([[0,      d[0],  d[1],  d[2]],
                      [-d[0],  0,     s[2], -s[1]],
                      [-d[1], -s[2],  0,     s[0]],
                      [-d[2],  s[1], -s[0],  0]])
        m += np.dot(k.T, k)


        l, vs = np.linalg.eigh(m)
        q = vs[:, np.argsort(l)[0]] # the eigenvector for the smallest eigenvalue
        if q[0] < 0:                # ensure that q[0] is non-negatie
            q = -q
        delta = np.sqrt(m[0,0]/8.)  # angular distance
        results.append((q, delta))
    return results

def screwframe_rotation_centers(calphas):
    """
    Return the successive rotation centers between 2 frenet basis

    Parameters
    ----------
    calphas: np.ndarray
    Coordinates of the successive C-alpha atoms

    Return
    ------
    np.ndarray
    Successive coordinates of the rotation center 
    (length: number of frenet basis - 1) 
    """
    fbs = frenet_bases(calphas)
    # frenet basis for first and last C-alpha are not computed, so these 2 atoms are discarded:
    calphas = calphas[1:-1]
    assert len(fbs[0]) == len(calphas)  # frenet basis for first and last C-alpha are not computed
    # quaternions that describe the rotation between the frenet bases
    q_delta = frame_rotation_and_distance(fbs)
    q = [q_delta[i][0] for i in range(len(q_delta))]
    q = np.array(q)
    t = calphas[1:] - calphas[:-1]
    assert len(q) == len(t) 
    # screw parameters from quaternions and translation
    r0_axis_angle_d = screw_motion(q, t)
    r0 = [i[0] for i in r0_axis_angle_d]
    r0 = np.array(r0)
    # r_screw is the point on the screw axis that is closest to the C atom
    r_screw = calphas[:-1] + r0
    return r_screw


