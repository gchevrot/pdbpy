import numpy as np

def norm(v):
    """
    Norm of a vector

    Parameter
    ---------
    v: numpy array
       vector

    Return
    ------
    numpy float
    norm of the vector v
    """
    return np.sqrt(v.dot(v))

def normalize(vectors):
    """
    Return the unit vector of the vectors (could be one or several vectors)

    Parameter
    ---------
    vector: numpy array
       vector

    Return
    ------
    numpy array
    The unit vector of vectors
    """
    if len(vectors.shape) == 1:
        unit_vectors = vectors/norm(vectors)
    elif len(vectors.shape) == 2:
        unit_vectors = [v/norm(v) for v in vectors]
        unit_vectors = np.array(unit_vectors)
    else:
        return NotImplemented
    return unit_vectors 

def screw_motion(quaternions, translations):
    """
    Compute screw parameters from quaternions and translation
    (adapted from the original code in screwframe.ap)

    Parameters
    ----------
    quaternions: np.ndarray
       all normalized quaternion with q[0]>0, describing the successive rotations between frenet bases
    translations: np.ndarray
       all translation vectors between C-alpha atoms

    Returns
    -------
    list of tuple (r0, axis, phi, d): 
        r0: point on the screw axis
        axis: normalized vector indicating the axis direction
        phi: angle of the rotation around the axis
        d: scalar displacement along the axis
    """
    results = []
    for q, t in zip(quaternions, translations):
        cosphi2 = q[0]
        if abs(cosphi2-1.) < 1.e-6:
            # No rotation, pure translation.
            # We arbitrarily choose the axis parallel to the translation.
            phi = 0.
            r0 = np.array([0., 0., 0.])
            d = norm(t)
            axis = t/d
        else:
            # From the definition of np.arccos, and with
            # cosphi2 > 0, we get 0 <= phi <= pi.
            phi = 2.*np.arccos(cosphi2)
            # We choose sinphi2 as the positive solution,
            # which is consistent with the range of phi.
            sinphi2 = np.sqrt(1.-cosphi2*cosphi2)
            axis = q[1:]/sinphi2
            # We choose the direction of the axis such that the
            # scalar displacement parameter is positive.
            d = np.dot(t, axis)
            if d < 0:
                d = -d
                axis = -axis
                phi = -phi
                sinphi2 = -sinphi2
            # r0 is the axis point closest to the reference point for the
            # rotation, relative to the reference point of the rotation.
            x = t-d*axis
            r0 = 0.5*(x - (cosphi2/sinphi2)*np.cross(axis, x))
        results.append((r0, axis, phi, d))
    return results

