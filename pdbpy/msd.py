import numpy as np

# MSD straightforward implementation
def msd(r):
    """
    Mean square displacement

    Parameters
    ----------
    r: numpy array. dimensions: (t, 3)
       Array containing the coordinates along t
    """
    shifts = np.arange(len(r))
    msds = np.zeros(shifts.size)    
    for i, shift in enumerate(shifts):
        diffs = r[:-shift if shift else None] - r[shift:]
        sqdist = np.square(diffs).sum(axis=1)
        msds[i] = sqdist.mean()
    return msds

# MSD FFT
# Algorithm comes from this paper - DOI:  http://dx.doi.org/10.1051/sfn/201112010 
# An implemenation has been proposed here: http://stackoverflow.com/questions/34222272/computing-mean-square-displacement-using-python-and-fft
def autocorrfft(x):
    N=len(x)
    F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
    #  power spectral density
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    # Autocorrelation in convention B
    res = (res[:N]).real
    # divide res(m) by (N-m) to obtain autocorrelation in convention A
    n = N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
    return res/n

def msd_fft(r):
    """
    Mean square displacement (using FFT)

    Parameters
    ----------
    r: numpy array. dimensions: (t, 3)
       Array containing the coordinates along t
    """
    N = len(r)
    D = np.square(r).sum(axis=1)
    D = np.append(D,0)
    S2 = sum([autocorrfft(r[:, i]) for i in range(r.shape[1])])
    Q = 2*D.sum()
    S1 = np.zeros(N)
    for m in range(N):
        Q = Q-D[m-1]-D[N-m]
        S1[m] = Q/(N-m)
    return S1-2*S2
