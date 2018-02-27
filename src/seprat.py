"""Functions for converting arrays of frequencies into separation
ratios.  The kinds of ratios are `r010`, `r02` and `r13`.

"""

import numpy as np


def get_ratios_from_arrays(f, n, l, ks, points=3):
    """Get ratios from arrays `f`, `n` and `l`, which contain the
    frequencies, radial orders and angular degrees of the modes.  The
    type of ratios to be computed is selected by the list of strings
    `ks`, which can be `r010`, `r02` or `r13`.  The optional argument
    `points` can be set to 3 or 5 (default 3) to select three- or
    five-point differences."""
    d = get_freq_dict(n, l, f)
    return get_ratios_from_dict(d, ks, points=points)


def get_ratios_from_dict(d, ks, points=3):
    """Get ratios from a dictionary of frequencies `d`, with entries of
    the form `d[n,l] = freq`, where `freq`, `n` and `l` are the
    frequency, radial order and angular degree of a mode.  The type of
    ratios to be computed is selected by the list of strings `ks`,
    which can be `r010`, `r02` or `r13`.  The optional argument
    `points` can be set to 3 or 5 (default 3) to select three- or
    five-point differences."""
    freqs = np.array([])
    ratios = np.array([])
    for k in ks:
        if k == 'r010':
            output = get_r010_from_dict(d, points=points)
            freqs = np.hstack((freqs, output[0]))
            ratios = np.hstack((ratios, output[1]))
        elif k == 'r02':
            output = get_r02_from_dict(d)
            freqs = np.hstack((freqs, output[0]))
            ratios = np.hstack((ratios, output[1]))
        elif k == 'r13':
            output = get_r13_from_dict(d)
            freqs = np.hstack((freqs, output[0]))
            ratios = np.hstack((ratios, output[1]))
        else:
            raise ValueError('ks must be r010, r02 or r13')
        
    return freqs, ratios


def get_r010_from_arrays(f, n, l, points=3):
    """Get r010 ratios from arrays `f`, `n` and `l`, which contain the
    frequencies, radial orders and angular degrees of the modes.  The
    optional argument `points` can be set to 3 or 5 (default 3) to
    select three- or five-point differences.    """
    d = get_freq_dict(n, l, f)
    return get_r010_from_dict(d, points=points)


def get_r02_from_arrays(f, n, l):
    """Get r02 ratios from arrays `f`, `n` and `l`, which contain the
    frequencies, radial orders and angular degrees of the modes."""
    d = get_freq_dict(n, l, f)
    return get_r02_from_dict(d)


def get_r13_from_arrays(f, n, l):
    """Get r02 ratios from arrays `f`, `n` and `l`, which contain the
    frequencies, radial orders and angular degrees of the modes."""
    d = get_freq_dict(n, l, f)
    return get_r13_from_dict(d)


def get_r010_from_dict(d, points=3):
    """Get r010 ratios from a dictionary of frequencies `d`, with entries
    of the form `d[n,l] = freq`, where `freq`, `n` and `l` are the
    frequency, radial order and angular degree of a mode.  The
    optional argument `points` can be set to 3 or 5 (default 3) to
    select three- or five-point differences."""
    smallest_l0 = np.min([n for (n,l) in d.keys() if l==0])
    smallest_l1 = np.min([n for (n,l) in d.keys() if l==1])
    if smallest_l0 <= smallest_l1:
        n = 1*smallest_l0
        l = 0
    else:
        n = 1*smallest_l1
        l = 1

    seqs = [[]]
    while n < np.max([k[0] for k in d.keys()])+1:
        try:
            seqs[-1].append(d[(n, l)])
        except KeyError:
            seqs.append([])
            
        if l == 0:
            l = 1
        else:
            l = 0
            n += 1

    if points == 3:
        freqs = np.hstack([seq[1:-1] for seq in seqs])
    elif points == 5:
        freqs = np.hstack([seq[2:-2] for seq in seqs])
    else:
        raise ValueError('points must be 3 or 5')

    ratios = np.hstack([get_r010(np.array(seq), points=points) for seq in seqs])
    return freqs, ratios


def get_r02_from_dict(d):
    """Get r02 ratios from a dictionary of frequencies `d`, with entries
    of the form `d[n,l] = freq`, where `freq`, `n` and `l` are the
    frequency, radial order and angular degree of a mode."""
    l1_keys = [k for k in d.keys() if k[1]==1]
    n = np.min([k[0] for k in l1_keys])

    seqs = [[d[(n, 1)]]]
    while n < np.max([k[0] for k in l1_keys])+1:
        try:
            if not d[(n,1)] in seqs[-1]:
                seqs[-1].append(d[n,1])
        except KeyError:
            n += 1
            continue

        try:
            seqs[-1].extend([d[n,2], d[n+1,0], d[n+1,1]])
        except KeyError:
            seqs.append([])

        n += 1

    seqs = [seq for seq in seqs if len(seq)>0]
    freqs = np.hstack([seq[2::3] for seq in seqs])
    ratios = np.hstack([get_r02(np.array(seq)) for seq in seqs])

    return freqs, ratios


def get_r13_from_dict(d):
    """Get r13 ratios from a dictionary of frequencies `d`, with entries
    of the form `d[n,l] = freq`, where `freq`, `n` and `l` are the
    frequency, radial order and angular degree of a mode."""
    l0_keys = [k for k in d.keys() if k[1]==0]
    n = np.min([k[0] for k in l0_keys])

    seqs = [[d[(n, 0)]]]
    while n < np.max([k[0] for k in l0_keys])+1:
        try:
            if not d[(n,0)] in seqs[-1]:
                seqs[-1].append(d[n,0])
        except KeyError:
            n += 1
            continue

        try:
            seqs[-1].extend([d[n-1,3], d[n,1], d[n+1,0]])
        except KeyError:
            seqs.append([])

        n += 1

    seqs = [seq for seq in seqs if len(seq)>2]
    freqs = np.hstack([seq[2::3] for seq in seqs])
    # cheat, because getting r13 is roughly the same as getting r02
    ratios = np.hstack([get_r02(np.array(seq)) for seq in seqs])

    return freqs, ratios


def get_r010(freqs, points=3):
    """Computes r010 ratios from a continuous sequence of radial and
    dipole mode frequencies in `freqs`.  The optional argument
    `points` can be set to 3 or 5 (default 3) to select three- or
    five-point differences."""

    if points == 3:
        return np.abs(freqs[1:-1]-(freqs[2:]+freqs[:-2])/2.)/(freqs[2:]-freqs[:-2])
    elif points == 5:
        return np.abs(freqs[2:-2]*3./4.-(freqs[1:-3]+freqs[3:-1])/2.+(freqs[:-4]+freqs[4:])/8.) \
            /(freqs[3:-1]-freqs[1:-3])
    else:
        raise ValueError('points must be equal to 3 or 5')



def get_r02(freqs):
    """Computes r02 (r13) ratios from a continuous sequence of mode
    frequencies `freqs` with angular degrees 1, 2, 0, 1, ... (0, 3, 1,
    0, ...)."""
    return (freqs[2::3]-freqs[1::3])/(freqs[3::3]-freqs[:-3:3])


def get_freq_dict(n, l, f):
    return {(int(ni), int(li)): fi for (ni, li, fi) in zip(n, l, f)}
