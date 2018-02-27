#!/usr/bin/env python

import numpy as np
from seprat import get_ratios_from_arrays

def get_ells(x):
    idiff = np.where(np.diff(x)<0.)[0].reshape((-1,1))
    i = np.arange(len(x)).reshape((1,-1))
    return np.sum(i > idiff, axis=0)

obs = np.array([2496.29, 2629.84, 2764.32, 2899.11, 3033.92, 3168.94,
                3303.73, 3439.17, 3575.23, 3711.6, 3848.2, 3984.79,
                2559.31, 2693.56, 2828.36, 2963.45, 3098.39, 3233.4,
                3368.88, 3504.5, 3640.95, 3777.21, 3913.93, 4051.51,
                2486.11, 2619.87, 2754.72, 2889.66, 3024.81, 3160.05,
                3295.05, 3430.88, 3566.95, 3703.7, 3840.9, 3977.52])

err = np.array([0.07, 0.06, 0.05, 0.04, 0.04, 0.05, 0.06, 0.1, 0.15,
                0.26, 0.5, 0.65, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05,
                0.08, 0.09, 0.16, 0.2, 0.36, 0.61, 0.13, 0.08, 0.07,
                0.06, 0.07, 0.07, 0.09, 0.14, 0.24, 0.32, 0.64, 0.87])
               
l = get_ells(obs)

# the radial orders don't have to be right:
# they just have to be relatively right
n = np.zeros_like(l, dtype=int)
Dnu = np.median(np.diff(obs[l==0]))  # might need improvement
n[0] = np.floor(obs[0]/Dnu)
n[l==0] = n[0] + np.around((obs[l==0]-obs[0])/Dnu)
n[l==1] = n[0] + np.around((obs[l==1] - 0.5*Dnu - obs[0])/Dnu)
n[l==2] = n[0] + np.around((obs[l==2] - obs[0])/Dnu) - 1
n[l==3] = n[0] + np.around((obs[l==3] - 1.5*Dnu - obs[0])/Dnu)

x, r = get_ratios_from_arrays(obs, n, l, ['r010', 'r02'], points=5)
randerr = np.random.randn(100000, len(obs))*err
R = np.vstack([get_ratios_from_arrays(obs + row, n, l,
                                         ['r010', 'r02'], points=5)[1]
               for row in randerr]).T

# print(np.cov(R))
C = np.cov(R)

try:
    tmp = np.linalg.cholesky(C)
    print('Cholesky test: covariance matrix is positive definite')
except LinAlgError:
    print('Cholesky test: covariance matrix is NOT positive definite')

if np.all(np.linalg.eigvals(C) > 0):
    print('Eigval test:   covariance matrix is positive definite')
else:
    print('Eigval test:   covariance matrix is NOT positive definite')

np.savetxt('ratio_means.txt', np.vstack((x, np.mean(R, axis=1))).T)
np.savetxt('ratio_cov.txt', C)
np.savetxt('ratio_invcov.txt', np.linalg.inv(C))
