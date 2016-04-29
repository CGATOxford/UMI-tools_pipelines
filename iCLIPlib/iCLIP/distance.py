''' This module contains a varety of diffrent functions for calculating
distance between two profiles'''

from utils import spread
import numpy as np


##################################################
def calcAverageDistance(profile1, profile2):
    ''' This function calculates the average distance of all
    pairwise distances in two profiles'''

    def _cartesian(x, y):
        return np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])

    positions = _cartesian(profile1.index.values, profile2.index.values)
    counts = _cartesian(profile1.values, profile2.values)
    counts = np.prod(counts, axis=1)
    distances = np.abs(positions[:, 0] - positions[:, 1])
    mean_distance = (distances.astype("float64") *
                     counts).sum() / np.sum(counts)

    return mean_distance


##################################################
def findMinDistance(profile1, profile2):
    '''Finds mean distance between each read in profile1
    and a read in profile2'''

    locations1 = profile1.index.values

    locations2 = profile2.index.values  # .astype("int16")

    mat1 = np.repeat(locations1, locations2.size).reshape(
        (locations1.size, locations2.size))
    mat2 = np.tile(locations2, locations1.size).reshape(
        (locations1.size, locations2.size))

    distances = np.abs(mat1-mat2).min(axis=1)

    return distances.mean()


##################################################
def corr_profile(profile1, profile2, nspread, profile2_ready=False):
    
    profile1 = profile1.reindex(
                    range(int(profile1.index.values.min())-1,
                          int(profile1.index.values.max())+1)).fillna(0)
    profile1 = spread(profile1, nspread, False)
        
    if not profile2_ready:
        profile2 = profile2.reindex(
                       range(int(profile2.index.values.min()),
                             int(profile2.index.values.max()))).fillna(0)
        profile2 = spread(profile2, nspread, False)
        
    return profile1.corr(profile2, method="spearman")
