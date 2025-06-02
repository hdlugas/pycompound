
##### Similarity Score Functions #####
# Note that the input for all similarity measures are two 1-d np arrays of the same length. 
# These 1-d arrays must be normalized to sum to 1 for the Shannon, Renyi, and Tsallis Entropy Similarity Measures.

import scipy.stats
import numpy as np


def S_cos(ints_a, ints_b):
    # Cosine Similarity Measure
    if np.sum(ints_a) == 0 or np.sum(ints_b) == 0:
        return(0)
    else:
        return np.dot(ints_a,ints_b) / (np.sqrt(sum(np.power(ints_a,2))) * np.sqrt(sum(np.power(ints_b,2))))


def ent_renyi(ints, q):
    # Computes the Renyi entropy of a probability distribution for a given positive entropy dimension q
    return np.log(sum(np.power(ints,q))) / (1-q)


def ent_tsallis(ints, q):
    # Computes the Tsallis entropy of a probability distribution for a given positive entropy dimension q
    return (sum(np.power(ints,q))-1) / (1-q)


def S_shannon(ints_a, ints_b):
    '''
    Shannon Entropy Similarity Measure

    This similarity function was presented by: 
    Li, Y.; Kind, T.; Folz, J.; Vaniya, A.; Mehta, S. S.; Fiehn, O.
    Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. 
    * Note that since scipy.stats.entropy normalizes the input vector to sum to 1, vec1 and vec1 need not be normalized when computing ent_ab
    '''

    ent_a = scipy.stats.entropy(ints_a)
    ent_b = scipy.stats.entropy(ints_b)
    ent_ab = scipy.stats.entropy(ints_a + ints_b)
    return(1 - (2 * ent_ab - ent_a - ent_b)/np.log(4))


def S_renyi(ints_a, ints_b, q):
    '''
    Renyi Entropy Similarity Measure
    * This is a novel similarity measure which generalizes the Shannon Entropy Similarity Measure
    * The Renyi Similarity Measure approaches the Shannon Entropy Similiarity Measure as q approaches 1
    * ints_a and ints_b must be normalized to sum to 1
    '''
    if q == 1:
        print('Warning: the Renyi Entropy Similarity Measure is equivalent to the Shannon Entropy Similarity Measure when the entropy dimension is 1')
        return S_shannon(ints_a, ints_b)
    else:
        ent_a = ent_renyi(ints_a, q)
        ent_b = ent_renyi(ints_b, q)
        ent_merg = ent_renyi(ints_a/2 + ints_b/2, q)
        N = (1/(1-q)) * (2*np.log(np.sum(np.power(ints_a/2,q))+np.sum(np.power(ints_b/2,q))) - np.log(np.sum(np.power(ints_a,q))) - np.log(np.sum(np.power(ints_b,q))))
        return 1 - (2 * ent_merg - ent_a - ent_b) / N


def S_tsallis(ints_a, ints_b, q):
    '''
    Tsallis Entropy Similarity Measure
    * This is a novel similarity measure which generalizes the Shannon Entropy Similarity Measure
    * The Tsallis Similarity Measure approaches the Shannon Entropy Similiarity Measure as q approaches 1
    * ints_a and ints_b must be normalized to sum to 1
    '''
    if q == 1:
        print('Warning: the Tsallis Entropy Similarity Measure is equivalent to the Shannon Entropy Similarity Measure when the entropy dimension is 1')
        return S_shannon(ints_a, ints_b)
    else:
        ent_a = ent_tsallis(ints_a, q)
        ent_b = ent_tsallis(ints_b, q)
        ent_merg = ent_tsallis(ints_a/2 + ints_b/2, q)
        N = np.sum(2*np.power(ints_a/2,q)+2*np.power(ints_b/2,q)-np.power(ints_a,q)-np.power(ints_b,q)) / (1-q)
        return 1 - (2 * ent_merg - ent_a - ent_b) / N


