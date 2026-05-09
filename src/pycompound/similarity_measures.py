
##### Similarity Score Functions #####
# Note that the input for all similarity measures is two 1-d np arrays of the same length. 
# These 1-d arrays must be normalized to sum to 1 for the Shannon, Rényi, and Tsallis Entropy Similarity Measures.

from .processing import *
import scipy.stats
import numpy as np
import sys



def S_cos(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute cosine similarity between two intensity vectors.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Cosine similarity score.
    """
    if np.sum(ints_a) == 0 or np.sum(ints_b) == 0:
        return(0)
    else:
        return np.dot(ints_a,ints_b) / (np.sqrt(sum(np.power(ints_a,2))) * np.sqrt(sum(np.power(ints_b,2))))


def ent_renyi(ints: np.ndarray, q: float) -> float:
    """Compute Rényi entropy for an intensity vector.

    Args:
        ints: Intensity vector.
        q: Rényi entropy dimension.

    Returns:
        Rényi entropy value.
    """
    return np.log(sum(np.power(ints,q))) / (1-q)


def ent_tsallis(ints: np.ndarray, q: float) -> float:
    """Compute Tsallis entropy for an intensity vector.

    Args:
        ints: Intensity vector.
        q: Tsallis entropy dimension.

    Returns:
        Tsallis entropy value.
    """
    return (sum(np.power(ints,q))-1) / (1-q)


def S_shannon(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Shannon entropy similarity between two intensity vectors.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Shannon entropy similarity score.
    """
    ent_a = scipy.stats.entropy(ints_a)
    ent_b = scipy.stats.entropy(ints_b)
    ent_ab = scipy.stats.entropy(ints_a + ints_b)
    return(1 - (2 * ent_ab - ent_a - ent_b)/np.log(4))


def S_renyi(ints_a: np.ndarray, ints_b: np.ndarray, q: float) -> float:
    """Compute Rényi entropy similarity between two intensity vectors.

    Args:
        ints_a: First normalized intensity vector.
        ints_b: Second normalized intensity vector.
        q: Rényi entropy dimension.

    Returns:
        Rényi entropy similarity score.
    """
    if q == 1:
        print('Warning: the Renyi Entropy Similarity Measure is equivalent to the Shannon Entropy Similarity Measure when the entropy dimension is 1')
        return S_shannon(ints_a, ints_b)
    else:
        ent_a = ent_renyi(ints_a, q)
        ent_b = ent_renyi(ints_b, q)
        ent_merg = ent_renyi(ints_a/2 + ints_b/2, q)
        N = (1/(1-q)) * (2*np.log(np.sum(np.power(ints_a/2,q))+np.sum(np.power(ints_b/2,q))) - np.log(np.sum(np.power(ints_a,q))) - np.log(np.sum(np.power(ints_b,q))))
        return 1 - (2 * ent_merg - ent_a - ent_b) / N


def S_tsallis(ints_a: np.ndarray, ints_b: np.ndarray, q: float) -> float:
    """Compute Tsallis entropy similarity between two intensity vectors.

    Args:
        ints_a: First normalized intensity vector.
        ints_b: Second normalized intensity vector.
        q: Tsallis entropy dimension.

    Returns:
        Tsallis entropy similarity score.
    """
    if q == 1:
        print('Warning: the Tsallis Entropy Similarity Measure is equivalent to the Shannon Entropy Similarity Measure when the entropy dimension is 1')
        return S_shannon(ints_a, ints_b)
    else:
        ent_a = ent_tsallis(ints_a, q)
        ent_b = ent_tsallis(ints_b, q)
        ent_merg = ent_tsallis(ints_a/2 + ints_b/2, q)
        N = np.sum(2*np.power(ints_a/2,q)+2*np.power(ints_b/2,q)-np.power(ints_a,q)-np.power(ints_b,q)) / (1-q)
        return 1 - (2 * ent_merg - ent_a - ent_b) / N


def S_mixture(ints_a: np.ndarray, ints_b: np.ndarray, weights: dict[str, float] = {'Cosine':0.25, 'Shannon':0.25, 'Renyi':0.25, 'Tsallis':0.25}, q: float = 1.1) -> float:
    """Compute a weighted mixture of supported similarity measures.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.
        weights: Weights for Cosine, Shannon, Renyi, and Tsallis similarities.
        q: Entropy dimension used for Rényi and Tsallis similarities.

    Returns:
        Weighted mixture similarity score.
    """
    if set(weights.keys()).issubset(set(['Cosine','Shannon','Renyi','Tsallis'])) is False:
        print('Error: the keys to the weight parameter dict of the function S_mixture must be one of the four: Cosine, Shannon, Renyi, Tsallis')
        sys.exit()

    similarity = 0
    for key, value in weights.items():
        if key == 'Cosine':
            similarity += value * S_cos(ints_a,ints_b)
        if key == 'Shannon':
            similarity += value * S_shannon(ints_a,ints_b)
        if key == 'Renyi':
            similarity += value * S_renyi(ints_a,ints_b,q)
        if key == 'Tsallis':
            similarity += value * S_tsallis(ints_a,ints_b,q)
    return similarity


def get_contingency_entries(ints_a: np.ndarray, ints_b: np.ndarray) -> list[int]:
    """Compute binary contingency entries for two intensity vectors.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        List containing counts [a, b, c], where a counts features present only in
        ints_a, b counts features present only in ints_b, and c counts features
        present in both.
    """
    a = 0
    b = 0
    c = 0

    for x, y in zip(ints_a, ints_b):
        if x != 0 and y != 0:
            c += 1
        elif x != 0 and y == 0:
            a += 1
        elif x == 0  and y != 0:
            b += 1
    return [a,b,c]


def S_jaccard(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Jaccard similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Jaccard similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = a + b + c
    if denom == 0:
        similarity = 0
    else:
        similarity = c / (a + b + c)
    return similarity


def S_dice(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Dice similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Dice similarity score.
    """

    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = a + b + 2 * c
    if denom == 0:
        similarity = 0
    else:
        similarity = 2 * c / denom
    return similarity


def S_3w_jaccard(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute three-weighted Jaccard similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Three-weighted Jaccard similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = a + b + 3 * c
    if denom == 0:
        similarity = 0
    else:
        similarity = 3 * c / denom
    return similarity


def S_sokal_sneath(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Sokal-Sneath similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Sokal-Sneath similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = 2 * a + 2 * b + c
    if denom == 0:
        similarity = 0
    else:
        similarity = c / denom
    return similarity


def S_binary_cosine(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute binary cosine similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Binary cosine similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = np.sqrt((a + c) * (b + c))
    if denom == 0:
        similarity = 0
    else:
        similarity = c / denom
    return similarity


def S_mountford(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Mountford similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Mountford similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = c * (a + b) + 2 * a * b
    if denom == 0:
        similarity = 1
    else:
        similarity = 2 * c / denom
    return similarity


def S_mcconnaughey(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute McConnaughey similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        McConnaughey similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = (a + c) * (b + c)
    if denom == 0:
        similarity = 0
    else:
        similarity = (c**2 - a * b) / denom
    return similarity


def S_driver_kroeber(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Driver-Kroeber similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Driver-Kroeber similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = 2 * (a + c) * (b + c)
    if denom == 0:
        similarity = 0
    else:
        similarity = c * (a + b + 2 * c) / denom
    return similarity


def S_simpson(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Simpson similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Simpson similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = min(a + c, b + c)
    if denom == 0:
        similarity = 0
    else:
        similarity = c / denom
    return similarity


def S_braun_banquet(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Braun-Banquet similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Braun-Banquet similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = max(a + c, b + c)
    if denom == 0:
        similarity = 0
    else:
        similarity = c / denom
    return similarity


def S_fager_mcgowan(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Fager-McGowan similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Fager-McGowan similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom1 = np.sqrt((a + c) * (b + c))
    denom2 = 2 * np.sqrt(max(a + c, b + c))
    if denom1 == 0 or denom2 == 0:
        similarity = 0
    else:
        similarity = c / denom1 - 1 / denom2
    return similarity


def S_kulczynski(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Kulczynski similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Kulczynski similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    denom = a + b
    if denom == 0:
        similarity = 1
    else:
        similarity = c / denom
    return similarity


def S_intersection(ints_a: np.ndarray, ints_b: np.ndarray) -> int:
    """Compute the number of shared nonzero entries.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Count of entries that are nonzero in both vectors.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    c = tmp[2]
    return c


def S_hamming(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Hamming-style similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Hamming-style similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    denom = a + b
    if denom == 0:
        similarity = 1
    else:
        similarity = 1 / denom
    return similarity


def S_hellinger(ints_a: np.ndarray, ints_b: np.ndarray) -> float:
    """Compute Hellinger-style similarity from binary presence patterns.

    Args:
        ints_a: First intensity vector.
        ints_b: Second intensity vector.

    Returns:
        Hellinger-style similarity score.
    """
    tmp = get_contingency_entries(ints_a, ints_b)
    a = tmp[0]
    b = tmp[1]
    c = tmp[2]
    similarity = 1 - np.sqrt((1 - c / np.sqrt((a + c) * (b + c))))
    return similarity


def get_similarity(similarity_measure: str, q_ints: np.ndarray, r_ints: np.ndarray, weights: dict[str, float] | None, q: float) -> float:
    """Dispatch to the requested similarity measure.

    Args:
        similarity_measure: Name of the similarity measure to compute.
        q_ints: Query intensity vector.
        r_ints: Reference intensity vector.
        weights: Weights used by the mixture similarity measure.
        q: Entropy dimension used by Rényi and Tsallis similarities.

    Returns:
        Requested similarity score.
    """

    if similarity_measure == 'cosine':
        similarity = S_cos(q_ints, r_ints)

    elif similarity_measure in ['shannon', 'renyi', 'tsallis']:
            q_ints = normalize(q_ints, method = 'standard')
            r_ints = normalize(r_ints, method = 'standard')
            if similarity_measure == 'shannon':
                similarity = S_shannon(q_ints, r_ints)
            elif similarity_measure == 'renyi':
                similarity = S_renyi(q_ints, r_ints, q)
            elif similarity_measure == 'tsallis':
                similarity = S_tsallis(q_ints, r_ints, q)

    elif similarity_measure == 'mixture':
        similarity = S_mixture(q_ints, r_ints, weights, q)

    elif similarity_measure == 'jaccard':
        similarity = S_jaccard(q_ints, r_ints)

    elif similarity_measure == 'dice':
        similarity = S_dice(q_ints, r_ints)

    elif similarity_measure == '3w_jaccard':
        similarity = S_3w_jaccard(q_ints, r_ints)

    elif similarity_measure == 'sokal_sneath':
        similarity = S_sokal_sneath(q_ints, r_ints)

    elif similarity_measure == 'binary_cosine':
        similarity = S_binary_cosine(q_ints, r_ints)

    elif similarity_measure == 'mountford':
        similarity = S_mountford(q_ints, r_ints)

    elif similarity_measure == 'mcconnaughey':
        similarity = S_mcconnaughey(q_ints, r_ints)

    elif similarity_measure == 'driver_kroeber':
        similarity = S_driver_kroeber(q_ints, r_ints)

    elif similarity_measure == 'simpson':
        similarity = S_simpson(q_ints, r_ints)

    elif similarity_measure == 'braun_banquet':
        similarity = S_braun_banquet(q_ints, r_ints)

    elif similarity_measure == 'fager_mcgowan':
        similarity = S_fager_mcgowan(q_ints, r_ints)

    elif similarity_measure == 'kulczynski':
        similarity = S_kulczynski(q_ints, r_ints)

    elif similarity_measure == 'intersection':
        similarity = S_intersection(q_ints, r_ints)

    elif similarity_measure == 'hamming':
        similarity = S_hamming(q_ints, r_ints)

    elif similarity_measure == 'hellinger':
        similarity = S_hellinger(q_ints, r_ints)

    return similarity

