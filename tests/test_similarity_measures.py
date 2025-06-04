
# this script runs some tests on the similarity measures to ensure they have the expected behavior when computing the similarity between (1) two completely dissimilar spectra and (2) two identical spectra
from similarity_measures import *
import numpy as np


intensities_A = np.array([0,0.2,0.1,0,0.5,0.2])
intensities_B = np.array([0.3,0,0,0.7,0,0])
print(f'\nIntensities of spectrum A: {intensities_A}')
print(f'Intensities of spectrum B: {intensities_B}')
print(f'Cosine similarity between spectrum A and spectrum B: {S_cos(intensities_A,intensities_B)}')
print(f'Shannon Entropy similarity between spectrum A and spectrum B: {S_shannon(intensities_A,intensities_B)}')
print(f'Renyi Entropy similarity between spectrum A and spectrum B with entropy dimension 1.5: {S_renyi(intensities_A,intensities_B,1.5)}')
print(f'Tsallis Entropy similarity between spectrum A and spectrum B with entropy dimension 1.5: {S_tsallis(intensities_A,intensities_B,1.5)}')


intensities_A = np.array([0.1,0.1,0.6,0,0,0.2])
intensities_B = intensities_A
print(f'\nIntensities of spectrum A: {intensities_A}')
print(f'Intensities of spectrum B: {intensities_B}')
print(f'Cosine similarity between spectrum A and spectrum B: {S_cos(intensities_A,intensities_B)}')
print(f'Shannon Entropy similarity between spectrum A and spectrum B: {S_shannon(intensities_A,intensities_B)}')
print(f'Renyi Entropy similarity between spectrum A and spectrum B with entropy dimension 1.5: {S_renyi(intensities_A,intensities_B,1.5)}')
print(f'Tsallis Entropy similarity between spectrum A and spectrum B with entropy dimension 1.5: {S_tsallis(intensities_A,intensities_B,1.5)}\n')


