'''
    JAquantum.py

    September 5, 2020
    José Alfredo de León

    This is a module for useful functions to compute quantum
    calculations.
'''


import numpy as np

# Pauli matrices:
s1 = np.array([[0,1],[1,0]])        # sigma_x
s2 = np.array([[0,-1j],[1j,0]])     # sigma_y
s3 = np.array([[1,0],[0,-1]])       # sigma_z

# Array with identity matrix and 3 Pauli matrices
sigma = np.array([np.identity(2),s1,s2,s3])

# Eigenbasis of sigma_z
ket0 = np.array([[0],[1]])
ket1 = np.array([[1],[0]])
zEigenbasis = np.array([ket0,ket1])

def Pauli(indices):
    '''
    Recursive function to calculate the tensor product of
    Pauli matrices.

    Examples:
    Pauli([1,2]) = sigma_x tensor sigma_y
    Pauli([2,3,0]) = sigma_y tensor sigma_z tensor identity
    '''
    if len(indices) == 0:
        return 1
    else:
        return np.kron(sigma[indices.pop(0)], Pauli(indices))


def TensorZkets(indices):
    '''
    Recursive function to calculate the tensor product of
    the elements from the eigenbasis of sigma_z.

    Examples:
    TensorZkets([1,0]) = ket1 tensor ket0
    TensorZkets([1,1,1,1]) = ket1 tensor ket1 tensor ket1 tensor ket1
    '''
    if len(indices) == 0:
        return 1
    else:
        return np.kron(zEigenbasis[indices.pop(0)], TensorZkets(indices))


def Comm(A, B):
    '''
    This function calculates the commutator of two matrices A and B.

    Example:
    Comm(A,B) = AB - BA
    '''
    return np.matmul(A, B) - np.matmul(B, A)
