###############################################################################
###############################################################################
#                                quantumCh.py
#
#       quantumCh.py is a Python module that contains quantum mechanics
#       and quantum channels related functions.
#
#       Collaborators:
#           1. Jose Alfredo de León, deleongarrido.jose@gmail.com
#           2. Alejandro Barillas
#           3. Amado Cabrera
#
#       Feel free to contact any of the collaboratos in case you have any question.
#
#
#       Last updated: 09.02.2021
###############################################################################
###############################################################################

import numpy as np
from itertools import product
from math import log

###############################################################################
#                    QUANTUM MECHANICS-RELATED FUNCTIONS
#
#      1. Dirac(vector, nomenclature): prints vector in Dirac notation using
#         computational (10) or plus-minus (+-) nomenclature.
#      2. CPTest: CPTest(E) tests the complete positivity of a matrix E.
###############################################################################

def Comm(A, B):
    '''
    Comm(A, B) computes the commutator of A and B.

    ### __Parameters:__
    A, B: 2D array_like
        Operators acting on a Hilbert space.

    ### __Returns:__
    out: 2D array_like
        AB - BA

    _Author: J.A. de León_
    '''
    return np.matmul(A, B) - np.matmul(B, A)

def Dirac(vector, nomenclature='10'):
    '''
    Dirac(vector, nomenclature) prints _vector_ in Dirac notation using
    computational (10) or plus-minus (+-) nomenclature for state vectors
    of a system of spin-1/2 particles.

    ### __Parameters:__
    vector: 1D array_like
        State vector of a system of spin-1/2 particles.

    nomenclature (kwarg): default '10' or '+-'
        Nomenclature used to print kets.

    ### __Returns:__
    out: string
        Vector in Dirac notation.

    ### __Example:__
    Dirac(np.array([1,-4,2,-3])) = |00> - 4|01> + 2|10> - 3|11>.

    Dirac(np.array([1,-4,2,-3]), nomenclature='+-') = |++> - 4|+-> + 2|-+> - 3|-->

    _Author: J.A. de León_
    '''
    particles = int(log(len(vector), 2))
    strings = [];
    for i in np.argwhere(vector).flatten():
        if len(np.binary_repr(i)) < particles:
          strings.append((particles-len(np.binary_repr(i)))*"0" + str(np.binary_repr(i)))
        else:
          strings.append(str(np.binary_repr(i)))

    PlusMinusNotation = []

    for string in strings:
        dummy = string.replace("0","+")
        PlusMinusNotation.append(dummy.replace("1", "-"))

    if nomenclature == '+-':
        stringsToUse = PlusMinusNotation
    else:
        stringsToUse = strings

    k = 0
    for i in range(len(strings)):
        if k != 0:
          if vector[int(strings[i], 2)] > 0:
            print(" + ", end="")
          elif vector[int(strings[i], 2)] < 0:
            print(" - ", end="")
        if vector[int(strings[i], 2)] == 1:
          coef = ""
        else:
          coef = str(vector[int(strings[i], 2)])
        if str(vector[int(strings[i], 2)])[0] == "-":
          coef = str(vector[int(strings[i], 2)])[1:]
          if vector[int(strings[i], 2)] == -1:
            coef = ""
        print(coef + "|" + stringsToUse[i] + ">", end="")
        if strings.index(strings[i]) == len(strings)-1:
          break
        k += 1
    return

def Pauli(indices):
    '''
    Pauli(indices) computes the tensor product of Pauli matrices with indices
    in _indices_.

    ### __Parameters:__
    indices: 1D array_like
        Indices of Pauli matrices.

    ### __Returns:__
    out: 2D array_like
        Tensor product of Pauli matrices.

    ### __Examples:__
    Pauli([1,2]) = sigma_x $tensor$ sigma_y

    Pauli([2,3,0]) = sigma_y $tensor$ sigma_z $tensor$ identity

    _Author: J.A. de León_
    '''
    # Pauli matrices:
    s1 = np.array([[0,1],[1,0]])        # sigma_x
    s2 = np.array([[0,-1j],[1j,0]])     # sigma_y
    s3 = np.array([[1,0],[0,-1]])       # sigma_z

    # Identity and Pauli matrices in an array
    sigma = np.array([np.identity(2),s1,s2,s3])

    if len(indices) == 0:
        return 1
    else:
        return np.kron(sigma[indices.pop(0)], Pauli(indices))

def PTest(A):
    '''
    PTest(A) returns True if A is a positive-semidefinite matrix.

    ### __Parameters:__
    A: 2D square array_like
        Matrix.
    ### __Returns:__
    out: boolean
        True if positive-semidefinite holds, false if not.

    _Author: J.A. de León_
    '''
    # Calculate the smallest eigenvalue of A
    smallest_eig = chop(np.linalg.eigh(A)[0])[0]

    # Compute the logical value of the positivity
    positivity = smallest_eig >= 0
    return positivity

def HSInnerP(A, B):
	'''
    HSInnerP(A, B) computes the Hilbert-Schmidt inner product of A and B.

    ### __Parameters:__
    A, B: 2D array_like
        Operators acting on Hilbert-Schmidt space.

    ### __Returns:__
    out: float
        HS inner product between A and B.

    _Author: A. Barillas_
    '''

    # Definition of HS inner product
	output = np.trace(np.matmul(np.transpose(np.conjugate(A)),B))
	return output

def Reshuffle(E):
    '''
    Reshuffle(A) reshuffles A, with A a superoperator acting on density matrices
    of a system of spin 1/2 particles.

    ### __Parameters__:
    A: 2D array_like
        Superoperator acting on density matrices of a system of spin 1/2 particles.

    ### __Returns__:
    out: 2D array_like
        Choi matrix of A, A^R.

    _Author: J.A. de León_
    '''
    # FALTA REVISAR ESTOS COMENTARIOS
    # Calculate number of qubits n and dim(total Hilbert space)
    n = int(log(E.shape[1], 4))
    dimHilbert = 2**n

    # Create a 2D square array to storage Choi matrix
    choiMatrix = np.zeros((4**n,4**n), dtype=complex)

    # Compute the reshuffling of A
    k = 0 # qbit number label (?)
    p = 0 # index of the rows after every row in A has been reshaped in a square matrix of size 2**n
    for i in range(4**n): # Go over each of the 4**n rows of D
        if i % dimHilbert == 0 and i != 0: # If it has completed a dimHilbert number of rows of A
          k = k + 1
          p = 0
        for j in range(dimHilbert): # Each iteration goes over 4**n elements in a row of D
          choiMatrix[i][dimHilbert*j:dimHilbert*(j+1)] = E[dimHilbert*k+j].reshape((dimHilbert,dimHilbert))[p]
        p = p + 1

    return choiMatrix

###############################################################################
#                EXCLUSIVELY QUANTUM CHANNELS RELATED FUNCTIONS

#      1. PauliChannel: PauliChannel(diagonal) computes the superoperator
#         of a Pauli quantum channel in computational basis, given its
#         diagonal elements oh the superoperator in the tensor product of
#         Pauli matrices basis.
#      2. CPTest: CPTest(E) tests the complete positivity of a matrix E.
###############################################################################


def PauliChannel(diagonal):
    '''
    PauliChannel(diagonal) computes the superoperator of a Pauli quantum channel
    in computational basis, given its diagonal elements oh the superoperator
    in the tensor product of Pauli matrices basis.

    ### __Parameters:__
    diagonal: 1D array_like
        Diagonal of superoperator in Pauli tensor products.

    ### __Returns:__
    PauliChannel: 2D array_like
        Superoperator in computational basis.

    _Author: J.A. de León_
    '''
    # Calculate number of qubits n
    # (Pauli Channel superoperator's dimension is 4**n)
    n = int(log(len(diagonal), 4))

    # Compute change of basis matrix from Pauli tensor products to computacional basis
    indices = list(product([0,1,2,3], repeat=n))
    pauliToComputational = []
    for i in range(len(indices)):
      pauliToComputational = np.append(pauliToComputational, Pauli(list(indices[i])).reshape(4**n))
    pauliToComputational = np.transpose(pauliToComputational.reshape((4**n,4**n)))

    # Compute Pauli channel in computacional basis by making a change of basis
    # (Pauli channels are diagonal in Pauli tensor products basis)
    PauliChannel = np.matmul(pauliToComputational, np.matmul(np.diagflat(diagonal), np.linalg.inv(pauliToComputational)))

    return PauliChannel

def CPTest(E):
    '''
    CPTest(E) tests the complete positivity of a matrix E.

    ### __Parameters__:
    E: 2D array_like
        Superoperator acting on density matrices.
    ### __Returns__:
    out: boolean
        True if CP holds, false if not.

    _Author: J.A. de León_
    '''
    # Get the smallest eigenvalue of Choi matrix E^R
    smallestEigval = min(Chop(np.linalg.eigvals(Reshuffle(E))))

    # Return true if smallest eigenvalue is greater or equal to zero, false otherwise
    if smallestEigval >= 0:
        return True
    else:
        return False

###############################################################################
#                TECHNICAL PYTHON PROGRAMMING RELATED FUNCTIONS
###############################################################################

def MatrixForm(A, limiter='p'):
    """
    MatrixForm(A) prints A in matrix form.

    ### __Parameters:__
    A: 2D square array_like
        Matrix.
    ### __Returns:__
    out: string
        A printed in matrix form.

    ### __Example:__
    MatrixForm(np.identity(2)) = ⎛ 1.0 0.0 ⎞
                                 ⎝ 0.0 1.0 ⎠

    _Author: A. Cabrera_
    """
    if limiter == 'p':
        ls, li = '⎛', '⎝'
        rs, ri = ' ⎞', ' ⎠'
    elif limiter == 'b':
        ls, li = '⎡', '⎣'
        rs, ri = ' ⎤', ' ⎦'
    elif limiter == 'v':
        ls, li = '⎢', '⎢'
        rs, ri = ' ⎥', ' ⎥'
    else:
        ls, li = '⎡', '⎣'
        rs, ri = ' ⎤', ' ⎦'
    max = 0
    for row in A:
        for ele in row:
            if len(str(ele)) > max:
                max = len(str(ele))
    print_str = ''
    for i in range(len(A)):
        if i == 0:
            print_str += f'{ls}'
        elif i==(len(A)-1):
            print_str += f'{li}'
        else:
            print_str += '⎢'
        for j in range(len(A[i])):
            print_str += ' {val:^{max}}'.format(val=A[i][j], max=max)
        if i == 0:
            print_str += f'{rs}\n'
        elif i==(len(A)-1):
            print_str += f'{ri}'
        else:
            print_str += ' ⎥\n'
    print(print_str, end='\n\n')
    return

# When diagonalizing with Python one can encounter eigenvalues of the order
# 10^{-15} that are zeros. Chop() appproximates those eigenvalues to 0.
def Chop(expr, *, max=1e-10):
    '''
    Chop() appproximates numbers of order _max_ to 0.

    ### __Parameters__:
    expr: float, array_like
        Number or array of numbers to chop.
    ### __Returns__:
    out: float, array_like
        Number or array of numbers chopped.

    _Author: J.A. de León_
    '''
    return [i if abs(i) > max else 0 for i in expr]
