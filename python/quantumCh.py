import numpy as np
from itertools import product, combinations
import sys

def combinations_of_components(n, num_of_corr):
	rIndices = list(range(4**n-1))

	positions = list(combinations(rIndices, num_of_corr))

	r = np.zeros((len(positions),4**n))

	for i in range(len(positions)):
		r[i][0] = 1
		for j in positions[i]:
			r[i][j+1] = 1
	return r


def diagonalMatrix(dim, diagEntries):
	# Constructs a diagonal square matrix.

	# dim : dimension of the matrix, int expected
	# diagEntries : list with diagonal entries, list expected

	if (type(dim) != int):
		sys.exit("Integer expected as first argument of diagonalMatrix function.")

	if (type(diagEntries) != list):
		sys.exit("List expected as second argument of diagonalMatrix function.")

	if (len(diagEntries) != dim):
		sys.exit("In diagonalMatrix function: length of diagonal-entries' list does not match dimension of the matrix.")

	A = np.zeros((dim, dim))

	for i in range(dim):
		A[i][i] = diagEntries[i]

	return A

def onesAndZeros(oneIndices, listLength):
	# Construct a list of 1's and 0's from tuple of indices of 1's.

	#if ( (type(oneIndices) != tuple) or (type(oneIndices) != list)):
	#	sys.exit("List or a tuple expected as argument of onesAndZeros.")

	List = [0]*(listLength)

	for i in oneIndices:
		List[i] = 1

	return List




'''
r = R(3, 3)

counter = 0

qubits_positions = np.zeros((3,3))

for i in range(3):
	for j in range(3):
		qubits_positions[i][j] = 4**i*(j+1)

for k in range(r.shape[0]):
	qubits_corr = np.zeros((3,3))
	norms = np.zeros(3)
	for i in range(3):
		for j in range(3):
			qubits_corr[i][j] = r[k][int(qubits_positions[i][j])]
		norms[i] = np.dot(qubits_corr[i],qubits_corr[i])
	for elem in list(product([0,1,2,3], repeat=3))[1:22]:
		#print(elem)
		if norms[0] == elem[0] and norms[1] == elem[1] and norms[2] == elem[2]:
			#print(r[k])
			counter = counter + 1
'''


# 1 2 3, 4 8 12, 16 32 48

def maps_in_Pauli_basis(n, param, num_of_corr):
    '''
    Constructs, in the Pauli basis, all the posible maps that leave invariant a certain qbit-system's number of components.

    Input:
    - n: number of qbits
    - param: number of maps to consider (usually depends on the number of components to leave invariant)
    - num_of_corr: number of components to leave invariant

    Returns an array with param number of all possible maps that leave invariant the num_of_corr n-qbit-system's number of components.
    '''
    # Storage all possible combinations of arrays of size 4**n with num_of_corr ones and the rest zeroes
    r = combinations_of_components(n, num_of_corr)

    # Create an array to storage the param number of maps (square matrices of size 4**n)
    maps_P = np.zeros((param, 4**n, 4**n))

    # Assign to every diagonal element of every map its correspondent element in each array of posible combinations of 1's and 0's
    for i in range(param):
        for j in range(4**n):
            maps_P[i][j][j] = r[i][j]

    # Return all possible maps in the Pauli basis
    return maps_P

def S(N,alpha,n):
	"""
	PAULI OPERATORS FUNCTION
	N: the number of qubits
	alpha: the Pauli matrix: 1->x, 2->y, 3->z
	n: label of the qubit
	"""

	# Pauli basis:
	sx = np.array([[0,1],[1,0]])
	sy = np.array([[0,-1],[1,0]])
	sz = np.array([[1,0],[0,-1]])
	sid = np.array([[1,0],[0,1]])
	s0 = np.array([[0,0],[0,0]])

	if alpha == 1: #if we define the x spin operator
		if n == 0: #if we consider the first site of the chain
			Sn = np.kron(sx,np.identity(2**(N-1)))
		else: #for the rest of the chain
			Sn = np.kron(np.kron(np.identity(2**(n-1)),sx),np.identity(2**(N-n)))
	elif alpha == 2: #if we define the y spin operator
		if n == 0: #if we consider the first site of the chain
			Sn = np.kron(sy,np.identity(2**(N-1)))
		else: #for the rest of the chain
			Sn = np.kron(np.kron(np.identity(2**(n-1)),sy),np.identity(2**(N-n)))
	elif alpha == 3: #if we define the z spin operator
		if n == 0: #if we consider the first site of the chain
			Sn = np.kron(sz,np.identity(2**(N-1)))
		else: #for the rest of the chain
			Sn = np.kron(np.kron(np.identity(2**(n-1)),sz),np.identity(2**(N-n)))
	elif alpha == 0:
		Sn = np.identity(2**N)
	return Sn

def change_of_basis_matrix_P_to_C(n):
	alpha = [0,1,2,3]
	R = np.array(list(product(alpha, repeat=(n))))

	# Create a 2D square array to storage in each row every vector of the Pauli basis
	change_of_basis_M = np.zeros((4**n, 4**n))

	# Actually asign to each row of change_of_basis_M every vector of the Pauli basis
	for i in range(4**n):
		dummy = np.identity(2**n)
		for j in range(n):
			dummy = dummy @ S(n,R[i][j],j+1)
		change_of_basis_M[i] = dummy.reshape(4**n)

	# Now that you have each column in every row of the change-of-basis matrix in change_of_basis_M transpose it to get the actual
	# change-of-basis matrix
	# M_cb_PaC = np.column_stack(change_of_basis_M)
	change_of_basis_M = change_of_basis_M.T

	return change_of_basis_M

def change_of_basis(A, change_of_basis_M):
	'''
	Makes a change of a basis of a matrix.

	Input:
	- A: matrix (2D square np.array)
	- change_of_basis_M: change-of-basis matrix (2D square np.array)

	Returns the matrix in the new basis.
	'''

	# Calculate de inverse of the change-of-basis matrix
	inverse_change_of_basis_M = np.linalg.inv(change_of_basis_M)

	# Compute de change of basis
	matrix_in_new_basis = np.matmul(change_of_basis_M, np.matmul(A, inverse_change_of_basis_M))

	return matrix_in_new_basis

def reshuffle(A, n):
	'''
	Computes the reshuffling of a square matrix of size 4**n.

	Input:
	- A: 2D squared np.array
	- n: number of qbits in the system

	Returns the dynamical matrix (also known as Choi matrix) of A.
	'''
	# Calculate the dimension of the density matrix of the qbit system
	rho_size = 2**n

	# Create a 2D square array to storage the dynamical matrix
	dynamical_Matrix = np.zeros((4**n,4**n))

	# Compute the reshuffling of A
	k = 0 # qbit number label (?)
	p = 0 # index of the rows after every row in A has been reshaped in a square matrix of size 2**n
	for i in range(4**n): # Go over each of the 4**n rows of D
		if i % rho_size == 0 and i != 0: # If it has completed a rho_size number of rows of A
			k = k + 1
			p = 0
		for j in range(rho_size): # Each iteration goes over 4**n elements in a row of D
			dynamical_Matrix[i][rho_size*j:rho_size*(j+1)] = A[rho_size*k+j].reshape((rho_size,rho_size))[p]
		p = p + 1

	return dynamical_Matrix

def chop(expr, *, max=1e-10):
	return [i if abs(i) > max else 0 for i in expr]

def positivity_test(A):
    '''
    Checks the positivity of a matrix.

    Input:
    - A: matrix

    Returns a boolean value. True if A is indeed positive and false if it's not.
    '''
    # Calculate the smallest eigenvalue of A
    smallest_eig = chop(np.linalg.eigh(A)[0])[0]

    # Compute the logical value of the positivity
    positivity = smallest_eig >= 0
    return positivity

def HSInnerP(A, B):
	'''
    Computes the Hilbert-Schimidt inner product.

    Input:
    - A, B: matrices of same dimension

    Returns Hilbert-Schimidt inner product of A and B.
    '''

	output = np.trace(np.matmul(np.transpose(np.conjugate(A)),B))
	return output
