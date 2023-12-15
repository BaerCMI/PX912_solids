# Program name: tools.py
# Description: Tools needed to answer later problems.
import sympy as sym
from sympy.matrices import Matrix, eye, zeros, ones, diag, MatMul


def deformation_gradient(phi, variables):    
    # Set up an empty matrix
    F = Matrix.zeros(len(variables))
    
    # Loop over transformed coordinates
    for i, x_i in enumerate(phi):
        # Loop over base coordinates
        for j, X_j in enumerate(variables):
    
            # Compute relevant derivative
            F[i,j] = sym.diff(x_i, X_j)
    
    return F

def infinitesimal_strain(phi, variables):
    
    F = deformation_gradient(phi, variables)
    
    # sym.eye is the identity matrix
    FmI = F - Matrix.eye(len(variables))
    
    # epsilon matrix as defined in notes
    epsilon = 0.5*(FmI + FmI.T)
    
    return epsilon