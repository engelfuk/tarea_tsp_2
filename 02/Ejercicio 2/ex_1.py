import numpy as np

def multMatMat(mat1,mat2):
    """
        Method to multiply two matrix of 4x4
    """
    if mat1.shape != (4,4) or mat2.shape != (4,4):
        raise Exception('Not valid: Matrixs must be 4x4')
    return mat1.dot(mat2)

def multMatVect(mat,vect):
    if mat.shape != (4,4):
        raise Exception('Not Valid: Matrix must be 4x4')
    if vect.shape != (4,1):
        raise Exception('Not Valid: Vector must be 4x1')
    return mat.dot(vect)
  




