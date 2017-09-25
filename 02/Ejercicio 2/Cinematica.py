import math as m
import numpy as np
from sympy import *
from ex_1 import multMatMat

class Cinematica(object):
    """
    Imprementacion directa de los parametros de Denavit-Hartenberg

    ---------------------
    Parametros
    ---------------------
    -N : Numero de grados de Libertad
    """
       
    def __init__(self,N):
        self.__N=N
        self.__A0n = np.eye(N)


    
    
    

    def compute_dh(self, tablaDH):
        """
        Devuelve la Matriz de Transformacion de A0n, dada una matriz con los parametros de DH
        ----------------------
        Tabla DH
        ----------------------
        Orden de los parametros
        theta_i = 0 # rad
        d_i = 1 
        alpha_i-1 = 0 #rad
        a_i-1 = 1
        
        """
        
         
        matAux = np.eye(self.__N)
        for x in range(0,self.__N):
            
            c_th,s_th = cos(tablaDH[x,0]),sin(tablaDH[x,0])
            d = tablaDH[x,1]
            c_al,s_al = cos(tablaDH[x,2]),sin(tablaDH[x,2])
            a = tablaDH[x,3]
            
            matT = np.array([
                [c_th, -c_al*s_th, s_al*s_th , a*c_th],
                [s_th, c_al*c_th , -s_al*c_th, a*s_th],
                [0   , s_al      , c_al      ,      d],
                [0   ,          0,          0,      1]
                ])
            
            matAux = multMatMat(matAux,matT)
            
        self.__A0n = trigsimp(matAux)
        return self.__A0n


    
        

def compute_jacobian(tabla_DH):
        """
        Calcula el Jacobiano de dada una matriz de transformacion
        -------------
        Parametros
        -------------
        - mat_trans : matrix de transformacion (4x4)      
        """

        

        th, th1, th2, th3, th4 ,al = symbols("th th1 th2 th3 th4 al")
        l1, l2, l3, l4 = symbols("l1 l2 l3 l4")
        
        A01 = mat_transformacion(0,1,tabla_DH)
        A02 = mat_transformacion(0,2,tabla_DH)
        A03 = mat_transformacion(0,3,tabla_DH)
        A04 = mat_transformacion(0,4,tabla_DH)

        matPos = A04[0:3,3]
        
        dif_th1 = np.array([diff(matPos[0],th1),diff(matPos[1],th1),diff(matPos[2],th1)])
        dif_th2 = np.array([diff(matPos[0],th2),diff(matPos[1],th2),diff(matPos[2],th2)])
        dif_th3 = np.array([diff(matPos[0],th3),diff(matPos[1],th3),diff(matPos[2],th3)])
        dif_th4 = np.array([diff(matPos[0],th4),diff(matPos[1],th4),diff(matPos[2],th4)])
        
        z00 = np.array([0,0,1])
        z01 = A01[0:3,2]
        z02 = A02[0:3,2]
        z03 = A03[0:3,2]
        
        
        jacobiano = np.array([
            [dif_th1[0], dif_th2[0], dif_th3[0], dif_th4[0]],
            [dif_th1[1], dif_th2[1], dif_th3[1], dif_th4[1]],
            [dif_th1[2], dif_th2[2], dif_th3[2], dif_th4[2]],
            [z00[0],z01[0],z02[0],z03[0]],
            [z00[1],z01[1],z02[1],z03[1]],
            [z00[2],z01[2],z02[2],z03[2]]
            ])
        
        return jacobiano
        
def mat_transformacion(sist_1, sist_2, tablaDH):
    """
    Calcula las matrices de transformacion Axy

    ---------------
    Parametros
    --------------
    - sist_1 = inicio
    - sist_2 = fin
    - tablaDH = matrix con los parametros de DH 
    """
    matAux = np.eye(4)
    for x in range(sist_1,sist_2):
            
        c_th,s_th = cos(tablaDH[x,0]),sin(tablaDH[x,0])
        d = tablaDH[x,1]
        c_al,s_al = cos(tablaDH[x,2]),sin(tablaDH[x,2])
        a = tablaDH[x,3]
            
        matT = np.array([
             [c_th, -c_al*s_th, s_al*s_th , a*c_th],
             [s_th, c_al*c_th , -s_al*c_th, a*s_th],
             [0   , s_al      , c_al      ,      d],
             [0   ,          0,          0,      1]
             ])
           
        matAux = multMatMat(matAux,matT)
    matAux = trigsimp(matAux)
    return matAux
        






    








