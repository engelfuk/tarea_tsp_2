from Cinematica import *

def main():
    
    th, th1, th2, th3, th4 ,al = symbols("th th1 th2 th3 th4 al")
    l1, l2, l3, l4 = symbols("l1 l2 l3 l4")
        
   

    tablaDH_1 = np.array([
        [th1,l3,0,l1],
        [th2,l4,0,l2],
        [th3,0,0,0],
        [0,th4,0,0]
        ])




    prueba = Cinematica(4)
    print("--------- Matriz de Transformacion Simbolica ----------------")
    print(prueba.compute_dh(tablaDH_1))
    print("--------- Jacobano ----------------")
    print(compute_jacobian(tablaDH_1))
    

main()










