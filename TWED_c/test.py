import numpy as np
import TWED

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    A=np.array([[0],[0],[1],[1],[2],[3],[5],[2],[0],[1],[-0.1]])
    tA=np.array(list(i for i in range(len(A))))
    B=np.array([[0],[1],[2],[2.5],[3],[3.5],[4],[4.5],[5.5],[2],[0],[0],[.25],[.05],[0]])
    tB=np.array(list(i for i in range(len(B))))
    C=np.array([[4],[4],[3],[3],[3],[3],[2],[5],[2],[.5],[.5],[.5]])
    tC=np.array(list(i for i in range(len(C))))
    nu=.1
    _lambda=.2
    degree=1
    print("twed(A,B,nu,lambda)=", TWED.distance(A, tA, B, tB, nu, _lambda, degree))
    print("twed(A,C,nu,lambda)=", TWED.distance(A, tA, C, tC, nu, _lambda, degree))
    print("twed(B,C,nu,lambda)=", TWED.distance(B, tB, C, tC, nu, _lambda, degree))
        
    plt.plot(A, label='A')
    plt.plot(B, label='B')
    plt.plot(C, label='C')
    plt.legend()
    plt.show()
