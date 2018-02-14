import numpy as np
import matplotlib.pyplot as plt
import Spline_Bezier as s_b
import bezier


def f(t):
    #funzione da interpolare
    return np.exp(-t) * np.cos(2*np.pi*t)

n = 10
t1 = np.arange(0.0, 5.0, 0.01) # [a,b] = [0,0.01,0.02,...,5] 
t2 = [] 
for i in t1:
    t2.append(f(i)) # f(x) per ogni x in [a,b]
time = np.arange(0.0,1.01,0.01) #tempo
Bezier_curva = s_b.Bezier(t1,n,f,time)
nodes2 = np.asfortranarray([ [0.0,0.0] for p in range(len(Bezier_curva.controlPoint_x))])
for i in range(len(nodes2)):
    nodes2[i] = [Bezier_curva.controlPoint_x[i],Bezier_curva.controlPoint_y[i]]
curve2 = bezier.Curve.from_nodes(nodes2)
curve2.plot(len(t1))
plt.plot(t1,f(t1),'o',Bezier_curva.result_x,Bezier_curva.result_y,'k')
plt.legend(['CurveModuleBezier','real','MyBezier'],loc='best')
plt.show()
