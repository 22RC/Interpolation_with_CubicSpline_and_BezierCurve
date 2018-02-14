import numpy as np
import matplotlib.pyplot as plt
import Spline_Bezier as s_b
import bezier


def f(t):
    #funzione da interpolare
    return np.exp(-t) * np.cos(2*np.pi*t)

x = np.arange(0.0, 5.0, 0.01) # [a,b] = [0,0.01,0.02,...,5] 
n = 10 # n° nodi di interpolazione
dist = (int(len(x)/n)) # distanza tra i nodi
y = [] 
for i in x:
    y.append(f(i)) # f(x) per ogni x in [a,b]
time = np.arange(0.0,1.01,0.01) #tempo
Bezier_curva = s_b.Bezier(x,dist,f,time)
nodes2 = np.asfortranarray([ [0.0,0.0] for p in range(len(Bezier_curva.controlPoint_x))])
for i in range(len(nodes2)):
    nodes2[i] = [Bezier_curva.controlPoint_x[i],Bezier_curva.controlPoint_y[i]]
curve2 = bezier.Curve.from_nodes(nodes2)
curve2.plot(500,color='b')
plt.plot(x,y,'o',Bezier_curva.result_x,Bezier_curva.result_y,'--')
plt.legend(['real','bèzier','MyBèzier'],loc='best')
plt.show()
