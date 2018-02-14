import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import Spline_Bezier as s_b

def f(t):
    #funzione da interpolare
    return np.exp(-t) * np.cos(2*np.pi*t)

def f_prime(t):
    #derivata
    return -2*np.pi*(np.exp(-t))*np.sin(2*np.pi*t) -((np.exp(-t))*(np.cos(2*np.pi*t)))

x = np.arange(0.0, 5.0, 0.01) # [a,b] = [0,0.01,0.02,...,5] 
n = 30 # n° nodi di interpolazione
dist = (int(len(x)/n)) # distanza tra i nodi
y = [] 
for i in x:
    y.append(f(i)) # f(x) per ogni x in [a,b]
f1_a , f1_b = f_prime(x[0]) , f_prime(x[len(x)-1]) # valori della derivata prima agli estremi
spline_completa = s_b.Spline(x,dist,f,target=['complete',[(1,f1_a),(1,f1_b)]]) # creazione della spline completa
x_knots = [k for k in spline_completa.partition] # ascisse dei nodi di interpolazione
y_knots = [f(k) for k in spline_completa.partition] # ordinate dei nodi di interpolazione 
tck = interpolate.CubicSpline(x_knots, y_knots,bc_type=((1,f1_a),(1,f1_b))) # rappresentazione Scipyspline s(x)
y_scipy = tck(x) # calcolo s(x) per ogni x in [a,b]
plt.plot(x,y,'k',color='b') #real function
plt.plot(x_knots,y_knots,'o',color='y') # k_not
plt.plot(x,y_scipy,'k--',color='g') # spline scipy
plt.plot(spline_completa.result_x,spline_completa.result_y,'k:',color='r') #Myspline
plt.legend(['real','knot','ScipyCompl','MysplineCompl'],loc='best')
plt.show()