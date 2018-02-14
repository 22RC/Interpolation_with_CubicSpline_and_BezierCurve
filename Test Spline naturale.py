import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import Spline_Bezier as s_b

def f(t):
    #funzione da interpolare
    return np.exp(-t) * np.cos(2*np.pi*t)


x = np.arange(0.0, 5.0, 0.01) # [a,b] = [0,0.01,0.02,...,5] 
n = 20 # n° nodi di interpolazione
dist = (int(len(x)/n)) # distanza tra i nodi
y = [] 
for i in x:
    y.append(f(i)) # f(x) per ogni x in [a,b]
spline_natural = s_b.Spline(x,dist,f) # creazione della spline naturale
x_knots = [k for k in spline_natural.partition] # ascisse dei nodi di interpolazione
y_knots = [f(k) for k in spline_natural.partition] # ordinate dei nodi di interpolazione 
tck = interpolate.CubicSpline(x_knots, y_knots,bc_type='natural') # Scipyspline s(x)
y_scipy = tck(x) # calcolo s(x) per ogni x in [a,b]
plt.plot(x,y,'k',color='b') #real function
plt.plot(x_knots,y_knots,'o',color='y') # k_not
plt.plot(x,y_scipy,'k--',color='g') # spline scipy
plt.plot(spline_natural.result_x,spline_natural.result_y,'k:',color='r') #Myspline
plt.legend(['real','knot','ScipyNatural','MysplineNatural'],loc='best')
plt.show()
