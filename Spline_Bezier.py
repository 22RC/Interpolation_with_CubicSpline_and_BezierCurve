import numpy as np
from scipy import linalg
import math

class Spline():
    def __init__(self,interval , n_Knot  , function  , target=('natural',None)):
        self.dist = n_Knot # distanza nodi
        self.f = function # funzione da approssimare 
        self.partition , self.index_n_Knot = self.partition_interval(interval,self.dist)     
        self.y_newData = [ self.f(x) for x in self.partition] # nodo (x_i,f(x_i))
        self.h = self.h_i(self.partition) # distanze
        self.n_knot = len(self.partition)-1 # numero di nodi
        if target[0] == 'natural' : # spline naturale
            self.m_i = self.momenti_naturalSpline(self.n_knot,self.y_newData,self.h)
        elif target[0] == 'complete': #spline completa
            if target[1] == None: 
                raise ReferenceError
            else:
                f1_a,f1_b = target[1][0][1] , target[1][1][1] #valore derivata agli estremi
                self.m_i = self.momenti_splineCompleta(self.n_knot+1,self.y_newData,self.h,f1_a,f1_b)
        elif target[0] == 'not-and-knot': # spline not-and-knot
            self.m_i = self.momenti_splineNot_and_knot(self.n_knot,self.y_newData,self.h)
        self.beta , self.alfa = self.beta_i_alfa_i(self.m_i,self.y_newData,self.n_knot,self.h)
        self.result_x , self.result_y = self.interpola_nodi(interval) #calcolo s(x) perogni x in [a,b]
   
    def partition_interval (self,interval,n):
        '''Partiziona l'intevallo in n+1 nodi equidistanti
        Args: interval = intervallo [a,b] da partizionare
              n = distanza tra i nodi
        Return: tavola dei nodi di interpolazione e indici dei nodi'''
        i = 0 
        part = []
        idx = []
        while i < len(interval):
            part.append(float(interval[i]))
            idx.append(i)
            i += n
        if part.count(interval[len(interval)-1]) == 0: # aggiungo l'estremo dell'intervallo se non
            part.append(float(interval[len(interval)-1])) # e' stato incluso nella tavola
            idx.append(len(interval)-1)
        return part , idx
    
    def h_i (self,x):
        '''Calcola il vettore che contiene le distanze tra i nodi x_i , x_i+1
           Args: x = vettore delle ascisse dei nodi di interpolazione
           Return: vettore delle distanze'''
        h_i = []
        for i in range(len(x)-1):
            h_i.append(x[i+1]-x[i])
        return h_i   

    def momenti_naturalSpline(self,n,ydata,h):
        '''Calcola i momenti relativi alle spline naturali
           Args: n = numero di nodi
                 ydata = vettore delle ordinate dei nodi di interpolazione
                 h = vettore delle distanze tra nodi
           Return: vettore dei momenti per le spline naturali'''
        M = [[0.0 for x in range(n-1)] for y in range(n-1)] #matrice (n-1)x(n-1)
        for i in range(n-1): 
            for j in range(n-1): 
                if i == j :
                    M[i][j] = 2*(h[i]+h[i+1])
                elif i+1 == j :
                    M[i][j] = h[j]
                elif i == j+1 :
                    M[i][j] = h[i]
        vect_f = [] 
        for i in range(1,n):
            vect_f.append(6*(((ydata[i+1]-ydata[i])/h[i])-((ydata[i]-ydata[i-1])/h[i-1])))
        solution = linalg.solve(M,vect_f)
        momenti = []
        for i in solution:
            momenti.append(i)
        momenti.insert(0,0.0) # aggiungo mu_0
        momenti.insert(n,0.0) # aggiungo mu_n
        return momenti
    
    def momenti_splineCompleta(self,n,ydata,h,f1a,f1b):
        '''Calcola i momenti relativi alle spline complete
           Args: n = numero di nodi
                 ydata = vettore delle ordinate dei nodi di interpolazione
                 partition = vettore delle ascisse dei nodi di interpolazione
                 h = vettore delle distanze tra nodi
                 f1a,f1b = valori della derivata agli estremi
           Return: vettore dei momenti per le spline complete'''        
        M = [[0.0 for x in range(n)] for y in range(n)] # matrice (n+1)x(n+1)
        for i in range(n): 
            for j in range(n): 
                if i == j :
                    if i != 0 and i!= n-1:
                        M[i][j] = 2*(float(h[i]+h[i-1]))
                    elif i == 0:
                        M[i][j] = 2*(float(h[i]))
                    else:
                        M[i][j] = 2*(float(h[i-1]))
                elif i+1 == j :
                    M[i][j] = float(h[j-1])
                elif i == j+1 :
                    M[i][j] = float(h[i-1])
        vect_f = []
        vect_f.append(6*(((ydata[1]-ydata[0])/h[0])-f1a))
        for i in range (1,n-1):
            vect_f.append(6*(((ydata[i+1]-ydata[i])/h[i])-((ydata[i]-ydata[i-1])/h[i-1])))
        vect_f.append(6*(f1b-((ydata[n-1]-ydata[n-2])/h[n-2])))
        solution = linalg.solve(M,vect_f)
        return solution    
    
    def momenti_splineNot_and_knot(self,n,ydata,h):
        '''Calcola i momenti relativi alle spline not-and-knot
           Args: n = numero di nodi
                 ydata = vettore delle ordinate dei nodi di interpolazione
                 h = vettore delle distanze tra nodi
           Return: vettore dei momenti per le spline not-and-knot'''        
        M = [[0.0 for x in range(n-1)] for y in range(n-1)] #matrice (n-1)x(n-1)
        for i in range(n-1): 
            for j in range(n-1): 
                if i == j :
                    if i != 0 and i!= n-1:
                        M[i][j] = 2*(float(h[i]+h[i+1]))
                    elif i == 0:
                        M[i][j] = 3*(float(h[i])) + 2*(float(h[i+1])) + (float(pow(h[i],2)/h[i+1]))
                    else:
                        M[i][j] = 3*(float(h[i])) + 2*(float(h[i-1])) + (float(pow(h[i],2)/h[i-1]))
                elif i+1 == j :
                    if i == 0:
                        M[i][j] = float(h[j]) - float(pow(h[i],2)/h[j])
                    else:
                        M[i][j] = float(h[j])
                elif i == j+1 :
                    if i == n-1:
                        M[i][j] = float(h[j]) - float(pow(h[i],2)/h[j])
                    else:
                        M[i][j] = float(h[i])
        vect_f = [] 
        for i in range(1,n):
            vect_f.append(6*(((ydata[i+1]-ydata[i])/h[i])-((ydata[i]-ydata[i-1])/h[i-1])))
        solution = linalg.solve(M,vect_f)
        momenti = []
        mu_0 = solution[0]*(1+(float(h[0]/h[1]))) - solution[1]*(float(h[0]/h[1]))
        momenti.append(mu_0)
        lenght = len(solution)-1
        for i in solution:
            momenti.append(i)
        mu_n = solution[lenght]*(1-(float(h[lenght-1]/h[lenght])))-solution[lenght-1]*(float(h[lenght-1]/h[lenght]))
        momenti.append(mu_n)
        return momenti        
    def beta_i_alfa_i(self,mu_i,y,n,h):
        '''Calcola il vettore degli alpha_i e beta_i
           Args: mu_i = momenti che definiscono la spline
              y = vettore delle ordinate dei nodi di intrpolazione
              n = numero di nodi di interpolazione 
              h = distanza tra x_i e x_i+1
        Return: vettore degli alpha_i beta_i'''
        beta_i , alfa_i = [] , []
        for i in range(n):
            beta_i.append(y[i]-(mu_i[i]*(h[i]**2)/6))
            alfa_i.append(((y[i+1]-y[i])/h[i]) - ((h[i]/6)*(mu_i[i+1]-mu_i[i])))
        return beta_i ,alfa_i    

    def index_(self,xdata , x):
        '''Trova l'indice dell' i-esima spline s_i(x)
        Args: xdata = tupla dei nodi di interpolazione in cui il secondo elemento
                      serve per la valutazione nei nodi di raccordo e per stabilire
                      a quale intervallo appartiene il punto x
              x = punto in [a,b]
        Return: indice della spline che servira' per la valutazione in quell'intervallo'''
        trovato = -1
        i = 0
        while trovato == -1 and i < len(xdata)-1:
            if x == xdata[i][0] and i == 0 :
                trovato = i
            elif x == xdata[i][0] and i != 0 :
                if xdata[i][1] and i != len(xdata)-1:
                    trovato = i+1
                else:
                    xdata[i] = [xdata[i][0],True]
                    trovato = i
            elif xdata[i][0] < x < xdata[i+1][0]:
                trovato = i
            i+=1
        if i == len(xdata)-1:
            trovato = i
        return trovato 
    
    def interp_tratti(self, mu_i , alfa , beta , xdata , x , n ,lenght , h) :
        '''Trova l'indice dell'i-esima spline e ne restituisce il valore calcolato
           nel punto x
           Args: mu_i = momenti
                 alfa , beta = coefficienti dell'interpolazione
                 xdata = tuple di nodi di interpolazione
                 x = punto in [a,b]
                 n = distanza tra nodi
                 lenght = lunghezza intervallo
                 h = distanza tra x_i e x_i+1
            Returns: y = s(x)'''
        idx  = self.index_(xdata , x )
        if idx == -1:
            return self.eval_spline(mu_i,alfa,beta,x,xdata,n,0,h)
        else:
            return self.eval_spline(mu_i,alfa,beta,x,xdata,n,idx,h)    
    
    def eval_spline(self,mu_i , alpha , beta , x , n_x , n , idx , h):
        '''Calcola la spline s_i(x)'''
        if idx == len(n_x) -1: # s_n-1(x_n)
            g = mu_i[idx]*((x-n_x[idx-1][0])**3)/(6*h[idx-1])
            g1 = mu_i[idx-1] *((x-n_x[idx][0])**3)/(6*h[idx-1])
            a = alpha[idx-1]*(x-n_x[idx-1][0])
            b = beta[idx-1]
        else:
            g = mu_i[idx+1]*((x-n_x[idx][0])**3)/(6*h[idx])
            g1 = mu_i[idx] *((x-n_x[idx+1][0])**3)/(6*h[idx])
            a = alpha[idx]*(x-n_x[idx][0])
            b = beta[idx]
        return g - g1 + a + b    
    
    def interpola_nodi(self,interval):
        '''Calcolo di s(x) perogni x in [a,b]'''
        y = [] #risultato interpolazione
        data = [] # ascisse interpolazione 
        for i in range(len(interval)):
            if self.partition.count(float(interval[i])) > 0:
                if i != 0 and i != len(interval)-1: #se mi trovo nel nodo di raccordo
                    data.append(float(interval[i])) #aggiungo due volte alla lista 
                    data.append(float(interval[i])) 
                else:
                    data.append(float(interval[i]))
            else:
                data.append(float(interval[i]))
        k_not = []
        for i in self.partition: 
            k_not.append((i,False))
        for i  in data:
            y.append(self.interp_tratti(self.m_i,self.alfa,self.beta,k_not,i,self.dist,len(data),self.h))
        for i in self.index_n_Knot:
            if i != len(interval) -1 and i != 0:
                data.pop(i+1)
                y.pop(i+1)
        return data , y
        


class Bezier():
    def __init__(self,interval , n_Knot  , function  , time):
        self.n = n_Knot
        self.time = time
        self.partition , self.index_n_Knot = self.partition_interval(interval,self.n)
        self.partition.append(interval[len(interval)-1])
        
        self.f = function
        self.Matrix = self.init_system()
        self.controlPoint_x , self.controlPoint_y = self.find_controlPoint()
        self.result_x , self.result_y = self.bezierCurve()
        
    def init_system (self):
        M = [[0.0 for x in range(len(self.partition))] 
             for y in range(len(self.partition))] #curva di bezie
        for i in range(len(self.partition)):           #tolgo x0 , xn
            for j in range(len(self.partition)):
                if (i == 0 and j == 0) or (i == len(self.partition)-1 and j == len(self.partition)-1):
                    M[i][j] = 1.0
                elif i == 0 or i == len(self.partition)-1:
                    M[i][j] = 0.0
                else:
                    M[i][j] = self.binomial(len(self.partition)-1,j,i) #matrice n-1 per costrire i punti di controllo 

        return M

    def partition_interval (self,interval,n):
        #partizione intervallo
        i = 0 
        part = []
        idx = []
        while i < len(interval):
            part.append(float(interval[i]))
            idx.append(i)
            i += n
        return part , idx    
    
    def find_controlPoint(self):
        x_b , y_b = [] , []
        for x in self.partition:
            x_b.append(float(x))
            y_b.append(float(self.f(x)))
        x_i = linalg.solve(self.Matrix,x_b)
        y_i = linalg.solve(self.Matrix,y_b)
        
        return x_i , y_i
        
    
    def fattoriale(self,n,limit=1):
        if n == 0 or n == 1:
            return 1
        result = 1
        while n > limit:
            result *= n
            n -= 1
        return result
    
    def berstain_polynomial (self, n_b , k , t ) :
        binomial = (float(self.fattoriale(n_b)) / (float(self.fattoriale(k)) * float(self.fattoriale(n_b-k))))
        return ((binomial * pow(t,n_b-k)) * pow((1-t),(k)))
    
    def binomial (self, n_b , j , i ):
        binomial = (self.fattoriale(n_b)) / (self.fattoriale(j) * self.fattoriale(n_b-j))
        k = i
        p = j
        l = (float(k) / float(n_b))
        c = pow(l,p)
        d =  pow((1-l),(n_b-j))     
        return binomial * c * d    
    
    def bezier (self,x_b , n_b  , y_b , t):
        somma_x = 0
        somma_y = 0
        for i in range(n_b):
            somma_x += (x_b[i] * self.berstain_polynomial(n_b-1,i,t))
            somma_y += (y_b[i] * self.berstain_polynomial(n_b-1,i,t))
        return somma_x , somma_y 
    
    def bezierCurve(self):
        y_bezier , r_x , r_y = [] , [] , []
        for i in range(len(self.time)):
            y_bezier.append(self.bezier(self.controlPoint_x,len(self.partition),self.controlPoint_y,self.time[i]))
        for i , j in y_bezier:
            r_x.append(float(i))
            r_y.append(float(j))
        
        return r_x , r_y
