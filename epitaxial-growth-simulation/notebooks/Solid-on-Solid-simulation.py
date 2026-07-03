#!/usr/bin/env python
# coding: utf-8

# In[6]:


# SOLID-ON-SOLID 
import numpy as np
def rand_sch(a,c,m,x):
    q = int(m/a)
    #q = round(m/a)
    #print("q:",q)
    r = m%a
    b = a*(x%q) - r*(x/q)
    if(b < 0):
        b += m
    if (b == 1 or b == 0):
        b = rand_sch(a,c,m,x)
    return b
def rand0(a,c,m,x):
    x = (a*x + c)%m
    return(x)
def rand_std(n1,n2):
    a1 = 40692
    c1 = 0
    m1 = 2**31 - 249
    a2 = 7**5
    c2 = 10
    m2 = 2**31 - 1
    u1 = ((n2-n1)%m1)/m1
    n2 = rand_sch(a2,c2,m2,n2)
    u2 = ((n2-n1)%m2)/m2
    n1 = rand0(a1,c1,m1,n1) 
    return(u1,u2,n1,n2)
def D(m,Y,X,transx,transy):
    donadores = []
    #NORTE,SUR,ESTE,OESTE
    N = [(transy + Y)%(Y+1),transx]
    S = [(transy + 1)%(Y+1),transx]
    E = [transy,(transx + X)%(X+1)]
    O = [transy,(transx + 1)%(X+1)]
    #print(m[N[0],N[1]])
    if (m[N[0],N[1]] > m[transy,transx]):
        donadores.append(N) 
    if (m[S[0],S[1]] > m[transy,transx]):
        donadores.append(S) 
    if (m[E[0],E[1]] > m[transy,transx]):
        donadores.append(E) 
    if (m[O[0],O[1]] > m[transy,transx]):
        donadores.append(O) 
        
    return donadores

def R(m,Y,X,transx,transy):
    receptores = []
    #NORTE,SUR,ESTE,OESTE
    N = [(transy + Y)%(Y+1),transx]
    S = [(transy + 1)%(Y+1),transx]
    E = [transy,(transx + X)%(X+1)]
    O = [transy,(transx + 1)%(X+1)]
    
    if (m[N[0],N[1]] < m[transy,transx]):
        receptores.append(N) 
    if (m[S[0],S[1]] < m[transy,transx]):
        receptores.append(S) 
    if (m[E[0],E[1]] < m[transy,transx]):
        receptores.append(E)
    if (m[O[0],O[1]] < m[transy,transx]):
        receptores.append(O) 
        
    return receptores

def kernel1(m,d,Y,X,y,x,idy,idx,n1,n2):
    # PASO 1 matrix de donadores (lista de listas)
    for j in range(y):
        transy = idy*y + j ##transformaciones en y y en x
        for i in range(x):
            transx = idx*x + i
            donadores = D(m,Y,X,transx,transy)
            #print(donadores)
            if(len(donadores) == 0):
                coordenadas_d = [-1,-1]
            else:
                (r1,r2,n1,n2) = rand_std(n1,n2)
                coordenadas_d = donadores[round(r1*100)%len(donadores)]
            d[transy,transx,0] = coordenadas_d[0]
            d[transy,transx,1] = coordenadas_d[1]
    return(d,n1,n2)
    #<SYNCH>
def kernel2(m,d,r,Y,X,y,x,idy,idx,n1,n2):
    # PASO 1 matrix de donadores (lista de listas)
    for j in range(y):
        transy = idy*y + j ##transformaciones en y y en x
        for i in range(x):
            transx = idx*x + i
            receptores = R(m,Y,X,transx,transy)
            #print(receptores)
            if(len(receptores) == 0):
                coordenadas_r = [-1,-1]
            else:
                refinal = []
                for k in receptores:
                    if(d[k[0],k[1],0] == transy and d[k[0],k[1],1] == transx):
                        refinal.append(k)
                #print(refinal)
                if(len(refinal) == 0):
                    coordenadas_r = [-1,-1]
                else:
                    (r1,r2,n1,n2) = rand_std(n1,n2)
                    coordenadas_r = refinal[round(r1*100)%len(refinal)]
            #print(coordenadas_r)
            r[transy,transx,0] = coordenadas_r[0]
            r[transy,transx,1] = coordenadas_r[1]
    return(r,n1,n2)
    #<SYNCH>
def kernel3(m,r,Y,X,y,x,idy,idx,n1,n2):
    #PASO 3 y 4 MONEDA AL AIRE
    K = 1.386*10**-23 #J/K
    Tprom = 4 #KELVIN
    B = 10**-23 #factor de normalización J/radio_atómico
    for j in range(y):
        transy = idy*y + j ##transformaciones en y y
        for i in range(x):
            transx = idx*x + i
            (r1,r2,n1,n2) = rand_std(n1,n2)
            #calculo de Pm r1 < Pm
            ##dz = abs(m[j,i]-m[rmatrix[j][i][0],rmatrix[j][i][1]])
            ##Pm = 
            if(round(r1) < 0.5 and (r[transy,transx,0] != -1)):
                m[transy,transx] -= 1
                m[r[transy,transx,0],r[transy,transx,1]] += 1
    #PASO 5 Pa: mas o menos 1%
            (r1,r2,n1,n2) = rand_std(n1,n2)
            Pa = 1/(X*Y)
            #print(Pa)
            if(r1 < Pa):
                m[transy,transx] += 1
                
    return(m,n1,n2)


# In[7]:


from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.pyplot as plt
import numpy as np
cdict = {'red': ((0.0, 0.0, 0.0),   
                 (1.0, 0.0, 0.0)),

        'green': ((0.0, 0.0, 0.0),
                  (1.0, 1.0, 1.0)),

        'blue':  ((0.0, 1.0, 1.0),
                  (1.0, 0.0, 0.0))
        }
blue_green = LinearSegmentedColormap('GreenRed', cdict, N=10, gamma=1.0)
colores = ['red','green','blue']


# In[8]:


#SEMILLAS
n1 = -4545829
n2 = -7**5
#CALENTAMIENTO DE ALGORITMO  r1,r2: números pseudoaleatorios(0,1)
for i in range(100):
    (r1,r2,n1,n2) = rand_std(n1,n2)
    
#INICIALIZACIÓN DE MATRICES Y DIMENSIONES
X = 30 ## dimensiones Y,X totales de la matriz
Y = 30
x = 5 ## dimensiones y,x por cada bloque
y = 5
m = np.zeros((Y,X),dtype=int) ##matriz principal [y,x] | m[0] = y | m[1] = x 
d = np.zeros((Y,X,2),dtype=int)##matriz de donadores d[posy][posx][y,x del donador]
r = np.zeros((Y,X,2),dtype=int)##matriz de receptores d[posy][posx][y,x del receptor]
iteraciones = 1000
matriz_salida = []
#print(r)
#matriz_salida.append(m.copy())
for n in range(iteraciones):
    if (n%(iteraciones/10) == 0):
        matriz_salida.append(m.copy())

    for j in range(int(X/x)):
        idy = j
        for i in range(int(Y/y)):
            idx = i
            #print([j,i])
            d,n1,n2 = kernel1(m,d,Y-1,X-1,y,x,idy,idx,n1,n2)
            r,n1,n2 = kernel2(m,d,r,Y-1,X-1,y,x,idy,idx,n1,n2)
            m,n1,n2 = kernel3(m,r,Y-1,X-1,y,x,idy,idx,n1,n2)
    
            #m = kernel(m,idt,x,y,n1,n2)


# In[9]:


import matplotlib.pyplot as plt
   
#print(matriz_salida[0])
#print("\n\n\n")
#print(matriz_salida[1])

for n in range(len(matriz_salida)):
    # set up the figure and axes
    fig = plt.figure(figsize=(7,7))
    ax1 = fig.add_subplot(111, projection='3d')
    # fake data
    _x = np.arange(X)
    _y = np.arange(Y)
    width = depth = 1
    ax1.set_zlim(0, X)
    #_xx, _yy = np.meshgrid(_x, _y)
    #x, y = _xx.ravel(), _yy.ravel()
    for j in range(Y):
        for i in range(X):
            top = matriz_salida[n][j,i]
            bottom = np.zeros_like(top)
            ax1.bar3d(i, j, bottom, width, depth, top, shade=True,color = blue_green(top*3))
    ax1.set_title('Iteración'+ str(n*iteraciones/10))
    #print('Iteración'+ str(n+1))
    plt.pause(0.000000001)

# In[ ]:




