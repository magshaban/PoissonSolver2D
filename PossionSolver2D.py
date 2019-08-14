# This is the main file.
# this code can be used to 
# Code by: Maged Shaban


import numpy as np 
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from math import sqrt 
 
from meshgen import *
from poltall import *

p,t,e = meshgen()

print(e.shape)

def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))


def gradients(x,y):
    area = PolyArea(x,y)
    b = np.array(([y[1]-y[2]],[y[2]-y[0]],[y[0]-y[1]]))/2/area
    c = np.array(([x[2]-x[1]],[x[0]-x[2]],[x[1]-x[0]]))/2/area
    return area,b,c


def a(x,y):
    return 1

# This function to evaluate the Stiffness Matrix a(u,v)
def stiffMat2D(p,t):
    
    np1 = p.shape[1]    # number of nodes
    nt = t.shape[1]     # number of elements 
    
    A = np.zeros((np1,np1)) # allocate the Mass Matrix 
    
    for k in range(nt):
        loc2glob = t[0:3,k]-1
        x = p[0,loc2glob]
        y = p[1,loc2glob]
        area,b,c = gradients(x,y)
        xc = np.mean(x)
        yc = np.mean(y)
        abar = a(xc,yc)
        
        AK = abar*(b@b.T + c@c.T)*area
        #assmble the local matrix to the global 
        A[np.ix_(loc2glob,loc2glob)] += AK   
         
    return A 

A =  stiffMat2D(p,t)
 

# to plot the stifness matrix 
fig, (ax1, ax2) = plt.subplots(1,2)
fig.suptitle('The Mass Matrix')
ax1.matshow(A)
ax2.spy(A)
plt.show()


# The given function
def func(x,y):
    return x+y**2


# This laod vector f(v) evaluation.
def LoadVec2D(p,t):

    np1 = p.shape[1]
    nt = t.shape[1]

    b = np.zeros((np1,1))

    for k in range(nt):
        loc2glob = t[0:3,k]-1
        x = p[0,loc2glob]
        y = p[1,loc2glob]     
        area = PolyArea(x,y)        
        bk =np.array(([func(x[0],y[0])],
                      [func(x[1],y[1])],
                      [func(x[2],y[2])])) /3*area
        b[loc2glob] += bk

    return b


def kappa(x,y):
    return 1


def gN(x,y):
    return x


def gD(x,y):
    return 1

    
# The following codes for the boundary condition. 
def RobinMat2D(p,e): 
    np1 = p.shape[1]
    ne = e.shape[1]
    
    R = np.zeros((np1,np1))
    
    for k in range(ne):
        loc2glob = t[0:2,k]-1
        x = p[0,loc2glob]
        y = p[1,loc2glob]
        length = sqrt((x[0]-x[1])**2+(y[0]-y[1])**2)
        xc = np.mean(x)
        yc = np.mean(y)
        k = kappa(xc,yc) # to be defined later
        RK = k/6*np.array([[2,1],[1,2]])*length
        R[np.ix_(loc2glob,loc2glob)] += RK
    
    return R


def RobinVec2D(p,e):
    np1 = p.shape[1]
    ne = e.shape[1]
    
    r = np.zeros((np1,1))   
    
    for k in range(ne):
        loc2glob = t[0:2,k]-1
        x = p[0,loc2glob]
        y = p[1,loc2glob]
        length = sqrt((x[0]-x[1])**2+(y[0]-y[1])**2)
        xc = np.mean(x)
        yc = np.mean(y)
        tmp = kappa(xc,yc)*gD(xc,yc) + gN(xc,yc) # to be defined later
        rk = tmp*np.array([[1],[1]])*length/2
        r[loc2glob] += rk

    return r


# (A+R)u = (b+r)
A = stiffMat2D(p,t)
R = RobinMat2D(p,e)
b = LoadVec2D(p,t)
r = RobinVec2D(p,e)

A_R = A + R
b_r = b + r

u = np.linalg.inv(A_R)@b_r
 
print(u)


##M = massmat(p,t)
#b = LoadVec2D(p,t)
#
#pf =  A@b
##
##
#print('The Mass Matrix(M) = \n',M)
#print('\nThe load vector (b)= \n',b)
#print('\nThe L2 Projection (pf) of the given function f(x) is: \n',pf)
#print('\n=======================================\n')
#print('The total number of nodes = ',p.shape[1])
#print('The total number of elements = ', t.shape[1])
#print('\n=======================================\n')
#
#delunay_tri(A)
poltall(x = p[0,:],y = p[1,:], z = u)
