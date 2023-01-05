# -*- coding: utf-8 -*-
"""
Disclamer
---------
This is not an original work. The code bellow is just a translation 
from matlab to Python of the original QPmR algorithm, for further details 
visit Prof. Tomas Vyhlidal site following the link bellow
http://www.cak.fs.cvut.cz/vyhlidal


References
----------
[1] Vyhlídal, T., & Zítek, P. (2003). 
Quasipolynomial mapping based rootfinder for analysis of time delay systems. 
IFAC Proceedings Volumes, 36(19), 227-232.

Example
-------
# s from Laplace domain, needed to definde the quasi-polynomial
s = sp.symbols("s")

# Some parameters
kp = 1
kd = 30
tau1 = 0.01

# Definition of the quasi-polynomial. This is literraly the quasi-polynomial of which you want to find solutions
QP=s**3+(2/5)*s**2-(13/20)*s+1/20+(kp+kd*(1-exp(-tau1*s))*(3-exp(-s*tau1))/(2*tau1))*(1/20-(1/100)*s**2)

# This is the region of interest on the complex given as: [minRealValue maxRealValue minImagValue maxImagValue] 
Region = [-200, 5000, -10, 10]

# The output of the QPmR function is an array containing all the roots found in the given region
r1 = QPmR(QP,np.array(Region),0.1*np.pi/2.5,0.000000001);

#Several examples were performed in order to verify that the results given by this function match those of the matlab function of Prof. Tomas Vyhlidal 
"""

from lib2to3.pygram import Symbols
import numpy as np
import sympy as sp
from sympy import exp
from matplotlib.pyplot import contour

import matplotlib.pyplot as plt
def QPmR(M,Reg,d,e):
    b_min = Reg[0] - 3*d
    b_max = Reg[1] + 3*d
    w_min = Reg[2] - 3*d
    w_max = Reg[3] + 3*d
    
    """
    m = M.shape
    
    if m[0]>1:
        M = poly(M,s)
    """

    P='No roots in the selected region'

    nW = round((w_max-w_min)/d + 1)
    nB = round((b_max-b_min)/d + 1)
    W = np.linspace((w_min), (w_max), nW)
    B = np.linspace((b_min), (b_max), nB)
    W = np.reshape(W, (1, nW))
    B = np.reshape(B, (1, nB))
    oB = np.ones((1, nW))
    oW = np.ones((1, nB))

    WW = np.matmul(np.transpose(oW), W)
    BB = np.matmul(np.transpose(oB), B)
    D = BB + np.transpose(WW)*1j

    s = sp.MatrixSymbol("s", D.shape[0], D.shape[1])
    func = sp.utilities.lambdify(s, M,'numpy') # returns a numpy-ready function
    MM = func(D)
    #MM = M.subs({s: sp.Matrix(D)})

    #MM = M.double(MM)

    R = MM.real
    I = MM.imag
 
    Cr = contour(B[0], W[0],R,[0])
    Cr = Cr.allsegs

    aux = np.empty((0, 2))
    for seg in Cr[0]:
        aux = np.append(aux, np.array(seg), axis=0)
    Cr = np.transpose(aux)

    pCr=Cr.shape
    pCr=pCr[1]
    Crr=Cr
    Cr=Cr[0,:] + Cr[1,:]*1j

    IIc=np.imag(func(Cr))

    Np=1
    ch=1
    Pp = np.empty((0))
    n = -1
    for no in range(2, pCr-1):
        if((IIc[no-2]*IIc[no])<0) and (ch==1):
            if (Reg[0] <= np.real(Cr[no-1])) and (Reg[1] >= np.real(Cr[no-1])) and (Reg[2] <= np.imag(Cr[no-1])) and (Reg[3] >= np.imag(Cr[no-1])):
                if (abs(Crr[0][no-1]-Crr[0][no]) < 2*d) and (abs(Crr[1][no-1]-Crr[1][no]) < 2*d):
                    Pp = np.append(Pp, np.array(Cr[no - 1]))
                    Np = Np + 1
                    ch = 0
        else:
            ch=1
    
    P=np.empty((0))

    #numeric calculation of poles
    if not (isinstance(P, str) and len(P) == 1):
        H = M
        s = sp.symbols("s")
        dH = sp.diff(M, s)
        dfunc = sp.utilities.lambdify(s, dH,'numpy') # returns a numpy-ready function
        n = Pp.shape
        n = n[0]
        for o in range(1,n+1):
            a = Pp[o-1]
            
            Newton_run = 1
            while Newton_run:
                a_k1 = a
                #s = a
                a = a - func(a) / dfunc(a)
                
                Newton_run = (abs(a - a_k1) > e)

            if (not np.isreal(a)) and (abs(np.imag(a)) < 0.1*e):
                a = np.real(a)
            
            P = np.append(P, a)

    Pmax = max(abs(P))

    PM = np.empty((0))

    for k in range(n):
        o = np.where(np.abs(P) == np.amin(np.abs(P)))
        o = o[0][0]
        PM = np.append(PM, P[o])
        P[o]=2*Pmax
    for elem in PM:
        print(elem)
    return PM
