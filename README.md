# QPmR
This is a python version of the QPmR algorithm developed in Matlab by Prof. Tomas Vyhlidal

# Update
A Python version of the most recent Matlab implementation of the QPmR with pip installation available can be find in: [link](https://github.com/LockeErasmus/qpmr)

<a href="https://github.com/LockeErasmus/qpmr">link</a>

## Disclamer
---------

This is not an original work. The code in qpmr.py is just a translation from matlab to Python of the original QPmR algorithm, for further details visit Prof. Tomas Vyhlidal site following the link bellow:

http://www.cak.fs.cvut.cz/vyhlidal

## Acknowledgments
----------
The translation was performed with the help of E.E. Arturo Tapia

Github: https://github.com/lohug1

## References
----------

[1] Vyhlídal, T., & Zítek, P. (2003). Quasipolynomial mapping based rootfinder for analysis of time delay systems. IFAC Proceedings Volumes, 36(19), 227-232.

## Example
---------
```
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

#Several examples were performed in order to verify that the results given by this function matched those of the matlab function of Prof. Tomas Vyhlidal 
```
