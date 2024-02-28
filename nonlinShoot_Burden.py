import numpy as np 
import matplotlib.pyplot as plt 

def method(a, b, alpha, beta, Nx, tol, Niter):
    h  = (b-a)/Nx
    k  = 0
    tk = (beta-alpha)/(b-a)
    
    w1 = np.zeros(Nx+1)
    w2 = np.zeros(Nx+1)

    while (k<=Niter):
        w1[0] = alpha
        w2[0] = tk
        u1    = 0
        u2    = 1
        
        for i in range(1, Nx+1):
            x   = a + (i-1)*h
            k11 = h*w2[i-1]
            k12 = h*f(x      , w1[i-1]         , w2[i-1]          )
            k21 = h*(w2[i-1] + 0.5*k12)
            k22 = h*f(x+0.5*h, w1[i-1]+0.5*k11 , w2[i-1] + 0.5*k12)

            k31 = h*(w2[i-1] + 0.5*k22)
            k32 = h*f(x+0.5*h, w1[i-1] + 0.5*k21, w2[i-1] + 0.5*k22)

            k41 = h*(w2[i-1] + k32)
            k42 = h*f(x+h, w1[i-1]+k31, w2[i-1]+k32)

            w1[i] = w1[i-1] + (k11 + 2*k21 + 2*k31 + k41)/6
            w2[i] = w2[i-1] + (k12 + 2*k22 + 2*k32 + k42)/6

            kt11  = h*u2
            kt12  = h*(fy(x, w1[i-1], w2[i-1])*u1 + fyt(x, w1[i-1], w2[i-1])*u2)

            kt21  = h*(u2 + 0.5*kt12)
            kt22  = h*(fy(x+0.5*h, w1[i-1], w2[i-1])*(u1+0.5*kt11) + fyt(x+0.5*h, w1[i-1], w2[i-1]*(u2+0.5*kt12)))

            kt31  = h*(u2 + 0.5*kt22)
            kt32  = h*(fy(x+0.5*h, w1[i-1], w2[i-1])*(u1 + 0.5*kt21) + fyt(x+0.5*h, w1[i-1], w2[i-1])*(u2+0.5*kt22)) 

            kt41  = h*(u2+kt32)
            kt42  = h*(fy(x+h, w1[i-1], w2[i-1])*(u1+kt31) + fyt(x+h, w1[i-1], w2[i-1])*(u2+kt32))

            u1    = u1 + (1/6)*(kt11 + 2*kt21 + 2*kt31 + kt41)
            u2    = u2 + (1/6)*(kt12 + 2*kt22 + 2*kt32 + kt42)

            test  = np.abs(w1[-1]-beta) 


            if test < tol: 
                x   = np.zeros(Nx+1)
                for i in range(0, Nx+1):
                    x[i] = a + i*h
                sol = np.hstack((x.reshape(Nx+1,1), w1.reshape(Nx+1,1), w2.reshape(Nx+1,1)))
                return sol
            
            tk = tk - (w1[-1]- beta)/u1

            k  += k 

#-------------------------------------------

def f(x, y, yp): 
    mu    = 0.3
    omega = 0.5
    return (1-omega**2)*y - 2*y**3 + 3*mu*y**5 # RHS

#---
def fy(xp, z, zp):
    mu    = 0.3
    omega = 0.5
    return  (1-omega**2) - 6*z**2 + 15*mu*z**4 #d(RHS)/d(psi)

#---
def fyt(xpp, zpp ,zppp):
    return 0    #d(RHS)/d(psi')

#---
a =  0   # interval 
b = 10

alpha = 1
beta  = 0 # we should give ""beta""" inside the method because we want 
          #    y'(b) = beta  

Nx    = 20
Nmax  = 5

tol   = 1e-5


x      = np.zeros(Nx)
h      = (b-a)/Nx
for i in range(0, Nx):
    x[i] = a + (i)*h

sol     = method(a,b,alpha,beta,Nx,tol,Nmax)
print(sol)