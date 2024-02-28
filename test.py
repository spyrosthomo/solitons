import numpy as np 
import matplotlib.pyplot as plt 

def shoot_nonlinear(a,b,alpha, beta, n, tol, M):

    w1 = np.zeros(n+1)  
    w2 = np.zeros(n+1)
    h = (b-a)/n
    k = 1
    TK = (beta - alpha)/(b - a)

    print("i""  x" "     " "W1""     " "W2")
    while k <= M:

        w1[0] = alpha
        w2[0] = TK
        u1    = 0
        u2    = 1

        for i in range(1,n+1):
            x = a + (i-1)*h     #step 5

            t = x + 0.5*(h)

            k11 = h*w2[i-1]     #step 6

            k12 = h*f(x,w1[i-1],w2[i-1])
            k21 = h*(w2[i-1] + (1/2)*k12)
            k22 = h*f(t, w1[i-1] + (1/2)*k11, w2[i-1] + (1/2)*k12)
            k31 = h*(w2[i-1] + (1/2)*k22)
            k32 = h*f(t, w1[i-1] + (1/2)*k21, w2[i-1] + (1/2)*k22)
            t   = x + h
            k41 = h*(w2[i-1]+k32)
            k42 = h*f(t, w1[i-1] + k31, w2[i-1] + k32)
            w1[i] = w1[i-1] + (k11 + 2*k21 + 2*k31 + k41)/6
            w2[i] = w2[i-1] + (k12 + 2*k22 + 2*k32 + k42)/6   
            kp11 = h*u2
            kp12 = h*(fy(x,w1[i-1],w2[i-1])*u1 + fyp(x,w1[i-1], w2[i-1])*u2)
            t    = x + 0.5*(h)
            kp21 = h*(u2 + (1/2)*kp12)
            kp22 = h*((fy(t, w1[i-1],w2[i-1])*(u1 + (1/2)*kp11)) + fyp(x+h/2, w1[i-1],w2[i-1])*(u2 +(1/2)*kp12))
            kp31 = h*(u2 + (1/2)*kp22)
            kp32 = h*((fy(t, w1[i-1],w2[i-1])*(u1 + (1/2)*kp21)) + fyp(x+h/2, w1[i-1],w2[i-1])*(u2 +(1/2)*kp22))
            t    = x + h
            kp41 = h*(u2 + kp32)
            kp42 = h*(fy(t, w1[i-1], w2[i-1])*(u1+kp31) + fyp(x + h, w1[i-1], w2[i-1])*(u2 + kp32))
            u1   = u1 + (1/6)*(kp11 + 2*kp21 + 2*kp31 + kp41)
            u2   = u2 + (1/6)*(kp12 + 2*kp22 + 2*kp32 + kp42)


        r = abs(w1[n] - beta)
        #print(r)
        if r < tol:
            for i in range(0,n+1):
                x = a + i*h
                #print("%.2f %.2f %.4f %.4f" %(i,x,w1[i],w2[i]))
            return w1

        TK = TK -(w1[n]-beta)/u1

        k = k+1


    print("Maximum number of iterations exceeded")   
    return w1



def f(x,y,yp):
    #fx = (1/8)*(32 + 2*x**3 -y*yp)
    mu    = 0.3
    omega = 0.5
    fx    =  (1-omega**2)*y - 2*y**3 + 3*mu*y**5 
    return fx

def fy(xp,z,zp):
    #fyy = -(1/8)*(zp)
    mu    = 0.3
    omega = 0.5
    fyy = (1-omega**2) - 6*z**2 + 15*mu*z**4
    return fyy

def fyp(xpp,zpp,zppp):
    #fypp = -(1/8)*(zpp)
    fypp = 0
    return fypp

a = 1        # start point
b = 10       # end point
alpha = 0.5    # boundary condition
beta = 0   # boundary condition
N = 2000        # number of subintervals
M = 10000        # maximum number of iterations
tol = 0.00001 # tolerance


x  = np.zeros(N+1)
h  = (b-a)/N
for i in range(0,N+1):
    x[i] = a + (i)*h

w1 = shoot_nonlinear(a,b,alpha,beta,N,tol,M)


# plot results 
fig = plt.figure()
ax  = fig.add_subplot(1,1,1)
ax.plot(x, w1)
plt.show()