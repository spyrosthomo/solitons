import numpy as np
import math
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

#-------------------------------------------------
## BVP
N=50
a=1
b=2
h=(b-a)/N
x=np.linspace(a, b,N+1)
fig1 = plt.figure(figsize=(10,4))
plt.plot(x,0*x,'o:',color='red')
plt.xlim((a,b))
plt.xlabel('x',fontsize=16)
plt.title('Illustration of discrete time points for h=%s'%(h),fontsize=32)
plt.show()

#-------------------------------------------------
U1=np.zeros(N+1)
U2=np.zeros(N+1)
Z1=np.zeros(N+1)
Z2=np.zeros(N+1)

beta=43/3
alpha=17
lambda_app=[-16.53]
U2[0]=alpha##-2.5           # swap 1<->2
U1[0]=lambda_app[0]
                     # swap 1<->2
Z2[0]=0
Z1[0]=1
#-------------------------------------------------
tol=0.0001
k=0
fig2 = plt.figure(figsize=(10,8))
gs   = fig2.add_gridspec(2,1)
ax1  = fig2.add_subplot(gs[0,0])
ax1.set_title('U1(x)')
ax1.plot(x,beta*np.ones(N+1))
ax2  = fig2.add_subplot(gs[1,0])
ax2.set_title('U2(x)')
ax2.plot(x,alpha*np.ones(N+1))
ax1.grid(True); ax2.grid(True)
while  k<20:
    k=k+1
    for i in range (0,N):
        U1[i+1]=U1[i]+h*(U2[i])
        U2[i+1]=U2[i]+h*(32+2*x[i]**3-U1[i]*U2[i])/8#(-2*U2[i]*U1[i])
    
        Z1[i+1]=Z1[i]+h*(Z2[i])
        Z2[i+1]=Z2[i]+h*(-U2[i]*Z1[i]+U1[i]*Z2[i])/8#(-2*U2[i]*Z1[i]-2*Z2[i]*U1[i])

    lambda_app.append(lambda_app[k-1]-(U1[N]-beta)/Z1[N])
    
    ax1.plot(x,U1,':o',label=r"$\lambda$ %s"%(k))
    #ax1.xlabel('x',fontsize=16)
    #ax1.ylabel('U1',fontsize=16)
    
    ax2.plot(x,U2,':o',label=r"$\lambda$ %s"%(k))
    #ax2.xlabel('x', fontsize=16)
    #ax2.ylabel('U2', fontsize=16)

    U1[0]=lambda_app[k]             # swap 1<->2    
    if abs(lambda_app[k]-lambda_app[k-1])<tol:
        break
       
plt.legend(loc='best')
#plt.title(r"Iterations of $\lambda_k$",fontsize=32)
#plt.tight_layout()
plt.show()
print(lambda_app)
#-------------------------------------------------/
fig3 = plt.figure(figsize=(10,8))
gs   = fig3.add_gridspec(2,1)
ax1  = fig3.add_subplot(gs[0,0])
ax1.set_title('U1(x)')
ax1.plot(x,beta*np.ones(N+1))
ax2  = fig3.add_subplot(gs[1,0])
ax2.set_title('U2(x)')
ax2.plot(x,alpha*np.ones(N+1))

ax1.grid(True)
ax1.plot(x,U1,'b:o')
ax2.grid(True)
ax2.plot(x,U2,'b:o')

#plt.title("Numerical Approximation of the non-linear Boundary Value Problem",fontsize=32)
ax1.set_title("U1(x)")
ax2.set_title("U2(x)")#, fontsize=16)
plt.show()
#-------------------------------------------------
fig = plt.figure(figsize=(10,4))
plt.grid(True)
plt.plot(lambda_app,'o')
#plt.title("Values of $\lambda$ for each interation ",fontsize=32)
plt.xlabel('Iterations (k)',fontsize=16)
plt.ylabel("lambda",fontsize=16)
plt.show()