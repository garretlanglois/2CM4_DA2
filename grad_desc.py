import numpy as np
import matplotlib.pyplot as plt
import subprocess

def run_code():


def f(r):
    x=r[0]
    y=r[1]
    return(1./((x-3.147128)**2+(y-2.73)**2+1)+.1*x+0.01*np.cos(x*10)+0.01*np.sin(y*10))

def grad(r,f0,delta=1e-3,f=f):
    dx=[delta,0]
    dy=[0,delta]
    return(np.array([f(r+dx)-f0, f(r+dy)-f0])/delta)

def dotprod(v1, v2):
    '''given two vectors, return their dot product'''
    return v1[0]*v2[0]+v1[1]*v2[1]

def relative_direction(v1, v2):
    '''given two vectors, find dot(v1hat,v2hat)'''
    return dotprod(v1,v2)/np.sqrt(dotprod(v1,v1)*dotprod(v2,v2))

scale = 2
'''Get initial grad & step'''
r0 = np.array([0,0]) #initial guess of 0,0,0
f0 = f(r0)
plt.plot(r0[0], r0[1], '*')
print(f'{r0=} {f0=}')
current_gradient = grad(r0,f0)

rused=[]
rused.append([r0[0],r0[1]])
rtried=rused.copy()
max_cosine=.2
for i in range(20):
    #plt.plot(r0[0], r0[1], '*')
    reldir=-1
    while(reldir<max_cosine):
        r1 = r0+scale*current_gradient
        f1 = f(r1)
        rtried.append([r1[0],r1[1]])

        nextgrad = grad(r1,f1)
        
        print(f'{current_gradient=}')
        print(f'{nextgrad=}')
        reldir = relative_direction(current_gradient,nextgrad)
        
        if(reldir<max_cosine):
            r2 = r1+scale*nextgrad #just for plotting where it would've gone
            plt.plot([r1[0],r2[0]], [r1[1],r2[1]], 'r-')
        
        print(f'{i}:f({r1})={f1} -> {reldir=} {scale=}')
        scale=scale/2
    scale=scale*3 #undo last halving if we were OK        
    current_gradient=nextgrad
    r0,f0=r1,f1 #advance points
    rused.append([r0[0],r0[1]])
    
rusedarray=np.array(rused)
rtriedarray=np.array(rtried)

plt.plot(rtriedarray[:,0], rtriedarray[:,1], '-*')
plt.plot(rusedarray[:,0], rusedarray[:,1], '-*')
