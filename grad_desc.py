import numpy as np
import matplotlib.pyplot as plt
import subprocess

flex_code = """"
TITLE 'New Problem'     { the problem identification }

COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }

VARIABLES        { system variables }

Temp(threshold=.01) !temperature in degC
SELECT         { method controls }
ngrid = 5
DEFINITIONS    { parameter definitions }

!Variables Params
T_air = 120 !Degrees celcius
h_convection = 200 !The convection heat transfer coefficient
P = 500 !The power of the microwave W

rho
k
cp
epsilon_m
T_ideal
Init_Temp

!Skillet values
k_skillet = 10 !W/m^2-K
rho_skillet = 3500 !kg/m^3
cp_skillet = 1300 !J/kg-K
thickness_skillet = .5*0.01 !m thick
width_skillet = 9 * 0.01 !m wide
length_skillet = 41 * 0.01 !m long

qdotvol_skillet = k_skillet * rho_skillet * cp_skillet !Volumetric Heating in the Skillet
qdotvol_microwave = 0


radius_crust = 0.04 !radius of the crust in meters
radius_filling = 0.03


INITIAL VALUES
	Temp = 4
    
EQUATIONS        { PDE's, one for each variable }

dt(rho*cp*Temp) = div(k_skillet*grad(Temp)) + qdotvol_skillet + h_convection*(T_air - Temp)


! CONSTRAINTS    { Integral constraints }

BOUNDARIES       { The domain definition }

  REGION 'crust'       { For each material region }
  !Crust values
  	k = 0.5 !W/m^2-K
	rho = 700 !kg/m^3
	cp = 2500 !J/kg-K
	epsilon_m = 0.05 !Percent * 0.01
	T_ideal= 75 !Degrees celcius
    START(radius_crust,0)   { Walk the domain boundary }
		arc(center=0,0) angle 180 load(temp) = 0
		line to (0,0) to close

	
    REGION 'filling'
        k = 1 !W/m^2-K
		rho = 1200 !kg/m^3
		epsilon_m = 0.4 ! Percent + 0.01
        cp = 4200 ! J/kg-K
		T_ideal = 75 !Degrees celcius
		start(radius_filling,0.005) !value(temp)=0
			arc(center=0,0.005) angle 180 load(temp) = 0 
			line to (0,0.005) to close

TIME 0 TO 5   { if time dependent }

PLOTS            { save result displays }
for t = 0 by endtime/5 to endtime
  CONTOUR(Temp) painted
!	vector(qdot) norm
	summary
END
"""

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
