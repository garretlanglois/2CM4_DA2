import numpy as np
import matplotlib.pyplot as plt
import subprocess


'''

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

flexfilename = "assignment_2_heatflow.pde"

def run_code():
    with open(flexfilename, 'w') as f:
        print(flex_code)
    
    subprocess.run(["C:\FkexPDE"])

'''

def f(r):
    x = r[0]
    y = r[1]
    z = r[2]
    return(1./((x-3.147128)**2+(y-2.73)**2+1)+.1*x+0.01*np.cos(x*10)+0.01*np.sin(y*10))

def grad_f(r, delta=0.01, f=f):
    x = r[0]
    y = r[1]
    z = r[2]
    grad = np.zeros(3)
    grad[0] = (f([x+delta, y, z])-f([x-delta, y, z]))/(2*delta)
    grad[1] = (f([x, y+delta, z])-f([x, y-delta, z]))/(2*delta)
    grad[2] = (f([x, y, z+delta])-f([x, y, z-delta]))/(2*delta)

    return grad

def grad_desc(r, alpha=0.01, tol=1e-6, max_iter=10000):
    for i in range(max_iter):
        r = r + alpha*grad_f(r)
        if np.linalg.norm(grad_f(r)) < tol:
            break
    return r

r = np.array([0, 0, 0])
r = grad_desc(r)
print(r)

