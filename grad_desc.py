import numpy as np
import matplotlib.pyplot as plt
import subprocess

flexcode = """
TITLE 'DA2'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
Temp(threshold=.01) !temperature in C
SELECT         { method controls }
ngrid = 19
DEFINITIONS    { parameter definitions }
!Variables Params
T_air = 120	!Degrees celcius
h_convection = 200 !The convection heat transfer coefficient
P_mw = %s !The power of the microwave W
P_skillet = %s !Power of Skillet

 
! vairable constants 
rho
k
cp
epsilon_m
T_ideal
Init_Temp
! heat equation stuff
qvol
qdot = -k*grad(Temp)
 
! total volume weighted absorbidity 
absorb_tot =  0.4*area_integral(epsilon_m,'crust') + 0.39*(area_integral(epsilon_m, "filling"))
! calculating total mass
rho_crust = 700
rho_filling = 1200
area_filling = 1/2*pi*0.035^2 !m^3
area_crust = (0.005*0.08+1/2*pi*0.04^2-1/2*pi*0.035^2) 
massf = rho_filling*area_filling*.39 !filling mass
massc = rho_crust*area_crust*.4 !crust mass
mass = massc+massf
! setting up tastiness
Ttol = 4
T_integral=.4*area_integral(rho*((Temp-T_ideal)/ttol)^4, 'crust')+0.39*area_integral(rho*((Temp-T_ideal)/ttol)^4, 'filling')
Tastiness=1/(1+T_integral/mass)
!Skillet values
thickness_skillet = .005 !m thick
width_skillet = 9 * 0.01 !m wide
length_skillet = 41 * 0.01 !m long
! other geometry values
radius_crust = 0.04 !radius in meters
radius_filling = 0.035

INITIAL VALUES
	Temp = Init_temp
EQUATIONS        
! heat equation
dt(rho*cp*Temp) = div(k*grad(Temp)) + qvol

! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  Region 'skillet'
 
  ! contants for skillet
  k = 10
  rho = 3500
  cp = 1300
  epsilon_m = 0
  T_Ideal = 0
  Init_Temp = 24
 
  ! volumetric heat generation over volume so it's per m^3  
  qvol =  P_skillet/(integral(1,"skillet")*0.41)
 
  START (-width_skillet/2, thickness_skillet)
    	line to (width_skillet/2, thickness_skillet)  load(temp) = 0
        line to (width_skillet/2, 0)
        line to (-width_skillet/2, 0)
        line to close
 

  REGION 'crust'       
  !Crust values
  	k = 0.5 !W/m^2-K
	rho = rho_crust !kg/m^3
	cp = 2500 !J/kg-K
	epsilon_m = 0.05 !5
	T_ideal= 75 !Degrees celcius
    Init_Temp = 4
    ! microwave generation
    qvol =  (P_mw*epsilon_m)/absorb_tot
    START(radius_crust,thickness_skillet) load(temp) = h_convection*(temp-t_air)  
		line to (radius_crust, thickness_skillet+0.005) 
        arc(center=0,thickness_skillet+0.005) angle 180 
    	line to (-radius_crust, thickness_skillet) load(temp) = 0 
        line to (radius_crust,thickness_skillet) to close

    REGION 'filling'
        k = 1 !W/m^2-K
		rho = rho_filling !kg/m^3
		epsilon_m = 0.4 ! 40
        cp = 4200 ! J/kg-K
		T_ideal = 75 !Degrees celcius
        Init_Temp = 4
        qvol =  (P_mw*epsilon_m)/absorb_tot
		start(radius_filling,0.005+thickness_skillet) !load(temp) = 0 
			arc(center=0,0.005+thickness_skillet) angle 180  
			line to (0,0.005+thickness_skillet) to close
 
 
TIME 0 TO  60{ if time dependent }
PLOTS            { save result displays }
for t = endtime
  CONTOUR (temp) painted
  vector(qdot) norm
SUMMARY
EXPORT FILE 'ELOut.txt'
report(tastiness*100) !accounts for percent
END
"""

max_tastiness = 0
max_P_mw = 0
max_P_skillet = 0

flexfilename = "assignment_2_heatflow.pde"

def run_code(r):
    global max_tastiness
    global max_P_mw
    global max_P_skillet
    print(r)

    P_mw = 0
    P_skillet = 0

    P_mw = r[0]
    P_skillet = r[1]

    with open(flexfilename, "w") as f:
        print(flexcode%(P_skillet, P_mw), file=f)
    
    subprocess.run(["C:\FlexPDE6student\FlexPDE6s.exe", "-S", flexfilename], timeout=5)

    tastiness = 0

    with open("ELOut_1.txt", 'r') as f:
        for i, line in enumerate(f):
            if i == 7:  # Since line counting starts at 0
                parts = line.split('=')
                if len(parts) > 1:
                    tastiness = float(parts[1].strip())

    if tastiness > max_tastiness:
        max_tastiness = tastiness
        max_P_mw = P_mw
        max_P_skillet = P_skillet

    return tastiness


'''
def f(r):
    x = r[0]
    y = r[1]
    return(1./((x-3.147128)**2+(y-2.73)**2+1)+.1*x+0.01*np.cos(x*10)+0.01*np.sin(y*10))
'''

def grad_f(r, delta=10, f=run_code):
    x = r[0]
    y = r[1]
    grad = np.zeros(2)
    grad[0] = (f([x+delta, y])-f([x-delta, y]))/(2*delta)
    grad[1] = (f([x, y+delta])-f([x, y-delta]))/(2*delta)

    return grad*100

def grad_desc(r, alpha=1000, tol=0.000001, max_iter=1000):
    for i in range(max_iter):
        gradient = grad_f(r)
        r = r + alpha*gradient
        if np.linalg.norm(grad_f(r)) < tol:
            print(np.linalg.norm(grad_f(r)))
            break
    return r


r = np.array([200, 3000])
r = grad_desc(r)
print(r)
tastiness_value = run_code(r)

#This accounts for weird issues in flex that lead to value fringing around points, 
#this helps to fix issues with discountinuity which gradient ascent can't handle nicely
#The gradient ascent is used to get close to the optimal value, then we use this to get the final value

if max_tastiness > tastiness_value:
    tastiness_value = max_tastiness
    r[0] = max_P_mw
    r[1] = max_P_skillet

print(tastiness_value, " this is the final tastiness value.")
print(r, " this is the final power values.")

