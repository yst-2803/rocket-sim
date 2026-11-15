#Works backwards from OpenRocket's trajectory to the drag coefficient it must
#have been using, so it can be compared with the constant Cd = 0.56 in my sim.
#
#From the equation of motion,   m*a = T - m*g - 0.5*rho*Cd*A*v^2
#so                             Cd = (T - m*g - m*a) / (0.5*rho*A*v^2)
#
#OpenRocket's acceleration a is found from the slope of its velocity data.
#Assumes OpenRocket used the same masses and thrust curve as rocket_sim.py.

import numpy as np

g = 9.81
A = np.pi * 0.0125**2
m0 = 0.0571
m_prop = 0.0108

#Thrust curve and the impulse delivered up to each data point

data = np.loadtxt("Estes_C6.csv",delimiter=",",skiprows=5)
t_data = np.insert(data[:,0],0,0)
T_data = np.insert(data[:,1],0,0)
impulse = np.concatenate(([0.0], np.cumsum(0.5*(T_data[1:]+T_data[:-1])*np.diff(t_data))))
I_total = impulse[-1]
t_burn = t_data[-1]

#OpenRocket export

or_data = np.loadtxt("openrocket_C6-7.csv",delimiter=",",comments="#")
or_t = or_data[:,0]
or_h = or_data[:,1]
or_v = or_data[:,2]

#ISA air density, which is what OpenRocket uses

def density(h):
    return 1.225 * (1 - 2.25577e-5*h)**4.25588

#Same thrust and mass model as rocket_sim.py

def thrust(t):
    if t >= t_burn:
        return 0.0
    return np.interp(t,t_data,T_data)

def mass(t):
    return m0 - m_prop*np.interp(t,t_data,impulse)/I_total

#Fit a quadratic to OpenRocket's velocity in a short window around time t.
#Its gradient is the acceleration, which is less noisy than differencing
#neighbouring points.

def effective_cd(t,window):
    near = abs(or_t - t) < window
    fit = np.polyfit(or_t[near],or_v[near],2)
    v = np.polyval(fit,t)
    a = np.polyval(np.polyder(fit),t)
    h = np.interp(t,or_t,or_h)
    m = mass(t)
    return v, (thrust(t) - m*g - m*a) / (0.5*density(h)*A*v**2)

print("Phase   t (s)   v (m/s)   effective Cd")
for t in [0.9,1.1,1.4,1.6]:
    v,cd = effective_cd(t,0.1)
    print("burn   ",t,"  ",round(v,1),"   ",round(cd,3))
for t in [2.0,2.5,3.0,4.0,5.0,6.0,7.0]:
    v,cd = effective_cd(t,0.15)
    print("coast  ",t,"  ",round(v,1),"   ",round(cd,3))
