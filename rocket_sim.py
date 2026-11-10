import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
# Load the csv file

data = np.loadtxt("Estes_C6.csv",delimiter=",",skiprows=5)

#Split into columns

t_data = data[:,0]
T_data = data[:,1]

#Add initial values (0,0)

t_data = np.insert(t_data,0,0)
T_data = np.insert(T_data,0,0)

#Calculate the trapezium rule for the data

I_total = np.trapezoid(T_data,t_data)
t_burn = t_data[-1]

#Constants

g = 9.81
rho = 1.225
Cd = 0.56
A = np.pi * 0.0125**2
k = 0.5 * rho * Cd * A
m0 = 0.0571
m_prop = 0.0108
dt = 0.01

#T(t)

def thrust(t):

    if t>= t_burn:
        return 0.0
    return np.interp(t,t_data, T_data)

#State model function

def derivs(t,y):
    h,v,m = y
    T = thrust(t)
    a = (T - m*g - k*v*abs(v))/m
    if h<=0 and a< 0:
        a = 0
    dh_dt = v
    dv_dt = a
    dm_dt = -(m_prop/I_total)*T
    return np.array([dh_dt,dv_dt,dm_dt])

#Variables for rk45 and euler to ref

y = np.array([0.0,0.0,m0])
t = 0.0
ts,hs,vs = [t], [y[0]], [y[1]]

#euler implementation

while t < 10:
    y = y + derivs(t,y)*dt
    t = t + dt
    ts.append(t)
    hs.append(y[0])
    vs.append(y[1])

#rk45 implementation

def apogee(t,y):
    if y[0] < 1.0:
        return 1.0
    return y[1]
apogee.terminal = True
apogee.direction = -1
sol = solve_ivp(derivs, t_span=(0,10),y0 =[0.0,0.0,m0],
events=apogee, max_step=0.01)
euler_apogee = max(hs)
rk_apogee = sol.y_events[0][0][0]

#Reading from the open rocket csv

or_data = np.loadtxt("openrocket_C6-7.csv",delimiter=",",comments="#")
or_t = or_data[:,0]
or_h = or_data[:,1]
or_v = or_data[:,2]
or_apogee = max(or_h)

#Prints and error

print("OpenRocket apogee:", or_apogee)
print("RK45 vs OpenRocket:", (rk_apogee - or_apogee) / or_apogee * 100, "%")
print("Euler vs OpenRocket:", (euler_apogee - or_apogee) / or_apogee * 100, "%")
print("RK45 apogee:", rk_apogee)
print("Euler apogee:",euler_apogee)
print("RK45 time of apogee:", sol.t_events[0][0])
print("Difference:", (euler_apogee - rk_apogee) / rk_apogee * 100, "%")


#Plotting the data points
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(ts, hs, label="My sim (Euler)")
ax1.plot(or_t, or_h, "--", label="OpenRocket")
ax1.set_ylabel("Altitude (m)")
ax1.legend()
ax2.plot(ts, vs, label="My sim (Euler)")
ax2.plot(or_t, or_v, "--", label="OpenRocket")
ax2.set_ylabel("Velocity (m/s)")
ax2.set_xlabel("Time (s)")
ax2.legend()
ax1.set_xlim(0, 8.5)
plt.savefig("validation_vs_openrocket.png", dpi=200, bbox_inches="tight")
plt.show()