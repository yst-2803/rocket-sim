import numpy as np
import matplotlib.pyplot as plt
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

y = np.array([0.0,0.0,m0])
t = 0.0
ts,hs,vs = [t], [y[0]], [y[1]]

while t < 10:
    y = y + derivs(t,y)*dt
    t = t + dt
    ts.append(t)
    hs.append(y[0])
    vs.append(y[1])

print("First time: ",t_data[0])
print("Number of points: ",len(t_data))
print("Total impulse: ",I_total)
print("Apogee:", max(hs))

#Plotting the data points

fig, (ax1, ax2) = plt.subplots(2,1,sharex = True)
ax1.plot(ts,hs)
ax1.set_ylabel("Altitude (m)")
ax2.set_xlabel("Time (s)")
ax2.plot(ts,vs)
ax2.set_ylabel("Velocity (m/s)")
plt.show()