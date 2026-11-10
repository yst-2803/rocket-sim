import numpy as np

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

print("First time: ",t_data[0])
print("Number of points: ",len(t_data))
print("Total impulse: ",I_total)