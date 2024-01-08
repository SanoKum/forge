import numpy as np
import matplotlib.pyplot as plt

#with open("log", "r") as f:
#	for line in f:
#		if "ros" in line:
#			print(line.rstrip("\n"))

with open("log", "r") as f:
	lines = f.read().splitlines()

nLoop = 3

time=[]
ros =[]
rok =[]

for line in lines:
	if "Total Time" in line:
		time.append(float(line.split(" ")[-1]))

	if "ros" in line:
		ros.append(float(line.split(" ")[-1]))

	if "rok" in line:
		rok.append(float(line.split(" ")[-1]))

ros2 = np.array(ros)
ros3 = ros2[nLoop-1::nLoop]
ros4 = (ros3-ros3[0])/ros3[0]

rok2 = np.array(rok)
rok3 = rok2[nLoop-1::nLoop]
rok4 = (rok3)/rok3[0]

# plot

fig = plt.figure(figsize=(12,6))

ax1 = fig.add_subplot(122)
ax1.plot(time[0:len(ros4)], ros4)
ax1.set_xlabel("Time [s]")
ax1.set_ylabel("(del ros/init ros) [-]")
ax1.set_xlim([0,200])
ax1.set_ylim([-0.02,0.005])

ax2 = fig.add_subplot(121)
ax2.plot(time[0:len(rok4)], rok4)
ax2.set_xlabel("Time [s]")
ax2.set_ylabel("(rok/init rok) [-]")
ax2.set_xlim([0,200])
ax2.set_ylim([0.5,2.5])

plt.show()



