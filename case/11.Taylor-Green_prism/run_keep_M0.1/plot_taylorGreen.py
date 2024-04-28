import numpy as np
import matplotlib.pyplot as plt

#with open("log", "r") as f:
#	for line in f:
#		if "ros" in line:
#			print(line.rstrip("\n"))

#4th
with open("log", "r") as f:
	lines = f.read().splitlines()

nLoop = 4

M0=0.4

time=[]
ros =[]
rok =[]

for line in lines:
	if "Total Time" in line:
		time.append(M0*float(line.split(" ")[-1]))

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
ax1.plot(time[0:len(ros4)], ros4, linewidth=1.5, color='g', label="4th")
ax1.set_xlabel("Time [s]")
ax1.set_ylabel("(del ros/init ros) [-]")
ax1.set_xlim([0,200])
ax1.set_xticks( np.arange(0, 201, 20))
ax1.set_ylim([-0.04,0.01])
ax1.hlines(0.0, 0, 200, color='black', linestyles='dotted')

ax2 = fig.add_subplot(121)
ax2.plot(time[0:len(rok4)], rok4, linewidth=1.5, color='g', label="4th")
ax2.set_xlabel("Time [s]")
ax2.set_ylabel("(rok/init rok) [-]")
ax2.set_xlim([0,200])
ax2.set_ylim([0.6,1.6])
ax2.set_xticks( np.arange(0, 201, 20))

ax2.hlines(1.0, 0, 200, color='black',linestyles='dotted')


##3rd
#with open("log_3rd", "r") as f:
#	lines = f.read().splitlines()
#
#nLoop = 3
#
#M0=0.4
#
#time=[]
#ros =[]
#rok =[]
#
#for line in lines:
#	if "Total Time" in line:
#		time.append(M0*float(line.split(" ")[-1]))
#
#	if "ros" in line:
#		ros.append(float(line.split(" ")[-1]))
#
#	if "rok" in line:
#		rok.append(float(line.split(" ")[-1]))
#
#ros2 = np.array(ros)
#ros3 = ros2[nLoop-1::nLoop]
#ros4 = (ros3-ros3[0])/ros3[0]
#
#rok2 = np.array(rok)
#rok3 = rok2[nLoop-1::nLoop]
#rok4 = (rok3)/rok3[0]
#
## plot
#
#
#ax1.plot(time[0:len(ros4)], ros4, linewidth=1.0, color='b', label="3rd TVD")
##ax1.set_xlabel("Time [s]")
##ax1.set_ylabel("(del ros/init ros) [-]")
##ax1.set_xlim([0,200])
##ax1.set_ylim([-0.04,0.01])
##ax1.hlines(0.0, 0, 200, color='black', linestyles='dotted')
#
#ax2.plot(time[0:len(rok4)], rok4, linewidth=1.0, color='b', label="3rd TVD")
##ax2.set_xlabel("Time [s]")
##ax2.set_ylabel("(rok/init rok) [-]")
##ax2.set_xlim([0,200])
##ax2.set_ylim([0.6,1.6])
#
##ax2.hlines(1.0, 0, 200, color='black',linestyles='dotted')
#
#
##4th redo
#with open("log_M0.1", "r") as f:
#	lines = f.read().splitlines()
#
#
#nLoop = 4
#
#M0=0.4
#
#time=[]
#ros =[]
#rok =[]
#
#for line in lines:
#	if "Total Time" in line:
#		time.append(M0*float(line.split(" ")[-1]))
#
#	if "ros" in line:
#		ros.append(float(line.split(" ")[-1]))
#
#	if "rok" in line:
#		rok.append(float(line.split(" ")[-1]))
#
#ros2 = np.array(ros)
#ros3 = ros2[nLoop-1::nLoop]
#ros4 = (ros3-ros3[0])/ros3[0]
#
#rok2 = np.array(rok)
#rok3 = rok2[nLoop-1::nLoop]
#rok4 = (rok3)/rok3[0]
#
## plot
#
#
#ax1.plot(time[0:len(ros4)], ros4, linewidth=1, color='r')
##ax1.set_xlabel("Time [s]")
##ax1.set_ylabel("(del ros/init ros) [-]")
##ax1.set_xlim([0,200])
##ax1.set_ylim([-0.04,0.01])
##ax1.hlines(0.0, 0, 200, color='black', linestyles='dotted')
#
#ax2.plot(time[0:len(rok4)], rok4, linewidth=1, color='r')
##ax2.set_xlabel("Time [s]")
##ax2.set_ylabel("(rok/init rok) [-]")
##ax2.set_xlim([0,200])
##ax2.set_ylim([0.6,1.6])
#
##ax2.hlines(1.0, 0, 200, color='black',linestyles='dotted')
#
#
#
#
plt.legend()
plt.show()



