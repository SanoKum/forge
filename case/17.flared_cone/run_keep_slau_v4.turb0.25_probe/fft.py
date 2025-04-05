import numpy as np # numpyをインポート
import matplotlib.pyplot as plt # matplotlibをインポート

dt = 5e-10*20 # [s]
f_s = 1.0/dt  # [1/s]
#N = int(f_s * t_fin) # サンプル数 [個]
#plot_probe_id_list = [0 , 1 , 2 , 3]
plot_probe_id_list = [0 , 2,  3]

xmax = 500   # kHz
ymax_p = 5.0 # kHz
ymax_v = 3.0  # kHz

startStep = 0

step_l = []
time_l = []
temp_summary = []
ux_summary   = []
uy_summary   = []
uz_summary   = []
pres_summary = []

temp_amp_summary = []
ux_amp_summary   = []
uy_amp_summary   = []
uz_amp_summary   = []
pres_amp_summary = []



for ip in plot_probe_id_list:
    
    temp_l = []
    ux_l   = []
    uy_l   = []
    uz_l   = []
    pres_l = []

    with open("point_probe_"+str(ip) + ".out") as f:
        irow = 0
        for line in f:
            if (irow == 0):
                irow += 1
                labels = line.split(",")
                continue

            value_l = line.split(",")
            if (int(value_l[0]) < startStep):
                continue

            step_l.append(int(value_l[0]))
            time_l.append(float(value_l[1]))
            temp_l.append(float(value_l[2]))
            pres_l.append(float(value_l[3]))
            ux_l.append(float(value_l[4]))
            uy_l.append(float(value_l[5]))
            uz_l.append(float(value_l[6]))

            irow += 1

    N = len(pres_l)

    ux_summary.append(ux_l)
    ux_fft = np.fft.fft(ux_l)
    freq = np.fft.fftfreq(N, d=dt)
    ux_amp = abs(ux_fft/(N/2))
    ux_amp_summary.append(ux_amp)
    
    uy_summary.append(uy_l)
    uy_fft = np.fft.fft(uy_l)
    freq = np.fft.fftfreq(N, d=dt)
    uy_amp = abs(uy_fft/(N/2))
    uy_amp_summary.append(uy_amp)

    uz_summary.append(uz_l)
    uz_fft = np.fft.fft(uz_l)
    freq = np.fft.fftfreq(N, d=dt)
    uz_amp = abs(uz_fft/(N/2))
    uz_amp_summary.append(uz_amp)
 
    pres_summary.append(pres_l)
    pres_fft = np.fft.fft(pres_l)
    freq = np.fft.fftfreq(N, d=dt)
    pres_amp = abs(pres_fft/(N/2))
    pres_amp_summary.append(pres_amp)
 

### Frequency vs Pressure ###
##################


fig = plt.figure()

ax1 = fig.add_subplot(1, 3, 1)

for ilp, ip in enumerate(plot_probe_id_list):
    ax1.plot(freq[1:int(N/2)]/1000, pres_amp_summary[ilp][1:int(N/2)], label="probe_"+str(ip)) 

ax1.legend()
ax1.set_xlabel("Frequency [kHz]")
ax1.set_ylabel("Pressure [Pa]")
ax1.set_xlim(0, xmax)
ax1.set_ylim(0, ymax_p)

### Frequency vs Ux ###
#######################


ax2 = fig.add_subplot(1, 3, 2)

for ilp, ip in enumerate(plot_probe_id_list):
    ax2.plot(freq[1:int(N/2)]/1000, ux_amp_summary[ilp][1:int(N/2)], label="probe_"+str(ip)) # A-f グラフのプロット

ax2.set_xlabel("Frequency [kHz]")
ax2.set_ylabel("Velocity x [m/s]")


ax2.legend()
ax2.set_xlim(0, xmax)
ax2.set_ylim(0, ymax_v)


### Time vs pressure ###
#######################

ax3 = fig.add_subplot(1, 3, 3)

for ilp, ip in enumerate(plot_probe_id_list):
    ax3.plot(time_l[1:N], pres_summary[ilp][1:N], label="probe_"+str(ip)) # A-f グラフのプロット

ax3.set_xlabel("Time [s]")
ax3.set_ylabel("Pressure [Pa]")


ax3.legend()
#ax3.set_xlim(0, xmax)


fig.tight_layout()


plt.show()

#
#
#        irow += 1
#
#
#
#### FFT: tの関数をfの関数にする ###
#y_fft = np.fft.fft(y) # 離散フーリエ変換
#freq = np.fft.fftfreq(N, d=dt) # 周波数を割り当てる（※後述）
#Amp = abs(y_fft/(N/2)) # 音の大きさ（振幅の大きさ）
