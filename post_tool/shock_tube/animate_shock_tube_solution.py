# coded by Kumpei Sano. 2022 05 04
import numpy as np
from scipy import optimize
#from matplotlib import pyplot
import matplotlib.pyplot as plt
import shock_tube_analytical_solution as shock_ana
import glob
import subprocess, sys

# Referance
# https://en.wikipedia.org/wiki/Sod_shock_tube
# https://physics.stackexchange.com/questions/423758/how-to-get-exact-solution-to-sod-shock-tube-test

if __name__ == '__main__':

    file_names_1st = sorted(glob.glob('./results_time_*.dat'))
    #file_names_MUSCL = sorted(glob.glob('./MUSCL/results_time_*.dat'))

    png_names = []

    i = 0
    #for fname1, fname2 in zip(file_names_1st, file_names_MUSCL):
    for fname1 in file_names_1st:

        cstep = str(i).zfill(6)
        #ctime = fname1[13:21]
        ctime = fname1.split("_")[-1][0:-4]
        time = float(ctime)
        result = shock_ana.shock_tube_analytical_solution(time=time)

        x       = result["x"]
        u       = result["u"]
        p       = result["p"]
        r       = result["r"]

        #result_cfd = np.genfromtxt("endVal.dat")
        result_cfd1 = np.genfromtxt(fname1)
        x_cfd       = result_cfd1[:,0]
        r_cfd       = result_cfd1[:,1]
        u_cfd       = result_cfd1[:,2]
        e_cfd       = result_cfd1[:,3]
        p_cfd       = result_cfd1[:,4]

        #result_cfd2 = np.genfromtxt(fname2)
        #x_cfd2       = result_cfd2[:,0]
        #r_cfd2       = result_cfd2[:,1]
        #u_cfd2       = result_cfd2[:,2]
        #e_cfd2       = result_cfd2[:,3]
        #p_cfd2       = result_cfd2[:,4]
        #int_eng_cfd2 = result_cfd2[:,5]


        # ------------
        # *** Plot ***
        # ------------
        fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(12,6.5), sharex=True )

        axes[0,0].plot(x, r, linewidth=1, label="Exact", color="black")
        axes[0,0].plot(x_cfd, r_cfd, label="FR Solver 1", linewidth=2, color="red")
        #axes[0,0].plot(x_cfd2, r_cfd2, label="Roe Solver MUSCL" , linestyle='None', marker="^", markersize=3, color="blue")
        axes[0,0].set_xlabel("x")
        axes[0,0].set_xlim([0,1])
        axes[0,0].set_ylabel("Density")
        axes[0,0].set_ylim([-0.05,1.05])
        axes[0,0].legend()

        axes[1,0].plot(x, p, linewidth=1, label="Exact", color="black")
        #axes[1,0].plot(x_cfd, p_cfd, label="Roe Solver")
        axes[1,0].plot(x_cfd, p_cfd, label="FR Solver 1", linewidth=2, color="red")
        #axes[1,0].plot(x_cfd2, p_cfd2, label="Roe Solver MUSCL" , linestyle='None', marker="^", markersize=3, color="blue")
        axes[1,0].set_xlabel("x")
        axes[1,0].set_xlim([0,1])
        axes[1,0].set_ylabel("Pressure")
        axes[1,0].set_ylim([-0.05,1.05])
        #axes[1,0].legend()

        axes[0,1].plot(x, u, linewidth=1, label="Exact", color="black")
        #axes[0,1].plot(x_cfd, u_cfd, label="JST Solver 1", linestyle='None', marker="o", markersize=3, color="red")
        axes[0,1].plot(x_cfd, u_cfd, label="FR Solver 1", linewidth=2, color="red")
        #axes[0,1].plot(x_cfd2, u_cfd2, label="Roe Solver MUSCL" , linestyle='None', marker="^", markersize=3, color="blue")
        axes[0,1].set_xlabel("x")
        axes[0,1].set_xlim([0,1])
        axes[0,1].set_ylabel("Velocity")
        axes[0,1].set_ylim([-0.05,1.05])
        #axes[0,1].legend()

        #axes[1,1].plot(x, int_eng, linewidth=1, label="Exact", color="black")
        #axes[1,1].plot(x_cfd , int_eng_cfd , label="JST Solver 1", linestyle='None', marker="o", markersize=3, color="red")
        ##axes[1,1].plot(x_cfd2, int_eng_cfd2, label="Roe Solver MUSCL"   , linestyle='None', marker="^", markersize=3, color="blue")
        #axes[1,1].set_xlabel("x")
        #axes[1,1].set_xlim([0,1])
        #axes[1,1].set_ylabel("Internal Energy")
        #axes[1,1].set_ylim([1.7,3.0 ])

        fig.tight_layout(rect=[0,0,1,0.96])

        fig.suptitle("Time = " + ctime)

        fig.savefig("step_"+cstep+".png")
        plt.close()

        png_names.append("step_"+cstep+".png")

        i = i + 1

    # ---------------------
    # *** save last png ***
    # ---------------------
    src=png_names[-1]
    dst="last.png"
    cmd='copy "%s" "%s"' % (src, dst)

    status = subprocess.call(cmd, shell=True)

    if status != 0:
        if status < 0:
            print("Killed by signal", status)
        else:
            print("Command failed with return code - ", status)
    else:
        print('Execution of %s passed!\n' % cmd)


    # ----------------------
    # *** make animation ***
    # ----------------------
    command = ['ffmpeg' , '-framerate', '10', '-i', 'step_%06d.png', '-vcodec', 'libx264', '-pix_fmt', 'yuv420p', '-r', '60', 'output.mp4']
    print(command)
    cp = subprocess.run(command)
    if cp.returncode != 0:
        print('ffmpeg failded', file=sys.stderr)
        sys.exit(1)

    command = ['ffmpeg' , '-loop', '1', '-i', png_names[-1], '-vcodec', 'libx264', '-pix_fmt', 'yuv420p', '-t', '5', '-r', '60', 'outlast.mp4']
    print(command)
    cp = subprocess.run(command)
    if cp.returncode != 0:
        print('ffmpeg failded', file=sys.stderr)
        sys.exit(1)

    with open('list.txt', mode='w') as f:
        f.write('file output.mp4\n')
        f.write('file outlast.mp4')

    command = ['ffmpeg' , '-f', 'concat', '-i', 'list.txt', '-c', 'copy', 'output2.mp4']
    print(command)
    cp = subprocess.run(command)
    if cp.returncode != 0:
        print('ffmpeg failded', file=sys.stderr)
        sys.exit(1)

    # ----------------------
    # *** delete results ***
    # ----------------------
    #file_names = sorted(glob.glob('results_time_*.dat'))

    command = ['rm'] + png_names 
    cp = subprocess.run(command)
    if cp.returncode != 0:
        print('rm failded', file=sys.stderr)
        sys.exit(1)



