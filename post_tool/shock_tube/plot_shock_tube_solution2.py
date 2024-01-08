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

    #for fname1, fname2 in zip(file_names_1st, file_names_MUSCL):
    fname1 = file_names_1st[-1]

    #ctime = fname1[13:21]
    ctime = fname1.split("_")[-1][0:-4]
    time = float(ctime)
    result = shock_ana.shock_tube_analytical_solution(time=time)

    x       = result["x"]
    u       = result["u"]
    p       = result["p"]
    r       = result["r"]
    int_eng = result["int_eng"]

    #result_cfd = np.genfromtxt("endVal.dat")
    result_cfd1 = np.genfromtxt(fname1)
    x_cfd       = result_cfd1[:,0]
    r_cfd       = result_cfd1[:,1]
    u_cfd       = result_cfd1[:,2]
    e_cfd       = result_cfd1[:,3]
    p_cfd       = result_cfd1[:,4]
    #int_eng_cfd = result_cfd1[:,5]
    switch_cfd  = result_cfd1[:,5]
    beta_cfd  = result_cfd1[:,6]

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
    fig, axes = plt.subplots(ncols=3, figsize=(16,6), sharex=True )

    axes[0].plot(x, r, linewidth=1, label="Exact", color="black")
    #axes[0,0].plot(x_cfd, r_cfd, label="FR Solver 1", linestyle='None', marker="o", markersize=2, color="red")
    axes[0].plot(x_cfd, r_cfd, label="FR Solver 1", linewidth=2, color="red")
    #axes[0,0].plot(x_cfd2, r_cfd2, label="Roe Solver MUSCL" , linestyle='None', marker="^", markersize=3, color="blue")
    axes[0].set_xlabel("x")
    axes[0].set_xlim([0,1])
    axes[0].set_ylabel("Density")
    #axes[0,0].set_ylim([-0.05,1.05])
    axes[0].legend()

    axes[1].plot(x, p, linewidth=1, label="Exact", color="black")
    #axes[1,0].plot(x_cfd, p_cfd, label="Roe Solver")
    #axes[1,0].plot(x_cfd, p_cfd, label="FR Solver 1", linestyle='None', marker="o", markersize=2, color="red")
    axes[1].plot(x_cfd, p_cfd, label="FR Solver 1", linewidth=2, color="red")
    #axes[1,0].plot(x_cfd2, p_cfd2, label="Roe Solver MUSCL" , linestyle='None', marker="^", markersize=3, color="blue")
    axes[1].set_xlabel("x")
    axes[1].set_xlim([0,1])
    axes[1].set_ylabel("Pressure")
    #axes[1,0].set_ylim([-0.05,1.05])
    #axes[1,0].legend()

    axes[2].plot(x, u, linewidth=1, label="Exact", color="black")
    #axes[2,0].plot(x_cfd, u_cfd, label="FR Solver 1", linestyle='None', marker="o", markersize=2, color="red")
    axes[2].plot(x_cfd, u_cfd, label="FR Solver 1", linewidth=2, color="red")
    #axes[0,1].plot(x_cfd2, u_cfd2, label="Roe Solver MUSCL" , linestyle='None', marker="^", markersize=3, color="blue")
    axes[2].set_xlabel("x")
    axes[2].set_xlim([0,1])
    axes[2].set_ylabel("Velocity")
    #axes[2,0].set_ylim([-0.05,1.05])
    #axes[0,1].legend()
#
#    #axes[1,1].plot(x, int_eng, linewidth=1, label="Exact", color="black")
#    axes[0,1].plot(x_cfd , switch_cfd, label="FR Solver 1", linestyle='None', marker="o", markersize=2, color="red")
#    #axes[1,1].plot(x_cfd2, int_eng_cfd2, label="Roe Solver MUSCL"   , linestyle='None', marker="^", markersize=3, color="blue")
#    axes[0,1].set_xlabel("x")
#    axes[0,1].set_xlim([0,1])
#    axes[0,1].set_ylabel("Internal Energy")
#    #axes[1,1].set_ylim([1.7,3.0 ])
#
#    #axes[1,1].plot(x, int_eng, linewidth=1, label="Exact", color="black")
#    axes[1,1].plot(x_cfd , beta_cfd, label="FR Solver 1", linestyle='None', marker="o", markersize=2, color="red")
#    #axes[1,1].plot(x_cfd2, int_eng_cfd2, label="Roe Solver MUSCL"   , linestyle='None', marker="^", markersize=3, color="blue")
#    axes[1,1].set_xlabel("x")
#    axes[1,1].set_xlim([0,1])
#    axes[1,1].set_ylabel("Internal Energy")
#    #axes[1,1].set_ylim([1.7,3.0 ])
#

    fig.tight_layout(rect=[0,0,1,0.96])

    fig.suptitle("Time = " + ctime)

    plt.show()
    fig.savefig("last.png")
    plt.close()


