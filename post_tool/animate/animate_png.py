# coded by Kumpei Sano. 2024 01 08
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import shock_tube_analytical_solution as shock_ana
import glob
import subprocess, sys

# Referance

if __name__ == '__main__':

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



