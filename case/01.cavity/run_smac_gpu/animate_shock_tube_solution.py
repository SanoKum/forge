# coded by Kumpei Sano. 2022 05 04
import numpy as np
from scipy import optimize
#from matplotlib import pyplot
import matplotlib.pyplot as plt
import glob
import subprocess, sys

# Referance
# https://en.wikipedia.org/wiki/Sod_shock_tube
# https://physics.stackexchange.com/questions/423758/how-to-get-exact-solution-to-sod-shock-tube-test

if __name__ == '__main__':


    # ----------------------
    # *** make animation ***
    # ----------------------
    #command = ['ffmpeg' , '-framerate', '10', '-i', 'cavity_flow.%04d.png', '-vcodec', 'libx264', '-pix_fmt', 'yuv420p', '-r', '60', 'output.mp4']
    command = ['ffmpeg' , '-framerate', '10', '-i', 'cavity_flow.%04d.png', '-r', '60', 'output.mp4']
    print(command)
    cp = subprocess.run(command)
    if cp.returncode != 0:
        print('ffmpeg failded', file=sys.stderr)
        sys.exit(1)

    #command = ['ffmpeg' , '-loop', '1', '-i', "cavity_flow.0019.png", '-vcodec', 'libx264', '-pix_fmt', 'yuv420p', '-t', '5', '-r', '60', 'outlast.mp4']
    command = ['ffmpeg' , '-loop', '1', '-i', "cavity_flow.0019.png", '-t', '5', '-r', '60', 'outlast.mp4']
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


