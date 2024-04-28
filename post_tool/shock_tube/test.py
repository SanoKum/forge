#import numpy as np
#
#u = [0.0000000000000000     ,  0.91597047605606596    ,   0.93179751333905836    ,   0.95635980309124391    ,   0.98667010404392763    ,    0.0000000000000000]
#u = np.array(u)
#
#xl = [-1.0 , -0.861136 , -0.339981 , 0.339981, 0.861136 , 1.0]
#xl = np.array(xl)
#
#x_in = 1.0
#nSolPt = 4
#
#fai = np.ones(nSolPt)
#
#for i in range(nSolPt):
#    for j in range(nSolPt):
#        fai[i] = fai[i]*(x_in - xl[j])