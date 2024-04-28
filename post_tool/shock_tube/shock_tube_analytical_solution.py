# coded by Kumpei Sano. 2022 05 04
import numpy as np
from scipy import optimize
#from matplotlib import pyplot
import matplotlib.pyplot as plt

# Referance
# https://en.wikipedia.org/wiki/Sod_shock_tube
# https://physics.stackexchange.com/questions/423758/how-to-get-exact-solution-to-sod-shock-tube-test

def shock_tube_analytical_solution(gamma=1.4, r1=1., p1=1., u1=0.,r5=0.125, p5=0.1, u5=0., x_mid=0.5, time=0.2):
    # ------------
    # *** calc ***
    # ------------
    c1 = np.sqrt(gamma*p1/r1)
    c5 = np.sqrt(gamma*p5/r5)

    zeta = (gamma-1)/(gamma+1)
    beta = (gamma-1)/(2*gamma)

    f = lambda x: (x-p5)*np.sqrt((1.0-zeta)/(r5*(x+zeta*p5))) \
                  -(p1**beta -x**beta)*np.sqrt((1-zeta**2)*p1**(1.0/gamma) /(zeta**2*r1))

    p3 = optimize.brentq(f, 0, 1)
    u3 = u5 + (p3-p5)/np.sqrt(r5*0.5*((gamma+1.0)*p3 + (gamma-1)*p5))
    r3 = r1*(p3/p1)**(1.0/gamma)

    p4 = p3
    u4 = u3
    r4 = r5*(p4+zeta*p5)/(p5+zeta*p4)


    u_contact = u4
    x_contact = u_contact*time + 0.5
    u_shock = u4*r4/(r4-r5)
    x_shock = u_shock*time + 0.5
    u_wave = -c1
    x_wave = u_wave*time + 0.5

    x_sep = (u3/2.0*(gamma+1) - c1)*time + x_mid

    if (time!=0):
        x_range2 = np.linspace(x_wave, x_sep, 10)
        u2 = 2.0/(gamma+1.0)*(c1 + (x_range2-x_mid)/time)
    else:
        x_range2 = x_mid
        u2 = u1
    r2 = r1*(1.0-(gamma-1.0)/2.0*u2/c1)**(2.0/(gamma-1))
    p2 = p1*(1.0-0.5*(gamma-1.0)*u2/c1)**(2.0*gamma/(gamma-1))


    # ---------------
    # *** summary ***
    # ---------------
    x = np.copy(x_range2)
    u = np.copy(u2)
    r = np.copy(r2)
    p = np.copy(p2)


    #add region 1 left
    x = np.insert(x, 0, 0.0)
    u = np.insert(u, 0, u1)
    r = np.insert(r, 0, r1)
    p = np.insert(p, 0, p1)

    #add region 3 right
    x = np.append(x, x_contact)
    u = np.append(u, u3)
    r = np.append(r, r3)
    p = np.append(p, p3)

    #add region 4 left
    x = np.append(x, x_contact)
    u = np.append(u, u4)
    r = np.append(r, r4)
    p = np.append(p, p4)

    #add region 4 right
    x = np.append(x, x_shock)
    u = np.append(u, u4)
    r = np.append(r, r4)
    p = np.append(p, p4)

    #add region 5 left
    x = np.append(x, x_shock)
    u = np.append(u, u5)
    r = np.append(r, r5)
    p = np.append(p, p5)

    #add region 5 right
    x = np.append(x, 1.0)
    u = np.append(u, u5)
    r = np.append(r, r5)
    p = np.append(p, p5)

    #set internal energy
    int_eng = p/(r*(gamma-1.0))

    result = {}

    result["x"] = x
    result["u"] = u
    result["r"] = r
    result["p"] = p
    result["int_eng"] = int_eng

    return result


if __name__ == '__main__':

    result = shock_tube_analytical_solution()

    x = result["x"]
    u = result["u"]
    p = result["p"]
    r = result["r"]

    # ------------
    # *** Plot ***
    # ------------

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(10,3), sharex=True, tight_layout=True)
    ax1.plot(x, r, linewidth=2)
    ax1.set_xlabel("x")
    ax1.set_ylabel("Density")

    ax2.plot(x, p, linewidth=2)
    ax2.set_xlabel("x")
    ax2.set_ylabel("Pressure")

    ax3.plot(x, u, linewidth=2)
    ax3.set_xlabel("x")
    ax3.set_ylabel("Velocity")

    plt.show()
