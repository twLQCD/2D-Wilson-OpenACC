import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    data = np.loadtxt('fDtime_LX64_LY64_B6.000000_M0.000000_tau1.000000_nHMCstep25.dat')
    timefD = np.sum(data)
    
    data = np.loadtxt('trajtime_LX64_LY64_B6.000000_M0.000000_tau1.000000_nHMCstep25.dat')
    timetraj = np.sum(data)
    
    data = np.loadtxt('data_LX64_LY64_B6.000000_M0.000000_tau1.000000_nHMCstep25.dat')
    totaltime = data[-1,1]
    
    data = np.loadtxt('invtime.dat')
    invtime = np.sum(data)
    
    data = np.loadtxt('mvpstime.dat')
    mvptime = np.sum(data)
    
    
    print("total time for HMC is: ")
    print(totaltime)
    print("time spent performing trajectories: ")
    print(timetraj)
    print("time spent calculating fermion force: ")
    print(timefD)
    print("time spent inverting Dirac operator")
    print(invtime)
    print("time spent performing matrix vector products: ")
    print(mvptime)
    print("% spent performing trajectories: ")
    print(timetraj/totaltime * 100)
    print("% spent calculation fermion force: ")
    print(timefD/totaltime * 100)
    print("% spent inverting")
    print(invtime/totaltime * 100)
    print("% spent doing mvps: ")
    print(mvptime/totaltime * 100)
    
    plt.bar(1,((timetraj/totaltime)*100), width = 0.5)
    plt.bar(2,((timefD/totaltime)*100), width = 0.5)
    plt.bar(3,((invtime/totaltime)*100), width = 0.5)
    plt.bar(4,((mvptime/totaltime)*100), width = 0.5)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
    plt.ylabel("% Time of HMC")
    plt.legend(('MD Trajectory','Fermion Force','Wilson-Dirac Inversion','Matrix-Vector Products'), loc='lower left')
    plt.show()