################################################################################
## The code defines the basic functions used by other files.                  ##
## Written by P.D.                                                            ##
################################################################################

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation  
import numpy as np
import scipy as sp
from tqdm import tqdm
from itertools import product

np.random.seed(42)

# Settings for continuos plots (old)
# N = 20
# R = 10
# num_points = 100

# Settings for scatter plots (new)
# N = 20
# R = 10
# num_points = 11

def Lattice(l=1,l0=[0,0],t0=0,N=3):
    '''Return an array of x coordinates and y coordinates that create a hexagonal lattice
      of N^2 points of length l, rotated by theta t0 and displaced by l0  '''
    mlist = np.linspace(-N,N,2*N+1)
    nlist = np.linspace(-N,N,2*N+1)
    m,n = np.meshgrid(mlist,nlist)
    a0 = np.array([(3**0.5)/2,1/2]) *l
    b0 = np.array([(3**0.5)/2,-1/2]) *l

    # apply rotation
    theta = np.radians(t0)
    c, s = np.cos(theta), np.sin(theta)
    Rot = np.array(((c, -s), (s, c)))
    a,b = Rot@a0 , Rot@b0 
    X = a[0]*m + b[0]*n + l0[0]
    Y = a[1]*m + b[1]*n + l0[1]
    return X,Y

def GaussLattice(x,y,X,Y,s=1):
    ''' A function with periodic Gaussians centered at X,Y lattice'''
    exponent = (np.subtract.outer(x,X))**2  + (np.subtract.outer(y,Y))**2
    # print(exponent.shape)
    explat = np.exp(exponent / (-2*(s**2)))     # lattice of exponentials
    return np.sum(explat, axis=(-2,-1))

def define_peak_variables(signal):
    ''' Define the peak of peaks and trough of troughs for further range calculation'''
    peaklocs, _ = sp.signal.find_peaks(signal)
    peaklocs2, _ = sp.signal.find_peaks(signal[peaklocs])
    troughlocs, _ = sp.signal.find_peaks(-signal)
    troughlocs2, _ = sp.signal.find_peaks(-signal[troughlocs])
    return peaklocs, peaklocs2, troughlocs, troughlocs2

def theta_var_data(t=0,dt=0):
    r = np.linspace(-R,R,num_points)
    theta = np.radians(t)
    x = r*np.cos(theta)
    y = r*np.sin(theta)

    X1,Y1 = Lattice(l=1,t0=0,N=N)
    X2,Y2 = Lattice(l=1,t0=-dt,N=N)

    Z1 = GaussLattice(x,y,X1,Y1,s=0.2)
    Z2 = GaussLattice(x,y,X2,Y2,s=0.2)
    signal = Z2-Z1
    peaklocs, peaklocs2, troughlocs, troughlocs2 = define_peak_variables(signal)
    
    return r[peaklocs][peaklocs2], r[troughlocs][troughlocs2]

def show_theta_var(t=0,dt=0,ax=0, conv=False, verbose=True):
    '''Plots the variation of grid response along the theta=t angle when 
    the grids have a difference in theta of dt'''
    r = np.linspace(-R,R,num_points)
    theta = np.radians(t)
    x = r*np.cos(theta)
    y = r*np.sin(theta)

    X1,Y1 = Lattice(l=1,t0=0,N=N)
    X2,Y2 = Lattice(l=1,t0=-dt,N=N)

    Z1 = GaussLattice(x,y,X1,Y1,s=0.2)
    Z2 = GaussLattice(x,y,X2,Y2,s=0.2)
    signal = Z2-Z1
    peaklocs, peaklocs2, troughlocs, troughlocs2 = define_peak_variables(signal)

    if conv:
        signal_conv = sp.ndimage.gaussian_filter(signal,sigma=8)
        if ax!=0:
            ax.plot(r,signal_conv)
            ax.set_title(r"along $\theta=$" +f"{t}" +r"$^\circ$at $\Delta\theta=$" +f"{dt}" +r"$^\circ$")
        else:
            plt.plot(r,signal_conv)
            plt.title(r"Smooth Diff of activity along $\theta=$" +f"{t}" +r"$^\circ$at $\Delta\theta=$" +f"{dt}" +r"$^\circ$")
        return r[peaklocs][peaklocs2], r[troughlocs][troughlocs2]

    if verbose:
        print(f"\n For {t} and {dt} Monotonicity of peaks break at {r[peaklocs][peaklocs2]}")
        print(f"Monotonicity of troughs break at {r[troughlocs][troughlocs2]} \n")

    if ax!=0:
        ax.plot(r,signal)
        ax.scatter(r[peaklocs][peaklocs2], signal[peaklocs][peaklocs2], color='r')
        ax.scatter(r[troughlocs][troughlocs2], signal[troughlocs][troughlocs2], color='g')
        ax.set_title(r"along $\theta=$" +f"{t:.2f}" +r"$^\circ$ with $\Delta\theta=$" +f"{dt:.2f}" +r"$^\circ$")
    else:
        plt.plot(r,signal)
        plt.scatter(r[peaklocs][peaklocs2], signal[peaklocs][peaklocs2], color='r')
        plt.scatter(r[troughlocs][troughlocs2], signal[troughlocs][troughlocs2], color='g')
        plt.title(r"along $\theta=$" +f"{t}" +r"$^\circ$ with $\Delta\theta=$" +f"{dt}" +r"$^\circ$ \n")
    return r[peaklocs][peaklocs2], r[troughlocs][troughlocs2]

def show_dist_var(t=0,dl=0,ax=0, conv=False, verbose=True):
    '''Plots the variation of grid response along the theta=t angle when 
    the grids have a difference in lambda of dl'''
    r = np.linspace(-R,R,num_points)
    theta = np.radians(t)
    x = r*np.cos(theta)
    y = r*np.sin(theta)

    X1,Y1 = Lattice(l=1,t0=0,N=N)
    X2,Y2 = Lattice(l=1+dl,t0=0,N=N)
    # X3,Y3 = Lattice(l=1,t0=10,N=10)

    Z1 = GaussLattice(x,y,X1,Y1,s=0.2)
    Z2 = GaussLattice(x,y,X2,Y2,s=0.2)
    signal = Z2-Z1
    peaklocs, peaklocs2, troughlocs, troughlocs2 = define_peak_variables(signal)

    if conv:
        signal_conv = sp.ndimage.gaussian_filter(signal,sigma=8)
        if ax!=0:
            ax.plot(r,signal_conv)
            ax.set_title(r"along $\theta=$" +f"{t}" +r"$^\circ$at $\Delta\lamda=$" +f"{dl}")
        else:
            plt.plot(r,signal_conv)
            plt.title(r"Smooth Diff of activity along $\theta=$" +f"{t}" +r"$^\circ$at $\Delta\lambda=$" +f"{dl}")
        return r[peaklocs][peaklocs2], r[troughlocs][troughlocs2]

    if verbose:
        print(f"\n For {t} and {dl} Monotonicity of peaks break at {r[peaklocs][peaklocs2]}")
        print(f"Monotonicity of troughs break at {r[troughlocs][troughlocs2]} \n")

    if ax!=0:
        ax.plot(r,signal)
        ax.scatter(r[peaklocs][peaklocs2], signal[peaklocs][peaklocs2], color='r')
        ax.scatter(r[troughlocs][troughlocs2], signal[troughlocs][troughlocs2], color='g')
        ax.set_title(r"along $\theta=$" +f"{t}" +r"$^\circ$ with $\Delta\lambda=$" +f"{dl}")
    else:
        plt.plot(r,signal)
        plt.scatter(r[peaklocs][peaklocs2], signal[peaklocs][peaklocs2], color='r')
        plt.scatter(r[troughlocs][troughlocs2], signal[troughlocs][troughlocs2], color='g')
        plt.title(r"along $\theta=$" +f"{t}" +r"$^\circ$ with $\Delta\lambda=$" +f"{dl}")
    return r[peaklocs][peaklocs2], r[troughlocs][troughlocs2]

def range_theta_data(ts=[0], dts=[0]):
    ranges = np.zeros((len(ts), len(dts), 2))
    # for i,t in tqdm(enumerate(ts), total=len(ts)):
    for i,j in tqdm(product(range(len(ts)), range(len(dts))), total=len(ts)*len(dts), leave=False):
        peaks, troughs = theta_var_data(ts[i], dts[j])
        if len(peaks[peaks>0.1]) != 0:
            ranges[i,j,0] = peaks[peaks>0.1][0]
        if len(troughs[troughs>0.1]) != 0:
            ranges[i,j,1] = troughs[troughs>0.1][0]
    return ranges

def range_theta_fullpage():
    ts = np.linspace(0,48,25)
    dts = np.linspace(0,20,50)
    fig, ax = plt.subplots(len(ts), len(dts), figsize=(5*len(ts),5*len(dts)))
    range_theta_dat = range_theta_data(ts, dts, ax)
    plt.close()

    num_rows = 5
    fig, ax = plt.subplots(num_rows, len(ts)//num_rows, figsize=(2*len(ts)/num_rows, 2.5*num_rows), sharey=True)
    fig.suptitle("Range along different directions")
    for i, ax in enumerate(ax.flatten()):
        # plt.title(f"Along angle theta={ts[i]} deg")
        # plt.plot(dts,range_theta_dat[i,:,0], color='r')
        # plt.plot(dts,range_theta_dat[i,:,1], color='g')
        # plt.xlabel('dt / diff of angle')
        # plt.ylabel("range")
        # plt.show()
        ax.set_title(f"Theta={ts[i]:.2f} deg")
        ax.plot(dts,range_theta_dat[i,:,0], color='r', label='Peak Range')
        ax.plot(dts,range_theta_dat[i,:,1], color='g', label='Trough Range')
        ax.set_xlabel('Diff of angle/ dt')
        ax.set_ylabel("Range")
    fig.set_tight_layout(True)
    plt.legend()
    plt.savefig('save.png')
    plt.show()

def range_theta_animator(t):
    ax.set_title(f"Theta={ts[t]} deg")
    # print(dts, range_theta_dat[0,:,0])
    # plt.plot( dts, range_theta_dat[0,:,0] )
    while len(plots)!=0:
        plot = plots.pop()
        # print(plot)
        plot[0].remove()

    plots.append( ax.plot( dts, range_theta_dat[t,:,0], color='r', label='Peak Range' ) )
    plots.append( ax.plot( dts, range_theta_dat[t,:,1], color='g', label='Trough Range' ) )
    plt.legend()

def single_theta_variation_animator(t):
    ax.clear()
    ax.set_title(f"Theta={ts[0]}, dt={dts[t]} deg")
    ax.set_ylim(-1,1)
    show_theta_var(t=ts[0],dt=dts[t],ax=ax, conv=False, verbose=False)

def fullpage_theta_variation_animator(t):
    for i, ax in enumerate(axs.flatten()):
        ax.clear()
        ax.set_ylim(-1,1)
        show_theta_var(t=ts[i],dt=dts[t],ax=ax, conv=False, verbose=False)

if __name__ == "__main__":
    # For fullpage_theta_variation_animator
    plots = []
    num = 100
    ts = np.linspace(0,55,12)
    dts = np.linspace(1,10,num)
    fig, axs = plt.subplots(len(ts)//4, 4, figsize=(19.2, 10.8))
    fig.set_tight_layout(True)

    range_theta_dat = range_theta_data(ts, dts)

    anim = FuncAnimation(fig, fullpage_theta_variation_animator,  
                               frames = num, interval = 200)  
   
    # saves the animation in our desktop 
    # plt.show()
    anim.save('fullpage.mkv', writer='ffmpeg') 

# if __name__ == "__main__":
#     # For single theta_variation_animator
#     plots = []
#     fig = plt.figure()
#     ax = fig.add_subplot(1,1,1)

#     num = 500
#     ts = [45]
#     dts = np.linspace(0,20,num)
#     range_theta_dat = range_theta_data(ts, dts)

#     anim = FuncAnimation(fig, single_theta_variation_animator,  
#                                frames = num, interval = 200)  
   
#     # saves the animation in our desktop 
#     # plt.show()
#     anim.save('45.mp4', writer='ffmpeg') 

# if __name__ == "__main__":
#     # For range_theta_animator
#     plots = []
#     fig = plt.figure()
#     ax = fig.add_subplot(1,1,1)
#     ax.set_xlabel('Diff of angle/ dt')
#     ax.set_ylabel("Range")

#     num = 500
#     ts = np.linspace(0,10,num)
#     dts = np.linspace(0,20,100)
#     range_theta_dat = range_theta_data(ts, dts)

#     anim = FuncAnimation(fig, range_theta_animator,  
#                                frames = num, interval = 200)  
   
#     # saves the animation in our desktop 
#     # plt.show()
#     anim.save('5.mp4', writer='ffmpeg') 
    
# if __name__ == "__main__":
#     range_theta_fullpage()

