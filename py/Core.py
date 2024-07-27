import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation  
import numpy as np
import scipy as sp
from tqdm import tqdm
from itertools import product

N = 30
R = 20
num_points = 500

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
    R = np.array(((c, -s), (s, c)))
    a,b = R@a0 + l0[0], R@b0 + l0[1]
    X = a[0]*m + b[0]*n
    Y = a[1]*m + b[1]*n
    return X,Y

def GaussLattice(x,y,X,Y,s=1):
    ''' A function with periodic Gaussians centered at X,Y lattice'''
    exponent = (np.subtract.outer(X,x))**2  + (np.subtract.outer(Y,y))**2
    # print(exponent.shape)
    explat = np.exp(exponent / (-2*(s**2)))     # lattice of exponentials
    return np.sum(explat, axis=(0,1))

def define_peak_variables(signal):
    # caller_frame = inspect.currentframe().f_back

    peaklocs, _ = sp.signal.find_peaks(signal)
    peaklocs2, _ = sp.signal.find_peaks(signal[peaklocs])
    troughlocs, _ = sp.signal.find_peaks(-signal)
    troughlocs2, _ = sp.signal.find_peaks(-signal[troughlocs])

    # caller_frame.f_locals['peaklocs'] = peaklocs
    # caller_frame.f_locals['peaklocs2'] = peaklocs2
    # caller_frame.f_locals['troughlocs'] = troughlocs
    # caller_frame.f_locals['troughlocs2'] = troughlocs2
    return peaklocs, peaklocs2, troughlocs, troughlocs2

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
        ax.set_title(r"along $\theta=$" +f"{t}" +r"$^\circ$ with $\Delta\theta=$" +f"{dt}" +r"$^\circ$")
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
            ax.set_title(r"along $\theta=$" +f"{t}" +r"$^\circ$at $\Delta\theta=$" +f"{dt}" +r"$^\circ$")
        else:
            plt.plot(r,signal_conv)
            plt.title(r"Smooth Diff of activity along $\theta=$" +f"{t}" +r"$^\circ$at $\Delta\theta=$" +f"{dt}" +r"$^\circ$")
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

def range_theta_data(ts=[0], dts=[0], ax=0):
    ranges = np.zeros((len(ts), len(dts), 2))
    # for i,t in tqdm(enumerate(ts), total=len(ts)):
    for i,j in tqdm(product(range(len(ts)), range(len(dts))), total=len(ts)*len(dts), leave=False):
        peaks, troughs = show_theta_var(ts[i], dts[j], ax[i,j], verbose=False)
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
        ax.set_title(f"Theta={ts[i]} deg")
        ax.plot(dts,range_theta_dat[i,:,0], color='r', label='Peak Range')
        ax.plot(dts,range_theta_dat[i,:,1], color='g', label='Trough Range')
        ax.set_xlabel('Diff of angle/ dt')
        ax.set_ylabel("Range")
    fig.set_tight_layout(True)
    plt.legend()
    plt.savefig('save.png')
    plt.show()

def range_theta_animator(t):
    dts = np.linspace(0,20,50)
    range_theta_dat = range_theta_data([t], dts, ax)
    ax.set_title(f"Theta={t} deg")
    line1.set_data( dts, range_theta_dat[0,:,0] )
    line2.set_data( dts, range_theta_dat[0,:,1] )
    return line1, line2

if __name__ == "__main__":
    fig, ax = plt.subplots(1)
    line1, = ax.plot([], [])
    line2, = ax.plot([], [])
    ax.set_xlabel('Diff of angle/ dt')
    ax.set_ylabel("Range")

    anim = FuncAnimation(fig, range_theta_animator,  
                               frames = 5, interval = 5)  
   
    # saves the animation in our desktop 
    # anim.save('growingCoil.mp4', writer = 'ffmpeg', fps = 30) 
    plt.show()
    


