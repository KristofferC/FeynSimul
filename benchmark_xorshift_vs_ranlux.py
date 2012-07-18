import numpy as np
import pylab as pl
from matplotlib import rc
from matplotlib.ticker import *
import matplotlib.pyplot as plt

xs = np.array([0.000835, 0.021259, 0.077566, 2.113641])   # xorshift single
xd = np.array([0.000833, 0.021227, 0.077535, 2.116926])   # xorshift double

r0s = np.array([0.002772, 0.013913, 0.124960, 1.234409])  # ranlux0 single
r0d = np.array([0.003888, 0.024675, 0.233399, 2.315815])  # ranlux0 double

r1s = np.array([0.003372, 0.019870, 0.184641, 1.833603])  # ranlux1 single
r1d = np.array([0.005065, 0.036373, 0.349675, 3.486245])  # ranlux1 double

r2s = np.array([0.005217, 0.037840, 0.364948, 3.635847])  # ranlux2 single
r2d = np.array([0.008705, 0.072788, 0.715689, 7.143060])  # ranlux2 double

r3s = np.array([0.008220, 0.068618, 0.674307, 6.730529])  # ranlux3 single
r3d = np.array([0.014878, 0.134723, 1.334791, 13.330784)] # ranlux3 double

r4s = np.array([0.012606, 0.113489, 1.124352, 11.224786]) # ranlux4 single
r4d = np.array([0.023882, 0.225516, 2.243392, 22.421897]) # ranlux4 double

logplot = True   # Log the y axis or not
bigimage = True  # 2 column or 1 column image
plotstep = 16    # Set this to how many points you want to skip
plotoffset = 15  # Offset for the first point to plot, set to one less than
                 # plotstep to make sure that you get global maxima
xr = 3           # Max for x axis (times 10**4)


if bigimage:
    msize = 40
    gtx470label = 'GPU: nVidia GeForce GTX 470'
    gtx560label = 'GPU: nVidia GeForce GTX 560 Ti'
    hd5850label = 'GPU: ATI Radeon HD 5850'
    gt9600label = 'GPU: nVidia GeForce 9600 GT'
    cpulabel = 'CPU: Intel Core i5-2500K CPU @ 4GHz'
    fig_width_pt = 512.0
    legendloc = 4
    fsize=16
else:
    msize = 20
    gtx470label = 'GTX470'
    gtx560label = 'GTX560'
    hd5850label = 'HD5850'
    gt9600label = 'GT9600'
    cpulabel = 'CPU'
    fig_width_pt = 246.0
    legendloc = 4
    fsize=10

# Calculate figure size using golden ratio method
inches_per_pt = 1.0/72.27                    # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0           # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt       # width in inches
fig_height = fig_width*golden_mean           # height in inches
fig_height = fig_height + 40 * inches_per_pt # Add space for title
fig_size =  [fig_width,fig_height]

# Set the RC LaTeX parameters
params = {'backend': 'ps',
         'axes.labelsize': fsize,
         'text.fontsize': 2*fsize,
         'legend.fontsize': 0.8*fsize,
         'xtick.labelsize': 0.8*fsize,
         'ytick.labelsize': 0.8*fsize,
         'text.usetex': True,
         'figure.figsize': fig_size}

pl.rcParams.update(params)

pl.figure(1)
if logplot:
    pl.semilogy()

# Dict with settings for different platforms
gpus = {'gtx470': {'Label': gtx470label,'Show': True, 'Marker': '^',
                   'MarkerSize': msize, 'Color': 'b', 'norm': None,
                   'cmap': None, 'vmin': None, 'vmax': None, 'Alpha': 1.0,
                   'LineWidths': 0.1, 'Verts': None},
        'gtx560ti': {'Label': gtx560label, 'Show': False, 'Marker': '>',
                     'MarkerSize': msize, 'Color': 'tomato', 'norm': None,
                     'cmap': None, 'vmin': None, 'vmax': None, 'Alpha': 0.5,
                     'LineWidths': 0.1, 'Verts': None},
        'hd5850': {'Label': hd5850label, 'Show': True, 'Marker': '<',
                   'MarkerSize': msize, 'Color': 'g', 'norm': None,
                   'cmap': None, 'vmin': None, 'vmax': None, 'Alpha': 0.5,
                   'LineWidths': 0.1, 'Verts': None},
        '9600gt': {'Label': gt9600label, 'Show': False, 'Marker': 'o',
                   'MarkerSize': msize, 'Color': 'g', 'norm': None,
                   'cmap': None, 'vmin': None, 'vmax': None, 'Alpha': 0.5,
                   'LineWidths': 0.1, 'Verts': None}}



gmax = 0.0
for gpu in gpus:
    print gpus[gpu]['Label'] + ': ' + str(gpus[gpu]['Show'])
    if(gpus[gpu]['Show']):
        data = np.loadtxt(str(gpu)+'_single_walker.dat')
        a = data[:,0].astype(np.uint32)
        a = a[plotoffset::plotstep]
        if logplot:
            b = data[:,1]*(10**9)
            #b = np.log10(data[:,1]*(10**9))
        else :
            b = data[:,1]*(10**9)
        if b.max() > gmax:
            gmax = b.max()
        b = b[plotoffset::plotstep]
        pl.scatter(a,b,
                   s=gpus[gpu]['MarkerSize'],
                   color=gpus[gpu]['Color'],
                   marker=gpus[gpu]['Marker'],
                   label=gpus[gpu]['Label'],
                   cmap=gpus[gpu]['cmap'],
                   norm=gpus[gpu]['norm'],
                   vmin=gpus[gpu]['vmin'],
                   vmax=gpus[gpu]['vmax'],
                   alpha=gpus[gpu]['Alpha'],
                   linewidths=gpus[gpu]['LineWidths'],
                   faceted=True,
                   verts=gpus[gpu]['Verts'],
                   hold=None)


# CPU results
cpu_val = 0.1*(10**9)

# CPU plot
pl.plot(np.linspace(0,32*2999,150),
        np.ones(150)*cpu_val,'r--',
        linewidth=3.0,
        label=cpulabel)

pl.xlabel(r'Number of Threads')
pl.ylabel(r'Metropolis steps [sec$^{-1}$]')

# Plot horizontal line at the maximum result
pl.axhline(gmax,alpha=0.1)

# Plot vertical lines on occupancy sweet spots
pl.axvline(x=448*32,alpha=0.25,linestyle='--')
pl.axvline(x=448*32*2,alpha=0.25,linestyle='--')
pl.axvline(x=448*32*3,alpha=0.25,linestyle='--')
pl.axvline(x=448*32*4,alpha=0.25,linestyle='--')
pl.axvline(x=448*32*5,alpha=0.25,linestyle='--')
pl.axvline(x=448*32*6,alpha=0.25,linestyle='--')

if logplot:
    pl.axis([0.0, float(xr)*10000.0, 2*10.0**7, 10.0**10])
else:
    pl.axis([0.0, float(xr)*10000.0, 0.0, 6.0*10**9])

# Fix the x axis to display 3.0*10^4 instead of 30000.0 etc
pl.xticks(np.arange(np.ceil(xr).astype(np.uint32)+1)*10**4,
          np.arange(np.ceil(xr).astype(np.uint32)+1).astype(np.float32))
plt.text(float(xr + 0.1)*10000.0,
         2*10.0**7,r'$\times10^4$',
         fontsize=0.8*fsize)

# Draw legend
pl.legend(loc=legendloc)

pl.show()
