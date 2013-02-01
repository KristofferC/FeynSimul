import numpy as np
import pylab as pl
from matplotlib import rc
from matplotlib.ticker import *
import matplotlib.pyplot as plt


logplot = True   # Log the x axis or not
showdouble = False

#xor = np.array([0.000835, 0.021259, 0.077566, 2.113641, 7.748644, 213.113910]) #xorshift single
xor = np.array([0.000835, 0.021259, 0.077566, 2.113641, 7.748644]) #xorshift single

ranlux = np.array([[0.002772, 0.013913, 0.124960, 1.234409, 12.359268],   #ranlux0 single
                   [0.003888, 0.024675, 0.233399, 2.315815, 23.163764],   #ranlux0 double
               
                   [0.003372, 0.019870, 0.184641, 1.833603, 18.334329],   #ranlux1 single
                   [0.005065, 0.036373, 0.349675, 3.486245, 34.851969],   #ranlux1 double
               
                   [0.005217, 0.037840, 0.364948, 3.635847, 36.360428],   #ranlux2 single
                   [0.008705, 0.072788, 0.715689, 7.143060, 71.408121],   #ranlux2 double
               
                   [0.008220, 0.068618, 0.674307, 6.730529, 67.230933],   #ranlux3 single
                   [0.014878, 0.134723, 1.334791, 13.330784, 133.259846],  #ranlux3 double
               
                   [0.012606, 0.113489, 1.124352, 11.224786, 112.292905],  #ranlux4 single
                   [0.023882, 0.225516, 2.243392, 22.421897, 224.170527]]) #ranlux4 double
"""ranlux = np.array([[0.002772, 0.013913, 0.124960, 1.234409, 12.359268, ],   #ranlux0 single
                   [0.003888, 0.024675, 0.233399, 2.315815, 23.163764, ],   #ranlux0 double
               
                   [0.003372, 0.019870, 0.184641, 1.833603, 18.334329,],   #ranlux1 single
                   [0.005065, 0.036373, 0.349675, 3.486245, 34.851969,],   #ranlux1 double
               
                   [0.005217, 0.037840, 0.364948, 3.635847, 36.360428,],   #ranlux2 single
                   [0.008705, 0.072788, 0.715689, 7.143060, 71.408121,],   #ranlux2 double
               
                   [0.008220, 0.068618, 0.674307, 6.730529, 67.230933,],   #ranlux3 single
                   [0.014878, 0.134723, 1.334791, 13.330784, 133.259846,],  #ranlux3 double
               
                   [0.012606, 0.113489, 1.124352, 11.224786, 112.292905, 1122.878422],  #ranlux4 single
                   [0.023882, 0.225516, 2.243392, 22.421897, 224.170527, ]]) #ranlux4 double
"""

# 10 000 000 * nbrOfThreads:
# xor:  213.113910 
# ran0: 
# ran1: 
# ran2: 
# ran3: 
# ran4: 1122.878422


nbrOfThreads = 448*32*4;
randsteps = np.array([1000,10000,100000,1000000,10000000])
#randsteps = np.array([1000,10000,100000,1000000,10000000,100000000])


msize = 10
fig_width_pt = 512.0
legendloc = 4
fsize=16

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
         'xtick.labelsize': 1*fsize,
         'ytick.labelsize': 1*fsize,
         'text.usetex': True,
         'figure.figsize': fig_size}

pl.rcParams.update(params)

pl.figure(1)
if logplot:
    pl.semilogx()

# xorshift plot
pl.plot(randsteps,
        randsteps/xor,'r-',
        linewidth=3.0,
        label=('xorshift, 32bit'),
        markersize=msize,
        marker='o')


if showdouble:
    linestyle = ('g--','b--','c--','k--','m--')
    markerstyle = ('^','<','>','v','*')
    for i in range(0,5):
        pl.plot(randsteps,
                randsteps/ranlux[i*2+1],
                linestyle[i],
                linewidth=2.0,
                label=('RANLUX, lux='+str(i)+', 64bit'),
                markersize=msize,
                marker=markerstyle[i])
else:
    linestyle = ('g-','b-','c-','k-','m-')
    markerstyle = ('^','<','>','v','*')
    for i in range(0,5):
        pl.plot(randsteps,
                randsteps/ranlux[i*2],
                linestyle[i],
                linewidth=2.0,
                label=('RANLUX, lux='+str(i)+', 32bit'),
                markersize=msize,
                marker=markerstyle[i])


if logplot:
    pl.title(r'logplot of xorshift vs RANLUX on 448*32*4 (=57344) threads')
else:
    pl.title(r'plot xorshift vs RANLUX on 448*32*4 (=57344) threads')
pl.xlabel(r'Number of randoms per thread')
pl.ylabel(r'Randoms on each thread per second [sec$^{-1}$]')

#        'ranlux4': {'Label': gt9600label, 'Show': False, 'Marker': 'o',
#                   'MarkerSize': msize, 'Color': 'g', 'norm': None,
#                   'cmap': None, 'vmin': None, 'vmax': None, 'Alpha': 0.5,
#                   'LineWidths': 0.1, 'Verts': None}}


"""
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
"""
# Draw legend
pl.legend(loc=legendloc)

pl.show()
