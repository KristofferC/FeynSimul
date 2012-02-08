Tutorial
========

This tutorial will go through all the functionality in FeynSimul.
First 

Import FeynSimul and a physical system
--------------------------------------

This imports the kernel module that contains the class PIMCKernel class that is
the heart of the computations. The kernel_args module contains an empty classed
called kernelArgs which is used to hold the arguments to the PIMCKernel class.
A physical system is also imported. In this case this is a simple harmonic
oscillator, that is a particle in a :math:`V(x) = kx/2` potential.

    >>> from FeynSimul.kernel import *
    >>> from FeynSimul.kernel_args import *
    >>> from FeynSimul.physical_systems.harm_osc import *


Setting the parameters in a run
--------------------------------

This creates the instance of the class that will hold the parameters for our
run. This class does not contain any methods and is just for convenience.

    >>> ka = KernelArgs()

   
We set the system to be a harmonic oscillator which is imported from the
harm_osc module. 

    >>> ka.system = HarmOsc()

    
This determines how many independent simulations to be run at the same time and
how many of these should be in a separate work group. Variation of these
parameter might increase (or decrease) performance.

    >>> ka.nbrOfWalkers = 448 * 2
    >>> ka.nbrOfWalkersPerWorkGroup = 4


N determines the number of discretization points to use in a path
and beta can be regarded as the inverse temperature.    

    >>> ka.N = 128
    >>> ka.beta = 11.0


    
This determines how many times the sampling algorithm should be used before a
measurement with the operator should be done and also how many operator runs
should be done before the kernel exits.

    >>> ka.metroStepsPerOperatorRun = 40
    >>> ka.operatorRuns = 150

   
Here we set what types of samplings to use. Bisection is the fastest by far.
The others are there for academic purposes. 

    >>> ka.enableBisection = True
    >>> ka.enablePathShift = False
    >>> ka.enableSingleNodeMove = False

    
 
These are parameters to the bisection sampling algorithm. S determines the
width of a bisection interval. enableParallelizePath determines if multiple
threads should work on the same path simultaneously. enableGlobalPath
determines where to store the whole path. If N is too large this might have to
be enabled. enableGlobalOldPath determines where to store the  Similarly if S
is too large this need to be enabled. Using global memory results in slower
computations.   

    >>> ka.S = 6
    >>> ka.enableParallelizePath = True
    >>> ka.enableGlobalPath = False
    >>> ka.enableGlobalOldPath = False

  
For now we only want to do calculations involving ground state
operators so we disable the other.

    >>> ka.enableOperator = True
    >>> ka.enableCorrelator = False
    >>> ka.enableBins = False


Here we set what operators to use. Our physical system already defines the
operator that gives the ground state energy so we use it here. The operators
need to be in a tuple. 

    >>> ka.operators = (ka.system.energyOp,)


When all kernel parameters are set the PIMCKernel class can be created.

    >>> kernel = PIMCKernel(ka)



Running the kernel and retrieving values
----------------------------------------


To run the simulation one simple uses the run method in the PIMCkernel object.

    >>> kernel.run()

We can now fetch the mean of the operators for each independent simulation and calculate for
example the mean of them using numpy.

    >>> import numpy as np
    >>> operatorValues = kernel.getOperators()
    >>> print np.mean(operatorValues)
    0.499803373124

This is fairly close to the analytical result 0.5.


Calculating the probability density
-----------------------------------

If we want the calculations to include the probability density of the particle
in its ground state we need to enable it and rebuild the kernel with the
updated settings

    >>> ka.enableBins = True

Some parameters need to be set when the kernel is run with enableBins active.
These determines what interval to calculate the density function in and how fine
resolution should be used. Points outside the interval will be discarded.

    >>> ka.xMin = -3.5
    >>> ka.xMax = 3.5
    >>> binResolutionPerDOF = 80


Now after we recreate the PIMCkernel object to update it with the new setttings
and run it, the probability density can be fetched.

    >>> kernel = PIMCkernel(ka)
    >>> kernel.run()
    >>> probDensity = kernel.getBinCounts()

We could now plot it and compare it with the analytical solution:

    
    >>> import pylab as pl
    >>> binSize = (ka.xMax - ka.xMin) / ka.binResolutionPerDOF
    >>> x = np.linspace(ka.xMin, ka.xMax - binSize,
                ka.binResolutionPerDOF) + 0.5 * binSize
    >>> pl.plot(x, probDensity, '*', label="Simulated")
    >>> pl.plot(x, 1/np.sqrt(np.pi) * np.exp(-x ** 2), label="Analytical")
    >>> pl.legend(loc="best")
    >>> pl.show()

This would give this figure where the calculated and analytical result agree
well:

.. image:: prob_dens.*


    
Using correlator function to get excited energy state
-----------------------------------------------------

We can get the energy for the first excited state by exploiting calculating the
autocorrelation of the values of the nodes in the path. Enabling the correlator
calculation will calculate the autocorrelation for lags up to N / 2.

   >>> ka.enableCorrelators = True
   >>> ka.correlators = ("x1",)
   >>> kernel = PIMCkernel(ka)
   >>> kernel.run()

We can now get the correlation like with the probability density

    >>> corrs = kernel.getCorrelator()

The negative log derivative of the correlation gives the difference in energy
between the ground state and first excited state. We can plot this

    >>> import pylab as pl
    >>> logDerCorr = -np.gradient(np.log(corrs[0]), ka.beta / ka.N)
    >>> pl.plot(logDerCorr, '*')
    >>> pl.show()

This would give the following figue in which we can read that the different in
energies are about 1.0 which is correct.

.. image:: corr.*

The modN function
-----------------



