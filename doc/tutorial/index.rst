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

    >>> from FeynSimul.kernel import PIMCKernel
    >>> from FeynSimul.kernel_args import KernelArgs
    >>> from FeynSimul.physical_systems.harm_osc import HarmOsc


Setting the parameters in a run
--------------------------------

What we first do is create an instance of the :class:`KernelArgs` which will
hold the parameters for the simulation. This class does not contain any methods
and is just for convenience. The class initiatior contains quite a lot of
parameters and what these do and which ones are needed to define can be read
about in the documentation for :ref:`kernelargs`.

    >>> phys_system = HarmOsc()
    >>> ka = KernelArgs(system = phys_system,
                        nbrOfWalkers = 32 * 28,
                        N = 128,
                        beta = 11.0,
                        S = 6,
                        enableOperator = True,
                        enableCorrelator = False,
                        metroStepsPerOperatorRun = 40,
                        enableBisection = True,
                        enablePathShift = False,
                        enableSingleNodeMove = False,
                        enableParallelizePath = True,
                        enableGlobalPath = False,
                        enableGlobalOldPath = False,
                        enableBins = False,
                        operatorRuns = 150,
                        nbrOfWalkersPerWorkGroup = 4,
                        operators = (phys_system.energyOp,))


When the kernel parameters are set the PIMCKernel class can be created with the
KernelArgs instance as argument.

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
    >>> ka.binResolutionPerDOF = 80


Now after we recreate the PIMCkernel object to update it with the new setttings
and run it, the probability density can be fetched.

    >>> kernel = PIMCKernel(ka)
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

   >>> ka.enableCorrelator = True
   >>> ka.correlators = ("x1",)
   >>> kernel = PIMCKernel(ka)
   >>> kernel.run()

We can now get the correlation like with the probability density. getCorrelator
returns a tuple with the correlator means and the standard error. For now we
are only interested in the means so we extract the first index.

    >>> corrs = kernel.getCorrelator()[0]

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



